#!/usr/bin/env python3
"""Bounded prospective causal test at physical approach-storage onset.

Select the first majority-source onset in each of 14 seeds using current/past
storage only.  Restart four arms: original/symmetrized internal alignment,
each with production angular noise or zero post-restart angular noise.  Stop
at the established 0.5-grid return of total approach storage to zero.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_nonreciprocal_alignment_intervention as nr
from analyze_approach_current_balance import ELL
from analyze_individual_heading_geometry import angular_apertures, global_required_turn, wrap
from analyze_initiation_propagation import RawReader


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO.parent / "BEM_bichamber_parallel_instance_no_hydro.py"
RAW_DEFAULT = Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
FEATURES = HERE / "outputs/initiation_propagation/cache/population_kinetics_features_dt0p5.npz"
OUT_DEFAULT = HERE / "outputs/prospective_entry_intervention"
N_STATES = 14
N_CLONES = 8
SAMPLE_STEPS = 500
MAX_STEPS = 20000
ZERO_TOL = 1e-14


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def load_features(path: Path) -> dict[int, dict[str, np.ndarray]]:
    with np.load(path) as z:
        arrays = {key: np.asarray(z[key]) for key in z.files}
    out = {}
    for seed in range(40):
        take = arrays["seed"] == seed; order = np.argsort(arrays["index"][take])
        out[seed] = {key: value[take][order] for key, value in arrays.items()
                     if key not in ("seed", "index")}
    return out


def select_states(path: Path) -> list[dict]:
    features = load_features(path); selected = []
    for seed in range(N_STATES):
        a = features[seed]; total = a["b_left"] + a["b_right"]
        candidates = np.flatnonzero((total > ZERO_TOL) & np.r_[False, total[:-1] <= ZERO_TOL])
        chosen = None
        for index in candidates:
            direction = 1 if a["b_left"][index] > a["b_right"][index] else -1
            q_past = float(a["q"][index - 1])
            majority_direction = 1 if q_past < 0 else (-1 if q_past > 0 else 0)
            if direction == majority_direction:
                chosen = (int(index), direction, q_past, float(a["b_left"][index]), float(a["b_right"][index]))
                break
        if chosen is None:
            raise RuntimeError(f"No prospective majority-source onset for seed {seed}")
        index, direction, q_past, bl, br = chosen
        native_step = round((200 + .5 * index) / nr.DT)
        selected.append(dict(
            sample_id=f"onset_s{seed}_i{index}", seed=seed, onset_index=index,
            native_step=native_step, direction=direction, source_sign=-direction,
            q_past=q_past, b_left_onset=bl, b_right_onset=br,
            selection_rule="first majority-source storage onset after transient; current and previous 0.5-grid samples only",
        ))
    return selected


def setup_origins(raw_dir: Path, step, para: dict, selected: list[dict]) -> list[dict]:
    zero = np.zeros((3, nr.N)); out = []
    for item in selected:
        reader = RawReader(nr.raw_path(raw_dir, item["seed"])); k = item["native_step"]
        previous = reader.frame("state", k - 1); state = reader.frame("state", k); target = reader.frame("state", k + 1); reader.close()
        _, _, _, memory, _ = step(previous, para, zero, None, item["source_sign"], "baseline")
        aperture = angular_apertures(state, item["source_sign"], float(nr.ELL))
        required_global = global_required_turn(aperture["minimum_turn_signed"], item["source_sign"])
        out.append(dict(
            item, state=state, dis_last=memory, archived_target=target,
            nright0=int(np.sum(state[0] >= 0)), original_indices=np.flatnonzero(nr.source_mask(state, item["source_sign"])),
            initial_heading=aperture["mean_angle"], initial_polarization=aperture["polarization"],
            minimum_turn_global=required_global, minimum_turn=abs(required_global),
            collective_noise_sigma=aperture["collective_angular_sigma"], fish_noise_sigma=.5,
        ))
    return out


def total_storage(state: np.ndarray) -> float:
    x, y = state[:2]; transverse = np.abs(y) <= ELL / 2
    left = (x < 0) & (x >= -ELL) & transverse
    right = (x >= 0) & (x <= ELL) & transverse
    return float((np.sum(ELL + x[left]) + np.sum(ELL - x[right])) / (nr.N * ELL))


def approach_current(state: np.ndarray, velocity: np.ndarray, direction: int) -> float:
    x, y = state[:2]
    inside = (direction * x < 0) & (-direction * x <= ELL) & (np.abs(y) <= ELL / 2)
    return float(np.sum(direction * velocity[0, inside]) / nr.N)


def original_heading(state: np.ndarray, indices: np.ndarray) -> float:
    return float(math.atan2(np.mean(np.sin(state[2, indices])), np.mean(np.cos(state[2, indices]))))


def audit(raw_dir: Path, mod, step, para: dict, origins: list[dict], adapter_sha: str) -> dict:
    rows = []; zero = np.zeros((3, nr.N))
    for item in origins:
        state = item["state"]; base, _, _, _, _ = step(state, para, zero, item["dis_last"], item["source_sign"], "baseline")
        increment = np.zeros((3, nr.N)); increment[2] = item["archived_target"][2] - base[2]
        predicted, velocity, _, _, _ = step(state, para, increment, item["dis_last"], item["source_sign"], "baseline")
        reader = RawReader(nr.raw_path(raw_dir, item["seed"])); archived_velocity = reader.frame("rdot", item["native_step"] + 1); reader.close()
        rows.append(dict(sample_id=item["sample_id"], max_state_error=float(np.max(np.abs(predicted - item["archived_target"]))),
                         max_velocity_error=float(np.max(np.abs(velocity - archived_velocity))), initial_storage=total_storage(state)))
    maximum = max(max(r["max_state_error"], r["max_velocity_error"]) for r in rows)
    return dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), adapter_sha256=adapter_sha,
                parity_passed=maximum < 1e-12, max_parity_error=maximum,
                all_initial_storage_positive=all(r["initial_storage"] > ZERO_TOL for r in rows), rows=rows)


def run_task(task: dict) -> list[dict]:
    rng = np.random.default_rng(task["noise_seed"])
    noise = .5 * math.sqrt(nr.DT) * rng.normal(size=(MAX_STEPS, nr.N))
    output = []
    for interaction in ("baseline", "symmetric_internal_alignment"):
        for noise_mode in ("production_noise", "zero_noise"):
            state = task["state"].copy(); memory = task["dis_last"].copy(); indices = task["original_indices"]
            theta_previous = original_heading(state, indices); unwrapped = 0.0; entered = False
            entry_time = float("nan"); heading_at_entry = float("nan")
            integrated_current = 0.0; forward = backward = 0; completed = False
            for n in range(MAX_STEPS):
                increment = np.zeros((3, nr.N))
                if noise_mode == "production_noise":
                    increment[2] = noise[n]
                previous = state
                state, velocity, _, memory, _ = nr._STEP(previous, nr._PARA, increment, memory, task["source_sign"], interaction)
                integrated_current += approach_current(previous, velocity, task["direction"]) * nr.DT
                f = (task["source_sign"] * previous[0] >= 0) & (task["source_sign"] * state[0] < 0)
                b = (task["source_sign"] * previous[0] < 0) & (task["source_sign"] * state[0] >= 0)
                forward += int(np.sum(f)); backward += int(np.sum(b))
                theta = original_heading(state, indices); unwrapped += float(wrap(theta - theta_previous)); theta_previous = theta
                if not entered and np.any(f):
                    entered = True; entry_time = (n + 1) * nr.DT; heading_at_entry = unwrapped
                if (n + 1) % SAMPLE_STEPS == 0 and total_storage(state) <= ZERO_TOL:
                    completed = True; steps = n + 1; break
            else:
                steps = MAX_STEPS
            duration = steps * nr.DT; nright = int(np.sum(state[0] >= 0))
            net = int(task["direction"] * (nright - task["nright0"]))
            required_sign = float(np.sign(task["minimum_turn_global"]))
            output.append(dict(
                sample_id=task["sample_id"], seed=task["seed"], clone=task["clone"],
                interaction=interaction, noise_mode=noise_mode, source_sign=task["source_sign"], direction=task["direction"],
                completed_excursion=int(completed), duration=duration, first_entry=int(entered), entry_time=entry_time,
                net_transfer=net, success5=int(net >= 5), forward_crossings=forward, backward_crossings=backward,
                mean_Ja=integrated_current / duration, integrated_Ja=integrated_current,
                heading_change_total=unwrapped, heading_change_to_entry=heading_at_entry,
                heading_toward_required_total=required_sign * unwrapped,
                heading_toward_required_to_entry=required_sign * heading_at_entry if entered else float("nan"),
                initial_minimum_turn=task["minimum_turn"], collective_noise_sigma=task["collective_noise_sigma"],
                fish_noise_sigma=task["fish_noise_sigma"], nright0=task["nright0"], nright_end=nright,
            ))
    return output


def simulate(origins: list[dict], adapter_text: str, output: Path, workers: int) -> list[dict]:
    tasks = []
    for item in origins:
        clean = {key: value for key, value in item.items() if key != "archived_target"}
        for clone in range(N_CLONES):
            tasks.append(dict(clean, clone=clone, noise_seed=202611230000 + item["seed"] * 100 + clone))
    partial = output / "branches_partial.csv"
    complete_file = output / "branches.csv"
    resume_file = partial if partial.exists() else complete_file
    existing = read_csv(resume_file) if resume_file.exists() else []
    complete = {(r["sample_id"], int(r["clone"])) for r in existing
                if sum(x["sample_id"] == r["sample_id"] and int(x["clone"]) == int(r["clone"]) for x in existing) == 4}
    rows = [r for r in existing if (r["sample_id"], int(r["clone"])) in complete]
    todo = [task for task in tasks if (task["sample_id"], task["clone"]) not in complete]
    start = time.monotonic()
    if todo:
        with mp.get_context("spawn").Pool(min(workers, len(todo)), initializer=nr.init_worker,
                                           initargs=(str(SOURCE), adapter_text)) as pool:
            for count, result in enumerate(pool.imap_unordered(run_task, todo), 1):
                rows.extend(result); write_csv(partial, rows)
                print(f"completed {count}/{len(todo)} four-arm blocks; elapsed {time.monotonic()-start:.1f}s", flush=True)
    rows.sort(key=lambda r: (r["sample_id"], int(r["clone"]), r["interaction"], r["noise_mode"]))
    write_csv(output / "branches.csv", rows)
    if partial.exists():
        partial.unlink()
    return rows


def cluster_effect(rows: list[dict], field: str, arm_a: tuple[str, str], arm_b: tuple[str, str],
                   nboot: int = 3000) -> dict:
    groups = defaultdict(dict)
    for row in rows:
        groups[(int(row["seed"]), row["sample_id"], int(row["clone"]))][(row["interaction"], row["noise_mode"])] = row
    paired = []
    for (seed, sample, clone), group in groups.items():
        if arm_a not in group or arm_b not in group:
            continue
        a = float(group[arm_a][field]); b = float(group[arm_b][field])
        if np.isfinite(a) and np.isfinite(b):
            paired.append(dict(seed=seed, value=b - a))
    by = defaultdict(list)
    for row in paired:
        by[row["seed"]].append(row["value"])
    seeds = sorted(by); seed_values = np.asarray([np.mean(by[s]) for s in seeds]); point = float(np.mean(seed_values))
    rng = np.random.default_rng(20261124); draw = rng.choice(seed_values, (nboot, len(seed_values)), replace=True).mean(axis=1)
    return dict(arm_a="/".join(arm_a), arm_b="/".join(arm_b), effect_b_minus_a=point,
                ci95=np.percentile(draw, [2.5, 97.5]).tolist(), n_pairs=len(paired), n_seeds=len(seeds))


def cluster_interaction(rows: list[dict], field: str, nboot: int = 3000) -> dict:
    """Seed-clustered interaction: symmetry effect with noise minus without it."""
    groups = defaultdict(dict)
    for row in rows:
        groups[(int(row["seed"]), row["sample_id"], int(row["clone"]))][
            (row["interaction"], row["noise_mode"])
        ] = row
    required = (
        ("baseline", "production_noise"),
        ("symmetric_internal_alignment", "production_noise"),
        ("baseline", "zero_noise"),
        ("symmetric_internal_alignment", "zero_noise"),
    )
    by = defaultdict(list)
    for (seed, _sample, _clone), group in groups.items():
        if not all(arm in group for arm in required):
            continue
        values = [float(group[arm][field]) for arm in required]
        if all(np.isfinite(values)):
            # (symmetric - baseline)_production - (symmetric - baseline)_zero
            by[seed].append((values[1] - values[0]) - (values[3] - values[2]))
    seeds = sorted(by)
    seed_values = np.asarray([np.mean(by[s]) for s in seeds])
    if not len(seed_values):
        return dict(effect_difference=float("nan"), ci95=[float("nan"), float("nan")],
                    n_pairs=0, n_seeds=0)
    rng = np.random.default_rng(20261126)
    draw = rng.choice(seed_values, (nboot, len(seed_values)), replace=True).mean(axis=1)
    return dict(effect_difference=float(np.mean(seed_values)),
                ci95=np.percentile(draw, [2.5, 97.5]).tolist(),
                n_pairs=sum(len(by[s]) for s in seeds), n_seeds=len(seeds))


def arm_summary(rows: list[dict]) -> list[dict]:
    output = []
    for interaction in ("baseline", "symmetric_internal_alignment"):
        for noise_mode in ("production_noise", "zero_noise"):
            selected = [r for r in rows if r["interaction"] == interaction and r["noise_mode"] == noise_mode]
            by = defaultdict(list)
            for row in selected:
                by[int(row["seed"])].append(row)
            for field in ("completed_excursion", "first_entry", "entry_time", "net_transfer", "success5", "mean_Ja", "integrated_Ja",
                          "heading_toward_required_total", "heading_toward_required_to_entry", "duration"):
                values = []
                for seed, part in by.items():
                    finite = [float(r[field]) for r in part if np.isfinite(float(r[field]))]
                    if finite:
                        values.append(float(np.mean(finite)))
                array = np.asarray(values); rng = np.random.default_rng(20261125)
                draw = rng.choice(array, (3000, len(array)), replace=True).mean(axis=1)
                output.append(dict(interaction=interaction, noise_mode=noise_mode, metric=field,
                                   estimate=float(np.mean(array)), ci95_low=float(np.percentile(draw, 2.5)),
                                   ci95_high=float(np.percentile(draw, 97.5)), n_rows=len(selected), n_seeds=len(array)))
    return output


def analyze(rows: list[dict], output: Path) -> dict:
    metrics = ("first_entry", "entry_time", "net_transfer", "success5", "mean_Ja", "integrated_Ja",
               "heading_toward_required_total", "heading_toward_required_to_entry", "duration")
    comparisons = []
    for noise_mode in ("production_noise", "zero_noise"):
        for field in metrics:
            comparisons.append(dict(contrast=f"symmetry_effect_at_{noise_mode}", metric=field,
                                    **cluster_effect(rows, field, ("baseline", noise_mode), ("symmetric_internal_alignment", noise_mode))))
    for interaction in ("baseline", "symmetric_internal_alignment"):
        for field in metrics:
            comparisons.append(dict(contrast=f"noise_suppression_effect_at_{interaction}", metric=field,
                                    **cluster_effect(rows, field, (interaction, "production_noise"), (interaction, "zero_noise"))))
    for field in metrics:
        comparisons.append(dict(contrast="symmetry_by_noise_interaction", metric=field,
                                **cluster_interaction(rows, field)))
    summaries = arm_summary(rows); write_csv(output / "arm_summary.csv", summaries); write_csv(output / "paired_effects.csv", comparisons)
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    labels = ["orig\nnoise", "sym\nnoise", "orig\nzero", "sym\nzero"]
    arms = [("baseline", "production_noise"), ("symmetric_internal_alignment", "production_noise"),
            ("baseline", "zero_noise"), ("symmetric_internal_alignment", "zero_noise")]
    summary_lookup = {(r["interaction"], r["noise_mode"], r["metric"]): r for r in summaries}
    panels = (("first_entry", "First-entry probability"),
              ("net_transfer", "Excursion-end net transfer"),
              ("integrated_Ja", "Integrated approach current"),
              ("heading_toward_required_total", "Collective turn toward initial aperture"))
    for ax, (field, title) in zip(axes.flat, panels):
        values = np.asarray([float(summary_lookup[(*arm, field)]["estimate"]) for arm in arms])
        low = np.asarray([float(summary_lookup[(*arm, field)]["ci95_low"]) for arm in arms])
        high = np.asarray([float(summary_lookup[(*arm, field)]["ci95_high"]) for arm in arms])
        ax.errorbar(np.arange(4), values, yerr=np.vstack((values - low, high - values)),
                    fmt="o", capsize=4, color="#286090")
        ax.set_xticks(np.arange(4), labels); ax.set_title(title); ax.axhline(0, color="k", ls=":", lw=.8)
    fig.tight_layout(); fig.savefig(output / "prospective_entry_intervention.png", dpi=220); plt.close(fig)
    initial = {(r["sample_id"], float(r["collective_noise_sigma"]), float(r["fish_noise_sigma"])) for r in rows}
    ratios = [fish / collective for _, collective, fish in initial]
    population_identity = all(int(r["net_transfer"]) == int(r["forward_crossings"]) - int(r["backward_crossings"])
                              for r in rows)
    zero_deterministic = True
    for sample in {r["sample_id"] for r in rows}:
        for interaction in ("baseline", "symmetric_internal_alignment"):
            part = [r for r in rows if r["sample_id"] == sample and r["interaction"] == interaction
                    and r["noise_mode"] == "zero_noise"]
            for field in ("first_entry", "entry_time", "net_transfer", "forward_crossings", "backward_crossings",
                          "duration", "mean_Ja", "heading_change_total"):
                values = np.asarray([float(r[field]) for r in part])
                finite = values[np.isfinite(values)]
                if len(finite) and np.ptp(finite) > 1e-12:
                    zero_deterministic = False
    return dict(n_rows=len(rows), n_states=len({r["sample_id"] for r in rows}), n_clones=N_CLONES,
                censored_rows=sum(not int(r["completed_excursion"]) for r in rows),
                median_collective_noise_sigma=float(np.median([x[1] for x in initial])),
                fish_noise_sigma=.5, median_fish_to_collective_noise_ratio=float(np.median(ratios)),
                population_crossing_identity_passed=population_identity,
                zero_noise_clone_determinism_passed=zero_deterministic,
                comparisons=comparisons)


def self_test() -> None:
    state = np.zeros((3, 4)); state[0] = [-1, -.5, .5, 1]; state[1] = [0, .1, 0, 3]
    assert total_storage(state) > 0
    state[0] = [-3, -3, 3, 3]; state[1] = [3, 3, 3, 3]
    assert total_storage(state) == 0
    print("PASS: exact storage-onset reset surface")


def main() -> None:
    p = argparse.ArgumentParser(); p.add_argument("--raw-dir", type=Path, default=RAW_DEFAULT)
    p.add_argument("--features", type=Path, default=FEATURES); p.add_argument("--output", type=Path, default=OUT_DEFAULT)
    p.add_argument("--workers", type=int, default=max(1, min(12, os.cpu_count() or 1)))
    p.add_argument("--stage", choices=("audit", "simulate", "all"), default="all"); p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test(); return
    args.output.mkdir(parents=True, exist_ok=True)
    mod = nr.load_gap20_production(SOURCE); para = nr.parameters(mod); step, text, adapter_sha = nr.extract_step(mod, SOURCE)
    selected = select_states(args.features); origins = setup_origins(args.raw_dir, step, para, selected)
    write_csv(args.output / "selection_manifest.csv", [{k: v for k, v in item.items() if k not in ("state", "dis_last", "archived_target", "original_indices")} for item in origins])
    audit_result = audit(args.raw_dir, mod, step, para, origins, adapter_sha)
    (args.output / "restart_audit.json").write_text(json.dumps(audit_result, indent=2) + "\n")
    if not audit_result["parity_passed"] or not audit_result["all_initial_storage_positive"]:
        raise RuntimeError("Prospective restart audit failed")
    print(json.dumps({k: audit_result[k] for k in ("parity_passed", "max_parity_error", "all_initial_storage_positive")}, indent=2))
    if args.stage == "audit":
        return
    rows = simulate(origins, text, args.output, args.workers); result = analyze(rows, args.output)
    metadata = dict(source_sha256=audit_result["source_sha256"], adapter_sha256=adapter_sha,
                    selection="first majority-source storage onset from current/past samples only",
                    arms="original/symmetric internal alignment x production/zero post-restart angular noise",
                    excursion_end="first 0.5-grid sample with total weighted approach storage exactly zero",
                    max_time=MAX_STEPS * nr.DT, paired_noise=True, new_full_trajectories=0, **result)
    (args.output / "summary.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({k: metadata[k] for k in ("n_rows", "n_states", "n_clones", "censored_rows", "median_collective_noise_sigma")}, indent=2))


if __name__ == "__main__":
    main()
