#!/usr/bin/env python3
"""Test a parameter-free coherent-school ballistic encounter prediction.

The predictor rigidly translates the actual source-school positions with the
instantaneous mean self-propulsion velocity.  It uses no outcome, fitted
coefficient, or added collective coordinate.  Existing storage-excursion
attempts, exact event budgets, torque decompositions, and matched-noise
intervention outputs are read only.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from analyze_approach_current_balance import ELL, N, load_gap20_production, raw_path
from analyze_initiation_propagation import RawReader, auc_score


DT_NATIVE = 0.001
ATTEMPT_DT = 0.5
ENTRY_HORIZON = 2.0
SUSTAINED_HORIZON = 5.0
TRANSFER_THRESHOLD = 5


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


def source_mask(state: np.ndarray, source_sign: int) -> np.ndarray:
    return state[0] >= 0 if source_sign > 0 else state[0] < 0


def segment_intersection(a: np.ndarray, b: np.ndarray,
                         c: np.ndarray, d: np.ndarray) -> bool:
    """Strict proper intersection, matching the production visibility rule."""
    def ccw(p, q, r):
        return (r[1] - p[1]) * (q[0] - p[0]) > (q[1] - p[1]) * (r[0] - p[0])
    return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)


def circle_intersection(a: np.ndarray, b: np.ndarray,
                        center: np.ndarray, radius: float) -> bool:
    delta = b - a; offset = a - center
    aa = float(delta @ delta)
    if aa == 0:
        return False
    bb = 2 * float(offset @ delta); cc = float(offset @ offset) - radius**2
    disc = bb * bb - 4 * aa * cc
    if disc < 0:
        return False
    root = math.sqrt(disc)
    t1 = (-bb - root) / (2 * aa); t2 = (-bb + root) / (2 * aa)
    # Tangency at the physical aperture endpoint is allowed.
    eps = 1e-10
    return (eps < t1 < 1 - eps) or (eps < t2 < 1 - eps)


def path_blocked(a: np.ndarray, b: np.ndarray, geometry: dict) -> bool:
    for wall in geometry["segments"]:
        if segment_intersection(a, b, wall[0], wall[1]):
            return True
    return any(circle_intersection(a, b, center, geometry["tip_radius"])
               for center in geometry["tip_centers"])


def production_geometry(mod) -> dict:
    return dict(
        gap=float(mod.gap),
        segments=np.asarray(mod.segments, dtype=float),
        tip_centers=np.asarray([mod.circle_center5, mod.circle_center6], dtype=float),
        tip_radius=float(mod.r),
        boundary_x=np.asarray(mod.x_coord, dtype=float),
        boundary_y=np.asarray(mod.y_coord, dtype=float),
    )


def ballistic_prediction(state: np.ndarray, source_sign: int,
                         horizon: float, geometry: dict) -> dict:
    """Rigid translation by v0 times the source mean-heading vector.

    v0=1 in production units.  The mean vector magnitude therefore supplies
    the actual coherent COM speed (apart from collision overrides), rather
    than imposing unit COM speed.
    """
    mask = source_mask(state, source_sign)
    theta = state[2, mask]
    velocity = np.asarray([np.mean(np.cos(theta)), np.mean(np.sin(theta))])
    polarization = float(np.linalg.norm(velocity))
    direction = -source_sign
    inward_speed = float(direction * velocity[0])
    indices = np.flatnonzero(mask)
    arrivals = []
    obstructed = 0
    if inward_speed > 0:
        distance = source_sign * state[0, indices]
        times = distance / inward_speed
        y_hit = state[1, indices] + velocity[1] * times
        candidate = (times > 0) & (times <= horizon + 1e-12) & (np.abs(y_hit) <= geometry["gap"] / 2 + 1e-12)
        for local in np.flatnonzero(candidate):
            start = state[:2, indices[local]]
            endpoint = np.asarray([0.0, y_hit[local]])
            if path_blocked(start, endpoint, geometry):
                obstructed += 1
            else:
                arrivals.append(float(times[local]))
    return dict(
        n_ballistic=len(arrivals),
        first_ballistic_time=min(arrivals) if arrivals else float("nan"),
        mean_ballistic_time=float(np.mean(arrivals)) if arrivals else float("nan"),
        polarization=polarization,
        mean_vx=float(velocity[0]), mean_vy=float(velocity[1]),
        inward_speed=inward_speed,
        n_source=int(mask.sum()), n_geometry_obstructed=obstructed,
    )


def approach_state(state: np.ndarray, stored_velocity: np.ndarray,
                   direction: int) -> dict:
    x, y = state[:2]
    inside = (direction * x < 0) & (-direction * x <= ELL) & (np.abs(y) <= ELL / 2)
    distance = -direction * x[inside]
    storage = float(np.sum(ELL - distance) / (N * ELL))
    current = float(np.sum(direction * stored_velocity[0, inside]) / N)
    return dict(Ja=current, b=storage, n_approach=int(inside.sum()))


def load_feature_cache(path: Path) -> dict[int, dict[str, np.ndarray]]:
    with np.load(path) as z:
        arrays = {key: np.asarray(z[key]) for key in z.files}
    out = {}
    for seed in range(40):
        take = arrays["seed"] == seed; order = np.argsort(arrays["index"][take])
        out[seed] = {key: value[take][order] for key, value in arrays.items()
                     if key not in ("seed", "index")}
    return out


def extract_attempts(raw_dir: Path, attempts_path: Path, cache_path: Path,
                     geometry: dict) -> tuple[list[dict], dict]:
    attempts = [r for r in read_csv(attempts_path) if int(r["majority_source"]) == 1]
    features = load_feature_cache(cache_path); by = defaultdict(list)
    for row in attempts:
        by[int(row["seed"])].append(row)
    output = []; obstructed = 0; source_count_mismatch = 0
    for count, seed in enumerate(sorted(by), 1):
        print(f"ballistic attempt snapshots {count}/{len(by)} seed={seed}", flush=True)
        reader = RawReader(raw_path(raw_dir, seed)); f = features[seed]
        try:
            for source in by[seed]:
                start = int(source["start_index"]); step = int(round(float(source["onset_upper"]) / DT_NATIVE))
                if start + round(SUSTAINED_HORIZON / ATTEMPT_DT) >= len(f["q"]):
                    continue
                state = reader.frame("state", step); velocity = reader.frame("rdot", step + 1)
                direction = int(source["direction"]); source_sign = -direction
                p2 = ballistic_prediction(state, source_sign, ENTRY_HORIZON, geometry)
                p5 = ballistic_prediction(state, source_sign, SUSTAINED_HORIZON, geometry)
                q0 = float(f["q"][start]); q2 = float(f["q"][start + round(ENTRY_HORIZON / ATTEMPT_DT)])
                q5 = float(f["q"][start + round(SUSTAINED_HORIZON / ATTEMPT_DT)])
                net2 = int(round(direction * N * (q2 - q0) / 2))
                net5 = int(round(direction * N * (q5 - q0) / 2))
                entry_time = float(source["first_entry_time_upper"])
                entry2 = int(source["first_entry"] == "1" and np.isfinite(entry_time) and
                             entry_time <= float(source["onset_upper"]) + ENTRY_HORIZON + 1e-12)
                ap = approach_state(state, velocity, direction)
                expected_source = int(round(N * (1 + source_sign * q0) / 2))
                source_count_mismatch += int(expected_source != p2["n_source"])
                obstructed += p5["n_geometry_obstructed"]
                output.append(dict(
                    seed=seed, attempt_id=int(source["attempt_id"]), onset_time=float(source["onset_upper"]),
                    direction=direction, source_sign=source_sign, q0=q0, abs_q0=abs(q0),
                    entry2=entry2, net2=net2, positive_transfer2=int(net2 > 0),
                    net5=net5, sustained5=int(net5 >= TRANSFER_THRESHOLD),
                    n_ballistic2=p2["n_ballistic"], n_ballistic5=p5["n_ballistic"],
                    first_ballistic_time2=p2["first_ballistic_time"],
                    ballistic_entry2=int(p2["n_ballistic"] > 0),
                    ballistic_sustained5=int(p5["n_ballistic"] >= TRANSFER_THRESHOLD),
                    polarization=p2["polarization"], inward_speed=p2["inward_speed"],
                    mean_vx=p2["mean_vx"], mean_vy=p2["mean_vy"], n_source=p2["n_source"],
                    **ap,
                ))
        finally:
            reader.close()
    audit = dict(n_majority_source_attempts_input=len(attempts), n_rows=len(output),
                 n_source_count_mismatch=source_count_mismatch,
                 n_solid_wall_obstructions=obstructed)
    return output, audit


def cluster_boot(rows: list[dict], metric, seed: int, nboot: int = 2000) -> dict:
    by = defaultdict(list)
    for row in rows:
        by[int(row["seed"])].append(row)
    seeds = sorted(by); point = float(metric(rows)); rng = np.random.default_rng(seed); draws = []
    for _ in range(nboot):
        sample = []
        for replicate, chosen in enumerate(rng.choice(seeds, len(seeds), replace=True)):
            for row in by[int(chosen)]:
                copied = dict(row); copied["seed"] = replicate; sample.append(copied)
        value = metric(sample)
        if np.isfinite(value):
            draws.append(float(value))
    return dict(estimate=point, ci95=np.percentile(draws, [2.5, 97.5]).tolist(),
                n=len(rows), n_seeds=len(seeds))


def auc(rows: list[dict], score: str, target: str) -> float:
    return float(auc_score(np.asarray([int(r[target]) for r in rows]),
                           np.asarray([float(r[score]) for r in rows])))


def corr(rows: list[dict], x: str, y: str) -> float:
    a = np.asarray([float(r[x]) for r in rows]); b = np.asarray([float(r[y]) for r in rows])
    return float(spearmanr(a, b).statistic) if np.std(a) and np.std(b) else float("nan")


def confusion(rows: list[dict], prediction: str, target: str) -> dict:
    p = np.asarray([int(r[prediction]) for r in rows], dtype=bool)
    y = np.asarray([int(r[target]) for r in rows], dtype=bool)
    sensitivity = float(np.mean(p[y])) if np.any(y) else float("nan")
    specificity = float(np.mean(~p[~y])) if np.any(~y) else float("nan")
    return dict(sensitivity=sensitivity, specificity=specificity,
                balanced_accuracy=(sensitivity + specificity) / 2,
                predicted_positive_fraction=float(np.mean(p)), observed_positive_fraction=float(np.mean(y)))


def summarize_attempts(rows: list[dict]) -> list[dict]:
    summaries = []
    for cohort, selected in (("train", [r for r in rows if r["seed"] < 30]),
                             ("heldout", [r for r in rows if r["seed"] >= 30])):
        for target, score, comparator in (("entry2", "n_ballistic2", "Ja"),
                                           ("positive_transfer2", "n_ballistic2", "Ja"),
                                           ("sustained5", "n_ballistic5", "Ja")):
            ab = cluster_boot(selected, lambda x, s=score, t=target: auc(x, s, t), 20261160)
            aj = cluster_boot(selected, lambda x, s=comparator, t=target: auc(x, s, t), 20261161)
            diff = cluster_boot(selected, lambda x, s=score, c=comparator, t=target: auc(x, s, t) - auc(x, c, t), 20261162)
            summaries.append(dict(cohort=cohort, target=target, metric="AUC", predictor="ballistic_count", **ab))
            summaries.append(dict(cohort=cohort, target=target, metric="AUC", predictor="instantaneous_Ja", **aj))
            summaries.append(dict(cohort=cohort, target=target, metric="AUC_difference", predictor="ballistic_minus_Ja", **diff))
        for target, score in (("net2", "n_ballistic2"), ("net5", "n_ballistic5")):
            cr = cluster_boot(selected, lambda x, s=score, t=target: corr(x, s, t), 20261163)
            bias = cluster_boot(selected, lambda x, s=score, t=target: float(np.mean([r[s] - r[t] for r in x])), 20261164)
            summaries.append(dict(cohort=cohort, target=target, metric="Spearman", predictor="ballistic_count", **cr))
            summaries.append(dict(cohort=cohort, target=target, metric="prediction_minus_actual", predictor="ballistic_count", **bias))
        for prediction, target in (("ballistic_entry2", "entry2"), ("ballistic_sustained5", "sustained5")):
            c = confusion(selected, prediction, target)
            boot = cluster_boot(selected, lambda x, p=prediction, t=target: confusion(x, p, t)["balanced_accuracy"], 20261165)
            summaries.append(dict(cohort=cohort, target=target, metric="parameter_free_classification", predictor=prediction,
                                  estimate=c["balanced_accuracy"], ci95=boot["ci95"], n=len(selected), n_seeds=len({r["seed"] for r in selected}), **c))
    return summaries


def aggregate_event_budgets(path: Path) -> dict[tuple[int, int], dict]:
    groups = defaultdict(list)
    for row in read_csv(path):
        if float(row["tau_start"]) >= -1 - 1e-9 and float(row["tau_end"]) <= 1 + 1e-9:
            groups[(int(row["seed"]), int(row["event_id"]))].append(row)
    out = {}
    for key, part in groups.items():
        def integral(field):
            return float(sum(float(r[field]) * float(r["duration"]) for r in part))
        out[key] = dict(
            actual_neck=float(sum(float(r["net_neck_count"]) for r in part)),
            supply=integral("plug_flux"), storage=integral("storage_contribution"),
            lateral=integral("lateral_contribution"), correction=integral("path_discretization_correction"),
        )
    return out


def extract_event_residuals(raw_dir: Path, cohort_path: Path, budget_path: Path,
                            torque_path: Path, geometry: dict) -> tuple[list[dict], dict]:
    budgets = aggregate_event_budgets(budget_path)
    torque = {(int(r["seed"]), int(r["event_id"])): r for r in read_csv(torque_path) if abs(float(r["lead"]) - 1) < 1e-9}
    events = read_csv(cohort_path); by = defaultdict(list)
    for row in events:
        by[int(row["seed"])].append(row)
    output = []; max_closure = 0.0
    for seed in sorted(by):
        reader = RawReader(raw_path(raw_dir, seed))
        try:
            for source in by[seed]:
                event_id = int(source["event_id"]); key = (seed, event_id)
                if key not in budgets or key not in torque:
                    continue
                step = int(source["first_native"]) - round(1 / DT_NATIVE)
                state = reader.frame("state", step); velocity = reader.frame("rdot", step + 1)
                source_sign = int(source["source_sign"]); direction = -source_sign
                p1 = ballistic_prediction(state, source_sign, 1.0, geometry)
                p2 = ballistic_prediction(state, source_sign, 2.0, geometry)
                p6 = ballistic_prediction(state, source_sign, 6.0, geometry)
                b = budgets[key]; reconstructed = b["supply"] + b["storage"] + b["lateral"] + b["correction"]
                max_closure = max(max_closure, abs(b["actual_neck"] - reconstructed))
                ap = approach_state(state, velocity, direction); t = torque[key]
                output.append(dict(
                    seed=seed, event_id=event_id, source_sign=source_sign,
                    success4=int(source["success4"]), abs_q_pre=float(source["abs_q_pre"]),
                    n_ballistic_to_entry=p1["n_ballistic"], n_ballistic_through_plus1=p2["n_ballistic"],
                    n_ballistic_through_plus5=p6["n_ballistic"], actual_neck_minus1_plus1=b["actual_neck"],
                    ballistic_residual=b["actual_neck"] - p2["n_ballistic"],
                    exact_supply=b["supply"], exact_storage=b["storage"], exact_lateral=b["lateral"], exact_correction=b["correction"],
                    nonreciprocal_alignment_turn=float(t["alignment_internal_alignment_visual_asymmetry"]),
                    total_social_alignment_turn=float(t["alignment_social_total"]), **ap,
                ))
        finally:
            reader.close()
    return output, dict(max_abs_integrated_budget_residual=max_closure, n_events=len(output))


def summarize_event_residuals(rows: list[dict]) -> list[dict]:
    output = []
    for cohort, part in (("train", [r for r in rows if r["seed"] < 30]),
                         ("heldout", [r for r in rows if r["seed"] >= 30])):
        for field in ("nonreciprocal_alignment_turn", "exact_supply", "exact_storage", "exact_lateral"):
            result = cluster_boot(part, lambda x, f=field: corr(x, "ballistic_residual", f), 20261170)
            output.append(dict(cohort=cohort, target="ballistic_residual", predictor=field, metric="Spearman", **result))
        for score in ("n_ballistic_to_entry", "n_ballistic_through_plus1", "Ja"):
            result = cluster_boot(part, lambda x, s=score: auc(x, s, "success4"), 20261171)
            output.append(dict(cohort=cohort, target="success4", predictor=score, metric="AUC", **result))
    return output


def intervention_test(raw_dir: Path, manifest_path: Path, branches_path: Path,
                      geometry: dict) -> tuple[list[dict], dict]:
    manifest = read_csv(manifest_path); predictors = {}
    by = defaultdict(list)
    for row in manifest:
        by[int(row["seed"])].append(row)
    for seed in sorted(by):
        reader = RawReader(raw_path(raw_dir, seed))
        try:
            for row in by[seed]:
                state = reader.frame("state", int(row["native_step"])); source_sign = int(row["source_sign"])
                p1 = ballistic_prediction(state, source_sign, 1.0, geometry)
                p5 = ballistic_prediction(state, source_sign, 5.0, geometry)
                predictors[row["sample_id"]] = dict(n_ballistic1=p1["n_ballistic"], n_ballistic5=p5["n_ballistic"])
        finally:
            reader.close()
    rows = []
    for row in read_csv(branches_path):
        rows.append(dict(sample_id=row["sample_id"], seed=int(row["seed"]), clone=int(row["clone"]), arm=row["arm"],
                         net_transfer=float(row["net_transfer"]), first_entry=int(row["first_entry"]), success4=int(row["success4"]),
                         turning_alignment_rate_1=float(row["turning_alignment_rate_1"]), **predictors[row["sample_id"]]))
    paired = defaultdict(dict)
    for row in rows:
        paired[(row["sample_id"], row["clone"])][row["arm"]] = row
    effects = []
    for pair in paired.values():
        a = pair["baseline"]; b = pair["symmetric_internal_alignment"]
        effects.append(dict(seed=a["seed"], sample_id=a["sample_id"], clone=a["clone"],
                            geometry_delta=b["n_ballistic5"] - a["n_ballistic5"],
                            net_transfer_delta=b["net_transfer"] - a["net_transfer"],
                            first_entry_delta=b["first_entry"] - a["first_entry"],
                            success4_delta=b["success4"] - a["success4"],
                            turning_delta=b["turning_alignment_rate_1"] - a["turning_alignment_rate_1"]))
    summary = {}
    for field in ("geometry_delta", "net_transfer_delta", "first_entry_delta", "success4_delta", "turning_delta"):
        summary[field] = cluster_boot(effects, lambda x, f=field: float(np.mean([r[f] for r in x])), 20261180)
    summary["max_abs_geometry_delta"] = max(abs(r["geometry_delta"]) for r in effects)
    summary["n_states"] = len(manifest); summary["n_paired_clones"] = len(effects)
    return rows, summary


def make_figure(attempts: list[dict], events: list[dict], intervention: list[dict],
                output: Path) -> None:
    held = [r for r in attempts if r["seed"] >= 30]
    ev = [r for r in events if r["seed"] >= 30]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].scatter([r["n_ballistic2"] for r in held], [r["net2"] for r in held],
                       c=[r["entry2"] for r in held], s=12, alpha=.55, cmap="coolwarm")
    lim = max(max(r["n_ballistic2"] for r in held), max(r["net2"] for r in held), 1)
    axes[0, 0].plot([0, lim], [0, lim], "k:")
    axes[0, 0].set(xlabel="ballistic fish by +2", ylabel="actual net fish at +2",
                   title="Held-out prospective attempts")
    axes[0, 1].scatter([r["n_ballistic5"] for r in held], [r["net5"] for r in held],
                       c=[r["sustained5"] for r in held], s=12, alpha=.55, cmap="coolwarm")
    axes[0, 1].set(xlabel="ballistic fish by +5", ylabel="actual net fish at +5",
                   title="Longer horizon: turning accumulates")
    axes[1, 0].scatter([r["nonreciprocal_alignment_turn"] for r in ev],
                       [r["ballistic_residual"] for r in ev], s=15, alpha=.6)
    axes[1, 0].axhline(0, color="k", ls=":")
    axes[1, 0].set(xlabel="lead-1 nonreciprocal alignment turn",
                   ylabel="actual minus ballistic fish (-1 to +1)",
                   title="Held-out event residual")
    by = defaultdict(dict)
    for row in intervention:
        by[(row["sample_id"], row["clone"])][row["arm"]] = row
    x = []; y = []
    for pair in by.values():
        x.append(pair["baseline"]["net_transfer"]); y.append(pair["symmetric_internal_alignment"]["net_transfer"])
    axes[1, 1].scatter(x, y, s=14, alpha=.55); lo = min(x + y); hi = max(x + y)
    axes[1, 1].plot([lo, hi], [lo, hi], "k:")
    axes[1, 1].set(xlabel="baseline transfer", ylabel="symmetrized transfer",
                   title="Matched noise; identical initial geometry")
    fig.tight_layout(); fig.savefig(output / "ballistic_encounter_transfer.png", dpi=220); plt.close(fig)


def self_test() -> None:
    geometry = dict(gap=2.0, segments=np.empty((0, 2, 2)), tip_centers=np.empty((0, 2)), tip_radius=.5)
    state = np.zeros((3, 4)); state[0] = [-2, -1, 1, 2]; state[1] = [0, .9, 0, 0]; state[2] = 0
    left = ballistic_prediction(state, -1, 2.1, geometry)
    reflected = state.copy(); reflected[0] *= -1; reflected[2] = np.pi
    right = ballistic_prediction(reflected, 1, 2.1, geometry)
    assert left["n_ballistic"] == right["n_ballistic"] == 2
    assert abs(left["first_ballistic_time"] - 1) < 1e-12
    wall = dict(gap=2.0, segments=np.asarray([[[ -.5, -.2], [-.5, .2]]]),
                tip_centers=np.empty((0, 2)), tip_radius=.5)
    assert ballistic_prediction(state, -1, 2.1, wall)["n_geometry_obstructed"] == 1
    print("PASS: finite-horizon rays, reflection, and physical obstruction")


def main() -> None:
    here = Path(__file__).resolve().parent; repo = here.parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data"))
    parser.add_argument("--source", type=Path, default=repo.parent / "BEM_bichamber_parallel_instance_no_hydro.py")
    parser.add_argument("--attempts", type=Path, default=here / "outputs/approach_renewal/attempt_episodes.csv")
    parser.add_argument("--feature-cache", type=Path, default=here / "outputs/initiation_propagation/cache/population_kinetics_features_dt0p5.npz")
    parser.add_argument("--cohort", type=Path, default=here / "outputs/collective_qab_pathway/event_outcomes.csv")
    parser.add_argument("--budgets", type=Path, default=here / "outputs/event_transport_budget/event_window_budgets.csv")
    parser.add_argument("--torques", type=Path, default=here / "outputs/social_torque_organization/pair_budget_samples.csv")
    parser.add_argument("--manifest", type=Path, default=here / "outputs/nonreciprocal_alignment_intervention/selection_manifest.csv")
    parser.add_argument("--branches", type=Path, default=here / "outputs/nonreciprocal_alignment_intervention/branches.csv")
    parser.add_argument("--output", type=Path, default=here / "outputs/ballistic_encounter_transfer")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    args.output.mkdir(parents=True, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        mod = load_gap20_production(args.source)
    geometry = production_geometry(mod)
    attempts, attempt_audit = extract_attempts(args.raw_dir, args.attempts, args.feature_cache, geometry)
    attempt_summary = summarize_attempts(attempts)
    events, budget_audit = extract_event_residuals(args.raw_dir, args.cohort, args.budgets, args.torques, geometry)
    event_summary = summarize_event_residuals(events)
    intervention, intervention_summary = intervention_test(args.raw_dir, args.manifest, args.branches, geometry)
    write_csv(args.output / "prospective_attempts.csv", attempts)
    write_csv(args.output / "heldout_validation.csv", attempt_summary)
    write_csv(args.output / "event_ballistic_budget.csv", events)
    write_csv(args.output / "event_residual_summary.csv", event_summary)
    write_csv(args.output / "intervention_geometry.csv", intervention)
    make_figure(attempts, events, intervention, args.output)
    summary = dict(
        configuration=dict(train_seeds="0-29", heldout_seeds="30-39", entry_horizon=ENTRY_HORIZON,
                           sustained_horizon=SUSTAINED_HORIZON, transfer_threshold=TRANSFER_THRESHOLD,
                           predictor="rigid translation of actual source positions by v0 times instantaneous mean-heading vector",
                           new_simulations=0, fitted_coefficients=0),
        audits=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                    gap=geometry["gap"], **attempt_audit, **budget_audit),
        attempt_validation=attempt_summary, event_residual_validation=event_summary,
        intervention=intervention_summary,
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["audits"], indent=2))


if __name__ == "__main__":
    main()
