#!/usr/bin/env python3
"""Separate individual-heading geometry from encounter-time reorientation.

At the established storage-excursion onset, compare coherent mean-heading
rays with independent fish-heading rays.  Derive the angular aperture and
minimum collective turn, audit complete-excursion reward, and use existing
source-exact torque samples for entries missed by both frozen models.
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

from analyze_approach_current_balance import load_gap20_production, raw_path
from analyze_ballistic_encounter_transfer import (
    cluster_boot, path_blocked, production_geometry, read_csv, source_mask, write_csv,
)
from analyze_initiation_propagation import RawReader, auc_score


DT = 0.001
HORIZON = 2.0
IN_NOISE = 0.5


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def global_required_turn(local_turn: float, source_sign: int) -> float:
    """Convert an inward-frame angular increment to the global frame."""
    return float(wrap(-source_sign * local_turn))


def mean_heading(theta: np.ndarray) -> tuple[float, float, np.ndarray]:
    vector = np.asarray([np.mean(np.cos(theta)), np.mean(np.sin(theta))])
    p = float(np.linalg.norm(vector)); angle = float(math.atan2(vector[1], vector[0]))
    return angle, p, vector


def angular_apertures(state: np.ndarray, source_sign: int, gap: float) -> dict:
    mask = source_mask(state, source_sign); indices = np.flatnonzero(mask)
    direction = -source_sign; distance = source_sign * state[0, indices]
    lower = np.arctan2(-gap / 2 - state[1, indices], distance)
    upper = np.arctan2(gap / 2 - state[1, indices], distance)
    theta_bar, p, vector = mean_heading(state[2, indices])
    phi_bar = float(math.atan2(vector[1], direction * vector[0]))
    # The inward-frame angle is periodic.  Outside the aperture, choose the
    # nearer wrapped endpoint rather than comparing unwrapped scalars.
    dlo = wrap(lower - phi_bar); dup = wrap(upper - phi_bar)
    nearest = np.where(np.abs(dlo) <= np.abs(dup), dlo, dup)
    signed = np.where((phi_bar >= lower) & (phi_bar <= upper), 0.0, nearest)
    pick = int(np.argmin(np.abs(signed)))
    weights = np.cos(state[2, indices] - theta_bar)
    variance_rate = IN_NOISE**2 * float(np.sum(weights**2)) / (len(indices)**2 * p**2)
    return dict(indices=indices, lower=lower, upper=upper,
                mean_angle=theta_bar, polarization=p, mean_vector=vector,
                minimum_turn_signed=float(signed[pick]),
                minimum_turn=abs(float(signed[pick])),
                median_aperture_width=float(np.median(upper - lower)),
                collective_angular_sigma=float(math.sqrt(variance_rate)),
                n_source=len(indices))


def count_rays(state: np.ndarray, source_sign: int, velocities: np.ndarray,
               geometry: dict, horizon: float = math.inf) -> tuple[int, int, list[float]]:
    indices = np.flatnonzero(source_mask(state, source_sign)); direction = -source_sign
    if velocities.ndim == 1:
        velocities = np.repeat(velocities[:, None], len(indices), axis=1)
    inward = direction * velocities[0]; distance = source_sign * state[0, indices]
    valid = inward > 0; times = np.full(len(indices), np.nan); times[valid] = distance[valid] / inward[valid]
    yhit = state[1, indices] + velocities[1] * times
    candidates = valid & (times > 0) & (times <= horizon + 1e-12) & (np.abs(yhit) <= geometry["gap"] / 2 + 1e-12)
    arrivals = []; blocked = 0
    for local in np.flatnonzero(candidates):
        a = state[:2, indices[local]]; b = np.asarray([0.0, yhit[local]])
        if path_blocked(a, b, geometry):
            blocked += 1
        else:
            arrivals.append(float(times[local]))
    return len(arrivals), blocked, arrivals


def geometry_metrics(state: np.ndarray, source_sign: int, geometry: dict) -> dict:
    aperture = angular_apertures(state, source_sign, geometry["gap"])
    indices = aperture["indices"]; common = aperture["mean_vector"]
    individual = np.vstack((np.cos(state[2, indices]), np.sin(state[2, indices])))
    coherent_any, cb_any, _ = count_rays(state, source_sign, common, geometry)
    individual_any, ib_any, _ = count_rays(state, source_sign, individual, geometry)
    coherent_2, cb_2, ct = count_rays(state, source_sign, common, geometry, HORIZON)
    individual_2, ib_2, it = count_rays(state, source_sign, individual, geometry, HORIZON)
    return dict(
        n_coherent_target=coherent_any, n_individual_target=individual_any,
        n_coherent_2=coherent_2, n_individual_2=individual_2,
        coherent_first_time=min(ct) if ct else float("nan"),
        individual_first_time=min(it) if it else float("nan"),
        coherent_blocked=cb_any + cb_2, individual_blocked=ib_any + ib_2,
        minimum_turn=aperture["minimum_turn"], minimum_turn_signed=aperture["minimum_turn_signed"],
        median_aperture_width=aperture["median_aperture_width"],
        collective_angular_sigma=aperture["collective_angular_sigma"],
        polarization=aperture["polarization"], n_source=aperture["n_source"],
        mean_heading_angle=aperture["mean_angle"],
    )


def extract_attempts(raw_dir: Path, attempts_path: Path, fixed5_path: Path,
                     geometry: dict) -> tuple[list[dict], dict]:
    original = [r for r in read_csv(attempts_path) if int(r["majority_source"]) == 1]
    fixed5 = {(int(r["seed"]), int(r["attempt_id"])): r for r in read_csv(fixed5_path)}
    by = defaultdict(list)
    for row in original:
        by[int(row["seed"])].append(row)
    output = []; source_mismatch = 0; blocked = 0
    for count, seed in enumerate(sorted(by), 1):
        print(f"individual/coherent geometry {count}/{len(by)} seed={seed}", flush=True)
        reader = RawReader(raw_path(raw_dir, seed)); sequence = sorted(by[seed], key=lambda r: int(r["attempt_id"]))
        try:
            for j, source in enumerate(sequence):
                key = (seed, int(source["attempt_id"])); step = int(round(float(source["onset_upper"]) / DT))
                state = reader.frame("state", step); direction = int(source["direction"]); source_sign = -direction
                metrics = geometry_metrics(state, source_sign, geometry)
                f5 = fixed5.get(key, {})
                q_onset = float(f5.get("q0", source["q_pre"]))
                n_right_onset = int(round(100 * (q_onset + 1) / 2))
                expected = n_right_onset if source_sign > 0 else 100 - n_right_onset
                source_mismatch += int(expected != metrics["n_source"])
                blocked += metrics["coherent_blocked"] + metrics["individual_blocked"]
                next_onset = float(sequence[j + 1]["onset_upper"]) if j + 1 < len(sequence) else float("nan")
                reset = float(source["reset_time"]); onset = float(source["onset_upper"])
                fixed_net = float(f5.get("net5", "nan"))
                multi = int(np.isfinite(next_onset) and reset < onset + 5 - 1e-12 and next_onset <= onset + 5 + 1e-12)
                entry_delay = float(source["first_entry_time_upper"]) - onset if int(source["first_entry"]) else float("nan")
                output.append(dict(
                    seed=seed, attempt_id=int(source["attempt_id"]), onset_time=onset, reset_time=reset,
                    direction=direction, source_sign=source_sign, q_pre=float(source["q_pre"]), q_onset=q_onset,
                    first_entry=int(source["first_entry"]), entry_delay=entry_delay,
                    excursion_reward=int(source["population_transfer"]), outcome=source["outcome"],
                    sustained_excursion=int(source["sustained"]), fixed5_net=fixed_net,
                    fixed5_sustained=int(np.isfinite(fixed_net) and fixed_net >= 5),
                    multiple_excursions_within5=multi, next_onset=next_onset, **metrics,
                ))
        finally:
            reader.close()
    return output, dict(n_attempts=len(output), n_source_count_mismatch=source_mismatch,
                        n_native_wall_obstructions=blocked)


def auc(rows: list[dict], score: str, target: str) -> float:
    return float(auc_score(np.asarray([int(r[target]) for r in rows]),
                           np.asarray([float(r[score]) for r in rows])))


def corr(rows: list[dict], x: str, y: str) -> float:
    a = np.asarray([float(r[x]) for r in rows]); b = np.asarray([float(r[y]) for r in rows])
    return float(spearmanr(a, b).statistic) if np.std(a) and np.std(b) else float("nan")


def binary_metrics(rows: list[dict], score: str, target: str) -> dict:
    p = np.asarray([float(r[score]) > 0 for r in rows]); y = np.asarray([bool(int(r[target])) for r in rows])
    sensitivity = float(np.mean(p[y])); specificity = float(np.mean(~p[~y]))
    return dict(sensitivity=sensitivity, specificity=specificity,
                balanced_accuracy=(sensitivity + specificity) / 2,
                predicted_fraction=float(np.mean(p)), observed_fraction=float(np.mean(y)))


def summarize_attempts(rows: list[dict]) -> tuple[list[dict], dict]:
    out = []
    for cohort, part in (("train", [r for r in rows if r["seed"] < 30]),
                         ("heldout", [r for r in rows if r["seed"] >= 30])):
        for score in ("n_coherent_target", "n_individual_target", "n_coherent_2", "n_individual_2"):
            result = cluster_boot(part, lambda x, s=score: auc(x, s, "first_entry"), 20261190)
            out.append(dict(cohort=cohort, target="first_entry", predictor=score, metric="AUC", **result))
            binary = binary_metrics(part, score, "first_entry")
            bci = cluster_boot(part, lambda x, s=score: binary_metrics(x, s, "first_entry")["balanced_accuracy"], 20261191)
            out.append(dict(cohort=cohort, target="first_entry", predictor=score, metric="binary_targeting",
                            estimate=binary["balanced_accuracy"], ci95=bci["ci95"], n=len(part), n_seeds=len({r["seed"] for r in part}), **binary))
        for score in ("n_coherent_target", "n_individual_target"):
            result = cluster_boot(part, lambda x, s=score: corr(x, s, "excursion_reward"), 20261192)
            out.append(dict(cohort=cohort, target="excursion_reward", predictor=score, metric="Spearman", **result))
            sustained = cluster_boot(part, lambda x, s=score: auc(x, s, "sustained_excursion"), 20261194)
            out.append(dict(cohort=cohort, target="sustained_excursion", predictor=score, metric="AUC", **sustained))
        success = [r for r in part if r["first_entry"]]
        for field in ("minimum_turn", "median_aperture_width", "collective_angular_sigma"):
            result = cluster_boot(success, lambda x, f=field: float(np.median([r[f] for r in x])), 20261193)
            out.append(dict(cohort=cohort, target="successful_entry", predictor=field, metric="median", **result))
    contamination = {}
    for cohort, part in (("all", rows), ("heldout", [r for r in rows if r["seed"] >= 30])):
        valid = [r for r in part if np.isfinite(r["fixed5_net"])]
        contamination[cohort] = dict(
            n=len(valid), multiple_excursions_fraction=float(np.mean([r["multiple_excursions_within5"] for r in valid])),
            fixed5_reward_mismatch_fraction=float(np.mean([r["fixed5_net"] != r["excursion_reward"] for r in valid])),
            sustained_label_disagreement=float(np.mean([r["fixed5_sustained"] != r["sustained_excursion"] for r in valid])),
            mean_fixed5_minus_excursion=float(np.mean([r["fixed5_net"] - r["excursion_reward"] for r in valid])),
        )
    return out, contamination


def associate_events(attempts: list[dict], cohort_path: Path) -> list[dict]:
    by = defaultdict(list)
    for row in attempts:
        by[int(row["seed"])].append(row)
    associated = []
    for source in read_csv(cohort_path):
        seed = int(source["seed"]); time = int(source["first_native"]) * DT; direction = -int(source["source_sign"])
        matches = [r for r in by[seed] if r["direction"] == direction and r["onset_time"] <= time <= r["reset_time"] + 1e-9]
        if len(matches) == 1:
            associated.append(dict(matches[0], event_id=int(source["event_id"]), first_native=int(source["first_native"]),
                                   success4=int(source["success4"]), event_time=time))
    return associated


def missed_entry_evolution(raw_dir: Path, associated: list[dict], kinematics_path: Path,
                           pair_path: Path) -> tuple[list[dict], list[dict], dict]:
    kin = {(int(r["seed"]), int(r["event_id"]), float(r["lead"])): r for r in read_csv(kinematics_path)}
    pair = {(int(r["seed"]), int(r["event_id"]), float(r["lead"])): r for r in read_csv(pair_path)}
    missed = [r for r in associated if r["first_entry"] and r["n_coherent_target"] == 0 and r["n_individual_target"] == 0]
    by = defaultdict(list)
    for row in missed:
        by[int(row["seed"])].append(row)
    evolution = []; torques = []; membership_mismatch = 0
    for seed in sorted(by):
        reader = RawReader(raw_path(raw_dir, seed))
        try:
            for source in by[seed]:
                start = int(round(source["onset_time"] / DT)); before = source["first_native"] - 1
                s0 = reader.frame("state", start); s1 = reader.frame("state", before)
                mask = source_mask(s0, source["source_sign"]); indices = np.flatnonzero(mask)
                _, _, v0 = mean_heading(s0[2, indices]); theta0 = math.atan2(v0[1], v0[0])
                _, _, v1 = mean_heading(s1[2, indices]); theta1 = math.atan2(v1[1], v1[0])
                membership_mismatch += int(np.any(source_mask(s1, source["source_sign"])[indices] == 0))
                change = float(wrap(theta1 - theta0))
                required_local = float(source["minimum_turn_signed"])
                required = global_required_turn(required_local, int(source["source_sign"]))
                duration = source["event_time"] - source["onset_time"]
                noise_rms = source["collective_angular_sigma"] * math.sqrt(max(duration, 0))
                sign = float(np.sign(required))
                evolution.append(dict(
                    seed=seed, event_id=source["event_id"], attempt_id=source["attempt_id"],
                    onset_time=source["onset_time"], event_time=source["event_time"], duration=duration,
                    required_turn_local=required_local, required_turn_signed=required, minimum_turn=abs(required), actual_heading_change=change,
                    actual_toward_required=sign * change, collective_noise_rms=noise_rms,
                    turn_to_noise_ratio=abs(required) / noise_rms if noise_rms > 0 else float("inf"),
                    reaches_required_direction=int(sign * change > 0),
                    reaches_required_magnitude=int(sign * change >= abs(required)),
                    seed_split="heldout" if seed >= 30 else "train",
                ))
                for lead in (1.0, .5, .2):
                    key = (seed, source["event_id"], lead)
                    if key not in kin or lead > duration + .01:
                        continue
                    k = kin[key]; p = pair.get(key, {})
                    torques.append(dict(
                        seed=seed, event_id=source["event_id"], attempt_id=source["attempt_id"], lead=lead,
                        required_turn_sign=sign,
                        wall_toward_required=sign * float(k["omega_wall"]),
                        social_toward_required=sign * float(k["omega_social"]),
                        deterministic_toward_required=sign * float(k["omega_deterministic"]),
                        nonreciprocal_gap_alignment=float(p.get("alignment_internal_alignment_visual_asymmetry", "nan")),
                        seed_split="heldout" if seed >= 30 else "train",
                    ))
        finally:
            reader.close()
    return evolution, torques, dict(n_associated=len(associated), n_missed_both=len(missed),
                                    n_membership_mismatch=membership_mismatch)


def summarize_missed(evolution: list[dict], torques: list[dict]) -> list[dict]:
    out = []
    for cohort in ("train", "heldout"):
        part = [r for r in evolution if r["seed_split"] == cohort]
        for field in ("minimum_turn", "actual_toward_required", "collective_noise_rms", "turn_to_noise_ratio"):
            result = cluster_boot(part, lambda x, f=field: float(np.median([r[f] for r in x])), 20261200)
            out.append(dict(cohort=cohort, lead="onset_to_entry", term=field, metric="median", **result))
        for field in ("reaches_required_direction", "reaches_required_magnitude"):
            result = cluster_boot(part, lambda x, f=field: float(np.mean([r[f] for r in x])), 20261201)
            out.append(dict(cohort=cohort, lead="onset_to_entry", term=field, metric="fraction", **result))
        for lead in (1.0, .5, .2):
            selected = [r for r in torques if r["seed_split"] == cohort and r["lead"] == lead]
            for field in ("wall_toward_required", "social_toward_required", "deterministic_toward_required", "nonreciprocal_gap_alignment"):
                result = cluster_boot(selected, lambda x, f=field: float(np.mean([r[f] for r in x])), 20261202 + round(10 * lead))
                out.append(dict(cohort=cohort, lead=lead, term=field, metric="mean_rate", **result))
    return out


def make_figure(attempts: list[dict], evolution: list[dict], torques: list[dict], output: Path) -> None:
    held = [r for r in attempts if r["seed"] >= 30]; evo = [r for r in evolution if r["seed_split"] == "heldout"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].hist([r["n_coherent_target"] for r in held], bins=np.arange(-.5, 15.5), alpha=.65, label="mean heading")
    axes[0, 0].hist([r["n_individual_target"] for r in held], bins=np.arange(-.5, 15.5), alpha=.55, label="individual headings")
    axes[0, 0].set(xlabel="fish geometrically targeting opening", ylabel="attempts", title="Held-out onset geometry"); axes[0, 0].legend(frameon=False)
    axes[0, 1].scatter([r["n_individual_target"] for r in held], [r["excursion_reward"] for r in held],
                       c=[r["first_entry"] for r in held], cmap="coolwarm", s=12, alpha=.55)
    axes[0, 1].set(xlabel="individual frozen rays", ylabel="excursion-end net fish", title="Complete-excursion reward")
    axes[1, 0].scatter([r["minimum_turn"] for r in evo], [r["actual_toward_required"] for r in evo], s=16, alpha=.65)
    limit = max([r["minimum_turn"] for r in evo] + [r["actual_toward_required"] for r in evo]); axes[1, 0].plot([0, limit], [0, limit], "k:")
    axes[1, 0].set(xlabel="minimum geometric turn at onset", ylabel="actual heading turn toward aperture", title="Held-out entries missed by both rays")
    data = []; labels = []
    for lead in (1.0, .5, .2):
        part = [r for r in torques if r["seed_split"] == "heldout" and r["lead"] == lead]
        data.extend([[r["wall_toward_required"] for r in part], [r["social_toward_required"] for r in part]])
        labels.extend([f"wall\n-{lead:g}", f"social\n-{lead:g}"])
    axes[1, 1].boxplot(data, tick_labels=labels, showfliers=False); axes[1, 1].axhline(0, color="k", ls=":")
    axes[1, 1].set(ylabel="collective angular rate toward required turn", title="Source-exact torque projections")
    fig.tight_layout(); fig.savefig(output / "individual_vs_coherent_geometry.png", dpi=220); plt.close(fig)


def self_test() -> None:
    state = np.zeros((3, 3)); state[0] = [-2, -1, 2]; state[1] = [0, .5, 0]; state[2] = [0, .4, 0]
    a = angular_apertures(state, -1, 2.0)
    reflected = state.copy(); reflected[0] *= -1; reflected[2] = wrap(np.pi - reflected[2])
    b = angular_apertures(reflected, 1, 2.0)
    assert abs(a["minimum_turn"] - b["minimum_turn"]) < 1e-12
    ga = global_required_turn(a["minimum_turn_signed"], -1)
    gb = global_required_turn(b["minimum_turn_signed"], 1)
    assert abs(float(wrap(ga + gb))) < 1e-12
    # A global heading increment and every global torque reverse under the
    # same reflection, so their projection toward the requirement is invariant.
    test_change = .17
    assert abs(np.sign(ga) * test_change - np.sign(gb) * (-test_change)) < 1e-12
    # A heading just across the +/-pi seam must choose the wrapped-nearest
    # aperture endpoint, not the endpoint nearest on the unwrapped line.
    seam = state.copy(); seam[2, :2] = math.pi - .01
    c = angular_apertures(seam, -1, 2.0)
    lo, hi = c["lower"], c["upper"]; phi = math.pi - .01
    brute = min((float(wrap(x - phi)) for x in np.r_[lo, hi]), key=abs)
    assert abs(float(wrap(c["minimum_turn_signed"] - brute))) < 1e-12
    assert a["collective_angular_sigma"] > 0
    print("PASS: reflected global turn, periodic aperture, and collective diffusion scale")


def main() -> None:
    here = Path(__file__).resolve().parent; repo = here.parents[1]
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, default=Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data"))
    p.add_argument("--source", type=Path, default=repo.parent / "BEM_bichamber_parallel_instance_no_hydro.py")
    p.add_argument("--attempts", type=Path, default=here / "outputs/approach_renewal/attempt_episodes.csv")
    p.add_argument("--fixed5", type=Path, default=here / "outputs/ballistic_encounter_transfer/prospective_attempts.csv")
    p.add_argument("--cohort", type=Path, default=here / "outputs/collective_qab_pathway/event_outcomes.csv")
    p.add_argument("--kinematics", type=Path, default=here / "outputs/geometric_first_passage_reduction/kinematic_samples.csv")
    p.add_argument("--pairs", type=Path, default=here / "outputs/social_torque_organization/pair_budget_samples.csv")
    p.add_argument("--output", type=Path, default=here / "outputs/individual_heading_geometry")
    p.add_argument("--self-test", action="store_true"); args = p.parse_args()
    if args.self_test:
        self_test(); return
    args.output.mkdir(parents=True, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        mod = load_gap20_production(args.source)
    geometry = production_geometry(mod)
    attempts, attempt_audit = extract_attempts(args.raw_dir, args.attempts, args.fixed5, geometry)
    validation, contamination = summarize_attempts(attempts)
    associated = associate_events(attempts, args.cohort)
    evolution, torques, missed_audit = missed_entry_evolution(args.raw_dir, associated, args.kinematics, args.pairs)
    missed_summary = summarize_missed(evolution, torques)
    write_csv(args.output / "prospective_geometry.csv", attempts)
    write_csv(args.output / "heldout_validation.csv", validation)
    write_csv(args.output / "missed_entry_evolution.csv", evolution)
    write_csv(args.output / "missed_entry_torques.csv", torques)
    write_csv(args.output / "missed_entry_summary.csv", missed_summary)
    make_figure(attempts, evolution, torques, args.output)
    summary = dict(
        configuration=dict(train_seeds="0-29", heldout_seeds="30-39", fixed_ray_horizon=HORIZON,
                           primary_reward="population at physical excursion reset minus population before onset",
                           angular_noise_amplitude=IN_NOISE, new_simulations=0, fitted_parameters=0),
        audits=dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(), gap=geometry["gap"],
                    **attempt_audit, **missed_audit),
        fixed5_contamination=contamination, validation=validation, missed_entry=missed_summary,
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["audits"], indent=2))


if __name__ == "__main__":
    main()
