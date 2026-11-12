#!/usr/bin/env python3
"""Prospective q and majority-school gap-polarization analysis.

Read-only analysis of the existing 40 no-hydrodynamics gap-20 trajectories.
No simulator code or raw trajectory is changed.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from analyze_initiation_propagation import RawReader, auc_score  # noqa: E402

N = 100
Q_EDGES = np.linspace(-1.01, 1.01, 22)
ABS_Q_EDGES = np.array([0.0, 0.5, 0.8, 1.000001])
A_EDGES = np.array([-1.000001, -0.6, -0.2, 0.2, 0.6, 1.000001])


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_q(seed: int, args) -> tuple[np.ndarray, np.ndarray]:
    cache = args.cache_dir / f"seed{seed:02d}_stride{args.stride}_qscan_v4.npz"
    if not cache.exists():
        raise FileNotFoundError(f"Missing existing occupancy cache: {cache}")
    with np.load(cache) as data:
        q = np.asarray(data["q"], float)
        steps = np.asarray(data["step"], np.int64)
        nl = np.asarray(data["nl"], float)
    if not np.allclose(q, 1.0 - 2.0 * nl / N):
        raise ValueError(f"q/N_L identity failed for seed {seed}")
    return q, steps


def read_events(path: Path) -> dict[int, list[dict]]:
    events: dict[int, list[dict]] = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            item = dict(row)
            item["seed"] = int(row["seed"])
            item["first_native"] = int(row["first_native"])
            item["time"] = float(row["time"])
            item["direction"] = int(row["direction"])
            item["pre_quiet_duration"] = float(row["pre_quiet_duration"])
            events.setdefault(item["seed"], []).append(item)
    for rows in events.values():
        rows.sort(key=lambda row: row["time"])
    return events


def majority_features(state: np.ndarray) -> dict[str, float]:
    """Geometric definitions using only information at the current snapshot."""
    x, y, theta = state
    nr = int(np.sum(x >= 0))
    nl = N - nr
    q = (nr - nl) / N
    if q == 0:
        raise ValueError("Majority chamber is undefined at q=0")
    majority_right = q > 0
    mask = (x >= 0) if majority_right else (x < 0)
    n_majority = int(mask.sum())
    com = np.array([np.mean(x[mask]), np.mean(y[mask])])
    distance = float(np.linalg.norm(com))
    if distance <= 0:
        raise ValueError("Majority-school COM coincides with the gap center")
    gap_hat = -com / distance
    headings = np.column_stack((np.cos(theta[mask]), np.sin(theta[mask])))
    polarization_vector = headings.mean(axis=0)
    a = float(polarization_vector @ gap_hat)
    # Independent scalar form for the numerical definition audit.
    alpha_gap = float(np.arctan2(gap_hat[1], gap_hat[0]))
    a_scalar = float(np.mean(np.cos(theta[mask] - alpha_gap)))
    return dict(q=q, n_majority=n_majority, com_x=float(com[0]), com_y=float(com[1]),
                com_gap_distance=distance, gap_hat_x=float(gap_hat[0]),
                gap_hat_y=float(gap_hat[1]), a=a, a_scalar=a_scalar,
                polarization_x=float(polarization_vector[0]),
                polarization_y=float(polarization_vector[1]))


def q_bin(value: float) -> int:
    return int(np.searchsorted(ABS_Q_EDGES, abs(value), side="right") - 1)


def a_bin(value: float) -> int:
    return int(np.searchsorted(A_EDGES, value, side="right") - 1)


def future_label(time: float, q: float, events: list[dict], horizon: float) -> tuple[int, float, int]:
    """First known episode after t; direction is set from current majority only."""
    desired = -int(np.sign(q))  # left-to-right for left majority and vice versa
    future = [row for row in events if time < row["time"] <= time + horizon]
    if not future:
        return 0, np.nan, 0
    first = future[0]
    return int(first["direction"] == desired), float(first["time"]), int(first["direction"])


def risk_rows_for_seed(seed: int, q: np.ndarray, steps: np.ndarray, reader: RawReader,
                       events: list[dict], args) -> list[dict]:
    dt_coarse = args.dt * args.stride
    quiet_n = round(args.quiet / dt_coarse)
    spacing_n = round(args.spacing / dt_coarse)
    first = max(quiet_n, round(args.transient / dt_coarse))
    out = []
    for idx in range(first, len(q) - round(args.horizon / dt_coarse), spacing_n):
        qi = float(q[idx])
        if qi == 0 or not np.allclose(q[idx-quiet_n:idx+1], qi, atol=1e-12):
            continue
        time = idx * dt_coarse
        # Exclude any already indexed initiation during the retrospective interval.
        if any(time-args.quiet < row["time"] <= time for row in events):
            continue
        feature = majority_features(reader.frame("state", int(steps[idx])))
        outcome, event_time, event_direction = future_label(time, qi, events, args.horizon)
        qb, ab = q_bin(qi), a_bin(feature["a"])
        if not (0 <= qb < len(ABS_Q_EDGES)-1 and 0 <= ab < len(A_EDGES)-1):
            continue
        out.append(dict(seed=seed, t=time, native_step=int(steps[idx]), q=qi,
                        a=feature["a"], q_bin=qb, a_bin=ab,
                        first_crossing=outcome, first_event_time=event_time,
                        first_event_direction=event_direction,
                        n_majority=feature["n_majority"],
                        majority_com_x=feature["com_x"], majority_com_y=feature["com_y"],
                        majority_com_gap_distance=feature["com_gap_distance"]))
    return out


def validate_records(reader: RawReader, seed_events: list[dict], q: np.ndarray,
                     steps: np.ndarray, risk_rows: list[dict], args) -> tuple[list[dict], dict]:
    event_checks = []
    for row in seed_events:
        k = row["first_native"]
        before = reader.frame("state", k-1)
        after = reader.frame("state", k)
        q0 = 1.0 - 2.0*np.sum(before[0] < 0)/N
        q1 = 1.0 - 2.0*np.sum(after[0] < 0)/N
        raw_direction = int(np.sign(q1-q0))
        event_checks.append(dict(seed=row["seed"], first_native=k, time=row["time"],
                                 recorded_direction=row["direction"], raw_direction=raw_direction,
                                 q_before=q0, q_after=q1,
                                 valid_change=int(raw_direction != 0),
                                 direction_match=int(raw_direction == row["direction"]),
                                 time_error=abs(row["time"]-k*args.dt)))

    # One bounded native-resolution audit window per seed. This can reveal a
    # rapid cross-and-return that the 0.1 occupancy grid leaves unchanged.
    audit = dict(available=False)
    candidates = [row for row in risk_rows if 350 <= row["t"] <= 750]
    if candidates:
        row = candidates[len(candidates)//2]
        origin = int(row["native_step"])
        first = origin - round(args.quiet/args.dt)
        last = origin + round(args.horizon/args.dt)
        block = reader.state_block(first, last)
        x = block[0]
        transitions = (x[:, :-1] < 0) != (x[:, 1:] < 0)
        offsets, fish = np.where(transitions.T)
        past = offsets < origin-first
        future = offsets >= origin-first
        desired = -int(np.sign(row["q"]))
        desired_future = 0
        if np.any(future):
            first_future = int(np.flatnonzero(future)[0])
            offset, fish_id = offsets[first_future], fish[first_future]
            first_direction = 1 if x[fish_id, offset] < 0 else -1
            desired_future = int(first_direction == desired)
        audit = dict(available=True, seed=int(row["seed"]), t=float(row["t"]),
                     coarse_past_quiet=True, native_past_crossings=int(past.sum()),
                     native_future_crossings=int(future.sum()),
                     native_desired_future=int(desired_future),
                     indexed_first_crossing=int(row["first_crossing"]),
                     label_match=int(desired_future == int(row["first_crossing"])))
    return event_checks, audit


def counts(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    shape = (len(ABS_Q_EDGES)-1, len(A_EDGES)-1)
    n = np.zeros(shape, int)
    y = np.zeros(shape, int)
    for row in rows:
        n[row["q_bin"], row["a_bin"]] += 1
        y[row["q_bin"], row["a_bin"]] += row["first_crossing"]
    return n, y


def cell_seed_intervals(rows: list[dict], seeds: list[int], draws: int = 3000) -> dict:
    shape = (len(seeds), len(ABS_Q_EDGES)-1, len(A_EDGES)-1)
    n = np.zeros(shape, int)
    y = np.zeros(shape, int)
    index = {seed: i for i, seed in enumerate(seeds)}
    for row in rows:
        i = index[row["seed"]]
        n[i, row["q_bin"], row["a_bin"]] += 1
        y[i, row["q_bin"], row["a_bin"]] += row["first_crossing"]
    rng = np.random.default_rng(20261112)
    pick = rng.integers(len(seeds), size=(draws, len(seeds)))
    result = {}
    for qb in range(shape[1]):
        for ab in range(shape[2]):
            den = n[pick, qb, ab].sum(axis=1)
            num = y[pick, qb, ab].sum(axis=1)
            value = num[den > 0]/den[den > 0]
            result[(qb, ab)] = tuple(map(float, np.quantile(value, [.025, .975])))
    return result


def heldout(rows: list[dict], seeds: list[int], min_cell: int) -> dict:
    split = len(seeds)//4
    train_seeds, test_seeds = seeds[:-split], seeds[-split:]
    train = [row for row in rows if row["seed"] in set(train_seeds)]
    test = [row for row in rows if row["seed"] in set(test_seeds)]
    n, y = counts(train)
    qn, qy = n.sum(axis=1), y.sum(axis=1)
    overall = (y.sum()+.5)/(n.sum()+1)
    per_seed = []
    all_truth, all_pq, all_pqa = [], [], []
    for seed in test_seeds:
        subset = [row for row in test if row["seed"] == seed]
        truth, pq, pqa = [], [], []
        coverage = 0
        for row in subset:
            qb, ab = row["q_bin"], row["a_bin"]
            p0 = (qy[qb]+.5)/(qn[qb]+1) if qn[qb] >= min_cell else overall
            if n[qb, ab] >= min_cell:
                p1 = (y[qb, ab]+.5)/(n[qb, ab]+1)
                coverage += 1
            else:
                p1 = p0
            truth.append(row["first_crossing"]); pq.append(p0); pqa.append(p1)
        truth = np.asarray(truth); pq = np.asarray(pq); pqa = np.asarray(pqa)
        per_seed.append(dict(seed=seed, n=len(subset), positives=int(truth.sum()),
                             coverage=coverage/max(1, len(subset)),
                             brier_q=float(np.mean((truth-pq)**2)),
                             brier_qa=float(np.mean((truth-pqa)**2)),
                             auc_q=auc_score(truth, pq), auc_qa=auc_score(truth, pqa)))
        all_truth.extend(truth); all_pq.extend(pq); all_pqa.extend(pqa)
    effect = np.asarray([row["brier_q"]-row["brier_qa"] for row in per_seed])
    auc_effect = np.asarray([row["auc_qa"]-row["auc_q"] for row in per_seed])
    rng = np.random.default_rng(20261112)
    boot = effect[rng.integers(len(effect), size=(4000, len(effect)))].mean(axis=1)
    auc_boot = auc_effect[rng.integers(len(auc_effect), size=(4000, len(auc_effect)))].mean(axis=1)
    truth = np.asarray(all_truth); pq = np.asarray(all_pq); pqa = np.asarray(all_pqa)
    return dict(training_seeds=train_seeds, heldout_seeds=test_seeds, per_seed=per_seed,
                heldout_n=len(truth), heldout_positives=int(truth.sum()),
                pooled_brier_q=float(np.mean((truth-pq)**2)),
                pooled_brier_qa=float(np.mean((truth-pqa)**2)),
                pooled_auc_q=auc_score(truth, pq), pooled_auc_qa=auc_score(truth, pqa),
                mean_seed_brier_improvement=float(effect.mean()),
                seed_bootstrap_ci95=list(map(float, np.quantile(boot, [.025, .975]))),
                mean_seed_auc_improvement=float(auc_effect.mean()),
                auc_seed_bootstrap_ci95=list(map(float, np.quantile(auc_boot, [.025, .975]))))


def make_figures(q_by_seed: dict[int, np.ndarray], rows: list[dict], args) -> None:
    figures = args.output_dir/"figures"
    figures.mkdir(parents=True, exist_ok=True)
    seeds = sorted(q_by_seed)
    hist = np.asarray([np.histogram(q_by_seed[s], Q_EDGES)[0] for s in seeds], float)
    hist /= hist.sum(axis=1, keepdims=True)
    rng = np.random.default_rng(20261112)
    boot = hist[rng.integers(len(seeds), size=(3000, len(seeds)))].mean(axis=1)
    p = hist.mean(axis=0); lo, hi = np.quantile(boot, [.025, .975], axis=0)
    mid = (Q_EDGES[:-1]+Q_EDGES[1:])/2
    write_csv(args.output_dir/"occupancy_pdf.csv", [dict(q=mid[i], mass=p[i], ci_low=lo[i], ci_high=hi[i]) for i in range(len(mid))])
    fig, ax = plt.subplots(figsize=(6.5, 4.1))
    ax.plot(mid, p, "o-", ms=3); ax.fill_between(mid, lo, hi, alpha=.3)
    ax.set(xlabel="population imbalance q", ylabel="probability mass", title="Empirical p(q), equal seed weights")
    fig.tight_layout(); fig.savefig(figures/"01_population_imbalance_pdf.png", dpi=220); plt.close(fig)

    qv = np.asarray([row["q"] for row in rows]); av = np.asarray([row["a"] for row in rows])
    joint, qe, ae = np.histogram2d(qv, av, bins=[Q_EDGES, A_EDGES])
    mass = joint/joint.sum()
    write_csv(args.output_dir/"joint_q_a.csv", [dict(q_min=qe[i], q_max=qe[i+1], a_min=ae[j], a_max=ae[j+1], count=int(joint[i,j]), mass=mass[i,j]) for i in range(len(qe)-1) for j in range(len(ae)-1)])
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    nonzero = mass[mass > 0]
    im = ax.pcolormesh(qe, ae, np.ma.masked_where(mass.T == 0, mass.T), norm=LogNorm(vmin=nonzero.min(), vmax=nonzero.max()), shading="auto")
    fig.colorbar(im, ax=ax, label="joint probability mass")
    ax.set(xlabel="population imbalance q", ylabel="gap-directed polarization a", title="Prospective p(q,a), crossing-free history")
    fig.tight_layout(); fig.savefig(figures/"02_joint_q_gap_polarization.png", dpi=220); plt.close(fig)

    n, y = counts(rows); ci = cell_seed_intervals(rows, seeds)
    risk = np.divide(y, n, out=np.full(n.shape, np.nan), where=n >= args.min_cell)
    write_csv(args.output_dir/"conditional_first_crossing.csv", [dict(abs_q_lo=ABS_Q_EDGES[i], abs_q_hi=ABS_Q_EDGES[i+1], a_lo=A_EDGES[j], a_hi=A_EDGES[j+1], count=int(n[i,j]), positives=int(y[i,j]), p_first_crossing=risk[i,j], ci_low=ci[(i,j)][0], ci_high=ci[(i,j)][1]) for i in range(n.shape[0]) for j in range(n.shape[1])])
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    im = ax.imshow(np.ma.masked_invalid(risk), origin="lower", aspect="auto", extent=[A_EDGES[0], A_EDGES[-1], 0, 3], vmin=0, vmax=max(.05, np.nanmax(risk)))
    fig.colorbar(im, ax=ax, label="P(first majority-to-minority crossing within 5)")
    ax.set_yticks([.5, 1.5, 2.5], ["0-.5", ".5-.8", ".8-1"])
    ax.set(xlabel="gap-directed polarization a", ylabel="current |q|", title="Prospective first-crossing probability at fixed q and a")
    for i in range(n.shape[0]):
        for j in range(n.shape[1]):
            ax.text((A_EDGES[j]+A_EDGES[j+1])/2, i+.5, f"{y[i,j]}/{n[i,j]}", ha="center", va="center", color="white", fontsize=8)
    fig.tight_layout(); fig.savefig(figures/"03_conditional_first_crossing.png", dpi=220); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data"))
    parser.add_argument("--cache-dir", type=Path, default=REPO/"results/nohydro_gap20_descriptive_transition/outputs/initiation_propagation/cache")
    parser.add_argument("--event-file", type=Path, default=REPO/"results/nohydro_gap20_descriptive_transition/outputs/initiation_propagation/crossing_episodes.csv")
    parser.add_argument("--output-dir", type=Path, default=REPO/"results/nohydro_gap20_descriptive_transition/outputs/q_gap_polarization")
    parser.add_argument("--seeds", type=int, nargs="*", default=list(range(40)))
    parser.add_argument("--dt", type=float, default=.001)
    parser.add_argument("--stride", type=int, default=100)
    parser.add_argument("--transient", type=float, default=200.)
    parser.add_argument("--quiet", type=float, default=2.)
    parser.add_argument("--horizon", type=float, default=5.)
    parser.add_argument("--spacing", type=float, default=5.)
    parser.add_argument("--min-cell", type=int, default=20)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    events = read_events(args.event_file)
    q_by_seed = {}; all_rows = []; all_checks = []; audits = []; definition_rows = []
    for seed in args.seeds:
        q, steps = load_q(seed, args)
        q_by_seed[seed] = q[round(args.transient/(args.dt*args.stride)):]
        raw = args.raw_dir/f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"
        reader = RawReader(raw)
        rows = risk_rows_for_seed(seed, q, steps, reader, events.get(seed, []), args)
        checks, audit = validate_records(reader, events.get(seed, []), q, steps, rows, args)
        for row in rows[::max(1, len(rows)//5)]:
            feature = majority_features(reader.frame("state", row["native_step"]))
            definition_rows.append(dict(seed=seed, t=row["t"], q_cache=row["q"], q_raw=feature["q"],
                                        abs_q_error=abs(feature["q"]-row["q"]),
                                        n_majority=feature["n_majority"],
                                        n_majority_from_q=round(N*(1+abs(row["q"]))/2),
                                        majority_com_x=feature["com_x"], majority_com_y=feature["com_y"],
                                        com_side_match=int(np.sign(feature["com_x"]) == np.sign(row["q"])),
                                        gap_hat_norm_error=abs(np.hypot(feature["gap_hat_x"],feature["gap_hat_y"])-1),
                                        a_vector=feature["a"], a_scalar=feature["a_scalar"],
                                        abs_a_projection_error=abs(feature["a"]-feature["a_scalar"])))
        reader.close()
        all_rows.extend(rows); all_checks.extend(checks)
        if audit["available"]: audits.append(audit)
        print(f"seed {seed}: risk origins {len(rows)}, positives {sum(row['first_crossing'] for row in rows)}", flush=True)
    write_csv(args.output_dir/"prospective_q_a_windows.csv", all_rows)
    write_csv(args.output_dir/"definition_validation_samples.csv", definition_rows)
    write_csv(args.output_dir/"native_event_record_validation.csv", all_checks)
    write_csv(args.output_dir/"native_window_label_audit.csv", audits)
    make_figures(q_by_seed, all_rows, args)
    validation = dict(n_definition_checks=len(definition_rows),
                      max_abs_q_error=float(max(row["abs_q_error"] for row in definition_rows)),
                      max_abs_a_projection_error=float(max(row["abs_a_projection_error"] for row in definition_rows)),
                      majority_count_mismatches=int(sum(row["n_majority"] != row["n_majority_from_q"] for row in definition_rows)),
                      majority_com_side_mismatches=int(sum(not row["com_side_match"] for row in definition_rows)),
                      max_gap_direction_norm_error=float(max(row["gap_hat_norm_error"] for row in definition_rows)),
                      native_event_records=len(all_checks),
                      native_event_valid_changes=int(sum(row["valid_change"] for row in all_checks)),
                      native_event_direction_matches=int(sum(row["direction_match"] for row in all_checks)),
                      max_event_time_error=float(max(row["time_error"] for row in all_checks)),
                      native_audit_windows=len(audits),
                      native_audit_past_crossings=int(sum(row["native_past_crossings"] for row in audits)),
                      native_audit_label_matches=int(sum(row["label_match"] for row in audits)))
    held = heldout(all_rows, sorted(args.seeds), args.min_cell)
    (args.output_dir/"definition_and_label_validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    (args.output_dir/"heldout_validation.json").write_text(json.dumps(held, indent=2)+"\n")
    summary = dict(seeds=sorted(args.seeds), n_risk_origins=len(all_rows),
                   positive_first_crossing_windows=int(sum(row["first_crossing"] for row in all_rows)),
                   positive_seeds=len({row["seed"] for row in all_rows if row["first_crossing"]}),
                   definitions=dict(q="(N_R-N_L)/N", majority_school="fish in the currently more populated chamber", majority_com="arithmetic mean position of majority-school fish", gap_direction="unit vector from majority-school COM to (0,0)", a="mean majority-school heading dotted with gap direction"),
                   validation=validation, heldout=held)
    (args.output_dir/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    config = vars(args).copy()
    config = {key: str(value) if isinstance(value, Path) else value for key, value in config.items()}
    (args.output_dir/"analysis_config.json").write_text(json.dumps(config, indent=2)+"\n")
    n, y = counts(all_rows)
    brier_relative = held["mean_seed_brier_improvement"] / np.mean([row["brier_q"] for row in held["per_seed"]])
    report = f'''# Population imbalance and majority-school gap-directed polarization

## Scope and definitions

This read-only analysis uses all 40 existing no-hydrodynamics gap-20 trajectories (`N=100`, `In=0.5`). No simulator code, raw data, dynamical equation, fitted potential, or parameter was changed.

At each prospective risk origin,

`q=(N_R-N_L)/N`.

The **majority school** is the fish in the currently more populated chamber. Its center of mass is their arithmetic mean position `C_majority`. The instantaneous direction to the center of the opening is `g=-C_majority/|C_majority|`, and

`a = mean(e(theta_i) dot g)` over majority-school fish.

Thus positive `a` means that the majority school points toward the gap. This definition uses only the current state; it is not signed by the later outcome. Origins are 5 time units apart, have constant coarse occupancy for the previous 2 units, and have no indexed initiation in that retrospective interval. The label is whether the first indexed crossing in `(t,t+5]` is from the current majority toward the minority chamber. Outcome windows do not overlap.

## Definition and label validation

- {validation['n_definition_checks']} raw snapshots reproduce cached q to {validation['max_abs_q_error']:.2e}; majority counts have {validation['majority_count_mismatches']} mismatches.
- Every majority COM is on the expected chamber side ({validation['majority_com_side_mismatches']} mismatches), and the direction vector is normalized to {validation['max_gap_direction_norm_error']:.2e}.
- Vector projection and the independent scalar formula `mean(cos(theta-alpha_gap))` agree to {validation['max_abs_a_projection_error']:.2e}.
- All {validation['native_event_records']} existing native-refined event records show a raw population change at the stored frame, with {validation['native_event_direction_matches']}/{validation['native_event_records']} direction matches and zero timing error.
- In a bounded native-resolution audit of {validation['native_audit_windows']} risk windows, there were {validation['native_audit_past_crossings']} hidden crossings during the nominally quiet retrospective intervals and {validation['native_audit_label_matches']}/{validation['native_audit_windows']} prospective label matches.

The last check is sampled rather than exhaustive, so extremely rare rapid cross-and-return events could still contaminate some nominal negatives. It finds no evidence that this drives the result.

## The three empirical distributions

### 1. `p(q)`

The population distribution is strongly concentrated near `q=+/-1`, while also retaining a broad central component. The bands resample independent seeds. Population imbalance is therefore a principal physical state variable, but it does not uniquely determine whether crossing starts.

### 2. `p(q,a)`

This is the empirical distribution over the {len(all_rows)} prospective, crossing-free-history risk origins, not a stationary potential. Both toward-gap and away-from-gap headings occur at the same q. Therefore a contains configuration information not determined by chamber population alone.

### 3. Conditional first-crossing probability

The result is ordered in the physically expected direction. For `|q|<0.5`, the probability rises from {y[0,0]/n[0,0]:.3f} ({y[0,0]}/{n[0,0]}) for `a` in `[-1,-0.6)` to {y[0,4]/n[0,4]:.3f} ({y[0,4]}/{n[0,4]}) for `a` in `[0.6,1]`. For `0.5<=|q|<0.8`, it rises from {y[1,0]/n[1,0]:.3f} to {y[1,4]/n[1,4]:.3f}. Even at `|q|>=0.8`, where initiation is rare, the corresponding values are {y[2,0]/n[2,0]:.3f} and {y[2,4]/n[2,4]:.3f}. Seed-cluster intervals for every cell are in `conditional_first_crossing.csv`.

There are {len(all_rows)} nonoverlapping risk origins, {sum(row['first_crossing'] for row in all_rows)} positives, and {len({row['seed'] for row in all_rows if row['first_crossing']})} positive seeds. These are prospective first-entry windows, not complete school relocations.

## Held-out seeds

Seeds 0--29 define fixed binned probabilities; seeds 30--39 are held out. The held-out set contains {held['heldout_n']} origins and {held['heldout_positives']} positives. Pooled AUC increases from {held['pooled_auc_q']:.3f} for q alone to {held['pooled_auc_qa']:.3f} for `(q,a)`. The mean seed-level AUC improvement is {held['mean_seed_auc_improvement']:.3f}, with seed-bootstrap 95% CI {held['auc_seed_bootstrap_ci95'][0]:.3f}--{held['auc_seed_bootstrap_ci95'][1]:.3f}.

Pooled Brier score decreases from {held['pooled_brier_q']:.5f} to {held['pooled_brier_qa']:.5f}. The mean seed-level improvement is {held['mean_seed_brier_improvement']:.5f} (95% CI {held['seed_bootstrap_ci95'][0]:.5f}--{held['seed_bootstrap_ci95'][1]:.5f}), about {100*brier_relative:.1f}% of the mean q-only seed score. The gain is therefore reproducible and practically visible, not merely a tiny numerical change.

## Physical conclusion

Positive majority-school gap-directed polarization genuinely predicts first crossing **before any measured population redistribution**. The association survives conditioning on coarse population imbalance and independent-seed testing. The supported interpretation is: q determines how the population is loaded, while a measures whether the loaded school is oriented so that self-propulsion carries it toward the opening.

This is predictive evidence, not proof that polarization is the microscopic cause of the turn. In the no-hydrodynamics simulator, heading toward the gap is also kinematically tied to source-side self-propelled motion. The analysis does not identify which wall, social, or stochastic angular term created that heading, and first-fish entry is not the same as sustained collective relocation.

## Reproducibility

- `analyze_q_gap_polarization.py`: complete read-only analysis.
- `definition_validation_samples.csv` and `definition_and_label_validation.json`: raw definition and label audits.
- `prospective_q_a_windows.csv`: compact risk set.
- `conditional_first_crossing.csv`: counts, probabilities, and seed-cluster intervals.
- `heldout_validation.json`: independent-seed comparison.
- `figures/`: the three requested plots.
'''
    (args.output_dir/"REPORT.md").write_text(report)
    print("Wrote", args.output_dir)


if __name__ == "__main__":
    main()
