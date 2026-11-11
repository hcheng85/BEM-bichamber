#!/usr/bin/env python3
"""Prospective initiation/propagation analysis for no-hydro gap-20 crossings.

Raw trajectories are read-only MATLAB-v5 memmaps.  No simulator code or raw
data are modified.  The analysis distinguishes pre-flux heading/current from
post-initiation neck transport.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr
from scipy.io import loadmat


HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "theory" / "circle_collective_first_stage"))
from matv5_memmap import inspect_numeric_arrays  # noqa: E402


class RawReader:
    """Direct fixed-offset reader for uncompressed MATLAB-v5 numeric arrays."""
    def __init__(self, path: Path):
        self.path = path
        self.info = inspect_numeric_arrays(path)
        self.fd = os.open(path, os.O_RDONLY)

    def frame(self, name: str, index: int) -> np.ndarray:
        info = self.info[name]
        frame_shape = info.shape[:-1]
        frame_bytes = int(np.prod(frame_shape)) * info.dtype.itemsize
        raw = os.pread(self.fd, frame_bytes, info.data_offset + index * frame_bytes)
        return np.frombuffer(raw, dtype=info.dtype).reshape(frame_shape, order="F")

    def state_block(self, first: int, last_inclusive: int) -> np.ndarray:
        info = self.info["state"]
        frame_bytes = int(np.prod(info.shape[:-1])) * info.dtype.itemsize
        count = last_inclusive - first + 1
        raw = os.pread(self.fd, count * frame_bytes, info.data_offset + first * frame_bytes)
        return np.frombuffer(raw, dtype=info.dtype).reshape((*info.shape[:-1], count), order="F")

    def close(self) -> None:
        os.close(self.fd)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, bool)
    score = np.asarray(score, float)
    ok = np.isfinite(score)
    y, score = y[ok], score[ok]
    n1, n0 = int(y.sum()), int((~y).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = rankdata(score)
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def extract_seed(path: Path, processed: Path | None, cache: Path, stride: int, gap: float) -> dict[str, np.ndarray]:
    if cache.exists():
        with np.load(cache) as data:
            return {key: data[key] for key in data.files}
    reader = RawReader(path)
    ntime = reader.info["state"].shape[2]
    steps = np.arange(0, ntime - 1, stride, dtype=np.int64)
    # Only chamber occupancy is required continuously.  At stride=100 this is
    # a 0.1-time-unit risk grid; full state/velocity snapshots are read later
    # only at event leads and matched controls.
    if processed is not None and processed.exists():
        native_nl = np.asarray(loadmat(processed, variable_names=["num_fish_left"])["num_fish_left"]).reshape(-1)
        nl = native_nl[:ntime - 1:stride].astype(np.int16)
    else:
        nl = np.empty(len(steps), np.int16)
        # Files without a processed occupancy index are read in contiguous
        # blocks.  This avoids one mechanical-disk seek per 0.1-time-unit
        # sample while keeping memory bounded (~24 MB per 10,000-frame block).
        block_frames = 10_000
        filled = 0
        for first in range(0, ntime - 1, block_frames):
            last = min(ntime - 2, first + block_frames - 1)
            block = reader.state_block(first, last)
            selected = steps[(steps >= first) & (steps <= last)] - first
            count = len(selected)
            nl[filled:filled + count] = np.sum(block[0, :, selected] < 0, axis=1)
            filled += count
        if filled != len(steps):
            raise RuntimeError(f"occupancy scan filled {filled} of {len(steps)} samples for {path}")
    reader.close()
    arrays = {"step": steps, "nl": nl, "q": 1.0 - 2.0 * nl / 100}
    np.savez_compressed(cache, **arrays)
    return arrays


def snapshot_features(reader: RawReader, idx_native: int, direction: int, gap: float) -> tuple[float, float, float, float, int]:
    """Polarization and exact currents at one stored microscopic state."""
    st = reader.frame("state", idx_native)
    vel = reader.frame("rdot", idx_native + 1)
    x, y = st[0], st[1]
    s = float(direction)
    source = s * x < 0
    approach = source & ((-s * x) <= gap) & (np.abs(y) <= gap / 2)
    nsource = int(source.sum())
    a = float(np.sum(s * np.cos(st[2]) * source) / nsource) if nsource else np.nan
    q = float(1.0 - 2.0 * np.sum(x < 0) / x.size)
    return q, a, float(np.sum(s * vel[0] * source) / x.size), float(np.sum(s * vel[0] * approach) / x.size), int(approach.sum())


def episodes(nl: np.ndarray, dt: float, quiet: float, transient: float, follow: float) -> list[dict]:
    changes = np.flatnonzero(np.diff(nl) != 0) + 1
    if not len(changes):
        return []
    max_gap = int(round(quiet / dt))
    groups: list[list[int]] = [[int(changes[0])]]
    for idx in changes[1:]:
        if int(idx) - groups[-1][-1] <= max_gap:
            groups[-1].append(int(idx))
        else:
            groups.append([int(idx)])
    out = []
    start_min = int(round(transient / dt))
    follow_n = int(round(follow / dt))
    for group_index, group in enumerate(groups):
        first, last = group[0], group[-1]
        if first < start_min or first + follow_n >= len(nl):
            continue
        dq_first = -2.0 * (int(nl[first]) - int(nl[first - 1])) / 100.0
        if dq_first == 0:
            continue
        s = int(np.sign(dq_first))
        q_pre = 1.0 - 2.0 * int(nl[first - 1]) / 100.0
        q5 = 1.0 - 2.0 * int(nl[first + follow_n]) / 100.0
        directional_fish = s * (q5 - q_pre) * 100 / 2
        out.append({
            "first": first,
            "last": last,
            "direction": s,
            "q_pre": q_pre,
            "duration": (last - first) * dt,
            "gross_changes": int(np.sum(np.abs(np.diff(nl[first - 1:last + 1])))),
            "net_fish_episode": float(s * (int(nl[first - 1]) - int(nl[last]))),
            "net_fish_5": float(directional_fish),
            "continuation": int(directional_fish > 0),
            "reversal": int(directional_fish <= 0),
            "pre_quiet_duration": (first - (groups[group_index - 1][-1] if group_index > 0 else 0)) * dt,
        })
    return out


def cluster_bootstrap_difference(rows: list[dict], value: str, rng: np.random.Generator, nboot: int = 4000) -> tuple[float, float, float]:
    by_seed: dict[int, list[float]] = {}
    for row in rows:
        if np.isfinite(float(row[value])):
            seed_key = "seed" if "seed" in row else "event_seed"
            by_seed.setdefault(int(row[seed_key]), []).append(float(row[value]))
    means = np.array([np.mean(v) for v in by_seed.values()])
    if not len(means):
        return np.nan, np.nan, np.nan
    draws = rng.choice(means, (nboot, len(means)), replace=True).mean(axis=1)
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return float(means.mean()), float(lo), float(hi)


def cluster_bootstrap_spearman(rows: list[dict], xname: str, yname: str, rng: np.random.Generator, nboot: int = 3000) -> tuple[float, float, float]:
    seeds = np.array(sorted({int(r["seed"]) for r in rows}))
    x = np.array([float(r[xname]) for r in rows])
    y = np.array([float(r[yname]) for r in rows])
    sid = np.array([int(r["seed"]) for r in rows])
    rho = float(spearmanr(x, y, nan_policy="omit").statistic)
    vals = []
    for _ in range(nboot):
        chosen = rng.choice(seeds, len(seeds), replace=True)
        xx, yy = [], []
        for seed in chosen:
            mask = sid == seed
            xx.extend(x[mask]); yy.extend(y[mask])
        val = spearmanr(xx, yy, nan_policy="omit").statistic
        if np.isfinite(val): vals.append(val)
    lo, hi = np.quantile(vals, [0.025, 0.975]) if vals else (np.nan, np.nan)
    return rho, float(lo), float(hi)


def cluster_bootstrap_pair_auc(rows: list[dict], name: str, rng: np.random.Generator, nboot: int = 3000) -> tuple[float, float, float]:
    seeds = np.array(sorted({int(r["event_seed"]) for r in rows}))
    def score(sampled: np.ndarray) -> float:
        ev, ct = [], []
        for seed in sampled:
            subset = [r for r in rows if int(r["event_seed"]) == seed]
            ev.extend(float(r[f"event_{name}"]) for r in subset)
            ct.extend(float(r[f"control_{name}"]) for r in subset)
        return auc_score(np.r_[np.ones(len(ev)), np.zeros(len(ct))], np.r_[ev, ct])
    estimate = score(seeds)
    draws = [score(rng.choice(seeds, len(seeds), replace=True)) for _ in range(nboot)]
    lo, hi = np.nanquantile(draws, [0.025, 0.975])
    return estimate, float(lo), float(hi)


def cluster_bootstrap_binary_auc(rows: list[dict], value: str, outcome: str, rng: np.random.Generator, nboot: int = 3000) -> tuple[float, float, float]:
    seeds = np.array(sorted({int(r["seed"]) for r in rows}))
    def score(sampled: np.ndarray) -> float:
        selected = []
        for seed in sampled:
            selected.extend(r for r in rows if int(r["seed"]) == seed)
        return auc_score(np.array([int(r[outcome]) for r in selected]), np.array([float(r[value]) for r in selected]))
    estimate = score(seeds)
    draws = [score(rng.choice(seeds, len(seeds), replace=True)) for _ in range(nboot)]
    lo, hi = np.nanquantile(draws, [0.025, 0.975])
    return estimate, float(lo), float(hi)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--processed-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/initiation_propagation"))
    parser.add_argument("--seeds", type=int, nargs="*", default=list(range(40)))
    parser.add_argument("--stride", type=int, default=100)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--gap", type=float, default=2.8635642126552705)
    parser.add_argument("--transient", type=float, default=200.0)
    parser.add_argument("--quiet", type=float, default=2.0)
    parser.add_argument("--follow", type=float, default=5.0)
    parser.add_argument("--index-only", action="store_true", help="build bounded occupancy caches and stop")
    args = parser.parse_args()
    out = args.output_dir
    cache_dir, fig_dir = out / "cache", out / "figures"
    cache_dir.mkdir(parents=True, exist_ok=True); fig_dir.mkdir(exist_ok=True)
    dt = args.dt * args.stride
    rng = np.random.default_rng(20261110)
    leads = (5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.2, 0.1)

    all_data: dict[int, dict[str, np.ndarray]] = {}
    raw_maps: dict[int, RawReader] = {}
    event_rows: list[dict] = []
    for seed in args.seeds:
        raw = args.raw_dir / f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"
        if not raw.exists():
            raise FileNotFoundError(raw)
        print(f"extract seed {seed}", flush=True)
        processed = None if args.processed_dir is None else args.processed_dir / f"event_histogram_If00_In05_num100_R50_gap20_n{seed}_BEM_no_hydro.mat"
        data = extract_seed(raw, processed, cache_dir / f"seed{seed:02d}_stride{args.stride}_qscan_v4.npz", args.stride, args.gap)
        all_data[seed] = data
        if args.index_only:
            continue
        reader = RawReader(raw)
        raw_maps[seed] = reader
        for event_id, ep in enumerate(episodes(data["nl"], dt, args.quiet, args.transient, args.follow)):
            first = ep["first"]
            lo_native = int(data["step"][first - 1])
            hi_native = int(data["step"][first])
            xlocal = reader.state_block(lo_native, hi_native)[0]
            nlocal = np.sum(xlocal < 0, axis=0)
            changed = np.flatnonzero(nlocal != nlocal[0])
            if not len(changed):
                continue
            first_native = lo_native + int(changed[0])
            row = {"seed": seed, "event_id": event_id, "time": first_native * args.dt, "first_native": first_native, **ep}
            for lead in leads:
                if ep["pre_quiet_duration"] < lead or first_native - int(round(lead / args.dt)) < int(round(args.transient / args.dt)):
                    continue
                idx_native = first_native - int(round(lead / args.dt))
                tag = str(lead).replace(".", "p")
                qv, av, js, ja, na = snapshot_features(reader, idx_native, ep["direction"], args.gap)
                row[f"q_m{tag}"] = qv; row[f"a_m{tag}"] = av
                row[f"j_source_m{tag}"] = js; row[f"j_approach_m{tag}"] = ja
                row[f"n_approach_m{tag}"] = na
                source_fraction = (1 - ep["direction"] * qv) / 2
                row[f"collision_current_residual_m{tag}"] = js - source_fraction * av
            for window in (0.2, 0.5, 1.0):
                end_native = first_native + int(round(window / args.dt))
                q0 = float(1 - 2*np.sum(reader.frame("state", first_native-1)[0] < 0)/100)
                q1 = float(1 - 2*np.sum(reader.frame("state", end_native)[0] < 0)/100)
                wtag = str(window).replace('.', 'p')
                current = ep["direction"] * (q1 - q0) / (2 * window)
                row[f"j_neck_0_{wtag}"] = current
                row[f"net_fish_remaining_after_{wtag}"] = ep["net_fish_5"] - current * window * 100
                row[f"duration_remaining_after_{wtag}"] = max(0.0, ep["duration"] - window)
                end_coarse = min(ep["last"], first + int(round(window / dt)))
                row[f"gross_changes_remaining_after_{wtag}"] = int(np.sum(np.abs(np.diff(data["nl"][end_coarse:ep["last"] + 1])))) if end_coarse < ep["last"] else 0
            event_rows.append(row)

    if args.index_only:
        print(f"indexed {len(all_data)} seeds", flush=True)
        return

    # Prospective controls: exact same q at each lead, different seed, no
    # occupancy change for 2 time units before and 5 after the control frame.
    matched_rows: list[dict] = []
    quiet_n, follow_n = int(round(args.quiet / dt)), int(round(args.follow / dt))
    for row in event_rows:
        for lead in leads:
            tag = str(lead).replace(".", "p")
            if f"q_m{tag}" not in row or row[f"q_m{tag}"] == "":
                continue
            q_target = float(row[f"q_m{tag}"])
            candidates = []
            for seed, data in all_data.items():
                if seed == int(row["seed"]):
                    continue
                valid = np.flatnonzero(np.isclose(data["q"], q_target, atol=1e-12))
                valid = valid[(valid > int(args.transient/dt) + quiet_n) & (valid + follow_n < len(data["q"]))]
                if len(valid):
                    # Sparse deterministic scan; require no flux around control.
                    order = rng.permutation(valid)
                    for idx in order[: min(200, len(order))]:
                        if data["nl"][idx-quiet_n] == data["nl"][idx+follow_n] and np.all(np.diff(data["nl"][idx-quiet_n:idx+follow_n+1]) == 0):
                            candidates.append((seed, int(idx)))
                            break
            if not candidates:
                continue
            control_seed, idx = candidates[int(rng.integers(len(candidates)))]
            direction = int(row["direction"])
            control = all_data[control_seed]
            control_q, control_a, control_js, control_ja, _ = snapshot_features(
                raw_maps[control_seed], int(control["step"][idx]), direction, args.gap
            )
            matched_rows.append({
                "event_seed": int(row["seed"]), "event_id": int(row["event_id"]),
                "lead": lead, "direction": direction, "q": q_target,
                "control_seed": control_seed, "control_time": idx * dt,
                "event_a": float(row[f"a_m{tag}"]),
                "control_a": control_a,
                "delta_a": float(row[f"a_m{tag}"]) - control_a,
                "event_j_source": float(row[f"j_source_m{tag}"]),
                "control_j_source": control_js,
                "delta_j_source": float(row[f"j_source_m{tag}"]) - control_js,
                "event_j_approach": float(row[f"j_approach_m{tag}"]),
                "control_j_approach": control_ja,
                "delta_j_approach": float(row[f"j_approach_m{tag}"]) - control_ja,
            })

    write_csv(out / "crossing_episodes.csv", event_rows)
    write_csv(out / "matched_precursor_pairs.csv", matched_rows)

    stats_rows: list[dict] = []
    for lead in leads:
        subset = [r for r in matched_rows if float(r["lead"]) == lead]
        for metric in ("delta_a", "delta_j_source", "delta_j_approach"):
            est, lo, hi = cluster_bootstrap_difference(subset, metric, rng)
            stats_rows.append({"analysis": "matched_precursor", "lead": lead, "metric": metric, "estimate": est, "ci95_low": lo, "ci95_high": hi, "n_pairs": len(subset), "n_event_seeds": len({r['event_seed'] for r in subset})})
        y = np.r_[np.ones(len(subset)), np.zeros(len(subset))]
        for name in ("a", "j_source", "j_approach"):
            estimate, lo, hi = cluster_bootstrap_pair_auc(subset, name, rng)
            stats_rows.append({"analysis": "matched_auc", "lead": lead, "metric": name, "estimate": estimate, "ci95_low": lo, "ci95_high": hi, "n_pairs": len(subset), "n_event_seeds": len({r['event_seed'] for r in subset})})
    for current in ("j_neck_0_0p2", "j_neck_0_0p5", "j_neck_0_1p0"):
        for outcome in ("net_fish_5", "duration", "gross_changes"):
            rho, lo, hi = cluster_bootstrap_spearman(event_rows, current, outcome, rng)
            stats_rows.append({"analysis": "propagation_spearman", "lead": "", "metric": f"{current}_vs_{outcome}", "estimate": rho, "ci95_low": lo, "ci95_high": hi, "n_pairs": len(event_rows), "n_event_seeds": len({r['seed'] for r in event_rows})})
        wtag = current.rsplit("_", 1)[-1]
        for outcome in (f"net_fish_remaining_after_{wtag}", f"duration_remaining_after_{wtag}", f"gross_changes_remaining_after_{wtag}"):
            rho, lo, hi = cluster_bootstrap_spearman(event_rows, current, outcome, rng)
            stats_rows.append({"analysis": "subsequent_propagation_spearman", "lead": "", "metric": f"{current}_vs_{outcome}", "estimate": rho, "ci95_low": lo, "ci95_high": hi, "n_pairs": len(event_rows), "n_event_seeds": len({r['seed'] for r in event_rows})})
        estimate, lo, hi = cluster_bootstrap_binary_auc(event_rows, current, "continuation", rng)
        stats_rows.append({"analysis": "continuation_auc", "lead": "", "metric": current, "estimate": estimate, "ci95_low": lo, "ci95_high": hi, "n_pairs": len(event_rows), "n_event_seeds": len({r['seed'] for r in event_rows})})
    write_csv(out / "statistical_summary.csv", stats_rows)

    # Compact uncertainty view: prospective initiation versus post-initiation propagation.
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    for name, label in (("a", "source polarization a"), ("j_approach", "approach current")):
        rows_auc = [r for r in stats_rows if r["analysis"] == "matched_auc" and r["metric"] == name]
        x = -np.array([float(r["lead"]) for r in rows_auc])
        y = np.array([float(r["estimate"]) for r in rows_auc])
        lo = np.array([float(r["ci95_low"]) for r in rows_auc]); hi = np.array([float(r["ci95_high"]) for r in rows_auc])
        order = np.argsort(x)
        axes[0].errorbar(x[order], y[order], yerr=[y[order]-lo[order], hi[order]-y[order]], marker="o", capsize=3, label=label)
    axes[0].axhline(.5, color="0.5", ls=":"); axes[0].set(xlabel="time from first crossing", ylabel="same-q event/control AUC", ylim=(.30, 1.0)); axes[0].legend(frameon=False)
    cont = [r for r in stats_rows if r["analysis"] == "continuation_auc"]
    windows = np.array([{"j_neck_0_0p2": 0.2, "j_neck_0_0p5": 0.5, "j_neck_0_1p0": 1.0}[r["metric"]] for r in cont])
    y = np.array([float(r["estimate"]) for r in cont]); lo = np.array([float(r["ci95_low"]) for r in cont]); hi = np.array([float(r["ci95_high"]) for r in cont])
    axes[1].errorbar(windows, y, yerr=[y-lo, hi-y], marker="o", capsize=3)
    axes[1].axhline(.5, color="0.5", ls=":"); axes[1].set(xlabel="post-initiation current window", ylabel="continuation/reversal AUC", ylim=(.45, 1.0))
    fig.tight_layout(); fig.savefig(fig_dir / "predictive_effects_with_seed_uncertainty.png", dpi=220); plt.close(fig)

    # Event/control precursor chronology.
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    leads_array = np.array(leads)
    for ax, name, label in zip(axes, ("a", "j_source", "j_approach"), ("source polarization a", "source current", "approach-zone current")):
        ev, ct = [], []
        for lead in leads_array:
            sub = [r for r in matched_rows if float(r["lead"]) == lead]
            ev.append(np.mean([float(r[f"event_{name}"]) for r in sub])); ct.append(np.mean([float(r[f"control_{name}"]) for r in sub]))
        ax.plot(-leads_array, ev, "o-", label="future crossing")
        ax.plot(-leads_array, ct, "o--", label="same-q no-flux control")
        ax.axvline(0, color="0.5", ls=":"); ax.set(xlabel="time from initiation", ylabel=label)
    axes[0].legend(frameon=False); fig.tight_layout(); fig.savefig(fig_dir / "matched_preflux_chronology.png", dpi=220); plt.close(fig)

    # Early flux versus later event outcomes.
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    x = np.array([r["j_neck_0_0p2"] for r in event_rows])
    for ax, outcome, label in zip(axes, ("net_fish_5", "duration", "gross_changes"), ("directional net fish at +5", "episode duration", "gross crossing count")):
        y = np.array([r[outcome] for r in event_rows])
        ax.scatter(x, y, s=12, alpha=.45); rho = spearmanr(x, y).statistic
        ax.set(xlabel="directional neck current, 0--0.2", ylabel=label, title=fr"Spearman $\rho={rho:.2f}$")
    fig.tight_layout(); fig.savefig(fig_dir / "early_current_vs_outcomes.png", dpi=220); plt.close(fig)

    # Representative time traces from independent seeds, chosen by largest net continuation.
    reps = sorted(event_rows, key=lambda r: r["net_fish_5"], reverse=True)
    selected = []
    for row in reps:
        if row["seed"] not in {r["seed"] for r in selected}:
            selected.append(row)
        if len(selected) == 6: break
    fig, axes = plt.subplots(len(selected), 1, figsize=(10, 1.9 * len(selected)), sharex=True)
    if len(selected) == 1: axes = [axes]
    for ax, row in zip(axes, selected):
        data = all_data[int(row["seed"])]; first = int(row["first"]); s = int(row["direction"])
        lo, hi = first-int(2/dt), first+int(5/dt); time = (np.arange(lo, hi)-first)*dt
        ax.plot(time, s*data["q"][lo:hi], label="direction-aligned q")
        leads_plot = -np.array([lead for lead in leads if f"a_m{str(lead).replace('.', 'p')}" in row and row[f"a_m{str(lead).replace('.', 'p')}"] != ""])
        avals = [row[f"a_m{str(abs(t)).replace('.', 'p')}"] for t in leads_plot]
        ax.plot(leads_plot, avals, "o-", label="source a", alpha=.8)
        qtrace = s * data["q"][lo:hi]
        ax.plot(time[:-1], np.diff(qtrace) / (2 * dt), label="neck current", alpha=.65)
        ax.axvline(0, color="k", ls=":"); ax.set_ylabel(f"seed {row['seed']}")
    axes[0].legend(frameon=False, ncol=3); axes[-1].set_xlabel("time from first crossing"); fig.tight_layout(); fig.savefig(fig_dir / "representative_event_traces.png", dpi=220); plt.close(fig)

    config = vars(args).copy(); config["raw_dir"] = str(config["raw_dir"]); config["output_dir"] = str(config["output_dir"])
    config["processed_dir"] = None if config["processed_dir"] is None else str(config["processed_dir"])
    config.update({
        "coarse_dt": dt,
        "source_polarization": "mean(s*cos(theta)) over fish in source chamber",
        "source_current": "s/N sum(v_x) over fish in source chamber; exact stored rdot at selected snapshots",
        "approach_zone": "source side, divider distance <= gap, |y| <= gap/2",
        "neck_current": "directional finite difference of q/2 after first crossing",
        "initiation": "first native-resolution centerline crossing after a >=2-time-unit no-crossing interval",
        "sample_scope": f"{len(args.seeds)} raw seeds; processed occupancy used when available and otherwise generated by bounded sequential raw-state indexing",
        "velocity_index": "rdot[:,:,k+1] advances state[:,:,k]",
    })
    (out / "analysis_configuration.json").write_text(json.dumps(config, indent=2) + "\n")
    residuals = {}
    for lead in ("2p0", "1p0", "0p5", "0p2", "0p1"):
        values = np.array([float(r[f"collision_current_residual_m{lead}"]) for r in event_rows])
        residuals[lead] = {"rms": float(np.sqrt(np.mean(values**2))), "mean": float(np.mean(values)), "max_abs": float(np.max(np.abs(values))), "fraction_nonzero": float(np.mean(np.abs(values) > 1e-12))}
    audit = {"n_episodes": len(event_rows), "n_event_seeds": len({r["seed"] for r in event_rows}), "n_continuations": int(sum(r["continuation"] for r in event_rows)), "n_reversals": int(sum(r["reversal"] for r in event_rows)), "matched_pairs_by_lead": {str(lead): sum(float(r["lead"]) == lead for r in matched_rows) for lead in leads}, "collision_override_current_residual": residuals}
    (out / "sample_and_identity_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    for reader in raw_maps.values():
        reader.close()


if __name__ == "__main__":
    main()
