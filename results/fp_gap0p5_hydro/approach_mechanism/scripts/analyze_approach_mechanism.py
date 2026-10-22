#!/usr/bin/env python3
"""Test pre-nucleation approach, geometric access, and recruitment in gap-50 BEM data."""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


DT = 0.01
BURN = 200.0
COMMIT_COUNT = 70
ALIGN_BEFORE = 20.0
ALIGN_AFTER = 10.0
MATCH_TIME = 10.0
GAP = 7.158910531638177
R = GAP


@dataclass
class Episode:
    seed: int
    origin: int                 # -1 left, +1 right
    full_start: int
    first_cross: int
    endpoint: int
    success: bool
    committed: int | None
    full_dwell: float


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_features(directory: Path) -> dict[int, dict[str, np.ndarray]]:
    result = {}
    for seed in range(20):
        path = directory / f"gap50_hydro_approach_seed{seed:02d}_stride10.npz"
        if not path.is_file():
            raise FileNotFoundError(path)
        with np.load(path) as archive:
            result[seed] = {key: np.asarray(archive[key]) for key in archive.files
                            if key not in ("metadata_json", "seed", "stride")}
    return result


def episodes_for_seed(seed: int, data: dict[str, np.ndarray]) -> list[Episode]:
    nl = data["left_n"]
    nr = data["right_n"]
    start = int(round(BURN / DT))
    result: list[Episode] = []
    i = start
    n = len(nl)
    while i < n - 1:
        while i < n and nl[i] != 100 and nr[i] != 100:
            i += 1
        if i >= n:
            break
        origin = -1 if nl[i] == 100 else 1
        full_start = i
        while full_start > start:
            full_before = nl[full_start - 1] == 100 if origin == -1 else nr[full_start - 1] == 100
            if not full_before:
                break
            full_start -= 1
        opp = nr if origin == -1 else nl
        while i < n and opp[i] == 0:
            i += 1
        if i >= n:
            break
        first = i
        committed = None
        success = False
        while i < n:
            if opp[i] >= COMMIT_COUNT:
                committed = i
                success = True
                break
            if opp[i] == 0:
                break
            i += 1
        endpoint = min(i, n - 1)
        result.append(Episode(seed, origin, full_start, first, endpoint, success, committed,
                              (first - full_start) * DT))
        i = endpoint + 1 if endpoint == first else endpoint
        if success:
            while i < n and nl[i] != 100 and nr[i] != 100:
                i += 1
    return result


def origin_series(data: dict[str, np.ndarray], origin: int) -> dict[str, np.ndarray]:
    side = "left" if origin == -1 else "right"
    s = 1.0 if origin == -1 else -1.0  # inward normal x component
    xc = data[f"{side}_xc"]
    yc = data[f"{side}_yc"]
    hx = data[f"{side}_hx"]
    hy = data[f"{side}_hy"]
    pol = data[f"{side}_pol"]
    front_x = data[f"{side}_xmax"] if origin == -1 else data[f"{side}_xmin"]
    d = -s * xc
    d_front = -s * front_x
    inward_cos = np.divide(s * hx, pol, out=np.full_like(pol, np.nan), where=pol > 1e-12)
    psi = np.arctan2(hy, s * hx)
    projected_y = yc + d * np.tan(psi)
    admissible = (inward_cos > 0) & (np.abs(projected_y) < GAP / 2)
    opp = data["right_n"] if origin == -1 else data["left_n"]
    return {
        "d": d, "d_front": d_front, "yc": yc, "hx_in": s * hx, "hy": hy,
        "pol": pol, "psi": psi, "projected_y": projected_y,
        "admissible": admissible, "opp": opp,
        "xstd": data[f"{side}_xstd"], "ystd": data[f"{side}_ystd"],
    }


def slope(time: np.ndarray, value: np.ndarray) -> float:
    valid = np.isfinite(value)
    if valid.sum() < 3:
        return np.nan
    return float(np.polyfit(time[valid], value[valid], 1)[0])


def finite_mean(values: np.ndarray, axis: int) -> np.ndarray:
    count = np.sum(np.isfinite(values), axis=axis)
    total = np.nansum(values, axis=axis)
    return np.divide(total, count, out=np.full_like(total, np.nan, dtype=float), where=count > 0)


def episode_row(ep: Episode, data: dict[str, np.ndarray]) -> dict[str, object]:
    f = origin_series(data, ep.origin)
    pre = int(round(MATCH_TIME / DT))
    has_pre10 = ep.first_cross - ep.full_start >= pre
    lo = ep.first_cross - pre if has_pre10 else ep.full_start
    t = (np.arange(lo, ep.first_cross + 1) - ep.first_cross) * DT
    d = f["d"][lo:ep.first_cross + 1]
    yc = f["yc"][lo:ep.first_cross + 1]
    admissible = f["admissible"][lo:ep.first_cross + 1]
    entry = np.flatnonzero(admissible & np.r_[True, ~admissible[:-1]])
    lead = -t[entry[-1]] if entry.size else np.nan
    first5 = np.flatnonzero(f["opp"][ep.first_cross:ep.endpoint + 1] >= 5)
    return {
        "seed": ep.seed, "origin": "left" if ep.origin == -1 else "right",
        "first_cross_time": ep.first_cross * DT,
        "success": int(ep.success),
        "commit_time": ep.committed * DT if ep.committed is not None else np.nan,
        "episode_duration": (ep.endpoint - ep.first_cross) * DT,
        "full_dwell_before": ep.full_dwell,
        "has_full_pre10": int(has_pre10),
        "d_at_minus10": f["d"][ep.first_cross - pre] if has_pre10 else np.nan,
        "yc_at_minus10": f["yc"][ep.first_cross - pre] if has_pre10 else np.nan,
        "cospsi_at_minus10": (f["hx_in"][ep.first_cross - pre] /
                            max(f["pol"][ep.first_cross - pre], 1e-12)) if has_pre10 else np.nan,
        "d_at_cross": f["d"][ep.first_cross],
        "front_distance_at_cross": f["d_front"][ep.first_cross - 1],
        "yc_at_cross": f["yc"][ep.first_cross],
        "cospsi_at_cross": f["hx_in"][ep.first_cross] / max(f["pol"][ep.first_cross], 1e-12),
        "polarization_at_cross": f["pol"][ep.first_cross],
        "xstd_at_cross": f["xstd"][ep.first_cross],
        "ystd_at_cross": f["ystd"][ep.first_cross],
        "pre_d_slope": slope(t, d) if has_pre10 else np.nan,
        "pre_abs_y_slope": slope(t, np.abs(yc)) if has_pre10 else np.nan,
        "admissible_at_cross": int(f["admissible"][ep.first_cross]),
        "admissible_any_pre10": int(admissible.any()) if has_pre10 else np.nan,
        "last_admissible_entry_lead": lead,
        "time_first_to_five": first5[0] * DT if first5.size else np.nan,
        "time_five_to_commit": ((ep.committed - ep.first_cross) * DT - first5[0] * DT)
                               if first5.size and ep.committed is not None else np.nan,
    }


def greedy_matches(rows: list[dict[str, object]]) -> dict[int, int]:
    # Require a complete 20-time-unit empty-opposite-chamber history so the
    # event-aligned comparison cannot include remnants of a preceding switch.
    success = [i for i, r in enumerate(rows) if r["success"] == 1 and r["full_dwell_before"] >= ALIGN_BEFORE]
    failure = [i for i, r in enumerate(rows) if r["success"] == 0 and r["full_dwell_before"] >= ALIGN_BEFORE]
    keys = ("d_at_minus10", "yc_at_minus10", "cospsi_at_minus10")
    matrix = np.array([[float(r[k]) for k in keys] for r in rows])
    scale = np.nanstd(matrix[failure], axis=0)
    scale[scale == 0] = 1
    available = set(failure)
    matches = {}
    for si in success:
        candidates = [fi for fi in available if rows[fi]["seed"] == rows[si]["seed"]
                      and rows[fi]["origin"] == rows[si]["origin"]]
        if not candidates:
            candidates = list(available)
        if not candidates:
            break
        dist = [np.nansum(((matrix[fi] - matrix[si]) / scale) ** 2) for fi in candidates]
        fi = candidates[int(np.argmin(dist))]
        matches[si] = fi
        available.remove(fi)
    return matches


def aligned_rows(episodes: list[Episode], features: dict[int, dict[str, np.ndarray]],
                 selected: list[int], group: str) -> list[dict[str, object]]:
    grid = np.arange(-ALIGN_BEFORE, ALIGN_AFTER + 1e-9, 0.1)
    offsets = np.rint(grid / DT).astype(int)
    by_seed: dict[int, list[np.ndarray]] = {}
    for idx in selected:
        ep = episodes[idx]
        ids = ep.first_cross + offsets
        if ids[0] < 0 or ids[-1] >= len(features[ep.seed]["time"]):
            continue
        f = origin_series(features[ep.seed], ep.origin)
        arr = np.vstack([f["d"][ids], f["d_front"][ids], f["yc"][ids],
                         f["hx_in"][ids] / np.maximum(f["pol"][ids], 1e-12),
                         f["pol"][ids], f["opp"][ids], f["admissible"][ids].astype(float)])
        # Do not display the previous split state before the chamber became
        # fully occupied, or later independent attempts after a failed event.
        arr[:, ids < ep.full_start] = np.nan
        if not ep.success:
            arr[:, ids > ep.endpoint] = np.nan
        by_seed.setdefault(ep.seed, []).append(arr)
    seed_means = np.stack([finite_mean(np.stack(values), axis=0) for values in by_seed.values()])
    rows = []
    for vi, name in enumerate(("wall_distance", "front_distance", "vertical_center", "cos_psi",
                               "polarization", "opposite_count", "admissible_fraction")):
        mean = np.nanmean(seed_means[:, vi], axis=0)
        if seed_means.shape[0] > 1:
            rng = np.random.default_rng(1200 + vi)
            boot = finite_mean(seed_means[rng.integers(0, seed_means.shape[0], (1000, seed_means.shape[0]))][:, :, vi], axis=1)
            low, high = np.nanpercentile(boot, [2.5, 97.5], axis=0)
        else:
            low = high = mean
        rows.extend({"group": group, "variable": name, "relative_time": t, "mean": m,
                     "ci_low": l, "ci_high": h, "independent_seeds": seed_means.shape[0]}
                    for t, m, l, h in zip(grid, mean, low, high))
    return rows


def phase_cycles(features: dict[int, dict[str, np.ndarray]], episodes: list[Episode]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    event_lookup = {seed: [ep for ep in episodes if ep.seed == seed] for seed in range(20)}
    cycle_rows: list[dict[str, object]] = []
    quality_rows: list[dict[str, object]] = []
    burn = int(round(BURN / DT))
    for seed, data in features.items():
        nl, nr = data["left_n"], data["right_n"]
        periods = []
        count_points = 0
        for origin in (-1, 1):
            full = nl == 100 if origin == -1 else nr == 100
            full[:burn] = False
            changes = np.flatnonzero(np.diff(np.r_[False, full, False]))
            for a, b in changes.reshape(-1, 2):
                if (b - a) * DT < 12:
                    continue
                f = origin_series(data, origin)
                ys = gaussian_filter1d(f["yc"][a:b], sigma=0.5 / DT, mode="nearest")
                for phase, signal in (("top", ys), ("bottom", -ys)):
                    peaks, _ = find_peaks(signal, prominence=0.75, distance=int(3 / DT))
                    if len(peaks) > 1:
                        periods.extend(np.diff(peaks) * DT)
                    for p0, p1 in zip(peaks[:-1], peaks[1:]):
                        i0, i1 = a + p0, a + p1
                        period = (i1 - i0) * DT
                        if not (3 <= period <= 40):
                            continue
                        count_points += 1
                        next_event = next((ep for ep in event_lookup[seed]
                                           if ep.origin == origin and i1 < ep.first_cross <= i1 + int(period / DT)), None)
                        cycle_rows.append({
                            "seed": seed, "origin": "left" if origin == -1 else "right", "phase": phase,
                            "time": i1 * DT, "period": period, "d_n": f["d"][i0],
                            "d_n1": f["d"][i1], "delta_d": f["d"][i1] - f["d"][i0],
                            "crossing_next_cycle": int(next_event is not None),
                            "successful_crossing_next_cycle": int(next_event.success) if next_event else 0,
                        })
        p = np.asarray(periods)
        quality_rows.append({"seed": seed, "phase_pairs": count_points, "median_period": np.nanmedian(p) if p.size else np.nan,
                             "period_cv": np.nanstd(p) / np.nanmean(p) if p.size else np.nan})
    return cycle_rows, quality_rows


def seed_bootstrap_difference(rows: list[dict[str, object]], variable: str,
                              group_key: str = "success") -> dict[str, float]:
    seeds = sorted({int(r["seed"]) for r in rows})
    def calc(chosen: list[int]) -> float:
        vals = {0: [], 1: []}
        for seed in chosen:
            for g in (0, 1):
                x = [float(r[variable]) for r in rows if int(r["seed"]) == seed and int(r[group_key]) == g
                     and np.isfinite(float(r[variable]))]
                if x:
                    vals[g].append(np.mean(x))
        return float(np.mean(vals[1]) - np.mean(vals[0])) if vals[0] and vals[1] else np.nan
    observed = calc(seeds)
    rng = np.random.default_rng(9421)
    boot = np.array([calc(list(rng.choice(seeds, len(seeds), replace=True))) for _ in range(2000)])
    return {"difference_success_minus_failure": observed,
            "ci_low": float(np.nanpercentile(boot, 2.5)), "ci_high": float(np.nanpercentile(boot, 97.5))}


def plot_aligned(rows: list[dict[str, object]], output: Path) -> None:
    variables = ("wall_distance", "vertical_center", "cos_psi", "opposite_count")
    labels = ("center distance to divider", "$Y_c$", "$\\cos\\psi$", "fish in opposite chamber")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True)
    for ax, variable, label in zip(axes.flat, variables, labels):
        for group, color in (("successful", "#0072B2"), ("matched_failure", "#D55E00")):
            sub = [r for r in rows if r["variable"] == variable and r["group"] == group]
            t = np.array([r["relative_time"] for r in sub])
            m = np.array([r["mean"] for r in sub])
            lo = np.array([r["ci_low"] for r in sub])
            hi = np.array([r["ci_high"] for r in sub])
            ax.plot(t, m, color=color, label=group.replace("_", " "))
            ax.fill_between(t, lo, hi, color=color, alpha=.18)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_ylabel(label)
        ax.grid(alpha=.2)
    axes[1, 0].set_xlabel("time from first-fish crossing")
    axes[1, 1].set_xlabel("time from first-fish crossing")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Pre-nucleation approach and post-nucleation recruitment")
    for ext in ("png", "pdf"):
        fig.savefig(output / f"02_event_aligned_histories.{ext}", dpi=220)
    plt.close(fig)


def plot_geometry(rows: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, success, title in zip(axes, (0, 1), ("failed nucleation", "successful nucleation")):
        sub = [r for r in rows if int(r["success"]) == success]
        d = np.array([r["d_at_cross"] for r in sub])
        psi = np.arccos(np.clip([r["cospsi_at_cross"] for r in sub], -1, 1))
        adm = np.array([r["admissible_at_cross"] for r in sub], bool)
        ax.scatter(d[~adm], psi[~adm], s=10, alpha=.35, color="#999999", label="outside COM criterion")
        ax.scatter(d[adm], psi[adm], s=12, alpha=.55, color="#009E73", label="inside COM criterion")
        ax.set_title(f"{title} (n={len(sub)})")
        ax.set_xlabel("center distance $d$")
        ax.set_ylabel("$|\\psi|$ [rad]")
        ax.set_ylim(0, np.pi)
        ax.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(output / f"03_geometric_condition.{ext}", dpi=220)
    plt.close(fig)


def plot_cycles(rows: list[dict[str, object]], quality: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for group, color, label in ((0, "#777777", "not successful next cycle"), (1, "#CC79A7", "successful next cycle")):
        vals = np.array([r["delta_d"] for r in rows if r["successful_crossing_next_cycle"] == group])
        axes[0].hist(vals, bins=np.linspace(-4, 4, 41), density=True, histtype="step", lw=2, color=color, label=f"{label} (n={len(vals)})")
    axes[0].axvline(0, color="k", lw=1)
    axes[0].set_xlabel("cycle-to-cycle $\\Delta d$")
    axes[0].set_ylabel("density")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=.2)
    axes[1].scatter([r["median_period"] for r in quality], [r["period_cv"] for r in quality], color="#0072B2")
    axes[1].axhline(.5, color="k", ls="--", lw=1, label="CV = 0.5")
    axes[1].set_xlabel("median detected vertical period")
    axes[1].set_ylabel("within-seed period CV")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=.2)
    for ext in ("png", "pdf"):
        fig.savefig(output / f"04_cycle_statistics.{ext}", dpi=220)
    plt.close(fig)


def plot_recruitment(rows: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    success = [r for r in rows if r["success"] == 1]
    fail = [r for r in rows if r["success"] == 0]
    axes[0].hist([r["episode_duration"] for r in fail], bins=np.linspace(0, 20, 41), color="#D55E00", alpha=.7)
    axes[0].set_xlabel("failed first-fish excursion duration")
    axes[0].set_ylabel("episodes")
    a = np.array([r["time_first_to_five"] for r in success], float)
    b = np.array([r["time_five_to_commit"] for r in success], float)
    valid = np.isfinite(a) & np.isfinite(b)
    axes[1].scatter(a[valid], b[valid], s=14, alpha=.5, color="#0072B2")
    axes[1].set_xlabel("first fish to five fish")
    axes[1].set_ylabel("five fish to 70-fish commitment")
    axes[1].grid(alpha=.2)
    for ext in ("png", "pdf"):
        fig.savefig(output / f"05_recruitment_timing.{ext}", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    (output / "tables").mkdir(parents=True, exist_ok=True)
    (output / "figures").mkdir(parents=True, exist_ok=True)
    features = load_features(args.feature_dir.expanduser().resolve())
    episodes = [ep for seed, data in features.items() for ep in episodes_for_seed(seed, data)]
    rows = [episode_row(ep, features[ep.seed]) for ep in episodes]
    matches = greedy_matches(rows)
    aligned = aligned_rows(episodes, features, list(matches), "successful")
    aligned += aligned_rows(episodes, features, list(matches.values()), "matched_failure")
    cycles, cycle_quality = phase_cycles(features, episodes)
    write_csv(output / "tables" / "nucleation_episodes.csv", rows)
    write_csv(output / "tables" / "event_aligned_histories.csv", aligned)
    write_csv(output / "tables" / "cycle_pairs.csv", cycles)
    write_csv(output / "tables" / "cycle_quality.csv", cycle_quality)
    write_csv(output / "tables" / "matched_controls.csv",
              [{"success_row": si, "failure_row": fi, "success_seed": rows[si]["seed"], "failure_seed": rows[fi]["seed"]}
               for si, fi in matches.items()])
    plot_aligned(aligned, output / "figures")
    plot_geometry(rows, output / "figures")
    plot_cycles(cycles, cycle_quality, output / "figures")
    plot_recruitment(rows, output / "figures")

    train = [r for r in rows if int(r["seed"]) < 15]
    test = [r for r in rows if int(r["seed"]) >= 15]
    def group_summary(group: list[dict[str, object]]) -> dict[str, object]:
        s = [r for r in group if r["success"] == 1]
        f = [r for r in group if r["success"] == 0]
        return {
            "episodes": len(group), "successful": len(s), "failed": len(f),
            "success_fraction": len(s) / len(group) if group else np.nan,
            "success_pre_d_slope_mean": float(np.nanmean([r["pre_d_slope"] for r in s])),
            "failure_pre_d_slope_mean": float(np.nanmean([r["pre_d_slope"] for r in f])),
            "success_pre10_inward_fraction": float(np.mean([r["pre_d_slope"] < 0 for r in s if np.isfinite(r["pre_d_slope"])])),
            "failure_pre10_inward_fraction": float(np.mean([r["pre_d_slope"] < 0 for r in f if np.isfinite(r["pre_d_slope"])])),
            "success_admissible_at_cross_fraction": float(np.mean([r["admissible_at_cross"] for r in s])),
            "failure_admissible_at_cross_fraction": float(np.mean([r["admissible_at_cross"] for r in f])),
            "success_admissible_pre10_fraction": float(np.nanmean([r["admissible_any_pre10"] for r in s])),
            "failure_admissible_pre10_fraction": float(np.nanmean([r["admissible_any_pre10"] for r in f])),
            "success_median_center_distance_at_cross": float(np.nanmedian([r["d_at_cross"] for r in s])),
            "failure_median_center_distance_at_cross": float(np.nanmedian([r["d_at_cross"] for r in f])),
            "success_median_front_distance_at_cross": float(np.nanmedian([r["front_distance_at_cross"] for r in s])),
            "failure_median_front_distance_at_cross": float(np.nanmedian([r["front_distance_at_cross"] for r in f])),
            "median_failed_excursion_duration": float(np.nanmedian([r["episode_duration"] for r in f])),
            "median_success_time_first_to_five": float(np.nanmedian([r["time_first_to_five"] for r in s])),
            "median_success_time_five_to_commit": float(np.nanmedian([r["time_five_to_commit"] for r in s])),
        }
    summary = {
        "geometry": {"R": R, "gap": GAP, "gap_over_R": 1.0, "gap_over_2R": 0.5},
        "definitions": {"burn_time": BURN, "first_crossing": "first opposite-side fish after N_origin=100",
                        "success": f"N_opposite reaches {COMMIT_COUNT} before returning to zero"},
        "train": group_summary(train), "heldout": group_summary(test),
        "seed_clustered_pre_d_slope_difference": seed_bootstrap_difference(rows, "pre_d_slope"),
        "seed_clustered_admissible_at_cross_difference": seed_bootstrap_difference(rows, "admissible_at_cross"),
        "seed_clustered_cospsi_at_cross_difference": seed_bootstrap_difference(rows, "cospsi_at_cross"),
        "seed_clustered_xstd_at_cross_difference": seed_bootstrap_difference(rows, "xstd_at_cross"),
        "cycle_pairs": len(cycles),
        "cycle_crossing_next": int(sum(r["crossing_next_cycle"] for r in cycles)),
        "cycle_delta_d_not_success_mean": float(np.nanmean([r["delta_d"] for r in cycles if not r["successful_crossing_next_cycle"]])),
        "cycle_delta_d_success_mean": float(np.nanmean([r["delta_d"] for r in cycles if r["successful_crossing_next_cycle"]])),
        "seed_clustered_cycle_delta_d_difference": seed_bootstrap_difference(
            cycles, "delta_d", "successful_crossing_next_cycle"),
        "median_period_cv": float(np.nanmedian([r["period_cv"] for r in cycle_quality])),
        "matched_pairs": len(matches),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    flat = []
    for split in ("train", "heldout"):
        flat.extend({"section": split, "metric": k, "value": v} for k, v in summary[split].items())
    flat.extend({"section": "overall", "metric": k, "value": v} for k, v in summary.items()
                if isinstance(v, (int, float)))
    for section, values in summary.items():
        if section.startswith("seed_clustered_") and isinstance(values, dict):
            flat.extend({"section": section, "metric": key, "value": value}
                        for key, value in values.items())
    write_csv(output / "summary.csv", flat)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
