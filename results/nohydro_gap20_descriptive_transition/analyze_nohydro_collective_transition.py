#!/usr/bin/env python3
"""Descriptive collective-crossing analysis for processed no-hydro data.

This script deliberately uses only the saved occupancy and accepted-event
products.  It does not infer polarization, geometry, or a physical mechanism
that is absent from those files.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


FILE_RE = re.compile(
    r"event_histogram_If00_In05_num100_R50_gap(?P<gap>\d+)_n(?P<seed>\d+)_BEM_no_hydro\.mat$"
)


def load_seed(path: Path) -> dict:
    match = FILE_RE.match(path.name)
    if match is None:
        raise ValueError(f"Unexpected filename: {path.name}")
    data = loadmat(
        path,
        variable_names=["num_fish_left", "event_occur", "event_size", "para"],
        squeeze_me=True,
        struct_as_record=False,
    )
    para = data["para"]
    n = int(para.num)
    dt = float(para.dt)
    nl = np.asarray(data["num_fish_left"], dtype=np.int16).reshape(-1)
    events = np.atleast_1d(data["event_occur"]).astype(np.int64)
    marks = np.atleast_1d(data["event_size"]).astype(np.int16)
    if len(events) != len(marks):
        raise ValueError(f"Event/mark mismatch in {path}")
    return {
        "gap": int(match.group("gap")),
        "seed": int(match.group("seed")),
        "n": n,
        "dt": dt,
        "nl": nl,
        "q": 1.0 - 2.0 * nl / n,
        "events": events,
        "marks": marks,
    }


def dominance_relocations(q: np.ndarray, dt: float, start: int, threshold: float) -> list[dict]:
    """Return first passages between the two threshold-defined dominance basins."""
    labels = np.where(q >= threshold, 1, np.where(q <= -threshold, -1, 0))
    occupied = np.flatnonzero(labels != 0)
    if occupied.size == 0:
        return []
    out = []
    last_sign = int(labels[occupied[0]])
    last_basin_index = int(occupied[0])
    last_arrival = int(occupied[0])
    for idx in occupied[1:]:
        sign = int(labels[idx])
        if sign == last_sign:
            last_basin_index = int(idx)
            continue
        if idx >= start and last_basin_index >= start:
            out.append(
                {
                    "origin": last_sign,
                    "destination": sign,
                    "departure_index": last_basin_index,
                    "arrival_index": int(idx),
                    "duration": (int(idx) - last_basin_index) * dt,
                    "time_since_previous_arrival": (int(idx) - last_arrival) * dt,
                }
            )
        last_sign = sign
        last_arrival = int(idx)
        last_basin_index = int(idx)
    return out


def seed_bootstrap(values: np.ndarray, rng: np.random.Generator, nboot: int = 5000) -> tuple[float, float]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan, np.nan
    draws = rng.choice(values, size=(nboot, values.size), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--gaps", type=int, nargs="+", default=[10, 20])
    parser.add_argument("--transient", type=float, default=200.0)
    parser.add_argument("--dominance-threshold", type=float, default=0.6)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = args.output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(20261109)

    ensembles: dict[int, list[dict]] = {}
    for gap in args.gaps:
        paths = sorted(
            args.processed_dir.glob(
                f"event_histogram_If00_In05_num100_R50_gap{gap}_n*_BEM_no_hydro.mat"
            ),
            key=lambda p: int(FILE_RE.match(p.name).group("seed")),
        )
        ensembles[gap] = [load_seed(path) for path in paths]

    per_seed = []
    relocation_rows = []
    event_rows = []
    summary_rows = []
    aligned_by_gap: dict[int, list[np.ndarray]] = {}

    for gap, seeds in ensembles.items():
        aligned_by_gap[gap] = []
        for record in seeds:
            start = int(round(args.transient / record["dt"]))
            duration = (len(record["q"]) - start) * record["dt"]
            keep = record["events"] >= start
            events = record["events"][keep]
            marks = record["marks"][keep]
            nonzero_marks = marks[marks != 0]
            relocs = dominance_relocations(
                record["q"], record["dt"], start, args.dominance_threshold
            )
            for relocation in relocs:
                relocation_rows.append({"gap": gap, "seed": record["seed"], **relocation})
            for event, mark in zip(events, marks):
                event_rows.append(
                    {
                        "gap": gap,
                        "seed": record["seed"],
                        "time": event * record["dt"],
                        "event_size": int(mark),
                        "abs_event_size": abs(int(mark)),
                        "q_at_marker": float(record["q"][event]),
                    }
                )
                half = int(round(20.0 / record["dt"]))
                stride = int(round(0.1 / record["dt"]))
                if event - half >= start and event + half < len(record["q"]):
                    aligned_by_gap[gap].append(record["q"][event - half : event + half + 1 : stride])
            per_seed.append(
                {
                    "gap": gap,
                    "seed": record["seed"],
                    "exposure": duration,
                    "accepted_events": len(events),
                    "nonzero_cluster_marks": len(nonzero_marks),
                    "event_rate": len(events) / duration,
                    "mean_abs_event_size": float(np.mean(np.abs(nonzero_marks))) if len(nonzero_marks) else np.nan,
                    "median_abs_event_size": float(np.median(np.abs(nonzero_marks))) if len(nonzero_marks) else np.nan,
                    "dominance_relocations": len(relocs),
                    "relocation_rate": len(relocs) / duration,
                    "fraction_abs_q_gt_0p8": float(np.mean(np.abs(record["q"][start::100]) > 0.8)),
                    "fraction_abs_q_lt_0p2": float(np.mean(np.abs(record["q"][start::100]) < 0.2)),
                }
            )

        gap_rows = [row for row in per_seed if row["gap"] == gap]
        gap_events = [row for row in event_rows if row["gap"] == gap]
        gap_relocs = [row for row in relocation_rows if row["gap"] == gap]
        for metric in [
            "event_rate",
            "mean_abs_event_size",
            "relocation_rate",
            "fraction_abs_q_gt_0p8",
            "fraction_abs_q_lt_0p2",
        ]:
            values = np.array([row[metric] for row in gap_rows], float)
            lo, hi = seed_bootstrap(values, rng)
            summary_rows.append(
                {
                    "gap": gap,
                    "metric": metric,
                    "estimate": float(np.nanmean(values)),
                    "ci95_low_seed_bootstrap": lo,
                    "ci95_high_seed_bootstrap": hi,
                    "n_seeds": len(gap_rows),
                }
            )
        summary_rows.extend(
            [
                {"gap": gap, "metric": "accepted_events_total", "estimate": len(gap_events), "ci95_low_seed_bootstrap": "", "ci95_high_seed_bootstrap": "", "n_seeds": len(gap_rows)},
                {"gap": gap, "metric": "dominance_relocations_total", "estimate": len(gap_relocs), "ci95_low_seed_bootstrap": "", "ci95_high_seed_bootstrap": "", "n_seeds": len(gap_rows)},
                {"gap": gap, "metric": "median_relocation_duration", "estimate": float(np.median([row["duration"] for row in gap_relocs])) if gap_relocs else np.nan, "ci95_low_seed_bootstrap": "", "ci95_high_seed_bootstrap": "", "n_seeds": len(gap_rows)},
            ]
        )

    write_csv(args.output_dir / "per_seed_summary.csv", per_seed)
    write_csv(args.output_dir / "accepted_cluster_events.csv", event_rows)
    write_csv(args.output_dir / "dominance_relocations.csv", relocation_rows)
    write_csv(args.output_dir / "summary.csv", summary_rows)

    # Stationary occupancy distributions and representative traces.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    bins = np.linspace(-1.01, 1.01, 52)
    for gap, seeds in ensembles.items():
        pooled = []
        for record in seeds:
            start = int(round(args.transient / record["dt"]))
            pooled.append(record["q"][start::100])
        axes[0].hist(np.concatenate(pooled), bins=bins, density=True, histtype="step", lw=2, label=f"gap {gap}")
        example = max(seeds, key=lambda r: len(dominance_relocations(r["q"], r["dt"], int(args.transient/r["dt"]), args.dominance_threshold)))
        step = int(round(0.1 / example["dt"]))
        t = np.arange(len(example["q"]))[::step] * example["dt"]
        axes[1].plot(t, example["q"][::step], lw=0.8, label=f"gap {gap}, seed {example['seed']}")
    axes[0].set(xlabel=r"chamber imbalance $q=(N_R-N_L)/N$", ylabel="density", title="Post-transient occupancy")
    axes[1].axhline(args.dominance_threshold, color="0.5", ls="--", lw=1)
    axes[1].axhline(-args.dominance_threshold, color="0.5", ls="--", lw=1)
    axes[1].set(xlabel="time", ylabel="q", title="Most-switching seed")
    for ax in axes: ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "occupancy_pdf_and_traces.png", dpi=220)
    plt.close(fig)

    # Accepted-event size and event-aligned q. Markers are cluster events, not full relocations.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for gap in args.gaps:
        sizes = np.array([row["abs_event_size"] for row in event_rows if row["gap"] == gap])
        sizes = sizes[sizes > 0]
        if sizes.size:
            axes[0].plot(np.arange(1, 101), [np.mean(sizes >= k) for k in range(1, 101)], label=f"gap {gap}")
        aligned = np.asarray(aligned_by_gap[gap])
        directional = []
        for record in ensembles[gap]:
            start = int(round(args.transient / record["dt"]))
            half = int(round(20.0 / record["dt"]))
            stride = int(round(0.1 / record["dt"]))
            for event, mark in zip(record["events"], record["marks"]):
                if mark != 0 and event - half >= start and event + half < len(record["q"]):
                    # Positive mark means N_L ultimately increased and q moved left.
                    direction = -np.sign(mark)
                    directional.append(direction * record["q"][event-half:event+half+1:stride])
        directional = np.asarray(directional)
        if directional.size:
            time = np.linspace(-20, 20, directional.shape[1])
            axes[1].plot(time, np.mean(directional, axis=0), lw=2, label=f"gap {gap}; n={len(directional)}")
    axes[0].set(xlabel="accepted net cluster size k", ylabel=r"$P(|\Delta N_L|\geq k)$", yscale="log", ylim=(5e-3, 1.1), title="Cluster-mark survival")
    axes[1].axvline(0, color="0.5", ls="--", lw=1)
    axes[1].set(xlabel="time from accepted marker", ylabel="direction-aligned q", title="Occupancy change around nonzero cluster marks")
    for ax in axes: ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "cluster_marks_and_alignment.png", dpi=220)
    plt.close(fig)

    # Per-seed counts show the actual independent replication and heterogeneity.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    width = 0.38
    seed_ids = np.arange(max(row["seed"] for row in per_seed) + 1)
    for offset, gap in zip([-width / 2, width / 2], args.gaps):
        rows = {row["seed"]: row for row in per_seed if row["gap"] == gap}
        axes[0].bar(seed_ids + offset, [rows[s]["accepted_events"] for s in seed_ids], width, label=f"gap {gap}")
        axes[1].bar(seed_ids + offset, [rows[s]["dominance_relocations"] for s in seed_ids], width, label=f"gap {gap}")
    axes[0].set(xlabel="independent seed", ylabel="accepted cluster events", title="Accepted events by seed")
    axes[1].set(xlabel="independent seed", ylabel="dominance relocations", title=fr"Full relocations, $|q|={args.dominance_threshold}$")
    for ax in axes: ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "seed_level_event_counts.png", dpi=220)
    plt.close(fig)

    config = {
        "gaps": args.gaps,
        "hydrodynamics": False,
        "noise": 0.5,
        "population": 100,
        "transient_time": args.transient,
        "dominance_threshold": args.dominance_threshold,
        "event_definition": "archived GetTransitEvent marker; 5-time-unit persistence check and 2-time-unit merge rule",
        "warning": "Processed files contain occupancy only; no position, heading, polarization, or cluster identity.",
    }
    (args.output_dir / "analysis_configuration.json").write_text(json.dumps(config, indent=2) + "\n")


if __name__ == "__main__":
    main()
