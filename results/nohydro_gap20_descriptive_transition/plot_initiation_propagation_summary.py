#!/usr/bin/env python3
"""Render the compact uncertainty figure from saved summary tables only."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--statistics", type=Path, required=True)
    parser.add_argument("--common-cohort", type=Path, help="optional fixed-event-cohort precursor table")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read(args.statistics)
    common = read(args.common_cohort) if args.common_cohort else None

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    for name, label in (("a", "source polarization a"), ("j_approach", "approach current")):
        if common is None:
            selected = [r for r in rows if r["analysis"] == "matched_auc" and r["metric"] == name]
            x = -np.array([float(r["lead"]) for r in selected])
            y = np.array([float(r["estimate"]) for r in selected])
            lo = np.array([float(r["ci95_low"]) for r in selected])
            hi = np.array([float(r["ci95_high"]) for r in selected])
        else:
            selected = [r for r in common if r["metric"] == name]
            x = -np.array([float(r["lead"]) for r in selected])
            y = np.array([float(r["auc"]) for r in selected])
            lo = np.array([float(r["auc_ci95_low"]) for r in selected])
            hi = np.array([float(r["auc_ci95_high"]) for r in selected])
        order = np.argsort(x)
        axes[0].errorbar(x[order], y[order], yerr=[y[order] - lo[order], hi[order] - y[order]], marker="o", capsize=3, label=label)
    axes[0].axhline(0.5, color="0.5", ls=":")
    cohort_label = "fixed-cohort" if common is not None else "available-cohort"
    axes[0].set(xlabel="time from first crossing", ylabel=f"same-q event/control AUC ({cohort_label})", ylim=(0.30, 1.0))
    axes[0].legend(frameon=False)

    selected = [r for r in rows if r["analysis"] == "continuation_auc"]
    window_map = {"j_neck_0_0p2": 0.2, "j_neck_0_0p5": 0.5, "j_neck_0_1p0": 1.0}
    x = np.array([window_map[r["metric"]] for r in selected])
    y = np.array([float(r["estimate"]) for r in selected])
    lo = np.array([float(r["ci95_low"]) for r in selected])
    hi = np.array([float(r["ci95_high"]) for r in selected])
    axes[1].errorbar(x, y, yerr=[y - lo, hi - y], marker="o", capsize=3)
    axes[1].axhline(0.5, color="0.5", ls=":")
    axes[1].set(xlabel="post-initiation current window", ylabel="continuation/reversal AUC", ylim=(0.45, 1.0))

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)


if __name__ == "__main__":
    main()
