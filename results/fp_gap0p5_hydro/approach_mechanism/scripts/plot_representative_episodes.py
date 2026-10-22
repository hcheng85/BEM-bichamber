#!/usr/bin/env python3
"""Plot one matched successful/failed nucleation pair using individual fish states."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--geometry-file", type=Path, required=True)
    args = parser.parse_args()
    root = args.analysis_dir.expanduser().resolve()
    episodes = list(csv.DictReader((root / "tables" / "nucleation_episodes.csv").open()))
    matches = list(csv.DictReader((root / "tables" / "matched_controls.csv").open()))
    candidates = []
    for pair in matches:
        s = episodes[int(pair["success_row"])]
        f = episodes[int(pair["failure_row"])]
        if s["seed"] == f["seed"] and float(f["episode_duration"]) >= 2.0:
            candidates.append((s, f))
    if not candidates:
        raise RuntimeError("no same-seed matched pair with a two-time-unit failed excursion")
    success_d = np.array([float(s["d_at_cross"]) for s, _ in candidates])
    success_c = np.array([float(s["cospsi_at_cross"]) for s, _ in candidates])
    score = ((success_d - np.median(success_d)) / max(np.std(success_d), 1e-9)) ** 2
    score += ((success_c - np.median(success_c)) / max(np.std(success_c), 1e-9)) ** 2
    success, failure = candidates[int(np.argmin(score))]
    seed = int(success["seed"])
    raw = args.raw_dir / f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
    loaded = sio.loadmat(raw, variable_names=["state"], squeeze_me=True)
    state = np.asarray(loaded["state"])
    geometry = sio.loadmat(args.geometry_file, variable_names=["Sx", "Sy"])
    sx = np.asarray(geometry["Sx"]).reshape(-1)
    sy = np.asarray(geometry["Sy"]).reshape(-1)
    relative = (-10.0, -2.0, 0.0, 2.0)
    fig, axes = plt.subplots(2, len(relative), figsize=(14, 7), sharex=True, sharey=True,
                             constrained_layout=True)
    selected = []
    for row, (label, event) in enumerate((("successful recruitment", success), ("failed nucleation", failure))):
        origin = -1 if event["origin"] == "left" else 1
        t0 = float(event["first_cross_time"])
        for col, rel in enumerate(relative):
            index = int(round((t0 + rel) / 0.001))
            x, y, theta = state[:, :, index]
            ax = axes[row, col]
            ax.scatter(sx, sy, s=.2, color="black", alpha=.65, rasterized=True)
            crossed = x >= 0 if origin == -1 else x < 0
            ax.scatter(x[~crossed], y[~crossed], s=10, color="#0072B2", alpha=.8)
            ax.scatter(x[crossed], y[crossed], s=13, color="#D55E00", alpha=.9)
            choose = np.arange(0, len(x), 5)
            ax.quiver(x[choose], y[choose], np.cos(theta[choose]), np.sin(theta[choose]),
                      angles="xy", scale_units="xy", scale=1.25, width=.004, color="#222222", alpha=.7)
            ax.axhline(7.158910531638177 / 2, color="#009E73", ls=":", lw=.8)
            ax.axhline(-7.158910531638177 / 2, color="#009E73", ls=":", lw=.8)
            ax.set_aspect("equal")
            ax.set_xlim(-7.5, 7.5)
            ax.set_ylim(-7.5, 7.5)
            if row == 0:
                ax.set_title(f"$t-t_A={rel:g}$")
            if col == 0:
                ax.set_ylabel(label + "\n$y$")
            if row == 1:
                ax.set_xlabel("$x$")
            selected.append({"episode": label, "seed": seed, "origin": event["origin"],
                             "first_cross_time": t0, "relative_time": rel, "raw_index": index,
                             "opposite_fish": int(crossed.sum())})
    fig.suptitle(f"Representative matched first-crossing episodes (seed {seed}); orange fish are across the divider")
    for ext in ("png", "pdf"):
        fig.savefig(root / "figures" / f"01_representative_episodes.{ext}", dpi=220)
    plt.close(fig)
    (root / "representative_episodes.json").write_text(json.dumps(selected, indent=2) + "\n")
    print(json.dumps({"seed": seed, "successful_time": float(success["first_cross_time"]),
                      "failed_time": float(failure["first_cross_time"])}, indent=2))


if __name__ == "__main__":
    main()
