#!/usr/bin/env python3
"""Read-only preliminary stationary-drift check from committed compact data.

This is a diagnostic for choosing the next *raw-trajectory* torque replay.
It does not validate the angular-noise convention or infer a stability line.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputs", type=Path, default=HERE / "outputs")
    parser.add_argument("--formation-time", type=float, default=20.0)
    parser.add_argument(
        "--stable-fraction", type=float, default=0.98,
        help="Minimum existing macrostate occupancy to label the whole interval",
    )
    parser.add_argument(
        "--noise-sigma", type=float,
        help="Angular Brownian amplitude, defaulting provisionally to config In",
    )
    args = parser.parse_args()

    params = json.loads((HERE / "simulation_parameters.json").read_text())
    n_fish = int(params["number_of_fish"])
    sigma = float(params["In"] if args.noise_sigma is None else args.noise_sigma)
    if sigma < 0:
        parser.error("--noise-sigma must be nonnegative")

    occupancy = csv_rows(args.outputs / "state_occupancy.csv")
    snapshots = csv_rows(args.outputs / "torque_snapshot_summary.csv")
    print(f"Provisional angular noise sigma={sigma:g}; source audit required")
    print(f"Full-series averaging begins at t={args.formation_time:g}")
    print(
        "seed  state       mean_Q  stationarity_det_target  "
        "replay_det  replay_total  replay_n  endpoint_dQ_dt"
    )
    grouped: dict[str, list[tuple[float, float]]] = {"schooling": [], "milling": []}
    for row in sorted(occupancy, key=lambda r: int(r["seed"])):
        seed = int(row["seed"])
        schooling = float(row["schooling_fraction"])
        milling = float(row["milling_fraction"])
        if schooling >= args.stable_fraction:
            state = "schooling"
        elif milling >= args.stable_fraction:
            state = "milling"
        else:
            state = "mixed"

        with np.load(args.outputs / f"metrics_seed{seed}.npz") as data:
            t = np.asarray(data["time"], dtype=float)
            q = np.asarray(data["Q"], dtype=float)
        valid = t >= args.formation_time
        if valid.sum() < 2:
            raise ValueError(f"Seed {seed}: fewer than two post-formation samples")
        qsel, tsel = q[valid], t[valid]
        mean_q = float(qsel.mean())
        stationary_target = sigma**2 * (mean_q - 1.0 / n_fish)
        endpoint_rate = float((qsel[-1] - qsel[0]) / (tsel[-1] - tsel[0]))
        replay = [s for s in snapshots if int(s["seed"]) == seed and s["macrostate"] == state]
        if replay:
            replay_det = float(np.mean([float(s["bQ_deterministic"]) for s in replay]))
            replay_total = float(np.mean([float(s["bQ_total"]) for s in replay]))
        else:
            replay_det = replay_total = float("nan")
        print(
            f"{seed:4d}  {state:10s}  {mean_q:7.5f}  {stationary_target:+23.6f}  "
            f"{replay_det:+10.6f}  {replay_total:+12.6f}  {len(replay):8d}  "
            f"{endpoint_rate:+15.7g}"
        )
        if state in grouped:
            grouped[state].append((mean_q, stationary_target))

    for state, entries in grouped.items():
        if entries:
            print(
                f"Mean across {len(entries)} persistent {state} seeds: "
                f"Q={np.mean([x[0] for x in entries]):.6f}, "
                f"stationary deterministic target={np.mean([x[1] for x in entries]):+.6f}"
            )
    print(
        "NOTE: replay snapshots are sparse, biased by their selection rule, "
        "and may not match time-uniform stationary averages."
    )
    print(
        "NOTE: the endpoint rate also contains stochastic martingale "
        "fluctuations; this is not a pointwise equality test."
    )


if __name__ == "__main__":
    main()
