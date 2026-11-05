#!/usr/bin/env python3
"""Chunked first-stage reduction of the ten circular full-BEM trajectories."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import numpy as np

from matv5_memmap import open_state_rdot


FILE_PATTERN = "If1e-02_Ip40_In05_R50_num100_wid100_n*_BEM.mat"


def seed_from_name(path: Path) -> int:
    match = re.search(r"_n(\d+)_BEM\.mat$", path.name)
    if not match:
        raise ValueError(path.name)
    return int(match.group(1))


def wrap(angle: np.ndarray) -> np.ndarray:
    return np.arctan2(np.sin(angle), np.cos(angle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def process_file(
    path: Path,
    outdir: Path,
    radius: float,
    dt: float,
    metric_stride: int,
    distribution_stride: int,
    chunk_frames: int,
) -> dict:
    seed = seed_from_name(path)
    state, rdot, info = open_state_rdot(path)
    if state.shape[0] != 3 or rdot.shape[0] != 2 or state.shape[1:] != rdot.shape[1:]:
        raise ValueError(f"Unexpected arrays: state={state.shape}, rdot={rdot.shape}")
    nfish, nt = state.shape[1], state.shape[2]
    metric_indices = np.arange(0, nt, metric_stride, dtype=np.int64)
    nm = len(metric_indices)
    fields = {
        name: np.empty(nm, dtype=np.float32)
        for name in (
            "P", "Q", "Pi", "M_heading", "M_velocity", "radial_mean",
            "radial_std", "com_radius", "mean_speed",
        )
    }
    radial_edges = np.linspace(0, 1.05, 107)
    heading_edges = np.linspace(-np.pi, np.pi, 145)
    radial_hist = np.zeros(len(radial_edges)-1, dtype=np.int64)
    heading_hist = np.zeros(len(heading_edges)-1, dtype=np.int64)
    safety = {
        "ray_layer_count": 0, "sticky_zero_velocity_count": 0,
        "outside_count": 0, "min_ray_distance": np.inf,
        "max_radius": 0.0, "first_ray_layer_step": -1,
        "first_sticky_step": -1,
    }

    for start in range(0, nt, chunk_frames):
        stop = min(nt, start + chunk_frames)
        # Copies remain bounded: 20 MB total at the default 5,000-frame chunk.
        s = np.array(state[:, :, start:stop], copy=True)
        v = np.array(rdot[:, :, start:stop], copy=True)
        x, y, theta = s[0], s[1], s[2]
        rr = np.hypot(x, y)
        speed = np.hypot(v[0], v[1])
        projection = x*np.cos(theta) + y*np.sin(theta)
        discriminant = np.maximum(projection**2 + radius**2 - rr**2, 0.0)
        ray_distance = -projection + np.sqrt(discriminant)
        ray_mask = ray_distance <= 0.1
        sticky_mask = (radius - rr <= 0.1000001) & (speed <= 1e-14)
        outside_mask = rr > radius + 1e-8
        safety["ray_layer_count"] += int(np.count_nonzero(ray_mask))
        safety["sticky_zero_velocity_count"] += int(np.count_nonzero(sticky_mask))
        safety["outside_count"] += int(np.count_nonzero(outside_mask))
        safety["min_ray_distance"] = min(safety["min_ray_distance"], float(np.min(ray_distance)))
        safety["max_radius"] = max(safety["max_radius"], float(np.max(rr)))
        if safety["first_ray_layer_step"] < 0 and np.any(ray_mask):
            safety["first_ray_layer_step"] = start + int(np.nonzero(np.any(ray_mask, axis=0))[0][0])
        if safety["first_sticky_step"] < 0 and np.any(sticky_mask):
            safety["first_sticky_step"] = start + int(np.nonzero(np.any(sticky_mask, axis=0))[0][0])

        # Aggregate distributions at a coarser, documented stride.
        local_dist = np.nonzero(np.arange(start, stop) % distribution_stride == 0)[0]
        if len(local_dist):
            alpha = wrap(theta[:, local_dist] - np.arctan2(y[:, local_dist], x[:, local_dist]))
            radial_hist += np.histogram((rr[:, local_dist]/radius).ravel(), radial_edges)[0]
            heading_hist += np.histogram(alpha.ravel(), heading_edges)[0]

        local = np.nonzero(np.arange(start, stop) % metric_stride == 0)[0]
        if not len(local):
            continue
        global_idx = np.arange(start, stop)[local]
        out_idx = global_idx // metric_stride
        xs, ys, ths = x[:, local], y[:, local], theta[:, local]
        # Index the two-dimensional component views so fish remains the first
        # axis (mixed basic/advanced indexing otherwise moves time in front).
        vxs, vys = v[0][:, local], v[1][:, local]
        c, sn = np.cos(ths), np.sin(ths)
        px, py = np.mean(c, axis=0), np.mean(sn, axis=0)
        P = np.hypot(px, py)
        rsel = np.hypot(xs, ys)
        fields["P"][out_idx] = P
        fields["Q"][out_idx] = P**2
        fields["radial_mean"][out_idx] = np.mean(rsel/radius, axis=0)
        fields["radial_std"][out_idx] = np.std(rsel/radius, axis=0)
        comx, comy = np.mean(xs, axis=0), np.mean(ys, axis=0)
        fields["com_radius"][out_idx] = np.hypot(comx, comy)/radius
        sp = np.hypot(vxs, vys)
        fields["mean_speed"][out_idx] = np.nanmean(sp, axis=0)

        # Original simulator definition: mean absolute normalized angular
        # momentum around the instantaneous school center, using actual rdot.
        rx, ry = xs-comx, ys-comy
        rn = np.hypot(rx, ry)
        with np.errstate(divide="ignore", invalid="ignore"):
            local_cross = (rx*vys - ry*vxs)/(rn*sp)
        fields["Pi"][out_idx] = np.nanmean(np.abs(local_cross), axis=0)

        # Signed milling about the fixed circle center. M_heading is the exact
        # variable used by the analytical moment equation; M_velocity is saved
        # as a diagnostic of hydrodynamic translation.
        with np.errstate(divide="ignore", invalid="ignore"):
            fields["M_heading"][out_idx] = np.nanmean((xs*sn-ys*c)/rsel, axis=0)
            fields["M_velocity"][out_idx] = np.nanmean((xs*vys-ys*vxs)/(rsel*sp), axis=0)

    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        outdir/f"metrics_seed{seed}.npz",
        step=metric_indices, time=metric_indices*dt,
        radial_edges=radial_edges, radial_hist=radial_hist,
        heading_edges=heading_edges, heading_hist=heading_hist,
        **fields,
    )
    row = {
        "seed": seed, "file": str(path), "file_bytes": path.stat().st_size,
        "state_shape": "x".join(map(str, state.shape)),
        "rdot_shape": "x".join(map(str, rdot.shape)),
        "state_data_offset": info["state"].data_offset,
        "rdot_data_offset": info["rdot"].data_offset,
        "metric_stride": metric_stride, "metric_dt": metric_stride*dt,
        "distribution_stride": distribution_stride,
        "n_metric_samples": nm, **safety,
    }
    (outdir/f"stream_audit_seed{seed}.json").write_text(json.dumps(row, indent=2)+"\n")
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent/"outputs")
    parser.add_argument("--radius", type=float, default=7.158910531638177)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--metric-stride", type=int, default=10)
    parser.add_argument("--distribution-stride", type=int, default=100)
    parser.add_argument("--chunk-frames", type=int, default=5000)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    files = sorted(args.data_dir.glob(FILE_PATTERN), key=seed_from_name)
    if args.limit is not None: files = files[:args.limit]
    if not files: raise FileNotFoundError(args.data_dir/FILE_PATTERN)
    rows = []
    for i, path in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {path.name}", flush=True)
        rows.append(process_file(path, args.outdir, args.radius, args.dt,
                                 args.metric_stride, args.distribution_stride,
                                 args.chunk_frames))
    write_csv(args.outdir/"streaming_inventory.csv", rows)
    print(json.dumps({"files":len(rows), "outdir":str(args.outdir)}, indent=2))


if __name__ == "__main__":
    main()
