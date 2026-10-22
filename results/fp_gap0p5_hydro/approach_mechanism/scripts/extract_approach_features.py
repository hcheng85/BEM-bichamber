#!/usr/bin/env python3
"""Extract compact school-level observables from the 20 gap-50 BEM trajectories.

The raw MATLAB files are read one realization at a time and are never changed.
Only stride-sampled aggregate quantities are saved; individual-fish trajectories
and raw fields are deliberately excluded from the output cache.
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import scipy.io as sio


def scalar(value: object) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def side_features(x: np.ndarray, y: np.ndarray, theta: np.ndarray, mask: np.ndarray) -> dict[str, np.ndarray]:
    n = mask.sum(axis=1)
    denom = np.maximum(n, 1)
    xf = np.where(mask, x, 0.0)
    yf = np.where(mask, y, 0.0)
    c = np.where(mask, np.cos(theta), 0.0)
    s = np.where(mask, np.sin(theta), 0.0)
    xc = xf.sum(axis=1) / denom
    yc = yf.sum(axis=1) / denom
    cx = c.sum(axis=1) / denom
    sy = s.sum(axis=1) / denom
    xvar = (np.where(mask, (x - xc[:, None]) ** 2, 0.0).sum(axis=1) / denom)
    yvar = (np.where(mask, (y - yc[:, None]) ** 2, 0.0).sum(axis=1) / denom)
    xmin = np.where(mask, x, np.inf).min(axis=1)
    xmax = np.where(mask, x, -np.inf).max(axis=1)
    ymin = np.where(mask, y, np.inf).min(axis=1)
    ymax = np.where(mask, y, -np.inf).max(axis=1)
    empty = n == 0
    arrays = {
        "n": n.astype(np.int16), "xc": xc, "yc": yc,
        "hx": cx, "hy": sy, "pol": np.hypot(cx, sy),
        "xstd": np.sqrt(xvar), "ystd": np.sqrt(yvar),
        "xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
    }
    for key, value in arrays.items():
        if key != "n":
            value[empty] = np.nan
    return arrays


def extract_seed(raw_dir: Path, out_dir: Path, seed: int, stride: int, overwrite: bool) -> dict[str, object]:
    raw = raw_dir / f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
    dest = out_dir / f"gap50_hydro_approach_seed{seed:02d}_stride{stride}.npz"
    if dest.exists() and not overwrite:
        with np.load(dest) as cached:
            return {"seed": seed, "status": "cached", "samples": int(cached["time"].size)}
    if not raw.is_file():
        raise FileNotFoundError(raw)
    shape = {name: shp for name, shp, _dtype in sio.whosmat(raw)}.get("state")
    if shape != (3, 100, 1_000_002):
        raise ValueError(f"{raw}: unexpected state shape {shape}")
    loaded = sio.loadmat(raw, variable_names=["state", "para"], squeeze_me=True, struct_as_record=False)
    state = np.asarray(loaded["state"])
    para = loaded["para"]
    metadata = {name: scalar(getattr(para, name)) for name in
                ("Ip", "In", "R", "dt", "Iw", "If", "num", "total_time", "gap", "radius_ratio", "gap_ratio")}
    required = {"Ip": 4, "In": 0.5, "dt": 0.001, "Iw": 0.94, "If": 0.01,
                "num": 100, "radius_ratio": 50, "gap_ratio": 50}
    for name, expected in required.items():
        if not np.isclose(metadata[name], expected):
            raise ValueError(f"seed {seed}: {name}={metadata[name]} != {expected}")
    indices = np.arange(0, state.shape[2], stride, dtype=np.int64)
    result: dict[str, np.ndarray] = {"time": indices * metadata["dt"]}
    fields = ("n", "xc", "yc", "hx", "hy", "pol", "xstd", "ystd", "xmin", "xmax", "ymin", "ymax")
    for side in ("left", "right"):
        for field in fields:
            dtype = np.int16 if field == "n" else np.float64
            result[f"{side}_{field}"] = np.empty(indices.size, dtype=dtype)
    chunk = 10_000
    for start in range(0, indices.size, chunk):
        stop = min(start + chunk, indices.size)
        idx = indices[start:stop]
        x = state[0, :, idx]
        y = state[1, :, idx]
        theta = state[2, :, idx]
        for side, mask in (("left", x < 0), ("right", x >= 0)):
            values = side_features(x, y, theta, mask)
            for field in fields:
                result[f"{side}_{field}"][start:stop] = values[field]
    result["q"] = (result["right_n"].astype(float) - result["left_n"].astype(float)) / metadata["num"]
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        dest, seed=np.int64(seed), stride=np.int64(stride),
        metadata_json=np.array(json.dumps(metadata, sort_keys=True)), **result,
    )
    del loaded, state
    gc.collect()
    return {"seed": seed, "status": "extracted", "samples": int(indices.size), "path": str(dest), **metadata}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--seeds", type=int, nargs="*", default=list(range(20)))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.stride < 1:
        raise SystemExit("stride must be positive")
    raw_dir = args.raw_dir.expanduser().resolve()
    out_dir = args.output_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for seed in args.seeds:
        row = extract_seed(raw_dir, out_dir, seed, args.stride, args.overwrite)
        manifest.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    (out_dir / "extraction_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
