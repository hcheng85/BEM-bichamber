#!/usr/bin/env python3
"""Extract compact continuous observables from the gap-50 full-BEM trajectories.

The raw MAT files are about 4 GB each.  This script reads one realization at a
time and stores only quantities needed for the Fokker--Planck closure tests.
It never modifies the raw or processed simulation files.
"""
from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import scipy.io as sio


RAW_DIR = Path("BEM_data")
PROCESSED_DIR = Path("Processed_data")
OUT_DIR = Path("features")
ELLS = (0.25, 0.5, 1.0)


def scalar(value: object) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def raw_path(seed: int) -> Path:
    return RAW_DIR / f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"


def processed_path(seed: int) -> Path:
    return PROCESSED_DIR / f"event_histogram_If1e-02_In05_num100_R50_gap50_n{seed}_BEM.mat"


def extract_seed(seed: int, stride: int, overwrite: bool) -> dict[str, object]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    destination = OUT_DIR / f"gap50_hydro_seed{seed:02d}_stride{stride}.npz"
    if destination.exists() and not overwrite:
        with np.load(destination) as cached:
            return {
                "seed": seed,
                "status": "cached",
                "samples": int(cached["time"].size),
                "path": str(destination),
            }

    raw = raw_path(seed)
    processed = processed_path(seed)
    if not raw.is_file() or not processed.is_file():
        raise FileNotFoundError(f"missing seed {seed}: raw={raw.exists()} processed={processed.exists()}")

    variables = {name: shape for name, shape, _dtype in sio.whosmat(raw)}
    expected = {"state": (3, 100, 1_000_002), "rdot": (2, 100, 1_000_002)}
    for name, shape in expected.items():
        if variables.get(name) != shape:
            raise ValueError(f"{raw}: {name} shape {variables.get(name)} != {shape}")

    data = sio.loadmat(raw, variable_names=["state", "rdot", "para"], squeeze_me=True, struct_as_record=False)
    state = np.asarray(data["state"])
    rdot = np.asarray(data["rdot"])
    para = data["para"]
    metadata = {name: scalar(getattr(para, name)) for name in ("Ip", "In", "R", "dt", "Iw", "If", "num", "total_time", "gap", "radius_ratio", "gap_ratio")}
    required = {"Ip": 4.0, "In": 0.5, "dt": 0.001, "Iw": 0.94, "If": 0.01, "num": 100.0, "total_time": 1000.0, "radius_ratio": 50.0, "gap_ratio": 50.0}
    for name, expected_value in required.items():
        if not np.isclose(metadata[name], expected_value):
            raise ValueError(f"seed {seed}: {name}={metadata[name]} != {expected_value}")

    indices = np.arange(0, state.shape[2], stride, dtype=np.int64)
    nsample = indices.size
    q_smooth = {ell: np.empty(nsample, dtype=np.float64) for ell in ELLS}
    # Advective time derivative of q_ell for ell=0.5.  Since translational
    # noise is absent, this is a physically defined neck-current observable.
    j_neck = np.empty(nsample, dtype=np.float64)
    neck_count = np.empty(nsample, dtype=np.int16)
    chunk = 10_000
    ell_current = 0.5
    gap_half = metadata["gap"] / 2.0
    for start in range(0, nsample, chunk):
        stop = min(start + chunk, nsample)
        idx = indices[start:stop]
        # NumPy moves the advanced time index to the leading axis here, so the
        # resulting arrays already have shape (sample, fish).
        x = state[0, :, idx]
        y = state[1, :, idx]
        vx = rdot[0, :, idx]
        for ell in ELLS:
            q_smooth[ell][start:stop] = np.mean(np.tanh(x / ell), axis=1)
        tanh_x = np.tanh(x / ell_current)
        j_neck[start:stop] = np.mean(vx * (1.0 - tanh_x * tanh_x) / ell_current, axis=1)
        neck_count[start:stop] = np.sum((np.abs(x) <= ell_current) & (np.abs(y) <= gap_half), axis=1)

    proc = sio.loadmat(processed, variable_names=["num_fish_left"], squeeze_me=True)
    nleft = np.asarray(proc["num_fish_left"], dtype=np.int16).reshape(-1)
    if nleft.size != state.shape[2]:
        raise ValueError(f"seed {seed}: processed length {nleft.size} != raw length {state.shape[2]}")
    q = 1.0 - 2.0 * nleft[indices].astype(np.float64) / metadata["num"]
    time = indices.astype(np.float64) * metadata["dt"]

    np.savez_compressed(
        destination,
        seed=np.int64(seed),
        stride=np.int64(stride),
        time=time,
        q=q,
        q_l025=q_smooth[0.25],
        q_l050=q_smooth[0.5],
        q_l100=q_smooth[1.0],
        j_neck=j_neck,
        neck_count=neck_count,
        metadata_json=np.array(json.dumps(metadata, sort_keys=True)),
    )
    del data, state, rdot, proc
    gc.collect()
    return {"seed": seed, "status": "extracted", "samples": nsample, "path": str(destination), **metadata}


def main() -> None:
    global RAW_DIR, PROCESSED_DIR, OUT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True, help="directory containing the full BEM MAT trajectories")
    parser.add_argument("--processed-dir", type=Path, required=True, help="directory containing processed event MAT files")
    parser.add_argument("--output-dir", type=Path, required=True, help="destination for compact NPZ feature caches")
    parser.add_argument("--stride", type=int, default=10, help="raw steps between cached samples")
    parser.add_argument("--seeds", type=int, nargs="*", default=list(range(20)))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    RAW_DIR = args.raw_dir.expanduser().resolve()
    PROCESSED_DIR = args.processed_dir.expanduser().resolve()
    OUT_DIR = args.output_dir.expanduser().resolve()
    if args.stride < 1:
        raise SystemExit("stride must be positive")
    summaries = []
    for seed in args.seeds:
        summary = extract_seed(seed, args.stride, args.overwrite)
        summaries.append(summary)
        print(json.dumps(summary, sort_keys=True), flush=True)
    (OUT_DIR / "extraction_manifest.json").write_text(json.dumps(summaries, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
