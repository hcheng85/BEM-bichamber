#!/usr/bin/env python3
"""Build the synchronized (m, P_L, P_R) table from validated stride-10 caches.

The raw 3.7 GB MAT files are only audited for presence and shape.  The values
come from the earlier source-audited approach and q/J extraction caches, which
were computed synchronously from those same raw trajectories.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.io as sio


DT = 0.01
BURN_INDEX = 20_000
PAST_LAG = 10  # 0.1 time unit, strictly retrospective
SEEDS = tuple(range(20))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--approach-cache", type=Path, required=True)
    p.add_argument("--qj-cache", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns: dict[str, list[np.ndarray]] = {name: [] for name in (
        "N_L", "N_R", "N_gap", "m", "P_L", "P_R", "valid_L", "valid_R",
        "dm_past", "heading_to_gap_L", "heading_to_gap_R", "net_heading_right",
        "J_neck", "q_l050", "global_P", "vertical_center",
    )}
    audits = []
    shared_time = None
    for seed in SEEDS:
        raw = args.raw_dir / f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        approach = args.approach_cache / f"gap50_hydro_approach_seed{seed:02d}_stride10.npz"
        qj = args.qj_cache / f"gap50_hydro_seed{seed:02d}_stride10.npz"
        if not raw.is_file() or not approach.is_file() or not qj.is_file():
            raise FileNotFoundError({"raw": str(raw), "approach": str(approach), "qj": str(qj)})
        shapes = {name: shape for name, shape, _ in sio.whosmat(raw)}
        if shapes.get("state") != (3, 100, 1_000_002) or shapes.get("rdot") != (2, 100, 1_000_002):
            raise ValueError(f"seed {seed}: unexpected raw shapes {shapes}")
        with np.load(approach) as a, np.load(qj) as q:
            time_a = np.asarray(a["time"], float)
            time_q = np.asarray(q["time"], float)
            if time_a.shape != (100_001,) or not np.array_equal(time_a, time_q):
                raise ValueError(f"seed {seed}: asynchronous cache time axes")
            nl = np.asarray(a["left_n"], np.int16)
            nr = np.asarray(a["right_n"], np.int16)
            m = np.asarray(a["q"], float)
            if np.any(nl + nr != 100):
                raise ValueError(f"seed {seed}: chamber assignment is not exhaustive")
            q_error = float(np.max(np.abs(m - np.asarray(q["q"], float))))
            count_error = float(np.max(np.abs(m - (nr - nl) / 100.0)))
            if q_error > 1e-12 or count_error > 1e-12:
                raise ValueError(f"seed {seed}: occupancy mismatch {q_error=} {count_error=}")

            sl = slice(BURN_INDEX, None)
            nls, nrs, ms = nl[sl], nr[sl], m[sl]
            pl = np.asarray(a["left_pol"], float)[sl]
            pr = np.asarray(a["right_pol"], float)[sl]
            valid_l, valid_r = nls > 0, nrs > 0
            if not (np.all(np.isnan(pl[~valid_l])) and np.all(np.isnan(pr[~valid_r]))):
                raise ValueError(f"seed {seed}: empty-chamber polarization is not NaN")
            h_l = np.asarray(a["left_hx"], float)[sl]       # left -> gap is +x
            h_r = -np.asarray(a["right_hx"], float)[sl]     # right -> gap is -x
            net_right = (nls * h_l - nrs * h_r) / 100.0
            dm = np.full(ms.shape, np.nan)
            dm[PAST_LAG:] = ms[PAST_LAG:] - ms[:-PAST_LAG]
            hx = (nls * np.asarray(a["left_hx"], float)[sl] +
                  nrs * np.asarray(a["right_hx"], float)[sl]) / 100.0
            hy = (nls * np.asarray(a["left_hy"], float)[sl] +
                  nrs * np.asarray(a["right_hy"], float)[sl]) / 100.0
            yc = (nls * np.asarray(a["left_yc"], float)[sl] +
                  nrs * np.asarray(a["right_yc"], float)[sl]) / 100.0
            values = {
                "N_L": nls.astype(np.uint8), "N_R": nrs.astype(np.uint8),
                "N_gap": np.zeros_like(nls, dtype=np.uint8), "m": ms.astype(np.float32),
                "P_L": pl.astype(np.float32), "P_R": pr.astype(np.float32),
                "valid_L": valid_l, "valid_R": valid_r, "dm_past": dm.astype(np.float32),
                "heading_to_gap_L": h_l.astype(np.float32),
                "heading_to_gap_R": h_r.astype(np.float32),
                "net_heading_right": net_right.astype(np.float32),
                "J_neck": np.asarray(q["j_neck"], float)[sl].astype(np.float32),
                "q_l050": np.asarray(q["q_l050"], float)[sl].astype(np.float32),
                "global_P": np.hypot(hx, hy).astype(np.float32),
                "vertical_center": yc.astype(np.float32),
            }
            for name, value in values.items():
                columns[name].append(value)
            this_time = time_a[sl]
            if shared_time is None:
                shared_time = this_time
            elif not np.array_equal(shared_time, this_time):
                raise ValueError(f"seed {seed}: post-burn time mismatch")
            audits.append({
                "seed": seed, "raw_path": str(raw), "raw_bytes": raw.stat().st_size,
                "samples_before_burn": len(time_a), "samples_retained": len(this_time),
                "time_start": float(this_time[0]), "time_end": float(this_time[-1]),
                "max_cache_q_difference": q_error, "max_count_identity_error": count_error,
                "empty_left_frames": int(np.sum(~valid_l)), "empty_right_frames": int(np.sum(~valid_r)),
            })
        print(json.dumps(audits[-1], sort_keys=True), flush=True)

    payload = {name: np.stack(values) for name, values in columns.items()}
    payload.update(seed=np.asarray(SEEDS, np.int8), time=np.asarray(shared_time, np.float64))
    np.savez_compressed(args.output, **payload)
    manifest = {
        "output": str(args.output), "output_bytes": args.output.stat().st_size,
        "shape": [20, int(len(shared_time))], "dt": DT, "burn_time": 200.0,
        "past_increment_lag": PAST_LAG * DT,
        "chamber_assignment": "left x<0; right x>=0; exhaustive; N_gap=0",
        "polarization": "magnitude of chamber-mean unit heading, undefined/NaN when empty",
        "heading_to_gap": "left +mean(cos theta); right -mean(cos theta)",
        "net_heading_right": "(N_L*h_L - N_R*h_R)/N",
        "raw_state": "state=[x,y,heading]; rdot is total translational velocity",
        "vertical_phase": "not stored in the source cache; vertical_center retained without inferred phase",
        "audits": audits,
    }
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
