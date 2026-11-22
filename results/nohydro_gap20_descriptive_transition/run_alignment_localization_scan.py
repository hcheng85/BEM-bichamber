#!/usr/bin/env python3
"""Bounded alignment-strength scan with compact population output.

The archived production source remains read-only.  Its exact one-step body is
extracted and supplied explicit noise.  New runs change only Ip, initialization,
or (for the independent control) set the entire social torque to zero.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np

from analyze_approach_current_balance import fast_check_voronoi, load_gap20_production
from analyze_initiation_propagation import RawReader
from run_bottleneck_commitment_branches import parameters, raw_path
from run_nonreciprocal_alignment_intervention import extract_step


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO.parent / "BEM_bichamber_parallel_instance_no_hydro.py"
RAW_DEFAULT = Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
OUT_DEFAULT = HERE / "outputs/alignment_localization_regime"
DT = 1e-3
N = 100
SAMPLE_STEPS = 100
IP_VALUES = (0.0, 1.0, 2.0, 4.0, 8.0)
INITIALIZATIONS = ("localized", "balanced")
N_SEEDS = 3
_MOD = _STEP = _ADAPTER_TEXT = None


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def initial_state(mod, seed: int, initialization: str) -> tuple[np.ndarray, np.random.RandomState]:
    """Production distributions with a local legacy RNG for exact reproducibility."""
    rng = np.random.RandomState(seed)
    alpha = rng.uniform(-np.pi, np.pi, N)
    radius = mod.R / 1.5
    if initialization == "localized":
        rpos = radius * np.sqrt(rng.uniform(0.2, 1.0, N))
        angle = np.pi * 0.4 * rng.uniform(-1.0, 1.0, N)
        position = np.vstack((-rpos * np.cos(angle), rpos * np.sin(angle)))
    elif initialization == "balanced":
        position = np.zeros((2, N))
        # Match the built-in bothside construction: independent paired draws,
        # first 50 on the right and last 50 reflected to the left.
        for i in range(N // 2):
            rpos = radius * math.sqrt(rng.uniform(0.2, 1.0)); angle = np.pi * 0.4 * rng.uniform(-1.0, 1.0)
            position[:, i] = (rpos * math.cos(angle), rpos * math.sin(angle))
        for i in range(N // 2, N):
            rpos = radius * math.sqrt(rng.uniform(0.2, 1.0)); angle = np.pi * 0.4 * rng.uniform(-1.0, 1.0)
            position[:, i] = (-rpos * math.cos(angle), rpos * math.sin(angle))
    else:
        raise ValueError(initialization)
    return np.vstack((position, alpha)), rng


def social_rule(mod, vn: np.ndarray, state: np.ndarray, para: dict, _source_sign: int, mode: str) -> np.ndarray:
    if mode == "social_off":
        return np.zeros(N)
    if mode != "baseline":
        raise ValueError(mode)
    # Production prints "no vision" whenever an isolated fish receives the
    # already-defined zero social torque.  Suppress only that diagnostic; the
    # returned torque and every branch of GetVision are unchanged.
    with np.errstate(invalid="ignore", divide="ignore"), contextlib.redirect_stdout(io.StringIO()):
        return np.asarray(mod.GetVision(vn, state, para)).reshape(-1)


def init_worker(source: str, adapter_text: str) -> None:
    global _MOD, _STEP, _ADAPTER_TEXT
    _MOD = load_gap20_production(Path(source))
    _MOD.CheckVoronoi = lambda vn, loc, pa: fast_check_voronoi(_MOD, vn, loc)
    env = dict(vars(_MOD))
    env["GetVision_intervention"] = lambda vn, st, pa, s, mode: social_rule(_MOD, vn, st, pa, s, mode)
    exec(adapter_text, env)
    _STEP = env["archived_step"]
    _ADAPTER_TEXT = adapter_text


def task_name(task: dict) -> str:
    condition = "independent" if task["condition"] == "independent" else f"Ip{task['ip']:g}"
    return f"{condition}_{task['initialization']}_seed{task['seed']}"


def run_task(task: dict) -> dict:
    state, rng = initial_state(_MOD, task["rng_seed"], task["initialization"])
    para = parameters(_MOD); para["Ip"] = float(task["ip"])
    mode = "social_off" if task["condition"] == "independent" else "baseline"
    steps = round(task["total_time"] / DT)
    q = np.empty(steps // SAMPLE_STEPS + 1, dtype=np.int8)
    q[0] = int(2 * np.sum(state[0] >= 0) - N)
    memory = None; sample = 1; zero = np.zeros((3, N))
    for step in range(steps):
        increment = zero.copy(); increment[2] = para["In"] * math.sqrt(DT) * rng.normal(size=N)
        state, _velocity, _theta_dot, memory, _audit = _STEP(state, para, increment, memory, -1, mode)
        if (step + 1) % SAMPLE_STEPS == 0:
            q[sample] = int(2 * np.sum(state[0] >= 0) - N); sample += 1
    return dict(task=task, q=q)


def tasks(total_time: float) -> list[dict]:
    rows = []
    for ip in IP_VALUES:
        for initialization in INITIALIZATIONS:
            # Reuse the 40 archived localized Ip=4 trajectories.
            if ip == 4.0 and initialization == "localized":
                continue
            for seed in range(N_SEEDS):
                rows.append(dict(condition="alignment", ip=ip, initialization=initialization,
                                 seed=seed, rng_seed=810000 + int(ip * 1000) + seed + (100 if initialization == "balanced" else 0),
                                 total_time=total_time))
    for initialization in INITIALIZATIONS:
        for seed in range(N_SEEDS):
            rows.append(dict(condition="independent", ip=0.0, initialization=initialization,
                             seed=seed, rng_seed=910000 + seed + (100 if initialization == "balanced" else 0),
                             total_time=total_time))
    return rows


def source_audit(raw_dir: Path, mod, step, adapter_sha: str) -> dict:
    para = parameters(mod); para["Ip"] = 4.0
    state, rng = initial_state(mod, 0, "localized")
    reader = RawReader(raw_path(raw_dir, 0))
    archived0 = reader.frame("state", 0); archived1 = reader.frame("state", 1); archived_v = reader.frame("rdot", 1); reader.close()
    increment = np.zeros((3, N)); increment[2] = para["In"] * math.sqrt(DT) * rng.normal(size=N)
    predicted, velocity, _, _, _ = step(state, para, increment, None, -1, "baseline")
    balanced, _ = initial_state(mod, 1, "balanced")
    independent_test = social_rule(mod, np.ones((N, N), dtype=bool), state, para, -1, "social_off")
    return dict(
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), adapter_sha256=adapter_sha,
        localized_initial_state_max_error=float(np.max(np.abs(state - archived0))),
        first_step_state_max_error=float(np.max(np.abs(predicted - archived1))),
        first_step_velocity_max_error=float(np.max(np.abs(velocity - archived_v))),
        balanced_initial_q=int(2 * np.sum(balanced[0] >= 0) - N),
        social_off_max_abs_torque=float(np.max(np.abs(independent_test))),
        independent_definition="all visual social torque disabled, including alignment and positional term; wall, noise, speed, clipping, and collision handling retained",
    )


def simulate(output: Path, adapter_text: str, run_tasks: list[dict], workers: int) -> None:
    cache = output / "cache"; cache.mkdir(parents=True, exist_ok=True)
    todo = [task for task in run_tasks if not (cache / f"{task_name(task)}.npz").exists()]
    start = time.monotonic()
    if todo:
        with mp.get_context("spawn").Pool(min(workers, len(todo)), initializer=init_worker,
                                           initargs=(str(SOURCE), adapter_text)) as pool:
            for count, result in enumerate(pool.imap_unordered(run_task, todo), 1):
                task = result["task"]
                np.savez_compressed(cache / f"{task_name(task)}.npz", q=result["q"])
                print(f"completed {count}/{len(todo)} {task_name(task)}; elapsed {time.monotonic()-start:.1f}s", flush=True)
    manifest = [] ; qrows = []
    for task in run_tasks:
        path = cache / f"{task_name(task)}.npz"
        if not path.exists():
            raise RuntimeError(f"Missing completed run {path}")
        with np.load(path) as data:
            qrows.append(np.asarray(data["q"], dtype=np.int8))
        manifest.append(dict(run_id=task_name(task), **task))
    qarray = np.stack(qrows)
    if not np.all(np.isfinite(qarray)) or np.any(np.abs(qarray) > N) or np.any((qarray + N) % 2):
        raise RuntimeError("Invalid chamber-count trajectory")
    expected_initial = np.asarray([-N if row["initialization"] == "localized" else 0 for row in manifest])
    if not np.array_equal(qarray[:, 0], expected_initial):
        raise RuntimeError("Initialization count audit failed")
    np.savez_compressed(output / "scan_q_trajectories.npz", q_count_difference=qarray,
                        time=np.arange(qarray.shape[1]) * SAMPLE_STEPS * DT)
    write_csv(output / "scan_manifest.csv", manifest)


def self_test() -> None:
    class Fake:
        R = 15.0
    for initialization, expected in (("localized", -100), ("balanced", 0)):
        state, _ = initial_state(Fake, 7, initialization)
        assert state.shape == (3, N)
        assert int(2 * np.sum(state[0] >= 0) - N) == expected
    print("PASS: localized and exactly balanced production-distribution initializations")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=RAW_DEFAULT)
    parser.add_argument("--output", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--total-time", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=max(1, min(20, os.cpu_count() or 1)))
    parser.add_argument("--stage", choices=("audit", "simulate", "all"), default="all")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    args.output.mkdir(parents=True, exist_ok=True)
    mod = load_gap20_production(SOURCE); para = parameters(mod); step, adapter_text, adapter_sha = extract_step(mod, SOURCE)
    env = dict(vars(mod)); env["GetVision_intervention"] = lambda vn, st, pa, s, mode: social_rule(mod, vn, st, pa, s, mode)
    exec(adapter_text, env); step = env["archived_step"]
    audit = source_audit(args.raw_dir, mod, step, adapter_sha)
    (args.output / "simulator_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    if max(audit["localized_initial_state_max_error"], audit["first_step_state_max_error"],
           audit["first_step_velocity_max_error"], audit["social_off_max_abs_torque"]) > 1e-12 or audit["balanced_initial_q"] != 0:
        raise RuntimeError("Simulator or initialization audit failed")
    print(json.dumps(audit, indent=2))
    if args.stage == "audit":
        return
    run_tasks = tasks(args.total_time)
    simulate(args.output, adapter_text, run_tasks, args.workers)
    metadata = dict(total_time=args.total_time, native_dt=DT, q_sample_dt=SAMPLE_STEPS * DT,
                    n_new_runs=len(run_tasks), seeds_per_initialization=N_SEEDS,
                    parameter_values=list(IP_VALUES), initializations=list(INITIALIZATIONS),
                    archived_ip4_localized_runs=40, production_source_unchanged=True,
                    fixed_parameters=dict(N=N, In=float(para["In"]), Iw=float(para["Iw"]), speed=1.0,
                                          radius_ratio=int(para["radius_ratio"]), gap_ratio=int(para["gap_ratio"]),
                                          dt=float(para["dt"]), positional_social_rule="unchanged production visual-position term"),
                    initialization_definitions=dict(localized="production distribution with all 100 fish in left chamber",
                                                    balanced="same chamber distribution with exactly 50 fish on each side; random headings"),
                    output="compact q count-difference trajectories only", **audit)
    (args.output / "simulation_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
