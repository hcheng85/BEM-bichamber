#!/usr/bin/env python3
"""Audit and run the bounded gap-20 microscopic restart pilot.

The production source and archived trajectories are read-only.  The native
timestep body is extracted with AST and changed only to accept an explicit
angular-noise increment and an optional cross-divider GetVision filter.
"""
from __future__ import annotations

import argparse
import ast
import contextlib
import copy
import csv
import hashlib
import io
import json
import multiprocessing as mp
import os
import sys
import time
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
from analyze_approach_current_balance import RawReader, fast_check_voronoi, load_gap20_production  # noqa: E402

N = 100
DT = 0.001
HORIZON_STEPS = 5000
N_CLONES = 8
N_SAMPLES = 14
EPSILON = 0.15
SOURCE = REPO.parent / "BEM_bichamber_parallel_instance_no_hydro.py"
RAW_DEFAULT = Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
OUT_DEFAULT = HERE / "outputs" / "bottleneck_commitment_branches"
RISK = HERE / "outputs" / "q_gap_polarization" / "prospective_q_a_windows.csv"
EVENTS = HERE / "outputs" / "event_aligned_continuation" / "event_aligned_samples.csv"

_MOD = None
_PARA = None
_STEP = None


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def parameters(mod) -> dict:
    return dict(Ip=4.0, In=0.5, R=mod.R, light=1.0, dt=DT, Iw=0.94, If=0.0,
                num=N, delta=0.01, total_time=1000, root="", n=-1,
                gap=mod.gap, radius_ratio=50, gap_ratio=20)


def raw_path(root: Path, seed: int) -> Path:
    return root / f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"


def social_with_filter(mod, voro: np.ndarray, state: np.ndarray, para: dict,
                       crossgap_off: bool) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"), contextlib.redirect_stdout(io.StringIO()):
        full = np.asarray(mod.GetVision(voro, state, para)).reshape(-1)
        if not crossgap_off:
            return full
        opposite = state[0, :, None] * state[0, None, :] < 0
        within = voro.copy(); within[opposite] = False
        within_renormalized = np.asarray(mod.GetVision(within, state, para)).reshape(-1)
    # GetVision is a visual-weighted numerator divided by sum_visual.  Restore
    # the original denominator so only opposite-chamber source terms vanish;
    # within-chamber terms retain their exact original weights.
    full_denominator = vision_denominator(mod, voro, state, para)
    within_denominator = vision_denominator(mod, within, state, para)
    return np.divide(within_renormalized * within_denominator, full_denominator,
                     out=np.zeros(N), where=full_denominator != 0)


def vision_denominator(mod, vn: np.ndarray, state: np.ndarray, para: dict) -> np.ndarray:
    """Reproduce GetVision's exact per-fish sum_visual without changing terms."""
    list_i = np.tile(np.arange(para["num"]), para["num"])
    neighbors = list_i[np.ndarray.flatten(vn)]
    counts = np.sum(vn, axis=0); maximum = int(np.amax(counts))
    neighbor_i = np.zeros((para["num"], maximum + 1))
    neighbor_i[np.arange(para["num"]), counts] = para["num"]
    neighbor_i = np.cumsum(neighbor_i[:, :-1], axis=1).astype(int)
    neighbor_i[neighbor_i == 0] = neighbors
    x = np.append(state[0], np.nan); y = np.append(state[1], np.nan)
    dx = x[neighbor_i] - state[0, :, None]; dy = y[neighbor_i] - state[1, :, None]
    theta = mod.getangle(state[2, :, None], dx, dy, np.sqrt(dx * dx + dy * dy))
    return np.nansum(1 + np.cos(theta), axis=1)


def extract_step(mod, source: Path):
    """Extract production statements, replacing only noise input and social filter."""
    tree = ast.parse(source.read_text())
    fn = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "step_Computation")
    loop = next(node for node in fn.body if isinstance(node, ast.For))
    first = next(i for i, node in enumerate(loop.body)
                 if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                 and node.targets[0].id == "CurState")
    body = copy.deepcopy(loop.body[first:])
    noise_replacements = social_replacements = 0
    for node in body:
        if not isinstance(node, ast.Assign) or not isinstance(node.targets[0], ast.Name):
            continue
        name = node.targets[0].id
        if name == "wNoise":
            node.value = ast.Name(id="noise_increment", ctx=ast.Load()); noise_replacements += 1
        elif name == "wVision":
            node.value = ast.Call(func=ast.Name(id="GetVision_intervention", ctx=ast.Load()),
                                  args=[ast.Name(id="voro", ctx=ast.Load()),
                                        ast.Name(id="CurState", ctx=ast.Load()),
                                        ast.Name(id="para", ctx=ast.Load()),
                                        ast.Name(id="crossgap_off", ctx=ast.Load())], keywords=[])
            social_replacements += 1
    if (noise_replacements, social_replacements) != (1, 1):
        raise RuntimeError("Could not identify unique production noise/social assignments")
    prefix = ast.parse('''def archived_step(current, para, noise_increment, dis_last, crossgap_off=False):
    state = np.empty((3, para["num"], 2))
    state[:,:,0] = current
    rDotList = np.empty((2, para["num"], 2))
    count = 0
''').body[0]
    suffix = ast.parse('''return state[:,:,count].copy(), rDot.copy(), thetaDot.copy(), dis_last.copy(), dict(collision=int(np.sum(np.linalg.norm(rDot,axis=0)==0)), safety_candidates=int(np.sum(dis>=safedis)))''').body
    prefix.body += body + suffix
    module = ast.fix_missing_locations(ast.Module(body=[prefix], type_ignores=[]))
    text = ast.unparse(module) + "\n"
    env = dict(vars(mod))
    env["CheckVoronoi"] = lambda vn, loc, pa: fast_check_voronoi(mod, vn, loc)
    env["GetVision_intervention"] = lambda vn, st, pa, off: social_with_filter(mod, vn, st, pa, off)
    exec(compile(module, str(source) + "::restart_adapter", "exec"), env)
    return env["archived_step"], text, hashlib.sha256(text.encode()).hexdigest()


def majority_polarization(state: np.ndarray, source_sign: int) -> float:
    mask = state[0] >= 0 if source_sign == 1 else state[0] < 0
    com = state[:2, mask].mean(axis=1); gap_hat = -com / np.linalg.norm(com)
    heading = np.vstack((np.cos(state[2, mask]), np.sin(state[2, mask]))).mean(axis=1)
    return float(heading @ gap_hat)


def nudge(state: np.ndarray, source_sign: int, epsilon: float) -> np.ndarray:
    out = state.copy(); mask = state[0] >= 0 if source_sign == 1 else state[0] < 0
    com = state[:2, mask].mean(axis=1); phi = np.arctan2(-com[1], -com[0])
    out[2, mask] += epsilon * np.sin(phi - out[2, mask])
    return out


def choose_pre() -> list[dict]:
    rows = [r for r in read_csv(RISK) if abs(float(r["q"])) > 0 and int(r["native_step"]) > 0]
    selected = []
    for seed in sorted({int(r["seed"]) for r in rows})[:N_SAMPLES]:
        candidates = [r for r in rows if int(r["seed"]) == seed]
        target_a = (0.35, 0.60, 0.85)[seed % 3]
        chosen = min(candidates, key=lambda r: (abs(float(r["a"]) - target_a)
                                                + 0.35 * abs(abs(float(r["q"])) - 0.8),
                                                int(r["native_step"])))
        selected.append(dict(phase="pre", sample_id=f"pre_s{seed}_k{chosen['native_step']}",
                             seed=seed, native_step=int(chosen["native_step"]), q=float(chosen["q"]),
                             a_table=float(chosen["a"]), selection_target_a=target_a,
                             selection_rule="nearest current-only target in (abs(q),a); future label unused"))
    return selected


def choose_post() -> list[dict]:
    rows = [r for r in read_csv(EVENTS) if not r.get("exclude", "")
            and int(r["first_majority_to_minority"]) == 1 and int(r["prelead_endpoints_match"]) == 1]
    selected = []
    for seed in sorted({int(r["seed"]) for r in rows})[:N_SAMPLES]:
        chosen = min((r for r in rows if int(r["seed"]) == seed), key=lambda r: int(r["first_native"]))
        selected.append(dict(phase="post", sample_id=f"post_s{seed}_e{chosen['event_id']}",
                             seed=seed, event_id=int(chosen["event_id"]),
                             native_step=int(chosen["first_native"]), q_pre=float(chosen["q_pre"]),
                             selection_rule="first eligible majority-to-minority event per seed; continuation unused"))
    return selected


def load_origin(root: Path, item: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    reader = RawReader(raw_path(root, int(item["seed"])))
    k = int(item["native_step"])
    state = np.asarray(reader.block("state", k, k)[:, :, 0], float)
    previous = np.asarray(reader.block("state", k - 1, k - 1)[:, :, 0], float)
    target = None
    if k + 1 <= 1_000_001:
        target = np.asarray(reader.block("state", k + 1, k + 1)[:, :, 0], float)
    reader.close()
    return state, previous, target


def setup_origins(root: Path, step, para: dict, pre: list[dict], post: list[dict]) -> list[dict]:
    zero = np.zeros((3, N)); origins = []
    for item in pre + post:
        state, previous, _ = load_origin(root, item)
        _, _, _, dis_last, _ = step(previous, para, zero, None, False)
        source_sign = int(np.sign(item["q"] if item["phase"] == "pre" else item["q_pre"]))
        a0 = majority_polarization(state, source_sign)
        record = dict(item, source_sign=source_sign, nright0=int(np.sum(state[0] >= 0)),
                      a_baseline=a0, state=state, dis_last=dis_last)
        if item["phase"] == "pre":
            toward = nudge(state, source_sign, EPSILON); away = nudge(state, source_sign, -EPSILON)
            record["a_toward"] = majority_polarization(toward, source_sign)
            record["a_away"] = majority_polarization(away, source_sign)
            if not record["a_toward"] >= a0 >= record["a_away"]:
                raise RuntimeError(f"Polarization monotonicity failed for {item['sample_id']}")
        origins.append(record)
    return origins


def audit(root: Path, mod, step, para: dict, origins: list[dict], source_sha: str, adapter_sha: str) -> dict:
    zero = np.zeros((3, N)); parity = []; coupling = []; fast = []
    original_check = mod.CheckVoronoi
    for item in origins[:20]:
        state, _, target = load_origin(root, item); k = int(item["native_step"])
        vn0 = mod.GetVoronoi(state[:2], para); slow = original_check(vn0.copy(), state[:2], para)
        quick = fast_check_voronoi(mod, vn0.copy(), state[:2])
        fast.append(dict(sample_id=item["sample_id"], equal=bool(np.array_equal(slow, quick))))
        full = np.asarray(social_with_filter(mod, slow, state, para, False)).reshape(-1)
        filtered = np.asarray(social_with_filter(mod, slow, state, para, True)).reshape(-1)
        cross = (state[0, :, None] * state[0, None, :] < 0) & slow
        coupling.append(dict(sample_id=item["sample_id"], phase=item["phase"],
                             admissible_crossgap_pairs=int(np.sum(np.triu(cross, 1))),
                             mean_abs_torque_change=float(np.mean(np.abs(full - filtered))),
                             max_abs_torque_change=float(np.max(np.abs(full - filtered)))))
        deterministic, _, _, _, _ = step(state, para, zero, item["dis_last"], False)
        increment = zero.copy(); increment[2] = target[2] - deterministic[2]
        predicted, velocity, _, _, flags = step(state, para, increment, item["dis_last"], False)
        reader = RawReader(raw_path(root, int(item["seed"])))
        archived_velocity = np.asarray(reader.block("rdot", k + 1, k + 1)[:, :, 0], float); reader.close()
        parity.append(dict(sample_id=item["sample_id"], seed=int(item["seed"]), native_step=k,
                           max_position_error=float(np.max(np.abs(predicted[:2] - target[:2]))),
                           max_heading_error=float(np.max(np.abs(predicted[2] - target[2]))),
                           max_velocity_error=float(np.max(np.abs(velocity - archived_velocity))),
                           inferred_noise_mean=float(np.mean(increment[2])),
                           inferred_noise_std=float(np.std(increment[2], ddof=1)), **flags))
    max_state = max(max(r["max_position_error"], r["max_heading_error"]) for r in parity)
    max_velocity = max(r["max_velocity_error"] for r in parity)
    cross_nonzero = any(r["max_abs_torque_change"] > 1e-14 for r in coupling)
    return dict(source_sha256=source_sha, adapter_sha256=adapter_sha,
                source_public_restart_support=False, source_public_exogenous_noise_support=False,
                geometry=dict(radius_ratio=50, gap_ratio=20, R=float(mod.R), gap=float(mod.gap)),
                parameters=dict(dt=DT, N=N, Ip=4.0, In=0.5, Iw=0.94, If=0.0,
                                wall_safety_distance=-0.1, wall_high_order=6),
                crossgap_intervention="remove only opposite-chamber GetVision numerator terms; retain original full-neighborhood visual denominator",
                parity=dict(passed=bool(max_state < 1e-12 and max_velocity < 1e-12),
                            n_steps=len(parity), max_state_error=max_state,
                            max_velocity_error=max_velocity, rows=parity),
                optimized_voronoi=dict(passed=all(r["equal"] for r in fast), n_frames=len(fast), rows=fast),
                crossgap_social=dict(nonzero=cross_nonzero, n_frames=len(coupling), rows=coupling,
                                     total_admissible_pairs=sum(r["admissible_crossgap_pairs"] for r in coupling),
                                     mean_abs_torque_change=float(np.mean([r["mean_abs_torque_change"] for r in coupling])),
                                     max_abs_torque_change=max(r["max_abs_torque_change"] for r in coupling)))


def init_worker(source: str, adapter_text: str):
    global _MOD, _PARA, _STEP
    _MOD = load_gap20_production(Path(source)); _MOD.CheckVoronoi = lambda vn, loc, para: fast_check_voronoi(_MOD, vn, loc)
    _PARA = parameters(_MOD)
    env = dict(vars(_MOD)); env["GetVision_intervention"] = lambda vn, st, pa, off: social_with_filter(_MOD, vn, st, pa, off)
    exec(adapter_text, env); _STEP = env["archived_step"]


def run_task(task: dict) -> list[dict]:
    phase = task["phase"]; clone = int(task["clone"]); state0 = task["state"]
    rng = np.random.default_rng(int(task["noise_seed"]))
    noises = 0.5 * np.sqrt(DT) * rng.normal(size=(HORIZON_STEPS, N))
    if phase == "pre":
        states = dict(baseline=state0.copy(), toward=nudge(state0, task["source_sign"], EPSILON),
                      away=nudge(state0, task["source_sign"], -EPSILON)); interventions = {k: False for k in states}
    else:
        states = {"full": state0.copy()}
        if task["crossgap_nonzero"]: states["crossgap_off"] = state0.copy()
        interventions = {"full": False, "crossgap_off": True}
    memory = {arm: task["dis_last"].copy() for arm in states}; entered = {arm: False for arm in states}
    for n in range(HORIZON_STEPS):
        increment = np.zeros((3, N)); increment[2] = noises[n]
        for arm in states:
            previous = states[arm]
            nxt, _, _, memory[arm], _ = _STEP(previous, _PARA, increment, memory[arm], interventions[arm])
            if phase == "pre" and np.any((task["source_sign"] * previous[0] >= 0)
                                          & (task["source_sign"] * nxt[0] < 0)):
                entered[arm] = True
            states[arm] = nxt
    rows = []
    for arm, state in states.items():
        rows.append(dict(phase=phase, sample_id=task["sample_id"], seed=task["seed"], clone=clone,
                         arm=arm, source_sign=task["source_sign"], nright0=task["nright0"],
                         nrightT=int(np.sum(state[0] >= 0)), first_entry=int(entered[arm]) if phase == "pre" else 1,
                         a_initial=majority_polarization((state0 if arm in ("baseline", "full", "crossgap_off")
                                                        else nudge(state0, task["source_sign"],
                                                                   EPSILON if arm == "toward" else -EPSILON)),
                                                       task["source_sign"]),
                         a_baseline=task["a_baseline"]))
    return rows


def simulate(origins: list[dict], audit_data: dict, adapter_text: str, output: Path, workers: int) -> list[dict]:
    tasks = []
    for item in origins:
        for clone in range(N_CLONES):
            phase_code = 1 if item["phase"] == "pre" else 2
            tasks.append(dict(item, clone=clone, noise_seed=202611130000 + phase_code * 100000
                              + int(item["seed"]) * 100 + clone,
                              crossgap_nonzero=audit_data["crossgap_social"]["nonzero"]))
    partial = output / "branches_partial.csv"; existing = read_csv(partial) if partial.exists() else []
    candidate_keys = {(r["phase"], r["sample_id"], int(r["clone"])) for r in existing}
    expected_arms = lambda t: 3 if t["phase"] == "pre" else (2 if t["crossgap_nonzero"] else 1)
    counts = {key: sum((r["phase"], r["sample_id"], int(r["clone"])) == key for r in existing)
              for key in candidate_keys}
    task_by_key = {(t["phase"], t["sample_id"], t["clone"]): t for t in tasks}
    complete = {key for key, count in counts.items()
                if key in task_by_key and count == expected_arms(task_by_key[key])}
    rows = [r for r in existing if (r["phase"], r["sample_id"], int(r["clone"])) in complete]
    todo = [t for t in tasks if (t["phase"], t["sample_id"], t["clone"]) not in complete]
    start = time.monotonic()
    if todo:
        context = mp.get_context("spawn")
        with context.Pool(min(workers, len(todo)), initializer=init_worker,
                          initargs=(str(SOURCE), adapter_text)) as pool:
            for i, result in enumerate(pool.imap_unordered(run_task, todo), 1):
                rows.extend(result); write_csv(partial, rows)
                print(f"completed {i}/{len(todo)} restart blocks; elapsed {time.monotonic()-start:.1f}s", flush=True)
    rows.sort(key=lambda r: (r["phase"], r["sample_id"], int(r["clone"]), r["arm"]))
    write_csv(output / "branches.csv", rows)
    if partial.exists(): partial.unlink()
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=RAW_DEFAULT)
    parser.add_argument("--output", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--workers", type=int, default=max(1, min(12, os.cpu_count() or 1)))
    parser.add_argument("--stage", choices=("audit", "all", "simulate"), default="all")
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    mod = load_gap20_production(SOURCE); para = parameters(mod)
    source_sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest(); step, text, adapter_sha = extract_step(mod, SOURCE)
    pre, post = choose_pre(), choose_post(); origins = setup_origins(args.raw_dir, step, para, pre, post)
    manifest = [{k: v for k, v in item.items() if k not in ("state", "dis_last")} for item in origins]
    write_csv(args.output / "selection_manifest.csv", manifest)
    audit_path = args.output / "simulator_audit.json"
    if args.stage in ("audit", "all"):
        result = audit(args.raw_dir, mod, step, para, origins, source_sha, adapter_sha)
        audit_path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        result = json.loads(audit_path.read_text())
    passed = result["parity"]["passed"] and result["optimized_voronoi"]["passed"]
    if not passed:
        raise RuntimeError("Restart/interaction audit failed; simulation stopped")
    print(json.dumps({"parity": result["parity"]["passed"],
                      "crossgap_nonzero": result["crossgap_social"]["nonzero"],
                      "crossgap_pairs": result["crossgap_social"]["total_admissible_pairs"]}, indent=2))
    if args.stage == "audit": return
    rows = simulate(origins, result, text, args.output, args.workers)
    metadata = dict(source_sha256=source_sha, adapter_sha256=adapter_sha, horizon_steps=HORIZON_STEPS,
                    horizon_time=HORIZON_STEPS * DT, clones=N_CLONES, pre_samples=len(pre), post_samples=len(post),
                    paired_noise="explicit PCG64 array shared by all arms of each sample/clone",
                    crossgap_arm_run=result["crossgap_social"]["nonzero"], rows=len(rows))
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
