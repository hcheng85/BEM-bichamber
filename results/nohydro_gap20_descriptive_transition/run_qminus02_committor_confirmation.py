#!/usr/bin/env python3
"""Independent fixed-q=-0.2 full-state committor confirmation experiment."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammaln
from scipy.stats import beta

import run_bottleneck_commitment_branches as bc
import run_full_state_committor as pilot


HERE = Path(__file__).resolve().parent
OUT_DEFAULT = HERE / "outputs" / "qminus02_committor_confirmation"
ORIGINAL_OUT = HERE / "outputs" / "full_state_committor"
HISTORICAL_DEFAULT = OUT_DEFAULT / "historical_completed_qminus02_shots.csv"
STATE_SEEDS = (8, 12, 14, 15)
CLONES = tuple(range(8, 24))
NOISE_BASE = 202611270000
HORIZONS = (600.0, 1200.0)
CHECKPOINT_INTERVAL = 25.0
FAMILY_ALPHA = 0.05


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def microscopic_hash(state: np.ndarray, memory: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(state).tobytes())
    digest.update(np.ascontiguousarray(memory).tobytes())
    return digest.hexdigest()


def noise_seed(seed: int, clone: int) -> int:
    return NOISE_BASE + 100 * seed + clone


def population(state: np.ndarray) -> tuple[int, float]:
    nright = int(np.count_nonzero(state[0] >= 0))
    if not 0 <= nright <= bc.N:
        raise RuntimeError("Population conservation failed")
    return nright, (2 * nright - bc.N) / bc.N


def checkpoint_path(output: Path, state_id: str, clone: int) -> Path:
    return output / "cache" / "checkpoints" / f"{state_id}_clone{clone:03d}.npz"


def save_checkpoint(path: Path, state: np.ndarray, memory: np.ndarray,
                    rng: np.random.Generator, step: int, outcome: str,
                    provenance: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(provenance, step=step, outcome=outcome,
                    rng_state=rng.bit_generator.state)
    temporary = path.with_suffix(".npz.tmp")
    with temporary.open("wb") as f:
        np.savez_compressed(f, state=state, memory=memory,
                            metadata=np.array(json.dumps(metadata)))
    os.replace(temporary, path)


def load_checkpoint(path: Path, expected: dict):
    with np.load(path, allow_pickle=False) as data:
        state = np.asarray(data["state"], float)
        memory = np.asarray(data["memory"], float)
        metadata = json.loads(str(data["metadata"].item()))
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise RuntimeError(f"Checkpoint provenance mismatch for {key}: {path}")
    if state.shape != (3, bc.N) or memory.shape != (bc.N,):
        raise RuntimeError(f"Invalid checkpoint shapes: {path}")
    if not np.all(np.isfinite(state)) or not np.all(np.isfinite(memory)):
        raise RuntimeError(f"Non-finite checkpoint: {path}")
    rng = np.random.default_rng()
    rng.bit_generator.state = metadata["rng_state"]
    return state, memory, rng, int(metadata["step"]), metadata["outcome"]


def integrate(state: np.ndarray, memory: np.ndarray, rng: np.random.Generator,
              start_step: int, stop_step: int):
    outcome = "censored"
    if stop_step <= start_step:
        return state, memory, start_step, outcome
    for step in range(start_step + 1, stop_step + 1):
        increment = np.zeros((3, bc.N))
        increment[2] = bc._PARA["In"] * math.sqrt(bc.DT) * rng.normal(size=bc.N)
        state, _, _, memory, _ = bc._STEP(
            state, bc._PARA, increment, memory, False)
        _, q = population(state)
        if q <= pilot.A_BOUND:
            outcome = "A"
            break
        if q >= pilot.B_BOUND:
            outcome = "B"
            break
    return state, memory, step, outcome


def run_one(task: dict) -> dict:
    path = Path(task["checkpoint"])
    expected_keys = ("experiment", "state_id", "archive_seed", "clone",
                     "noise_seed", "state_hash", "source_sha", "adapter_sha")
    expected = {key: task[key] for key in expected_keys}
    target_step = int(task["target_step"])
    interval_steps = int(task["interval_steps"])
    if path.exists():
        state, memory, rng, step, outcome = load_checkpoint(path, expected)
    else:
        state = task["state"].copy()
        memory = task["memory"].copy()
        rng = np.random.default_rng(task["noise_seed"])
        step, outcome = 0, "censored"
        save_checkpoint(path, state, memory, rng, step, outcome, expected)
    while outcome == "censored" and step < target_step:
        stop = min(step + interval_steps, target_step)
        state, memory, step, outcome = integrate(state, memory, rng, step, stop)
        save_checkpoint(path, state, memory, rng, step, outcome, expected)
    nright, q = population(state)
    if outcome == "A" and q > pilot.A_BOUND:
        raise RuntimeError("Invalid A terminal state")
    if outcome == "B" and q < pilot.B_BOUND:
        raise RuntimeError("Invalid B terminal state")
    if outcome == "censored" and step < target_step:
        raise RuntimeError("Unresolved shot stopped before requested horizon")
    return dict(state_id=task["state_id"], seed=task["archive_seed"],
                clone=task["clone"], noise_seed=task["noise_seed"],
                checkpoint_step=step, checkpoint_time=step * bc.DT,
                checkpoint_outcome=outcome, nright=nright)


def state_at_horizon(checkpoint: Path, expected: dict, horizon_step: int) -> dict:
    _, _, _, step, outcome = load_checkpoint(checkpoint, expected)
    if step < horizon_step and outcome == "censored":
        raise RuntimeError(f"Checkpoint has not reached horizon: {checkpoint}")
    if outcome in ("A", "B") and step <= horizon_step:
        result = outcome
        hit_step = step
    else:
        result = "censored"
        hit_step = ""
    return dict(outcome=result, hit_step=hit_step,
                hit_time=(hit_step * bc.DT if hit_step != "" else ""))


def logchoose(n: int, k: int) -> float:
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def exact_equal_binomial(counts: list[int], n_each: int = 16) -> dict:
    """Two-sided Fisher-Freeman-Halton test for four equal binomials."""
    total = int(sum(counts))
    observed_logp = sum(logchoose(n_each, k) for k in counts) - logchoose(4*n_each, total)
    probability = 0.0
    tables = 0
    for a in range(n_each + 1):
        for b in range(n_each + 1):
            for c in range(n_each + 1):
                d = total - a - b - c
                if not 0 <= d <= n_each:
                    continue
                logp = (logchoose(n_each, a) + logchoose(n_each, b) +
                        logchoose(n_each, c) + logchoose(n_each, d) -
                        logchoose(4*n_each, total))
                if logp <= observed_logp + 1e-12:
                    probability += math.exp(logp)
                tables += 1
    return dict(test="exact conditional 4x2 equal-probability test",
                counts=counts, total=total, p_value=min(1.0, probability),
                enumerated_tables=tables)


def pearson_stat(table: np.ndarray) -> float:
    pooled = table.sum(axis=0) / table.sum()
    expected = table.sum(axis=1)[:, None] * pooled[None, :]
    mask = expected > 0
    return float(np.sum((table[mask] - expected[mask])**2 / expected[mask]))


def permutation_three_category(table: np.ndarray, nperm: int = 100000) -> dict:
    observed = pearson_stat(table)
    labels = np.repeat(np.arange(3), table.sum(axis=0))
    sizes = table.sum(axis=1)
    rng = np.random.default_rng(2026112701)
    exceed = 0
    for _ in range(nperm):
        shuffled = rng.permutation(labels)
        cursor = 0
        permuted = []
        for size in sizes:
            group = shuffled[cursor:cursor+size]
            permuted.append(np.bincount(group, minlength=3))
            cursor += size
        exceed += pearson_stat(np.asarray(permuted)) >= observed - 1e-12
    return dict(test="prespecified permutation Pearson test for A/B/censored",
                statistic=observed, permutations=nperm,
                p_value=(exceed + 1)/(nperm + 1))


def clopper(k: int, n: int, alpha: float = .05) -> list[float]:
    lower = 0.0 if k == 0 else float(beta.ppf(alpha/2, k, n-k+1))
    upper = 1.0 if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))
    return [lower, upper]


def simultaneous_committor_bounds(b: int, c: int, n: int) -> list[float]:
    """Bonferroni 95% bounds for all two-sided state bounds (8 tails)."""
    tail = FAMILY_ALPHA / 8
    lower = 0.0 if b == 0 else float(beta.ppf(tail, b, n-b+1))
    upper_successes = b + c
    upper = 1.0 if upper_successes == n else float(
        beta.ppf(1-tail, upper_successes+1, n-upper_successes))
    return [lower, upper]


def summarize(shots: list[dict], horizon: float) -> tuple[list[dict], dict]:
    grouped = defaultdict(list)
    for row in shots:
        grouped[row["state_id"]].append(row)
    rows = []
    for state_id in sorted(grouped, key=lambda x: int(x.split("_s")[1].split("_")[0])):
        group = grouped[state_id]
        n = len(group)
        a = sum(row["outcome"] == "A" for row in group)
        b = sum(row["outcome"] == "B" for row in group)
        c = n-a-b
        simultaneous = simultaneous_committor_bounds(b, c, n)
        rows.append(dict(
            state_id=state_id, seed=group[0]["seed"], horizon=horizon,
            shots=n, A=a, B=b, censored=c,
            F_A=a/n, F_B=b/n, survival=c/n,
            empirical_committor_lower=b/n,
            empirical_committor_upper=(b+c)/n,
            simultaneous_confidence_lower=simultaneous[0],
            simultaneous_confidence_upper=simultaneous[1],
            eventual_committor_estimate=(b/n if c == 0 else ""),
            binomial_ci95=(json.dumps(clopper(b, n)) if c == 0 else "not applicable while censored")))
    table = np.asarray([[row["A"], row["B"], row["censored"]] for row in rows], int)
    if table[:, 2].sum() == 0:
        primary = exact_equal_binomial(table[:, 1].tolist())
        primary["estimand"] = "equality of eventual B-before-A committors"
    else:
        primary = permutation_three_category(table)
        primary["estimand"] = "equality of finite-horizon A/B/unresolved distributions"
    tests = dict(
        horizon=horizon,
        primary=primary,
        B_by_horizon=exact_equal_binomial(table[:, 1].tolist()),
        unresolved_by_horizon=exact_equal_binomial(table[:, 2].tolist()),
        eventual_test_performed=bool(table[:, 2].sum() == 0),
        interpretation=(
            "All shots resolved: the primary exact test targets equality of eventual B-before-A committors."
            if table[:, 2].sum() == 0 else
            "Resolved-only committors are not estimated while shots are censored; the primary test targets finite-horizon cumulative incidence and survival."))
    return rows, tests


def completed_original(path: Path, state_ids: set[str]) -> list[dict]:
    records = [row for row in read_csv(path) if row["state_id"] in state_ids]
    if len(records) != 32:
        raise RuntimeError(f"Expected 32 completed historical shots, got {len(records)}")
    for state_id in state_ids:
        group = [row for row in records if row["state_id"] == state_id]
        if sorted(int(row["clone"]) for row in group) != list(range(8)):
            raise RuntimeError(f"Historical clone IDs changed for {state_id}")
        if any(row["outcome"] not in ("A", "B") for row in group):
            raise RuntimeError(f"Historical shots not completed for {state_id}")
    return [dict(state_id=row["state_id"], seed=int(row["seed"]),
                 clone=int(row["clone"]), outcome=row["outcome"],
                 hit_time=float(row["hit_time"])) for row in records]


def combined_summary(original: list[dict], new: list[dict], horizon: float) -> list[dict]:
    rows = []
    for state_id in sorted({row["state_id"] for row in original},
                           key=lambda x: int(x.split("_s")[1].split("_")[0])):
        old = [row for row in original if row["state_id"] == state_id]
        fresh = [row for row in new if row["state_id"] == state_id]
        group = old + fresh
        n=len(group); a=sum(r["outcome"]=="A" for r in group)
        b=sum(r["outcome"]=="B" for r in group); c=n-a-b
        rows.append(dict(state_id=state_id,seed=old[0]["seed"],horizon=horizon,
                         original_shots=len(old),confirmatory_shots=len(fresh),shots=n,
                         A=a,B=b,censored=c,empirical_lower=b/n,
                         empirical_upper=(b+c)/n,
                         combined_estimate=(b/n if c==0 else ""),
                         binomial_ci95=(json.dumps(clopper(b,n)) if c==0 else "not applicable while censored"),
                         interpretation="descriptive combination; not independent confirmatory evidence"))
    return rows


def make_figure(all_state_rows: list[dict], combined: list[dict], path: Path) -> None:
    horizons=sorted({row["horizon"] for row in all_state_rows})
    fig, axes=plt.subplots(1,len(horizons)+1,figsize=(5*(len(horizons)+1),4),sharey=True)
    for ax,horizon in zip(axes,horizons):
        rows=[row for row in all_state_rows if row["horizon"]==horizon]
        x=np.arange(len(rows));bottom=np.zeros(len(rows))
        for key,color in (("A","#4477AA"),("B","#CC6677"),("censored","#BBBBBB")):
            values=np.asarray([row[key]/row["shots"] for row in rows])
            ax.bar(x,values,bottom=bottom,label=key,color=color);bottom+=values
        ax.set(title=f"New shots at t={horizon:g}",xticks=x,
               xticklabels=[f"seed {row['seed']}" for row in rows],ylim=(0,1),ylabel="fraction")
        ax.legend(frameon=False)
    ax=axes[-1];x=np.arange(len(combined))
    lower=np.asarray([row["empirical_lower"] for row in combined]);upper=np.asarray([row["empirical_upper"] for row in combined])
    if all(row["combined_estimate"] != "" for row in combined):
        point=np.asarray([row["combined_estimate"] for row in combined],float)
        intervals=np.asarray([json.loads(row["binomial_ci95"]) for row in combined],float)
        ax.errorbar(x,point,yerr=np.vstack((point-intervals[:,0],intervals[:,1]-point)),
                    fmt="o",color="#CC6677",capsize=4)
        title="Combined 24-shot estimates\n(descriptive; 95% binomial CI)"
    else:
        ax.vlines(x,lower,upper,color="k",lw=3);ax.scatter(x,lower,color="#CC6677",zorder=3)
        title="Combined 24-shot bounds\n(descriptive)"
    ax.set(title=title,xticks=x,
           xticklabels=[f"seed {row['seed']}" for row in combined],ylim=(0,1),ylabel="B-before-A probability")
    fig.tight_layout();fig.savefig(path,dpi=220);plt.close(fig)


def prepare(args):
    inventory, selected = pilot.inventory_and_select(args.cache_q)
    del inventory
    selected = [row for row in selected if np.isclose(row["q_level"], -.2)]
    if tuple(sorted(row["seed"] for row in selected)) != STATE_SEEDS:
        raise RuntimeError("The four committed q=-0.2 configurations changed")
    mod=bc.load_gap20_production(bc.SOURCE);para=bc.parameters(mod)
    step,adapter_text,adapter_sha=bc.extract_step(mod,bc.SOURCE)
    origins=pilot.load_origins(args.raw_dir,selected,step,para)
    audit=pilot.audit(args.raw_dir,mod,step,para,origins,adapter_sha)
    previous=json.loads((ORIGINAL_OUT/"restart_audit.json").read_text())
    if not audit["native_replay_passed"] or not audit["reflection_one_step_passed"]:
        raise RuntimeError("Native replay or reflection audit failed")
    source_sha=sha256(bc.SOURCE)
    if source_sha != previous["source_sha256"] or adapter_sha != previous["adapter_sha256"]:
        raise RuntimeError("Production source or restart adapter differs from original pilot")
    original_pilot=read_csv(ORIGINAL_OUT/"shots.csv")
    old_seeds={int(r["noise_seed"]) for r in original_pilot}
    new_seeds={noise_seed(seed,clone) for seed in STATE_SEEDS for clone in CLONES}
    if old_seeds & new_seeds or len(new_seeds)!=64:
        raise RuntimeError("Noise seeds are not independent/nonoverlapping")
    return origins,adapter_text,adapter_sha,source_sha,audit,new_seeds


def analyze(args,origins,adapter_sha,source_sha,new_seeds):
    state_ids={row["state_id"] for row in origins}
    historical=completed_original(args.original_completed,state_ids)
    write_csv(args.output/"historical_completed_qminus02_shots.csv",historical)
    all_shots=[];state_rows=[];tests=[]
    requested=[h for h in HORIZONS if h<=args.target_time]
    for horizon in requested:
        horizon_step=round(horizon/bc.DT);shots=[]
        for item in origins:
            state_hash=microscopic_hash(item["state"],item["dis_last"])
            for clone in CLONES:
                expected=dict(experiment="qminus02_committor_confirmation_v1",
                              state_id=item["state_id"],archive_seed=item["seed"],clone=clone,
                              noise_seed=noise_seed(item["seed"],clone),state_hash=state_hash,
                              source_sha=source_sha,adapter_sha=adapter_sha)
                record=state_at_horizon(checkpoint_path(args.output,item["state_id"],clone),expected,horizon_step)
                shots.append(dict(state_id=item["state_id"],seed=item["seed"],clone=clone,
                                  noise_seed=expected["noise_seed"],horizon=horizon,**record))
        rows,test=summarize(shots,horizon);all_shots.extend(shots);state_rows.extend(rows);tests.append(test)
        write_csv(args.output/f"confirmatory_shots_t{int(horizon)}.csv",shots)
        write_csv(args.output/f"confirmatory_state_summary_t{int(horizon)}.csv",rows)
        (args.output/f"confirmatory_tests_t{int(horizon)}.json").write_text(json.dumps(test,indent=2)+"\n")
    latest=max(requested);latest_shots=[r for r in all_shots if r["horizon"]==latest]
    combined=combined_summary(historical,latest_shots,latest)
    write_csv(args.output/"combined_24shot_summary.csv",combined)
    make_figure(state_rows,combined,args.output/"committor_confirmation.png")
    manifest=dict(experiment="qminus02_committor_confirmation_v1",states=sorted(state_ids),
                  clones=list(CLONES),noise_seed_rule="202611270000 + 100*archive_seed + clone",
                  n_unique_noise_seeds=len(new_seeds),basin_A="q <= -0.6",basin_B="q >= +0.6",
                  dt=bc.DT,horizons=requested,checkpoint_interval=args.checkpoint_interval,
                  production_source_sha256=source_sha,restart_adapter_sha256=adapter_sha,
                  original_pilot_shots_sha256=sha256(ORIGINAL_OUT/"shots.csv"),
                  original_completed_shots_sha256=sha256(args.original_completed),
                  confirmatory_tests=tests,
                  separation="Only clones 8-23 enter confirmatory tests; combined 24-shot summaries are descriptive.")
    (args.output/"manifest_and_results.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest


def self_test() -> None:
    assert len({noise_seed(seed,clone) for seed in STATE_SEEDS for clone in CLONES})==64
    strong=exact_equal_binomial([0,16,0,16]);equal=exact_equal_binomial([8,8,8,8])
    assert strong["p_value"]<1e-8 and abs(equal["p_value"]-1)<1e-12
    bounds=simultaneous_committor_bounds(2,3,16)
    assert 0<=bounds[0]<=2/16<=(2+3)/16<=bounds[1]<=1
    rng=np.random.default_rng(17);rng.normal(size=30)
    other=np.random.default_rng();other.bit_generator.state=rng.bit_generator.state
    assert np.array_equal(rng.normal(size=50),other.normal(size=50))
    print("PASS: noise uniqueness, exact homogeneity test, censor bounds, RNG continuation")


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--raw-dir",type=Path,default=pilot.RAW_DEFAULT)
    parser.add_argument("--cache-q",type=Path,default=pilot.CACHE_Q)
    parser.add_argument("--output",type=Path,default=OUT_DEFAULT)
    parser.add_argument("--original-completed",type=Path,default=HISTORICAL_DEFAULT)
    parser.add_argument("--target-time",type=float,choices=HORIZONS,default=600.)
    parser.add_argument("--checkpoint-interval",type=float,default=CHECKPOINT_INTERVAL)
    parser.add_argument("--workers",type=int,default=9)
    parser.add_argument("--stage",choices=("audit","shoot","analyze"),default="shoot")
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:self_test();return
    if args.workers<1 or args.checkpoint_interval<=0:raise ValueError("Invalid workers/checkpoint interval")
    if abs(round(args.checkpoint_interval/bc.DT)*bc.DT-args.checkpoint_interval)>1e-12:
        raise ValueError("Checkpoint interval must be a multiple of dt")
    args.output.mkdir(parents=True,exist_ok=True)
    origins,adapter_text,adapter_sha,source_sha,audit,new_seeds=prepare(args)
    audit_output=dict(native_replay_passed=audit["native_replay_passed"],
                      max_native_replay_error=audit["max_native_replay_error"],
                      reflection_one_step_passed=audit["reflection_one_step_passed"],
                      max_reflection_error=audit["max_reflection_error"],
                      previous_step_memory=audit["previous_step_memory"],states=[x["state_id"] for x in origins])
    (args.output/"restart_audit.json").write_text(json.dumps(audit_output,indent=2)+"\n")
    if args.stage=="audit":
        print(json.dumps(audit_output,indent=2));return
    target_step=round(args.target_time/bc.DT);interval_steps=round(args.checkpoint_interval/bc.DT)
    if args.stage=="shoot":
        jobs=[]
        # Round-robin states so an interrupted run retains balanced coverage.
        for clone in CLONES:
            for item in origins:
                state_hash=microscopic_hash(item["state"],item["dis_last"])
                jobs.append(dict(experiment="qminus02_committor_confirmation_v1",
                                 state_id=item["state_id"],archive_seed=item["seed"],clone=clone,
                                 noise_seed=noise_seed(item["seed"],clone),state_hash=state_hash,
                                 source_sha=source_sha,adapter_sha=adapter_sha,
                                 state=item["state"],memory=item["dis_last"],target_step=target_step,
                                 interval_steps=interval_steps,
                                 checkpoint=str(checkpoint_path(args.output,item["state_id"],clone))))
        start=time.monotonic()
        with mp.get_context("spawn").Pool(min(args.workers,len(jobs)),initializer=bc.init_worker,
                                          initargs=(str(bc.SOURCE),adapter_text)) as pool:
            for index,result in enumerate(pool.imap_unordered(run_one,jobs),1):
                print(f"completed {index}/{len(jobs)} {result['state_id']} clone {result['clone']} "
                      f"{result['checkpoint_outcome']} t={result['checkpoint_time']:.3f}; "
                      f"elapsed={time.monotonic()-start:.1f}s",flush=True)
    manifest=analyze(args,origins,adapter_sha,source_sha,new_seeds)
    print(json.dumps(dict(horizons=manifest["horizons"],tests=manifest["confirmatory_tests"]),indent=2))


if __name__=="__main__":main()
