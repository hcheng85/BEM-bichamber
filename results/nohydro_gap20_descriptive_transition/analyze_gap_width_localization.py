#!/usr/bin/env python3
"""Macroscopic localization and transport across archived opening widths.

Only chamber counts are extracted from existing trajectories.  The production
simulators and microscopic trajectories remain unchanged.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_initiation_propagation import RawReader


HERE = Path(__file__).resolve().parent
RAW_DEFAULT = Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
OUT_DEFAULT = HERE / "outputs/gap_width_localization"
LEGACY_GAP20_CACHE = HERE / "outputs/initiation_propagation/cache"
SOURCE_NOHYDRO = HERE.parents[1].parent / "BEM_bichamber_parallel_instance_no_hydro.py"
SOURCE_HYDRO = HERE.parents[1].parent / "BEM_bichamber_parallel_instance.py"
GAPS = (10, 20, 30, 40, 50, 60, 70, 80, 85)
N = 100
NATIVE_DT = 0.001
STRIDE = 100
DT = NATIVE_DT * STRIDE
TRANSIENT = 200.0
STRONG = 0.8
BASIN = 0.6
LAG_TIME = 10.0
FILE_RE = re.compile(
    r"If(?P<If>[^_]+)_Ip(?P<Ip>\d+)_In(?P<In>\d+)_R(?P<R>\d+)_num(?P<N>\d+)_wid(?P<gap>\d+)_n(?P<seed>\d+)_BEM(?P<nohydro>_no_hydro)?\.mat$"
)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def inventory(raw_dir: Path) -> list[dict]:
    rows = []
    for path in raw_dir.glob("*.mat"):
        match = FILE_RE.fullmatch(path.name)
        if not match:
            continue
        item = match.groupdict()
        rows.append(dict(path=str(path), filename=path.name,
                         dynamics="no_hydrodynamics" if item["nohydro"] else "hydrodynamics",
                         If=item["If"], Ip=int(item["Ip"]) / 10, In=int(item["In"]) / 10,
                         radius_ratio=int(item["R"]), N=int(item["N"]), gap_ratio=int(item["gap"]),
                         seed=int(item["seed"])))
    return sorted(rows, key=lambda r: (r["dynamics"], r["gap_ratio"], r["seed"]))


def compatible(rows: list[dict], dynamics: str) -> list[dict]:
    return [r for r in rows if r["dynamics"] == dynamics and r["Ip"] == 4 and r["In"] == .5
            and r["radius_ratio"] == 50 and r["N"] == N and r["gap_ratio"] in GAPS]


def cache_path(output: Path, item: dict) -> Path:
    tag = "nohydro" if item["dynamics"] == "no_hydrodynamics" else "hydro"
    return output / "cache" / f"{tag}_gap{item['gap_ratio']:02d}_seed{item['seed']:02d}_stride{STRIDE}.npz"


def extract_one(payload: tuple[dict, str]) -> dict:
    item, output_text = payload
    output = Path(output_text); cache = cache_path(output, item)
    if cache.exists():
        with np.load(cache) as data:
            q = np.asarray(data["q_count_difference"], np.int8)
            step = np.asarray(data["step"], np.int64)
        return dict(item=item, cache=str(cache), n_samples=len(q), q0=int(q[0]), cached=True)
    if item["dynamics"] == "no_hydrodynamics" and item["gap_ratio"] == 20:
        legacy = LEGACY_GAP20_CACHE / f"seed{item['seed']:02d}_stride100_qscan_v4.npz"
        if legacy.exists():
            with np.load(legacy) as data:
                step = np.asarray(data["step"], np.int64)
                q = np.rint(np.asarray(data["q"], float) * N).astype(np.int8)
            if q[0] != -N or np.any((q + N) % 2):
                raise RuntimeError(f"Legacy population-count audit failed for {legacy}")
            cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cache, step=step, q_count_difference=q)
            return dict(item=item, cache=str(cache), n_samples=len(q), q0=int(q[0]), cached=False,
                        reused_from=str(legacy))
    path = Path(item["path"]); reader = RawReader(path); info = reader.info["state"]
    if info.shape[0:2] != (3, N):
        raise RuntimeError(f"Unexpected state shape {info.shape} in {path}")
    step = np.arange(0, info.shape[-1] - 1, STRIDE, dtype=np.int64)
    state = np.memmap(path, dtype=info.dtype, mode="r", offset=info.data_offset, shape=info.shape, order="F")
    # Advanced indexing returns (time, fish); only sampled pages are touched.
    x = np.asarray(state[0, :, step])
    q = (2 * np.sum(x >= 0, axis=1) - N).astype(np.int8)
    del x, state; reader.close()
    if q[0] != -N or np.any(np.abs(q) > N) or np.any((q + N) % 2):
        raise RuntimeError(f"Population-count audit failed for {path}")
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, step=step, q_count_difference=q)
    return dict(item=item, cache=str(cache), n_samples=len(q), q0=int(q[0]), cached=False)


def extract(rows: list[dict], output: Path, workers: int) -> list[dict]:
    payload = [(item, str(output)) for item in rows]
    audit = []
    with mp.get_context("spawn").Pool(min(workers, len(payload))) as pool:
        for index, result in enumerate(pool.imap_unordered(extract_one, payload), 1):
            audit.append(result)
            print(f"extracted {index}/{len(payload)} {Path(result['item']['path']).name}", flush=True)
    return audit


def load_runs(rows: list[dict], output: Path) -> list[dict]:
    runs = []
    for item in rows:
        with np.load(cache_path(output, item)) as data:
            step = np.asarray(data["step"], np.int64)
            q = np.asarray(data["q_count_difference"], float) / N
        runs.append(dict(**item, time=step * NATIVE_DT, q=q))
    return runs


def mask_time(run: dict, lo: float, hi: float) -> np.ndarray:
    return (run["time"] >= lo - 1e-12) & (run["time"] <= hi + 1e-12)


def seed_boot(values: np.ndarray, seed: int, nboot: int = 5000) -> tuple[float, list[float]]:
    values = np.asarray(values, float); values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    draw = values[rng.integers(len(values), size=(nboot, len(values)))].mean(axis=1)
    return float(values.mean()), np.percentile(draw, [2.5, 97.5]).tolist()


def relocation_events(q: np.ndarray, dt: float = DT, bound: float = BASIN) -> list[dict]:
    current = 0; last_source = None; events = []
    for index, value in enumerate(q):
        basin = 1 if value >= bound else (-1 if value <= -bound else 0)
        if basin == current and basin:
            last_source = index
        elif basin and current and basin != current:
            events.append(dict(source_basin=current, target_basin=basin, departure_index=last_source,
                               arrival_index=index, duration=(index - last_source) * dt))
            current = basin; last_source = index
        elif basin and not current:
            current = basin; last_source = index
    return events


def distribution(runs: list[dict], lo: float, hi: float) -> np.ndarray:
    per_seed = []
    for run in runs:
        q = run["q"][mask_time(run, lo, hi)]
        count = np.bincount(np.rint((q + 1) * N / 2).astype(int), minlength=N + 1).astype(float)
        per_seed.append(count / count.sum())
    return np.mean(per_seed, axis=0)


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(.5 * np.abs(a - b).sum())


def per_seed_metrics(run: dict, lo: float, hi: float) -> dict:
    keep = mask_time(run, lo, hi); q = run["q"][keep]; lag = round(LAG_TIME / DT)
    q0, q1 = q[:-lag], q[lag:]; delta = np.abs(q1 - q0)
    out = dict(mean_q=float(np.mean(q)), mean_q2=float(np.mean(q * q)),
               strong=float(np.mean(np.abs(q) >= STRONG)),
               left_strong=float(np.mean(q <= -STRONG)), right_strong=float(np.mean(q >= STRONG)),
               mobility=float(np.mean(delta)))
    for label, low, high in (("central", 0., .4), ("intermediate", .4, .8), ("strong", .8, 1.01)):
        select = (np.abs(q0) >= low) & (np.abs(q0) < high)
        out[f"mobility_{label}"] = float(np.mean(delta[select])) if np.any(select) else np.nan
        out[f"mobility_{label}_samples"] = int(np.sum(select))
    return out


def summarize(runs: list[dict], lo: float, hi: float, window: str) -> list[dict]:
    rows = []
    grouped = defaultdict(list)
    for run in runs:
        grouped[(run["dynamics"], run["gap_ratio"])].append(run)
    metrics = ("mean_q", "mean_q2", "strong", "left_strong", "right_strong", "mobility",
               "mobility_central", "mobility_intermediate", "mobility_strong")
    for group_index, ((dynamics, gap), part) in enumerate(sorted(grouped.items()), 1):
        values = [per_seed_metrics(run, lo, hi) for run in part]
        row = dict(dynamics=dynamics, gap_ratio=gap, gap_over_radius=gap / 50,
                   window=window, n_seeds=len(part), samples_per_seed=int(mask_time(part[0], lo, hi).sum()))
        for metric_index, metric in enumerate(metrics):
            estimate, ci = seed_boot(np.asarray([v[metric] for v in values]), 20261201 + 100 * group_index + metric_index)
            row[metric] = estimate; row[metric + "_ci95"] = json.dumps(ci)
            if metric.startswith("mobility_"):
                row[metric + "_n_seeds"] = int(np.sum(np.isfinite([v[metric] for v in values])))
                row[metric + "_samples"] = int(np.sum([v.get(metric + "_samples", 0) for v in values]))
        row["q_distribution_tv_early_late"] = tv(distribution(part, 200, 600), distribution(part, 600, 1000)) if lo == 600 else ""
        rows.append(row)
    return rows


def relocations(runs: list[dict], lo: float = 200., hi: float = 1000.) -> tuple[list[dict], list[dict]]:
    events = []; summary = []
    grouped = defaultdict(list)
    for run in runs: grouped[(run["dynamics"], run["gap_ratio"])].append(run)
    for (dynamics, gap), part in sorted(grouped.items()):
        group_events = []
        for run in part:
            q = run["q"][mask_time(run, lo, hi)]
            for event_id, event in enumerate(relocation_events(q), 1):
                row = dict(dynamics=dynamics, gap_ratio=gap, seed=run["seed"], event_id=event_id,
                           departure_time=lo + event["departure_index"] * DT,
                           arrival_time=lo + event["arrival_index"] * DT, **event)
                events.append(row); group_events.append(row)
        summary.append(dict(dynamics=dynamics, gap_ratio=gap, n_seeds=len(part),
                            exposure_time=len(part) * (hi - lo), relocation_count=len(group_events),
                            seeds_with_relocation=len({r["seed"] for r in group_events}),
                            seeds_without_relocation=len(part) - len({r["seed"] for r in group_events}),
                            right_censored_terminal_intervals=len(part),
                            relocation_duration_median=float(np.median([r["duration"] for r in group_events])) if group_events else ""))
    return events, summary


def paired_contrasts(runs: list[dict], lo: float = 600., hi: float = 1000.) -> list[dict]:
    by_key={(r["dynamics"],r["gap_ratio"],r["seed"]):per_seed_metrics(r,lo,hi) for r in runs}
    requests=[("localization_peak_vs_narrow",50,10,"mean_q2"),
              ("localization_peak_vs_wide",50,85,"mean_q2"),
              ("strong_peak_vs_narrow",50,10,"strong"),
              ("strong_peak_vs_wide",50,85,"strong"),
              ("wide_vs_narrow_redistribution",85,10,"mobility"),
              ("wide_vs_narrow_central_q_band",85,10,"mobility_central"),
              ("wide_vs_narrow_intermediate_q_band",85,10,"mobility_intermediate"),
              ("wide_vs_narrow_strong_q_band",85,10,"mobility_strong")]
    rows=[]
    for dynamics in sorted({r["dynamics"] for r in runs}):
        for index,(name,gap_a,gap_b,metric) in enumerate(requests):
            seeds=sorted({r["seed"] for r in runs if r["dynamics"]==dynamics and r["gap_ratio"]==gap_a}
                         & {r["seed"] for r in runs if r["dynamics"]==dynamics and r["gap_ratio"]==gap_b})
            diff=np.asarray([by_key[(dynamics,gap_a,s)][metric]-by_key[(dynamics,gap_b,s)][metric] for s in seeds])
            estimate,ci=seed_boot(diff,20261280+index+(100 if dynamics=="hydrodynamics" else 0))
            rows.append(dict(dynamics=dynamics,contrast=name,gap_a=gap_a,gap_b=gap_b,metric=metric,
                             definition="gap_a minus gap_b, paired by archive seed",n_paired_seeds=len(seeds),
                             difference=estimate,difference_ci95=json.dumps(ci)))
    return rows


def figure(runs: list[dict], late: list[dict], relocation_summary: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5)); colors = {"no_hydrodynamics":"#2166ac", "hydrodynamics":"#b2182b"}
    for dynamics in sorted({r["dynamics"] for r in late}):
        part = sorted((r for r in late if r["dynamics"] == dynamics), key=lambda r:r["gap_ratio"])
        x = np.asarray([r["gap_ratio"] for r in part]); label = dynamics.replace("_", " ")
        for metric, ax in (("mean_q2", axes[0,0]), ("strong", axes[0,1]), ("mobility", axes[1,0])):
            y=np.asarray([r[metric] for r in part]); ci=np.asarray([json.loads(r[metric+"_ci95"]) for r in part])
            ax.errorbar(x,y,yerr=np.vstack((y-ci[:,0],ci[:,1]-y)),fmt="o-",capsize=3,color=colors[dynamics],label=label)
    axes[0,0].axhline(1/N,color="k",ls=":",label="independent $1/N$")
    axes[0,0].set(ylabel=r"$\langle q^2\rangle$",title="Macroscopic localization")
    axes[0,1].set(ylabel=r"$P(|q|\geq0.8)$",title="Strong-localization occupation")
    axes[1,0].set(ylabel=r"$\langle|q(t+10)-q(t)|\rangle$",title="Population redistribution")
    for ax in axes.flat[:3]: ax.set(xlabel="gap ratio (gap / body length)",xticks=GAPS); ax.legend(frameon=False,fontsize=8)
    ax=axes[1,1]; selected=(10,20,40,60,85); qvalues=np.arange(-N,N+1,2)/N
    nohydro=[r for r in runs if r["dynamics"]=="no_hydrodynamics"]
    for gap in selected:
        part=[r for r in nohydro if r["gap_ratio"]==gap]
        ax.plot(qvalues,distribution(part,600,1000),label=f"gap {gap}")
    ax.set(xlabel="q",ylabel="equal-seed P(q)",title="Late finite-time occupation, no hydrodynamics"); ax.legend(frameon=False,fontsize=8,ncol=2)
    fig.tight_layout(); fig.savefig(path,dpi=220); plt.close(fig)


def source_and_initialization_audit(rows: list[dict], output: Path) -> dict:
    selected = compatible(rows, "no_hydrodynamics")
    counts = defaultdict(int)
    for row in selected: counts[row["gap_ratio"]] += 1
    first_states = {}
    for row in selected:
        if row["seed"] == 0:
            reader=RawReader(Path(row["path"])); first_states[row["gap_ratio"]]=reader.frame("state",0); reader.close()
    reference=first_states[20]
    max_initial_error=max(float(np.max(np.abs(state-reference))) for state in first_states.values())
    audit=dict(nohydro_source=str(SOURCE_NOHYDRO), nohydro_source_sha256=hashlib.sha256(SOURCE_NOHYDRO.read_bytes()).hexdigest(),
               compatible_nohydro_counts={str(k):counts[k] for k in GAPS}, common_initial_state_seed0_max_error=max_initial_error,
               initialization="identical production all-left state across widths for a common seed",
               native_dt=NATIVE_DT, q_sample_dt=DT, common_duration=1000.0,
               physical_gap_over_radius={str(k):k/50 for k in GAPS})
    if SOURCE_HYDRO.exists():
        hydro = compatible(rows, "hydrodynamics"); hydro_counts = defaultdict(int)
        nohydro_lookup = {(r["gap_ratio"], r["seed"]):r for r in selected}
        initial_errors = []
        for row in hydro:
            hydro_counts[row["gap_ratio"]] += 1
            reference = nohydro_lookup.get((row["gap_ratio"], row["seed"]))
            if reference is None:
                continue
            reader=RawReader(Path(row["path"])); hydro_state=reader.frame("state",0); reader.close()
            reader=RawReader(Path(reference["path"])); nohydro_state=reader.frame("state",0); reader.close()
            initial_errors.append(float(np.max(np.abs(hydro_state-nohydro_state))))
        audit.update(hydro_source=str(SOURCE_HYDRO), hydro_source_sha256=hashlib.sha256(SOURCE_HYDRO.read_bytes()).hexdigest(),
                     compatible_hydro_counts={str(k):hydro_counts[k] for k in GAPS},
                     hydro_If=sorted({r["If"] for r in hydro}),
                     hydro_nohydro_common_initial_state_max_error=max(initial_errors),
                     hydro_comparison="same Ip, In, N, geometry, duration, seed initialization; hydrodynamic coupling If=1e-02 versus If=00")
    return audit


def self_test() -> None:
    q=np.r_[np.full(5,-.8),[-.4,0,.4],np.full(5,.8),[.4,0,-.4],np.full(5,-.8)]
    assert len(relocation_events(q)) == 2
    fake=dict(time=np.arange(301)*DT,q=np.linspace(-1,1,301))
    out=per_seed_metrics(fake,0,30); assert out["mean_q2"] > 0 and out["mobility"] > 0
    print("PASS: relocation state machine and fixed-lag macroscopic metrics")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--raw-dir",type=Path,default=RAW_DEFAULT)
    parser.add_argument("--output",type=Path,default=OUT_DEFAULT); parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--include-hydro",action="store_true"); parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test: self_test(); return
    args.output.mkdir(parents=True,exist_ok=True)
    rows=inventory(args.raw_dir); audit=source_and_initialization_audit(rows,args.output)
    selected=compatible(rows,"no_hydrodynamics")
    if args.include_hydro: selected += compatible(rows,"hydrodynamics")
    expected={10:40,20:40,30:20,40:20,50:20,60:20,70:20,80:20,85:20}
    actual=defaultdict(int)
    for r in compatible(rows,"no_hydrodynamics"): actual[r["gap_ratio"]]+=1
    if dict(actual) != expected: raise RuntimeError((dict(actual),expected))
    extraction=extract(selected,args.output,args.workers); runs=load_runs(selected,args.output)
    windows=[(0.,200.,"0-200"),(200.,600.,"200-600"),(600.,1000.,"600-1000"),(200.,1000.,"200-1000")]
    summaries=[]
    for lo,hi,name in windows: summaries.extend(summarize(runs,lo,hi,name))
    events, relocation_summary=relocations(runs)
    contrasts=paired_contrasts(runs)
    late=[r for r in summaries if r["window"]=="600-1000"]
    extraction_rows=[dict(**r["item"],cache=r["cache"],n_samples=r["n_samples"],q0=r["q0"],
                          cached=r["cached"],reused_from=r.get("reused_from","")) for r in extraction]
    write_csv(args.output/"archive_inventory.csv",rows); write_csv(args.output/"extraction_audit.csv",extraction_rows)
    write_csv(args.output/"trajectory_manifest.csv",selected)
    np.savez_compressed(args.output/"population_imbalance_trajectories.npz",
                        time=runs[0]["time"],
                        q_count_difference=np.stack([np.rint(r["q"]*N).astype(np.int8) for r in runs]))
    write_csv(args.output/"gap_width_macroscopic_summary.csv",summaries)
    write_csv(args.output/"complete_relocations.csv",events); write_csv(args.output/"relocation_censoring_summary.csv",relocation_summary)
    write_csv(args.output/"geometry_contrasts.csv",contrasts)
    figure(runs,late,relocation_summary,args.output/"gap_width_localization.png")
    audit["include_hydrodynamics"]=args.include_hydro; audit["n_analyzed_runs"]=len(runs)
    (args.output/"source_and_data_audit.json").write_text(json.dumps(audit,indent=2)+"\n")


if __name__ == "__main__": main()
