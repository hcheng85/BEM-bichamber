#!/usr/bin/env python3
"""Exact approach-region budgets for continued and failed crossing attempts.

This reuses the validated piecewise-linear conservation accounting without
modifying the simulator.  The source chamber is fixed from the pre-crossing
event and never reassigned after fish begin to transfer.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_approach_current_balance import DT, ELL, N, RawReader, raw_path
from analyze_conservation_correction import segment_budget


Q_EDGES = np.asarray([0.0, 0.5, 0.8, 1.000001])
WINDOW = 0.2
START = -3.0
STOP = 2.0
TERMS = (
    "plug_flux",
    "storage_contribution",
    "lateral_contribution",
    "neck_flux",
    "path_discretization_correction",
)
PERIODS = (
    ("early_pre", -3.0, -1.0),
    ("late_pre", -1.0, 0.0),
    ("immediate_post", 0.0, 1.0),
    ("later_post", 1.0, 2.0),
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def q_bin(q: float) -> int:
    idx = int(np.searchsorted(Q_EDGES, q, side="right") - 1)
    if not 0 <= idx < len(Q_EDGES) - 1:
        raise ValueError(f"Initial |q| outside declared strata: {q}")
    return idx


def load_cohort(path: Path) -> list[dict]:
    rows = []
    for row in read_csv(path):
        source_sign = int(row["source_sign"])
        rows.append(dict(seed=int(row["seed"]), event_id=int(row["event_id"]),
                         first_native=int(row["first_native"]),
                         source_sign=source_sign, direction=-source_sign,
                         success4=int(row["success4"]),
                         additional_fish=int(float(row["additional_fish"])),
                         abs_q_pre=float(row["abs_q_pre"]),
                         q_bin=q_bin(float(row["abs_q_pre"]))))
    keys={(r["seed"],r["event_id"]) for r in rows}
    if len(keys) != len(rows):
        raise RuntimeError("Duplicate event in fixed cohort")
    return rows


def extract_budgets(raw_dir: Path, cohort: list[dict]) -> list[dict]:
    """Measure nonoverlapping native-step budgets relative to first entry.

    A bin labelled by its center (lo+hi)/2 contains segments with relative end
    times in (lo,hi].  Consequently the first crossing lies in (-0.2,0].
    """
    edges = np.arange(START, STOP + WINDOW / 2, WINDOW)
    rows=[]
    by_seed=defaultdict(list)
    for event in cohort:
        by_seed[event["seed"]].append(event)
    for count,seed in enumerate(sorted(by_seed),1):
        print(f"exact event budgets: seed {count}/{len(by_seed)} ({seed})",flush=True)
        reader=RawReader(raw_path(raw_dir,seed))
        try:
            total=reader.info["state"].shape[-1]
            for event in by_seed[seed]:
                anchor=event["first_native"]
                first=anchor+int(round(START/DT))
                last=anchor+int(round(STOP/DT))
                if first < 0 or last+1 >= total:
                    raise RuntimeError(f"Event window out of range: {event}")
                states=reader.block("state",first,last)
                velocities=reader.block("rdot",first+1,last)
                anchor_offset=anchor-first
                pre=states[:,:,anchor_offset-1];at=states[:,:,anchor_offset]
                q_pre=float(2*np.count_nonzero(pre[0]>=0)/N-1)
                if abs(q_pre-event["source_sign"]*event["abs_q_pre"])>1e-10:
                    raise RuntimeError(f"Pre-event population mismatch: {event}, q={q_pre}")
                nr_pre=int(np.count_nonzero(pre[0]>=0));nr_at=int(np.count_nonzero(at[0]>=0))
                if -event["source_sign"]*(nr_at-nr_pre)<=0:
                    raise RuntimeError(f"Indexed first crossing has wrong direction: {event}")
                for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                    offset=int(round((lo-START)/DT))
                    steps=int(round((hi-lo)/DT))
                    budget=segment_budget(states[:,:,offset:offset+steps+1],
                                          velocities[:,:,offset:offset+steps],
                                          event["direction"])
                    rows.append({**event,"tau_start":float(lo),"tau_end":float(hi),
                                 "tau_center":float((lo+hi)/2),
                                 "crossing_bin":int(lo < 0 <= hi),
                                 **budget,
                                 "storage_contribution":-budget["weighted_accumulation"],
                                 "lateral_contribution":-budget["weighted_lateral_flux"],
                                 "supply_plus_lateral":budget["plug_flux"]-budget["weighted_lateral_flux"],
                                 "lateral_minus_supply":-budget["weighted_lateral_flux"]-budget["plug_flux"]})
        finally:
            reader.close()
    return rows


def target_weights(cohort: list[dict]) -> np.ndarray:
    counts=np.bincount([r["q_bin"] for r in cohort],minlength=len(Q_EDGES)-1)
    return counts/counts.sum()


def seed_cell_means(rows: list[dict], term: str) -> dict[tuple[int,int,int],float]:
    grouped=defaultdict(list)
    for row in rows:
        grouped[(int(row["success4"]),int(row["q_bin"]),int(row["seed"]))].append(float(row[term]))
    return {key:float(np.mean(values)) for key,values in grouped.items()}


def standardized(cell: dict, outcome: int, weights: np.ndarray,
                 multiplicity: dict[int,int] | None = None) -> float:
    total=0.0
    for qb,w in enumerate(weights):
        if w == 0:
            continue
        values=[];sample_weights=[]
        for (status,cell_q,seed),value in cell.items():
            if status==outcome and cell_q==qb:
                mult=1 if multiplicity is None else multiplicity.get(seed,0)
                if mult:
                    values.append(value);sample_weights.append(mult)
        if not values:
            return float("nan")
        total += float(w)*float(np.average(values,weights=sample_weights))
    return total


def bootstrap_contrast(rows: list[dict], term: str, weights: np.ndarray,
                       nboot: int=2000, rng_seed: int=20261118) -> dict:
    cell=seed_cell_means(rows,term)
    continued=standardized(cell,1,weights);failed=standardized(cell,0,weights)
    seeds=np.asarray(sorted({int(r["seed"]) for r in rows}))
    rng=np.random.default_rng(rng_seed)
    draws=[]
    for _ in range(nboot):
        sample=rng.choice(seeds,len(seeds),replace=True)
        unique,counts=np.unique(sample,return_counts=True)
        mult=dict(zip(unique.tolist(),counts.tolist()))
        c=standardized(cell,1,weights,mult);f=standardized(cell,0,weights,mult)
        if np.isfinite(c) and np.isfinite(f):
            draws.append((c,f,c-f))
    arr=np.asarray(draws)
    if len(arr) < nboot*.9:
        raise RuntimeError(f"Insufficient supported seed-bootstrap draws: {len(arr)}/{nboot}")
    intervals=np.percentile(arr,[2.5,97.5],axis=0)
    return dict(continued=continued,failed=failed,difference=continued-failed,
                continued_ci95=intervals[:,0].tolist(),failed_ci95=intervals[:,1].tolist(),
                difference_ci95=intervals[:,2].tolist(),bootstrap_draws=len(arr))


def summarize_time(rows: list[dict], weights: np.ndarray) -> list[dict]:
    out=[]
    for tau in sorted({float(r["tau_center"]) for r in rows}):
        part=[r for r in rows if np.isclose(float(r["tau_center"]),tau)]
        for term in TERMS:
            stats=bootstrap_contrast(part,term,weights,rng_seed=20261118+round((tau-START)*100)+TERMS.index(term))
            out.append(dict(tau_center=tau,term=term,**stats))
    return out


def summarize_periods(rows: list[dict], weights: np.ndarray) -> list[dict]:
    out=[]
    period_terms=TERMS+("supply_plus_lateral","lateral_minus_supply")
    for period,lo,hi in PERIODS:
        # Each event contributes its mean rate over the period before seed averaging.
        part=[r for r in rows if float(r["tau_center"])>lo and float(r["tau_center"])<=hi]
        for term in period_terms:
            stats=bootstrap_contrast(part,term,weights,rng_seed=20262000+PERIODS.index((period,lo,hi))*20+period_terms.index(term))
            out.append(dict(period=period,tau_start=lo,tau_end=hi,term=term,**stats))
    return out


def summarize_q_strata(rows: list[dict]) -> list[dict]:
    out=[]
    for period,lo,hi in PERIODS:
        if period not in ("late_pre","immediate_post"):
            continue
        part=[r for r in rows if float(r["tau_center"])>lo and float(r["tau_center"])<=hi]
        for qb in range(len(Q_EDGES)-1):
            weights=np.zeros(len(Q_EDGES)-1);weights[qb]=1
            for term in ("plug_flux","storage_contribution","lateral_contribution","neck_flux"):
                stats=bootstrap_contrast(part,term,weights,
                                         rng_seed=20263000+PERIODS.index((period,lo,hi))*20+qb*5+TERMS.index(term))
                out.append(dict(period=period,q_bin=qb,q_low=Q_EDGES[qb],q_high=Q_EDGES[qb+1],
                                term=term,**stats))
    return out


def make_figure(summary: list[dict], output: Path) -> None:
    labels={"plug_flux":"directed supply  $N J_a/\\ell$",
            "storage_contribution":"storage release  $-A_W$",
            "lateral_contribution":"lateral contribution  $-L_W$",
            "neck_flux":"net neck flux  $F_{neck}$"}
    fig,axes=plt.subplots(2,2,figsize=(9.2,6.8),sharex=True)
    for ax,term in zip(axes.ravel(),labels):
        part=sorted((r for r in summary if r["term"]==term),key=lambda r:r["tau_center"])
        x=np.asarray([r["tau_center"] for r in part])
        for status,color,name in (("continued","tab:blue","continued (>=4 additional)"),
                                  ("failed","tab:orange","failed")):
            y=np.asarray([r[status] for r in part]);ci=np.asarray([r[f"{status}_ci95"] for r in part])
            ax.plot(x,y,color=color,lw=1.8,label=name);ax.fill_between(x,ci[:,0],ci[:,1],color=color,alpha=.18,lw=0)
        ax.axhline(0,color=".55",lw=.7);ax.axvline(0,color="k",ls=":",lw=1)
        ax.set_title(labels[term]);ax.set_ylabel("fish / time")
    for ax in axes[-1]: ax.set_xlabel("time from first crossing")
    axes[0,0].legend(frameon=False,fontsize=8,loc="upper left")
    fig.suptitle("Exact source-fixed approach budgets, standardized to common initial |q|",fontsize=11)
    fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    rows=[]
    for seed in range(6):
        for qb in range(3):
            for outcome in (0,1):
                rows.append(dict(seed=seed,q_bin=qb,success4=outcome,value=10*qb+outcome))
    result=bootstrap_contrast(rows,"value",np.asarray([.2,.3,.5]),nboot=200)
    assert abs(result["difference"]-1)<1e-12
    # One crossing segment must close even when the stored velocity is arbitrary;
    # the separately retained path correction repairs that convention difference.
    states=np.zeros((3,100,2));states[0,:,:]=-2*ELL
    states[0,0,0]=-.1;states[0,0,1]=.1
    velocities=np.zeros((3,100,1))
    budget=segment_budget(states,velocities,1)
    assert abs(budget["closure_residual"])<1e-10 and budget["net_neck_count"]==1
    print("PASS: fixed-source conservative budget and q-standardized seed contrast")


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--raw-dir",type=Path)
    parser.add_argument("--cohort",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        self_test();return
    if not all((args.raw_dir,args.cohort,args.output)):
        parser.error("--raw-dir, --cohort, and --output are required")
    args.output.mkdir(parents=True,exist_ok=True)
    cohort=load_cohort(args.cohort)
    if len(cohort)!=337:
        raise RuntimeError(f"Expected fixed 337-event cohort, found {len(cohort)}")
    rows=extract_budgets(args.raw_dir,cohort)
    weights=target_weights(cohort)
    time_summary=summarize_time(rows,weights)
    period_summary=summarize_periods(rows,weights)
    stratum_summary=summarize_q_strata(rows)
    write_csv(args.output/"event_window_budgets.csv",rows)
    write_csv(args.output/"time_resolved_summary.csv",time_summary)
    write_csv(args.output/"period_contrasts.csv",period_summary)
    write_csv(args.output/"q_stratum_sensitivity.csv",stratum_summary)
    make_figure(time_summary,args.output/"physical_budget.png")
    closure=np.asarray([r["closure_residual"] for r in rows])
    path=np.asarray([r["path_discretization_correction"] for r in rows])
    result=dict(n_events=len(cohort),n_continued=sum(r["success4"] for r in cohort),
                n_failed=sum(not r["success4"] for r in cohort),
                n_original_seeds=len({r["seed"] for r in cohort}),
                q_edges=Q_EDGES.tolist(),q_standardization_weights=weights.tolist(),
                window_duration=WINDOW,time_range=[START,STOP],
                first_crossing_assignment="segment is in (-0.2,0] bin",
                identity="neck_flux = plug_flux - weighted_accumulation - weighted_lateral_flux + path_discretization_correction",
                maximum_absolute_closure_residual=float(np.max(np.abs(closure))),
                rms_closure_residual=float(np.sqrt(np.mean(closure**2))),
                maximum_absolute_path_correction=float(np.max(np.abs(path))),
                mean_absolute_path_correction=float(np.mean(np.abs(path))))
    (args.output/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),cohort=str(args.cohort),dt=DT,ell=ELL,
                source_direction="fixed from pre-event majority for every window",
                uncertainty="cluster bootstrap over original seeds; seed-equal within outcome and q stratum",
                standardization="direct standardization to pooled cohort frequencies in prior q-a-b |q_pre| strata",
                outcome="at least four additional net fish over five time units")
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":
    main()
