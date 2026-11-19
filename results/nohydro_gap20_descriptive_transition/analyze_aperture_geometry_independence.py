#!/usr/bin/env python3
"""Isolate aperture geometry from removal of negative approach current.

Uses the existing prospective cohorts and held-out split.  No outcome-fitted
predictor is constructed: AUCs are direct ranks, optionally restricted to
pre-existing (q,a,b) cells and outcome-blind inward-supply strata.
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

from analyze_aperture_capture_current import capture_geometry, fit_b_cuts, qab_cell
from analyze_approach_current_balance import raw_path
from analyze_initiation_propagation import RawReader, auc_score


METRICS=("Ja","Jin","Jcap")
CONDITIONS=("unconditional","q_a_b","q_a_b_and_Jin")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def enrich_jin(raw_dir: Path,source: list[dict]) -> list[dict]:
    by_seed=defaultdict(list)
    for row in source:by_seed[int(row["seed"])].append(row)
    out=[]
    for seed in sorted(by_seed):
        reader=RawReader(raw_path(raw_dir,seed))
        try:
            for row in by_seed[seed]:
                step=int(row["step"]);state=reader.frame("state",step)
                cap=capture_geometry(state,int(row["source_sign"]))
                cached=float(row["Jcap"])
                if abs(cap["j_capture"]-cached)>1e-12:
                    raise RuntimeError(f"Jcap parity failure seed={seed} step={step}: {cap['j_capture']} vs {cached}")
                out.append(dict(seed=seed,step=step,outcome=int(row["outcome"]),
                                abs_q=float(row["abs_q"]),a=float(row["a"]),b=float(row["b"]),
                                Ja=float(row["Ja"]),Jin=cap["j_inward"],Jcap=cached))
        finally:reader.close()
    return out


def fit_supply_strata(train: list[dict],bcuts: dict,quantiles=(.25,.5,.75)) -> tuple[dict,list[dict]]:
    groups=defaultdict(list)
    for row in train:
        cell=qab_cell(row,bcuts)
        if cell is not None:groups[cell].append(float(row["Jin"]))
    cuts={};audit=[]
    for cell,values in sorted(groups.items()):
        values=np.asarray(values);q=np.quantile(values,quantiles)
        cuts[cell]=q
        audit.append(dict(q_bin=cell[0],a_bin=cell[1],b_bin=cell[2],train_n=len(values),
                          jin_q25=float(q[0]),jin_q50=float(q[1]),jin_q75=float(q[2]),
                          train_zero_fraction=float(np.mean(values==0))))
    return cuts,audit


def cell_label(row: dict,condition: str,bcuts: dict,supply_cuts: dict) -> tuple|str|None:
    if condition=="unconditional":return "all"
    base=qab_cell(row,bcuts)
    if base is None:return None
    if condition=="q_a_b":return base
    if base not in supply_cuts:return None
    supply_bin=int(np.searchsorted(supply_cuts[base],float(row["Jin"]),side="right"))
    return base+(supply_bin,)


def conditional_auc(rows: list[dict],metric: str,condition: str,bcuts: dict,supply_cuts: dict) -> tuple[float,int,int]:
    groups=defaultdict(list)
    for row in rows:
        cell=cell_label(row,condition,bcuts,supply_cuts)
        if cell is not None:groups[cell].append(row)
    concordant=0.0;pairs=0;supported=0
    for part in groups.values():
        pos=np.asarray([float(r[metric]) for r in part if r["outcome"]==1]);neg=np.asarray([float(r[metric]) for r in part if r["outcome"]==0])
        if not len(pos) or not len(neg):continue
        diff=pos[:,None]-neg[None,:]
        concordant+=float(np.sum(diff>0)+.5*np.sum(diff==0));pairs+=diff.size;supported+=1
    return (concordant/pairs if pairs else float("nan")),pairs,supported


def bootstrap_comparison(rows: list[dict],condition: str,bcuts: dict,supply_cuts: dict,nboot: int=2000,seed: int=20261118) -> tuple[list[dict],list[dict]]:
    seeds=sorted({r["seed"] for r in rows});by={s:[r for r in rows if r["seed"]==s] for s in seeds}
    observed={m:conditional_auc(rows,m,condition,bcuts,supply_cuts)[0] for m in METRICS}
    _,pairs,cells=conditional_auc(rows,"Jcap",condition,bcuts,supply_cuts)
    rng=np.random.default_rng(seed);draws={m:[] for m in METRICS};diffs={"Jin_minus_Ja":[],"Jcap_minus_Jin":[]}
    for _ in range(nboot):
        sample=[]
        for s in rng.choice(seeds,len(seeds),replace=True):sample.extend(by[int(s)])
        values={m:conditional_auc(sample,m,condition,bcuts,supply_cuts)[0] for m in METRICS}
        if all(np.isfinite(list(values.values()))):
            for m in METRICS:draws[m].append(values[m])
            diffs["Jin_minus_Ja"].append(values["Jin"]-values["Ja"])
            diffs["Jcap_minus_Jin"].append(values["Jcap"]-values["Jin"])
    metric_rows=[dict(condition=condition,metric=m,auc=observed[m],ci95=np.percentile(draws[m],[2.5,97.5]).tolist(),
                      heldout_n=len(rows),heldout_seeds=len(seeds),concordant_pairs=pairs,supported_cells=cells) for m in METRICS]
    contrast_rows=[]
    for name,(first,second) in {"Jin_minus_Ja":("Jin","Ja"),"Jcap_minus_Jin":("Jcap","Jin")}.items():
        contrast_rows.append(dict(condition=condition,contrast=name,difference=observed[first]-observed[second],
                                  ci95=np.percentile(diffs[name],[2.5,97.5]).tolist(),bootstrap_draws=len(diffs[name])))
    return metric_rows,contrast_rows


def supply_balance(rows: list[dict],bcuts: dict,supply_cuts: dict) -> dict:
    """Outcome difference in Jin after joint-cell restriction; diagnostic only."""
    groups=defaultdict(list)
    for row in rows:
        cell=cell_label(row,"q_a_b_and_Jin",bcuts,supply_cuts)
        if cell is not None:groups[cell].append(row)
    weighted=0.0;total=0
    for part in groups.values():
        pos=[r["Jin"] for r in part if r["outcome"]];neg=[r["Jin"] for r in part if not r["outcome"]]
        if pos and neg:
            pairs=len(pos)*len(neg);weighted+=pairs*(float(np.mean(pos))-float(np.mean(neg)));total+=pairs
    return dict(pair_weighted_positive_minus_negative_Jin=weighted/total if total else None,concordant_pairs=total)


def make_figure(rows: list[dict],output: Path) -> None:
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),sharey=True);colors={"Ja":"#888888","Jin":"#E69F00","Jcap":"#0072B2"}
    for ax,(cohort,title) in zip(axes,(('initiation','future first entry'),('transfer','four-fish transfer'))):
        part=[r for r in rows if r["cohort"]==cohort];x=np.arange(len(CONDITIONS));width=.24
        for j,metric in enumerate(METRICS):
            data=[next(r for r in part if r["condition"]==c and r["metric"]==metric) for c in CONDITIONS]
            y=np.asarray([r["auc"] for r in data]);ci=np.asarray([r["ci95"] for r in data]).T
            ax.bar(x+(j-1)*width,y,width,label=metric,color=colors[metric],yerr=np.vstack([y-ci[0],ci[1]-y]),capsize=2,error_kw={"lw":.8})
        ax.axhline(.5,color=".35",ls=":");ax.set_xticks(x,["all","same q,a,b","same q,a,b\nand Jin stratum"]);ax.set(title=title,ylabel="held-out conditional AUC",ylim=(.45,.82));ax.legend(frameon=False)
    fig.suptitle("Does aperture geometry add information beyond positive inward supply?")
    fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    rows=[]
    for seed in range(4):
        for outcome in (0,1):
            rows.append(dict(seed=seed,outcome=outcome,abs_q=.4,a=0.,b=.1,Ja=outcome,Jin=outcome,Jcap=outcome))
    bcuts=fit_b_cuts(rows);scuts,audit=fit_supply_strata(rows,bcuts)
    assert audit and conditional_auc(rows,"Jcap","unconditional",bcuts,scuts)[0]==1
    print("PASS: Jin strata and paired conditional concordance")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--initiation",type=Path);p.add_argument("--transfer",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.initiation,args.transfer,args.output)):p.error("all paths required")
    args.output.mkdir(parents=True,exist_ok=True)
    all_metrics=[];all_contrasts=[];all_cuts=[];heldout_rows=[];summary={}
    for cohort,path in (("initiation",args.initiation),("transfer",args.transfer)):
        rows=enrich_jin(args.raw_dir,read_csv(path));train=[r for r in rows if r["seed"]<30];test=[r for r in rows if r["seed"]>=30]
        bcuts=fit_b_cuts(train);scuts,cut_audit=fit_supply_strata(train,bcuts)
        for r in cut_audit:r["cohort"]=cohort
        all_cuts.extend(cut_audit)
        for condition in CONDITIONS:
            metrics,contrasts=bootstrap_comparison(test,condition,bcuts,scuts,seed=20261118+len(all_metrics))
            for r in metrics:r["cohort"]=cohort
            for r in contrasts:r["cohort"]=cohort
            all_metrics.extend(metrics);all_contrasts.extend(contrasts)
        for r in test:
            r["cohort"]=cohort;r["qab_cell"]=str(qab_cell(r,bcuts));r["supply_cell"]=str(cell_label(r,"q_a_b_and_Jin",bcuts,scuts))
        heldout_rows.extend(test)
        summary[cohort]=dict(n_all=len(rows),n_train=len(train),n_heldout=len(test),heldout_positives=sum(r["outcome"] for r in test),
                             supply_balance=supply_balance(test,bcuts,scuts))
    write_csv(args.output/"heldout_auc.csv",all_metrics);write_csv(args.output/"paired_auc_contrasts.csv",all_contrasts)
    write_csv(args.output/"training_supply_strata.csv",all_cuts);write_csv(args.output/"heldout_validation_rows.csv",heldout_rows)
    make_figure(all_metrics,args.output/"geometry_independence.png")
    summary["conditions"]={"unconditional":"direct heldout ranks","q_a_b":"pairs only within prior q-a-b cells","q_a_b_and_Jin":"pairs only within prior q-a-b cells and training-only Jin quartiles"}
    summary["metrics"]={"Ja":"signed stored-velocity approach current","Jin":"positive inward heading current over all approach fish","Jcap":"Jin restricted by exact straight-ray aperture condition"}
    summary["auc"]=all_metrics;summary["contrasts"]=all_contrasts
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),initiation_input=str(args.initiation),transfer_input=str(args.transfer),
                train_seeds="0--29",heldout_seeds="30--39",inward_supply_strata="within-q-a-b training quartiles",
                supply_quantiles=[.25,.5,.75],bootstrap_unit="original trajectory seed",bootstrap_draws=2000,
                fitted_predictors=0,new_simulations=0)
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
