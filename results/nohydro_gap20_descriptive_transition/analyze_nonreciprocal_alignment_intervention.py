#!/usr/bin/env python3
"""Analyze the bounded paired nonreciprocal-alignment intervention."""
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
from scipy.stats import binomtest


METRICS=("turning_alignment_rate_1","turning_alignment_rate_5","geometry_alignment_rate_1","geometry_alignment_rate_5","mean_Ja_1","mean_Ja_5","mean_b_1","mean_b_5","mean_polarization_1","mean_polarization_5","mean_gap_alignment_1","mean_gap_alignment_5","net_transfer","success4","first_entry")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)


def analyze(rows: list[dict],nboot=10000) -> tuple[list[dict],list[dict],dict]:
    paired=defaultdict(dict)
    for r in rows:paired[(int(r["seed"]),int(r["clone"]))][r["arm"]]=r
    if len(paired)!=112 or any(set(v)!={"baseline","symmetric_internal_alignment"} for v in paired.values()):raise RuntimeError("Expected 112 complete paired blocks")
    seed_effects=[]
    for seed in sorted({k[0] for k in paired}):
        blocks=[v for (s,_),v in paired.items() if s==seed];row={"seed":seed,"n_clones":len(blocks)}
        for metric in METRICS:row[metric]=float(np.mean([float(v["symmetric_internal_alignment"][metric])-float(v["baseline"][metric]) for v in blocks]))
        seed_effects.append(row)
    rng=np.random.default_rng(2026112301);effects=[]
    for metric in METRICS:
        values=np.asarray([r[metric] for r in seed_effects]);draw=values[rng.integers(len(values),size=(nboot,len(values)))].mean(axis=1);effects.append(dict(metric=metric,symmetric_minus_baseline=float(values.mean()),ci95_low=float(np.percentile(draw,2.5)),ci95_high=float(np.percentile(draw,97.5)),n_seeds=len(values),n_pairs=len(paired)))
    arm_summary={}
    for arm in ("baseline","symmetric_internal_alignment"):
        arm_summary[arm]={}
        for metric in ("net_transfer","success4","first_entry","mean_Ja_1","mean_b_1","mean_polarization_1","mean_gap_alignment_1","turning_alignment_rate_1"):
            values=np.asarray([np.mean([float(v[arm][metric]) for (s,_),v in paired.items() if s==seed]) for seed in sorted({k[0] for k in paired})]);draw=values[rng.integers(len(values),size=(nboot,len(values)))].mean(axis=1);arm_summary[arm][metric]=dict(mean=float(values.mean()),ci95=[float(np.percentile(draw,2.5)),float(np.percentile(draw,97.5))])
    switches=defaultdict(int)
    for v in paired.values():switches[(int(v["baseline"]["success4"]),int(v["symmetric_internal_alignment"]["success4"]))]+=1
    lost=switches[(1,0)];gained=switches[(0,1)];summary=dict(n_states=14,n_clones=8,n_paired_blocks=112,arm_summary=arm_summary,
        success_switches={"both_failed":switches[(0,0)],"gained":gained,"lost":lost,"both_succeeded":switches[(1,1)],"discordant_exact_p":float(binomtest(lost,lost+gained,.5).pvalue)},effects={r["metric"]:r for r in effects})
    return effects,seed_effects,summary


def figure(effects: list[dict],seed_effects: list[dict],summary: dict,output: Path) -> None:
    lookup={r["metric"]:r for r in effects};fig,axes=plt.subplots(2,2,figsize=(10,7.5));axes=axes.ravel()
    terms=("turning_alignment_rate_1","mean_Ja_1","mean_b_1","mean_gap_alignment_1","mean_polarization_1");labels=("turning","approach current","storage","gap alignment","polarization")
    y=np.asarray([lookup[x]["symmetric_minus_baseline"] for x in terms]);lo=np.asarray([lookup[x]["ci95_low"] for x in terms]);hi=np.asarray([lookup[x]["ci95_high"] for x in terms]);axes[0].bar(range(len(terms)),y,yerr=[y-lo,hi-y],capsize=3);axes[0].axhline(0,color=".5",ls=":");axes[0].set_xticks(range(len(terms)),labels,rotation=35,ha="right");axes[0].set(title="First-unit paired effects",ylabel="symmetric - baseline")
    arms=("baseline","symmetric_internal_alignment");labels_arm=("baseline","symmetric")
    y=np.asarray([summary["arm_summary"][x]["net_transfer"]["mean"] for x in arms]);lo=np.asarray([summary["arm_summary"][x]["net_transfer"]["ci95"][0] for x in arms]);hi=np.asarray([summary["arm_summary"][x]["net_transfer"]["ci95"][1] for x in arms]);axes[1].bar(labels_arm,y,yerr=[y-lo,hi-y],capsize=3);axes[1].set(title="Five-unit net transfer",ylabel="net fish toward target")
    width=.35;x=np.arange(2)
    for j,arm in enumerate(arms):axes[2].bar(x+(j-.5)*width,[summary["arm_summary"][arm][m]["mean"] for m in ("success4","first_entry")],width,label=labels_arm[j])
    axes[2].set_xticks(x,("P(>=4 net)","P(first entry)"));axes[2].set(ylim=(0,1),title="Threshold outcomes");axes[2].legend(frameon=False)
    vals=np.asarray([r["net_transfer"] for r in seed_effects]);axes[3].bar(range(len(vals)),vals);axes[3].axhline(0,color=".5",ls=":");axes[3].set(xlabel="original seed",ylabel="net-transfer effect",title="Seed-level paired transfer")
    fig.tight_layout();fig.savefig(output,dpi=220);plt.close(fig)


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--branches",type=Path,required=True);p.add_argument("--audit",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    effects,seeds,summary=analyze(read_csv(a.branches));instant=json.loads(a.audit.read_text())["instantaneous_counterfactual"];instant={k:v for k,v in instant.items() if k!="rows"};instant["fixed_cohort_lead1"]={k:v for k,v in instant["fixed_cohort_lead1"].items() if k!="rows"};summary["instantaneous_audit"]=instant;write_csv(a.output/"paired_effects.csv",effects);write_csv(a.output/"per_seed_effects.csv",seeds);figure(effects,seeds,summary,a.output/"causal_intervention_effects.png");(a.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":main()
