#!/usr/bin/env python3
"""Recompute periodic angular metrics, paired contrasts, and final figures."""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAGS=(-2.0,-1.0,-0.2)
CONDS=("success","failure","control")
COLORS={"success":"#009E73","failure":"#D55E00","control":"#777777"}
TORQUES=("wall","alignment","attraction","direct_hydro","boundary_hydro")

def read(path): return list(csv.DictReader(path.open()))
def write(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)

def finite_mean(values):
    v=np.asarray(values,float); return float(np.nanmean(v)) if np.any(np.isfinite(v)) else np.nan

def main():
    ap=argparse.ArgumentParser();ap.add_argument("analysis_dir",type=Path);a=ap.parse_args();root=a.analysis_dir
    rows=read(root/"tables"/"leading_fish_torques.csv")
    for r in rows:
        psi=float(r["local_heading"]);d=float(r["distance_to_midline"]);y=float(r["y"])
        gap=7.158910531638177;tip_r=.5
        cy_lo=-gap/2-tip_r;cy_hi=gap/2+tip_r
        lo=np.arctan2(cy_lo-y,d)+np.arcsin(np.clip(tip_r/np.hypot(d,cy_lo-y),0,1))
        hi=np.arctan2(cy_hi-y,d)-np.arcsin(np.clip(tip_r/np.hypot(d,cy_hi-y),0,1))
        r["admissible_low"]=lo;r["admissible_high"]=hi
        if lo<=psi<=hi: err=0.0
        else:
            dlo=np.arctan2(np.sin(psi-lo),np.cos(psi-lo));dhi=np.arctan2(np.sin(psi-hi),np.cos(psi-hi))
            err=dlo if abs(dlo)<=abs(dhi) else dhi
        r["signed_angular_error"]=err;r["abs_angular_error"]=abs(err);r["admissible"]=int(err==0)
        s=1.0 if r["origin"]=="left" else -1.0;sg=0 if err==0 else np.sign(err)
        for name in TORQUES:
            r[f"{name}_closing_rate"]=np.nan if not sg else -sg*s*float(r[f"{name}_torque"])
        r["deterministic_closing_rate"]=np.nan if not sg else -sg*s*float(r["deterministic_total"])
        r["noise_closing_increment"]=np.nan if not sg else -sg*s*float(r["inferred_noise_increment"])
    write(root/"tables"/"leading_fish_torques.csv",rows)
    metrics=["abs_angular_error","admissible"]+[f"{x}_closing_rate" for x in TORQUES]+["deterministic_closing_rate","noise_closing_increment"]
    seed_rows=[]
    for seed in range(15,20):
      for cond in CONDS:
       for lag in LAGS:
        q=[r for r in rows if int(r["seed"])==seed and r["condition"]==cond and float(r["lag"])==lag]
        out=dict(seed=seed,condition=cond,lag=lag,n_leading=len(q),outside_fish=sum(int(r["admissible"])==0 for r in q))
        for m in metrics: out[m]=finite_mean([r[m] for r in q])
        seed_rows.append(out)
    write(root/"tables"/"seed_level_summary.csv",seed_rows)
    rng=np.random.default_rng(20261021);group=[]
    for cond in CONDS:
      for lag in LAGS:
       q=[r for r in seed_rows if r["condition"]==cond and r["lag"]==lag]
       for m in metrics:
        v=np.array([r[m] for r in q],float);v=v[np.isfinite(v)]
        if len(v):
            boot=np.mean(v[rng.integers(0,len(v),(10000,len(v)))],axis=1)
            group.append(dict(condition=cond,lag=lag,metric=m,mean=np.mean(v),median=np.median(v),
                seed_bootstrap_ci_low=np.percentile(boot,2.5),seed_bootstrap_ci_high=np.percentile(boot,97.5),independent_seeds=len(v)))
    write(root/"tables"/"group_summary.csv",group)
    contrasts=[]
    for lag in LAGS:
      for ref in ("failure","control"):
       for m in metrics:
        diffs=[]
        for seed in range(15,20):
            x=next(r for r in seed_rows if r["seed"]==seed and r["condition"]=="success" and r["lag"]==lag)
            y=next(r for r in seed_rows if r["seed"]==seed and r["condition"]==ref and r["lag"]==lag)
            d=float(x[m])-float(y[m])
            if np.isfinite(d): diffs.append(d)
        v=np.asarray(diffs);boot=np.mean(v[rng.integers(0,len(v),(20000,len(v)))],axis=1)
        contrasts.append(dict(lag=lag,metric=m,contrast=f"success_minus_{ref}",mean_difference=np.mean(v),
            median_difference=np.median(v),seed_bootstrap_ci_low=np.percentile(boot,2.5),
            seed_bootstrap_ci_high=np.percentile(boot,97.5),positive_seeds=np.sum(v>0),independent_seeds=len(v)))
    write(root/"tables"/"paired_contrasts.csv",contrasts)

    fig,axes=plt.subplots(1,3,figsize=(13,4.2),constrained_layout=True)
    for cond in CONDS:
      q=[r for r in seed_rows if r["condition"]==cond]
      for key,ax in (("abs_angular_error",axes[0]),("admissible",axes[1]),("deterministic_closing_rate",axes[2])):
        means=[finite_mean([r[key] for r in q if r["lag"]==lag]) for lag in LAGS]
        ax.plot([-x for x in LAGS],means,marker="o",label=cond,color=COLORS[cond])
    for ax in axes:ax.grid(alpha=.2);ax.set_xlabel("time before anchor A");ax.invert_xaxis()
    axes[0].set_ylabel("shortest angular distance to opening (rad)");axes[1].set_ylabel("fraction geometrically admissible")
    axes[2].axhline(0,color="k",lw=.8);axes[2].set_ylabel("deterministic rate toward opening (rad/time)");axes[0].legend(frameon=False)
    for ext in ("png","pdf"):fig.savefig(root/"figures"/f"01_pre_event_angular_dynamics.{ext}",dpi=220)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True);xs=np.arange(3);width=.05
    for j,name in enumerate(TORQUES):
        vals=[finite_mean([r[f"{name}_closing_rate"] for r in seed_rows if r["condition"]=="success" and r["lag"]==lag]) for lag in LAGS]
        axes[0].bar(xs+(j-2)*width,vals,width,label=name.replace("_"," "))
    axes[0].axhline(0,color="k",lw=.8);axes[0].set_xticks(xs,[str(x) for x in LAGS]);axes[0].set_xlabel("lag from A")
    axes[0].set_ylabel("success: rate toward opening (rad/time)");axes[0].legend(frameon=False,fontsize=8)
    for cond in CONDS:
        q=[r for r in rows if r["condition"]==cond and int(r["admissible"])==0]
        xx=np.array([float(r["signed_angular_error"]) for r in q]);yy=np.array([float(r["wall_closing_rate"]) for r in q])
        axes[1].scatter(xx,yy,s=12,alpha=.22,color=COLORS[cond])
        edges=np.linspace(-np.pi,np.pi,13);centers=(edges[:-1]+edges[1:])/2
        med=np.array([np.median(yy[(xx>=a)&(xx<b)]) if np.sum((xx>=a)&(xx<b))>=4 else np.nan
                      for a,b in zip(edges[:-1],edges[1:])])
        axes[1].plot(centers,med,marker="o",ms=3,lw=1.8,color=COLORS[cond],label=cond)
    axes[1].axhline(0,color="k",lw=.8);axes[1].axvline(0,color="k",lw=.5);axes[1].grid(alpha=.2)
    axes[1].set_xlabel("shortest signed angular error (rad)");axes[1].set_ylabel("wall rate toward opening");axes[1].legend(frameon=False,fontsize=8)
    for ext in ("png","pdf"):fig.savefig(root/"figures"/f"02_torque_decomposition_and_restoring.{ext}",dpi=220)
    plt.close(fig)
    print(f"updated {len(rows)} fish rows, {len(seed_rows)} seed summaries, {len(contrasts)} contrasts")
if __name__=="__main__":main()
