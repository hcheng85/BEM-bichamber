#!/usr/bin/env python3
"""Seed-cluster summaries and figures for the fixed-leading-fish analysis."""
from __future__ import annotations
import csv, gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(sys.argv[1]).resolve()
TAB=ROOT/"tables"; FIG=ROOT/"figures"; FIG.mkdir(exist_ok=True)
RNG=np.random.default_rng(20261022)

def ci(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    b=np.mean(v[RNG.integers(0,len(v),(30000,len(v)))],axis=1)
    return float(np.mean(v)),float(np.percentile(b,2.5)),float(np.percentile(b,97.5))

def write(name,rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields: fields.append(k)
    with (TAB/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

f=pd.read_csv(TAB/"fixed_leading_fish.csv")
h=pd.read_csv(TAB/"tracked_fish_history.csv.gz")
for suffix in ("full","to_entry"):
    for k in ("wall","alignment","attraction","direct_hydro","boundary_hydro_residual"):
        f[f"{k}_closing_{suffix}"]=-f.initial_eta_sign*f[f"{k}_impulse_{suffix}"]

metrics=[]
def paired_metric(name,col,subset=None):
    q=f if subset is None else f[subset]
    vals=[]; event=[]; control=[]
    for seed in range(20):
        a=q[(q.seed==seed)&(q.condition=="event")][col].mean()
        b=q[(q.seed==seed)&(q.condition=="control")][col].mean()
        vals.append(a-b);event.append(a);control.append(b)
    m,lo,hi=ci(vals)
    metrics.append(dict(metric=name,event_mean=float(np.nanmean(event)),control_mean=float(np.nanmean(control)),
        event_minus_control=m,seed_ci_low=lo,seed_ci_high=hi,independent_seeds=20))
    return np.asarray(vals)

paired_metric("initial_abs_angular_error_rad","initial_abs_eta")
paired_metric("gate_entry_probability","gate_entered")
paired_metric("spatial_crossing_probability_fixed_ten","reached_divider_through_opening")
valid=f.rng_valid==1
paired_metric("deterministic_closing_impulse_full","deterministic_closing_full",valid)
paired_metric("stochastic_closing_increment_full","stochastic_closing_full",valid)
for k in ("wall","alignment","attraction","direct_hydro","boundary_hydro_residual"):
    paired_metric(f"{k}_closing_impulse_full",f"{k}_closing_full",valid)
f["combined_hydro_closing_full"]=f.direct_hydro_closing_full+f.boundary_hydro_residual_closing_full
paired_metric("combined_hydro_closing_impulse_full","combined_hydro_closing_full",valid)

# Coverage and turnover, clustered by seed.
e=f[f.condition=="event"].groupby(["seed","pair_id"]).event_first_crosser_covered.first().reset_index()
coverage=[e[e.seed==s].event_first_crosser_covered.mean() for s in range(20)]
cov_mean,cov_lo,cov_hi=ci(coverage); metrics.append(dict(metric="first_crosser_coverage_by_fixed_ten",event_mean=cov_mean,control_mean="",
    event_minus_control="",seed_ci_low=cov_lo,seed_ci_high=cov_hi,independent_seeds=20))
for t in (1.0,2.0):
    q=h[np.isclose(h.time_from_A_minus_2,t)].drop_duplicates(["seed","pair_id","condition"])
    vals=[];ev=[];co=[]
    for s in range(20):
        a=q[(q.seed==s)&(q.condition=="event")].leading_retention.mean(); b=q[(q.seed==s)&(q.condition=="control")].leading_retention.mean()
        vals.append(a-b);ev.append(a);co.append(b)
    mm,ll,hh=ci(vals); metrics.append(dict(metric=f"leading_ten_retention_t{t:g}",event_mean=np.mean(ev),control_mean=np.mean(co),
        event_minus_control=mm,seed_ci_low=ll,seed_ci_high=hh,independent_seeds=20))
write("mechanism_summary.csv",metrics)
with (ROOT/"summary.csv").open("w",newline="") as fh:
    fields=[]
    for r in metrics:
        for k in r:
            if k not in fields: fields.append(k)
    ww=csv.DictWriter(fh,fieldnames=fields,lineterminator="\n");ww.writeheader();ww.writerows(metrics)

# Entry budgets use fish that begin outside and enter, but aggregate within seed.
entr=f[(f.rng_valid==1)&(f.initially_admissible==0)&(f.gate_entered==1)].copy()
budget_cols=("deterministic_closing_to_entry","stochastic_closing_to_entry","geometry_closing_to_entry")
budget=[]
for cond in ("event","control"):
    for col in budget_cols:
        vals=[entr[(entr.seed==s)&(entr.condition==cond)][col].mean() for s in range(20)]
        m,lo,hi=ci(vals);budget.append(dict(condition=cond,component=col,mean=m,seed_ci_low=lo,seed_ci_high=hi,
            independent_seeds=int(np.sum(np.isfinite(vals)))) )
write("gate_entry_budget.csv",budget)

source_budget=[]
for cond in ("event","control"):
    for k in ("wall","alignment","attraction","direct_hydro","boundary_hydro_residual"):
        col=f"{k}_closing_to_entry"; vals=[entr[(entr.seed==s)&(entr.condition==cond)][col].mean() for s in range(20)]
        m,lo,hi=ci(vals);source_budget.append(dict(condition=cond,component=k,mean=m,seed_ci_low=lo,seed_ci_high=hi,independent_seeds=20))
write("gate_entry_source_budget.csv",source_budget)

# Eventual first crosser versus the same A-2 rank in its matched control.
cross=f[f.is_eventual_first_crosser==1].copy()
control=f[f.condition=="control"].set_index(["seed","pair_id","leading_rank"])
matched=[]
for metric in ("initial_abs_eta","initially_admissible","gate_entered",
               "deterministic_closing_full","stochastic_closing_full"):
    vals=[]; neval=0
    for seed in range(20):
        local=[]
        for _,r in cross[cross.seed==seed].iterrows():
            c=control.loc[(r.seed,r.pair_id,r.leading_rank)]
            if metric in ("deterministic_closing_full","stochastic_closing_full") and not (r.rng_valid and c.rng_valid): continue
            local.append(float(r[metric])-float(c[metric])); neval+=1
        vals.append(np.mean(local) if local else np.nan)
    m,lo,hi=ci(vals);matched.append(dict(metric=metric,event_minus_same_rank_control=m,
        seed_ci_low=lo,seed_ci_high=hi,independent_seeds=int(np.sum(np.isfinite(vals))),fish_pairs=neval))
write("first_crosser_matched_rank_contrasts.csv",matched)

# Time histories: average fish within pair, pairs within seed, then seeds.
time_rows=[]
for t in sorted(h.time_from_A_minus_2.unique()):
    for cond in ("event","control"):
        q=h[(h.condition==cond)&np.isclose(h.time_from_A_minus_2,t)]
        for metric,col in (("abs_angular_error","signed_angular_error"),("admissible","admissible"),("leading_retention","leading_retention")):
            qq=q.copy(); qq["value"]=np.abs(qq[col]) if metric=="abs_angular_error" else qq[col]
            pair=qq.groupby(["seed","pair_id"]).value.mean().reset_index()
            vals=[pair[pair.seed==s].value.mean() for s in range(20)]
            m,lo,hi=ci(vals);time_rows.append(dict(time=t,condition=cond,metric=metric,mean=m,seed_ci_low=lo,seed_ci_high=hi))
write("time_resolved_summary.csv",time_rows)

plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
fig,axs=plt.subplots(2,2,figsize=(9,6.5),constrained_layout=True)
labels={"event":"imminent first crossing","control":"matched no crossing"}; colors={"event":"#b2182b","control":"#2166ac"}
for ax,metric,ylabel in ((axs[0,0],"abs_angular_error",r"mean $|\eta_i|$ (rad)"),(axs[0,1],"admissible","fraction inside angular gate"),(axs[1,0],"leading_retention","fraction of fixed ten still leading")):
    for cond in ("event","control"):
        q=pd.DataFrame(time_rows);q=q[(q.metric==metric)&(q.condition==cond)]
        ax.plot(q.time-2,q["mean"],color=colors[cond],label=labels[cond]);ax.fill_between(q.time-2,q.seed_ci_low,q.seed_ci_high,color=colors[cond],alpha=.16)
    ax.set_xlabel("time relative to A");ax.set_ylabel(ylabel);ax.axvline(0,color="k",lw=.7,ls="--")
axs[0,0].legend(frameon=False)
ranks=f[(f.condition=="event")&(f.is_eventual_first_crosser==1)].leading_rank
axs[1,1].hist(ranks,bins=np.arange(.5,11.6,1),color="#762a83",alpha=.8)
axs[1,1].set(xlabel="A-2 rank of eventual first crosser",ylabel="events captured")
axs[1,1].text(.98,.95,f"coverage = {cov_mean:.3f}\n95% seed CI [{cov_lo:.3f}, {cov_hi:.3f}]",transform=axs[1,1].transAxes,ha="right",va="top")
for ext in ("png","pdf"): fig.savefig(FIG/f"01_fixed_identity_time_course.{ext}",dpi=220)
plt.close(fig)

# Source-resolved and gate-entry budgets.
fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
source_names=["wall","alignment","attraction","combined hydro","total deterministic","stochastic"]
source_metrics=["wall_closing_impulse_full","alignment_closing_impulse_full","attraction_closing_impulse_full",
                "combined_hydro_closing_impulse_full","deterministic_closing_impulse_full","stochastic_closing_increment_full"]
mm=pd.DataFrame(metrics).set_index("metric").loc[source_metrics]
x=np.arange(len(mm));axs[0].errorbar(x,mm.event_minus_control,yerr=[mm.event_minus_control-mm.seed_ci_low,mm.seed_ci_high-mm.event_minus_control],fmt="o",color="#762a83",capsize=3)
axs[0].axhline(0,color="k",lw=.7);axs[0].set_xticks(x,source_names,rotation=35,ha="right");axs[0].set_ylabel("event - control closing impulse (rad)");axs[0].set_title("Net contrasts; seed 95% CI")
hydro_names=["direct hydro","BEM + override\nresidual"]
hydro_metrics=["direct_hydro_closing_impulse_full","boundary_hydro_residual_closing_impulse_full"]
hm=pd.DataFrame(metrics).set_index("metric").loc[hydro_metrics];x=np.arange(2)
axs[1].errorbar(x,hm.event_minus_control,yerr=[hm.event_minus_control-hm.seed_ci_low,hm.seed_ci_high-hm.event_minus_control],fmt="o",color="#762a83",capsize=3)
axs[1].axhline(0,color="k",lw=.7);axs[1].set_xticks(x,hydro_names,rotation=20,ha="right");axs[1].set_ylabel("event - control impulse (rad)");axs[1].set_title("Large hydrodynamic cancellation")
bb=pd.DataFrame(budget); comps=list(budget_cols); x=np.arange(3);w=.36
for j,cond in enumerate(("event","control")):
    q=bb[bb.condition==cond].set_index("component").loc[comps]
    axs[2].bar(x+(j-.5)*w,q["mean"],w,label=labels[cond],color=colors[cond],alpha=.8)
axs[2].axhline(0,color="k",lw=.7);axs[2].set_xticks(x,["deterministic torque","angular noise","moving corridor"],rotation=20,ha="right")
axs[2].set_ylabel("closing contribution before first gate entry (rad)");axs[2].set_title("Fish initially outside gate");axs[2].legend(frameon=False)
for ext in ("png","pdf"): fig.savefig(FIG/f"02_angular_budget_and_controls.{ext}",dpi=220)
plt.close(fig)
