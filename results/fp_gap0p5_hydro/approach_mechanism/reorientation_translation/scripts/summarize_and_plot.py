#!/usr/bin/env python3
"""Create compact summary and figures from validation tables."""
from __future__ import annotations
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read(path):
    with path.open() as f:return list(csv.DictReader(f))


def write(path,rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("analysis_dir",type=Path);args=ap.parse_args()
    out=args.analysis_dir.resolve();tab=out/"tables";(out/"figures").mkdir(exist_ok=True)
    geom=read(tab/"geometric_decomposition_summary.csv")
    models=read(tab/"coordinate_model_pooled.csv")
    kin=read(tab/"segment_kinematics_group.csv");kc=read(tab/"segment_kinematics_contrasts.csv")
    risk=read(tab/"distance_matched_risk_set.csv");riskm=read(tab/"distance_matched_risk_metrics.csv")

    fig,axes=plt.subplots(1,2,figsize=(10,4.2),constrained_layout=True)
    components=[("delta_log_Lambda",r"$\Delta\log\Lambda$","#0072B2"),
                ("delta_log_sigma_perp",r"$\Delta\log\sigma_\perp$","#009E73"),
                ("minus_delta_log_distance",r"$-\Delta\log d$","#D55E00")]
    for ax,lead in zip(axes,(-2.,-1.)):
        q=[r for r in geom if float(r["lead"])==lead and r["component"]!="distance_fraction_of_mean_delta_log_Lambda"]
        vals=[];lo=[];hi=[]
        for key,_,_ in components:
            r=next(x for x in q if x["component"]==key);v=float(r["mean"]);vals.append(v)
            lo.append(v-float(r["seed_bootstrap_ci_low"]));hi.append(float(r["seed_bootstrap_ci_high"])-v)
        ax.bar(range(3),vals,color=[x[2] for x in components],yerr=np.array([lo,hi]),capsize=4)
        ax.set_xticks(range(3),[x[1] for x in components]);ax.axhline(0,color="k",lw=.8);ax.set_title(f"lead {abs(lead):g}")
        ax.set_ylabel("paired event - control")
    for ext in ("png","pdf"):fig.savefig(out/"figures"/f"01_Lambda_decomposition.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.4),constrained_layout=True)
    order=["psi","psi_Lambda","psi_distance","psi_sigma_perp","psi_distance_sigma_perp"]
    labels=[r"$\psi$",r"$\psi,\Lambda$",r"$\psi,d$",r"$\psi,\sigma_\perp$",r"$\psi,d,\sigma_\perp$"]
    for ax,lead in zip(axes,(-2.,-1.)):
        q={r["model"]:float(r["auc"]) for r in models if float(r["lead"])==lead}
        ax.bar(range(len(order)),[q[k] for k in order],color=["#777777","#56B4E9","#D55E00","#009E73","#CC79A7"])
        ax.set_xticks(range(len(order)),labels,rotation=25);ax.axhline(.5,color="k",lw=.8);ax.set_ylim(.35,.9)
        ax.set_ylabel("held-out pooled AUC");ax.set_title(f"lead {abs(lead):g}")
    for ext in ("png","pdf"):fig.savefig(out/"figures"/f"02_coordinate_comparison.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
    keys=["actual_delta_distance","self_propulsion_delta_distance","direct_hydro_delta_distance","boundary_hydro_effective_delta_distance"]
    labels=["actual","self propulsion","direct hydro","effective BEM"]
    colors=["#333333","#0072B2","#D55E00","#009E73"]
    for cond,x0 in (("event",0),("control",1)):
        q={r["metric"]:float(r["mean_of_seed_means"]) for r in kin if r["condition"]==cond}
        axes[0].bar(np.arange(4)+(x0-.5)*.35,[q[k] for k in keys],.35,label=cond,color=colors,alpha=1 if cond=="event" else .45)
    axes[0].set_xticks(np.arange(4),labels,rotation=20);axes[0].axhline(0,color="k",lw=.8);axes[0].set_ylabel(r"change from $A-2$ to $A-1$")
    axes[0].legend(frameon=False)
    q={r["metric"]:float(r["event_minus_control"]) for r in kc}
    vals=[q[k] for k in keys]+[q["actual_delta_sigma_perp"]]
    axes[1].bar(range(5),vals,color=colors+["#CC79A7"]);axes[1].set_xticks(range(5),labels+[r"$\sigma_\perp$ change"],rotation=20)
    axes[1].axhline(0,color="k",lw=.8);axes[1].set_ylabel("event - control");axes[1].set_title("paired kinematic contrast")
    for ext in ("png","pdf"):fig.savefig(out/"figures"/f"03_inward_kinematics.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10,4.2),constrained_layout=True)
    for cond,color in (("event","#0072B2"),("control","#999999")):
        q=[r for r in risk if r["condition"]==cond]
        axes[0].hist([float(r["distance"]) for r in q],bins=25,density=True,histtype="step",lw=2,color=color,label=cond)
        axes[1].hist([float(r["cos_psi"]) for r in q],bins=25,density=True,histtype="step",lw=2,color=color,label=cond)
    axes[0].set_xlabel(r"distance $d$ at risk time");axes[1].set_xlabel(r"$\cos\psi$ at risk time")
    for ax in axes:ax.set_ylabel("density");ax.legend(frameon=False)
    auc=float(next(r for r in riskm if r["sample"]=="pooled")["auc"]);axes[1].set_title(f"heading-only held-out AUC = {auc:.3f}")
    for ext in ("png","pdf"):fig.savefig(out/"figures"/f"04_distance_matched_risk_set.{ext}",dpi=220)
    plt.close(fig)

    gs={r["component"]:r for r in geom if float(r["lead"])==-1}
    pm={r["model"]:r for r in models if float(r["lead"])==-1}
    cg={r["metric"]:r for r in kc}
    rp=read(tab/"distance_matched_paired_contrasts.csv")
    rs=read(tab/"distance_matched_seed_summary.csv")[0]
    summary=[]
    def add(k,v,n=""):summary.append(dict(metric=k,value=v,note=n))
    add("lead1_delta_log_Lambda",gs["delta_log_Lambda"]["mean"],f"95% seed CI {gs['delta_log_Lambda']['seed_bootstrap_ci_low']} to {gs['delta_log_Lambda']['seed_bootstrap_ci_high']}")
    add("lead1_delta_log_sigma_perp",gs["delta_log_sigma_perp"]["mean"],f"95% seed CI {gs['delta_log_sigma_perp']['seed_bootstrap_ci_low']} to {gs['delta_log_sigma_perp']['seed_bootstrap_ci_high']}")
    add("lead1_minus_delta_log_distance",gs["minus_delta_log_distance"]["mean"],f"95% seed CI {gs['minus_delta_log_distance']['seed_bootstrap_ci_low']} to {gs['minus_delta_log_distance']['seed_bootstrap_ci_high']}")
    add("lead1_distance_fraction",gs["distance_fraction_of_mean_delta_log_Lambda"]["mean"])
    for k in order:add(f"lead1_auc_{k}",pm[k]["auc"],"held-out seeds 15-19")
    for k in keys+["actual_delta_sigma_perp"]:add(f"segment_event_minus_control_{k}",cg[k]["event_minus_control"],f"95% seed CI {cg[k]['seed_bootstrap_ci_low']} to {cg[k]['seed_bootstrap_ci_high']}")
    add("distance_matched_heading_mean_seed_auc",rs["mean_seed_auc"],f"95% seed CI {rs['seed_bootstrap_ci_low']} to {rs['seed_bootstrap_ci_high']}")
    for r in rp:add(f"distance_matched_{r['metric']}_event_minus_control",r["event_minus_control"],f"95% seed CI {r['seed_bootstrap_ci_low']} to {r['seed_bootstrap_ci_high']}")
    write(out/"summary.csv",summary)

if __name__=="__main__":main()
