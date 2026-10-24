#!/usr/bin/env python3
"""Seed-cluster summaries and decisive figures for the rounded-tip test."""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment

CLASSES = ("deflection", "pinning", "transmission", "repelling")
RNG = np.random.default_rng(2601023)


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def f(row, key):
    try: return float(row[key])
    except (ValueError, KeyError): return np.nan


def cluster_ci(rows, statistic, nboot=10000):
    byseed = defaultdict(list)
    for r in rows: byseed[int(r["seed"])].append(r)
    seeds = sorted(byseed)
    point = statistic(rows)
    vals = np.empty(nboot)
    for b in range(nboot):
        sampled = RNG.choice(seeds, len(seeds), replace=True)
        q = []
        for newseed, seed in enumerate(sampled):
            for r in byseed[seed]:
                z = dict(r); z["seed"] = newseed; q.append(z)
        vals[b] = statistic(q)
    vals = vals[np.isfinite(vals)]
    return point, np.percentile(vals, 2.5), np.percentile(vals, 97.5), len(seeds)


def clustered_proportion(rows, predicate, nboot=10000):
    """Fast whole-seed bootstrap for a binary proportion."""
    seeds = sorted({int(r["seed"]) for r in rows})
    num = np.array([sum(predicate(r) for r in rows if int(r["seed"]) == s) for s in seeds], float)
    den = np.array([sum(1 for r in rows if int(r["seed"]) == s) for s in seeds], float)
    point = num.sum()/den.sum()
    take = RNG.integers(0, len(seeds), (nboot, len(seeds)))
    vals = num[take].sum(axis=1)/den[take].sum(axis=1)
    return point, *np.percentile(vals, [2.5,97.5]), len(seeds)


def class_frequencies(primary):
    rows = []
    groupings = [
        ("overall", lambda r: "all"),
        ("tip", lambda r: r["tip"]),
        ("condition", lambda r: r["condition"]),
        ("offset", lambda r: r["observation_offset"]),
        ("phase", lambda r: r["collective_vertical_phase"]),
        ("tip_phase", lambda r: r["tip"]+"_"+r["collective_vertical_phase"]),
    ]
    valid = [r for r in primary if r["classification"] != "invalid"]
    dvals = np.array([f(r, "divider_distance") for r in valid])
    q1, q2 = np.quantile(dvals, [1/3, 2/3])
    for r in valid:
        d = f(r, "divider_distance")
        r["wall_distance_bin"] = "near" if d <= q1 else ("middle" if d <= q2 else "far")
    groupings.append(("wall_distance", lambda r: r["wall_distance_bin"]))
    for grouping, fun in groupings:
        labels = sorted({fun(r) for r in valid})
        for label in labels:
            subset = [r for r in valid if fun(r) == label]
            for cls in CLASSES:
                point, lo, hi, ns = clustered_proportion(
                    subset, lambda x, c=cls: x["classification"] == c)
                rows.append(dict(grouping=grouping, group=label, classification=cls,
                                 n=len(subset), pooled_frequency=point,
                                 seed_bootstrap_ci_low=lo, seed_bootstrap_ci_high=hi,
                                 independent_seeds=ns, divider_tertile_q1=q1,
                                 divider_tertile_q2=q2, Iw=0.94))
    return rows


def branch_summary(primary):
    valid = [r for r in primary if r["classification"] != "invalid"]
    rows = []
    for tip in ("upper", "lower"):
        q = [r for r in valid if r["tip"] == tip]
        sign = 1 if tip == "upper" else -1
        def qwall(r, which):
            side = 1 if r["origin"] == "left" else -1
            return sign*side*f(r, which+"_wall")
        for which in ("out", "in"):
            a = np.array([qwall(r, which) for r in q])
            rows.append(dict(tip=tip, side=which, n=len(a),
                q_directed_wall_mean=np.mean(a), q_directed_wall_median=np.median(a),
                q_directed_wall_q05=np.quantile(a,.05), q_directed_wall_q95=np.quantile(a,.95),
                expected_local_torque_sign=("negative" if tip == "upper" else "positive")))
        jumps = np.array([qwall(r,"out")-qwall(r,"in") for r in q])
        rows.append(dict(tip=tip, side="jump_out_minus_in", n=len(jumps),
            q_directed_wall_mean=np.mean(jumps), q_directed_wall_median=np.median(jumps),
            q_directed_wall_q05=np.quantile(jumps,.05), q_directed_wall_q95=np.quantile(jumps,.95),
            expected_local_torque_sign="toward corridor when negative"))
    pairs = Counter((r["tip"], r.get("out_ray_branch",""), r.get("in_ray_branch","")) for r in valid)
    for (tip, outb, inb), n in pairs.items():
        rows.append(dict(tip=tip, side="branch_pair", n=n, out_branch=outb, in_branch=inb,
                         frequency=n/sum(v for (t,_,_),v in pairs.items() if t == tip)))
    return rows


def epsilon_sensitivity(allrows):
    keyed = defaultdict(dict)
    for r in allrows:
        key = (r["seed"],r["pair_id"],r["condition"],r["observation_offset"],r["tip"])
        keyed[key][r["epsilon"]] = r["classification"]
    pairs = Counter((v.get("0.001"), v.get("0.005")) for v in keyed.values())
    total = sum(pairs.values())
    return [dict(class_epsilon_1e3=a, class_epsilon_5e3=b, n=n, frequency=n/total)
            for (a,b),n in sorted(pairs.items())]


def seed_mean_ci(values):
    seeds = sorted(values)
    a = np.array([values[s] for s in seeds], float)
    point = np.nanmean(a)
    boot = np.nanmean(a[RNG.integers(0,len(a),(20000,len(a)))], axis=1)
    return point, *np.nanpercentile(boot,[2.5,97.5]), len(a), a


def association(primary):
    """Descriptive prospective association; no fitted surrogate model."""
    eligible = [r for r in primary if r["classification"] in ("deflection","pinning","transmission")
                and r["initially_outside"] == "1"]
    rows = []
    for offset in (0.0,1.5,1.8):
        for horizon in (0.2,0.5):
            tag = str(horizon).replace(".","p")
            for outcome in ("angular_entry", "focal_crossing", "any_crossing"):
                key = f"{outcome}_{tag}"
                q = [r for r in eligible if abs(f(r,"observation_offset")-offset)<1e-9]
                positive = [r for r in q if r["classification"] in ("pinning","transmission")]
                negative = [r for r in q if r["classification"] == "deflection"]
                rp = np.mean([int(r[key]) for r in positive]) if positive else np.nan
                rn = np.mean([int(r[key]) for r in negative]) if negative else np.nan
                perseed = {}
                for seed in range(20):
                    pp=[int(r[key]) for r in positive if int(r["seed"])==seed]
                    nn=[int(r[key]) for r in negative if int(r["seed"])==seed]
                    if pp and nn: perseed[seed]=np.mean(pp)-np.mean(nn)
                if perseed:
                    rd,lo,hi,ns,arr=seed_mean_ci(perseed)
                else: rd=lo=hi=np.nan;ns=0;arr=[]
                rows.append(dict(analysis="unadjusted_class_association",offset=offset,
                    horizon=horizon,outcome=outcome,n_positive=len(positive),n_deflection=len(negative),
                    positive_rate=rp,deflection_rate=rn,seed_mean_risk_difference=rd,
                    seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,independent_seeds=ns,
                    seed_differences=";".join(f"{v:.8g}" for v in arr)))
    return rows


def geometry_matched_class_association(primary):
    """Match access classes within seed/origin/tip on eta and divider distance."""
    eligible = [r for r in primary if r["classification"] in ("deflection","pinning","transmission")
                and r["initially_outside"] == "1"]
    pairrows=[]; summaries=[]
    for offset in (0.0,1.5,1.8):
      qoff=[r for r in eligible if abs(f(r,"observation_offset")-offset)<1e-9]
      eta_scale=np.std([f(r,"initial_abs_angular_error") for r in qoff]) or 1
      d_scale=np.std([f(r,"divider_distance") for r in qoff]) or 1
      for seed in range(20):
       for origin in ("left","right"):
        for tip in ("upper","lower"):
         q=[r for r in qoff if int(r["seed"])==seed and r["origin"]==origin and r["tip"]==tip]
         pos=[r for r in q if r["classification"] in ("pinning","transmission")]
         neg=[r for r in q if r["classification"]=="deflection"]
         if not pos or not neg: continue
         zp=np.array([[f(r,"initial_abs_angular_error")/eta_scale,f(r,"divider_distance")/d_scale] for r in pos])
         zn=np.array([[f(r,"initial_abs_angular_error")/eta_scale,f(r,"divider_distance")/d_scale] for r in neg])
         cost=np.sqrt(np.sum((zp[:,None,:]-zn[None,:,:])**2,axis=2))
         ii,jj=linear_sum_assignment(cost)
         for i,j in zip(ii,jj):
          a,b=pos[i],neg[j]
          row=dict(seed=seed,origin=origin,tip=tip,offset=offset,
                   standardized_match_distance=cost[i,j],positive_class=a["classification"])
          for horizon in (.2,.5):
           tag=str(horizon).replace(".","p")
           for outcome in ("angular_entry","focal_crossing","any_crossing"):
            key=f"{outcome}_{tag}"
            row[f"delta_{key}"]=int(a[key])-int(b[key])
          pairrows.append(row)
      pairs=[r for r in pairrows if r["offset"]==offset]
      for horizon in (.2,.5):
       tag=str(horizon).replace(".","p")
       for outcome in ("angular_entry","focal_crossing","any_crossing"):
        key=f"delta_{outcome}_{tag}"; perseed={}
        for seed in range(20):
         vals=[r[key] for r in pairs if r["seed"]==seed]
         if vals: perseed[seed]=np.mean(vals)
        point,lo,hi,ns,arr=seed_mean_ci(perseed)
        summaries.append(dict(analysis="geometry_matched_class_association",offset=offset,
            horizon=horizon,outcome=outcome,n_pairs=len(pairs),mean_standardized_match_distance=np.mean([r["standardized_match_distance"] for r in pairs]),
            seed_mean_risk_difference=point,seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,
            independent_seeds=ns,seed_differences=";".join(f"{v:.8g}" for v in arr)))
    return pairrows,summaries


def event_control_matching(primary):
    """Prospective event/control rematching on eta and divider distance only."""
    eligible=[r for r in primary if r["classification"] in ("deflection","pinning","transmission")
              and r["initially_outside"]=="1"]
    pairs=[]; summaries=[]
    for offset in (0.0,1.5,1.8):
      qoff=[r for r in eligible if abs(f(r,"observation_offset")-offset)<1e-9]
      es=np.std([f(r,"initial_abs_angular_error") for r in qoff]) or 1
      ds=np.std([f(r,"divider_distance") for r in qoff]) or 1
      for seed in range(20):
       for origin in ("left","right"):
        for tip in ("upper","lower"):
         q=[r for r in qoff if int(r["seed"])==seed and r["origin"]==origin and r["tip"]==tip]
         ev=[r for r in q if r["condition"]=="event"]; co=[r for r in q if r["condition"]=="control"]
         if not ev or not co: continue
         ze=np.array([[f(r,"initial_abs_angular_error")/es,f(r,"divider_distance")/ds] for r in ev])
         zc=np.array([[f(r,"initial_abs_angular_error")/es,f(r,"divider_distance")/ds] for r in co])
         cost=np.sqrt(np.sum((ze[:,None,:]-zc[None,:,:])**2,axis=2));ii,jj=linear_sum_assignment(cost)
         for i,j in zip(ii,jj):
          a,b=ev[i],co[j]
          pairs.append(dict(seed=seed,origin=origin,tip=tip,offset=offset,
            standardized_match_distance=cost[i,j],
            delta_access_class=int(a["classification"] in ("pinning","transmission"))-int(b["classification"] in ("pinning","transmission")),
            delta_pinning=int(a["classification"]=="pinning")-int(b["classification"]=="pinning"),
            event_class=a["classification"],control_class=b["classification"]))
      pp=[r for r in pairs if r["offset"]==offset]
      for metric in ("delta_access_class","delta_pinning"):
       perseed={}
       for seed in range(20):
        vals=[r[metric] for r in pp if r["seed"]==seed]
        if vals:perseed[seed]=np.mean(vals)
       point,lo,hi,ns,arr=seed_mean_ci(perseed)
       summaries.append(dict(analysis="event_control_geometry_matched",offset=offset,metric=metric,
          n_pairs=len(pp),mean_standardized_match_distance=np.mean([r["standardized_match_distance"] for r in pp]),
          seed_mean_difference=point,seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,
          independent_seeds=ns,seed_differences=";".join(f"{v:.8g}" for v in arr)))
    return pairs,summaries


def figures(root, primary, freq, assoc, branch, ident, sweeps):
    figdir=root/"figures";figdir.mkdir(exist_ok=True)
    # Geometry and exact planar identity.
    fig,axs=plt.subplots(1,4,figsize=(18,4.2))
    z=np.array([np.log(f(ident[0],"perpendicular_distance")/f(r,"perpendicular_distance")) for r in ident])
    theta=np.array([f(r,"theta") for r in ident]);pred=np.array([f(r,"identity_prediction") for r in ident])
    axs[0].plot(z,theta,label="exact wall-only integration",lw=2);axs[0].plot(z,pred,"--",label=r"$\theta_0+I_w\log(d_0/d)$")
    axs[0].set(xlabel=r"$\log(d_0/d)$",ylabel=r"$\theta$",title="Planar-wall identity");axs[0].legend(fontsize=8)
    colors={"upper_tip":"#d95f02","lower_tip":"#1b9e77","outer_circle":"#7570b3",
            "upper_left_smooth":"#66a61e","upper_right_smooth":"#e6ab02",
            "lower_left_smooth":"#a6761d","lower_right_smooth":"#e7298a"}
    cases=[c for c in ("lower_archived","upper_archived","closed_archived") if any(r["case"]==c for r in sweeps)]
    for ax,case in zip(axs[1:],cases):
      q=[r for r in sweeps if r["case"]==case and r["tip"]==( "upper" if "upper" in case else "lower")]
      for branch_name in sorted({r["ray_branch"] for r in q}):
       w=[r for r in q if r["ray_branch"]==branch_name]
       ax.plot([f(r,"delta") for r in w],[f(r,"wall_torque") for r in w],".",ms=4,label=branch_name,color=colors.get(branch_name))
      ax.axvline(0,color="k",lw=1);ax.set(xlabel=r"$\theta-\alpha_{tip}$",ylabel="global wall torque",title=case.replace("_"," "))
      ax.legend(fontsize=6)
    fig.tight_layout();fig.savefig(figdir/"01_exact_geometry_and_wall_identity.png",dpi=220);fig.savefig(figdir/"01_exact_geometry_and_wall_identity.pdf")
    plt.close(fig)

    # One-sided drift plane and class frequencies.
    fig,axs=plt.subplots(1,2,figsize=(11,4.4))
    valid=[r for r in primary if r["classification"]!="invalid"]
    cmap={"deflection":"#d95f02","pinning":"#1b9e77","transmission":"#7570b3","repelling":"#666666"}
    for cls in CLASSES:
      q=[r for r in valid if r["classification"]==cls]
      axs[0].scatter([f(r,"mu_out") for r in q],[f(r,"mu_in") for r in q],s=6,alpha=.28,label=f"{cls} ({len(q)})",color=cmap[cls])
    lim=30;axs[0].axhline(0,color="k",lw=.8);axs[0].axvline(0,color="k",lw=.8);axs[0].set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel=r"$\mu_{out}$",ylabel=r"$\mu_{in}$",title="Total drift relative to moving tangent");axs[0].legend(fontsize=7)
    labels=[];bottom=np.zeros(6);x=np.arange(6)
    groups=[]
    for off in (0.0,1.5,1.8):
      for cond in ("event","control"):
       groups.append([r for r in valid if abs(f(r,"observation_offset")-off)<1e-9 and r["condition"]==cond]);labels.append(f"A-{2-off:.1f}\n{cond[0].upper()}")
    for cls in ("deflection","pinning","transmission"):
      vals=np.array([np.mean([r["classification"]==cls for r in g]) for g in groups])
      axs[1].bar(x,vals,bottom=bottom,label=cls,color=cmap[cls]);bottom+=vals
    axs[1].set(xticks=x,xticklabels=labels,ylabel="frequency among valid tip evaluations",title="Drift classes by lead and outcome stratum");axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(figdir/"02_one_sided_drift_classes.png",dpi=220);fig.savefig(figdir/"02_one_sided_drift_classes.pdf");plt.close(fig)

    # Prospective risk differences.
    fig,axs=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
    for ax,h in zip(axs,(.2,.5)):
      q=[r for r in assoc if r["analysis"]=="unadjusted_class_association" and r["horizon"]==h]
      xpos=np.arange(len(q)); y=np.array([r["seed_mean_risk_difference"] for r in q]);lo=np.array([r["seed_bootstrap_ci_low"] for r in q]);hi=np.array([r["seed_bootstrap_ci_high"] for r in q])
      ax.errorbar(xpos,y,yerr=[y-lo,hi-y],fmt="o",capsize=3,color="#1b9e77");ax.axhline(0,color="k",lw=.8)
      ax.set(xticks=xpos,xticklabels=[f"A-{2-r['offset']:.1f}\n{r['outcome'].replace('_',' ')}" for r in q],title=f"future window {h}",ylabel="access-class minus deflection risk");ax.tick_params(axis="x",labelrotation=45)
    fig.tight_layout();fig.savefig(figdir/"03_prospective_association.png",dpi=220);fig.savefig(figdir/"03_prospective_association.pdf");plt.close(fig)


def geometry_verification(sweeps):
    rows=[]
    for case in sorted({r["case"] for r in sweeps}):
      for tip in ("upper","lower"):
       q=[r for r in sweeps if r["case"]==case and r["tip"]==tip]
       if not q:continue
       outside_delta=1e-3 if tip=="upper" else -1e-3
       inside_delta=-outside_delta
       out=min(q,key=lambda r:abs(f(r,"delta")-outside_delta));inn=min(q,key=lambda r:abs(f(r,"delta")-inside_delta))
       near=min(q,key=lambda r:abs(f(r,"delta")))
       rows.append(dict(case=case,tip=tip,corridor_valid=near["corridor_valid"],
          corridor_width=near["corridor_width"],analytic_tangent=near["tangent_angle"],
          sweep_resolution=5e-4,outside_branch=out["ray_branch"],inside_branch=inn["ray_branch"],
          outside_wall_torque=out["wall_torque"],inside_wall_torque=inn["wall_torque"],
          branch_changed=int(out["ray_branch"]!=inn["ray_branch"])))
    return rows


def trace_summary(trace):
    rows=[]
    keys=sorted({(r["seed"],r["pair_id"],r["tip"],r["reference_class"]) for r in trace})
    for key in keys:
      q=[r for r in trace if (r["seed"],r["pair_id"],r["tip"],r["reference_class"])==key]
      pre=[r for r in q if f(r,"time_from_crossing")<0]
      branches=[r["selected_ray_branch"] for r in pre]
      rows.append(dict(seed=key[0],pair_id=key[1],tip=key[2],reference_class=key[3],n_pre_crossing=len(pre),
        reference_class_fraction=np.mean([r["instantaneous_tangent_class"]==key[3] for r in pre]),
        pinning_fraction=np.mean([r["instantaneous_tangent_class"]=="pinning" for r in pre]),
        transmission_fraction=np.mean([r["instantaneous_tangent_class"]=="transmission" for r in pre]),
        deflection_fraction=np.mean([r["instantaneous_tangent_class"]=="deflection" for r in pre]),
        ray_branch_switches=sum(a!=b for a,b in zip(branches[:-1],branches[1:])),
        q_min=min(f(r,"q") for r in pre),q_max=max(f(r,"q") for r in pre)))
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--analysis-dir",type=Path,required=True)
    ap.add_argument("--legacy-observed-gate-velocity",action="store_true",
                    help="summarize the preserved legacy observed-velocity replay")
    args=ap.parse_args()
    root=args.analysis_dir.resolve();tab=root/"tables"
    legacy=args.legacy_observed_gate_velocity
    allrows=read(tab/("tangent_snapshot_replay.csv" if legacy else "tangent_snapshot_replay_corrected.csv"))
    primary=[r for r in allrows if r["epsilon"]=="0.001"]
    freq=class_frequencies(primary);branch=branch_summary(primary);eps=epsilon_sensitivity(allrows)
    assoc=association(primary);matched_pairs,matched=geometry_matched_class_association(primary)
    ec_pairs,ec_summary=event_control_matching(primary)
    prefix="" if legacy else "corrected_"
    write(tab/f"{prefix}class_frequencies.csv",freq);write(tab/f"{prefix}branch_jump_summary.csv",branch)
    write(tab/f"{prefix}epsilon_sensitivity.csv",eps);write(tab/f"{prefix}prospective_association.csv",assoc+matched)
    write(tab/f"{prefix}geometry_matched_class_pairs.csv",matched_pairs);write(tab/f"{prefix}event_control_geometry_matches.csv",ec_pairs)
    write(tab/f"{prefix}event_control_matched_summary.csv",ec_summary)
    ident=read(tab/"planar_wall_identity.csv");sweeps=read(tab/"ray_branch_sweeps.csv")
    write(tab/"geometry_verification.csv",geometry_verification(sweeps))
    tracepath=tab/"representative_traces.csv"
    if tracepath.exists():write(tab/"representative_trace_summary.csv",trace_summary(read(tracepath)))
    if not legacy: figures(root,primary,freq,assoc,branch,ident,sweeps)
    valid=[r for r in primary if r["classification"]!="invalid"]
    overview=[dict(metric="n_prospective_snapshots",value=len(primary)//2),
      dict(metric="n_valid_tip_evaluations",value=len(valid)),
      dict(metric="pinning_frequency",value=np.mean([r["classification"]=="pinning" for r in valid])),
      dict(metric="epsilon_class_agreement",value=np.mean([x["class_epsilon_1e3"]==x["class_epsilon_5e3"] for x in eps for _ in range(int(x["n"]))])),
      dict(metric="planar_identity_max_abs_error",value=max(abs(f(r,"identity_error")) for r in ident if r["branch"]=="left_vertical")),
      dict(metric="closed_or_invalid_tip_evaluations",value=sum(r["classification"]=="invalid" for r in primary)),
      dict(metric="unique_Iw_values",value=len({r["Iw"] for r in primary}))]
    write(root/("summary_legacy_observed_gate_velocity.csv" if legacy else "summary.csv"),overview)


if __name__=="__main__":main()
