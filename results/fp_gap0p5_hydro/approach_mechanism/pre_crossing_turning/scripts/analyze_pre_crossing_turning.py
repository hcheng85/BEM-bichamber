#!/usr/bin/env python3
"""Balanced, strictly pre-event angular-mechanism analysis for gap-50 BEM data."""
from __future__ import annotations

import argparse, csv, json, sys
from pathlib import Path
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from replay_dynamics import Parameters, analytic_geometry, replay_snapshot

LAGS = (-2.0, -1.0, -0.2)
N_LEAD = 10


def read_csv(path): return list(csv.DictReader(path.open()))


def write_csv(path, rows):
    if not rows: return
    with path.open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)


def cache_features(z, idx, origin):
    side = "left" if origin == "left" else "right"
    inward = 1 if origin == "left" else -1
    xc=float(z[f"{side}_xc"][idx]); yc=float(z[f"{side}_yc"][idx])
    edge=float(z[f"{side}_xmax"][idx] if origin == "left" else z[f"{side}_xmin"][idx])
    d=-inward*edge
    h=.1
    j=max(0,idx-10); k=min(len(z["time"])-1,idx+10)
    vy=(float(z[f"{side}_yc"][k])-float(z[f"{side}_yc"][j]))/(float(z["time"][k])-float(z["time"][j]))
    return np.array([yc, vy, d, float(z[f"{side}_xstd"][idx]), float(z[f"{side}_ystd"][idx])])


def select_triplets(episodes, cache_dir):
    rows=[]
    for seed in range(15,20):
        with np.load(cache_dir/f"gap50_hydro_approach_seed{seed:02d}_stride10.npz") as archive:
            z={key:archive[key] for key in archive.files}
        candidates=[r for r in episodes if int(r["seed"])==seed and float(r["full_dwell_before"])>=2.2
                    and float(r["first_cross_time"])+float(r["episode_duration"])<999.9]
        best=None
        for origin in ("left","right"):
            succ=[r for r in candidates if r["origin"]==origin and r["success"]=="1"]
            fail=[r for r in candidates if r["origin"]==origin and r["success"]=="0"]
            if len(succ)<2 or len(fail)<2: continue
            allf=[]
            for r in succ+fail:
                idx=int(round((float(r["first_cross_time"])-2)/.01)); allf.append(cache_features(z,idx,origin))
            scale=np.nanstd(np.vstack(allf),axis=0); scale[scale<1e-6]=1
            costs=[]
            for a in succ:
                fa=cache_features(z,int(round((float(a["first_cross_time"])-2)/.01)),origin)
                for b in fail:
                    fb=cache_features(z,int(round((float(b["first_cross_time"])-2)/.01)),origin)
                    score=float(np.linalg.norm((fa-fb)/scale))
                    costs.append((score,a,b,(fa+fb)/2))
            chosen=[]; used_s=set(); used_f=set()
            for item in sorted(costs,key=lambda x:x[0]):
                sa=float(item[1]["first_cross_time"]); fb=float(item[2]["first_cross_time"])
                if sa not in used_s and fb not in used_f:
                    chosen.append(item); used_s.add(sa); used_f.add(fb)
                if len(chosen)==2: break
            objective=float(np.mean([x[0] for x in chosen]))
            if len(chosen)==2 and (best is None or objective<best[0]):
                best=(objective,origin,chosen,scale)
        if best is None: raise RuntimeError(f"no balanced event pair for seed {seed}")
        objective,origin,pairs,scale=best
        n=z["left_n"] if origin=="left" else z["right_n"]
        exact=n==100; time=z["time"]
        # Prospective control anchor: full occupancy from -2.2 through +4, and
        # not near any observed A time. Candidate stride 0.1 time units.
        valid=np.ones(len(time),bool)
        pre=int(round(2.2/.01)); post=int(round(4/.01))
        valid[:pre]=False; valid[-post:]=False
        bad=np.flatnonzero(~exact)
        diff=np.zeros(len(time)+1,int)
        for b in bad:
            lo=max(0,b-post); hi=min(len(time),b+pre+1); diff[lo]+=1; diff[hi]-=1
        valid &= np.cumsum(diff[:-1])==0
        for e in [r for r in episodes if int(r["seed"])==seed]:
            valid &= np.abs(time-float(e["first_cross_time"]))>5
        inds=np.flatnonzero(valid & ((np.arange(len(time))%10)==0))
        used_controls=[]
        for replicate,(score,succ,fail,target) in enumerate(pairs,1):
            eligible=[i for i in inds if all(abs(float(time[i])-float(time[j]))>5 for j in used_controls)]
            control=min(eligible,key=lambda i: np.linalg.norm((cache_features(z,i-int(round(2/.01)),origin)-target)/scale))
            used_controls.append(control)
            # Note: cache feature target is at A-2; for pseudo-anchor t0 use t0-2.
            for label,event in (("success",succ),("failure",fail)):
                feat=cache_features(z,int(round((float(event["first_cross_time"])-2)/.01)),origin)
                rows.append(dict(seed=seed,replicate=replicate,condition=label,origin=origin,
                                 anchor_time=float(event["first_cross_time"]),pair_distance=score,
                                 success=int(event["success"]),control_match_distance="",
                                 match_y_center=feat[0],match_vertical_com_velocity=feat[1],
                                 match_front_distance=feat[2],match_x_extent=feat[3],match_y_extent=feat[4]))
            cf=cache_features(z,control-int(round(2/.01)),origin)
            rows.append(dict(seed=seed,replicate=replicate,condition="control",origin=origin,
                             anchor_time=float(time[control]),pair_distance=score,success="",
                             control_match_distance=float(np.linalg.norm((cf-target)/scale)),
                             match_y_center=cf[0],match_vertical_com_velocity=cf[1],
                             match_front_distance=cf[2],match_x_extent=cf[3],match_y_extent=cf[4]))
    return rows


def angular_geometry(state, origin, gap, tip_r=0.5):
    s=1.0 if origin=="left" else -1.0
    d=-s*state[0]
    psi=np.arctan2(state[1]*0+s*np.sin(state[2]), s*np.cos(state[2]))
    # Exact limiting rays are tangents from the fish to the rounded tip
    # circles, not rays to their inner endpoints on x=0.
    cy_lo=-gap/2-tip_r; cy_hi=gap/2+tip_r
    alpha_lo=np.arctan2(cy_lo-state[1],d); alpha_hi=np.arctan2(cy_hi-state[1],d)
    beta_lo=np.arcsin(np.clip(tip_r/np.hypot(d,cy_lo-state[1]),0,1))
    beta_hi=np.arcsin(np.clip(tip_r/np.hypot(d,cy_hi-state[1]),0,1))
    lo=alpha_lo+beta_lo; hi=alpha_hi-beta_hi
    # The interval lies in the forward half-plane. For a backward-pointing fish,
    # choose the shorter periodic rotation to either endpoint (never > pi).
    inside=(psi>=lo)&(psi<=hi)
    dlo=np.arctan2(np.sin(psi-lo),np.cos(psi-lo))
    dhi=np.arctan2(np.sin(psi-hi),np.cos(psi-hi))
    err=np.where(inside,0.0,np.where(np.abs(dlo)<=np.abs(dhi),dlo,dhi))
    return s,d,psi,lo,hi,err


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--analysis-dir",type=Path,required=True)
    ap.add_argument("--raw-dir",type=Path,required=True)
    ap.add_argument("--cache-dir",type=Path,required=True)
    ap.add_argument("--geometry-file",type=Path,required=True)
    ap.add_argument("--lu-file",type=Path,required=True)
    ap.add_argument("--piv-file",type=Path,required=True)
    args=ap.parse_args(); out=args.analysis_dir.resolve(); (out/"tables").mkdir(parents=True,exist_ok=True); (out/"figures").mkdir(exist_ok=True)
    p=Parameters(); geom=analytic_geometry(p)
    episodes=read_csv(out.parent/"tables"/"nucleation_episodes.csv")
    selected=select_triplets(episodes,args.cache_dir); write_csv(out/"tables"/"selected_triplets.csv",selected)
    g=sio.loadmat(args.geometry_file); normal=np.asarray(g["normal_vec"])
    boundary=dict(x=np.asarray(g["Sx"]).reshape(-1),y=np.asarray(g["Sy"]).reshape(-1),
                  weights=np.asarray(g["weights"]).reshape(-1),nx=normal[:,0],ny=normal[:,1])
    lu=np.load(args.lu_file,mmap_mode="r"); piv=np.load(args.piv_file)
    fish_rows=[]; validation=[]
    for seed in range(15,20):
        anchors=[r for r in selected if int(r["seed"])==seed]
        raw=args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state_all=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
        for a in anchors:
            origin=a["origin"]
            for lag in LAGS:
                t=float(a["anchor_time"])+lag; idx=int(round(t/p.dt))
                state=state_all[:,:,idx]; nxt=state_all[:,:,idx+1]
                if origin=="left" and np.any(state[0]>=0): raise RuntimeError("pre-event left snapshot not full")
                if origin=="right" and np.any(state[0]<0): raise RuntimeError("pre-event right snapshot not full")
                rep=replay_snapshot(state,nxt,boundary,lu,piv,p,geom)
                s,d,psi,lo,hi,err=angular_geometry(state,origin,p.gap)
                lead=np.argsort(d)[:N_LEAD]
                validation.append(dict(seed=seed,condition=a["condition"],lag=lag,
                    all_fish_noise_z_mean=float(np.mean(rep["noise_z"])),
                    all_fish_noise_z_std=float(np.std(rep["noise_z"],ddof=1)),
                    leading_noise_z_mean=float(np.mean(rep["noise_z"][lead])),
                    leading_noise_z_std=float(np.std(rep["noise_z"][lead],ddof=1)),
                    overridden_all=int(np.sum(rep["overridden"])),
                    overridden_leading=int(np.sum(rep["overridden"][lead])),
                    rms_reconstruction_error=float(np.sqrt(np.mean((rep["observed_increment"]-rep["deterministic_total"]*p.dt-rep["noise_increment"])**2)))))
                for rank,i in enumerate(lead,1):
                    sg=0 if err[i]==0 else np.sign(err[i])
                    def closing(omega): return float(-sg*s*omega[i]) if sg else np.nan
                    row=dict(seed=seed,condition=a["condition"],origin=origin,lag=lag,
                             replicate=int(a["replicate"]),anchor_time=a["anchor_time"],fish=int(i),leading_rank=rank,
                             x=float(state[0,i]),y=float(state[1,i]),distance_to_midline=float(d[i]),
                             local_heading=float(psi[i]),admissible_low=float(lo[i]),admissible_high=float(hi[i]),
                             signed_angular_error=float(err[i]),abs_angular_error=float(abs(err[i])),
                             admissible=int(err[i]==0),wall_distance=float(rep["wall_distance"][i]),
                             wall_torque=float(rep["terms"]["wall"][i]),
                             alignment_torque=float(rep["terms"]["alignment"][i]),
                             attraction_torque=float(rep["terms"]["attraction"][i]),
                             direct_hydro_torque=float(rep["terms"]["direct_hydro_effective"][i]),
                             boundary_hydro_torque=float(rep["terms"]["boundary_hydro_effective"][i]),
                             deterministic_total=float(rep["deterministic_total"][i]),
                             observed_heading_increment=float(rep["observed_increment"][i]),
                             inferred_noise_increment=float(rep["noise_increment"][i]),noise_z=float(rep["noise_z"][i]),
                             overridden=int(rep["overridden"][i]))
                    for key,col in (("wall","wall_closing_rate"),("alignment","alignment_closing_rate"),
                                    ("attraction","attraction_closing_rate"),("direct_hydro_effective","direct_hydro_closing_rate"),
                                    ("boundary_hydro_effective","boundary_hydro_closing_rate")):
                        row[col]=closing(rep["terms"][key])
                    row["deterministic_closing_rate"]=closing({"x":rep["deterministic_total"]} if False else rep["deterministic_total"])
                    row["noise_closing_increment"]=float(-sg*s*rep["noise_increment"][i]) if sg else np.nan
                    fish_rows.append(row)
        del state_all
    write_csv(out/"tables"/"leading_fish_torques.csv",fish_rows); write_csv(out/"tables"/"replay_validation.csv",validation)

    metrics=["abs_angular_error","admissible","wall_closing_rate","alignment_closing_rate","attraction_closing_rate",
             "direct_hydro_closing_rate","boundary_hydro_closing_rate","deterministic_closing_rate","noise_closing_increment"]
    seed_rows=[]
    for seed in range(15,20):
      for cond in ("success","failure","control"):
       for lag in LAGS:
        sub=[r for r in fish_rows if r["seed"]==seed and r["condition"]==cond and r["lag"]==lag]
        row=dict(seed=seed,condition=cond,lag=lag,n_leading=len(sub))
        for m in metrics:
            v=np.array([r[m] for r in sub],float); row[m]=float(np.nanmean(v)) if np.any(np.isfinite(v)) else np.nan
        seed_rows.append(row)
    write_csv(out/"tables"/"seed_level_summary.csv",seed_rows)
    group=[]
    rng=np.random.default_rng(20261021)
    for cond in ("success","failure","control"):
      for lag in LAGS:
       sub=[r for r in seed_rows if r["condition"]==cond and r["lag"]==lag]
       for m in metrics:
        v=np.array([r[m] for r in sub],float); v=v[np.isfinite(v)]
        if len(v):
            boot=np.mean(v[rng.integers(0,len(v),(10000,len(v)))],axis=1)
            group.append(dict(condition=cond,lag=lag,metric=m,mean=float(np.mean(v)),median=float(np.median(v)),
                              seed_bootstrap_ci_low=float(np.percentile(boot,2.5)),seed_bootstrap_ci_high=float(np.percentile(boot,97.5)),independent_seeds=len(v)))
    write_csv(out/"tables"/"group_summary.csv",group)

    # Figure 1: mean leading-edge geometry and torque decomposition, paired by seed.
    colors={"success":"#009E73","failure":"#D55E00","control":"#777777"}
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),constrained_layout=True)
    for cond in colors:
        sub=[r for r in seed_rows if r["condition"]==cond]
        for key,ax,ylabel in (("abs_angular_error",axes[0],"angular distance to opening (rad)"),
                             ("admissible",axes[1],"fraction of leading edge admissible"),
                             ("deterministic_closing_rate",axes[2],"deterministic closing rate (rad/time)")):
            means=[np.nanmean([r[key] for r in sub if r["lag"]==lag]) for lag in LAGS]
            ax.plot([-x for x in LAGS],means,marker="o",label=cond,color=colors[cond])
    for ax in axes: ax.grid(alpha=.2); ax.set_xlabel("time before anchor A (time units)"); ax.invert_xaxis()
    axes[0].set_ylabel("angular distance to opening (rad)"); axes[1].set_ylabel("fraction admissible")
    axes[2].axhline(0,color="k",lw=.8); axes[2].set_ylabel("closing rate; positive = toward opening")
    axes[0].legend(frameon=False)
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"01_pre_event_angular_dynamics.{ext}",dpi=220)
    plt.close(fig)

    # Figure 2: physical decomposition at each lag and restoring behavior.
    terms=[("wall_closing_rate","wall"),("alignment_closing_rate","alignment"),("attraction_closing_rate","attraction"),
           ("direct_hydro_closing_rate","direct hydro"),("boundary_hydro_closing_rate","boundary hydro")]
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    width=.05; xs=np.arange(len(LAGS),dtype=float)
    for j,(key,label) in enumerate(terms):
        vals=[]
        for lag in LAGS:
            vals.append(np.nanmean([r[key] for r in seed_rows if r["condition"]=="success" and r["lag"]==lag]))
        axes[0].bar(xs+(j-2)*width,vals,width,label=label)
    axes[0].axhline(0,color="k",lw=.8); axes[0].set_xticks(xs,[str(x) for x in LAGS]); axes[0].set_xlabel("lag from A")
    axes[0].set_ylabel("success: closing rate (rad/time)"); axes[0].legend(frameon=False,fontsize=8)
    outside=[r for r in fish_rows if np.isfinite(r["wall_closing_rate"])]
    for cond in colors:
        q=[r for r in outside if r["condition"]==cond]
        axes[1].scatter([r["signed_angular_error"] for r in q],[r["wall_closing_rate"] for r in q],s=14,alpha=.45,color=colors[cond],label=cond)
    axes[1].axhline(0,color="k",lw=.8); axes[1].axvline(0,color="k",lw=.5)
    axes[1].set_xlabel("signed angular error (rad)"); axes[1].set_ylabel("wall contribution toward opening")
    axes[1].legend(frameon=False,fontsize=8); axes[1].grid(alpha=.2)
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"02_torque_decomposition_and_restoring.{ext}",dpi=220)
    plt.close(fig)

    config=dict(lags=list(LAGS),leading_edge_fish=N_LEAD,seeds=list(range(15,20)),selection="two within-seed success/failure/control triplets; match chamber, y center, vertical COM velocity, front distance, x/y extent; never heading; exclude end-censored episodes",parameters=p.__dict__)
    (out/"configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    print(json.dumps(dict(selected=len(selected),fish_rows=len(fish_rows),validation=len(validation)),indent=2))

if __name__=="__main__": main()
