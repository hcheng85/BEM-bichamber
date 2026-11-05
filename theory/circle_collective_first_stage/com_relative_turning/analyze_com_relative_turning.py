#!/usr/bin/env python3
"""Source-exact COM-relative turning analysis for circular confinement.

The production simulator and large trajectory files are immutable.  MATLAB-v5
arrays are opened read-only and only bounded snapshots are replayed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.ndimage import uniform_filter1d


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))

from matv5_memmap import open_state_rdot  # noqa: E402
from process_collective import FILE_PATTERN, seed_from_name  # noqa: E402
from replay_torques import load_production, social_split  # noqa: E402


RADIUS = 7.158910531638177
DT = 0.001
N_FISH = 100
IN = 0.5
COM_SINGULAR_CUTOFF = 1e-8 * RADIUS


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def wrap(x: np.ndarray) -> np.ndarray:
    return np.arctan2(np.sin(x), np.cos(x))


def t_block_stats(times: np.ndarray, values: np.ndarray, nblocks: int=10) -> dict[str, float]:
    edges=np.linspace(float(np.min(times)),float(np.max(times))+1e-12,nblocks+1)
    means=[]
    for i in range(nblocks):
        use=(times>=edges[i])&(times<(edges[i+1]) if i<nblocks-1 else times<=edges[i+1])
        if np.any(use): means.append(float(np.mean(values[use])))
    means=np.asarray(means); mean=float(np.mean(values))
    se=float(np.std(means,ddof=1)/np.sqrt(len(means)))
    critical=float(stats.t.ppf(.975,len(means)-1))
    return {"mean":mean,"block_se":se,"ci_low":mean-critical*se,
            "ci_high":mean+critical*se,"n_blocks":len(means)}


def state_com_quantities(state: np.ndarray, velocity: np.ndarray | None=None,
                         cutoff: float=COM_SINGULAR_CUTOFF) -> dict[str, np.ndarray | float | int]:
    x,y,theta=state; cx=float(np.mean(x)); cy=float(np.mean(y))
    yx=x-cx; yy=y-cy; rr2=yx*yx+yy*yy; rr=np.sqrt(rr2)
    valid=rr>cutoff; phi=np.arctan2(yy,yx); alpha=wrap(theta-phi)
    result: dict[str,np.ndarray|float|int]={
        "yx":yx,"yy":yy,"rr":rr,"valid":valid,"phi":phi,"alpha":alpha,
        "M_C":float(np.mean(np.sin(alpha[valid]))),
        "n_valid":int(np.count_nonzero(valid)),
        "n_excluded":int(np.count_nonzero(~valid)),
        "n_within_1e3R":int(np.count_nonzero(rr<=1e-3*RADIUS)),
        "n_within_1e2R":int(np.count_nonzero(rr<=1e-2*RADIUS)),
        "min_r_com_over_R":float(np.min(rr)/RADIUS),
        "com_x":cx,"com_y":cy,
    }
    if velocity is not None:
        vx,vy=velocity; vbarx=float(np.mean(vx)); vbary=float(np.mean(vy))
        omega=np.full(N_FISH,np.nan)
        omega[valid]=(yx[valid]*(vy[valid]-vbary)-yy[valid]*(vx[valid]-vbarx))/rr2[valid]
        result.update(omega_C=omega,vbar_x=vbarx,vbar_y=vbary)
    return result


class COMReplay:
    def __init__(self,data_dir:Path,source:Path,metrics_dir:Path):
        self.data_dir=data_dir; self.source=source; self.metrics_dir=metrics_dir
        self.mod=load_production(source)
        self.para=dict(Ip=4.0,In=IN,R=self.mod.R,light=1.0,dt=DT,Iw=.94,
            If=.01,num=N_FISH,delta=.01,total_time=1000,gap=self.mod.gap,
            radius_ratio=50,gap_ratio=50)
        files=sorted(data_dir.glob(FILE_PATTERN),key=seed_from_name)
        if len(files)!=10: raise FileNotFoundError(f"Expected 10 trajectories, found {len(files)}")
        self.paths={seed_from_name(f):f for f in files}; self.maps={}; self.metrics={}; self.cache={}

    def arrays(self,seed:int):
        if seed not in self.maps: self.maps[seed]=open_state_rdot(self.paths[seed])[:2]
        return self.maps[seed]

    def metric(self,seed:int):
        if seed not in self.metrics: self.metrics[seed]=np.load(self.metrics_dir/f"metrics_seed{seed}.npz")
        return self.metrics[seed]

    def evaluate(self,seed:int,step:int) -> tuple[dict[str,Any],dict[str,np.ndarray]]:
        key=(seed,step)
        if key in self.cache: return self.cache[key]
        state_mm,rdot_mm=self.arrays(seed)
        state=np.asarray(state_mm[:,:,step],float); x,y,theta=state
        c=np.cos(theta); sn=np.sin(theta); propulsion=np.vstack((c,sn))
        direct_u,direct_w=self.mod.GetHydro(state,self.para)
        bem_u,bem_w=self.mod.GetHydroBEM(state,self.para)
        vn=self.mod.GetVoronoi(state[:2],self.para)
        alignment,attraction=social_split(self.mod,vn,state,self.para["Ip"])
        social_direct=np.asarray(self.mod.GetVision(vn,state,self.para)).ravel()
        wall=np.asarray(self.mod.GetWallAvoid(state,self.para)).ravel()
        direct_w=np.asarray(direct_w).ravel(); bem_w=np.asarray(bem_w).ravel()

        trial_velocity=propulsion+direct_u+bem_u
        collision=np.asarray(self.mod.Check_jumping(state,trial_velocity,self.para),bool)
        extreme=np.any(np.abs(trial_velocity)>500,axis=0); angular_override=collision|extreme
        effective_velocity=trial_velocity.copy(); effective_velocity[:,angular_override]=propulsion[:,angular_override]
        direct_eff=direct_w.copy(); bem_eff=bem_w.copy()
        direct_eff[angular_override]=0; bem_eff[angular_override]=0
        components={"alignment":alignment,"attraction":attraction,"wall":wall,
                    "direct_hydro":direct_eff,"bem_hydro":bem_eff}
        raw_total=sum(components.values()); omega_eff=np.clip(raw_total,-np.pi/DT,np.pi/DT)
        components["clip_residual"]=omega_eff-raw_total

        # Reproduce the later positional retry.  It changes translation only.
        next_trial=state[:2]+effective_velocity*DT
        dis_new=np.hypot(next_trial[0],next_trial[1])-self.mod.R
        dis_old=np.hypot(x,y)-self.mod.R
        sticky=(dis_new-self.mod.safedis>=0)&(dis_old<dis_new)
        effective_velocity[:,sticky]=0
        stored_velocity=np.asarray(rdot_mm[:,:,step+1],float)

        cq=state_com_quantities(state,effective_velocity); valid=np.asarray(cq["valid"])
        alpha=np.asarray(cq["alpha"]); omega_c=np.asarray(cq["omega_C"]); cosalpha=np.cos(alpha)
        nvalid=int(cq["n_valid"])
        source={name:float(np.sum(cosalpha[valid]*value[valid])/nvalid)
                for name,value in components.items()}
        geometry=float(-np.sum(cosalpha[valid]*omega_c[valid])/nvalid)
        noise=float(-.5*IN**2*float(cq["M_C"]))
        deterministic=float(sum(source.values())+geometry)
        total=deterministic+noise
        angular_total=sum(source.values())
        omega_valid=omega_eff[valid]; geom_valid=omega_c[valid]
        omega_corr=float(np.corrcoef(omega_valid,geom_valid)[0,1])

        # Exact conditional finite-step expectation: positions are deterministic
        # after the stored retry; angular increments are independent Gaussians.
        next_position=np.asarray(state_mm[:2,:,step+1],float)
        ncx,ncy=np.mean(next_position,axis=1); nyx=next_position[0]-ncx; nyy=next_position[1]-ncy
        nrr=np.hypot(nyx,nyy); next_valid=valid&(nrr>COM_SINGULAR_CUTOFF)
        common=int(np.count_nonzero(next_valid))
        next_phi=np.arctan2(nyy,nyx)
        alpha_det=theta+omega_eff*DT-next_phi
        decay=np.exp(-.5*IN**2*DT)
        expected_next=float(np.mean(decay*np.sin(alpha_det[next_valid])))
        observed_theta=np.asarray(state_mm[2,:,step+1],float)
        observed_next=float(np.mean(np.sin((observed_theta-next_phi)[next_valid])))
        current_common=float(np.mean(np.sin(alpha[next_valid])))
        finite_expected=(expected_next-current_common)/DT
        finite_observed=(observed_next-current_common)/DT
        var_each=.5*(1-np.exp(-2*IN**2*DT)*np.cos(2*alpha_det[next_valid]))-(decay*np.sin(alpha_det[next_valid]))**2
        martingale_variance=float(np.sum(var_each)/(common*common))
        martingale_increment=observed_next-expected_next

        metric=self.metric(seed); mi=min(int(step//10),len(metric["P"])-1)
        row={"seed":seed,"step":step,"time":step*DT,"P":float(metric["P"][mi]),
            "Pi":float(metric["Pi"][mi]),"M_C":float(cq["M_C"]),
            "bMC_alignment":source["alignment"],"bMC_attraction":source["attraction"],
            "bMC_wall":source["wall"],"bMC_direct_hydro":source["direct_hydro"],
            "bMC_bem_hydro":source["bem_hydro"],"bMC_clip_residual":source["clip_residual"],
            "bMC_angular_total":angular_total,"bMC_geometry":geometry,
            "bMC_deterministic":deterministic,"bMC_noise_ito":noise,
            "bMC_total":total,"bMC_finite_step_expected":finite_expected,
            "bMC_finite_step_observed":finite_observed,
            "martingale_increment":martingale_increment,"martingale_variance":martingale_variance,
            "martingale_z":martingale_increment/np.sqrt(max(martingale_variance,1e-30)),
            "n_com_valid":nvalid,"n_com_common_next":common,"n_com_excluded":int(cq["n_excluded"]),
            "n_within_1e3R":int(cq["n_within_1e3R"]),"n_within_1e2R":int(cq["n_within_1e2R"]),
            "min_r_com_over_R":float(cq["min_r_com_over_R"]),
            "omega_effective_rms":float(np.sqrt(np.mean(omega_valid**2))),
            "omega_com_geometry_rms":float(np.sqrt(np.mean(geom_valid**2))),
            "omega_minus_geometry_rms":float(np.sqrt(np.mean((omega_valid-geom_valid)**2))),
            "omega_geometry_correlation":omega_corr,
            "social_split_max_error":float(np.max(np.abs(alignment+attraction-social_direct))),
            "velocity_replay_max_error":float(np.max(np.abs(effective_velocity-stored_velocity))),
            "n_angular_clip":int(np.count_nonzero(components["clip_residual"])),
            "n_collision_override":int(np.count_nonzero(collision)),
            "n_extreme_override":int(np.count_nonzero(extreme)),"n_sticky_override":int(np.count_nonzero(sticky))}
        detail={"alpha":alpha,"omega_C":omega_c,"omega_effective":omega_eff,"components":components,
                "effective_velocity":effective_velocity,"stored_velocity":stored_velocity}
        self.cache[key]=(row,detail); return row,detail


def stream_com_series(replay:COMReplay,outdir:Path,stride:int=100) -> dict[int,dict[str,np.ndarray]]:
    series={}
    for seed in range(10):
        state,_=replay.arrays(seed); steps=np.arange(20000,state.shape[2]-1,stride,dtype=np.int64)
        mc=np.empty(len(steps),np.float32); minimum=np.empty(len(steps),np.float32)
        near1=np.empty(len(steps),np.uint8); near2=np.empty(len(steps),np.uint8)
        for start in range(0,len(steps),1000):
            ss=steps[start:start+1000]; batch=np.asarray(state[:,:,ss],float)
            if batch.shape[0]!=3: batch=np.moveaxis(batch,-1,0)
            x,y,theta=batch; cx=np.mean(x,axis=0); cy=np.mean(y,axis=0)
            yx=x-cx; yy=y-cy; rr=np.hypot(yx,yy); alpha=theta-np.arctan2(yy,yx)
            sl=slice(start,start+len(ss)); mc[sl]=np.mean(np.sin(alpha),axis=0)
            minimum[sl]=np.min(rr,axis=0)/RADIUS
            near1[sl]=np.sum(rr<=1e-3*RADIUS,axis=0); near2[sl]=np.sum(rr<=1e-2*RADIUS,axis=0)
        series[seed]={"step":steps,"time":steps*DT,"M_C":mc,"min_r_com_over_R":minimum,
                      "n_within_1e3R":near1,"n_within_1e2R":near2}
        print(f"streamed COM series seed {seed}: {len(steps)} frames",flush=True)
    np.savez_compressed(outdir/"com_relative_series.npz",**{
        f"{key}_seed{seed}":value for seed,data in series.items() for key,value in data.items()})
    return series


def interval_steps(t0:float,t1:float,n:int) -> np.ndarray:
    return np.unique(np.linspace(round(t0/DT),round(t1/DT),n,dtype=np.int64))


def stationary_analysis(replay:COMReplay,series:dict[int,dict[str,np.ndarray]],outdir:Path):
    intervals=[*(dict(seed=s,t0=20.,t1=999.9,state="schooling",n=100) for s in range(7)),
        dict(seed=7,t0=20.,t1=999.9,state="milling",n=200),
        dict(seed=8,t0=20.,t1=999.9,state="milling",n=200),
        dict(seed=9,t0=20.,t1=280.,state="schooling_pre_transition",n=150),
        dict(seed=9,t0=340.,t1=999.9,state="milling_post_transition",n=200)]
    snapshots=[]; summaries=[]
    fields=["bMC_alignment","bMC_attraction","bMC_wall","bMC_direct_hydro","bMC_bem_hydro",
            "bMC_clip_residual","bMC_geometry","bMC_deterministic","bMC_noise_ito","bMC_total",
            "bMC_angular_total","bMC_finite_step_expected","bMC_finite_step_observed",
            "omega_effective_rms","omega_com_geometry_rms","omega_minus_geometry_rms","omega_geometry_correlation"]
    for interval in intervals:
        seed=interval["seed"]; steps=interval_steps(interval["t0"],interval["t1"],interval["n"])
        rows=[]
        for step in steps:
            row,_=replay.evaluate(seed,int(step)); row={"interval":interval["state"],**row}; rows.append(row); snapshots.append(row)
        times=np.array([r["time"] for r in rows]); mc=np.array([r["M_C"] for r in rows])
        chirality=1 if np.mean(mc)>=0 else -1
        ss=series[seed]; use=(ss["time"]>=interval["t0"])&(ss["time"]<=interval["t1"])
        endpoint=float((ss["M_C"][use][-1]-ss["M_C"][use][0])/(ss["time"][use][-1]-ss["time"][use][0]))
        summary={"seed":seed,"interval":interval["state"],"t0":interval["t0"],"t1":interval["t1"],
            "n_replayed":len(rows),"chirality":chirality,"mean_M_C_stream":float(np.mean(ss["M_C"][use])),
            "mean_abs_M_C_stream":float(np.mean(np.abs(ss["M_C"][use]))),"endpoint_dM_C_dt":endpoint,
            "min_r_com_over_R":float(np.min(ss["min_r_com_over_R"][use])),
            "frames_with_fish_within_1e3R":int(np.count_nonzero(ss["n_within_1e3R"][use])),
            "frames_with_fish_within_1e2R":int(np.count_nonzero(ss["n_within_1e2R"][use]))}
        for field in fields:
            values=np.array([r[field] for r in rows]); st=t_block_stats(times,values)
            summary[field]=st["mean"]
            if field in ("bMC_total","bMC_finite_step_expected","bMC_finite_step_observed"):
                summary[field+"_block_se"]=st["block_se"]; summary[field+"_ci_low"]=st["ci_low"]
                summary[field+"_ci_high"]=st["ci_high"]
        summary["total_ci_contains_zero"]=int(summary["bMC_total_ci_low"]<=0<=summary["bMC_total_ci_high"])
        summary["endpoint_inside_expected_ci"]=int(summary["bMC_finite_step_expected_ci_low"]<=endpoint<=summary["bMC_finite_step_expected_ci_high"])
        summary["aligned_alignment"]=chirality*summary["bMC_alignment"]
        summary["aligned_attraction"]=chirality*summary["bMC_attraction"]
        summary["aligned_wall"]=chirality*summary["bMC_wall"]
        summary["aligned_direct_hydro"]=chirality*summary["bMC_direct_hydro"]
        summary["aligned_bem_hydro"]=chirality*summary["bMC_bem_hydro"]
        summary["aligned_geometry"]=chirality*summary["bMC_geometry"]
        summary["aligned_noise_ito"]=chirality*summary["bMC_noise_ito"]
        summaries.append(summary)
        print(f"stationary replay seed {seed} {interval['state']}: {len(rows)}",flush=True)
    write_csv(outdir/"stationary_snapshot_replay.csv",snapshots)
    write_csv(outdir/"stationary_balance_by_seed.csv",summaries)
    return snapshots,summaries


def transition_analysis(replay:COMReplay,outdir:Path):
    rows=[]
    for step in np.arange(300000,330001,100,dtype=np.int64):
        row,_=replay.evaluate(9,int(step)); rows.append(row)
    write_csv(outdir/"seed9_transition_com_drift.csv",rows)
    windows=[("baseline",300,306),("circulation_loss",306,311.5),("polarization_collapse",311.5,315),("early_milling",315,325)]
    fields=["M_C","P","Pi","bMC_alignment","bMC_attraction","bMC_wall","bMC_direct_hydro",
            "bMC_bem_hydro","bMC_clip_residual","bMC_geometry","bMC_noise_ito","bMC_total",
            "bMC_angular_total","bMC_finite_step_expected","omega_effective_rms",
            "omega_com_geometry_rms","omega_minus_geometry_rms","omega_geometry_correlation"]
    window_rows=[]
    for name,t0,t1 in windows:
        use=[r for r in rows if t0<=r["time"]<t1]
        window_rows.append({"window":name,"t0":t0,"t1":t1,"n":len(use),
            **{field:float(np.mean([r[field] for r in use])) for field in fields}})
    write_csv(outdir/"seed9_transition_window_means.csv",window_rows)

    t=np.array([r["time"] for r in rows]); width=11
    def smooth(field): return uniform_filter1d(np.array([r[field] for r in rows],float),width,mode="nearest")
    fig,axs=plt.subplots(4,1,figsize=(11,12),sharex=True)
    axs[0].plot(t,smooth("P"),label="P"); axs[0].plot(t,smooth("Pi"),label=r"$\Pi$")
    axs[0].plot(t,smooth("M_C"),label=r"$M_C$"); axs[0].legend(ncol=3); axs[0].set_ylabel("order parameter")
    for field,label in (("bMC_alignment","alignment"),("bMC_attraction","attraction"),("bMC_wall","wall")):
        axs[1].plot(t,smooth(field),label=label,lw=.9)
    axs[1].legend(ncol=3); axs[1].set_ylabel(r"turning part of $b_{M_C}$")
    for field,label in (("bMC_direct_hydro","direct hydro"),("bMC_bem_hydro","BEM hydro"),("bMC_geometry","geometry")):
        axs[2].plot(t,smooth(field),label=label,lw=.9)
    axs[2].legend(ncol=3); axs[2].set_ylabel(r"hydro / geometry")
    axs[3].plot(t,smooth("bMC_deterministic"),label="deterministic")
    axs[3].plot(t,smooth("bMC_noise_ito"),label="Ito noise decay")
    axs[3].plot(t,smooth("bMC_total"),label="total generator")
    axs[3].axhline(0,color="black",lw=.5); axs[3].legend(ncol=3); axs[3].set_ylabel(r"$b_{M_C}$"); axs[3].set_xlabel("time")
    for ax in axs: ax.axvline(312.41,color="black",ls="--",lw=.8); ax.grid(alpha=.2)
    fig.suptitle("Seed 9: source-exact COM-relative turning balance (1.1-time smoothing)")
    fig.tight_layout(); fig.savefig(outdir/"seed9_com_relative_turning.png",dpi=220); plt.close(fig)
    return rows,window_rows


def schooling_control_windows(replay:COMReplay,outdir:Path,event_windows:list[dict]):
    """Compare the 5.5-time pre-collapse window with prespecified controls."""
    fields=["bMC_alignment","bMC_attraction","bMC_wall","bMC_direct_hydro",
            "bMC_bem_hydro","bMC_geometry","bMC_noise_ito","bMC_total"]
    controls=[]
    for seed in range(7):
        for t0 in (100.,500.,900.):
            steps=np.arange(round(t0/DT),round((t0+5.5)/DT),100,dtype=np.int64)
            rows=[replay.evaluate(seed,int(step))[0] for step in steps]
            chirality=1 if np.mean([r["M_C"] for r in rows])>=0 else -1
            controls.append({"seed":seed,"t0":t0,"t1":t0+5.5,"n":len(rows),"chirality":chirality,
                "mean_M_C":float(np.mean([r["M_C"] for r in rows])),
                **{field:float(np.mean([r[field] for r in rows])) for field in fields},
                **{"aligned_"+field[4:]:chirality*float(np.mean([r[field] for r in rows])) for field in fields}})
    write_csv(outdir/"schooling_control_window_means.csv",controls)
    event=next(r for r in event_windows if r["window"]=="circulation_loss")
    event_chirality=1
    comparison=[]
    for field in fields:
        key="aligned_"+field[4:]; values=np.array([r[key] for r in controls])
        event_value=event_chirality*event[field]
        seed_means=np.array([np.mean([r[key] for r in controls if r["seed"]==seed]) for seed in range(7)])
        se=float(np.std(seed_means,ddof=1)/np.sqrt(len(seed_means)))
        crit=float(stats.t.ppf(.975,len(seed_means)-1)); mean=float(np.mean(seed_means))
        comparison.append({"component":field,"event_aligned_mean":event_value,
            "control_window_mean":float(np.mean(values)),"control_window_q025":float(np.quantile(values,.025)),
            "control_window_q975":float(np.quantile(values,.975)),"control_seed_mean":mean,
            "control_seed_ci_low":mean-crit*se,"control_seed_ci_high":mean+crit*se,
            "event_outside_control_window_95pct":int(event_value<np.quantile(values,.025) or event_value>np.quantile(values,.975))})
    write_csv(outdir/"precollapse_control_comparison.csv",comparison)

    fig,ax=plt.subplots(figsize=(10,5)); x=np.arange(len(fields))
    values=[]; lo=[]; hi=[]; event_values=[]
    for field in fields:
        key="aligned_"+field[4:]; v=np.array([r[key] for r in controls]); values.append(np.mean(v))
        lo.append(np.quantile(v,.025)); hi.append(np.quantile(v,.975)); event_values.append(event[field])
    values=np.array(values); lo=np.array(lo); hi=np.array(hi)
    ax.errorbar(x,values,yerr=[values-lo,hi-values],fmt="o",capsize=4,label="21 schooling control windows")
    ax.scatter(x,event_values,marker="D",color="#e45756",label="seed 9, t=306..311.5")
    ax.axhline(0,color="black",lw=.6); ax.set_xticks(x,[f[4:].replace("_"," ") for f in fields],rotation=25,ha="right")
    ax.set_ylabel("chirality-aligned contribution to $b_{M_C}$"); ax.legend(); ax.grid(axis="y",alpha=.2)
    fig.tight_layout(); fig.savefig(outdir/"precollapse_vs_schooling_controls.png",dpi=220); plt.close(fig)
    return controls,comparison


def validation_audit(source:Path,snapshots:list[dict],transition:list[dict],outdir:Path):
    rows=snapshots+transition; z=np.array([r["martingale_z"] for r in rows])
    audit={"source":str(source),"source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "velocity_index":"rdot[:,:,step+1] advances state[:,:,step]",
        "n_replayed_snapshots":len(rows),"n_fish_updates":len(rows)*N_FISH,
        "maximum_velocity_replay_error":float(max(r["velocity_replay_max_error"] for r in rows)),
        "maximum_social_split_error":float(max(r["social_split_max_error"] for r in rows)),
        "override_totals":{key:int(sum(r[key] for r in rows)) for key in
            ("n_angular_clip","n_collision_override","n_extreme_override","n_sticky_override")},
        "com_excluded_total":int(sum(r["n_com_excluded"] for r in rows)),
        "fish_snapshots_within_1e3R":int(sum(r["n_within_1e3R"] for r in rows)),
        "fish_snapshots_within_1e2R":int(sum(r["n_within_1e2R"] for r in rows)),
        "minimum_r_com_over_R":float(min(r["min_r_com_over_R"] for r in rows)),
        "finite_step_martingale_z_mean":float(np.mean(z)),"finite_step_martingale_z_variance":float(np.var(z,ddof=1)),
        "finite_step_martingale_z_skew":float(stats.skew(z)),"finite_step_martingale_z_excess_kurtosis":float(stats.kurtosis(z)),
        "theory":"d alpha_i^C = (Omega_i-omega_i^C)dt + In dW_i; b_MC = mean[cos(alpha_i^C)(Omega_i-omega_i^C)] - In^2 M_C/2",
        "com_cutoff":COM_SINGULAR_CUTOFF,"com_cutoff_over_R":COM_SINGULAR_CUTOFF/RADIUS}
    (outdir/"replay_validation.json").write_text(json.dumps(audit,indent=2)+"\n")
    return audit


def summary_figures(summaries:list[dict],outdir:Path):
    components=[("aligned_alignment","alignment"),("aligned_attraction","attraction"),("aligned_wall","wall"),
        ("aligned_direct_hydro","direct hydro"),("aligned_bem_hydro","BEM hydro"),("aligned_geometry","geometry"),
        ("aligned_noise_ito","noise decay")]
    school=[r for r in summaries if "schooling" in r["interval"]]
    mill=[r for r in summaries if "milling" in r["interval"]]
    group_rows=[]
    for rows,label in ((school,"schooling"),(mill,"milling")):
        for key,_ in components:
            values=np.array([r[key] for r in rows]); mean=float(np.mean(values))
            se=float(np.std(values,ddof=1)/np.sqrt(len(values))); crit=float(stats.t.ppf(.975,len(values)-1))
            group_rows.append({"macrostate":label,"quantity":key,"n_seed_intervals":len(values),
                "mean":mean,"seed_se":se,"ci_low":mean-crit*se,"ci_high":mean+crit*se})
        values=np.array([r["mean_abs_M_C_stream"] for r in rows]); mean=float(np.mean(values))
        se=float(np.std(values,ddof=1)/np.sqrt(len(values))); crit=float(stats.t.ppf(.975,len(values)-1))
        group_rows.append({"macrostate":label,"quantity":"mean_abs_M_C_stream","n_seed_intervals":len(values),
            "mean":mean,"seed_se":se,"ci_low":mean-crit*se,"ci_high":mean+crit*se})
    write_csv(outdir/"stationary_group_comparison.csv",group_rows)
    fig,ax=plt.subplots(figsize=(10,5)); x=np.arange(len(components)); width=.36
    for offset,rows,label,color in ((-.5,school,"schooling","#4c78a8"),(.5,mill,"milling","#e45756")):
        mean=np.array([np.mean([r[k] for r in rows]) for k,_ in components])
        se=np.array([np.std([r[k] for r in rows],ddof=1)/np.sqrt(len(rows)) for k,_ in components])
        ax.bar(x+offset*width,mean,width,yerr=se,capsize=3,label=label,color=color,alpha=.85)
    ax.axhline(0,color="black",lw=.6); ax.set_xticks(x,[label for _,label in components],rotation=25,ha="right")
    ax.set_ylabel("chirality-aligned contribution to $b_{M_C}$"); ax.legend(); ax.grid(axis="y",alpha=.2)
    fig.tight_layout(); fig.savefig(outdir/"stationary_com_turning_decomposition.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(7,5))
    for rows,label,color in ((school,"schooling","#4c78a8"),(mill,"milling","#e45756")):
        ax.scatter([r["mean_abs_M_C_stream"] for r in rows],[r["bMC_total"] for r in rows],label=label,color=color,s=45)
        for r in rows: ax.annotate(str(r["seed"]),(r["mean_abs_M_C_stream"],r["bMC_total"]),fontsize=7)
    ax.axhline(0,color="black",lw=.7); ax.set(xlabel=r"interval mean $|M_C|$",ylabel=r"mean exact generator drift $b_{M_C}$")
    ax.grid(alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(outdir/"stationary_balance_by_seed.png",dpi=220); plt.close(fig)


def main():
    p=argparse.ArgumentParser(); p.add_argument("--data-dir",type=Path,required=True)
    p.add_argument("--source",type=Path,default=PARENT.parents[1].parent/"BEM_circular_parallel_instance.py")
    p.add_argument("--outdir",type=Path,default=HERE/"outputs")
    args=p.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    params=json.loads((PARENT/"simulation_parameters.json").read_text())
    source=args.source.resolve(); source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    if source_hash!=params["source_sha256"]: raise RuntimeError(f"Source hash mismatch: {source_hash}")
    replay=COMReplay(args.data_dir,source,PARENT/"outputs")
    series=stream_com_series(replay,args.outdir)
    snapshots,summaries=stationary_analysis(replay,series,args.outdir)
    transition,windows=transition_analysis(replay,args.outdir)
    controls,control_comparison=schooling_control_windows(replay,args.outdir,windows)
    control_steps=sum(int(r["n"]) for r in controls)
    audit=validation_audit(source,snapshots+[
        {**replay.evaluate(seed,int(step))[0]} for seed in range(7)
        for t0 in (100.,500.,900.) for step in np.arange(round(t0/DT),round((t0+5.5)/DT),100,dtype=np.int64)
    ],transition,args.outdir)
    summary_figures(summaries,args.outdir)
    config={"radius":RADIUS,"dt":DT,"N":N_FISH,"In":IN,"data_pattern":FILE_PATTERN,
        "stream_stride_steps":100,"stream_dt":.1,"stationary_school_snapshots_per_seed":100,
        "stationary_milling_snapshots_per_seed":200,"seed9_pre_snapshots":150,"seed9_post_snapshots":200,
        "transition_window":[300,330],"transition_stride_steps":100,"transition_dt":.1,
        "schooling_control_windows":{"seeds":list(range(7)),"starts":[100,500,900],"duration":5.5,"dt":.1},
        "uncertainty":"10 contiguous time blocks with Student-t 95% interval",
        "com_singular_cutoff":COM_SINGULAR_CUTOFF,"com_singular_cutoff_over_R":COM_SINGULAR_CUTOFF/RADIUS,
        "source_sha256":source_hash}
    (args.outdir/"analysis_config.json").write_text(json.dumps(config,indent=2)+"\n")
    print(json.dumps({"outdir":str(args.outdir),"snapshots":len(snapshots)+len(transition),
        "control_snapshots":control_steps,"max_velocity_error":audit["maximum_velocity_replay_error"]},indent=2),flush=True)


if __name__=="__main__": main()
