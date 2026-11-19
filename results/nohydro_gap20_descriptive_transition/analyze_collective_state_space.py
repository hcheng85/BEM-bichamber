#!/usr/bin/env python3
"""Physical collective state-space audit for no-hydro gap-20 trajectories.

This script tests metastability, slow coordinates, memory, and the population
jump character before any effective potential or rate model is fitted.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from analyze_approach_current_balance import ELL, N, raw_path
from analyze_initiation_propagation import RawReader


DT=.001
TRANSIENT=200.
SAMPLE_DT=1.
Q_BINS=np.asarray([0.,.2,.4,.6,.8,1.000001])
J_BINS=np.asarray([-1.000001,-.25,0.,.25,1.000001])
PI_BINS=np.asarray([0.,.5,1.000001])
B_BINS=np.asarray([0.,.025,.075,.15,1.000001])
LAGS=(1,2,5,10,20)


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");writer.writeheader();writer.writerows(rows)


def chamber_order(state: np.ndarray,mask: np.ndarray) -> tuple[float,float,float,float]:
    if not np.any(mask):return (float("nan"),)*4
    theta=state[2,mask];px=float(np.mean(np.cos(theta)));py=float(np.mean(np.sin(theta)));polar=float(np.hypot(px,py))
    if mask.sum()<2:return px,py,polar,float("nan")
    pos=state[:2,mask];rel=pos-np.mean(pos,axis=1,keepdims=True);radius=np.hypot(rel[0],rel[1]);valid=radius>1e-12
    mill=float(np.mean((rel[0,valid]*np.sin(theta[valid])-rel[1,valid]*np.cos(theta[valid]))/radius[valid])) if np.any(valid) else float("nan")
    return px,py,polar,mill


def features(state: np.ndarray) -> dict:
    x,y,theta=state;right=x>=0;left=~right;q=float(2*np.sum(right)/N-1);e=np.vstack((np.cos(theta),np.sin(theta)))
    px=float(np.mean(e[0]));py=float(np.mean(e[1]));polar=float(np.hypot(px,py));com=np.mean(state[:2],axis=1);rel=state[:2]-com[:,None];radius=np.hypot(rel[0],rel[1]);valid=radius>1e-12
    milling=float(np.mean((rel[0,valid]*e[1,valid]-rel[1,valid]*e[0,valid])/radius[valid]));pi=abs(milling)
    lpx,lpy,lp,lm=chamber_order(state,left);rpx,rpy,rp,rm=chamber_order(state,right)
    right_approach=right&(x<=ELL)&(np.abs(y)<=ELL/2);left_approach=left&((-x)<=ELL)&(np.abs(y)<=ELL/2)
    b_right=float(np.sum(np.clip((ELL-x[right_approach])/ELL,0,1))/N);b_left=float(np.sum(np.clip((ELL+x[left_approach])/ELL,0,1))/N)
    if q==0:
        source_j=target_j=b=float("nan")
    else:
        sign=float(np.sign(q));source=right if sign>0 else left;target=~source;direction=-sign
        source_j=float(direction*np.mean(e[0,source]));target_j=float(direction*np.mean(e[0,target])) if np.any(target) else float("nan")
        u=sign*x;approach=source&(u<=ELL)&(np.abs(y)<=ELL/2);b=float(np.sum(np.clip((ELL-u[approach])/ELL,0,1))/N)
    return dict(q=q,px=px,py=py,P=polar,m=milling,Pi=pi,left_px=lpx,left_P=lp,left_m=lm,right_px=rpx,right_P=rp,right_m=rm,
                source_j=source_j,target_j=target_j,b=b,b_left=b_left,b_right=b_right,b_difference=b_right-b_left,b_total=b_right+b_left,n_left=int(left.sum()),n_right=int(right.sum()))


def load_q_cache(cache_dir: Path,seed: int) -> tuple[np.ndarray,np.ndarray]:
    path=cache_dir/f"seed{seed:02d}_stride100_qscan_v4.npz"
    with np.load(path) as data:return np.asarray(data["step"],int),np.asarray(data["q"],float)


def sample_all(raw_dir: Path,cache_dir: Path,seeds: list[int]) -> tuple[dict[int,list[dict]],dict[int,np.ndarray],dict]:
    series={};dense={};native_delta=[]
    for count,seed in enumerate(seeds,1):
        print(f"collective state snapshots {count}/{len(seeds)} seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed));ntime=reader.info["state"].shape[-1]
        steps=np.arange(round(TRANSIENT/DT),ntime-1,round(SAMPLE_DT/DT),dtype=int);rows=[]
        try:
            for step in steps:
                st=reader.frame("state",int(step));nxt=reader.frame("state",int(step+1));f=features(st);qnext=features(nxt)["q"]
                rows.append(dict(seed=seed,step=int(step),t=step*DT,**f));native_delta.append(qnext-f["q"])
        finally:reader.close()
        series[seed]=rows;cache_steps,q=load_q_cache(cache_dir,seed);dense[seed]=q[cache_steps>=round(TRANSIENT/DT)]
    return series,dense,dict(native_delta=np.asarray(native_delta))


def equal_seed_occupancy(dense: dict[int,np.ndarray]) -> dict:
    values=np.asarray([[np.mean(q>=.8),np.mean(q<=-.8),np.mean(abs(q)<.5)] for q in dense.values()]);rng=np.random.default_rng(20261121);pick=rng.integers(len(values),size=(3000,len(values)));boot=values[pick].mean(axis=1)
    return dict(right_basin=float(values[:,0].mean()),left_basin=float(values[:,1].mean()),central=float(values[:,2].mean()),
                right_ci95=np.percentile(boot[:,0],[2.5,97.5]).tolist(),left_ci95=np.percentile(boot[:,1],[2.5,97.5]).tolist(),central_ci95=np.percentile(boot[:,2],[2.5,97.5]).tolist())


def symmetry_audit(series: dict[int,list[dict]]) -> dict:
    rows=[r for rr in series.values() for r in rr];q=np.asarray([r["q"] for r in rows]);px=np.asarray([r["px"] for r in rows]);m=np.asarray([r["m"] for r in rows])
    return dict(mean_q=float(q.mean()),mean_px=float(px.mean()),mean_signed_milling=float(m.mean()),
                q_odd_moment=float(np.mean(q**3)),px_odd_moment=float(np.mean(px**3)),corr_q_px=float(np.corrcoef(q,px)[0,1]))


def order_parameter_summary(series: dict[int,list[dict]]) -> dict:
    rows=[r for rr in series.values() for r in rr];p=np.asarray([r["P"] for r in rows]);pi=np.asarray([r["Pi"] for r in rows]);source_p=[];target_p=[];source_m=[];target_m=[];source_px=[];target_px=[]
    for r in rows:
        if r["q"]==0:continue
        if r["q"]>0:sp,tp,sm,tm,sx,tx=r["right_P"],r["left_P"],r["right_m"],r["left_m"],-r["right_px"],-r["left_px"]
        else:sp,tp,sm,tm,sx,tx=r["left_P"],r["right_P"],r["left_m"],r["right_m"],r["left_px"],r["right_px"]
        source_p.append(sp);target_p.append(tp);source_m.append(sm);target_m.append(tm);source_px.append(sx);target_px.append(tx)
    def finite_summary(values):
        a=np.asarray(values,float);a=a[np.isfinite(a)];return dict(mean=float(a.mean()),median=float(np.median(a)),q10_q90=np.percentile(a,[10,90]).tolist(),n=len(a))
    return dict(global_P=finite_summary(p),global_Pi=finite_summary(pi),corr_P_Pi=float(np.corrcoef(p,pi)[0,1]),
                source_P=finite_summary(source_p),target_P=finite_summary(target_p),source_signed_axial=finite_summary(source_px),target_signed_axial=finite_summary(target_px),
                source_signed_milling=finite_summary(source_m),target_signed_milling=finite_summary(target_m))


def correlations(series: dict[int,list[dict]]) -> list[dict]:
    out=[]
    for lag in LAGS:
        for field in ("q","px","P","Pi","source_j","b"):
            vals=[]
            for rows in series.values():
                x=np.asarray([r[field] for r in rows],float);ok=np.isfinite(x[:-lag])&np.isfinite(x[lag:])
                if ok.sum()>3:vals.append(float(np.corrcoef(x[:-lag][ok],x[lag:][ok])[0,1]))
            out.append(dict(lag=lag,coordinate=field,mean_seed_autocorrelation=float(np.nanmean(vals)),n_seeds=len(vals)))
    return out


def lagged_singular_values(series: dict[int,list[dict]],fields: tuple[str,...],lag: int=5) -> list[float]:
    x=[];y=[]
    for rows in series.values():
        a=np.asarray([[r[f] for f in fields] for r in rows],float);ok=np.all(np.isfinite(a[:-lag]),axis=1)&np.all(np.isfinite(a[lag:]),axis=1);x.append(a[:-lag][ok]);y.append(a[lag:][ok])
    x=np.vstack(x);y=np.vstack(y);mean=x.mean(axis=0);std=x.std(axis=0);x=(x-mean)/std;y=(y-mean)/std;x-=x.mean(axis=0);y-=y.mean(axis=0)
    c00=x.T@x/len(x)+1e-8*np.eye(x.shape[1]);c11=y.T@y/len(y)+1e-8*np.eye(y.shape[1]);c01=x.T@y/len(x)
    def invsqrt(c):
        value,vec=np.linalg.eigh(c);return (vec/np.sqrt(np.maximum(value,1e-10)))@vec.T
    return np.linalg.svd(invsqrt(c00)@c01@invsqrt(c11),compute_uv=False).tolist()


def transition_events(q: np.ndarray,bound=.6) -> list[tuple[int,int]]:
    previous=0;last_index=None;events=[]
    for i,value in enumerate(q):
        state=1 if value>=bound else (-1 if value<=-bound else 0)
        if state:
            if previous and state!=previous:
                lo=last_index if last_index is not None else 0;segment=q[lo:i+1];center=lo+int(np.argmin(np.abs(segment)));events.append((center,state))
            previous=state;last_index=i
    return events


def transition_paths(series: dict[int,list[dict]],dense: dict[int,np.ndarray]) -> tuple[list[dict],dict]:
    out=[];counts={}
    for seed,q in dense.items():
        events=transition_events(q);counts[seed]=len(events);by_time={round(r["t"]):r for r in series[seed]}
        for eid,(index,target) in enumerate(events):
            center=TRANSIENT+.1*index
            for tau in range(-20,21):
                row=by_time.get(round(center+tau))
                if row is None:continue
                out.append(dict(seed=seed,event_id=eid,target_sign=target,tau=tau,q_oriented=target*row["q"],px_oriented=target*row["px"],
                                q=row["q"],px=row["px"],P=row["P"],Pi=row["Pi"],source_j=row["source_j"],b=row["b"]))
    return out,dict(n_full_relocations=sum(counts.values()),n_seeds=sum(v>0 for v in counts.values()),by_seed=counts)


def landscape(series: dict[int,list[dict]],outdir: Path) -> tuple[list[dict],dict]:
    rows=[r for rr in series.values() for r in rr];q=np.asarray([r["q"] for r in rows]);px=np.asarray([r["px"] for r in rows]);seed=np.asarray([r["seed"] for r in rows])
    q_aug=np.r_[q,-q];px_aug=np.r_[px,-px];seed_aug=np.r_[seed,seed];qe=np.linspace(-1.025,1.025,42);pe=np.linspace(-1.025,1.025,34);hist,_,_=np.histogram2d(q_aug,px_aug,bins=(qe,pe));records=[]
    maxmass=hist.max();supported=0
    for i in range(len(qe)-1):
        for j in range(len(pe)-1):
            cell=(q_aug>=qe[i])&(q_aug<qe[i+1])&(px_aug>=pe[j])&(px_aug<pe[j+1]);n=int(cell.sum());ns=int(len(np.unique(seed_aug[cell])));ok=n>=20 and ns>=5;supported+=ok
            records.append(dict(q_center=(qe[i]+qe[i+1])/2,px_center=(pe[j]+pe[j+1])/2,count=n,n_seeds=ns,supported=int(ok),probability_mass=n/len(q_aug),statistical_neglog=(-math.log(n/maxmass) if ok else "")))
    return records,dict(q_edges=qe,px_edges=pe,hist=hist,supported_cells=supported,total_cells=hist.size)


def bin_index(value: float,edges: np.ndarray) -> int:
    return int(np.clip(np.searchsorted(edges,value,side="right")-1,0,len(edges)-2))


def memory_rows(series: dict[int,list[dict]],future_lag=5,past_lag=2) -> list[dict]:
    out=[]
    for seed,rows in series.items():
        for i in range(past_lag,len(rows)-future_lag):
            r=rows[i];q=r["q"]
            if q==0 or not all(np.isfinite([r["source_j"],r["b"]])):continue
            sign=np.sign(q);past=-sign*(q-rows[i-past_lag]["q"]);future=-sign*(rows[i+future_lag]["q"]-q);j=-sign*r["px"]
            out.append(dict(seed=seed,q_abs=abs(q),j_global=j,j_source=r["source_j"],Pi=r["Pi"],b=r["b"],
                            b_contrast=sign*r["b_difference"],b_total=r["b_total"],past_progress=past,future_progress=future))
    return out


def candidate_key(r: dict,name: str) -> tuple:
    key=[bin_index(r["q_abs"],Q_BINS)]
    if name in ("q_global_axial","q_global_axial_milling","q_global_axial_storage","q_global_axial_storage_contrast","q_global_axial_chamber_storage"):key.append(bin_index(r["j_global"],J_BINS))
    if name in ("q_source_axial","q_source_axial_storage"):key.append(bin_index(r["j_source"],J_BINS))
    if name=="q_global_axial_milling":key.append(bin_index(r["Pi"],PI_BINS))
    if name in ("q_storage","q_global_axial_storage","q_source_axial_storage"):key.append(bin_index(r["b"],B_BINS))
    if name=="q_global_axial_storage_contrast":key.append(bin_index(r["b_contrast"],np.asarray([-1.000001,-.075,0.,.075,1.000001])))
    if name=="q_global_axial_chamber_storage":
        key.extend((bin_index(r["b_contrast"],np.asarray([-1.000001,-.075,0.,.075,1.000001])),bin_index(r["b_total"],np.asarray([0.,.05,.15,.3,2.000001]))))
    return tuple(key)


def coordinate_sufficiency(rows: list[dict]) -> list[dict]:
    train=[r for r in rows if r["seed"]<30];test=[r for r in rows if r["seed"]>=30];out=[];rng=np.random.default_rng(20261122)
    for name in ("q","q_global_axial","q_source_axial","q_storage","q_global_axial_milling","q_global_axial_storage","q_global_axial_storage_contrast","q_global_axial_chamber_storage","q_source_axial_storage"):
        by=defaultdict(list)
        for r in train:by[candidate_key(r,name)].append(r["future_progress"])
        means={k:float(np.mean(v)) for k,v in by.items() if len(v)>=30};supported=[r for r in test if candidate_key(r,name) in means]
        for r in supported:r["_residual"]=r["future_progress"]-means[candidate_key(r,name)]
        seed_stats=[]
        for seed in sorted({r["seed"] for r in supported}):
            part=[r for r in supported if r["seed"]==seed]
            if len(part)>10:seed_stats.append((seed,float(spearmanr([r["past_progress"] for r in part],[r["_residual"] for r in part]).statistic),float(np.mean([r["_residual"]**2 for r in part]))))
        rho=np.asarray([x[1] for x in seed_stats]);mse=np.asarray([x[2] for x in seed_stats]);pick=rng.integers(len(rho),size=(3000,len(rho)))
        out.append(dict(candidate=name,n_train=len(train),n_test=len(test),n_supported=len(supported),support_fraction=len(supported)/len(test),n_cells=len(means),
                        residual_history_spearman=float(rho.mean()),residual_history_ci95=json.dumps(np.percentile(rho[pick].mean(axis=1),[2.5,97.5]).tolist()),
                        future_progress_mse=float(mse.mean()),mse_ci95=json.dumps(np.percentile(mse[pick].mean(axis=1),[2.5,97.5]).tolist()),n_test_seeds=len(seed_stats)))
    return out


def jump_diagnostics(dense: dict[int,np.ndarray],native_delta: np.ndarray) -> list[dict]:
    out=[]
    for lag,dt in ((0,.001),(1,.1),(5,.5),(10,1.),(20,2.)):
        if lag==0:d=native_delta
        else:d=np.concatenate([q[lag:]-q[:-lag] for q in dense.values()])
        second=float(np.mean(d*d));fourth=float(np.mean(d**4));out.append(dict(lag=dt,n=len(d),fraction_zero=float(np.mean(abs(d)<1e-12)),fraction_single_lattice_jump=float(np.mean(np.isclose(abs(d),.02))),
                   fraction_multiple_jump=float(np.mean(abs(d)>.0200001)),mean_abs_increment=float(np.mean(abs(d))),second_kramers_moyal=second/(2*dt),
                   fourth_to_gaussian=float(fourth/(3*second*second)) if second else float("nan")))
    return out


def make_figure(series: dict[int,list[dict]],land: dict,paths: list[dict],corr: list[dict],jumps: list[dict],output: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10,7.8));qe,pe,h=land["q_edges"],land["px_edges"],land["hist"];masked=np.ma.masked_where(h.T<20,-np.log(np.maximum(h.T,1)/h.max()));im=axes[0,0].pcolormesh(qe,pe,masked,shading="auto",cmap="magma_r");fig.colorbar(im,ax=axes[0,0],label="-log empirical mass + const")
    for key in sorted({(r["seed"],r["event_id"]) for r in paths}):
        part=sorted((r for r in paths if (r["seed"],r["event_id"])==key),key=lambda r:r["tau"]);axes[0,0].plot([r["q_oriented"] for r in part],[r["px_oriented"] for r in part],color="cyan",alpha=.55,lw=1)
    axes[0,0].set(xlabel="oriented population imbalance q",ylabel="oriented global axial polarization",title="reflection-pooled occupancy and 9 relocations")
    rows=[r for rr in series.values() for r in rr];axes[0,1].hexbin([r["P"] for r in rows],[r["Pi"] for r in rows],gridsize=32,mincnt=5,bins="log",cmap="viridis");axes[0,1].set(xlabel="global polar order P",ylabel="global milling order Pi",title="collective-motion order parameters")
    ax=axes[1,0]
    for field in ("q","px","P","Pi","source_j","b"):
        part=[r for r in corr if r["coordinate"]==field];ax.plot([r["lag"] for r in part],[r["mean_seed_autocorrelation"] for r in part],"o-",label=field)
    ax.axhline(0,color="k",lw=.6);ax.set(xlabel="lag",ylabel="seed-mean autocorrelation",title="time dependence");ax.legend(frameon=False,ncol=2,fontsize=8)
    ax=axes[1,1];x=[r["lag"] for r in jumps];ax.plot(x,[r["fraction_zero"] for r in jumps],"o-",label="P(delta q=0)");ax.plot(x,[r["fraction_multiple_jump"] for r in jumps],"s-",label="P(|delta q|>0.02)");ax.set_xscale("log");ax.set(xlabel="sampling lag",ylabel="fraction",title="population-lattice jump character");ax.legend(frameon=False)
    fig.suptitle("Minimal collective state audit: occupancy is not yet a dynamical potential",fontsize=11);fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    theta=np.zeros(N);x=np.r_[np.full(60,1.),np.full(40,-1.)];state=np.vstack((x,np.linspace(-1,1,N),theta));f=features(state);ref=state.copy();ref[0]*=-1;ref[2]=np.pi-ref[2];g=features(ref)
    assert abs(f["q"]+g["q"])<1e-12 and abs(f["px"]+g["px"])<1e-12 and abs(f["P"]-g["P"])<1e-12
    assert abs(f["b_difference"]+g["b_difference"])<1e-12 and abs(f["b_total"]-g["b_total"])<1e-12
    q=np.r_[np.full(5,-.8),np.linspace(-.5,.5,11),np.full(5,.8)];assert len(transition_events(q))==1
    print("PASS: reflection-covariant orders, physical storage bounds, and relocation indexing")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--cache-dir",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.cache_dir,args.output)):p.error("--raw-dir, --cache-dir, --output required")
    args.output.mkdir(parents=True,exist_ok=True);seeds=list(range(40));series,dense,audit=sample_all(args.raw_dir,args.cache_dir,seeds);occupancy=equal_seed_occupancy(dense);symmetry=symmetry_audit(series);orders=order_parameter_summary(series);corr=correlations(series)
    paths,relocations=transition_paths(series,dense);land_rows,land=landscape(series,args.output);memory=coordinate_sufficiency(memory_rows(series));jumps=jump_diagnostics(dense,audit["native_delta"])
    slow_sets=[]
    for fields in (("q",),("q","px"),("q","px","P","Pi"),("q","px","source_j","b","Pi")):
        slow_sets.append(dict(coordinates="+".join(fields),lag=5,singular_values=lagged_singular_values(series,fields,5)))
    write_csv(args.output/"stationary_landscape.csv",land_rows);write_csv(args.output/"time_correlations.csv",corr);write_csv(args.output/"coordinate_sufficiency.csv",memory);write_csv(args.output/"transition_paths.csv",paths);write_csv(args.output/"population_jump_diagnostics.csv",jumps)
    make_figure(series,land,paths,corr,jumps,args.output/"collective_state_diagram.png")
    summary=dict(configuration=dict(seeds=seeds,native_dt=DT,stationary_start=TRANSIENT,snapshot_dt=SAMPLE_DT,n_snapshots=sum(map(len,series.values())),reflection_pooling="(q,px,m) with (-q,-px,-m)"),
                 occupancy=occupancy,symmetry=symmetry,order_parameters=orders,relocations=relocations,landscape_support=dict(supported_cells=land["supported_cells"],total_cells=land["total_cells"],minimum_count=20,minimum_seeds=5),
                 slow_coordinate_singular_values=slow_sets,coordinate_sufficiency=memory,population_jump_diagnostics=jumps,
                 warning="-log stationary mass is an occupancy landscape, not a drift-diffusion potential or finite barrier")
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n");(args.output/"analysis_configuration.json").write_text(json.dumps(dict(raw_dir=str(args.raw_dir),cache_dir=str(args.cache_dir),new_simulations=0,fitted_potential=0,nonlinear_latent_coordinate=0),indent=2)+"\n")


if __name__=="__main__":main()
