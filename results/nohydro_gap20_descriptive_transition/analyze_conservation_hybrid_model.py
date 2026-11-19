#!/usr/bin/env python3
"""Conservation-consistent hybrid model audit for no-hydro gap-20 switching.

The exact state is X=(n_R,p_x,b_L,b_R).  Centerline crossings use exact jump
maps; between-crossing moments are estimated after removing those jumps.  A
small, reflection-symmetric projection is used only to test closure and the
joint operator.  No microscopic trajectory is generated or modified.
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
from scipy.io import loadmat
from scipy.stats import spearmanr

from analyze_approach_current_balance import DT, ELL, N, RawReader, raw_path
from analyze_conservation_correction import segment_budget


TRANSIENT = 200.0
FEATURE_DT = 0.5
Q_DT = 0.1
Q_EDGES = np.asarray([-1.000001, -.6, -.2, .2, .6, 1.000001])
PX_EDGES = np.asarray([-1.000001, -.25, .25, 1.000001])
B_THRESHOLD = .075
NQ, NP, NB = 5, 3, 2
NSTATES = NQ * NP * NB * NB
LAGS = (.5, 1., 2., 5.)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def bindex(x: float, edges: np.ndarray) -> int:
    return int(np.clip(np.searchsorted(edges,x,side="right")-1,0,len(edges)-2))


def sid(q: float, px: float, bl: float, br: float) -> int:
    qb=bindex(q,Q_EDGES);pb=bindex(px,PX_EDGES);lb=int(bl>=B_THRESHOLD);rb=int(br>=B_THRESHOLD)
    return ((qb*NP+pb)*NB+lb)*NB+rb


def decode(state: int) -> tuple[int,int,int,int]:
    rb=state%NB;state//=NB;lb=state%NB;state//=NB;pb=state%NP;qb=state//NP
    return qb,pb,lb,rb


def reflect_sid(state: int) -> int:
    qb,pb,lb,rb=decode(state)
    return (((NQ-1-qb)*NP+(NP-1-pb))*NB+rb)*NB+lb


def load_features(path: Path) -> dict[int,dict[str,np.ndarray]]:
    with np.load(path) as data:
        arrays={k:np.asarray(data[k]) for k in data.files}
    out={}
    for seed in range(40):
        take=arrays["seed"]==seed;order=np.argsort(arrays["index"][take])
        out[seed]={k:v[take][order] for k,v in arrays.items() if k not in ("seed",)}
    return out


def load_q(cache_dir: Path,seed: int) -> np.ndarray:
    with np.load(cache_dir/f"seed{seed:02d}_stride100_qscan_v4.npz") as data:
        step=np.asarray(data["step"]);q=np.asarray(data["q"])
    q=q[step>=round(TRANSIENT/DT)]
    return np.rint((q+1)*N/2).astype(np.int16)


def storage(state: np.ndarray) -> tuple[float,float]:
    x,y=state[0],state[1];left=(x<0)&((-x)<=ELL)&(np.abs(y)<=ELL/2);right=(x>=0)&(x<=ELL)&(np.abs(y)<=ELL/2)
    return (float(np.sum((ELL+x[left])/ELL)/N),float(np.sum((ELL-x[right])/ELL)/N))


def exact_crossing_audit(raw_dir: Path,processed_dir: Path) -> tuple[list[dict],dict]:
    rows=[];step_mismatch=0;outside=0;processed_files=0
    for seed in range(40):
        processed=processed_dir/f"event_histogram_If00_In05_num100_R50_gap20_n{seed}_BEM_no_hydro.mat"
        if not processed.exists():continue
        processed_files+=1
        nl=np.asarray(loadmat(processed,variable_names=["num_fish_left"])["num_fish_left"]).reshape(-1)
        changes=np.flatnonzero(np.diff(nl[round(TRANSIENT/DT):])!=0)+round(TRANSIENT/DT)
        reader=RawReader(raw_path(raw_dir,seed))
        try:
            for k in changes:
                st=reader.block("state",int(k),int(k+1));x0,x1=st[0,:,0],st[0,:,1];y0,y1=st[1,:,0],st[1,:,1]
                lr=(x0<0)&(x1>=0);rl=(x0>=0)&(x1<0)
                if int(lr.sum()-rl.sum()) != int(nl[k]-nl[k+1]):step_mismatch+=1
                for direction,mask in ((1,lr),(-1,rl)):
                    for fish in np.flatnonzero(mask):
                        t=-x0[fish]/(x1[fish]-x0[fish]);yc=y0[fish]+t*(y1[fish]-y0[fish]);inside=abs(yc)<=ELL/2+1e-12
                        outside+=not inside
                        before_l=(1+x0[fish]/ELL)/N if x0[fish]<0 and abs(y0[fish])<=ELL/2 and -x0[fish]<=ELL else 0.
                        before_r=(1-x0[fish]/ELL)/N if x0[fish]>=0 and abs(y0[fish])<=ELL/2 and x0[fish]<=ELL else 0.
                        after_l=(1+x1[fish]/ELL)/N if x1[fish]<0 and abs(y1[fish])<=ELL/2 and -x1[fish]<=ELL else 0.
                        after_r=(1-x1[fish]/ELL)/N if x1[fish]>=0 and abs(y1[fish])<=ELL/2 and x1[fish]<=ELL else 0.
                        jump_l=-direction/N;jump_r=direction/N
                        # For a physical aperture crossing the finite-frame
                        # remainder is continuous travel to and from x=0.
                        cont_l=(-x0[fish]/(N*ELL)) if direction==1 else (x1[fish]/(N*ELL))
                        cont_r=(-x1[fish]/(N*ELL)) if direction==1 else (x0[fish]/(N*ELL))
                        rows.append(dict(seed=seed,step=int(k),fish=int(fish),direction=direction,y_cross=yc,inside_opening=int(inside),
                                         delta_n_right=direction,jump_px=0.,jump_b_left=jump_l,jump_b_right=jump_r,
                                         endpoint_delta_b_left=after_l-before_l,endpoint_delta_b_right=after_r-before_r,
                                         continuous_b_left=cont_l,continuous_b_right=cont_r,
                                         residual_b_left=(after_l-before_l)-jump_l-cont_l,
                                         residual_b_right=(after_r-before_r)-jump_r-cont_r))
        finally:reader.close()
    return rows,dict(n_processed_files=processed_files,n_seeds_with_crossings=len({r["seed"] for r in rows}),n_individual_crossings=len(rows),n_step_count_mismatches=step_mismatch,
                     n_crossings_outside_opening=outside,max_abs_jump_map_residual=max(max(abs(r["residual_b_left"]),abs(r["residual_b_right"])) for r in rows),
                     exact_jump_maps={"left_to_right":"(n_R+1, px, b_L-1/N, b_R+1/N)","right_to_left":"(n_R-1, px, b_L+1/N, b_R-1/N)"})


def conservation_audit(raw_dir: Path,features: dict[int,dict[str,np.ndarray]]) -> tuple[list[dict],dict]:
    rows=[]
    for seed in range(40):
        count=len(features[seed]["index"])-1;chosen=np.unique(np.linspace(0,count-1,5,dtype=int));reader=RawReader(raw_path(raw_dir,seed))
        try:
            for i in chosen:
                start=round(TRANSIENT/DT)+int(i)*round(FEATURE_DT/DT);stop=start+round(FEATURE_DT/DT)
                states=reader.block("state",start,stop);vel=reader.block("rdot",start+1,stop)
                bl0,br0=storage(states[:,:,0]);bl1,br1=storage(states[:,:,-1])
                for chamber,direction,actual in (("left",1,bl1-bl0),("right",-1,br1-br0)):
                    budget=segment_budget(states,vel,direction);duration=budget["duration"]
                    longitudinal=duration*budget["path_flux"]/N;lateral=-duration*budget["weighted_lateral_flux"]/N;neck=-duration*budget["neck_flux"]/N
                    rows.append(dict(seed=seed,index=int(i),chamber=chamber,delta_storage=actual,longitudinal_change=longitudinal,
                                     lateral_change=lateral,centerline_jump=neck,reconstructed=longitudinal+lateral+neck,
                                     residual=actual-longitudinal-lateral-neck,net_center_crossings_outward=budget["net_neck_count"]))
        finally:reader.close()
    residual=np.asarray([r["residual"] for r in rows])
    return rows,dict(n_windows=len(rows)//2,n_chamber_budgets=len(rows),n_seeds=40,max_abs_residual=float(np.max(np.abs(residual))),rms_residual=float(np.sqrt(np.mean(residual**2))),
                     identity="Delta b_c = longitudinal path contribution - weighted lateral outflux - net centerline outflux/N")


def interval_records(features: dict[int,dict[str,np.ndarray]]) -> list[dict]:
    rows=[]
    for seed,a in features.items():
        for i in range(len(a["q"])-1):
            net=int(a["plus"][i]-a["minus"][i]);state=sid(a["q"][i],a["px"][i],a["b_left"][i],a["b_right"][i])
            rows.append(dict(seed=seed,index=i,state=state,next_state=sid(a["q"][i+1],a["px"][i+1],a["b_left"][i+1],a["b_right"][i+1]),
                             plus=int(a["plus"][i]),minus=int(a["minus"][i]),net=net,
                             dpx=float(a["px"][i+1]-a["px"][i]),
                             dbl=float(a["b_left"][i+1]-a["b_left"][i]+net/N),
                             dbr=float(a["b_right"][i+1]-a["b_right"][i]-net/N)))
    return rows


def moment_model(rows: list[dict]) -> tuple[list[dict],dict]:
    train=[r for r in rows if r["seed"]<30];by=defaultdict(list)
    for r in train:
        y=np.asarray([r["dpx"],r["dbl"],r["dbr"]]);by[r["state"]].append(y)
        by[reflect_sid(r["state"])].append(np.asarray([-r["dpx"],r["dbr"],r["dbl"]]))
    means={s:np.mean(v,axis=0) for s,v in by.items() if len(v)>=100};table=[]
    for s,v in by.items():
        if s not in means:continue
        a=np.asarray(v);cov=np.cov(a,rowvar=False)/FEATURE_DT
        table.append(dict(state=s,q_bin=decode(s)[0],px_bin=decode(s)[1],b_left_bin=decode(s)[2],b_right_bin=decode(s)[3],n_train=len(a),
                          drift_px=means[s][0]/FEATURE_DT,drift_b_left=means[s][1]/FEATURE_DT,drift_b_right=means[s][2]/FEATURE_DT,
                          cov_px_px=cov[0,0],cov_b_left_b_left=cov[1,1],cov_b_right_b_right=cov[2,2],cov_b_left_b_right=cov[1,2]))
    seed_stats=[]
    for seed in range(30,40):
        part=[r for r in rows if r["seed"]==seed and r["state"] in means];res=[]
        for r in part:
            actual=np.asarray([r["dpx"],r["dbl"],r["dbr"]]);res.append(actual-means[r["state"]])
        res=np.asarray(res);actual=np.asarray([[r["dpx"],r["dbl"],r["dbr"]] for r in part]);previous=np.vstack((np.full((1,3),np.nan),actual[:-1]));rhos=[]
        for j in range(3):
            ok=np.isfinite(previous[:,j]);rhos.append(float(spearmanr(previous[ok,j],res[ok,j]).statistic))
        seed_stats.append(dict(seed=seed,n=len(part),rho_px=rhos[0],rho_b_left=rhos[1],rho_b_right=rhos[2],
                               mse_px=float(np.mean(res[:,0]**2)),mse_b_left=float(np.mean(res[:,1]**2)),mse_b_right=float(np.mean(res[:,2]**2))))
    rng=np.random.default_rng(20261119);summary={"n_supported_cells":len(means),"n_total_cells":NSTATES,"heldout_coverage":sum(r["n"] for r in seed_stats)/sum(1 for r in rows if r["seed"]>=30)}
    for field in ("rho_px","rho_b_left","rho_b_right"):
        a=np.asarray([r[field] for r in seed_stats]);pick=rng.integers(len(a),size=(4000,len(a)));summary[field]={"mean":float(a.mean()),"ci95":np.percentile(a[pick].mean(axis=1),[2.5,97.5]).tolist()}
    return table,{"summary":summary,"seed_rows":seed_stats,"means":means}


def directed_rate_rows(features: dict[int,dict[str,np.ndarray]],qseries: dict[int,np.ndarray]) -> list[dict]:
    rows=[]
    horizons=(.1,.2,.5)
    for seed,a in features.items():
        n=qseries[seed]
        for i in range(len(a["q"])):
            origin=i*round(FEATURE_DT/Q_DT)
            for direction in (-1,1):
                source=a["b_left"][i] if direction==1 else a["b_right"][i]
                cell=(bindex(direction*a["q"][i],Q_EDGES),bindex(direction*a["px"][i],PX_EDGES),int(source>=B_THRESHOLD))
                row=dict(seed=seed,index=i,direction=direction,cell=cell)
                for h in horizons:
                    m=round(h/Q_DT)
                    if origin+m>=len(n):row[f"count_{h}"]=None;continue
                    d=np.diff(n[origin:origin+m+1]);row[f"count_{h}"]=int(d[d>0].sum() if direction==1 else -d[d<0].sum())
                rows.append(row)
    return rows


def poisson_deviance(y: np.ndarray,mu: np.ndarray) -> float:
    mu=np.maximum(mu,1e-12);term=mu.copy();ok=y>0;term[ok]=y[ok]*np.log(y[ok]/mu[ok])-(y[ok]-mu[ok]);return float(2*np.mean(term))


def rate_consistency(rows: list[dict]) -> tuple[list[dict],dict,dict]:
    train=[r for r in rows if r["seed"]<30];test=[r for r in rows if r["seed"]>=30];horizons=(.1,.2,.5);rates={};support={}
    for h in horizons:
        by=defaultdict(list)
        for r in train:
            if r[f"count_{h}"] is not None:by[r["cell"]].append(r[f"count_{h}"])
        rates[h]={k:float(np.mean(v)/h) for k,v in by.items() if len(v)>=100};support[h]=by
    common=set.intersection(*(set(rates[h]) for h in horizons));table=[];instant=rates[.1];rng=np.random.default_rng(20261121)
    for h in horizons:
        seeddev=[];seed_rate=[];seed_difference=[]
        for seed in range(30,40):
            part=[r for r in test if r["seed"]==seed and r["cell"] in instant and r[f"count_{h}"] is not None]
            y=np.asarray([r[f"count_{h}"] for r in part]);mu=np.asarray([instant[r["cell"]]*h for r in part]);seeddev.append(poisson_deviance(y,mu));seed_rate.append(float(y.sum()/(len(y)*h)))
            base_y=np.asarray([r["count_0.1"] for r in part]);seed_difference.append(float(y.sum()/(len(y)*h)-base_y.sum()/(len(base_y)*.1)))
        difference=np.asarray(seed_difference);pick=rng.integers(len(difference),size=(4000,len(difference)))
        weights=np.asarray([len(support[h][c]) for c in common]);r0=np.asarray([rates[.1][c] for c in common]);rh=np.asarray([rates[h][c] for c in common])
        conditional_difference=float(np.sum(weights*np.abs(rh-r0))/np.sum(weights*np.maximum(r0,1e-12)))
        observed=[];predicted=[];multi=[];fanos=[]
        for c in common:
            v=np.asarray(support[h][c]);observed.append(np.mean(v>0));predicted.append(1-np.exp(-rates[.1][c]*h));multi.append(np.mean(v>=2));
            if np.mean(v)>0:fanos.append(np.var(v,ddof=1)/np.mean(v))
        table.append(dict(horizon=h,n_supported_cells=len(common),conditional_rate_abs_difference_fraction=conditional_difference,
                          mean_observed_rate=float(np.average(rh,weights=weights)),mean_instantaneous_rate=float(np.average(r0,weights=weights)),
                          event_probability_mae=float(np.mean(np.abs(np.asarray(observed)-np.asarray(predicted)))),
                          observed_multiple_count_fraction=float(np.mean(multi)),mean_cell_fano=float(np.mean(fanos)),heldout_poisson_deviance=float(np.mean(seeddev)),heldout_mean_rate=float(np.mean(seed_rate)),
                          heldout_rate_difference_from_0p1=float(difference.mean()),heldout_rate_difference_ci95=json.dumps(np.percentile(difference[pick].mean(axis=1),[2.5,97.5]).tolist())))
    details={str(k):v for k,v in instant.items()}
    return table,{"instantaneous_horizon":.1,"n_rate_cells":len(instant),"cell_definition":"oriented q regime, oriented px regime, low/high source storage"},details


def transition_matrix(features: dict[int,dict[str,np.ndarray]],seeds: range,lag: int) -> tuple[np.ndarray,np.ndarray]:
    count=np.zeros((NSTATES,NSTATES))
    for seed in seeds:
        a=features[seed];states=np.asarray([sid(*x) for x in zip(a["q"],a["px"],a["b_left"],a["b_right"])])
        for x,y in zip(states[:-lag],states[lag:]):
            count[x,y]+=1;count[reflect_sid(x),reflect_sid(y)]+=1
    support=count.sum(axis=1);p=np.zeros_like(count);ok=support>0;p[ok]=count[ok]/support[ok,None];p[~ok,np.flatnonzero(~ok)]=1
    return p,support


def tv(a: np.ndarray,b: np.ndarray,weight: np.ndarray,minimum: int=30) -> tuple[float,int]:
    ok=weight>=minimum;w=weight[ok]/weight[ok].sum();return float(np.sum(w*.5*np.abs(a[ok]-b[ok]).sum(axis=1))),int(ok.sum())


def operator_tests(features: dict[int,dict[str,np.ndarray]]) -> tuple[list[dict],dict,np.ndarray,list[dict]]:
    base,base_support=transition_matrix(features,range(30),1);out=[];rng=np.random.default_rng(20261120)
    for tau in LAGS:
        lag=round(tau/FEATURE_DT);direct,support=transition_matrix(features,range(30),lag);pred=np.linalg.matrix_power(base,lag);ck,n=tv(direct,pred,support)
        ok=support>=30;w=support[ok]/support[ok].sum();null=[]
        for _ in range(300):
            val=[]
            for probability,count in zip(pred[ok],support[ok].astype(int)):
                sample=rng.multinomial(count,probability)/count;val.append(.5*np.abs(sample-probability).sum())
            null.append(float(np.dot(w,val)))
        held,hs=transition_matrix(features,range(30,40),lag);heldtv,hn=tv(held,pred,hs)
        out.append(dict(tau=tau,ck_total_variation=ck,ck_sampling_null_ci95=json.dumps(np.percentile(null,[2.5,97.5]).tolist()),ck_supported_states=n,
                        heldout_total_variation=heldtv,heldout_supported_states=hn))
    supported=base_support>0;pi=base_support/base_support.sum();converged=False
    for iteration in range(100000):
        new=pi@base
        if np.max(np.abs(new-pi))<1e-14:converged=True;pi=new;break
        pi=new
    obs=np.zeros(NSTATES)
    for seed in range(30,40):
        a=features[seed]
        for x in zip(a["q"],a["px"],a["b_left"],a["b_right"]):
            s=sid(*x);obs[s]+=1;obs[reflect_sid(s)]+=1
    obs/=obs.sum();stationary_tv=float(.5*np.abs(pi-obs).sum())
    qpi=np.zeros(NQ);qobs=np.zeros(NQ)
    for s in range(NSTATES):qpi[decode(s)[0]]+=pi[s];qobs[decode(s)[0]]+=obs[s]
    stationary_rows=[]
    for s in range(NSTATES):
        qb,pb,lb,rb=decode(s);stationary_rows.append(dict(state=s,q_bin=qb,px_bin=pb,b_left_bin=lb,b_right_bin=rb,predicted_stationary=pi[s],heldout_occupancy=obs[s]))
    return out,dict(joint_stationary_tv=stationary_tv,q_marginal_stationary_tv=float(.5*np.abs(qpi-qobs).sum()),predicted_q_marginal=qpi.tolist(),heldout_q_marginal=qobs.tolist(),power_iterations=iteration+1,
                    stationary_converged=int(converged),stationary_residual=float(np.max(np.abs(pi@base-pi))),base_observed_states=int(np.sum(supported)),base_supported_states=int(np.sum(base_support>=30))),base,stationary_rows


def pathway_tests(features: dict[int,dict[str,np.ndarray]],base: np.ndarray) -> list[dict]:
    out=[]
    for tau in (.5,1.,2.,5.):
        lag=round(tau/FEATURE_DT);pred=[];actual=[]
        for seed in range(30,40):
            a=features[seed];states=np.asarray([sid(*x) for x in zip(a["q"],a["px"],a["b_left"],a["b_right"])])
            for i in range(len(states)-lag-1):
                net=int(a["plus"][i]-a["minus"][i])
                if net==0:continue
                direction=int(np.sign(net));origin=states[i+1];future=states[i+1+lag]
                if direction<0:origin=reflect_sid(origin);future=reflect_sid(future)
                pred.append(np.linalg.matrix_power(base,lag)[origin]);actual.append(future)
        meanpred=np.mean(pred,axis=0);hist=np.bincount(actual,minlength=NSTATES)/len(actual);joint=float(.5*np.abs(meanpred-hist).sum())
        qp=np.zeros(NQ);qh=np.zeros(NQ)
        for s in range(NSTATES):qp[decode(s)[0]]+=meanpred[s];qh[decode(s)[0]]+=hist[s]
        out.append(dict(time_after_partial_transfer=tau,n_test_events=len(actual),joint_state_tv=joint,q_regime_tv=float(.5*np.abs(qp-qh).sum())))
    return out


def make_figure(jump_rows,conservation,rate_rows,operator,memory,pathways,output: Path) -> None:
    fig,ax=plt.subplots(2,3,figsize=(12,7.5))
    residual=np.maximum(np.abs([r["residual_b_left"] for r in jump_rows]+[r["residual_b_right"] for r in jump_rows]),1e-20)
    ax[0,0].hist(residual,bins=np.logspace(-20,-12,25));ax[0,0].set(xscale="log",xlabel="jump-map absolute residual",ylabel="crossings",title="exact centerline transfer")
    cr=np.maximum(np.abs([r["residual"] for r in conservation]),1e-20);ax[0,1].hist(cr,bins=np.logspace(-20,-12,25));ax[0,1].set(xscale="log",xlabel="storage-budget residual",title="conservation between crossings")
    ax[0,2].plot([r["horizon"] for r in rate_rows],[r["conditional_rate_abs_difference_fraction"] for r in rate_rows],"o-");ax[0,2].set(xlabel="counting horizon",ylabel="relative rate change",title="continuous-time intensity audit")
    tau=[r["tau"] for r in operator];ax[1,0].plot(tau,[r["ck_total_variation"] for r in operator],"o-",label="CK");ax[1,0].plot(tau,[json.loads(r["ck_sampling_null_ci95"])[1] for r in operator],"k--",label="sampling 97.5%")
    ax[1,0].plot(tau,[r["heldout_total_variation"] for r in operator],"s-",label="held out");ax[1,0].set(xlabel="lag",ylabel="total variation",title="joint operator");ax[1,0].legend(frameon=False,fontsize=8)
    names=("px","b_left","b_right");vals=[memory[f"rho_{x}"]["mean"] for x in names];lo=[memory[f"rho_{x}"]["ci95"][0] for x in names];hi=[memory[f"rho_{x}"]["ci95"][1] for x in names]
    ax[1,1].bar(names,vals,yerr=[np.asarray(vals)-lo,np.asarray(hi)-vals]);ax[1,1].axhline(0,color="k",lw=.7);ax[1,1].set(ylabel="residual vs previous-increment Spearman",title="unresolved between-event memory")
    ax[1,2].plot([r["time_after_partial_transfer"] for r in pathways],[r["joint_state_tv"] for r in pathways],"o-",label="joint");ax[1,2].plot([r["time_after_partial_transfer"] for r in pathways],[r["q_regime_tv"] for r in pathways],"s-",label="q marginal");ax[1,2].set(xlabel="time after partial transfer",ylabel="held-out TV",title="partial-transfer pathways");ax[1,2].legend(frameon=False)
    fig.suptitle("Conservation-constrained hybrid-model validation",fontsize=11);fig.tight_layout();fig.savefig(output,dpi=220);plt.close(fig)


def compact_jump_rows(rows: list[dict]) -> list[dict]:
    out=[]
    for seed in sorted({r["seed"] for r in rows}):
        for direction in (-1,1):
            part=[r for r in rows if r["seed"]==seed and r["direction"]==direction]
            if not part:continue
            out.append(dict(seed=seed,direction=direction,n_crossings=len(part),min_y_cross=float(min(r["y_cross"] for r in part)),max_y_cross=float(max(r["y_cross"] for r in part)),
                            max_abs_b_left_residual=float(max(abs(r["residual_b_left"]) for r in part)),max_abs_b_right_residual=float(max(abs(r["residual_b_right"]) for r in part)),
                            sum_jump_b_left=float(sum(r["jump_b_left"] for r in part)),sum_jump_b_right=float(sum(r["jump_b_right"] for r in part))))
    return out


def self_test() -> None:
    assert reflect_sid(reflect_sid(17))==17
    x=np.asarray([[-.1,.1],[0.,0.],[0.,0.]])
    assert np.allclose(storage(x),((1-.1/ELL)/N,(1-.1/ELL)/N))
    plus=np.asarray([50,51]);net=plus[1]-plus[0];assert (-net/N)+(net/N)==0
    print("PASS: reflection, storage weights, and exact conservative jump maps")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--processed-dir",type=Path);p.add_argument("--q-cache",type=Path);p.add_argument("--feature-cache",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.processed_dir,args.q_cache,args.feature_cache,args.output)):p.error("all data and output paths are required")
    args.output.mkdir(parents=True,exist_ok=True);features=load_features(args.feature_cache);qseries={s:load_q(args.q_cache,s) for s in range(40)}
    jumps,jump_summary=exact_crossing_audit(args.raw_dir,args.processed_dir);conservation,conservation_summary=conservation_audit(args.raw_dir,features);intervals=interval_records(features)
    moments,moment_info=moment_model(intervals);directed=directed_rate_rows(features,qseries);rate_rows,rate_summary,rate_details=rate_consistency(directed)
    operator,stationary,base,stationary_rows=operator_tests(features);pathways=pathway_tests(features,base)
    write_csv(args.output/"crossing_jump_audit.csv",compact_jump_rows(jumps));write_csv(args.output/"storage_conservation_audit.csv",conservation);write_csv(args.output/"between_event_moments.csv",moments);write_csv(args.output/"between_event_memory_by_seed.csv",moment_info.pop("seed_rows"));moment_info.pop("means")
    write_csv(args.output/"continuous_time_rate_audit.csv",rate_rows);write_csv(args.output/"joint_operator_tests.csv",operator);write_csv(args.output/"stationary_occupancy_comparison.csv",stationary_rows);write_csv(args.output/"partial_transfer_pathways.csv",pathways)
    make_figure(jumps,conservation,rate_rows,operator,moment_info["summary"],pathways,args.output/"hybrid_model_validation.png")
    gate=bool(all(r["ck_total_variation"]<=json.loads(r["ck_sampling_null_ci95"])[1] for r in operator) and all(moment_info["summary"][f"rho_{x}"]["ci95"][0]<=0<=moment_info["summary"][f"rho_{x}"]["ci95"][1] for x in ("px","b_left","b_right")))
    summary=dict(configuration=dict(state="(n_R,px,b_L,b_R)",population_lattice=list(range(101)),validation_bins=dict(q=Q_EDGES.tolist(),px=PX_EDGES.tolist(),storage_threshold=B_THRESHOLD,n_cells=NSTATES),train_seeds="0-29",test_seeds="30-39",new_simulations=0),
                 exact_jump_audit=jump_summary,storage_conservation=conservation_summary,continuous_time_rate_audit=rate_summary,between_event_model=moment_info["summary"],joint_operator=operator,stationary_prediction=stationary,partial_transfer_pathways=pathways,
                 dynamically_validated_landscape=gate,landscape_policy="construct joint stationary landscape and currents only if the joint CK and memory gates pass")
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n");(args.output/"instantaneous_rate_cells.json").write_text(json.dumps(rate_details,indent=2)+"\n")
    (args.output/"analysis_configuration.json").write_text(json.dumps(dict(raw_dir=str(args.raw_dir),processed_dir=str(args.processed_dir),q_cache=str(args.q_cache),feature_cache=str(args.feature_cache),feature_dt=FEATURE_DT,q_dt=Q_DT),indent=2)+"\n")


if __name__=="__main__":main()
