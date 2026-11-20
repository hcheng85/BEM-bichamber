#!/usr/bin/env python3
"""Test adiabatic elimination of directed relative-heading dynamics.

Uses the actual visual graph and normalized production coefficients.  No drift,
potential, or empirical closure is fit.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_approach_current_balance import DT, N, fast_check_voronoi, load_gap20_production, replay_torques, raw_path
from analyze_geometric_first_passage_reduction import geometry, wrap
from analyze_initiation_propagation import RawReader


IP=4.;LEADS=(1.,.8,.6,.4,.2,.1)


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def cohort(path: Path) -> list[dict]:
    out=[]
    for r in read_csv(path):out.append(dict(seed=int(r["seed"]),event_id=int(r["event_id"]),first_native=int(r["first_native"]),source_sign=int(r["source_sign"]),success4=int(r["success4"]),abs_q_pre=float(r["abs_q_pre"])))
    if len(out)!=337:raise RuntimeError("Expected fixed 337-event cohort")
    return out


def source_mask(state: np.ndarray,sign: int):return state[0]>=0 if sign>0 else state[0]<0


def directed_snapshot(mod,para: dict,state: np.ndarray,sign: int) -> tuple[dict,np.ndarray,np.ndarray]:
    mask=source_mask(state,sign);idx=np.flatnonzero(mask);ns=len(idx);local={g:i for i,g in enumerate(idx)};z=geometry(state,mask);vn=fast_check_voronoi(mod,mod.GetVoronoi(state[:2],para),state[:2]);theta=state[2];e=np.column_stack((np.cos(theta),np.sin(theta)));visual=np.zeros((N,N));align=np.zeros((N,N));position=np.zeros((N,N))
    for i,j in zip(*np.nonzero(vn)):
        rho=state[:2,j]-state[:2,i];dist=float(np.linalg.norm(rho));visual[i,j]=1+float(e[i]@(rho/dist));align[i,j]=IP*math.sin(theta[j]-theta[i]);position[i,j]=e[i,0]*rho[1]-e[i,1]*rho[0]
    denom=visual.sum(axis=1);alpha=np.divide(visual,denom[:,None],out=np.zeros_like(visual),where=denom[:,None]>0)
    b=np.zeros((ns,ns));internal_exact=np.zeros(ns);internal_linear=np.zeros(ns);external=np.zeros(ns);pos=np.zeros(ns)
    delta=np.asarray([float(wrap(theta[i]-z["theta"])) for i in idx])
    for gi in idx:
        i=local[gi]
        for gj in np.flatnonzero(vn[gi]):
            a=alpha[gi,gj]
            if mask[gj]:
                j=local[gj];b[i,j]=a;internal_exact[i]+=a*align[gi,gj];internal_linear[i]+=a*IP*(delta[j]-delta[i])
            else:external[i]+=a*align[gi,gj]
            pos[i]+=a*position[gi,gj]
    lap=np.diag(b.sum(axis=1))-b;eigvals=np.linalg.eigvals(lap);order=np.argsort(np.abs(eigvals));zero_count=int(np.sum(np.abs(eigvals)<1e-9));positive=[v.real for v in eigvals if v.real>1e-9];gap=float(min(positive)) if positive else np.nan
    vals,vecs=np.linalg.eig(lap.T);pi=np.real(vecs[:,np.argmin(np.abs(vals))]);pi*=1 if pi.sum()>=0 else -1;pi/=pi.sum()
    wall,social,clip,_=replay_torques(mod,state,para);forcing=np.asarray([wall[i]+pos[k]+external[k]+clip[i] for k,i in enumerate(idx)])
    effective=np.asarray([wall[i]+social[i]+clip[i] for i in idx]);c=np.cos(delta)/(ns*z["p"]);actual=float(c@effective);linear_actual=float(c@(internal_linear+forcing));pi_rate=float(pi@forcing);uniform_rate=float(np.mean(forcing))
    qforcing=forcing-pi_rate;aug=lap+np.ones((ns,1))*pi[None,:]
    try:quasi=np.linalg.solve(aug,qforcing)/IP
    except np.linalg.LinAlgError:quasi=np.full(ns,np.nan)
    delta_pi=delta-float(pi@delta);ok=np.isfinite(quasi);corr=float(np.corrcoef(delta_pi[ok],quasi[ok])[0,1]) if ok.sum()>2 and np.std(quasi[ok])>0 else np.nan
    row=dict(n_source=ns,polarization=z["p"],psi=z["psi"],spectral_gap=gap,relaxation_time=1/(IP*gap) if gap>0 else np.inf,zero_eigenvalues=zero_count,pi_effective_number=1/float(pi@pi),pi_l1_from_uniform=float(np.sum(np.abs(pi-1/ns))),
        internal_exact_collective=float(c@internal_exact),internal_linear_collective=float(c@internal_linear),linearization_error=float(c@(internal_exact-internal_linear)),actual_collective_rate=actual,linearized_collective_rate=linear_actual,
        directed_adiabatic_rate=pi_rate,uniform_forcing_rate=uniform_rate,quasi_relative_correlation=corr,quasi_relative_rmse=float(np.sqrt(np.mean((quasi-delta_pi)**2))),relative_heading_sd=float(np.std(delta_pi)),forcing_sd=float(np.std(forcing)),
        predicted_sym_minus_baseline_alignment=float(-math.sin(z["psi"])*(uniform_rate-pi_rate)))
    return row,pi,vn[np.ix_(idx,idx)]


def extract(mod,para: dict,raw_dir: Path,events: list[dict]) -> tuple[list[dict],dict]:
    rows=[];by=defaultdict(list);lookup={}
    for ev in events:by[ev["seed"]].append(ev)
    for count,seed in enumerate(sorted(by),1):
        print(f"directed elimination audit {count}/{len(by)} seed={seed}",flush=True);rd=RawReader(raw_path(raw_dir,seed))
        try:
            for ev in by[seed]:
                for lead in LEADS:
                    state=rd.frame("state",ev["first_native"]-round(lead/DT));row,pi,graph=directed_snapshot(mod,para,state,ev["source_sign"]);item=dict(**ev,lead=lead,**row);rows.append(item);lookup[(ev["seed"],ev["event_id"],lead)]=(item,pi,graph,state)
        finally:rd.close()
    path=[]
    for ev in events:
        seq=[lookup[(ev["seed"],ev["event_id"],lead)] for lead in LEADS];times=np.asarray([1-lead for lead in LEADS]);rates=np.asarray([x[0]["directed_adiabatic_rate"] for x in seq]);uniform=np.asarray([x[0]["uniform_forcing_rate"] for x in seq]);theta0=geometry(seq[0][3],source_mask(seq[0][3],ev["source_sign"]))["theta"];theta1=geometry(seq[-1][3],source_mask(seq[-1][3],ev["source_sign"]))["theta"];actual=float(wrap(theta1-theta0));pi0,pi1=seq[0][1],seq[-2][1];g0,g1=seq[0][2],seq[-2][2]
        path.append(dict(**ev,actual_heading_change_1_to_0p1=actual,directed_integrated_response=float(np.trapezoid(rates,times)),uniform_integrated_response=float(np.trapezoid(uniform,times)),frozen_directed_0p5=.5*rates[0],frozen_uniform_0p5=.5*uniform[0],
            pi_l1_change_1_to_0p2=float(np.sum(np.abs(pi1-pi0))),graph_edge_jaccard_1_to_0p2=float(np.sum(g0&g1)/np.sum(g0|g1))))
    audit=dict(n_snapshots=len(rows),n_paths=len(path),max_abs_linearization_error=max(abs(r["linearization_error"]) for r in rows),disconnected_snapshots=sum(r["zero_eigenvalues"]!=1 for r in rows))
    return rows,path,audit


def cluster_boot(rows: list[dict],metric,seed=2026112401,nboot=3000) -> tuple[float,list[float]]:
    by=defaultdict(list)
    for r in rows:by[r["seed"]].append(r)
    seeds=list(by);point=float(metric(rows));rng=np.random.default_rng(seed);draw=[]
    for _ in range(nboot):
        sample=[]
        for rep,s in enumerate(rng.choice(seeds,len(seeds),replace=True)):
            for r in by[int(s)]:z=dict(r);z["seed"]=rep;sample.append(z)
        v=metric(sample)
        if np.isfinite(v):draw.append(v)
    return point,np.percentile(draw,[2.5,97.5]).tolist()


def prediction_stats(rows: list[dict],actual: str,pred: str,intervals: bool=False) -> dict:
    def values(items):
        a=np.asarray([r[actual] for r in items]);p=np.asarray([r[pred] for r in items])
        return float(np.corrcoef(a,p)[0,1]),float(np.sqrt(np.mean((a-p)**2)))
    a=np.asarray([r[actual] for r in rows]);p=np.asarray([r[pred] for r in rows]);corr,rmse=values(rows)
    out=dict(correlation=corr,rmse=rmse,actual_sd=float(np.std(a)),bias=float(np.mean(p-a)),n=len(rows))
    if intervals:
        out["correlation_seed_interval"]=cluster_boot(rows,lambda x:values(x)[0])[1]
        out["rmse_seed_interval"]=cluster_boot(rows,lambda x:values(x)[1],seed=2026112402)[1]
    return out


def summarize(samples: list[dict],paths: list[dict],branch_effects: Path|None) -> dict:
    result={"timescale":{},"instantaneous":{},"integrated":{}}
    for name,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
        s=[r for r in samples if pred(r["seed"])];lead1=[r for r in s if r["lead"]==1.];p=[r for r in paths if pred(r["seed"])]
        tau=np.asarray([r["relaxation_time"] for r in s]);result["timescale"][name]=dict(median_relaxation_time=float(np.median(tau)),q10_q90=np.quantile(tau,[.1,.9]).tolist(),fraction_below_0p1=float(np.mean(tau<.1)),fraction_below_0p5=float(np.mean(tau<.5)),fraction_below_1=float(np.mean(tau<1)),median_pi_l1_change=float(np.median([r["pi_l1_change_1_to_0p2"] for r in p])),median_edge_jaccard=float(np.median([r["graph_edge_jaccard_1_to_0p2"] for r in p])),median_quasi_relative_correlation=float(np.nanmedian([r["quasi_relative_correlation"] for r in s])),median_quasi_rmse_over_sd=float(np.nanmedian([r["quasi_relative_rmse"]/r["relative_heading_sd"] for r in s])))
        with_intervals=name=="heldout"
        result["instantaneous"][name]={"directed":prediction_stats(lead1,"actual_collective_rate","directed_adiabatic_rate",with_intervals),"uniform":prediction_stats(lead1,"actual_collective_rate","uniform_forcing_rate",with_intervals),"linearized_exact_state":prediction_stats(lead1,"actual_collective_rate","linearized_collective_rate",with_intervals)}
        result["integrated"][name]={"directed_path_quadrature":prediction_stats(p,"actual_heading_change_1_to_0p1","directed_integrated_response",with_intervals),"uniform_path_quadrature":prediction_stats(p,"actual_heading_change_1_to_0p1","uniform_integrated_response",with_intervals)}
    if branch_effects:
        effects={int(r["seed"]):float(r["turning_alignment_rate_1"]) for r in read_csv(branch_effects)}
        manifest=branch_effects.with_name("selection_manifest.csv")
        selected={(int(r["seed"]),int(r["event_id"])) for r in read_csv(manifest)}
        initial=[]
        for r in samples:
            if r["lead"]==1. and (r["seed"],r["event_id"]) in selected:initial.append(dict(seed=r["seed"],predicted=r["predicted_sym_minus_baseline_alignment"],observed=effects[r["seed"]]))
        stats=prediction_stats(initial,"observed","predicted",True)
        pred=np.asarray([x["predicted"] for x in initial]);actual=np.asarray([x["observed"] for x in initial]);stats.update(predicted_mean=float(np.mean(pred)),observed_mean=float(np.mean(actual)),sign_agreement=float(np.mean(np.sign(pred)==np.sign(actual))))
        result["causal_intervention_prediction"]=stats
    return result


def figures(samples,paths,summary,out: Path) -> None:
    fig,axes=plt.subplots(1,3,figsize=(12,3.8));held=[r for r in samples if r["seed"]>=30 and r["lead"]==1.]
    axes[0].hist([r["relaxation_time"] for r in held],bins=20);axes[0].axvline(.5,color="r",ls=":",label="0.5 approach scale");axes[0].set(xlabel="directed relative-heading relaxation time",title="Held-out spectral times");axes[0].legend(frameon=False)
    a=np.asarray([r["actual_collective_rate"] for r in held]);p=np.asarray([r["directed_adiabatic_rate"] for r in held]);axes[1].scatter(a,p,s=12,alpha=.6);lim=[min(a.min(),p.min()),max(a.max(),p.max())];axes[1].plot(lim,lim,"k:");axes[1].set(xlabel="actual instantaneous turn",ylabel="pi-weighted forcing",title="Adiabatic instantaneous test")
    hp=[r for r in paths if r["seed"]>=30];a=np.asarray([r["actual_heading_change_1_to_0p1"] for r in hp]);p=np.asarray([r["directed_integrated_response"] for r in hp]);axes[2].scatter(a,p,s=12,alpha=.6);lim=[min(a.min(),p.min()),max(a.max(),p.max())];axes[2].plot(lim,lim,"k:");axes[2].set(xlabel="actual heading change",ylabel="integrated pi-weighted forcing",title="Approach-integrated test");fig.tight_layout();fig.savefig(out/"directed_elimination_validation.png",dpi=220);plt.close(fig)


def self_test() -> None:
    b=np.asarray([[0.,.8,.2],[.1,0.,.9],[.4,.6,0.]]);l=np.diag(b.sum(axis=1))-b;vals,vecs=np.linalg.eig(l.T);pi=np.real(vecs[:,np.argmin(abs(vals))]);pi*=np.sign(pi.sum());pi/=pi.sum();assert np.max(abs(l@np.ones(3)))<1e-14 and np.max(abs(pi@l))<1e-14 and np.max(abs(pi-np.ones(3)/3))>1e-2
    print("PASS: directed Laplacian has distinct right and left null vectors")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--cohort",type=Path);p.add_argument("--source",type=Path);p.add_argument("--branch-effects",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.raw_dir,a.cohort,a.source,a.output)):p.error("all primary paths required")
    a.output.mkdir(parents=True,exist_ok=True);events=cohort(a.cohort)
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(a.source)
    para=dict(Ip=IP,In=.5,R=mod.R,light=1.,dt=DT,Iw=.94,If=0.,num=N,delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
    samples,paths,audit=extract(mod,para,a.raw_dir,events);summary=summarize(samples,paths,a.branch_effects);config=dict(leads=list(LEADS),n_events=len(events),train_seeds="0-29",heldout_seeds="30-39",new_simulations=0,fitted_terms=0,production_source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest())
    write_csv(a.output/"directed_operator_samples.csv",samples);write_csv(a.output/"integrated_path_tests.csv",paths);figures(samples,paths,summary,a.output);(a.output/"summary.json").write_text(json.dumps(dict(configuration=config,audits=audit,results=summary),indent=2)+"\n")


if __name__=="__main__":main()
