#!/usr/bin/env python3
"""Parameter-free straight-line aperture capture and persistence audit.

The simulator and archived trajectories are immutable inputs.  The aperture
test uses only instantaneous position, heading, and the physical gap geometry.
Statistical models are avoided: evaluation uses ranks, exact-q pairs, and
conditional concordance inside the previously declared (q,a,b) cells.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtr

from analyze_approach_current_balance import DT, ELL, N, raw_path
from analyze_collective_qab_pathway import A_EDGES, Q_EDGES, features
from analyze_initiation_propagation import RawReader, auc_score


HALF_WIDTH=ELL/2
IN_NOISE=0.5
EVENT_LEADS=(3.0,2.0,1.0,0.5,0.2)
Q_STANDARD_EDGES=np.asarray([0.0,0.5,0.8,1.000001])
PERSISTENCE_EDGES=np.asarray([0.0,0.25,0.5,1.0,2.0,4.1])


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def wrap(angle: np.ndarray|float) -> np.ndarray|float:
    return np.arctan2(np.sin(angle),np.cos(angle))


def capture_geometry(state: np.ndarray,source_sign: int) -> dict:
    """Exact ray-opening intersection for the source-side approach rectangle.

    source_sign is +1 for the right source and -1 for the left source.  Local
    forward direction s=-source_sign makes z=s*x increase toward the divider.
    A ray is captured iff u=s*cos(theta)>0 and its divider intercept lies in
    |y_hit|<=ELL/2.  Because both endpoints are inside the lateral interval,
    the straight segment cannot exit laterally before reaching the divider.
    """
    s=-int(source_sign);x,y,theta=state
    z=s*x;u=s*np.cos(theta);vy=np.sin(theta)
    inside=(z<0)&(z>=-ELL)&(np.abs(y)<=HALF_WIDTH)
    inward=inside&(u>0)
    tau=np.full(N,np.nan);y_hit=np.full(N,np.nan)
    tau[inward]=-z[inward]/u[inward]
    y_hit[inward]=y[inward]+vy[inward]*tau[inward]
    captured=inward&(np.abs(y_hit)<=HALF_WIDTH)
    j_inward=float(np.sum(u[inward])/N)
    j_capture=float(np.sum(u[captured])/N)
    return dict(mask=captured,inside=inside,u=u,tau=tau,y_hit=y_hit,
                j_inward=j_inward,j_capture=j_capture,n_capture=int(captured.sum()),
                capture_weight=float(np.sum(u[captured])),
                mean_capture_time=float(np.average(tau[captured],weights=u[captured])) if captured.any() else float("nan"))


def snapshot_metrics(state: np.ndarray,velocity: np.ndarray,source_sign: int) -> dict:
    cap=capture_geometry(state,source_sign);s=-source_sign
    ja=float(np.sum(s*velocity[0]*cap["inside"])/N)
    f=features(state,source_sign)
    return dict(abs_q=abs(float(f["q"])),a=float(f["a"]),b=float(f["b"]),
                Ja=ja,Jin=cap["j_inward"],Jcap=cap["j_capture"],n_capture=cap["n_capture"],
                mean_capture_time=cap["mean_capture_time"]),cap


def load_frames(raw_dir: Path,requests: dict[int,set[int]]) -> dict[tuple[int,int],tuple[np.ndarray,np.ndarray]]:
    out={}
    for seed in sorted(requests):
        reader=RawReader(raw_path(raw_dir,seed))
        try:
            for step in sorted(requests[seed]):
                out[(seed,step)]=(reader.frame("state",step),reader.frame("rdot",step+1))
        finally:reader.close()
    return out


def qbin(q: float) -> int|None:
    i=int(np.searchsorted(Q_EDGES,abs(q),side="right")-1)
    return i if 0<=i<len(Q_EDGES)-1 else None


def abin(a: float) -> int|None:
    i=int(np.searchsorted(A_EDGES,a,side="right")-1)
    return i if 0<=i<len(A_EDGES)-1 else None


def fit_b_cuts(train: list[dict]) -> dict[tuple[int,int],float]:
    groups=defaultdict(list)
    for r in train:
        key=(qbin(r["abs_q"]),abin(r["a"]))
        if None not in key:groups[key].append(float(r["b"]))
    return {k:float(np.median(v)) for k,v in groups.items()}


def qab_cell(row: dict,cuts: dict) -> tuple[int,int,int]|None:
    qa=(qbin(row["abs_q"]),abin(row["a"]))
    if None in qa or qa not in cuts:return None
    return qa+(int(float(row["b"])>cuts[qa]),)


def conditional_auc(rows: list[dict],score: str,outcome: str,cuts: dict) -> tuple[float,int]:
    groups=defaultdict(list)
    for r in rows:
        cell=qab_cell(r,cuts)
        if cell is not None:groups[cell].append(r)
    concordant=0.0;pairs=0
    for part in groups.values():
        pos=np.asarray([float(r[score]) for r in part if int(r[outcome])==1])
        neg=np.asarray([float(r[score]) for r in part if int(r[outcome])==0])
        if not len(pos) or not len(neg):continue
        diff=pos[:,None]-neg[None,:]
        concordant+=float(np.sum(diff>0)+0.5*np.sum(diff==0));pairs+=diff.size
    return (concordant/pairs if pairs else float("nan")),pairs


def cluster_auc(rows: list[dict],score: str,outcome: str,nboot: int=1200,
                conditional_cuts: dict|None=None,seed: int=20261118) -> dict:
    seeds=sorted({int(r["seed"]) for r in rows});by={s:[r for r in rows if int(r["seed"])==s] for s in seeds}
    def evaluate(rr):
        if conditional_cuts is not None:return conditional_auc(rr,score,outcome,conditional_cuts)[0]
        return auc_score(np.asarray([int(r[outcome]) for r in rr]),np.asarray([float(r[score]) for r in rr]))
    observed=evaluate(rows);rng=np.random.default_rng(seed);draws=[]
    for _ in range(nboot):
        rr=[]
        for s in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[int(s)])
        value=evaluate(rr)
        if np.isfinite(value):draws.append(value)
    return dict(auc=float(observed),ci95=np.percentile(draws,[2.5,97.5]).tolist(),n=len(rows),n_seeds=len(seeds))


def paired_auc_difference(rows: list[dict],outcome: str,cuts: dict|None=None,nboot: int=1500) -> dict:
    seeds=sorted({int(r["seed"]) for r in rows});by={s:[r for r in rows if int(r["seed"])==s] for s in seeds}
    def evaluate(rr,metric):
        if cuts is None:return auc_score(np.asarray([int(r[outcome]) for r in rr]),np.asarray([float(r[metric]) for r in rr]))
        return conditional_auc(rr,metric,outcome,cuts)[0]
    obs=evaluate(rows,"Jcap")-evaluate(rows,"Ja");rng=np.random.default_rng(20261119);draws=[]
    for _ in range(nboot):
        rr=[]
        for s in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[int(s)])
        d=evaluate(rr,"Jcap")-evaluate(rr,"Ja")
        if np.isfinite(d):draws.append(d)
    return dict(auc_difference_Jcap_minus_Ja=float(obs),ci95=np.percentile(draws,[2.5,97.5]).tolist())


def paired_metric_auc_difference(rows: list[dict],outcome: str,first: str,second: str,nboot: int=1500) -> dict:
    seeds=sorted({int(r["seed"]) for r in rows});by={s:[r for r in rows if int(r["seed"])==s] for s in seeds}
    def evaluate(rr,metric):return auc_score(np.asarray([int(r[outcome]) for r in rr]),np.asarray([float(r[metric]) for r in rr]))
    obs=evaluate(rows,first)-evaluate(rows,second);rng=np.random.default_rng(20261122);draws=[]
    for _ in range(nboot):
        rr=[]
        for s in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[int(s)])
        draws.append(evaluate(rr,first)-evaluate(rr,second))
    return dict(first=first,second=second,auc_difference=float(obs),ci95=np.percentile(draws,[2.5,97.5]).tolist())


def enrich_prospective(raw_dir: Path,rows: list[dict],kind: str) -> list[dict]:
    requests=defaultdict(set)
    for r in rows:requests[int(r["seed"])].add(int(float(r.get("native_step",r.get("step")))))
    frames=load_frames(raw_dir,requests);out=[]
    for source in rows:
        seed=int(source["seed"]);step=int(float(source.get("native_step",source.get("step"))))
        q=float(source["q"]);source_sign=1 if q>0 else -1
        metrics,_=snapshot_metrics(*frames[(seed,step)],source_sign)
        row=dict(seed=seed,step=step,t=float(source["t"]),source_sign=source_sign,kind=kind,**metrics)
        if kind=="initiation":row["outcome"]=int(source["first_crossing"])
        else:
            row["outcome"]=int(source["future4"]);row["future_net"]=int(float(source["future_net"]))
        out.append(row)
    return out


def matched_rows(raw_dir: Path,pairs: list[dict],event_steps: dict[tuple[int,int],int]) -> list[dict]:
    selected=[r for r in pairs if (int(r["event_seed"]),int(r["event_id"])) in event_steps and float(r["lead"]) in (1.0,2.0)]
    requests=defaultdict(set);meta=[]
    for r in selected:
        es=int(r["event_seed"]);cs=int(r["control_seed"]);lead=float(r["lead"])
        event_step=event_steps[(es,int(r["event_id"]))]-int(round(lead/DT));control_step=int(round(float(r["control_time"])/DT))
        requests[es].add(event_step);requests[cs].add(control_step);meta.append((r,event_step,control_step))
    frames=load_frames(raw_dir,requests);out=[]
    for r,ek,ck in meta:
        q=float(r["q"]);source_sign=1 if q>0 else -1
        ev,_=snapshot_metrics(*frames[(int(r["event_seed"]),ek)],source_sign)
        ct,_=snapshot_metrics(*frames[(int(r["control_seed"]),ck)],source_sign)
        row=dict(seed=int(r["event_seed"]),event_id=int(r["event_id"]),lead=float(r["lead"]),q=q,control_seed=int(r["control_seed"]))
        for metric in ("a","b","Ja","Jcap","n_capture"):
            row[f"event_{metric}"]=ev[metric];row[f"control_{metric}"]=ct[metric];row[f"delta_{metric}"]=ev[metric]-ct[metric]
        out.append(row)
    return out


def event_rows(raw_dir: Path,outcomes: list[dict]) -> list[dict]:
    requests=defaultdict(set);events=[]
    for r in outcomes:
        event=dict(seed=int(r["seed"]),event_id=int(r["event_id"]),first_native=int(r["first_native"]),
                   source_sign=int(r["source_sign"]),success4=int(r["success4"]),abs_q_pre=float(r["abs_q_pre"]))
        events.append(event)
        for lead in EVENT_LEADS:requests[event["seed"]].add(event["first_native"]-round(lead/DT))
    frames=load_frames(raw_dir,requests);out=[]
    for ev in events:
        for lead in EVENT_LEADS:
            step=ev["first_native"]-round(lead/DT);metrics,_=snapshot_metrics(*frames[(ev["seed"],step)],ev["source_sign"])
            out.append({**ev,"lead":lead,**metrics})
    return out


def standardized_difference(rows: list[dict],metric: str,nboot: int=1500) -> dict:
    weights=np.bincount([int(np.searchsorted(Q_STANDARD_EDGES,r["abs_q_pre"],side="right")-1) for r in rows],minlength=3).astype(float)
    weights/=weights.sum();groups=defaultdict(list)
    for r in rows:
        qb=int(np.searchsorted(Q_STANDARD_EDGES,r["abs_q_pre"],side="right")-1)
        groups[(r["success4"],qb,r["seed"])].append(float(r[metric]))
    cell={k:float(np.mean(v)) for k,v in groups.items()};seeds=sorted({r["seed"] for r in rows})
    def mean(status,mult=None):
        ans=0.0
        for qb,w in enumerate(weights):
            vals=[];ww=[]
            for (s,q,seed),value in cell.items():
                m=1 if mult is None else mult.get(seed,0)
                if s==status and q==qb and m:vals.append(value);ww.append(m)
            if not vals:return float("nan")
            ans+=w*float(np.average(vals,weights=ww))
        return ans
    c=mean(1);f=mean(0);rng=np.random.default_rng(20261120);draws=[]
    for _ in range(nboot):
        unique,count=np.unique(rng.choice(seeds,len(seeds),replace=True),return_counts=True);mult=dict(zip(unique,count))
        cc=mean(1,mult);ff=mean(0,mult)
        if np.isfinite(cc) and np.isfinite(ff):draws.append((cc,ff,cc-ff))
    arr=np.asarray(draws);ci=np.percentile(arr,[2.5,97.5],axis=0)
    return dict(continued=c,failed=f,difference=c-f,continued_ci95=ci[:,0].tolist(),failed_ci95=ci[:,1].tolist(),difference_ci95=ci[:,2].tolist())


def pair_difference(rows: list[dict],metric: str,nboot: int=2000) -> dict:
    by=defaultdict(list)
    for r in rows:by[int(r["seed"])].append(float(r[f"delta_{metric}"]))
    vals=np.asarray([np.mean(v) for v in by.values()]);rng=np.random.default_rng(20261121)
    draws=rng.choice(vals,(nboot,len(vals)),replace=True).mean(axis=1)
    return dict(mean=float(vals.mean()),ci95=np.percentile(draws,[2.5,97.5]).tolist(),positive_seed_fraction=float(np.mean(vals>0)),n_seeds=len(vals))


def persistence_audit(raw_dir: Path,risk_source: list[dict]) -> tuple[list[dict],list[dict]]:
    """Endpoint audit at each ballistic arrival time plus first-fish recall."""
    by_seed=defaultdict(list)
    for r in risk_source:by_seed[int(r["seed"])].append(r)
    fish_rows=[];origin_rows=[]
    for count,seed in enumerate(sorted(by_seed),1):
        print(f"persistence audit {count}/{len(by_seed)} seed={seed}",flush=True)
        reader=RawReader(raw_path(raw_dir,seed));endpoint_cache={}
        try:
            for source in by_seed[seed]:
                step=int(source["native_step"]);q=float(source["q"]);source_sign=1 if q>0 else -1
                state=reader.frame("state",step);cap=capture_geometry(state,source_sign)
                indices=np.flatnonzero(cap["mask"]);captured_at_origin=set(indices.tolist())
                for fish in indices:
                    eta=int(round(float(cap["tau"][fish])/DT));end_step=step+eta
                    if end_step not in endpoint_cache:endpoint_cache[end_step]=reader.frame("state",end_step)
                    end=endpoint_cache[end_step];dtheta=abs(float(wrap(end[2,fish]-state[2,fish])))
                    error=float(np.hypot(end[0,fish],end[1,fish]-cap["y_hit"][fish]))
                    direction=-source_sign
                    fish_rows.append(dict(seed=seed,step=step,fish=int(fish),tau=float(cap["tau"][fish]),
                                          abs_heading_change=dtheta,heading_retention=math.cos(dtheta),
                                          endpoint_error=error,noise_rms_angle=IN_NOISE*math.sqrt(float(cap["tau"][fish])),
                                          reached_target_by_eta=int(direction*end[0,fish]>=-1.1*DT and abs(end[1,fish])<=HALF_WIDTH),
                                          aperture_margin=float(HALF_WIDTH-abs(cap["y_hit"][fish]))))
                recall=float("nan");crossing_fish=-1
                if int(source["first_crossing"]):
                    event_step=int(round(float(source["first_event_time"])/DT))
                    pre=reader.frame("state",event_step-1);at=reader.frame("state",event_step)
                    crossed=np.flatnonzero((source_sign*pre[0]>0)&(source_sign*at[0]<=0)&(np.abs(at[1])<=HALF_WIDTH+1e-10))
                    if len(crossed):
                        crossing_fish=int(crossed[0]);recall=float(crossing_fish in captured_at_origin)
                origin_rows.append(dict(seed=seed,step=step,n_capture=len(indices),first_crossing=int(source["first_crossing"]),
                                        time_to_event=float(source["first_event_time"])-float(source["t"]) if int(source["first_crossing"]) else float("nan"),
                                        crossing_fish=crossing_fish,eventual_first_fish_captured=recall))
        finally:reader.close()
    return fish_rows,origin_rows


def summarize_persistence(fish: list[dict],origins: list[dict]) -> tuple[list[dict],list[dict]]:
    rows=[]
    for lo,hi in zip(PERSISTENCE_EDGES[:-1],PERSISTENCE_EDGES[1:]):
        part=[r for r in fish if lo<=r["tau"]<hi]
        if not part:continue
        rows.append(dict(tau_low=lo,tau_high=hi,n=len(part),n_seeds=len({r["seed"] for r in part}),
                         median_abs_heading_change=float(np.median([r["abs_heading_change"] for r in part])),
                         median_noise_rms_angle=float(np.median([r["noise_rms_angle"] for r in part])),
                         mean_heading_retention=float(np.mean([r["heading_retention"] for r in part])),
                         median_endpoint_error=float(np.median([r["endpoint_error"] for r in part])),
                         target_fraction_by_eta=float(np.mean([r["reached_target_by_eta"] for r in part]))))
    recall=[]
    for lo,hi in zip((0,.25,.5,1,2),(0.25,.5,1,2,5.0001)):
        part=[r for r in origins if r["first_crossing"] and lo<r["time_to_event"]<=hi and np.isfinite(r["eventual_first_fish_captured"])]
        if part:recall.append(dict(lead_low=lo,lead_high=hi,n=len(part),n_seeds=len({r["seed"] for r in part}),
                                   first_fish_capture_recall=float(np.mean([r["eventual_first_fish_captured"] for r in part]))))
    return rows,recall


def finite_persistence_intercept(z: float,y: float,theta_local: float,omega: float,noise: float=IN_NOISE) -> dict:
    """First-order microscopic correction for locally constant angular drift.

    For dtheta=omega*dt+noise*dW and ballistic arrival tau=-z/cos(theta),
    delta y_hit=(1/cos(theta))*integral_0^tau delta_theta(t)dt.
    No coefficient is fitted.
    """
    u=math.cos(theta_local)
    if z>=0 or u<=0:raise ValueError("Requires source-side z<0 and inward heading")
    tau=-z/u;ballistic=y-z*math.tan(theta_local)
    mean_shift=omega*tau*tau/(2*u);sigma=noise*math.sqrt(tau**3/(3*u*u))
    mean=ballistic+mean_shift
    probability=float(ndtr((HALF_WIDTH-mean)/sigma)-ndtr((-HALF_WIDTH-mean)/sigma)) if sigma>0 else float(abs(mean)<=HALF_WIDTH)
    return dict(tau=tau,ballistic_y_hit=ballistic,mean_y_hit=mean,sigma_y_hit=sigma,capture_probability=probability)


def make_figure(init_stats: list[dict],transfer_stats: list[dict],event_summary: list[dict],persistence: list[dict],output: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10.2,7.6));ax=axes[0,0]
    # Geometry schematic in local coordinates.
    ax.add_patch(plt.Rectangle((-ELL,-HALF_WIDTH),ELL,ELL,fill=False,color=".25"));ax.plot([0,0],[-HALF_WIDTH,HALF_WIDTH],color="tab:green",lw=5,label="opening")
    starts=[(-2.35,-.75,.30,"tab:blue"),(-1.8,.55,-.18,"tab:blue"),(-2.0,-.2,.72,"tab:red")]
    for z,y,ang,color in starts:
        u=np.cos(ang);yh=y+(-z)*np.sin(ang)/u;ax.arrow(z,y,-z,yh-y,width=.015,head_width=.12,length_includes_head=True,color=color,alpha=.85)
    ax.set(xlim=(-ELL-.15,.25),ylim=(-HALF_WIDTH-.25,HALF_WIDTH+.25),xlabel="local source coordinate z",ylabel="y",title="parameter-free ray capture");ax.set_aspect("equal")

    ax=axes[0,1];metrics=["abs_q","a","b","Ja","Jcap"];x=np.arange(len(metrics));w=.36
    init=[next(r for r in init_stats if r["metric"]==m) for m in metrics];transfer=[next(r for r in transfer_stats if r["metric"]==m) for m in metrics]
    for offset,data,label in ((-w/2,init,"future first entry"),(w/2,transfer,">=4 net transfer")):
        y=np.asarray([r["auc"] for r in data]);ci=np.asarray([r["ci95"] for r in data]).T
        ax.bar(x+offset,y,w,label=label,yerr=np.vstack([y-ci[0],ci[1]-y]),capsize=2,error_kw={"lw":.8})
    ax.axhline(.5,color=".4",ls=":");ax.set_xticks(x,["-|q|","a","b","Ja","Jcap"]);ax.set(ylim=(.35,.9),ylabel="held-out AUC",title="prospective independent seeds");ax.legend(frameon=False,fontsize=8)

    ax=axes[1,0]
    for metric,color,ls in (("Jcap","tab:blue","-"),("Ja","tab:orange","--")):
        for status,label in ((1,"continued"),(0,"failed")):
            part=sorted([r for r in event_summary if r["metric"]==metric and r["status"]==status],key=lambda r:r["lead"])
            xx=np.asarray([-r["lead"] for r in part]);yy=np.asarray([r["mean"] for r in part]);ci=np.asarray([r["ci95"] for r in part])
            ax.plot(xx,yy,color=color,ls=ls,marker="o" if status else "s",label=f"{metric} {label}");ax.fill_between(xx,ci[:,0],ci[:,1],color=color,alpha=.08)
    ax.axvline(0,color=".4",ls=":");ax.set(xlabel="time from first crossing",ylabel="current per fish",title="q-standardized crossing attempts");ax.legend(frameon=False,fontsize=7,ncol=2)

    ax=axes[1,1];mid=np.asarray([(r["tau_low"]+r["tau_high"])/2 for r in persistence]);angle=np.asarray([r["median_abs_heading_change"] for r in persistence]);error=np.asarray([r["median_endpoint_error"] for r in persistence])
    ax.plot(mid,angle,"o-",color="tab:purple",label="median |heading change|");ax.set(xlabel="ballistic time to opening",ylabel="angle (rad)",title="finite-persistence audit")
    twin=ax.twinx();twin.plot(mid,error,"s--",color="tab:green",label="median endpoint error");twin.set_ylabel("distance")
    lines=ax.lines+twin.lines;ax.legend(lines,[l.get_label() for l in lines],frameon=False,fontsize=7,loc="center right")
    fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    st=np.zeros((3,N));st[0]=2*ELL;st[1]=0;st[2]=np.pi
    st[0,0]=ELL/2;st[1,0]=0;st[2,0]=np.pi
    st[0,1]=ELL/2;st[1,1]=HALF_WIDTH*.9;st[2,1]=np.pi/2+.2
    cap=capture_geometry(st,1)
    assert cap["mask"][0] and not cap["mask"][1] and np.isclose(cap["tau"][0],ELL/2)
    assert cap["j_capture"]<=cap["j_inward"] and cap["j_inward"]>0
    corr=finite_persistence_intercept(-1,0,0,0)
    assert corr["tau"]==1 and corr["mean_y_hit"]==0 and np.isclose(corr["sigma_y_hit"],IN_NOISE/math.sqrt(3))
    print("PASS: aperture ray condition and finite-persistence limit")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--risk-file",type=Path);p.add_argument("--transfer-file",type=Path);p.add_argument("--event-outcomes",type=Path);p.add_argument("--pairs",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:self_test();return
    required=(args.raw_dir,args.risk_file,args.transfer_file,args.event_outcomes,args.pairs,args.output)
    if not all(required):p.error("all input and output paths are required")
    args.output.mkdir(parents=True,exist_ok=True)
    risk_source=read_csv(args.risk_file);transfer_source=read_csv(args.transfer_file);outcomes=read_csv(args.event_outcomes)
    initiation=enrich_prospective(args.raw_dir,risk_source,"initiation")
    transfer=enrich_prospective(args.raw_dir,transfer_source,"transfer")
    event=event_rows(args.raw_dir,outcomes);event_steps={(int(r["seed"]),int(r["event_id"])):int(r["first_native"]) for r in outcomes}
    matched=matched_rows(args.raw_dir,read_csv(args.pairs),event_steps)
    fish,origin=persistence_audit(args.raw_dir,risk_source);persistence,recall=summarize_persistence(fish,origin)

    # Prior coordinates define cells using training seeds only; test seeds remain untouched.
    init_cuts=fit_b_cuts([r for r in initiation if r["seed"]<30]);transfer_cuts=fit_b_cuts([r for r in transfer if r["seed"]<30])
    init_test=[r for r in initiation if r["seed"]>=30];transfer_test=[r for r in transfer if r["seed"]>=30]
    for rows in (initiation,transfer):
        for r in rows:r["minus_abs_q"]=-r["abs_q"]
    stats=[]
    for label,rows,cuts in (("initiation",init_test,init_cuts),("transfer",transfer_test,transfer_cuts)):
        for metric,field in (("abs_q","minus_abs_q"),("a","a"),("b","b"),("Ja","Ja"),("Jcap","Jcap")):
            result=cluster_auc(rows,field,"outcome",seed=20261118+len(stats));stats.append(dict(cohort=label,metric=metric,conditional=0,**result))
        for metric in ("Ja","Jcap"):
            result=cluster_auc(rows,metric,"outcome",conditional_cuts=cuts,seed=20261218+len(stats));auc,pairs=conditional_auc(rows,metric,"outcome",cuts)
            stats.append(dict(cohort=label,metric=metric,conditional=1,conditional_pairs=pairs,**result))
    comparisons={
        "initiation_overall":paired_auc_difference(init_test,"outcome"),
        "initiation_conditional_qab":paired_auc_difference(init_test,"outcome",init_cuts),
        "transfer_overall":paired_auc_difference(transfer_test,"outcome"),
        "transfer_conditional_qab":paired_auc_difference(transfer_test,"outcome",transfer_cuts),
        "initiation_Jcap_minus_a":paired_metric_auc_difference(init_test,"outcome","Jcap","a"),
        "initiation_Jcap_minus_b":paired_metric_auc_difference(init_test,"outcome","Jcap","b"),
        "transfer_Jcap_minus_a":paired_metric_auc_difference(transfer_test,"outcome","Jcap","a"),
        "transfer_Jcap_minus_b":paired_metric_auc_difference(transfer_test,"outcome","Jcap","b"),
    }
    event_summary=[]
    for lead in EVENT_LEADS:
        part=[r for r in event if r["lead"]==lead]
        for metric in ("a","b","Ja","Jcap","n_capture"):
            d=standardized_difference(part,metric)
            event_summary.extend([dict(lead=lead,metric=metric,status=1,mean=d["continued"],ci95=d["continued_ci95"]),dict(lead=lead,metric=metric,status=0,mean=d["failed"],ci95=d["failed_ci95"]),dict(lead=lead,metric=metric,status=-1,mean=d["difference"],ci95=d["difference_ci95"])])
    pair_summary=[]
    for lead in (2.0,1.0):
        part=[r for r in matched if r["lead"]==lead]
        for metric in ("a","b","Ja","Jcap","n_capture"):pair_summary.append(dict(lead=lead,metric=metric,**pair_difference(part,metric)))

    write_csv(args.output/"prospective_initiation.csv",initiation);write_csv(args.output/"prospective_transfer.csv",transfer)
    write_csv(args.output/"matched_entry_controls.csv",matched);write_csv(args.output/"event_continuation_paths.csv",event)
    write_csv(args.output/"predictive_auc.csv",stats);write_csv(args.output/"matched_and_event_summary.csv",pair_summary+event_summary)
    write_csv(args.output/"persistence_fish.csv",fish);write_csv(args.output/"persistence_summary.csv",persistence);write_csv(args.output/"first_fish_recall.csv",recall)
    init_plot=[r for r in stats if r["cohort"]=="initiation" and not r["conditional"]];transfer_plot=[r for r in stats if r["cohort"]=="transfer" and not r["conditional"]]
    make_figure(init_plot,transfer_plot,event_summary,persistence,args.output/"aperture_capture_summary.png")
    summary=dict(geometry=dict(ell=ELL,half_width=HALF_WIDTH,condition="inside approach, u=s*cos(theta)>0, abs(y-z*sin(theta)/u)<=ell/2",current="Jcap=N^-1 sum captured max(s*cos(theta),0)"),
                 cohorts=dict(initiation=len(initiation),initiation_positive=sum(r["outcome"] for r in initiation),transfer=len(transfer),transfer_positive=sum(r["outcome"] for r in transfer),events=len(outcomes),matched_pairs=len(matched)),
                 heldout_auc_comparisons=comparisons,persistence=persistence,first_fish_recall=recall,
                 finite_persistence_formula=dict(microscopic_law="dtheta=clip(omega_wall+omega_social)dt+In*dW; In=0.5",mean_shift="omega*tau^2/(2u)",transverse_sigma="In*sqrt(tau^3/(3u^2))",capture_probability="Phi((h-mu)/sigma)-Phi((-h-mu)/sigma)",fitted_coefficients=0))
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    production=Path(__file__).resolve().parents[3]/"BEM_bichamber_parallel_instance_no_hydro.py"
    config=dict(raw_dir=str(args.raw_dir),risk_file=str(args.risk_file),transfer_file=str(args.transfer_file),
                event_outcomes=str(args.event_outcomes),matched_pairs=str(args.pairs),dt=DT,N=N,ell=ELL,
                train_seeds="0--29 (q-a-b storage cuts only)",heldout_seeds="30--39",
                production_source=str(production),production_source_sha256=hashlib.sha256(production.read_bytes()).hexdigest(),
                simulations_run=0,fitted_physical_coefficients=0)
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
