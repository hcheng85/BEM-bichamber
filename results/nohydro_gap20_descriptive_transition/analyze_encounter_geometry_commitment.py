#!/usr/bin/env python3
"""Encounter geometry versus genuine source-school reorientation.

Compares fixed-axis and source-COM-to-opening alignment on the established
337-event cohort, decomposes their pre-crossing changes exactly, and defines
prospective COM-entry attempts using only current geometry and past population
history.  Existing trajectories are read only.
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
from scipy.stats import rankdata

from analyze_approach_current_balance import ELL
from analyze_initiation_propagation import RawReader


DT=.001
GRID_DT=.5
LEADS=(5.,3.,2.,1.,.5,.2,.1)
Q_EDGES=np.asarray([0.,.25,.5,.75,.9,1.000001])
RHO_EDGES=np.asarray([0.,ELL,2*ELL,np.inf])


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def raw_path(raw_dir: Path,seed: int) -> Path:
    return raw_dir/f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"


def geometry(state: np.ndarray,source_sign: int|None=None) -> dict:
    x,y,theta=np.asarray(state);q=float(2*np.sum(x>=0)/len(x)-1)
    invalid=lambda:dict(q=q,valid=0,source_sign=0,n_source=0,com_x=np.nan,com_y=np.nan,h_x=np.nan,h_y=np.nan,polarization=np.nan,gap_x=np.nan,gap_y=np.nan,rho=np.nan,axial_distance=np.nan,lateral_distance=np.nan,alignment_fixed=np.nan,alignment_com=np.nan,in_approach=0)
    if source_sign is None:
        if abs(q)<1e-12:return invalid()
        source_sign=1 if q>0 else -1
    source=(x>=0) if source_sign>0 else (x<0)
    if source.sum()<3:return invalid()
    com=np.asarray([x[source].mean(),y[source].mean()]);h=np.asarray([np.cos(theta[source]).mean(),np.sin(theta[source]).mean()]);rho=float(np.linalg.norm(com))
    if rho<1e-12:return invalid()
    gap=-com/rho;u=float(source_sign*com[0]);ay=float(abs(com[1]));af=float(-source_sign*h[0]);ag=float(h@gap)
    return dict(q=q,valid=1,source_sign=source_sign,n_source=int(source.sum()),com_x=float(com[0]),com_y=float(com[1]),h_x=float(h[0]),h_y=float(h[1]),
                polarization=float(np.linalg.norm(h)),gap_x=float(gap[0]),gap_y=float(gap[1]),rho=rho,axial_distance=u,lateral_distance=ay,
                alignment_fixed=af,alignment_com=ag,in_approach=int(u<=ELL and ay<=ELL/2))


def bindex(x: float,edges: np.ndarray) -> int:
    return max(0,min(int(np.searchsorted(edges,x,side="right")-1),len(edges)-2))


def auc(labels,scores) -> float:
    y=np.asarray(labels,bool);s=np.asarray(scores,float);ok=np.isfinite(s);y=y[ok];s=s[ok];n1=int(y.sum());n0=len(y)-n1
    if not n1 or not n0:return np.nan
    return float((rankdata(s)[y].sum()-n1*(n1+1)/2)/(n1*n0))


def stratified_auc(rows: list[dict],feature: str,outcome: str,proximity=True) -> float:
    correct=0.;pairs=0;groups=defaultdict(list)
    for r in rows:
        key=(bindex(abs(float(r["q"])),Q_EDGES),bindex(float(r["rho"]),RHO_EDGES) if proximity else 0);groups[key].append(r)
    for part in groups.values():
        n1=sum(int(r[outcome]) for r in part);n0=len(part)-n1
        if n1 and n0:
            value=auc([int(r[outcome]) for r in part],[float(r[feature]) for r in part]);correct+=value*n1*n0;pairs+=n1*n0
    return float(correct/pairs) if pairs else np.nan


def seed_boot(rows: list[dict],metric,nboot=1500,seed=20261120) -> tuple[float,list[float]]:
    by=defaultdict(list)
    for r in rows:by[int(r["seed"])].append(r)
    seeds=list(by);point=float(metric(rows));rng=np.random.default_rng(seed);draw=[]
    for _ in range(nboot):
        sample=[]
        for replicate,s in enumerate(rng.choice(seeds,len(seeds),replace=True)):
            for r in by[int(s)]:copy=dict(r);copy["seed"]=replicate;sample.append(copy)
        value=metric(sample)
        if np.isfinite(value):draw.append(value)
    return point,np.percentile(draw,[2.5,97.5]).tolist() if draw else [np.nan,np.nan]


def equal_seed_mean(rows: list[dict],field: str) -> float:
    by=defaultdict(list)
    for r in rows:
        if np.isfinite(float(r[field])):by[int(r["seed"])].append(float(r[field]))
    return float(np.mean([np.mean(v) for v in by.values()])) if by else np.nan


def event_geometry(raw_dir: Path,outcomes: list[dict]) -> list[dict]:
    readers={};rows=[]
    try:
        for ev in outcomes:
            seed=int(ev["seed"])
            if seed not in readers:readers[seed]=RawReader(raw_path(raw_dir,seed))
            first=int(ev["first_native"]);sign=int(ev["source_sign"])
            for lead in LEADS:
                f=geometry(readers[seed].frame("state",first-round(lead/DT)),sign)
                if abs(f["q"]-sign*float(ev["abs_q_pre"]))>1e-12:raise RuntimeError(f"event q mismatch seed={seed} event={ev['event_id']}")
                rows.append(dict(seed=seed,event_id=int(ev["event_id"]),lead=lead,success4=int(ev["success4"]),additional_fish=int(ev["additional_fish"]),abs_q_pre=float(ev["abs_q_pre"]),**f))
    finally:
        for r in readers.values():r.close()
    return rows


def alignment_prediction(rows: list[dict]) -> list[dict]:
    out=[]
    for lead in LEADS:
        for cohort,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            part=[r for r in rows if r["lead"]==lead and pred(int(r["seed"]))]
            values={}
            for field in ("alignment_fixed","alignment_com"):
                metric=lambda x,f=field:stratified_auc(x,f,"success4",True);point,ci=seed_boot(part,metric,seed=20262000+round(10*lead));values[field]=(point,ci)
                out.append(dict(lead=lead,cohort=cohort,metric=field,q_proximity_stratified_auc=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(part),n_success=sum(r["success4"] for r in part),n_seeds=len({r["seed"] for r in part})))
            metric=lambda x:stratified_auc(x,"alignment_com","success4",True)-stratified_auc(x,"alignment_fixed","success4",True)
            point,ci=seed_boot(part,metric,seed=20262100+round(10*lead));out.append(dict(lead=lead,cohort=cohort,metric="com_minus_fixed",q_proximity_stratified_auc=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(part),n_success=sum(r["success4"] for r in part),n_seeds=len({r["seed"] for r in part})))
    return out


def q_rho_standardized_difference(rows: list[dict],field: str) -> float:
    groups=defaultdict(list)
    for r in rows:groups[(bindex(abs(float(r["q"])),Q_EDGES),bindex(float(r["rho"]),RHO_EDGES))].append(r)
    diffs=[];weights=[]
    for part in groups.values():
        yes=[r for r in part if int(r["success4"])];no=[r for r in part if not int(r["success4"])]
        if yes and no:diffs.append(equal_seed_mean(yes,field)-equal_seed_mean(no,field));weights.append(len(part))
    return float(np.average(diffs,weights=weights)) if diffs else np.nan


def exact_decomposition(rows: list[dict]) -> tuple[list[dict],list[dict]]:
    lookup={(r["seed"],r["event_id"],r["lead"]):r for r in rows};samples=[];summary=[]
    for target in (1.,.5,.2,.1):
        for r0 in rows:
            if r0["lead"]!=2.:continue
            r1=lookup[(r0["seed"],r0["event_id"],target)];h0=np.asarray([r0["h_x"],r0["h_y"]]);h1=np.asarray([r1["h_x"],r1["h_y"]]);g0=np.asarray([r0["gap_x"],r0["gap_y"]]);g1=np.asarray([r1["gap_x"],r1["gap_y"]])
            gbar=(g0+g1)/2;heading=float((h1-h0)@gbar);geometry_term=float(((h0+h1)/2)@(g1-g0));total=float(r1["alignment_com"]-r0["alignment_com"])
            p0=float(np.linalg.norm(h0));p1=float(np.linalg.norm(h1));m0=h0/p0;m1=h1/p1
            coherence=float((p1-p0)*((m0+m1)/2)@gbar);heading_direction=float(((p0+p1)/2)*(m1-m0)@gbar)
            samples.append(dict(seed=r0["seed"],event_id=r0["event_id"],target_lead=target,success4=r0["success4"],q=r0["q"],rho=r0["rho"],heading_contribution=heading,geometry_contribution=geometry_term,
                                coherence_contribution=coherence,heading_direction_contribution=heading_direction,heading_split_residual=heading-coherence-heading_direction,
                                total_com_alignment_change=total,fixed_alignment_change=r1["alignment_fixed"]-r0["alignment_fixed"],rho_change=r1["rho"]-r0["rho"],axial_change=r1["axial_distance"]-r0["axial_distance"],residual=total-heading-geometry_term))
    for target in (1.,.5,.2,.1):
        part=[r for r in samples if r["target_lead"]==target]
        for field in ("heading_contribution","coherence_contribution","heading_direction_contribution","geometry_contribution","total_com_alignment_change","fixed_alignment_change","rho_change","axial_change"):
            metric=lambda x,f=field:q_rho_standardized_difference(x,f);point,ci=seed_boot(part,metric,seed=20262200+round(target*10))
            summary.append(dict(baseline_lead=2.,target_lead=target,term=field,success_minus_failure=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(part),n_seeds=len({r["seed"] for r in part})))
    return samples,summary


def load_q(cache_dir: Path,seed: int) -> tuple[np.ndarray,np.ndarray]:
    with np.load(cache_dir/f"seed{seed:02d}_stride100_qscan_v4.npz") as data:steps=np.asarray(data["step"],int);q=np.asarray(data["q"])
    return steps,np.rint((q+1)*50).astype(np.int16)


def sample_grid(raw_dir: Path,q_cache: Path,cache: Path) -> dict[int,list[dict]]:
    if cache.exists():
        with np.load(cache) as data:
            arrays={k:np.asarray(data[k]) for k in data.files}
        out={}
        for seed in range(40):
            take=arrays["seed"]==seed;out[seed]=[{k:arrays[k][i].item() for k in arrays if k!="seed"} for i in np.flatnonzero(take)]
        return out
    allrows=[];out={}
    for seed in range(40):
        print(f"source geometry grid {seed+1}/40 seed={seed}",flush=True);rd=RawReader(raw_path(raw_dir,seed));steps=np.arange(round(200/DT),rd.info["state"].shape[-1]-1,round(GRID_DT/DT),dtype=int);rows=[]
        try:
            for index,step in enumerate(steps):
                f=geometry(rd.frame("state",int(step)));rows.append(dict(index=index,step=int(step),t=step*DT,**f))
        finally:rd.close()
        out[seed]=rows
        for r in rows:allrows.append(dict(seed=seed,**r))
    arrays={k:np.asarray([r[k] for r in allrows]) for k in allrows[0]};cache.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(cache,**arrays)
    return out


def prospective_attempts(grid: dict[int,list[dict]],q_cache: Path,episodes: list[dict],outcomes: list[dict]) -> tuple[list[dict],list[dict]]:
    by_episode=defaultdict(list)
    for e in episodes:by_episode[int(e["seed"])].append(e)
    for value in by_episode.values():value.sort(key=lambda e:int(e["first_native"]))
    outcome={(int(e["seed"]),int(e["event_id"])):e for e in outcomes};attempts=[];risk=[]
    for seed,rows in grid.items():
        steps,nright=load_q(q_cache,seed);step_to_index={int(s):i for i,s in enumerate(steps)};eligible_time=0.;next_allowed=-1
        for i in range(1,len(rows)):
            r=rows[i];previous=rows[i-1]
            if not r.get("valid") or r["step"] not in step_to_index:continue
            qi=step_to_index[r["step"]];quiet=qi>=50 and np.all(nright[qi-50:qi+1]==nright[qi]);future_ok=r["step"]+round(7/DT)<steps[-1]
            if quiet and future_ok and i>=next_allowed:eligible_time+=GRID_DT
            onset=quiet and future_ok and i>=next_allowed and int(r["in_approach"])==1 and int(previous.get("in_approach",0))==0
            if not onset:continue
            sign=int(r["source_sign"]);candidate=[e for e in by_episode[seed] if r["step"]<=int(e["first_native"])<=r["step"]+round(2/DT)]
            first=min(candidate,key=lambda e:int(e["first_native"])) if candidate else None;entry=bool(first and int(first["direction"])==-sign);matched=outcome.get((seed,int(first["event_id"]))) if entry else None;commit=int(matched["success4"]) if matched else 0
            failure="committed" if commit else ("entry_without_commitment" if entry else ("opposite_crossing" if first else "no_entry"))
            attempts.append(dict(seed=seed,attempt_id=sum(a["seed"]==seed for a in attempts),step=r["step"],t=r["t"],q=r["q"],rho=r["rho"],axial_distance=r["axial_distance"],lateral_distance=r["lateral_distance"],
                                 alignment_fixed=r["alignment_fixed"],alignment_com=r["alignment_com"],polarization=r["polarization"],com_inward_speed=(previous["rho"]-r["rho"])/GRID_DT,
                                 entry=int(entry),commit=commit,outcome=failure,crossing_delay=((int(first["first_native"])-r["step"])*DT if first else np.nan),event_id=(int(first["event_id"]) if first else "")))
            next_allowed=i+round(2/GRID_DT)
        risk.append(dict(seed=seed,eligible_time=eligible_time,n_attempts=sum(a["seed"]==seed for a in attempts),n_entries=sum(a["seed"]==seed and a["entry"] for a in attempts),n_commits=sum(a["seed"]==seed and a["commit"] for a in attempts)))
    return attempts,risk


def attempt_prediction(rows: list[dict]) -> list[dict]:
    out=[]
    for cohort,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
        part=[r for r in rows if pred(int(r["seed"]))]
        for outcome in ("entry","commit"):
            for field in ("alignment_fixed","alignment_com","com_inward_speed"):
                metric=lambda x,f=field,o=outcome:stratified_auc(x,f,o,True);point,ci=seed_boot(part,metric,seed=20262300)
                out.append(dict(cohort=cohort,outcome=outcome,metric=field,q_proximity_stratified_auc=point,ci95_low=ci[0],ci95_high=ci[1],n_attempts=len(part),n_positive=sum(int(r[outcome]) for r in part),n_seeds=len({r["seed"] for r in part})))
            metric=lambda x,o=outcome:stratified_auc(x,"alignment_com",o,True)-stratified_auc(x,"alignment_fixed",o,True)
            point,ci=seed_boot(part,metric,seed=20262301)
            out.append(dict(cohort=cohort,outcome=outcome,metric="com_minus_fixed",q_proximity_stratified_auc=point,ci95_low=ci[0],ci95_high=ci[1],n_attempts=len(part),n_positive=sum(int(r[outcome]) for r in part),n_seeds=len({r["seed"] for r in part})))
    return out


def renewal_summary(attempts: list[dict],risk: list[dict]) -> dict:
    active=[r for r in risk if r["eligible_time"]>0];rates=np.asarray([r["n_attempts"]/r["eligible_time"] for r in active]);entry=np.asarray([r["n_entries"]/r["n_attempts"] for r in active if r["n_attempts"]]);commit=np.asarray([r["n_commits"]/r["n_attempts"] for r in active if r["n_attempts"]]);product=np.asarray([r["n_commits"]/r["eligible_time"] for r in active])
    intervals=[];interval_pairs=[];commit_pairs=[]
    for seed in range(40):
        sequence=sorted((r for r in attempts if r["seed"]==seed),key=lambda r:r["t"]);t=[r["t"] for r in sequence];dt=np.diff(t);intervals.extend(dt)
        interval_pairs.extend(dict(seed=seed,x=float(dt[i]),y=float(dt[i+1])) for i in range(len(dt)-1))
        commit_pairs.extend(dict(seed=seed,x=float(sequence[i]["commit"]),y=float(sequence[i+1]["commit"])) for i in range(len(sequence)-1))
    def boot(values,seed):
        values=np.asarray(values,float);rng=np.random.default_rng(seed);pick=rng.integers(len(values),size=(4000,len(values)));return dict(mean=float(values.mean()),ci95=np.percentile(values[pick].mean(axis=1),[2.5,97.5]).tolist(),n=len(values))
    def seed_boot_corr(rows,seed):
        by=defaultdict(list)
        for r in rows:by[r["seed"]].append(r)
        def corr(sample):
            x=np.asarray([r["x"] for r in sample]);y=np.asarray([r["y"] for r in sample]);return float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else np.nan
        seeds=list(by);point=corr(rows);rng=np.random.default_rng(seed);draw=[]
        for _ in range(4000):
            sample=[]
            for s in rng.choice(seeds,len(seeds),replace=True):sample.extend(by[int(s)])
            value=corr(sample)
            if np.isfinite(value):draw.append(value)
        return dict(correlation=point,ci95=np.percentile(draw,[2.5,97.5]).tolist(),n_pairs=len(rows))
    return dict(attempt_rate_per_time=boot(rates,20262401),entry_probability_equal_seed=boot(entry,20262402),commit_probability_equal_seed=boot(commit,20262403),committed_rate_per_time=boot(product,20262404),
                n_attempts=len(attempts),n_entries=sum(r["entry"] for r in attempts),n_commits=sum(r["commit"] for r in attempts),n_no_target_entry=sum(not r["entry"] for r in attempts),
                n_opposite_crossings=sum(r["outcome"]=="opposite_crossing" for r in attempts),
                interattempt_interval_mean=float(np.mean(intervals)) if intervals else None,interattempt_interval_cv=float(np.std(intervals)/np.mean(intervals)) if intervals else None,
                lag1_interattempt_interval=seed_boot_corr(interval_pairs,20262405),lag1_commitment=seed_boot_corr(commit_pairs,20262406),
                identity="committed rate = attempt rate x per-attempt commitment only for consistently weighted aggregate estimates; equal-seed summaries need not multiply exactly")


def make_figures(events,prediction,decomposition,attempts,attempt_pred,output: Path) -> None:
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,field,title in zip(axes,("alignment_fixed","alignment_com","rho"),("fixed-axis alignment","COM-to-opening alignment","COM distance to opening")):
        for success,color,label in ((0,"#B66A50","<4 additional"),(1,"#327C91",">=4 additional")):
            x=[];y=[]
            for lead in LEADS:
                part=[r for r in events if r["lead"]==lead and r["success4"]==success];x.append(-lead);y.append(equal_seed_mean(part,field))
            ax.plot(x,y,"o-",color=color,label=label)
        ax.axvline(0,color="k",ls=":",lw=.7);ax.set(xlabel="time from first crossing",ylabel=title,title=title)
    axes[0].legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(output/"encounter_geometry_chronology.png",dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,term,title in zip(axes,("heading_direction_contribution","geometry_contribution","total_com_alignment_change"),("mean-heading rotation","COM-geometry contribution","total COM-alignment change")):
        part=[r for r in decomposition if r["term"]==term];x=[-r["target_lead"] for r in part];y=np.asarray([r["success_minus_failure"] for r in part]);lo=np.asarray([r["ci95_low"] for r in part]);hi=np.asarray([r["ci95_high"] for r in part]);ax.errorbar(x,y,yerr=[y-lo,hi-y],marker="o",capsize=2);ax.axhline(0,color="0.5",ls=":");ax.set(xlabel="target time",ylabel="continued - failed change",title=title)
    fig.tight_layout();fig.savefig(output/"turning_vs_geometry_decomposition.png",dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(8.5,3.8))
    for ax,outcome in zip(axes,("entry","commit")):
        part=[r for r in attempt_pred if r["cohort"]=="heldout" and r["outcome"]==outcome and r["metric"]!="com_minus_fixed"];names=[r["metric"].replace("alignment_","") for r in part];y=np.asarray([r["q_proximity_stratified_auc"] for r in part]);lo=np.asarray([r["ci95_low"] for r in part]);hi=np.asarray([r["ci95_high"] for r in part]);ax.bar(names,y,yerr=[y-lo,hi-y],capsize=3);ax.axhline(.5,color="0.5",ls=":");ax.set(ylim=(0,1),ylabel="held-out stratified AUC",title=f"prospective attempt: {outcome}")
    fig.tight_layout();fig.savefig(output/"prospective_attempt_validation.png",dpi=220);plt.close(fig)


def self_test() -> None:
    x=np.r_[np.full(6,2.),np.full(4,-2.)];y=np.r_[np.linspace(-1,1,6),np.linspace(-1,1,4)];theta=np.zeros(10);st=np.vstack((x,y,theta));a=geometry(st,1);mirror=st.copy();mirror[0]*=-1;mirror[2]=np.pi-mirror[2];b=geometry(mirror,-1)
    for key in ("rho","axial_distance","lateral_distance","alignment_fixed","alignment_com"):assert abs(a[key]-b[key])<1e-12
    h0=np.asarray([.3,.4]);h1=np.asarray([.1,.7]);g0=np.asarray([1.,0.]);g1=np.asarray([0.,1.]);head=(h1-h0)@((g0+g1)/2);geo=((h0+h1)/2)@(g1-g0);assert abs((h1@g1-h0@g0)-head-geo)<1e-12
    print("PASS: reflection-covariant encounter geometry and exact finite-change identity")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--outcomes",type=Path);p.add_argument("--episodes",type=Path);p.add_argument("--q-cache",type=Path);p.add_argument("--geometry-cache",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.outcomes,args.episodes,args.q_cache,args.geometry_cache,args.output)):p.error("all paths required")
    args.output.mkdir(parents=True,exist_ok=True);outcomes=read_csv(args.outcomes);episodes=read_csv(args.episodes);events=event_geometry(args.raw_dir,outcomes);prediction=alignment_prediction(events);samples,decomposition=exact_decomposition(events)
    grid=sample_grid(args.raw_dir,args.q_cache,args.geometry_cache);attempts,risk=prospective_attempts(grid,args.q_cache,episodes,outcomes);attempt_pred=attempt_prediction(attempts);renewal=renewal_summary(attempts,risk)
    write_csv(args.output/"event_geometry.csv",events);write_csv(args.output/"alignment_prediction.csv",prediction);write_csv(args.output/"turning_geometry_samples.csv",samples);write_csv(args.output/"turning_geometry_decomposition.csv",decomposition);write_csv(args.output/"prospective_attempts.csv",attempts);write_csv(args.output/"attempt_risk_by_seed.csv",risk);write_csv(args.output/"attempt_prediction.csv",attempt_pred)
    make_figures(events,prediction,decomposition,attempts,attempt_pred,args.output)
    summary=dict(configuration=dict(event_cohort=337,event_definition=">=4 additional net fish after first majority-to-minority crossing",grid_dt=GRID_DT,attempt_definition="first source-COM entry into conservation-law approach rectangle after >=5 crossing-free time units; nonoverlapping 2-unit outcome windows",q_edges=Q_EDGES.tolist(),rho_edges=RHO_EDGES.tolist(),train_seeds="0-29",heldout_seeds="30-39",new_simulations=0),
                 alignment_prediction=prediction,turning_geometry_decomposition=decomposition,prospective_attempt_prediction=attempt_pred,renewal_summary=renewal,
                 audits=dict(max_abs_decomposition_residual=max(abs(r["residual"]) for r in samples),max_abs_heading_split_residual=max(abs(r["heading_split_residual"]) for r in samples),strict_pre_crossing_events=all(r["lead"]>0 for r in events),attempt_definition_uses_future=0))
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":main()
