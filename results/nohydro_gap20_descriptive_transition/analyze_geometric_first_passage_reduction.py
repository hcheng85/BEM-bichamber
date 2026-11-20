#!/usr/bin/env python3
"""Audit a source-school COM/relative-bearing first-passage reduction.

Derives and verifies the exact crossing-free kinematics, evaluates the
coherent-aggregate velocity approximation, projects the unchanged production
wall/social torques onto the mean-heading angle, and audits COM-entry against
fish-level occupancy of the established approach rectangle.  No model is fit.
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

from analyze_approach_current_balance import DT, ELL, load_gap20_production, raw_path
from analyze_collective_polarization_evolution import D_THETA, IN_NOISE, replay_extended, source_mask
from analyze_initiation_propagation import RawReader


LEADS=(3.,2.,1.,.5,.2,.1)
Q_EDGES=np.asarray([0.,.5,.8,1.000001])
D_EDGES=np.asarray([0.,ELL,2*ELL,np.inf])
C_EDGES=np.asarray([-1.,0.,.5,math.cos(math.pi/6),1.000001])


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def wrap(x):return np.arctan2(np.sin(x),np.cos(x))


def bindex(x: float,edges: np.ndarray) -> int:
    return max(0,min(int(np.searchsorted(edges,x,side="right")-1),len(edges)-2))


def geometry(state: np.ndarray,mask: np.ndarray) -> dict:
    r=np.mean(state[:2,mask],axis=1);d=float(np.linalg.norm(r));g=-r/d
    e=np.column_stack((np.cos(state[2]),np.sin(state[2])));h=np.mean(e[mask],axis=0);p=float(np.linalg.norm(h));m=h/p
    theta=float(math.atan2(m[1],m[0]));chi=float(math.atan2(g[1],g[0]));psi=float(wrap(theta-chi));gp=np.asarray([-g[1],g[0]])
    return dict(r=r,d=d,g=g,gp=gp,e=e,h=h,p=p,m=m,theta=theta,chi=chi,psi=psi,cospsi=float(math.cos(psi)))


def angular_terms(state: np.ndarray,mask: np.ndarray,geo: dict,wall_raw: np.ndarray,social: np.ndarray,wall_clip: np.ndarray,total_clip: np.ndarray) -> dict:
    n=int(mask.sum());weights=mask*np.cos(state[2]-geo["theta"])/(n*geo["p"])
    wall=float(np.sum(weights*wall_raw));soc=float(np.sum(weights*social));clip=float(np.sum(weights*(wall_clip+total_clip)))
    # Ito correction for arg(h), with independent sigma_i=(In/n)e_i_perp.
    ep=np.column_stack((-np.sin(state[2]),np.cos(state[2])))
    mp=np.asarray([-geo["m"][1],geo["m"][0]])
    ito=-IN_NOISE**2/(n*n*geo["p"]**2)*float(np.sum(mask*(ep@mp)*(ep@geo["m"])))
    variance=IN_NOISE**2/(n*n*geo["p"]**2)*float(np.sum(mask*(ep@mp)**2))
    return dict(omega_wall=wall,omega_social=soc,omega_clip=clip,omega_ito=ito,omega_deterministic=wall+soc+clip+ito,omega_noise_variance=variance)


def com_wall_approximation(mod,para: dict,geo: dict) -> float:
    probe=np.asarray([geo["r"][0],geo["r"][1],geo["theta"]]);angle,distance=mod.GetInterc(probe);phi=float(wrap(geo["theta"]-angle))
    value=para["Iw"]*np.sign(phi)/distance
    if abs(distance)<=abs(mod.safedis):value*=1+1/distance**(mod.high_order-1)
    return float(np.clip(value,-np.pi/DT,np.pi/DT))


def event_rows(cohort_path: Path) -> list[dict]:
    rows=[]
    for r in read_csv(cohort_path):
        rows.append(dict(seed=int(r["seed"]),event_id=int(r["event_id"]),first_native=int(r["first_native"]),source_sign=int(r["source_sign"]),
                         direction=-int(r["source_sign"]),success4=int(r["success4"]),abs_q_pre=float(r["abs_q_pre"])))
    if len(rows)!=337:raise RuntimeError("Expected established 337-event cohort")
    return rows


def extract_kinematics(mod,para,raw_dir: Path,cohort: list[dict]) -> tuple[list[dict],dict]:
    rows=[];by=defaultdict(list)
    for e in cohort:by[e["seed"]].append(e)
    max_position_residual=0.;max_psi_identity=0.
    for count,seed in enumerate(sorted(by),1):
        print(f"geometric/torque audit {count}/{len(by)} seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed))
        try:
            for ev in by[seed]:
                for lead in LEADS:
                    k=ev["first_native"]-round(lead/DT);state=reader.frame("state",k);nxt=reader.frame("state",k+1);velocity=reader.frame("rdot",k+1)
                    mask=source_mask(state,ev["direction"])
                    if np.any(mask!=source_mask(nxt,ev["direction"])):raise RuntimeError("Unexpected pre-crossing membership change")
                    z=geometry(state,mask);z1=geometry(nxt,mask);v=np.mean(velocity[:,mask],axis=1)
                    max_position_residual=max(max_position_residual,float(np.max(np.abs(np.mean(nxt[:2,mask]-state[:2,mask],axis=1)-DT*v))))
                    wall,social,wclip,tclip=replay_extended(mod,state,para);turn=angular_terms(state,mask,z,wall,social,wclip,tclip)
                    centered=state[:2,mask]-z["r"][:,None];school_rms=float(np.sqrt(np.mean(np.sum(centered**2,axis=0))));wall_com=com_wall_approximation(mod,para,z)
                    alignment_turn={"alignment_wall":-math.sin(z["psi"])*turn["omega_wall"],"alignment_social":-math.sin(z["psi"])*turn["omega_social"],
                                    "alignment_clip":-math.sin(z["psi"])*turn["omega_clip"],"alignment_geometry":-math.sin(z["psi"])*float(z["gp"]@v/z["d"])}
                    d_rate=float((z1["d"]-z["d"])/DT);theta_rate=float(wrap(z1["theta"]-z["theta"])/DT);bearing_rate=float(wrap(z1["chi"]-z["chi"])/DT);psi_rate=float(wrap(z1["psi"]-z["psi"])/DT)
                    max_psi_identity=max(max_psi_identity,abs(psi_rate-theta_rate+bearing_rate))
                    heading_speed=float(v@z["m"]);vnorm=float(np.linalg.norm(v));hres=float(np.linalg.norm(v-z["h"]));align=float(v@z["m"]/vnorm) if vnorm else np.nan
                    future={}
                    if lead>=.2:
                        j=k+round(.1/DT);zf=geometry(reader.frame("state",j),mask);future=dict(theta_rate_0p1=float(wrap(zf["theta"]-z["theta"])/.1),psi_rate_0p1=float(wrap(zf["psi"]-z["psi"])/.1),d_rate_0p1=(zf["d"]-z["d"])/.1)
                    rows.append(dict(**ev,lead=lead,n_source=int(mask.sum()),com_x=z["r"][0],com_y=z["r"][1],mean_heading_angle=z["theta"],d=z["d"],psi=z["psi"],cospsi=z["cospsi"],polarization=z["p"],school_rms_radius=school_rms,
                                     vx=v[0],vy=v[1],vnorm=vnorm,heading_speed=heading_speed,velocity_heading_cos=align,velocity_minus_heading=hres,
                                     d_rate=d_rate,d_rate_continuous=-float(z["g"]@v),d_rate_error=d_rate+float(z["g"]@v),
                                     theta_rate=theta_rate,bearing_rate=bearing_rate,bearing_rate_continuous=-float(z["gp"]@v/z["d"]),
                                     psi_rate=psi_rate,psi_geometry_rate=float(z["gp"]@v/z["d"]),psi_identity_residual=psi_rate-theta_rate+bearing_rate,
                                     realized_collective_noise=theta_rate-turn["omega_deterministic"],wall_com_approximation=wall_com,wall_com_error=turn["omega_wall"]-wall_com,**turn,**alignment_turn,**future))
        finally:reader.close()
    return rows,dict(max_abs_position_update_residual=max_position_residual,max_abs_psi_rate_identity=max_psi_identity)


def leading_edge_audit(raw_dir: Path,cohort: list[dict],attempts_path: Path) -> list[dict]:
    matched={(int(r["seed"]),int(float(r["event_id"]))) for r in read_csv(attempts_path) if r["event_id"] not in ("","nan")}
    rows=[];by=defaultdict(list)
    for e in cohort:by[e["seed"]].append(e)
    offsets=np.arange(-round(2/DT),0,round(.1/DT),dtype=int)
    for seed,events in by.items():
        reader=RawReader(raw_path(raw_dir,seed))
        try:
            for ev in events:
                n_inside=[];com_inside=[]
                for off in offsets:
                    state=reader.frame("state",ev["first_native"]+int(off));mask=source_mask(state,ev["direction"]);x,y=state[0],state[1];s=ev["source_sign"]
                    inside=mask&(s*x<=ELL)&(np.abs(y)<=ELL/2);n_inside.append(int(inside.sum()))
                    com=np.mean(state[:2,mask],axis=1);com_inside.append(int(s*com[0]<=ELL and abs(com[1])<=ELL/2))
                lead_onsets=[-offsets[i]*DT for i in range(len(offsets)) if n_inside[i]>0 and (i==0 or n_inside[i-1]==0)]
                com_onsets=[-offsets[i]*DT for i in range(len(offsets)) if com_inside[i] and (i==0 or not com_inside[i-1])]
                rows.append(dict(**ev,leading_present_last2=int(any(n_inside)),com_present_last2=int(any(com_inside)),max_leading_count=max(n_inside),
                                 latest_leading_onset_lead=min(lead_onsets) if lead_onsets else np.nan,latest_com_onset_lead=min(com_onsets) if com_onsets else np.nan,
                                 matched_by_prospective_com_attempt=int((seed,ev["event_id"]) in matched)))
        finally:reader.close()
    return rows


def cell_contrast(rows: list[dict],field: str) -> float:
    cells=defaultdict(list)
    for r in rows:
        key=(bindex(r["abs_q_pre"],Q_EDGES),bindex(r["d"],D_EDGES),bindex(r["cospsi"],C_EDGES));cells[key].append(r)
    diff=[];weight=[]
    for part in cells.values():
        a=[r[field] for r in part if r["success4"]];b=[r[field] for r in part if not r["success4"]]
        if a and b:diff.append(float(np.mean(a)-np.mean(b)));weight.append(len(part))
    return float(np.average(diff,weights=weight)) if diff else np.nan


def seed_boot(rows: list[dict],metric,seed: int,nboot=1500) -> tuple[float,list[float]]:
    by=defaultdict(list)
    for r in rows:by[r["seed"]].append(r)
    seeds=list(by);point=float(metric(rows));rng=np.random.default_rng(seed);draw=[]
    for _ in range(nboot):
        sample=[]
        for replicate,s in enumerate(rng.choice(seeds,len(seeds),replace=True)):
            for row in by[int(s)]:z=dict(row);z["seed"]=replicate;sample.append(z)
        value=metric(sample)
        if np.isfinite(value):draw.append(value)
    return point,np.percentile(draw,[2.5,97.5]).tolist()


def summarize(kin: list[dict],leading: list[dict]) -> tuple[list[dict],list[dict],dict]:
    velocity=[];turning=[]
    for lead in LEADS:
        part=[r for r in kin if r["lead"]==lead]
        for cohort,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            selected=[r for r in part if pred(r["seed"])]
            for field in ("velocity_heading_cos","velocity_minus_heading","d_rate","psi_rate"):
                metric=lambda x,f=field:float(np.mean([r[f] for r in x]));point,ci=seed_boot(selected,metric,20261140+round(lead*10))
                velocity.append(dict(lead=lead,cohort=cohort,metric=field,mean=point,ci95_low=ci[0],ci95_high=ci[1],n=len(selected),n_seeds=len({r["seed"] for r in selected})))
        for cohort,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            selected=[r for r in part if pred(r["seed"])]
            for field in ("alignment_wall","alignment_social","alignment_clip","alignment_geometry","omega_deterministic","psi_geometry_rate","d_rate"):
                point,ci=seed_boot(selected,lambda x,f=field:cell_contrast(x,f),20261150+round(lead*10))
                turning.append(dict(lead=lead,cohort=cohort,term=field,continued_minus_failed_geometry_controlled=point,ci95_low=ci[0],ci95_high=ci[1],n=len(selected),n_seeds=len({r["seed"] for r in selected})))
    held=[r for r in leading if r["seed"]>=30];allr=leading
    def rate(rows,field):
        point,ci=seed_boot(rows,lambda x:float(np.mean([r[field] for r in x])),20261170);return dict(mean=point,ci95=ci,n=len(rows),n_seeds=len({r["seed"] for r in rows}))
    closure={"0.001":{"d_rate_rmse":float(np.sqrt(np.mean([r["d_rate_error"]**2 for r in kin]))),"bearing_rate_rmse":float(np.sqrt(np.mean([(r["bearing_rate"]-r["bearing_rate_continuous"])**2 for r in kin])))}}
    for cohort,pred_seed in (("all",lambda s:True),("train",lambda s:s<30),("heldout",lambda s:s>=30)):
        valid=[r for r in kin if isinstance(r.get("psi_rate_0p1"),(int,float,np.integer,np.floating)) and np.isfinite(r["psi_rate_0p1"]) and pred_seed(r["seed"])]
        pp=np.asarray([r["omega_deterministic"]+r["psi_geometry_rate"] for r in valid]);pa=np.asarray([r["psi_rate_0p1"] for r in valid]);dp=np.asarray([r["d_rate_continuous"] for r in valid]);da=np.asarray([r["d_rate_0p1"] for r in valid])
        closure[f"0.1_{cohort}"]={"psi_rate_rmse":float(np.sqrt(np.mean((pa-pp)**2))),"psi_rate_correlation":float(np.corrcoef(pa,pp)[0,1]),"psi_actual_sd":float(np.std(pa)),
                                    "d_rate_rmse":float(np.sqrt(np.mean((da-dp)**2))),"d_rate_correlation":float(np.corrcoef(da,dp)[0,1]),"d_actual_sd":float(np.std(da)),"n":len(valid)}
        selected=[r for r in kin if pred_seed(r["seed"])];actual=np.asarray([r["omega_wall"] for r in selected]);compact=np.asarray([r["wall_com_approximation"] for r in selected])
        closure[f"compact_wall_{cohort}"]={"rmse":float(np.sqrt(np.mean((actual-compact)**2))),"correlation":float(np.corrcoef(actual,compact)[0,1]),"actual_sd":float(np.std(actual)),
                                             "median_school_rms_radius":float(np.median([r["school_rms_radius"] for r in selected])),"median_com_distance":float(np.median([r["d"] for r in selected])),"n":len(selected)}
        dactual=np.asarray([r["d_rate_continuous"] for r in selected]);unit_speed=-np.asarray([r["cospsi"] for r in selected])
        closure[f"unit_speed_translation_{cohort}"]={"rmse":float(np.sqrt(np.mean((dactual-unit_speed)**2))),"bias_predicted_minus_actual":float(np.mean(unit_speed-dactual)),"actual_sd":float(np.std(dactual)),"n":len(selected)}
    audit=dict(all=dict(leading_presence=rate(allr,"leading_present_last2"),com_presence=rate(allr,"com_present_last2"),prospective_com_match=rate(allr,"matched_by_prospective_com_attempt")),
               heldout=dict(leading_presence=rate(held,"leading_present_last2"),com_presence=rate(held,"com_present_last2"),prospective_com_match=rate(held,"matched_by_prospective_com_attempt")),closure=closure)
    return velocity,turning,audit


def figures(kin,turning,leading,out: Path) -> None:
    fig,axes=plt.subplots(1,3,figsize=(12,3.7))
    for success,color,label in ((0,"#B66A50","failed"),(1,"#327C91","continued")):
        for ax,field,title in zip(axes,("d","cospsi","velocity_minus_heading"),("COM distance d","relative alignment cos psi","translation residual |V - h|")):
            x=[];y=[]
            for lead in LEADS:
                p=[r[field] for r in kin if r["lead"]==lead and r["success4"]==success];x.append(-lead);y.append(np.mean(p))
            ax.plot(x,y,"o-",color=color,label=label);ax.set(xlabel="time from first crossing",title=title)
    axes[0].legend(frameon=False);fig.tight_layout();fig.savefig(out/"geometric_kinematics.png",dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));terms=("alignment_social","alignment_wall","alignment_geometry")
    x=np.arange(len(LEADS));width=.24
    for j,term in enumerate(terms):
        rows=[r for r in turning if r["term"]==term and r["cohort"]=="all"];y=np.asarray([r["continued_minus_failed_geometry_controlled"] for r in rows]);lo=np.asarray([r["ci95_low"] for r in rows]);hi=np.asarray([r["ci95_high"] for r in rows]);ax.bar(x+(j-1)*width,y,width,yerr=[y-lo,hi-y],capsize=2,label=term)
    ax.axhline(0,color="0.4",ls=":");ax.set_xticks(x,[str(v) for v in LEADS]);ax.set(xlabel="lead before crossing",ylabel="continued - failed angular rate",title="Residual turning after controlling q, d, and psi");ax.legend(frameon=False);fig.tight_layout();fig.savefig(out/"geometry_controlled_turning.png",dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4));names=("leading_present_last2","com_present_last2","matched_by_prospective_com_attempt");vals=[np.mean([r[n] for r in leading]) for n in names];ax.bar(("any source fish\nin rectangle","COM in rectangle","prospective COM\nattempt matched"),vals);ax.set(ylim=(0,1),ylabel="fraction of 337 crossing episodes",title="COM-entry sensitivity audit");fig.tight_layout();fig.savefig(out/"com_entry_sensitivity.png",dpi=220);plt.close(fig)


def self_test() -> None:
    state=np.asarray([[2.,2.2,1.8],[1.,1.1,.9],[2.4,2.5,2.3]]);mask=np.ones(3,bool);z=geometry(state,mask);v=np.asarray([.3,-.2]);eps=1e-7
    moved=state.copy();moved[:2]+=eps*v[:,None];z1=geometry(moved,mask)
    assert abs((z1["d"]-z["d"])/eps+z["g"]@v)<1e-7
    assert abs(wrap(z1["chi"]-z["chi"])/eps+z["gp"]@v/z["d"])<1e-7
    theta_rate=.4;rot=state.copy();rot[2]+=eps*theta_rate;rot[:2]+=eps*v[:,None];zr=geometry(rot,mask)
    assert abs(wrap(zr["psi"]-z["psi"])/eps-(theta_rate+z["gp"]@v/z["d"]))<1e-7
    print("PASS: exact COM-distance and relative-bearing kinematics")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--cohort",type=Path);p.add_argument("--attempts",type=Path);p.add_argument("--source",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.raw_dir,a.cohort,a.attempts,a.source,a.output)):p.error("all input/output paths required")
    a.output.mkdir(parents=True,exist_ok=True);cohort=event_rows(a.cohort)
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(a.source)
    para=dict(Ip=4.,In=IN_NOISE,R=mod.R,light=1.,dt=DT,Iw=.94,If=0.,num=100,delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
    kin,exact_audit=extract_kinematics(mod,para,a.raw_dir,cohort);leading=leading_edge_audit(a.raw_dir,cohort,a.attempts)
    velocity,turning,audit=summarize(kin,leading);audit.update(exact_audit)
    write_csv(a.output/"kinematic_samples.csv",kin);write_csv(a.output/"velocity_summary.csv",velocity);write_csv(a.output/"geometry_controlled_turning.csv",turning);write_csv(a.output/"leading_edge_audit.csv",leading)
    figures(kin,turning,leading,a.output)
    config=dict(cohort_size=len(cohort),leads=list(LEADS),approach_length=ELL,approach_width=ELL,train_seeds="0-29",heldout_seeds="30-39",new_simulations=0,fitted_turning_law=0,
                production_source=str(a.source),production_source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest())
    (a.output/"summary.json").write_text(json.dumps(dict(configuration=config,audits=audit,velocity_summary=velocity,geometry_controlled_turning=turning),indent=2)+"\n")


if __name__=="__main__":main()
