#!/usr/bin/env python3
"""Exact source-school polarization evolution from production torques.

The stochastic identity is evaluated at prospective snapshots from the fixed
337-event cohort.  Every chamber-crossing step is evaluated separately so
membership jumps are not folded into a smooth drift.  No drift is fitted.
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

from analyze_approach_current_balance import DT, ELL, N, RawReader, load_gap20_production, raw_path, replay_torques


IN_NOISE=.5
D_THETA=.5*IN_NOISE**2
START,STOP,SAMPLE_DT=-2.,1.,.1
PERIODS=(("early_pre",-2.,-1.),("late_pre",-1.,0.),("immediate_post",0.,1.))
Q_EDGES=np.asarray([0.,.5,.8,1.000001])
DRIFT_TERMS=("wall","social","clip","noise_relaxation","geometry_motion","realized_noise","delta_a")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");writer.writeheader();writer.writerows(rows)


def source_mask(state: np.ndarray,direction: int) -> np.ndarray:
    return state[0]<0 if direction==1 else state[0]>=0


def geometry(state: np.ndarray,mask: np.ndarray) -> dict:
    if not np.any(mask):raise RuntimeError("Original source chamber is empty")
    com=np.mean(state[:2,mask],axis=1);radius=float(np.linalg.norm(com))
    if radius<1e-12:raise RuntimeError("Source COM at gap: reference undefined")
    g=-com/radius;h=np.column_stack((np.cos(state[2]),np.sin(state[2])));p=np.mean(h[mask],axis=0)
    return dict(mask=mask,n=int(mask.sum()),com=com,radius=radius,g=g,p=p,a=float(p@g))


def raw_wall_and_clips(mod,state: np.ndarray,para: dict,returned_wall: np.ndarray,final_clip: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Expose both saturation operations in the unchanged production rule."""
    raw=np.zeros(N);limit=np.pi/DT
    for i in range(N):
        angle,distance=mod.GetInterc(state[:,i]);phi=math.atan2(math.sin(state[2,i]-angle),math.cos(state[2,i]-angle))
        value=para["Iw"]*np.sign(phi)/distance
        if abs(distance)<=abs(mod.safedis):value*=1+1/distance**(mod.high_order-1)
        raw[i]=value
    internal=np.clip(raw,-limit,limit)-raw
    reconstructed=raw+internal
    if np.max(np.abs(reconstructed-returned_wall))>2e-10:raise RuntimeError("Wall saturation replay mismatch")
    return raw,internal,final_clip


def exact_step(state: np.ndarray,nxt: np.ndarray,direction: int,wall_raw: np.ndarray,
               social: np.ndarray,wall_clip: np.ndarray,total_clip: np.ndarray) -> dict:
    """Exact finite-step Ito/source/membership split in production update order."""
    old_mask=source_mask(state,direction);new_mask=source_mask(nxt,direction)
    before=geometry(state,old_mask);after=geometry(nxt,new_mask)
    effective=wall_raw+social+wall_clip+total_clip;delta=effective*DT
    alpha=math.atan2(before["g"][1],before["g"][0]);phi=state[2]-alpha
    slope=-np.sin(phi+.5*delta)*np.sinc(delta/(2*np.pi));weight=old_mask/before["n"]
    terms={
        "wall":float(np.sum(weight*slope*wall_raw*DT)),
        "social":float(np.sum(weight*slope*social*DT)),
        "wall_clip":float(np.sum(weight*slope*wall_clip*DT)),
        "total_clip":float(np.sum(weight*slope*total_clip*DT)),
    }
    terms["clip"]=terms["wall_clip"]+terms["total_clip"]
    deterministic=np.cos(phi+delta);factor=math.exp(-D_THETA*DT)
    projected_next=np.cos(nxt[2]-alpha)
    terms["noise_relaxation"]=float(np.sum(weight*(factor-1)*deterministic))
    terms["realized_noise"]=float(np.sum(weight*(projected_next-factor*deterministic)))
    post_h=np.column_stack((np.cos(nxt[2]),np.sin(nxt[2])))
    p_old=np.mean(post_h[old_mask],axis=0)
    old_next=geometry(np.vstack((nxt[:2],nxt[2])),old_mask)
    terms["geometry_motion"]=float(p_old@(old_next["g"]-before["g"]))
    p_new=np.mean(post_h[new_mask],axis=0)
    terms["membership_heading"]=float((p_new-p_old)@old_next["g"])
    terms["membership_geometry"]=float(p_new@(after["g"]-old_next["g"]))
    terms["membership"]=terms["membership_heading"]+terms["membership_geometry"]
    terms["delta_a"]=after["a"]-before["a"]
    terms["reconstruction"]=sum(terms[k] for k in ("wall","social","clip","noise_relaxation","realized_noise","geometry_motion","membership"))
    terms["residual"]=terms["delta_a"]-terms["reconstruction"]
    terms.update(a_before=before["a"],a_after=after["a"],n_before=before["n"],n_after=after["n"],
                 n_membership_changes=int(np.count_nonzero(old_mask!=new_mask)),net_source_loss=before["n"]-after["n"])
    return terms


def load_cohort(path: Path) -> list[dict]:
    rows=[]
    for r in read_csv(path):
        q=float(r["abs_q_pre"]);qb=int(np.clip(np.searchsorted(Q_EDGES,q,side="right")-1,0,len(Q_EDGES)-2))
        rows.append(dict(seed=int(r["seed"]),event_id=int(r["event_id"]),first_native=int(r["first_native"]),
                         source_sign=int(r["source_sign"]),direction=-int(r["source_sign"]),success4=int(r["success4"]),
                         additional_fish=int(float(r["additional_fish"])),abs_q_pre=q,q_bin=qb))
    if len(rows)!=337 or len({(r["seed"],r["event_id"]) for r in rows})!=337:raise RuntimeError("Expected unique fixed 337-event cohort")
    return rows


def replay_extended(mod,state: np.ndarray,para: dict) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    wall,social,total_clip,_=replay_torques(mod,state,para)
    wall_raw,wall_clip,total_clip=raw_wall_and_clips(mod,state,para,wall,total_clip)
    return wall_raw,social,wall_clip,total_clip


def extract(mod,para: dict,raw_dir: Path,cohort: list[dict]) -> tuple[list[dict],list[dict],dict]:
    samples=[];jumps=[];max_residual=0.;wall_clip_count=0;total_clip_count=0
    by_seed=defaultdict(list)
    for event in cohort:by_seed[event["seed"]].append(event)
    sample_offsets=np.rint(np.arange(START,STOP+SAMPLE_DT/2,SAMPLE_DT)/DT).astype(int)
    for count,seed in enumerate(sorted(by_seed),1):
        print(f"collective polarization replay {count}/{len(by_seed)} seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed))
        try:
            for event in by_seed[seed]:
                anchor=event["first_native"];lo=anchor+int(round(START/DT));hi=anchor+int(round(STOP/DT))+1
                block=reader.block("state",lo,hi)
                old=source_mask(block[:,:,0],event["direction"])
                changes=[]
                for k in range(block.shape[2]-1):
                    new=source_mask(block[:,:,k+1],event["direction"])
                    if np.any(old!=new):changes.append(k)
                    old=new
                indices=sorted(set((sample_offsets+anchor-lo).tolist()+changes));cache={}
                for k in indices:
                    st=block[:,:,k];nxt=block[:,:,k+1];wall_raw,social,wall_clip,total_clip=replay_extended(mod,st,para)
                    z=exact_step(st,nxt,event["direction"],wall_raw,social,wall_clip,total_clip);cache[k]=z
                    max_residual=max(max_residual,abs(z["residual"]));wall_clip_count+=int(np.count_nonzero(wall_clip));total_clip_count+=int(np.count_nonzero(total_clip))
                for offset in sample_offsets:
                    k=int(offset+anchor-lo);z=cache[k];tau=offset*DT
                    row={**event,"tau":tau,"crossing_step":int(z["n_membership_changes"]>0)}
                    for term in DRIFT_TERMS:row[term]=z[term]/DT
                    row.update(a=z["a_before"],n_source=z["n_before"],wall_clip=z["wall_clip"]/DT,total_clip=z["total_clip"]/DT,
                               deterministic_angular=(z["wall"]+z["social"]+z["clip"]+z["noise_relaxation"])/DT,
                               smooth_total=(z["wall"]+z["social"]+z["clip"]+z["noise_relaxation"]+z["geometry_motion"])/DT,
                               closure_residual=z["residual"])
                    samples.append(row)
                for k in changes:
                    z=cache[k];tau=(lo+k+1-anchor)*DT
                    jumps.append({**event,"tau":tau,**{name:z[name] for name in ("membership_heading","membership_geometry","membership","net_source_loss","n_membership_changes","residual")}})
        finally:reader.close()
    audit=dict(max_abs_step_closure=max_residual,wall_internal_clip_count=wall_clip_count,total_torque_clip_count=total_clip_count,
               evaluated_snapshot_steps=len(samples),evaluated_crossing_steps=len(jumps),n_events=len(cohort),n_seeds=len(by_seed))
    return samples,jumps,audit


def target_weights(cohort: list[dict]) -> np.ndarray:
    count=np.bincount([r["q_bin"] for r in cohort],minlength=len(Q_EDGES)-1);return count/count.sum()


def cell_means(rows: list[dict],value: str) -> dict[tuple[int,int,int],float]:
    by=defaultdict(list)
    for r in rows:by[(int(r["success4"]),int(r["q_bin"]),int(r["seed"]))].append(float(r[value]))
    return {k:float(np.mean(v)) for k,v in by.items()}


def standardized(cell: dict,outcome: int,weights: np.ndarray,multiplicity: dict[int,int]|None=None) -> float:
    answer=0.
    for qb,w in enumerate(weights):
        vals=[];ww=[]
        for (status,q,seed),value in cell.items():
            mult=1 if multiplicity is None else multiplicity.get(seed,0)
            if status==outcome and q==qb and mult:vals.append(value);ww.append(mult)
        if w and not vals:return float("nan")
        if w:answer+=w*float(np.average(vals,weights=ww))
    return answer


def contrast(rows: list[dict],value: str,weights: np.ndarray,nboot=1500,seed=20261118) -> dict:
    cell=cell_means(rows,value);c=standardized(cell,1,weights);f=standardized(cell,0,weights);seeds=np.asarray(sorted({int(r["seed"]) for r in rows}));rng=np.random.default_rng(seed);draw=[]
    for _ in range(nboot):
        unique,count=np.unique(rng.choice(seeds,len(seeds),replace=True),return_counts=True);mult=dict(zip(unique.tolist(),count.tolist()))
        cc=standardized(cell,1,weights,mult);ff=standardized(cell,0,weights,mult)
        if np.isfinite(cc) and np.isfinite(ff):draw.append((cc,ff,cc-ff))
    arr=np.asarray(draw)
    if len(arr)<.85*nboot:raise RuntimeError(f"Insufficient bootstrap support {len(arr)}/{nboot}")
    ci=np.percentile(arr,[2.5,97.5],axis=0)
    return dict(continued=c,failed=f,difference=c-f,continued_ci95=json.dumps(ci[:,0].tolist()),failed_ci95=json.dumps(ci[:,1].tolist()),difference_ci95=json.dumps(ci[:,2].tolist()),bootstrap_draws=len(arr))


def summarize_time(samples: list[dict],weights: np.ndarray) -> list[dict]:
    out=[]
    for j,tau in enumerate(sorted({r["tau"] for r in samples})):
        part=[r for r in samples if r["tau"]==tau and not r["crossing_step"]]
        for i,term in enumerate(("a","wall","social","clip","noise_relaxation","geometry_motion","realized_noise","smooth_total")):
            out.append(dict(tau=tau,term=term,**contrast(part,term,weights,seed=20262000+50*j+i)))
    return out


def summarize_periods(samples: list[dict],weights: np.ndarray) -> list[dict]:
    out=[]
    for j,(name,lo,hi) in enumerate(PERIODS):
        part=[r for r in samples if lo<=r["tau"]<hi and not r["crossing_step"]]
        for i,term in enumerate(("wall","social","clip","noise_relaxation","geometry_motion","realized_noise","smooth_total")):
            out.append(dict(period=name,tau_start=lo,tau_end=hi,term=term,**contrast(part,term,weights,seed=20263000+50*j+i)))
    return out


def summarize_jumps(jumps: list[dict],cohort: list[dict],weights: np.ndarray) -> list[dict]:
    out=[]
    for j,(name,lo,hi) in enumerate(PERIODS):
        by={(e["seed"],e["event_id"]):{**e,"membership":0.,"membership_heading":0.,"membership_geometry":0.,"net_source_loss":0.} for e in cohort}
        for r in jumps:
            if lo<r["tau"]<=hi:
                target=by[(r["seed"],r["event_id"])]
                for term in ("membership","membership_heading","membership_geometry","net_source_loss"):target[term]+=float(r[term])/(hi-lo)
        part=list(by.values())
        for i,term in enumerate(("membership","membership_heading","membership_geometry","net_source_loss")):
            out.append(dict(period=name,tau_start=lo,tau_end=hi,term=term,**contrast(part,term,weights,seed=20264000+50*j+i)))
    return out


def make_figure(time_rows: list[dict],period_rows: list[dict],jump_rows: list[dict],output: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10,7.2));colors={"continued":"tab:blue","failed":"tab:orange"}
    ax=axes[0,0]
    for status in colors:
        part=sorted((r for r in time_rows if r["term"]=="a"),key=lambda r:r["tau"]);y=np.asarray([r[status] for r in part]);ci=np.asarray([json.loads(r[f"{status}_ci95"]) for r in part]);x=np.asarray([r["tau"] for r in part]);ax.plot(x,y,color=colors[status],label=status);ax.fill_between(x,ci[:,0],ci[:,1],color=colors[status],alpha=.15)
    ax.set(title="source-school gap polarization",ylabel="a");ax.legend(frameon=False)
    ax=axes[0,1]
    for term,color in (("social","tab:green"),("wall","tab:red"),("geometry_motion","tab:purple"),("noise_relaxation","tab:gray")):
        part=sorted((r for r in time_rows if r["term"]==term),key=lambda r:r["tau"]);ax.plot([r["tau"] for r in part],[r["difference"] for r in part],label=term,color=color)
    ax.axhline(0,color="k",ls=":");ax.set(title="continued - failed drift",ylabel="polarization / time");ax.legend(frameon=False,fontsize=8)
    ax=axes[1,0];terms=("wall","social","clip","noise_relaxation","geometry_motion","realized_noise");x=np.arange(len(terms));width=.24
    for j,(period,_,_) in enumerate(PERIODS):
        vals=[next(r["difference"] for r in period_rows if r["period"]==period and r["term"]==term) for term in terms];ax.bar(x+(j-1)*width,vals,width,label=period)
    ax.axhline(0,color="k",ls=":");ax.set_xticks(x,terms,rotation=35,ha="right");ax.set(title="period-mean continued - failed",ylabel="polarization / time");ax.legend(frameon=False,fontsize=8)
    ax=axes[1,1];terms=("membership_heading","membership_geometry");x=np.arange(len(PERIODS));width=.34
    for j,term in enumerate(terms):
        vals=[next(r["difference"] for r in jump_rows if r["period"]==p[0] and r["term"]==term) for p in PERIODS];ax.bar(x+(j-.5)*width,vals,width,label=term)
    ax.axhline(0,color="k",ls=":");ax.set_xticks(x,[p[0] for p in PERIODS],rotation=20);ax.set(title="crossing membership jumps",ylabel="polarization / time");ax.legend(frameon=False,fontsize=8)
    for ax in axes[0]:ax.axvline(0,color="k",ls=":",lw=.8);ax.set_xlabel("time from first crossing")
    fig.suptitle("Exact production-source polarization terms; outcome contrasts are retrospective",fontsize=11);fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    rng=np.random.default_rng(20261120);x=np.r_[rng.uniform(-3,-.1,60),rng.uniform(.1,3,40)];state=np.vstack((x,rng.uniform(-1,1,N),rng.uniform(-np.pi,np.pi,N)));nxt=state.copy()
    wall=rng.normal(0,.3,N);social=rng.normal(0,.4,N);wall_clip=np.zeros(N);total_clip=np.zeros(N);nxt[2]+=DT*(wall+social)+rng.normal(0,IN_NOISE*np.sqrt(DT),N);nxt[:2]+=DT*np.vstack((np.cos(state[2]),np.sin(state[2])));nxt[0,59]=.01
    z=exact_step(state,nxt,1,wall,social,wall_clip,total_clip)
    assert abs(z["residual"])<2e-12 and z["n_membership_changes"]==1
    stationary=exact_step(state,state,1,np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N));assert abs(stationary["residual"])<1e-14
    print("PASS: exact smooth and crossing-jump polarization identities")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--cohort",type=Path);p.add_argument("--source",type=Path);p.add_argument("--transport-periods",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.cohort,args.source,args.output)):p.error("--raw-dir, --cohort, --source, and --output required")
    args.output.mkdir(parents=True,exist_ok=True);cohort=load_cohort(args.cohort);weights=target_weights(cohort)
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(args.source)
    para=dict(Ip=4.,In=IN_NOISE,R=mod.R,light=1.,dt=DT,Iw=.94,If=0.,num=N,delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
    samples,jumps,audit=extract(mod,para,args.raw_dir,cohort);time_rows=summarize_time(samples,weights);period_rows=summarize_periods(samples,weights);jump_rows=summarize_jumps(jumps,cohort,weights)
    write_csv(args.output/"time_resolved_terms.csv",time_rows);write_csv(args.output/"period_term_contrasts.csv",period_rows);write_csv(args.output/"crossing_jump_contrasts.csv",jump_rows)
    make_figure(time_rows,period_rows,jump_rows,args.output/"collective_polarization_mechanism.png")
    transport=[]
    if args.transport_periods:
        wanted={"plug_flux","storage_contribution","lateral_contribution","neck_flux"}
        transport=[r for r in read_csv(args.transport_periods) if r["period"] in {p[0] for p in PERIODS} and r["term"] in wanted]
        write_csv(args.output/"linked_transport_budget.csv",transport)
    audit.update(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),q_stratum_weights=weights.tolist(),
                 equation="wall + social + wall_internal_clip + final_total_clip + Ito_noise_relaxation + realized_noise + fixed_membership_COM_geometry + crossing_membership_jump")
    (args.output/"summary.json").write_text(json.dumps(dict(audit=audit,period_terms=period_rows,crossing_jumps=jump_rows,linked_transport=transport),indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),cohort=str(args.cohort),source=str(args.source),transport_periods=str(args.transport_periods) if args.transport_periods else None,
                window=[START,STOP],snapshot_spacing=SAMPLE_DT,n_events=len(cohort),n_seeds=len({r['seed'] for r in cohort}),new_simulations=0,fitted_terms=0,
                interpretation="source identities exact; success/failure contrasts retrospective and q-standardized")
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
