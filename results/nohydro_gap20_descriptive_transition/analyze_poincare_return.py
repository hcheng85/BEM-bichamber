#!/usr/bin/env python3
"""Test a physical chamber-orbit Poincare return hypothesis.

The section is fixed a priori: the outward ray from the exact area centroid of
each production chamber.  Its state is only COM radius and mean-heading angle
relative to the direction of section crossing.  Population is bookkeeping.
No nonlinear map, potential, or new microscopic simulation is introduced.
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

from analyze_approach_current_balance import load_gap20_production, raw_path
from analyze_initiation_propagation import RawReader


DT=.5
N=100
TRANSIENT=200.
OUTCOMES=("no_entry","returned","partial","sustained")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");writer.writeheader();writer.writerows(rows)


def numeric(row: dict[str,str]) -> dict:
    out={}
    for key,value in row.items():
        try:
            x=float(value);out[key]=int(x) if x.is_integer() else x
        except (ValueError,TypeError):out[key]=value
    return out


def wrap(x):return np.arctan2(np.sin(x),np.cos(x))


def clip_half(poly: list[tuple[float,float]],side: int) -> list[tuple[float,float]]:
    out=[]
    for a,b in zip(poly,[*poly[1:],poly[0]]):
        inside_a=side*a[0]>=0;inside_b=side*b[0]>=0
        if inside_a:out.append(a)
        if inside_a!=inside_b:
            f=-a[0]/(b[0]-a[0]);out.append((0.,a[1]+f*(b[1]-a[1])))
    return out


def polygon_centroid(poly: list[tuple[float,float]]) -> tuple[float,float,float]:
    x=np.asarray([p[0] for p in poly]);y=np.asarray([p[1] for p in poly]);cross=x*np.roll(y,-1)-np.roll(x,-1)*y;area=cross.sum()/2
    cx=((x+np.roll(x,-1))*cross).sum()/(6*area);cy=((y+np.roll(y,-1))*cross).sum()/(6*area)
    return float(area),float(cx),float(cy)


def production_geometry(source: Path) -> dict:
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(source)
    poly=list(zip(np.asarray(mod.x_coord,float),np.asarray(mod.y_coord,float)));left=polygon_centroid(clip_half(poly,-1));right=polygon_centroid(clip_half(poly,1))
    if not np.allclose([left[0],left[2],right[0],right[2]],[right[0],0,left[0],0],atol=1e-10):raise RuntimeError("production chamber symmetry audit failed")
    return dict(R=float(mod.R),gap=float(mod.gap),left_area=left[0],right_area=right[0],center_abs=right[1],left_center=[left[1],left[2]],right_center=[right[1],right[2]],boundary_x=np.asarray(mod.x_coord),boundary_y=np.asarray(mod.y_coord))


def sample_orbit(raw_dir: Path,center_abs: float,cache: Path|None) -> list[dict]:
    if cache and cache.exists():
        with np.load(cache) as z:
            fields=list(z.files);return [{key:z[key][i].item() for key in fields} for i in range(len(z[fields[0]]))]
    rows=[]
    for seed in range(40):
        print(f"chamber orbit snapshots {seed+1}/40 seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed));steps=np.arange(round(TRANSIENT/.001),reader.info["state"].shape[-1]-1,round(DT/.001),dtype=int)
        try:
            for index,step in enumerate(steps):
                state=reader.frame("state",int(step));x,y,theta=state;q=float(2*np.sum(x>=0)/N-1)
                for side in (-1,1):
                    mask=x<0 if side<0 else x>=0;n=int(mask.sum())
                    if n<3:continue
                    comx=float(x[mask].mean());comy=float(y[mask].mean());hx=float(np.cos(theta[mask]).mean());hy=float(np.sin(theta[mask]).mean());p=float(np.hypot(hx,hy))
                    xi=float(-side*(comx-side*center_abs));eta=comy;hxi=float(-side*hx);heta=hy
                    rows.append(dict(seed=seed,index=index,step=int(step),t=float(step*.001),side=side,n=n,q=q,com_x=comx,com_y=comy,xi=xi,eta=eta,r=float(np.hypot(xi,eta)),phi=float(math.atan2(eta,xi)),hxi=hxi,heta=heta,polarization=p))
        finally:reader.close()
    if cache:
        arrays={key:np.asarray([r[key] for r in rows]) for key in rows[0]};cache.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(cache,**arrays)
    return rows


def grouped(rows: list[dict]) -> dict[tuple[int,int],list[dict]]:
    by=defaultdict(list)
    for r in rows:by[(int(r["seed"]),int(r["side"]))].append(r)
    for seq in by.values():seq.sort(key=lambda r:r["index"])
    return by


def section_events(orbit: list[dict]) -> tuple[list[dict],list[dict]]:
    by=grouped(orbit);events=[];segments=[]
    for (seed,side),seq in by.items():
        local=[]
        for k,(a,b) in enumerate(zip(seq[:-1],seq[1:])):
            if b["index"]!=a["index"]+1 or b["n"]!=a["n"] or a["eta"]==b["eta"] or a["eta"]*b["eta"]>0:continue
            f=-a["eta"]/(b["eta"]-a["eta"]);xi=a["xi"]+f*(b["xi"]-a["xi"])
            if xi>=0:continue
            deta=b["eta"]-a["eta"];chirality=-1 if deta>0 else 1 # sign of angular phase velocity at phi=pi
            hxi=a["hxi"]+f*(b["hxi"]-a["hxi"]);heta=a["heta"]+f*(b["heta"]-a["heta"]);alpha=float(math.atan2(chirality*hxi,-chirality*heta))
            local.append(dict(seed=seed,side=side,event_id=len(local),left_index=k,right_index=k+1,t=float(a["t"]+f*(b["t"]-a["t"])),chirality=chirality,r=float(-xi),alpha=alpha,n=int(a["n"]),q=float(a["q"]+f*(b["q"]-a["q"])),polarization=float(a["polarization"]+f*(b["polarization"]-a["polarization"]))))
        for previous,current in zip(local[:-1],local[1:]):
            lo,hi=int(previous["right_index"]),int(current["left_index"]);path=seq[lo:hi+1];dphi=np.asarray([float(wrap(v["phi"]-u["phi"])) for u,v in zip(path[:-1],path[1:])]);membership=sum(v["n"]!=u["n"] for u,v in zip(path[:-1],path[1:]));same_chirality=int(previous["chirality"]==current["chirality"])
            chirality=int(previous["chirality"]);tangent=[]
            for row in path:
                tangent_angle=row["phi"]+chirality*math.pi/2;tangent.append((row["hxi"]*math.cos(tangent_angle)+row["heta"]*math.sin(tangent_angle))/max(row["polarization"],1e-12))
            segments.append(dict(seed=seed,side=side,event_id=previous["event_id"],next_event_id=current["event_id"],chirality=chirality,next_chirality=current["chirality"],same_chirality=same_chirality,t=previous["t"],next_t=current["t"],return_time=current["t"]-previous["t"],r=previous["r"],alpha=previous["alpha"],next_r=current["r"],next_alpha=current["alpha"],n=previous["n"],next_n=current["n"],delta_n=current["n"]-previous["n"],membership_changes=membership,phase_advance=float(dphi.sum()),directional_phase_fraction=float(np.mean(np.sign(dphi)==chirality)) if len(dphi) else np.nan,mean_tangent_alignment=float(np.mean(tangent)),radial_sd=float(np.std([r["r"] for r in path]))))
        events.extend(local)
    return events,segments


def attempt_pairs(attempt_path: Path) -> tuple[list[dict],list[dict]]:
    attempts=[numeric(r) for r in read_csv(attempt_path)];by=defaultdict(list)
    for r in attempts:by[int(r["seed"])].append(r)
    pairs=[]
    for seed,seq in by.items():
        seq.sort(key=lambda r:r["attempt_id"])
        for a,b in zip(seq[:-1],seq[1:]):
            if a["outcome"]==b["outcome"]=="no_entry" and a["gross_crossings"]==b["gross_crossings"]==0:
                pairs.append(dict(seed=seed,attempt_id=a["attempt_id"],next_attempt_id=b["attempt_id"],direction=a["direction"],next_direction=b["direction"],same_direction=int(a["direction"]==b["direction"]),onset_time=a["onset_upper"],next_onset_time=b["onset_upper"],onset_interval=b["onset_upper"]-a["onset_upper"]))
    return attempts,pairs


def onset_phase_audit(pairs: list[dict],orbit: list[dict]) -> list[dict]:
    by=grouped(orbit);out=[]
    for pair in pairs:
        if not pair["same_direction"]:continue
        side=-int(pair["direction"]);seq=by[(int(pair["seed"]),side)];times=np.asarray([r["t"] for r in seq]);lo=int(np.searchsorted(times,pair["onset_time"]));hi=int(np.searchsorted(times,pair["next_onset_time"]));path=seq[max(0,lo):min(len(seq),hi+1)]
        if len(path)<2:continue
        dphi=np.asarray([float(wrap(v["phi"]-u["phi"])) for u,v in zip(path[:-1],path[1:])]);changes=sum(v["n"]!=u["n"] for u,v in zip(path[:-1],path[1:]))
        out.append(dict(**pair,side=side,phase_advance=float(dphi.sum()),absolute_windings=float(abs(dphi.sum())/(2*math.pi)),path_length=float(np.sum(np.abs(dphi))),directional_coherence=float(abs(dphi.sum())/max(np.sum(np.abs(dphi)),1e-12)),membership_changes=changes))
    return out


def associate_encounters(events: list[dict],attempts: list[dict]) -> list[dict]:
    egroup=defaultdict(list);agroup=defaultdict(list)
    for e in events:egroup[(int(e["seed"]),int(e["side"]))].append(e)
    for a in attempts:agroup[(int(a["seed"]),-int(a["direction"]))].append(a)
    rows=[]
    for key,seq in egroup.items():
        seq.sort(key=lambda r:r["t"]);alist=sorted(agroup[key],key=lambda r:r["onset_upper"])
        for event,nxt in zip(seq[:-1],seq[1:]):
            candidates=[a for a in alist if event["t"]<a["onset_upper"]<=nxt["t"]];a=candidates[0] if candidates else None
            rows.append(dict(seed=event["seed"],side=event["side"],event_id=event["event_id"],t=event["t"],r=event["r"],alpha=event["alpha"],chirality=event["chirality"],n=event["n"],next_section_time=nxt["t"],encounter=int(a is not None),encounter_delay=(a["onset_upper"]-event["t"] if a else np.nan),attempt_id=(a["attempt_id"] if a else ""),outcome=(a["outcome"] if a else "none"),population_transfer=(a["population_transfer"] if a else 0),first_entry=(a["first_entry"] if a else 0),sustained=(a["sustained"] if a else 0)))
    return rows


def post_outcome_cycles(attempts: list[dict]) -> dict:
    result={}
    for cohort,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30)):
        result[cohort]={}
        for outcome in OUTCOMES:
            part=[r for r in attempts if pred(int(r["seed"])) and r["outcome"]==outcome and np.isfinite(float(r["cycle_time"]))]
            result[cohort][outcome]=dict(n=len(part),median_recovery_wait=float(np.median([r["recovery_wait"] for r in part])),median_onset_interval=float(np.median([r["cycle_time"] for r in part])))
    return result


def cluster_boot(rows: list[dict],metric,seed=20261210,nboot=4000) -> dict:
    by=defaultdict(list)
    for r in rows:by[int(r["seed"])].append(r)
    seeds=list(by);point=float(metric(rows));rng=np.random.default_rng(seed);draw=[]
    for _ in range(nboot):
        sample=[]
        for rep,s in enumerate(rng.choice(seeds,len(seeds),replace=True)):
            for r in by[int(s)]:z=dict(r);z["seed"]=rep;sample.append(z)
        value=metric(sample)
        if np.isfinite(value):draw.append(float(value))
    return dict(estimate=point,ci95=np.percentile(draw,[2.5,97.5]).tolist(),n=len(rows),n_seeds=len(seeds))


def fit_affine(rows: list[dict],inputs: tuple[str,...],target: tuple[str,...]):
    x=np.asarray([[1.,*[r[k] for k in inputs]] for r in rows]);y=np.asarray([[r[k] for k in target] for r in rows]);return np.linalg.lstsq(x,y,rcond=None)[0]


def predict(rows: list[dict],coef,inputs):
    return np.asarray([[1.,*[r[k] for k in inputs]] for r in rows])@coef


def map_validation(segments: list[dict],encounters: list[dict]) -> tuple[dict,list[dict]]:
    # Population is recorded but not supplied to the map.  Including changed-
    # membership returns is the direct test requested by the bookkeeping-only
    # role of q; the fixed-membership count is reported separately.
    clean=[r for r in segments if r["same_chirality"] and abs(r["alpha"])<math.pi*.95]
    train=[r for r in clean if r["seed"]<30];test=[r for r in clean if r["seed"]>=30];state=("r","alpha");target=("next_r","next_alpha")
    coef=fit_affine(train,state,target);mean=np.mean([[r[k] for k in target] for r in train],axis=0);scale=np.std([[r[k] for k in target] for r in train],axis=0);pmap=predict(test,coef,state);actual=np.asarray([[r[k] for k in target] for r in test]);pconst=np.tile(mean,(len(test),1))
    by=defaultdict(list)
    for i,r in enumerate(test):by[r["seed"]].append(i)
    score=[]
    for seed,idx in by.items():
        i=np.asarray(idx);m0=float(np.mean(((actual[i]-pconst[i])/scale)**2));m1=float(np.mean(((actual[i]-pmap[i])/scale)**2));score.append(dict(seed=seed,target="section_state",n=len(i),constant_mse=m0,map_mse=m1,improvement=m0-m1))
    # Physical return time: current section state versus one constant period.
    ct=fit_affine(train,state,("return_time",));mean_t=float(np.mean([r["return_time"] for r in train]));pt=predict(test,ct,state).reshape(-1);yt=np.asarray([r["return_time"] for r in test])
    for seed,idx in by.items():
        i=np.asarray(idx);m0=float(np.mean((yt[i]-mean_t)**2));m1=float(np.mean((yt[i]-pt[i])**2));score.append(dict(seed=seed,target="return_time",n=len(i),constant_mse=m0,map_mse=m1,improvement=m0-m1))
    # Previous-return memory diagnostic, only where a preceding clean return exists.
    lookup={(r["seed"],r["side"],r["next_event_id"]):r for r in clean};memory=[]
    for r in clean:
        previous=lookup.get((r["seed"],r["side"],r["event_id"]))
        if previous and previous["chirality"]==r["chirality"]:
            z=dict(r);z["previous_r"]=previous["r"];z["previous_alpha"]=previous["alpha"];memory.append(z)
    mtrain=[r for r in memory if r["seed"]<30];mtest=[r for r in memory if r["seed"]>=30];coef_current=fit_affine(mtrain,state,target);coef_history=fit_affine(mtrain,("r","alpha","previous_r","previous_alpha"),target);ma=np.asarray([[r[k] for k in target] for r in mtest]);pc=predict(mtest,coef_current,state);ph=predict(mtest,coef_history,("r","alpha","previous_r","previous_alpha"));mscale=np.std([[r[k] for k in target] for r in mtrain],axis=0)
    mby=defaultdict(list)
    for i,r in enumerate(mtest):mby[r["seed"]].append(i)
    for seed,idx in mby.items():
        i=np.asarray(idx);score.append(dict(seed=seed,target="history_state",n=len(i),constant_mse=float(np.mean(((ma[i]-pc[i])/mscale)**2)),map_mse=float(np.mean(((ma[i]-ph[i])/mscale)**2)),improvement=float(np.mean(((ma[i]-pc[i])/mscale)**2-((ma[i]-ph[i])/mscale)**2))))
    # Section state as predictor of the subsequent encounter and reward.
    etrain=[r for r in encounters if r["seed"]<30];etest=[r for r in encounters if r["seed"]>=30]
    for field in ("encounter","population_transfer","sustained"):
        use_train=etrain if field=="encounter" else [r for r in etrain if r["encounter"]]
        use_test=etest if field=="encounter" else [r for r in etest if r["encounter"]]
        ec=fit_affine(use_train,state,(field,));base=float(np.mean([r[field] for r in use_train]));yp=np.asarray([r[field] for r in use_test]);p0=np.full(len(use_test),base);p1=predict(use_test,ec,state).reshape(-1);eby=defaultdict(list)
        for i,r in enumerate(use_test):eby[r["seed"]].append(i)
        for seed,idx in eby.items():
            i=np.asarray(idx);score.append(dict(seed=seed,target=field,n=len(i),constant_mse=float(np.mean((yp[i]-p0[i])**2)),map_mse=float(np.mean((yp[i]-p1[i])**2)),improvement=float(np.mean((yp[i]-p0[i])**2-(yp[i]-p1[i])**2))))
    summaries=[]
    for j,name in enumerate(("section_state","return_time","history_state","encounter","population_transfer","sustained")):
        part=[r for r in score if r["target"]==name];v=np.asarray([r["improvement"] for r in part]);rng=np.random.default_rng(20261230+j);pick=rng.integers(len(v),size=(10000,len(v)))
        summaries.append(dict(target=name,n=sum(r["n"] for r in part),n_seeds=len(part),constant_mse=float(np.mean([r["constant_mse"] for r in part])),map_mse=float(np.mean([r["map_mse"] for r in part])),improvement=float(v.mean()),improvement_ci95=np.percentile(v[pick].mean(axis=1),[2.5,97.5]).tolist()))
    eig=np.linalg.eigvals(coef[1:,:].T)
    model=dict(n_clean_train=len(train),n_clean_test=len(test),n_no_sampled_count_change_train=sum(r["membership_changes"]==0 for r in train),n_no_sampled_count_change_test=sum(r["membership_changes"]==0 for r in test),fixed_point=np.linalg.solve(np.eye(2)-coef[1:,:].T,coef[0,:]).tolist(),linear_map=coef.tolist(),eigenvalues=[[float(z.real),float(z.imag)] for z in eig],constant_period=mean_t,return_time_coefficients=ct.reshape(-1).tolist(),heldout=summaries,n_memory_train=len(mtrain),n_memory_test=len(mtest))
    return model,score


def shuffle_control(segments: list[dict],n=500) -> dict:
    clean=[r for r in segments if r["same_chirality"] and abs(r["alpha"])<math.pi*.95];train=[r for r in clean if r["seed"]<30];test=[r for r in clean if r["seed"]>=30];state=("r","alpha");target=("next_r","next_alpha");scale=np.std([[r[k] for k in target] for r in train],axis=0);actual=np.asarray([[r[k] for k in target] for r in test]);mean=np.mean([[r[k] for k in target] for r in train],axis=0);base=float(np.mean(((actual-mean)/scale)**2));coef=fit_affine(train,state,target);observed=base-float(np.mean(((actual-predict(test,coef,state))/scale)**2));rng=np.random.default_rng(20261250);draw=[]
    y=np.asarray([[r[k] for k in target] for r in train])
    for _ in range(n):
        shuffled=y[rng.permutation(len(y))];x=np.asarray([[1.,r["r"],r["alpha"]] for r in train]);c=np.linalg.lstsq(x,shuffled,rcond=None)[0];draw.append(base-float(np.mean(((actual-predict(test,c,state))/scale)**2)))
    return dict(observed_improvement=observed,shuffled_improvement_ci95=np.percentile(draw,[2.5,97.5]).tolist(),fraction_shuffled_at_least_observed=float(np.mean(np.asarray(draw)>=observed)),n_shuffle=n)


def summaries(events,segments,onset_phase,pairs,encounters) -> dict:
    same=[r for r in segments if r["same_chirality"] and r["membership_changes"]==0];train=[r for r in same if r["seed"]<30];held=[r for r in same if r["seed"]>=30]
    def med(rows,field):return cluster_boot(rows,lambda x:float(np.median([r[field] for r in x])),20261300)
    ptrain=[r for r in pairs if r["seed"]<30];pheld=[r for r in pairs if r["seed"]>=30];oph=[r for r in onset_phase if r["seed"]>=30]
    return dict(
        n_section_events=len(events),
        n_same_chirality_no_sampled_count_change_returns=len(same),
        chirality_reversal_fraction=float(np.mean([not r["same_chirality"] for r in segments])),
        no_sampled_count_change_return_period_train=med(train,"return_time"),
        no_sampled_count_change_return_period_heldout=med(held,"return_time"),
        no_entry_onset_pairs=dict(
            all=len(pairs),train=len(ptrain),heldout=len(pheld),
            median_all=float(np.median([r["onset_interval"] for r in pairs])),
            median_train=float(np.median([r["onset_interval"] for r in ptrain])),
            median_heldout=float(np.median([r["onset_interval"] for r in pheld]))),
        same_direction_onset_phase_heldout=dict(
            n=len(oph),
            median_absolute_windings=float(np.median([r["absolute_windings"] for r in oph])),
            median_directional_coherence=float(np.median([r["directional_coherence"] for r in oph])),
            fraction_0p35_to_0p65_windings=float(np.mean([.35<=r["absolute_windings"]<=.65 for r in oph])),
            fraction_0p8_to_1p2_windings=float(np.mean([.8<=r["absolute_windings"]<=1.2 for r in oph]))),
        heldout_orbit_shape=dict(
            median_phase_advance_per_return=float(np.median([abs(r["phase_advance"])/(2*math.pi) for r in held])),
            median_directional_phase_fraction=float(np.median([r["directional_phase_fraction"] for r in held])),
            median_tangent_alignment=float(np.median([r["mean_tangent_alignment"] for r in held])),
            median_radial_sd=float(np.median([r["radial_sd"] for r in held]))),
        encounter_after_section=dict(
            n=len(encounters),fraction=float(np.mean([r["encounter"] for r in encounters])),
            heldout_fraction=float(np.mean([r["encounter"] for r in encounters if r["seed"]>=30]))))


def figures(orbit,events,segments,onset_phase,model,out: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10,8));sample=[r for r in orbit if r["seed"]==30 and r["side"]==1 and 300<=r["t"]<=380];axes[0,0].plot([r["xi"] for r in sample],[r["eta"] for r in sample],".-",ms=2);axes[0,0].axhline(0,color="k",lw=.6);axes[0,0].plot([-6,0],[0,0],"r",lw=3,label="outward section");axes[0,0].set(aspect="equal",xlabel="local coordinate toward gap",ylabel="vertical coordinate",title="Held-out source-COM trajectory");axes[0,0].legend(frameon=False)
    held_period=[r for r in segments if r["seed"]>=30 and r["same_chirality"] and r["membership_changes"]==0];held_map=[r for r in segments if r["seed"]>=30 and r["same_chirality"] and abs(r["alpha"])<math.pi*.95];axes[0,1].hist([r["return_time"] for r in held_period],bins=25,alpha=.7,label="returns without sampled n change");axes[0,1].hist([r["onset_interval"] for r in onset_phase if r["seed"]>=30],bins=25,alpha=.6,label="same-side no-entry onsets");axes[0,1].set(xlabel="time",ylabel="count",title="Physical period and approach recurrence");axes[0,1].legend(frameon=False)
    axes[1,0].scatter([r["r"] for r in held_map],[r["next_r"] for r in held_map],s=10,alpha=.5);lo=min([r["r"] for r in held_map]+[r["next_r"] for r in held_map]);hi=max([r["r"] for r in held_map]+[r["next_r"] for r in held_map]);axes[1,0].plot([lo,hi],[lo,hi],"k:");axes[1,0].set(xlabel="current section radius",ylabel="next section radius",title="Held-out radial return pairs")
    axes[1,1].scatter([r["alpha"] for r in held_map],[r["next_alpha"] for r in held_map],s=10,alpha=.5);axes[1,1].axhline(0,color="k",ls=":");axes[1,1].axvline(0,color="k",ls=":");axes[1,1].set(xlabel="current relative heading",ylabel="next relative heading",title="Held-out heading return pairs");fig.tight_layout();fig.savefig(out/"poincare_return_validation.png",dpi=220);plt.close(fig)


def self_test() -> None:
    square=[(-1,-1),(1,-1),(1,1),(-1,1)];left=polygon_centroid(clip_half(square,-1));right=polygon_centroid(clip_half(square,1));assert np.allclose(left,(2,-.5,0)) and np.allclose(right,(2,.5,0));assert abs(wrap(3.2)-(-3.083185307179586))<1e-12;print("PASS: exact half-domain centroids, wrapped phase, and reflection orientation")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--source",type=Path);p.add_argument("--attempts",type=Path);p.add_argument("--cache",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.raw_dir,a.source,a.attempts,a.output)):p.error("--raw-dir, --source, --attempts and --output required")
    a.output.mkdir(parents=True,exist_ok=True);geom=production_geometry(a.source);orbit=sample_orbit(a.raw_dir,geom["center_abs"],a.cache);attempts,pairs=attempt_pairs(a.attempts);events,segments=section_events(orbit);phase=onset_phase_audit(pairs,orbit);encounters=associate_encounters(events,attempts);model,scores=map_validation(segments,encounters);shuffle=shuffle_control(segments);physical=summaries(events,segments,phase,pairs,encounters);physical["post_outcome_cycles"]=post_outcome_cycles(attempts);figures(orbit,events,segments,phase,model,a.output)
    write_csv(a.output/"section_events.csv",events);write_csv(a.output/"section_return_pairs.csv",segments);write_csv(a.output/"no_entry_onset_phase.csv",phase);write_csv(a.output/"section_encounter_rewards.csv",encounters);write_csv(a.output/"heldout_map_scores_by_seed.csv",scores)
    summary=dict(configuration=dict(section="outward ray from exact chamber area centroid; both chiralities oriented by crossing direction",state="(COM radius, mean-heading angle relative to oriented tangent)",train_seeds="0-29",heldout_seeds="30-39",population_role="bookkeeping only",new_simulations=0,nonlinear_map=0),production_geometry={k:v for k,v in geom.items() if not isinstance(v,np.ndarray)},audits=dict(source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),n_orbit_samples=len(orbit),n_no_entry_pairs=len(pairs),expected_no_entry_pairs=1195),physical_orbit=physical,linear_return_test=model,shuffled_return_control=shuffle)
    (a.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":main()
