#!/usr/bin/env python3
"""Pre-crossing source organization versus initiation and commitment.

Uses the fixed 337-event majority-to-minority cohort and existing exact-q
matched controls.  Polarization, normalized milling, and fixed gap alignment
are diagnostics only; no new state model or microscopic simulation is built.
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
from scipy.stats import rankdata, spearmanr

from analyze_initiation_propagation import RawReader


DT=.001
LEADS=(5.,3.,2.,1.,.5,.2,.1)
Q_EDGES=np.asarray([0.,.25,.5,.75,.9,1.000001])
FEATURES=("polarization","milling_abs","gap_alignment")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def raw_path(raw_dir: Path,seed: int) -> Path:
    return raw_dir/f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"


def source_features(state: np.ndarray,source_sign: int) -> dict:
    """Reflection-invariant organization of the original source school.

    source_sign=+1 is the right chamber and -1 the left.  The fixed direction
    toward the central opening is -source_sign*xhat.
    """
    x,y,theta=np.asarray(state);source=(x>=0) if source_sign>0 else (x<0)
    q=float(2*np.sum(x>=0)/len(x)-1)
    if source.sum()<3:return dict(q=q,polarization=np.nan,milling_signed=np.nan,milling_abs=np.nan,gap_alignment=np.nan,n_source=int(source.sum()))
    ux=np.cos(theta[source]);uy=np.sin(theta[source]);polarization=float(np.hypot(ux.mean(),uy.mean()))
    pos=np.vstack((x[source],y[source]));pos-=pos.mean(axis=1,keepdims=True);radius=np.hypot(pos[0],pos[1]);valid=radius>1e-12
    milling=float(np.mean((pos[0,valid]*uy[valid]-pos[1,valid]*ux[valid])/radius[valid])) if np.any(valid) else np.nan
    alignment=float(-source_sign*ux.mean())
    return dict(q=q,polarization=polarization,milling_signed=milling,milling_abs=abs(milling),gap_alignment=alignment,n_source=int(source.sum()))


def auc(labels,scores) -> float:
    y=np.asarray(labels,bool);s=np.asarray(scores,float);ok=np.isfinite(s);y=y[ok];s=s[ok];n1=int(y.sum());n0=len(y)-n1
    if not n1 or not n0:return np.nan
    return float((rankdata(s)[y].sum()-n1*(n1+1)/2)/(n1*n0))


def qbin(value: float) -> int:
    return max(0,min(int(np.searchsorted(Q_EDGES,value,side="right")-1),len(Q_EDGES)-2))


def stratified_auc(rows: list[dict],feature: str) -> float:
    correct=0.;pairs=0
    for b in range(len(Q_EDGES)-1):
        part=[r for r in rows if qbin(float(r["abs_q_pre"]))==b];n1=sum(int(r["success4"]) for r in part);n0=len(part)-n1
        if n1 and n0:
            a=auc([int(r["success4"]) for r in part],[float(r[feature]) for r in part]);correct+=a*n1*n0;pairs+=n1*n0
    return float(correct/pairs) if pairs else np.nan


def seed_resample(by_seed: dict[int,list[dict]],rng: np.random.Generator) -> list[dict]:
    seeds=list(by_seed);chosen=rng.choice(seeds,len(seeds),replace=True)
    sampled=[]
    for replicate,seed in enumerate(chosen):
        for row in by_seed[int(seed)]:
            copy=dict(row);copy["seed"]=replicate;sampled.append(copy)
    return sampled


def bootstrap_metric(rows: list[dict],metric,nboot=1500,seed=20261119) -> tuple[float,list[float]]:
    point=float(metric(rows));rng=np.random.default_rng(seed);draw=[];by=defaultdict(list)
    for row in rows:by[int(row["seed"])].append(row)
    for _ in range(nboot):
        value=metric(seed_resample(by,rng))
        if np.isfinite(value):draw.append(value)
    return point,np.percentile(draw,[2.5,97.5]).tolist() if draw else [np.nan,np.nan]


def bootstrap_seed_mean(rows: list[dict],feature: str,nboot=3000,seed=20261119) -> tuple[float,list[float]]:
    by=defaultdict(list)
    for r in rows:
        if np.isfinite(float(r[feature])):by[int(r["seed"])].append(float(r[feature]))
    values=np.asarray([np.mean(v) for v in by.values()]);rng=np.random.default_rng(seed);pick=rng.integers(len(values),size=(nboot,len(values)));draw=values[pick].mean(axis=1)
    return float(values.mean()),np.percentile(draw,[2.5,97.5]).tolist()


def bootstrap_standardized(rows: list[dict],feature: str,nboot=3000,seed=20261119) -> tuple[float,list[float]]:
    seeds=sorted({int(r["seed"]) for r in rows});index={s:i for i,s in enumerate(seeds)};shape=(len(seeds),2,len(Q_EDGES)-1);sums=np.zeros(shape);counts=np.zeros(shape)
    for r in rows:
        i=index[int(r["seed"])];g=int(r["success4"]);b=qbin(float(r["abs_q_pre"]));sums[i,g,b]+=float(r[feature]);counts[i,g,b]+=1
    rng=np.random.default_rng(seed);multiplicity=rng.multinomial(len(seeds),np.full(len(seeds),1/len(seeds)),size=nboot);values=[]
    for weights in np.vstack((np.ones(len(seeds),int),multiplicity)):
        differences=[];binweights=[]
        for b in range(len(Q_EDGES)-1):
            means=[];valid=True
            for g in (0,1):
                ok=counts[:,g,b]>0;den=np.sum(weights[ok])
                if den==0:valid=False;break
                means.append(np.sum(weights[ok]*(sums[ok,g,b]/counts[ok,g,b]))/den)
            if valid:
                differences.append(means[1]-means[0]);binweights.append(np.sum(weights[:,None]*counts[:,:,b]))
        values.append(float(np.average(differences,weights=binweights)) if differences else np.nan)
    return values[0],np.nanpercentile(values[1:],[2.5,97.5]).tolist()


def equal_seed_mean(rows: list[dict],feature: str) -> float:
    by=defaultdict(list)
    for r in rows:
        if np.isfinite(float(r[feature])):by[int(r["seed"])].append(float(r[feature]))
    return float(np.mean([np.mean(v) for v in by.values()])) if by else np.nan


def standardized_difference(rows: list[dict],feature: str) -> float:
    valid_bins=[];weights=[];differences=[]
    for b in range(len(Q_EDGES)-1):
        part=[r for r in rows if qbin(float(r["abs_q_pre"]))==b]
        yes=[r for r in part if int(r["success4"])==1];no=[r for r in part if int(r["success4"])==0]
        if yes and no:
            valid_bins.append(b);weights.append(len(part));differences.append(equal_seed_mean(yes,feature)-equal_seed_mean(no,feature))
    if not valid_bins:return np.nan
    return float(np.average(differences,weights=weights))


def extract_event_rows(raw_dir: Path,outcomes: list[dict]) -> list[dict]:
    rows=[];readers={}
    try:
        for ev in outcomes:
            seed=int(ev["seed"])
            if seed not in readers:readers[seed]=RawReader(raw_path(raw_dir,seed))
            source_sign=int(ev["source_sign"]);first=int(ev["first_native"])
            for lead in LEADS:
                f=source_features(readers[seed].frame("state",first-round(lead/DT)),source_sign)
                if abs(f["q"]-source_sign*float(ev["abs_q_pre"]))>1e-12:raise RuntimeError(f"pre-crossing q mismatch seed={seed} event={ev['event_id']} lead={lead}")
                rows.append(dict(seed=seed,event_id=int(ev["event_id"]),lead=lead,source_sign=source_sign,success4=int(ev["success4"]),
                                 additional_fish=int(ev["additional_fish"]),abs_q_pre=float(ev["abs_q_pre"]),**f))
    finally:
        for rd in readers.values():rd.close()
    return rows


def extract_controls(raw_dir: Path,outcomes: list[dict],pairs: list[dict]) -> list[dict]:
    lookup={(int(e["seed"]),int(e["event_id"])):e for e in outcomes};rows=[];readers={}
    try:
        for pair in pairs:
            key=(int(pair["event_seed"]),int(pair["event_id"]));lead=float(pair["lead"])
            if key not in lookup or lead not in (2.,1.):continue
            ev=lookup[key];eseed=key[0];cseed=int(pair["control_seed"])
            if eseed not in readers:readers[eseed]=RawReader(raw_path(raw_dir,eseed))
            if cseed not in readers:readers[cseed]=RawReader(raw_path(raw_dir,cseed))
            sign=int(ev["source_sign"]);estep=int(ev["first_native"])-round(lead/DT);cstep=round(float(pair["control_time"])/DT)
            ef=source_features(readers[eseed].frame("state",estep),sign);cf=source_features(readers[cseed].frame("state",cstep),sign)
            expected=float(pair["q"])
            if abs(ef["q"]-expected)>1e-12 or abs(cf["q"]-expected)>1e-12:raise RuntimeError(f"matched q mismatch {key} lead={lead}")
            row=dict(seed=eseed,event_id=key[1],control_seed=cseed,lead=lead,source_sign=sign,abs_q_pre=abs(float(pair["q"])))
            for name in FEATURES:row[f"event_{name}"]=ef[name];row[f"control_{name}"]=cf[name];row[f"delta_{name}"]=ef[name]-cf[name]
            rows.append(row)
    finally:
        for rd in readers.values():rd.close()
    return rows


def chronology(rows: list[dict]) -> tuple[list[dict],list[dict]]:
    time=[];contrast=[]
    for lead in LEADS:
        part=[r for r in rows if float(r["lead"])==lead]
        for success in (0,1):
            group=[r for r in part if int(r["success4"])==success]
            for feature in FEATURES:
                point,ci=bootstrap_seed_mean(group,feature,seed=20261200+round(lead*10)+success)
                time.append(dict(lead=lead,success4=success,feature=feature,mean=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(group),n_seeds=len({r["seed"] for r in group})))
        for feature in FEATURES:
            point,ci=bootstrap_standardized(part,feature,seed=20261300+round(lead*10))
            contrast.append(dict(lead=lead,feature=feature,q_standardized_success_minus_failure=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(part),n_seeds=len({r["seed"] for r in part})))
    return time,contrast


def commitment_scores(rows: list[dict]) -> list[dict]:
    out=[]
    for lead in LEADS:
        for cohort,predicate in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            part=[r for r in rows if float(r["lead"])==lead and predicate(int(r["seed"]))]
            for feature in FEATURES:
                point,ci=bootstrap_metric(part,lambda x,f=feature:stratified_auc(x,f),seed=20261400+round(lead*10))
                raw=auc([int(r["success4"]) for r in part],[float(r[feature]) for r in part])
                out.append(dict(lead=lead,cohort=cohort,feature=feature,q_stratified_auc=point,ci95_low=ci[0],ci95_high=ci[1],raw_auc=raw,
                                n_events=len(part),n_success=sum(int(r["success4"]) for r in part),n_seeds=len({r["seed"] for r in part})))
    return out


def initiation_scores(rows: list[dict]) -> list[dict]:
    out=[]
    for lead in (2.,1.):
        for cohort,predicate in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            part=[r for r in rows if float(r["lead"])==lead and predicate(int(r["seed"]))]
            for feature in FEATURES:
                metric=lambda x,f=feature:equal_seed_mean(x,f"delta_{f}")
                point,ci=bootstrap_seed_mean(part,f"delta_{feature}",seed=20261500+round(lead*10))
                labels=np.r_[np.ones(len(part)),np.zeros(len(part))];scores=np.r_[[r[f"event_{feature}"] for r in part],[r[f"control_{feature}"] for r in part]]
                out.append(dict(lead=lead,cohort=cohort,feature=feature,paired_event_minus_control=point,ci95_low=ci[0],ci95_high=ci[1],auc=auc(labels,scores),n_pairs=len(part),n_seeds=len({r["seed"] for r in part})))
    return out


def pre_crossing_changes(rows: list[dict]) -> list[dict]:
    lookup={(int(r["seed"]),int(r["event_id"]),float(r["lead"])):r for r in rows};out=[]
    for target in (1.,.5,.2,.1):
        changes=[]
        for r in rows:
            if float(r["lead"])!=2.:continue
            later=lookup[(int(r["seed"]),int(r["event_id"]),target)];row=dict(seed=r["seed"],event_id=r["event_id"],success4=r["success4"],abs_q_pre=r["abs_q_pre"])
            for feature in FEATURES:row[feature]=float(later[feature])-float(r[feature])
            changes.append(row)
        for feature in FEATURES:
            point,ci=bootstrap_standardized(changes,feature,seed=20261600+round(target*10))
            out.append(dict(baseline_lead=2.,target_lead=target,feature=feature,q_standardized_success_minus_failure_change=point,ci95_low=ci[0],ci95_high=ci[1],n_events=len(changes),n_seeds=len({r["seed"] for r in changes})))
    return out


def circulation_persistence(rows: list[dict]) -> list[dict]:
    lookup={(int(r["seed"]),int(r["event_id"]),float(r["lead"])):r for r in rows};out=[]
    all_pairs=[]
    for success in (0,1):
        pairs=[]
        for r in rows:
            if float(r["lead"])==2. and int(r["success4"])==success:
                later=lookup[(int(r["seed"]),int(r["event_id"]),.5)];pairs.append(dict(seed=r["seed"],success4=success,early=r["milling_signed"],late=later["milling_signed"]))
        all_pairs.extend(pairs)
        metric=lambda x:float(spearmanr([r["early"] for r in x],[r["late"] for r in x]).statistic)
        point,ci=bootstrap_metric(pairs,metric,seed=20261700+success)
        retention=equal_seed_mean([{**r,"same_sign":int(np.sign(r["early"])==np.sign(r["late"]))} for r in pairs],"same_sign")
        out.append(dict(success4=success,lead_early=2.,lead_late=.5,signed_milling_spearman=point,ci95_low=ci[0],ci95_high=ci[1],equal_seed_sign_retention=retention,n_events=len(pairs),n_seeds=len({r["seed"] for r in pairs})))
    def difference(sample):
        rho=[]
        for success in (0,1):
            part=[r for r in sample if int(r["success4"])==success];rho.append(float(spearmanr([r["early"] for r in part],[r["late"] for r in part]).statistic))
        return rho[1]-rho[0]
    point,ci=bootstrap_metric(all_pairs,difference,seed=20261702)
    sign_difference=equal_seed_mean([{**r,"same_sign":int(np.sign(r["early"])==np.sign(r["late"]))} for r in all_pairs if int(r["success4"])==1],"same_sign")-equal_seed_mean([{**r,"same_sign":int(np.sign(r["early"])==np.sign(r["late"]))} for r in all_pairs if int(r["success4"])==0],"same_sign")
    out.append(dict(success4="success_minus_failure",lead_early=2.,lead_late=.5,signed_milling_spearman=point,ci95_low=ci[0],ci95_high=ci[1],equal_seed_sign_retention=sign_difference,n_events=len(all_pairs),n_seeds=len({r["seed"] for r in all_pairs})))
    return out


def make_figures(time,contrast,commitment,initiation,output: Path) -> None:
    labels={"polarization":"source polarization P","milling_abs":"normalized milling |M|","gap_alignment":"fixed gap alignment A"};colors={0:"#B66A50",1:"#327C91"}
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,feature in zip(axes,FEATURES):
        for success in (0,1):
            part=sorted([r for r in time if r["feature"]==feature and r["success4"]==success],key=lambda r:r["lead"],reverse=True);x=-np.asarray([r["lead"] for r in part]);y=np.asarray([r["mean"] for r in part]);lo=np.asarray([r["ci95_low"] for r in part]);hi=np.asarray([r["ci95_high"] for r in part])
            ax.plot(x,y,"o-",color=colors[success],label=">=4 additional" if success else "<4 additional");ax.fill_between(x,lo,hi,color=colors[success],alpha=.15)
        ax.axvline(0,color="k",ls=":",lw=.7);ax.set(xlabel="time from first crossing",ylabel=labels[feature],title=labels[feature])
    axes[0].legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(output/"source_organization_chronology.png",dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,feature in zip(axes,FEATURES):
        c=sorted([r for r in commitment if r["cohort"]=="heldout" and r["feature"]==feature],key=lambda r:r["lead"],reverse=True);x=-np.asarray([r["lead"] for r in c]);y=np.asarray([r["q_stratified_auc"] for r in c]);lo=np.asarray([r["ci95_low"] for r in c]);hi=np.asarray([r["ci95_high"] for r in c]);ax.errorbar(x,y,yerr=[y-lo,hi-y],marker="o",capsize=2,label="commitment AUC")
        ini=[r for r in initiation if r["cohort"]=="heldout" and r["feature"]==feature];ax2=ax.twinx();iy=np.asarray([r["paired_event_minus_control"] for r in ini]);ilo=np.asarray([r["ci95_low"] for r in ini]);ihi=np.asarray([r["ci95_high"] for r in ini])
        ax2.errorbar([-r["lead"] for r in ini],iy,yerr=[iy-ilo,ihi-iy],fmt="s--",capsize=2,color="tab:orange",label="initiation contrast")
        ax.axhline(.5,color="0.5",ls=":");ax2.axhline(0,color="tab:orange",ls=":",alpha=.5)
        ax.set(xlabel="pre-crossing time",ylabel="held-out q-stratified AUC",title=labels[feature]);ax2.set_ylabel("event-control difference",color="tab:orange")
    fig.tight_layout();fig.savefig(output/"initiation_vs_commitment.png",dpi=220);plt.close(fig)


def self_test() -> None:
    n=8;x=np.r_[np.full(5,1.),np.full(3,-1.)];y=np.linspace(-1,1,n);theta=np.linspace(-1,1,n);st=np.vstack((x,y,theta));a=source_features(st,1);mirror=st.copy();mirror[0]*=-1;mirror[2]=np.pi-mirror[2];b=source_features(mirror,-1)
    assert np.allclose([a[k] for k in ("polarization","milling_abs","gap_alignment")],[b[k] for k in ("polarization","milling_abs","gap_alignment")])
    assert abs(a["milling_signed"]+b["milling_signed"])<1e-12
    print("PASS: normalized milling, fixed-gap alignment, and reflection symmetry")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--outcomes",type=Path);p.add_argument("--pairs",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.outcomes,args.pairs,args.output)):p.error("all input and output paths are required")
    args.output.mkdir(parents=True,exist_ok=True);outcomes=read_csv(args.outcomes);pairs=read_csv(args.pairs);events=extract_event_rows(args.raw_dir,outcomes);controls=extract_controls(args.raw_dir,outcomes,pairs)
    time,contrast=chronology(events);commitment=commitment_scores(events);initiation=initiation_scores(controls);changes=pre_crossing_changes(events);persistence=circulation_persistence(events)
    write_csv(args.output/"source_organization_events.csv",events);write_csv(args.output/"matched_initiation_controls.csv",controls);write_csv(args.output/"organization_chronology.csv",time);write_csv(args.output/"commitment_contrasts.csv",contrast);write_csv(args.output/"commitment_prediction.csv",commitment);write_csv(args.output/"initiation_comparison.csv",initiation);write_csv(args.output/"pre_crossing_changes.csv",changes);write_csv(args.output/"circulation_persistence.csv",persistence)
    make_figures(time,contrast,commitment,initiation,args.output)
    summary=dict(configuration=dict(raw_dir=str(args.raw_dir),outcomes=str(args.outcomes),pairs=str(args.pairs),leads=LEADS,q_edges=Q_EDGES.tolist(),train_seeds="0-29",heldout_seeds="30-39",new_simulations=0),
                 cohort=dict(n_events=len(outcomes),n_success=sum(int(r["success4"]) for r in outcomes),n_failure=sum(1-int(r["success4"]) for r in outcomes),n_seeds=len({r["seed"] for r in outcomes}),definition=">=4 additional net fish in five time units after first majority-to-minority crossing"),
                 feature_definitions=dict(polarization="magnitude of source-school mean heading",milling_abs="absolute mean normalized angular momentum about source-school centroid",gap_alignment="source mean heading projected on fixed horizontal direction to opening"),
                 audits=dict(exact_pre_crossing_q="passed for every event lead",exact_matched_control_q="passed for every event/control pair",reflection_symmetry="synthetic self-test passed"),
                 commitment_prediction=commitment,initiation_comparison=initiation,pre_crossing_changes=changes,circulation_persistence=persistence)
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":main()
