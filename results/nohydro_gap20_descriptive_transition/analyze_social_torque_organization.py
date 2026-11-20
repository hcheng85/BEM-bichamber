#!/usr/bin/env python3
"""Exact pair organization of the production visual-social torque.

Decomposes the source collective turning rate into reciprocal-cancellation,
heading-mismatch, directed-visibility, and external-neighbor terms.  Leading,
convex-hull, and approach-region receiver splits are diagnostics only.  No
turning law or transition model is fit.
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
from scipy.spatial import ConvexHull
from scipy.stats import spearmanr

from analyze_approach_current_balance import DT, ELL, N, fast_check_voronoi, load_gap20_production, raw_path
from analyze_collective_polarization_evolution import IN_NOISE, source_mask
from analyze_geometric_first_passage_reduction import C_EDGES, D_EDGES, Q_EDGES, bindex, geometry, seed_boot
from analyze_initiation_propagation import RawReader


LEADS=(2.,1.,.5,.2)
COMPONENTS=("internal_alignment_visual_asymmetry","internal_alignment_normalization_asymmetry","internal_position_heading_mismatch","internal_position_visibility_asymmetry","external_alignment","external_position")


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def cross(a,b):return float(a[0]*b[1]-a[1]*b[0])


def cohort(path: Path) -> list[dict]:
    out=[]
    for r in read_csv(path):out.append(dict(seed=int(r["seed"]),event_id=int(r["event_id"]),first_native=int(r["first_native"]),source_sign=int(r["source_sign"]),direction=-int(r["source_sign"]),success4=int(r["success4"]),abs_q_pre=float(r["abs_q_pre"])))
    if len(out)!=337:raise RuntimeError("Expected fixed 337-event cohort")
    return out


def social_pair_budget(mod,para: dict,state: np.ndarray,mask: np.ndarray,z: dict) -> tuple[dict,dict]:
    vn=fast_check_voronoi(mod,mod.GetVoronoi(state[:2],para),state[:2]);x=state[:2].T;theta=state[2];e=np.column_stack((np.cos(theta),np.sin(theta)));ns=int(mask.sum())
    visual=np.zeros((N,N));align=np.zeros((N,N));position=np.zeros((N,N))
    for i,j in zip(*np.nonzero(vn)):
        rho=x[j]-x[i];dist=float(np.linalg.norm(rho));unit=rho/dist;visual[i,j]=1+float(e[i]@unit);align[i,j]=para["Ip"]*math.sin(theta[j]-theta[i]);position[i,j]=cross(e[i],rho)
    normalizer=visual.sum(axis=1);perfish=np.divide(np.sum((align+position)*visual,axis=1),normalizer,out=np.zeros(N),where=normalizer>0)
    with np.errstate(invalid="ignore",divide="ignore"),contextlib.redirect_stdout(io.StringIO()):native=np.asarray(mod.GetVision(vn,state,para)).reshape(-1)
    coeff=np.zeros(N);coeff[mask]=np.cos(theta[mask]-z["theta"])/(ns*z["p"])
    exact=float(coeff@native);components={k:0. for k in COMPONENTS};internal_l1=0.;internal_pair_l1=0.
    for i in np.flatnonzero(mask):
        for j in np.flatnonzero(vn[i]):
            if mask[j] and j<i:continue
            ki=coeff[i]*visual[i,j]/normalizer[i] if normalizer[i]>0 else 0.
            if mask[j]:
                kj=coeff[j]*visual[j,i]/normalizer[j] if normalizer[j]>0 else 0.;rho=x[j]-x[i];av=align[i,j]
                ai=coeff[i]/normalizer[i] if normalizer[i]>0 else 0.;aj=coeff[j]/normalizer[j] if normalizer[j]>0 else 0.
                tav=av*.5*(ai+aj)*(visual[i,j]-visual[j,i]);tan=av*.5*(visual[i,j]+visual[j,i])*(ai-aj)
                th=.5*(ki+kj)*cross(e[i]-e[j],rho);tv=.5*(ki-kj)*cross(e[i]+e[j],rho)
                components["internal_alignment_visual_asymmetry"]+=tav;components["internal_alignment_normalization_asymmetry"]+=tan;components["internal_position_heading_mismatch"]+=th;components["internal_position_visibility_asymmetry"]+=tv
                tij=ki*(align[i,j]+position[i,j]);tji=kj*(align[j,i]+position[j,i]);internal_l1+=abs(tij)+abs(tji);internal_pair_l1+=abs(tij+tji)
            else:
                components["external_alignment"]+=ki*align[i,j];components["external_position"]+=ki*position[i,j]
    receiver=coeff*perfish;rel=(x-z["r"])@z["g"];leading=mask&(rel>=0);trailing=mask&~leading
    hull=np.zeros(N,bool)
    src=np.flatnonzero(mask)
    if len(src)>=3:hull[src[ConvexHull(x[src]).vertices]]=True
    approach=mask&(np.sign(z["r"][0])*state[0]>=0)&(np.sign(z["r"][0])*state[0]<=ELL)&(np.abs(state[1])<=ELL/2)
    # Coherent-heading spatial approximation on the exact graph and positions.
    coh_visual=np.zeros((N,N));coh_align=np.zeros((N,N));coh_position=np.zeros((N,N))
    for i in src:
        for j in np.flatnonzero(vn[i]):
            rho=x[j]-x[i];dist=float(np.linalg.norm(rho));coh_visual[i,j]=1+float(z["m"]@(rho/dist));coh_position[i,j]=cross(z["m"],rho);coh_align[i,j]=0. if mask[j] else para["Ip"]*math.sin(theta[j]-z["theta"])
    coh_z=coh_visual.sum(axis=1);coherent=float(sum(np.sum((coh_align[i]+coh_position[i])*coh_visual[i])/coh_z[i] for i in src if coh_z[i]>0)/ns)
    rec=dict(leading_receiver=float(receiver[leading].sum()),trailing_receiver=float(receiver[trailing].sum()),hull_receiver=float(receiver[hull].sum()),interior_receiver=float(receiver[mask&~hull].sum()),
             approach_receiver=float(receiver[approach].sum()),nonapproach_receiver=float(receiver[mask&~approach].sum()),leading_fraction=float(leading.sum()/ns),hull_fraction=float(hull.sum()/ns),approach_fraction=float(approach.sum()/ns))
    audit=dict(max_perfish_social_residual=float(np.max(np.abs(perfish-native))),collective_social=exact,reconstruction_residual=exact-sum(components.values()),
               internal_pair_cancellation=1-internal_pair_l1/internal_l1 if internal_l1 else np.nan,internal_directed_l1=internal_l1,coherent_spatial_approximation=coherent,coherent_spatial_error=exact-coherent,
               n_visible_directed=int(vn.sum()),n_external_directed=int(np.sum(vn[mask][:,~mask])))
    return {**components,**rec},audit


def storage_current(state: np.ndarray,velocity: np.ndarray,direction: int) -> tuple[float,float,int]:
    s=float(direction);x,y=state[0],state[1];inside=(s*x<0)&((-s*x)<=ELL)&(np.abs(y)<=ELL/2);u=-s*x
    b=float(np.sum((ELL-u[inside]))/(N*ELL));ja=float(np.sum(s*velocity[0,inside])/N)
    return b,ja,int(inside.sum())


def extract(mod,para: dict,raw_dir: Path,events: list[dict]) -> tuple[list[dict],dict]:
    rows=[];by=defaultdict(list);max_native=0.;max_recon=0.
    for ev in events:by[ev["seed"]].append(ev)
    for count,seed in enumerate(sorted(by),1):
        print(f"social pair audit {count}/{len(by)} seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed))
        try:
            for ev in by[seed]:
                for lead in LEADS:
                    k=ev["first_native"]-round(lead/DT);state=reader.frame("state",k);velocity=reader.frame("rdot",k+1);mask=source_mask(state,ev["direction"]);z=geometry(state,mask)
                    terms,audit=social_pair_budget(mod,para,state,mask,z);max_native=max(max_native,audit["max_perfish_social_residual"]);max_recon=max(max_recon,abs(audit["reconstruction_residual"]));b,ja,napp=storage_current(state,velocity,ev["direction"])
                    align_factor=-math.sin(z["psi"]);aligned={"alignment_"+name:align_factor*terms[name] for name in COMPONENTS}
                    for name in ("leading_receiver","trailing_receiver","hull_receiver","interior_receiver","approach_receiver","nonapproach_receiver"):aligned["alignment_"+name]=align_factor*terms[name]
                    aligned["alignment_social_total"]=align_factor*audit["collective_social"];aligned["alignment_coherent_spatial_approximation"]=align_factor*audit["coherent_spatial_approximation"]
                    rows.append(dict(**ev,lead=lead,d=z["d"],psi=z["psi"],cospsi=z["cospsi"],polarization=z["p"],b=b,Ja=ja,n_approach=napp,**terms,**aligned,**audit))
        finally:reader.close()
    return rows,dict(max_abs_native_social_reconstruction=max_native,max_abs_collective_pair_reconstruction=max_recon,n_snapshots=len(rows))


def cell_contrast(rows: list[dict],field: str) -> float:
    cells=defaultdict(list)
    for r in rows:cells[(bindex(r["abs_q_pre"],Q_EDGES),bindex(r["d"],D_EDGES),bindex(r["cospsi"],C_EDGES))].append(r)
    values=[];weights=[]
    for part in cells.values():
        yes=[r[field] for r in part if r["success4"]];no=[r[field] for r in part if not r["success4"]]
        if yes and no:values.append(float(np.mean(yes)-np.mean(no)));weights.append(len(part))
    return float(np.average(values,weights=weights)) if values else np.nan


def summaries(rows: list[dict]) -> tuple[list[dict],list[dict],dict]:
    contrasts=[];diagnostics=[]
    receiver_fields=["alignment_"+x for x in ("leading_receiver","trailing_receiver","hull_receiver","interior_receiver","approach_receiver","nonapproach_receiver")]
    fields=["alignment_"+x for x in COMPONENTS]+receiver_fields+["alignment_social_total","alignment_coherent_spatial_approximation","Ja","b"]
    for lead in LEADS:
        part=[r for r in rows if r["lead"]==lead]
        for cohort_name,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
            selected=[r for r in part if pred(r["seed"])]
            for field in fields:
                point,ci=seed_boot(selected,lambda x,f=field:cell_contrast(x,f),20261180+round(lead*10))
                contrasts.append(dict(lead=lead,cohort=cohort_name,term=field,continued_minus_failed=point,ci95_low=ci[0],ci95_high=ci[1],n=len(selected),n_seeds=len({r["seed"] for r in selected})))
    for cohort_name,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
        selected=[r for r in rows if pred(r["seed"])]
        for field in ("internal_pair_cancellation","leading_fraction","hull_fraction","approach_fraction"):
            point,ci=seed_boot(selected,lambda x,f=field:float(np.mean([r[f] for r in x])),20261191)
            diagnostics.append(dict(cohort=cohort_name,metric=field,value=point,ci95_low=ci[0],ci95_high=ci[1],n=len(selected)))
        exact=np.asarray([r["collective_social"] for r in selected]);approx=np.asarray([r["coherent_spatial_approximation"] for r in selected])
        diagnostics.append(dict(cohort=cohort_name,metric="coherent_spatial_approximation",value=float(np.corrcoef(exact,approx)[0,1]),rmse=float(np.sqrt(np.mean((exact-approx)**2))),actual_sd=float(np.std(exact)),n=len(selected)))
        for split,fraction in (("leading_receiver","leading_fraction"),("hull_receiver","hull_fraction"),("approach_receiver","approach_fraction")):
            share=np.asarray([abs(r[split])/(abs(r[split])+abs(r[{"leading_receiver":"trailing_receiver","hull_receiver":"interior_receiver","approach_receiver":"nonapproach_receiver"}[split]])+1e-15) for r in selected]);frac=np.asarray([r[fraction] for r in selected])
            diagnostics.append(dict(cohort=cohort_name,metric=split+"_absolute_share",value=float(np.mean(share)),population_fraction=float(np.mean(frac)),enrichment=float(np.mean(share)/np.mean(frac)),n=len(selected)))
    # Prespecified supply link: lead-1 social alignment versus Ja change to lead 0.5.
    lookup={(r["seed"],r["event_id"],r["lead"]):r for r in rows};link=[]
    for r in rows:
        if r["lead"]==1.:
            later=lookup[(r["seed"],r["event_id"],.5)];z=dict(r);z["delta_Ja_1_to_0p5"]=later["Ja"]-r["Ja"];link.append(z)
    def residual_spearman(part):
        cells=defaultdict(list)
        for r in part:cells[(bindex(r["abs_q_pre"],Q_EDGES),bindex(r["d"],D_EDGES),bindex(r["cospsi"],C_EDGES))].append(r)
        x=[];y=[]
        for cell in cells.values():
            if len(cell)>=3:
                xx=np.asarray([r["alignment_social_total"] for r in cell]);yy=np.asarray([r["delta_Ja_1_to_0p5"] for r in cell]);x.extend(xx-xx.mean());y.extend(yy-yy.mean())
        return float(spearmanr(x,y).statistic) if len(x)>3 else np.nan
    supply={}
    for name,pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
        selected=[r for r in link if pred(r["seed"])];point,ci=seed_boot(selected,residual_spearman,20261201);supply[name]=dict(residual_spearman=point,ci95=ci,n=len(selected),n_seeds=len({r["seed"] for r in selected}))
    return contrasts,diagnostics,supply


def figures(rows,contrasts,diagnostics,out: Path) -> None:
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    names=("internal_alignment_visual_asymmetry","internal_alignment_normalization_asymmetry","internal_position_heading_mismatch","internal_position_visibility_asymmetry","external_alignment","external_position")
    means=[np.mean([abs(r[n]) for r in rows]) for n in names];axes[0].bar(range(len(names)),means);axes[0].set_xticks(range(len(names)),("align visual\nasym","align norm\nasym","heading\nmismatch","position\nnonrec","external\nalign","external\nposition"),rotation=20);axes[0].set(ylabel="mean absolute angular rate",title="Exact social-torque magnitudes")
    metrics={r["metric"]:r for r in diagnostics if r["cohort"]=="all"};keys=("leading_receiver_absolute_share","hull_receiver_absolute_share","approach_receiver_absolute_share");axes[1].bar(("leading half","convex hull","approach region"),[metrics[k]["value"] for k in keys],label="absolute torque share");axes[1].plot(("leading half","convex hull","approach region"),[metrics[k]["population_fraction"] for k in keys],"ko",label="population share");axes[1].set(ylim=(0,1),title="Where social turning is received");axes[1].legend(frameon=False);fig.tight_layout();fig.savefig(out/"social_pair_organization.png",dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4));terms=("alignment_internal_alignment_visual_asymmetry","alignment_internal_alignment_normalization_asymmetry","alignment_internal_position_visibility_asymmetry")
    x=np.arange(len(LEADS));width=.25
    for j,term in enumerate(terms):
        part=[r for r in contrasts if r["cohort"]=="all" and r["term"]==term];y=np.asarray([r["continued_minus_failed"] for r in part]);lo=np.asarray([r["ci95_low"] for r in part]);hi=np.asarray([r["ci95_high"] for r in part]);ax.bar(x+(j-1)*width,y,width,yerr=[y-lo,hi-y],capsize=2,label=term.replace("alignment_", ""))
    ax.axhline(0,color=".5",ls=":");ax.set_xticks(x,[str(x) for x in LEADS]);ax.set(xlabel="lead before crossing",ylabel="continued - failed alignment rate",title="Geometry-controlled social organization");ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(out/"social_mechanism_contrasts.png",dpi=220);plt.close(fig)


def self_test() -> None:
    # Algebraic pair identity for arbitrary directed coefficients and headings.
    ei=np.asarray([.8,.6]);ej=np.asarray([-.2,math.sqrt(.96)]);rho=np.asarray([1.3,-.7]);ki=.4;kj=.1;A=.9
    direct=ki*(A+cross(ei,rho))+kj*(-A-cross(ej,rho));split=A*(ki-kj)+.5*(ki+kj)*cross(ei-ej,rho)+.5*(ki-kj)*cross(ei+ej,rho)
    assert abs(direct-split)<1e-14
    print("PASS: exact directed-pair social-torque identity")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--cohort",type=Path);p.add_argument("--source",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    if not all((a.raw_dir,a.cohort,a.source,a.output)):p.error("all paths required")
    a.output.mkdir(parents=True,exist_ok=True);events=cohort(a.cohort)
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(a.source)
    para=dict(Ip=4.,In=IN_NOISE,R=mod.R,light=1.,dt=DT,Iw=.94,If=0.,num=N,delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
    rows,audit=extract(mod,para,a.raw_dir,events);contrasts,diagnostics,supply=summaries(rows);write_csv(a.output/"pair_budget_samples.csv",rows);write_csv(a.output/"component_contrasts.csv",contrasts);write_csv(a.output/"organization_diagnostics.csv",diagnostics);figures(rows,contrasts,diagnostics,a.output)
    config=dict(cohort_size=len(events),leads=list(LEADS),train_seeds="0-29",heldout_seeds="30-39",new_simulations=0,fitted_terms=0,production_source=str(a.source),production_source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest())
    (a.output/"summary.json").write_text(json.dumps(dict(configuration=config,audits=audit,organization_diagnostics=diagnostics,component_contrasts=contrasts,supply_link=supply),indent=2)+"\n")


if __name__=="__main__":main()
