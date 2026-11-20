#!/usr/bin/env python3
"""Bounded matched-noise intervention on internal visual-alignment reciprocity.

The production source is read-only.  Internal source-pair alignment weights
alpha_ij=w_ij/Z_i are replaced by (alpha_ij+alpha_ji)/2.  Positional social
terms, external pairs, visual denominators, walls, clipping, collisions, and
noise are untouched.  The exact collective-k counterfactual is audited but is
not simulated because it requires singular global projection feedback.
"""
from __future__ import annotations

import argparse
import ast
import contextlib
import copy
import csv
import hashlib
import io
import json
import math
import multiprocessing as mp
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_approach_current_balance import ELL, fast_check_voronoi, load_gap20_production
from analyze_geometric_first_passage_reduction import geometry
from analyze_initiation_propagation import RawReader
from analyze_social_torque_organization import cross
from run_bottleneck_commitment_branches import DT, N, parameters, raw_path


HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
SOURCE=REPO.parent/"BEM_bichamber_parallel_instance_no_hydro.py"
RAW_DEFAULT=Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
COHORT=HERE/"outputs"/"collective_qab_pathway"/"event_outcomes.csv"
OUT_DEFAULT=HERE/"outputs"/"nonreciprocal_alignment_intervention"
N_STATES=14;N_CLONES=8;HORIZON_STEPS=5000;LEAD_STEPS=1000
_MOD=None;_PARA=None;_STEP=None


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def source_mask(state: np.ndarray,source_sign: int) -> np.ndarray:
    return state[0]>=0 if source_sign>0 else state[0]<0


def visible_arrays(vn: np.ndarray,state: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    theta=state[2];e=np.column_stack((np.cos(theta),np.sin(theta)));visual=np.zeros((N,N));align=np.zeros((N,N))
    for i,j in zip(*np.nonzero(vn)):
        rho=state[:2,j]-state[:2,i];dist=float(np.linalg.norm(rho));visual[i,j]=1+float(e[i]@(rho/dist));align[i,j]=4.*math.sin(theta[j]-theta[i])
    return visual,align,visual.sum(axis=1)


def social_intervention(mod,vn: np.ndarray,state: np.ndarray,para: dict,source_sign: int,mode: str) -> np.ndarray:
    with np.errstate(invalid="ignore",divide="ignore"),contextlib.redirect_stdout(io.StringIO()):full=np.asarray(mod.GetVision(vn,state,para)).reshape(-1)
    if mode=="baseline":return full
    if mode!="symmetric_internal_alignment":raise ValueError(mode)
    visual,align,z=visible_arrays(vn,state);mask=source_mask(state,source_sign);out=full.copy()
    for i in np.flatnonzero(mask):
        for j in np.flatnonzero(vn[i]&mask):
            if j<=i:continue
            ai=visual[i,j]/z[i] if z[i]>0 else 0.;aj=visual[j,i]/z[j] if z[j]>0 else 0.;abar=.5*(ai+aj);a=align[i,j]
            out[i]+=(abar-ai)*a;out[j]+=(aj-abar)*a
    return out


def extract_step(mod,source: Path):
    tree=ast.parse(source.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="step_Computation");loop=next(n for n in fn.body if isinstance(n,ast.For))
    first=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=="CurState");body=copy.deepcopy(loop.body[first:]);nr=sr=0
    for node in body:
        if not isinstance(node,ast.Assign) or not isinstance(node.targets[0],ast.Name):continue
        name=node.targets[0].id
        if name=="wNoise":node.value=ast.Name(id="noise_increment",ctx=ast.Load());nr+=1
        elif name=="wVision":
            node.value=ast.Call(func=ast.Name(id="GetVision_intervention",ctx=ast.Load()),args=[ast.Name(id="voro",ctx=ast.Load()),ast.Name(id="CurState",ctx=ast.Load()),ast.Name(id="para",ctx=ast.Load()),ast.Name(id="source_sign",ctx=ast.Load()),ast.Name(id="mode",ctx=ast.Load())],keywords=[]);sr+=1
    if (nr,sr)!=(1,1):raise RuntimeError("Unique production noise/social statements not found")
    prefix=ast.parse('''def archived_step(current, para, noise_increment, dis_last, source_sign, mode="baseline"):
    state=np.empty((3,para["num"],2));state[:,:,0]=current
    rDotList=np.empty((2,para["num"],2));count=0
''').body[0]
    suffix=ast.parse('''return state[:,:,count].copy(),rDot.copy(),thetaDot.copy(),dis_last.copy(),dict(collision=int(np.sum(np.linalg.norm(rDot,axis=0)==0)),safety_candidates=int(np.sum(dis>=safedis)))''').body
    prefix.body+=body+suffix;module=ast.fix_missing_locations(ast.Module(body=[prefix],type_ignores=[]));text=ast.unparse(module)+"\n";env=dict(vars(mod));env["CheckVoronoi"]=lambda vn,loc,pa:fast_check_voronoi(mod,vn,loc);env["GetVision_intervention"]=lambda vn,st,pa,s,m:social_intervention(mod,vn,st,pa,s,m);exec(compile(module,str(source)+"::reciprocity_adapter","exec"),env)
    return env["archived_step"],text,hashlib.sha256(text.encode()).hexdigest()


def select_states() -> list[dict]:
    rows=read_csv(COHORT);selected=[]
    for seed in sorted({int(r["seed"]) for r in rows})[:N_STATES]:
        r=min((x for x in rows if int(x["seed"])==seed),key=lambda x:int(x["first_native"]));k=int(r["first_native"])-LEAD_STEPS
        selected.append(dict(sample_id=f"lead1_s{seed}_e{r['event_id']}",seed=seed,event_id=int(r["event_id"]),native_step=k,first_native=int(r["first_native"]),source_sign=int(r["source_sign"]),success4_archived=int(r["success4"]),selection_rule="first fixed-cohort crossing per seed at lead 1; archived continuation label unused"))
    return selected


def load_state(root: Path,item: dict) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    rd=RawReader(raw_path(root,item["seed"]));k=item["native_step"];prev=rd.frame("state",k-1);state=rd.frame("state",k);target=rd.frame("state",k+1);rd.close();return state,prev,target


def setup(root: Path,step,para,items: list[dict]) -> list[dict]:
    zero=np.zeros((3,N));out=[]
    for item in items:
        state,prev,_=load_state(root,item);_,_,_,memory,_=step(prev,para,zero,None,item["source_sign"],"baseline");out.append(dict(item,state=state,dis_last=memory,nright0=int(np.sum(state[0]>=0))))
    return out


def collective_pair_audit(mod,para,state: np.ndarray,source_sign: int) -> dict:
    vn=fast_check_voronoi(mod,mod.GetVoronoi(state[:2],para),state[:2]);visual,align,z=visible_arrays(vn,state);mask=source_mask(state,source_sign);g=geometry(state,mask);c=np.zeros(N);c[mask]=np.cos(state[2,mask]-g["theta"])/(mask.sum()*g["p"])
    original=exact_cf=micro_cf=0.;max_exact_pair=0.;modified=social_intervention(mod,vn,state,para,source_sign,"symmetric_internal_alignment")
    baseline=social_intervention(mod,vn,state,para,source_sign,"baseline")
    for i in np.flatnonzero(mask):
        for j in np.flatnonzero(vn[i]&mask):
            if j<=i:continue
            ai=visual[i,j]/z[i] if z[i]>0 else 0.;aj=visual[j,i]/z[j] if z[j]>0 else 0.;ki=c[i]*ai;kj=c[j]*aj;a=align[i,j];original+=a*(ki-kj)
            kbar=.5*(ki+kj);pair_exact=a*(kbar-kbar);exact_cf+=pair_exact;max_exact_pair=max(max_exact_pair,abs(pair_exact));abar=.5*(ai+aj);micro_cf+=a*abar*(c[i]-c[j])
    observed_delta=float(c@(modified-baseline));return dict(original_internal_collective_alignment=original,exact_k_symmetric_collective_alignment=exact_cf,microscopic_alpha_symmetric_collective_alignment=micro_cf,
        microscopic_observed_delta=observed_delta,delta_reconstruction_residual=observed_delta-(micro_cf-original),max_abs_exact_pair_after_symmetry=max_exact_pair,mean_abs_fish_torque_change=float(np.mean(np.abs(modified-baseline))),max_abs_fish_torque_change=float(np.max(np.abs(modified-baseline))))


def audit(root: Path,mod,step,para,origins: list[dict],source_sha: str,adapter_sha: str) -> dict:
    zero=np.zeros((3,N));parity=[];pairs=[]
    for item in origins:
        state,_,target=load_state(root,item);base,_,_,_,_=step(state,para,zero,item["dis_last"],item["source_sign"],"baseline");noise=np.zeros((3,N));noise[2]=target[2]-base[2];pred,vel,_,_,_=step(state,para,noise,item["dis_last"],item["source_sign"],"baseline");rd=RawReader(raw_path(root,item["seed"]));archvel=rd.frame("rdot",item["native_step"]+1);rd.close();parity.append(dict(sample_id=item["sample_id"],max_state_error=float(np.max(np.abs(pred-target))),max_velocity_error=float(np.max(np.abs(vel-archvel)))));pairs.append(dict(sample_id=item["sample_id"],**collective_pair_audit(mod,para,state,item["source_sign"])))
    snapshot_pairs=[];events=read_csv(COHORT);by=defaultdict(list)
    for r in events:by[int(r["seed"])].append(r)
    for seed,items in by.items():
        rd=RawReader(raw_path(root,seed))
        try:
            for r in items:
                state=rd.frame("state",int(r["first_native"])-LEAD_STEPS);snapshot_pairs.append(dict(seed=seed,event_id=int(r["event_id"]),**collective_pair_audit(mod,para,state,int(r["source_sign"]))))
        finally:rd.close()
    return dict(source_sha256=source_sha,adapter_sha256=adapter_sha,parity=dict(passed=max(max(x["max_state_error"],x["max_velocity_error"]) for x in parity)<1e-12,rows=parity),
        instantaneous_counterfactual=dict(exact_k_symmetry="replace k_ij,k_ji by their mean; collective internal pair alignment cancels identically",microscopic_intervention="replace alpha_ij=w_ij/Z_i and alpha_ji by their mean for current internal source pairs",unavoidable_difference="microscopic alpha symmetry leaves residual proportional to alpha_bar*(c_i-c_j); exact k symmetry would require global c_i-dependent division and can be singular",rows=pairs,
        max_abs_exact_collective_after_symmetry=max(abs(x["exact_k_symmetric_collective_alignment"]) for x in pairs),max_abs_delta_reconstruction=max(abs(x["delta_reconstruction_residual"]) for x in pairs),
        mean_abs_original=float(np.mean([abs(x["original_internal_collective_alignment"]) for x in pairs])),mean_abs_microscopic_residual=float(np.mean([abs(x["microscopic_alpha_symmetric_collective_alignment"]) for x in pairs])),
        fixed_cohort_lead1=dict(n_snapshots=len(snapshot_pairs),rows=snapshot_pairs,max_abs_exact_collective_after_symmetry=max(abs(x["exact_k_symmetric_collective_alignment"]) for x in snapshot_pairs),max_abs_delta_reconstruction=max(abs(x["delta_reconstruction_residual"]) for x in snapshot_pairs),
            mean_abs_original=float(np.mean([abs(x["original_internal_collective_alignment"]) for x in snapshot_pairs])),mean_abs_microscopic_residual=float(np.mean([abs(x["microscopic_alpha_symmetric_collective_alignment"]) for x in snapshot_pairs])))),
        preserved_terms=["internal positional social numerator","external alignment and positional social numerator","all original visual denominators","Delaunay/wall visibility","wall torque input and both production clipping rules","collision handling rule","angular noise"])


def init_worker(source: str,adapter_text: str):
    global _MOD,_PARA,_STEP
    _MOD=load_gap20_production(Path(source));_MOD.CheckVoronoi=lambda vn,loc,pa:fast_check_voronoi(_MOD,vn,loc);_PARA=parameters(_MOD);env=dict(vars(_MOD));env["GetVision_intervention"]=lambda vn,st,pa,s,m:social_intervention(_MOD,vn,st,pa,s,m);exec(adapter_text,env);_STEP=env["archived_step"]


def order(state: np.ndarray,source_sign: int) -> tuple[float,float,float,float]:
    mask=source_mask(state,source_sign);z=geometry(state,mask);return z["p"],float(z["h"]@z["g"]),z["theta"],z["psi"]


def storage_current(state: np.ndarray,velocity: np.ndarray,source_sign: int) -> tuple[float,float]:
    direction=-float(source_sign);x,y=state[0],state[1];inside=(direction*x<0)&((-direction*x)<=ELL)&(np.abs(y)<=ELL/2);u=-direction*x;b=float(np.sum(ELL-u[inside])/(N*ELL));ja=float(np.sum(direction*velocity[0,inside])/N);return b,ja


def run_task(task: dict) -> list[dict]:
    rng=np.random.default_rng(task["noise_seed"]);noises=.5*np.sqrt(DT)*rng.normal(size=(HORIZON_STEPS,N));rows=[]
    for mode in ("baseline","symmetric_internal_alignment"):
        state=task["state"].copy();memory=task["dis_last"].copy();p0,a0,_,_=order(state,task["source_sign"]);entered=False;sum_ja=sum_b=sum_p=sum_a=turn=geom=0.;sum_ja1=sum_b1=sum_p1=sum_a1=turn1=geom1=0.;valid_turn=0;valid_turn1=0
        for n in range(HORIZON_STEPS):
            inc=np.zeros((3,N));inc[2]=noises[n];previous=state;oldmask=source_mask(previous,task["source_sign"]);p,a,theta,psi=order(previous,task["source_sign"]);state,velocity,_,memory,_=_STEP(previous,_PARA,inc,memory,task["source_sign"],mode);b,ja=storage_current(previous,velocity,task["source_sign"]);sum_ja+=ja;sum_b+=b;sum_p+=p;sum_a+=a
            if n<1000:sum_ja1+=ja;sum_b1+=b;sum_p1+=p;sum_a1+=a
            if np.any((task["source_sign"]*previous[0]>=0)&(task["source_sign"]*state[0]<0)):entered=True
            newmask=source_mask(state,task["source_sign"])
            if np.array_equal(oldmask,newmask):
                _,_,theta1,psi1=order(state,task["source_sign"]);dtheta=float(np.arctan2(np.sin(theta1-theta),np.cos(theta1-theta)));dpsi=float(np.arctan2(np.sin(psi1-psi),np.cos(psi1-psi)));t=-math.sin(psi)*dtheta;gg=-math.sin(psi)*(dpsi-dtheta);turn+=t;geom+=gg;valid_turn+=1
                if n<1000:turn1+=t;geom1+=gg;valid_turn1+=1
        pT,aT,_,_=order(state,task["source_sign"]);nrightT=int(np.sum(state[0]>=0));net=int(task["source_sign"]*(task["nright0"]-nrightT))
        rows.append(dict(sample_id=task["sample_id"],seed=task["seed"],event_id=task["event_id"],clone=task["clone"],arm=mode,source_sign=task["source_sign"],nright0=task["nright0"],nrightT=nrightT,net_transfer=net,success4=int(net>=4),first_entry=int(entered),
            polarization_initial=p0,polarization_final=pT,gap_alignment_initial=a0,gap_alignment_final=aT,mean_Ja_1=sum_ja1/1000,mean_Ja_5=sum_ja/HORIZON_STEPS,mean_b_1=sum_b1/1000,mean_b_5=sum_b/HORIZON_STEPS,mean_polarization_1=sum_p1/1000,mean_polarization_5=sum_p/HORIZON_STEPS,mean_gap_alignment_1=sum_a1/1000,mean_gap_alignment_5=sum_a/HORIZON_STEPS,
            turning_alignment_rate_1=turn1/(valid_turn1*DT),turning_alignment_rate_5=turn/(valid_turn*DT),geometry_alignment_rate_1=geom1/(valid_turn1*DT),geometry_alignment_rate_5=geom/(valid_turn*DT),turning_valid_fraction=valid_turn/HORIZON_STEPS))
    return rows


def simulate(origins: list[dict],adapter_text: str,output: Path,workers: int) -> list[dict]:
    tasks=[]
    for item in origins:
        for clone in range(N_CLONES):tasks.append(dict(item,clone=clone,noise_seed=202611220000+item["seed"]*100+clone))
    partial=output/"branches_partial.csv";existing=read_csv(partial) if partial.exists() else [];complete={(r["sample_id"],int(r["clone"])) for r in existing if sum(x["sample_id"]==r["sample_id"] and int(x["clone"])==int(r["clone"]) for x in existing)==2};rows=[r for r in existing if (r["sample_id"],int(r["clone"])) in complete];todo=[t for t in tasks if (t["sample_id"],t["clone"]) not in complete];start=time.monotonic()
    if todo:
        with mp.get_context("spawn").Pool(min(workers,len(todo)),initializer=init_worker,initargs=(str(SOURCE),adapter_text)) as pool:
            for i,result in enumerate(pool.imap_unordered(run_task,todo),1):rows.extend(result);write_csv(partial,rows);print(f"completed {i}/{len(todo)} paired blocks; elapsed {time.monotonic()-start:.1f}s",flush=True)
    rows.sort(key=lambda r:(r["sample_id"],int(r["clone"]),r["arm"]));write_csv(output/"branches.csv",rows)
    if partial.exists():partial.unlink()
    return rows


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path,default=RAW_DEFAULT);p.add_argument("--output",type=Path,default=OUT_DEFAULT);p.add_argument("--workers",type=int,default=max(1,min(12,os.cpu_count() or 1)));p.add_argument("--stage",choices=("audit","simulate","all"),default="all");a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    mod=load_gap20_production(SOURCE);para=parameters(mod);step,text,adapter_sha=extract_step(mod,SOURCE);items=select_states();origins=setup(a.raw_dir,step,para,items);write_csv(a.output/"selection_manifest.csv",[{k:v for k,v in x.items() if k not in ("state","dis_last")} for x in origins]);audit_path=a.output/"intervention_audit.json"
    if a.stage in ("audit","all"):result=audit(a.raw_dir,mod,step,para,origins,hashlib.sha256(SOURCE.read_bytes()).hexdigest(),adapter_sha);audit_path.write_text(json.dumps(result,indent=2)+"\n")
    else:result=json.loads(audit_path.read_text())
    if not result["parity"]["passed"]:raise RuntimeError("Baseline restart parity failed")
    print(json.dumps({"parity":True,"exact_cf_max":result["instantaneous_counterfactual"]["max_abs_exact_collective_after_symmetry"],"micro_residual":result["instantaneous_counterfactual"]["mean_abs_microscopic_residual"],"original":result["instantaneous_counterfactual"]["mean_abs_original"]},indent=2))
    if a.stage=="audit":return
    rows=simulate(origins,text,a.output,a.workers);metadata=dict(source_sha256=result["source_sha256"],adapter_sha256=adapter_sha,n_states=len(origins),clones=N_CLONES,horizon_steps=HORIZON_STEPS,horizon_time=HORIZON_STEPS*DT,paired_noise=True,rows=len(rows),selection="crossing-conditioned lead-1 states; archived continuation labels unused")
    (a.output/"run_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")


if __name__=="__main__":main()
