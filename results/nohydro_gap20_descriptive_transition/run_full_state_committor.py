#!/usr/bin/env python3
"""Bounded full-state forward-committor shooting experiment at gap 20."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import time
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta, spearmanr

import run_bottleneck_commitment_branches as bc
from analyze_approach_current_balance import RawReader


HERE=Path(__file__).resolve().parent
RAW_DEFAULT=Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
CACHE_Q=HERE/"outputs/initiation_propagation/cache"
OUT_DEFAULT=HERE/"outputs/full_state_committor"
LEVELS=(-0.4,-0.2,0.0,0.2)
STATES_PER_LEVEL=4
PILOT_SHOTS=8
EXPANDED_SHOTS=32
HORIZON=100.0
HORIZON_STEPS=round(HORIZON/bc.DT)
A_BOUND=-0.6
B_BOUND=0.6


def write_csv(path:Path,rows:list[dict])->None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def read_json(path:Path)->dict:
    return json.loads(path.read_text())


def wrap(x):
    return (np.asarray(x)+np.pi)%(2*np.pi)-np.pi


def inventory_and_select(cache_dir:Path)->tuple[list[dict],list[dict]]:
    inventory=[]; candidates={level:[] for level in LEVELS}
    for level in LEVELS:
        for seed in range(40):
            with np.load(cache_dir/f"seed{seed:02d}_stride100_qscan_v4.npz") as z:
                step=np.asarray(z["step"],int);q=np.asarray(z["q"],float)
            take=np.flatnonzero((step*bc.DT>=200)&(step*bc.DT<=900)&np.isclose(q,level,rtol=0,atol=1e-12))
            if len(take):
                candidates[level].append(dict(seed=seed,n_frames=len(take),native_step=int(step[take[0]]),time=float(step[take[0]]*bc.DT)))
        inventory.append(dict(q_level=level,eligible_seed_count=len(candidates[level]),
                              eligible_frame_count=sum(x["n_frames"] for x in candidates[level]),
                              time_window="200-900",exact_q_required=True))
    selected=[];used=set()
    for level in LEVELS:
        ordered=sorted(candidates[level],key=lambda x:(x["seed"] in used,x["seed"],x["native_step"]))
        chosen=[]
        for item in ordered:
            if item["seed"] in used:continue
            chosen.append(item);used.add(item["seed"])
            if len(chosen)==STATES_PER_LEVEL:break
        if len(chosen)<STATES_PER_LEVEL:
            raise RuntimeError(f"Only {len(chosen)} globally distinct seeds at q={level}")
        for rank,item in enumerate(chosen):
            selected.append(dict(state_id=f"q{level:+.1f}_s{item['seed']}_k{item['native_step']}",q_level=level,
                                 seed=item["seed"],native_step=item["native_step"],time=item["time"],selection_rank=rank,
                                 selection_rule="earliest exact-q frame after t=200 in first unused eligible seed; archived future unused"))
    return inventory,selected


def left_source_orientation(state:np.ndarray)->float:
    mask=state[0]<0
    com=state[:2,mask].mean(axis=1);gap_hat=-com/np.linalg.norm(com)
    heading=np.array([np.mean(np.cos(state[2,mask])),np.mean(np.sin(state[2,mask]))])
    return float(heading@gap_hat)


def load_origins(raw_dir:Path,selected:list[dict],step,para:dict)->list[dict]:
    zero=np.zeros((3,bc.N));out=[]
    for item in selected:
        reader=RawReader(bc.raw_path(raw_dir,item["seed"]));k=item["native_step"]
        previous=reader.block("state",k-1,k-1)[:,:,0].copy();state=reader.block("state",k,k)[:,:,0].copy()
        target=reader.block("state",k+1,k+1)[:,:,0].copy();reader.close()
        _,_,_,memory,_=step(previous,para,zero,None,False)
        q=float((2*np.sum(state[0]>=0)-bc.N)/bc.N)
        if not (A_BOUND<q<B_BOUND) or not np.isclose(q,item["q_level"],atol=1e-12):raise RuntimeError((item,q))
        out.append(dict(item,state=state,previous=previous,target=target,dis_last=memory,
                        nright=int(np.sum(state[0]>=0)),a_left=left_source_orientation(state)))
    return out


def transform_reflection(state:np.ndarray)->np.ndarray:
    out=state.copy();out[0]*=-1;out[2]=np.pi-out[2];return out


def audit(raw_dir:Path,mod,step,para:dict,origins:list[dict],adapter_sha:str)->dict:
    zero=np.zeros((3,bc.N));parity=[];reflection=[]
    for index,item in enumerate(origins):
        deterministic,_,_,_,_=step(item["state"],para,zero,item["dis_last"],False)
        increment=zero.copy();increment[2]=item["target"][2]-deterministic[2]
        predicted,velocity,_,_,flags=step(item["state"],para,increment,item["dis_last"],False)
        reader=RawReader(bc.raw_path(raw_dir,item["seed"]));archived_velocity=reader.block("rdot",item["native_step"]+1,item["native_step"]+1)[:,:,0];reader.close()
        parity.append(dict(state_id=item["state_id"],max_state_error=float(np.max(np.abs(predicted-item["target"]))),
                           max_velocity_error=float(np.max(np.abs(velocity-archived_velocity))),
                           inferred_noise_std=float(np.std(increment[2],ddof=1)),**flags))
        rng=np.random.default_rng(2026112300+index);noise=zero.copy();noise[2]=para["In"]*math.sqrt(bc.DT)*rng.normal(size=bc.N)
        original,v0,_,m0,_=step(item["state"],para,noise,item["dis_last"],False)
        reflected_state=transform_reflection(item["state"]);reflected_noise=zero.copy();reflected_noise[2]=-noise[2]
        reflected,v1,_,m1,_=step(reflected_state,para,reflected_noise,item["dis_last"],False)
        expected=transform_reflection(original);expected_velocity=v0.copy();expected_velocity[0]*=-1
        reflection.append(dict(state_id=item["state_id"],position_error=float(np.max(np.abs(reflected[:2]-expected[:2]))),
                               heading_vector_error=float(max(np.max(np.abs(np.cos(reflected[2])-np.cos(expected[2]))),
                                                              np.max(np.abs(np.sin(reflected[2])-np.sin(expected[2]))))),
                               velocity_error=float(np.max(np.abs(v1-expected_velocity))),memory_error=float(np.max(np.abs(m1-m0)))))
    max_parity=max(max(x["max_state_error"],x["max_velocity_error"]) for x in parity)
    max_reflection=max(max(x["position_error"],x["heading_vector_error"],x["velocity_error"],x["memory_error"]) for x in reflection)
    return dict(source_sha256=hashlib.sha256(bc.SOURCE.read_bytes()).hexdigest(),adapter_sha256=adapter_sha,
                native_replay_passed=max_parity<1e-12,max_native_replay_error=max_parity,parity_rows=parity,
                reflection_one_step_passed=max_reflection<1e-10,max_reflection_error=max_reflection,reflection_rows=reflection,
                previous_step_memory="reconstructed with validated production step from archived k-1 state; included in native k-to-k+1 replay",
                basin_A=f"q <= {A_BOUND}",basin_B=f"q >= {B_BOUND}",strict_interior_verified=True,
                first_hitting_check="basin tested after every native dt=0.001 step; simulation terminates immediately on first A or B hit",
                safety_horizon=HORIZON,censoring="unresolved at safety horizon; never classified as A")


def init_worker(source:str,adapter_text:str):
    bc.init_worker(source,adapter_text)


def run_shot(task:dict)->dict:
    state=task["state"].copy();memory=task["dis_last"].copy();rng=np.random.default_rng(task["noise_seed"])
    q0=(2*np.sum(state[0]>=0)-bc.N)/bc.N
    if not A_BOUND<q0<B_BOUND:raise RuntimeError("Shot did not start strictly between basins")
    outcome="censored";q_hit=float("nan");population_ok=True
    for n in range(HORIZON_STEPS):
        increment=np.zeros((3,bc.N));increment[2]=bc._PARA["In"]*math.sqrt(bc.DT)*rng.normal(size=bc.N)
        state,_,_,memory,_=bc._STEP(state,bc._PARA,increment,memory,False)
        nright=int(np.sum(state[0]>=0));population_ok &= 0<=nright<=bc.N
        q=(2*nright-bc.N)/bc.N
        if q<=A_BOUND:
            outcome="A";q_hit=float(q);break
        if q>=B_BOUND:
            outcome="B";q_hit=float(q);break
    steps=n+1
    return dict(state_id=task["state_id"],seed=task["seed"],q_level=task["q_level"],clone=task["clone"],
                noise_seed=task["noise_seed"],outcome=outcome,steps=steps,hit_time=steps*bc.DT,q_hit=q_hit,
                nright_end=int(np.sum(state[0]>=0)),population_conserved=int(population_ok),horizon=HORIZON)


def cache_file(output:Path,state_id:str,clone:int)->Path:
    return output/"cache/shots"/f"{state_id}_clone{clone:03d}.json"


def shoot(origins:list[dict],selected_ids:set[str],clones:range,adapter_text:str,output:Path,workers:int)->None:
    tasks=[]
    # Round-robin clone ordering obtains broad microscopic-state coverage
    # before adding repeated shots to any one state.
    for clone in clones:
        for item in origins:
            if item["state_id"] not in selected_ids:continue
            path=cache_file(output,item["state_id"],clone)
            if path.exists():continue
            tasks.append(dict(state_id=item["state_id"],seed=item["seed"],q_level=item["q_level"],clone=clone,
                              noise_seed=202611230000+item["seed"]*1000+clone,state=item["state"],dis_last=item["dis_last"]))
    if not tasks:return
    start=time.monotonic()
    with mp.get_context("spawn").Pool(min(workers,len(tasks)),initializer=init_worker,initargs=(str(bc.SOURCE),adapter_text)) as pool:
        for count,result in enumerate(pool.imap_unordered(run_shot,tasks),1):
            path=cache_file(output,result["state_id"],result["clone"]);path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(result,indent=2)+"\n")
            print(f"completed {count}/{len(tasks)} {result['state_id']} clone {result['clone']} {result['outcome']} t={result['hit_time']:.3f}; elapsed {time.monotonic()-start:.1f}s",flush=True)


def load_shots(output:Path)->list[dict]:
    return [read_json(path) for path in sorted((output/"cache/shots").glob("*.json"))]


def clopper(k:int,n:int,alpha:float=.05)->list[float]:
    if n==0:return [float("nan"),float("nan")]
    lo=0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1));hi=1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))
    return [lo,hi]


def state_summaries(origins:list[dict],shots:list[dict])->list[dict]:
    by=defaultdict(list)
    for row in shots:by[row["state_id"]].append(row)
    rows=[]
    for item in origins:
        part=by[item["state_id"]];a=sum(r["outcome"]=="A" for r in part);b=sum(r["outcome"]=="B" for r in part);c=sum(r["outcome"]=="censored" for r in part);resolved=a+b;n=len(part)
        rows.append(dict(state_id=item["state_id"],seed=item["seed"],q=item["q_level"],time=item["time"],a_left=item["a_left"],
                         n_shots=n,n_A=a,n_B=b,n_censored=c,n_resolved=resolved,resolved_fraction=resolved/n if n else 0,
                         committor_resolved=b/resolved if resolved else float("nan"),binomial_ci95=json.dumps(clopper(b,resolved)),
                         censoring_lower_bound=b/n if n else float("nan"),censoring_upper_bound=(b+c)/n if n else float("nan"),
                         mean_resolved_hit_time=float(np.mean([r["hit_time"] for r in part if r["outcome"]!="censored"])) if resolved else float("nan")))
    return rows


def variance_components(rows:list[dict],nboot:int=5000)->list[dict]:
    output=[]
    for level in LEVELS:
        part=[r for r in rows if np.isclose(r["q"],level) and r["n_resolved"]>=2]
        p=np.asarray([r["committor_resolved"] for r in part]);n=np.asarray([r["n_resolved"] for r in part],float)
        observed=float(np.var(p,ddof=1)) if len(p)>1 else float("nan")
        shooting=float(np.mean(p*(1-p)/(n-1))) if len(p) else float("nan")
        corrected=max(0.,observed-shooting) if np.isfinite(observed) else float("nan")
        rng=np.random.default_rng(202611231+int((level+1)*10));draw=[]
        if len(part)>1:
            for _ in range(nboot):
                idx=rng.integers(len(part),size=len(part));pp=p[idx];nn=n[idx]
                draw.append(max(0.,float(np.var(pp,ddof=1)-np.mean(pp*(1-pp)/(nn-1)))))
        output.append(dict(q=level,n_states=len(part),observed_state_variance=observed,estimated_shooting_variance=shooting,
                           corrected_committor_variance=corrected,corrected_committor_sd=math.sqrt(corrected) if np.isfinite(corrected) else float("nan"),
                           corrected_variance_bootstrap_ci95=json.dumps(np.percentile(draw,[2.5,97.5]).tolist()) if draw else json.dumps([float("nan")]*2),
                           method="sample variance of state estimates minus mean unbiased binomial proportion variance"))
    return output


def level_summaries(rows:list[dict])->list[dict]:
    """Equal-state summaries; censor bounds do not treat unresolved shots as A."""
    output=[]
    for level in LEVELS:
        part=[r for r in rows if np.isclose(r["q"],level)]
        output.append(dict(
            q=level,n_states=len(part),n_shots=sum(r["n_shots"] for r in part),
            n_A=sum(r["n_A"] for r in part),n_B=sum(r["n_B"] for r in part),
            n_censored=sum(r["n_censored"] for r in part),
            resolved_fraction=float(np.mean([r["resolved_fraction"] for r in part])),
            mean_state_lower_bound=float(np.mean([r["censoring_lower_bound"] for r in part])),
            mean_state_upper_bound=float(np.mean([r["censoring_upper_bound"] for r in part])),
            interpretation="bounds on the equal-state mean full committor at the safety horizon"))
    return output


def choose_expansion(rows:list[dict],variance:list[dict])->list[str]:
    eligible_levels=[]
    for row in variance:
        part=[r for r in rows if np.isclose(r["q"],row["q"]) and r["resolved_fraction"]>=.75]
        informative=sum(.1<r["committor_resolved"]<.9 for r in part)
        if len(part)>=3 and informative>=2:eligible_levels.append((row["corrected_committor_variance"],row["q"]))
    if not eligible_levels:return []
    level=max(eligible_levels)[1]
    return [r["state_id"] for r in rows if np.isclose(r["q"],level) and r["resolved_fraction"]>=.75]


def orientation_test(rows:list[dict],variance:list[dict],nperm:int=10000)->dict:
    # The additional coordinate is tested only after fixed-q heterogeneity is
    # resolved across at least three adequately sampled microscopic states.
    eligible=[]
    for v in variance:
        ci=json.loads(v["corrected_variance_bootstrap_ci95"])
        part=[r for r in rows if np.isclose(r["q"],v["q"]) and r["resolved_fraction"]>=.75]
        if len(part)>=3 and np.isfinite(ci[0]) and ci[0]>.01:eligible.append(v["q"])
    gate=bool(eligible)
    valid=[r for r in rows if r["n_resolved"]>=2 and np.isfinite(r["committor_resolved"])]
    if not gate:return dict(tested=False,reason="fixed-q spread was not resolved: no level had >=3 states with >=75% resolved shots and a bootstrap lower bound above 0.01")
    levels=sorted({r["q"] for r in valid});x=[];y=[]
    for level in levels:
        part=[r for r in valid if np.isclose(r["q"],level)];mx=np.mean([r["a_left"] for r in part]);my=np.mean([r["committor_resolved"] for r in part])
        x.extend(r["a_left"]-mx for r in part);y.extend(r["committor_resolved"]-my for r in part)
    rho=float(spearmanr(x,y).statistic);rng=np.random.default_rng(202611232);null=[]
    for _ in range(nperm):
        shuffled=[]
        for level in levels:
            vals=[r["a_left"] for r in valid if np.isclose(r["q"],level)];shuffled.extend(rng.permutation(vals)-np.mean(vals))
        null.append(float(spearmanr(shuffled,y).statistic))
    p=(1+np.sum(np.abs(null)>=abs(rho)))/(nperm+1)
    return dict(tested=True,variable="left/source-school mean heading dotted with its COM-to-gap direction",
                statistic="Spearman correlation after centering both variables within exact q level",rho=rho,permutation_p_two_sided=float(p),
                n_states=len(valid),interpretation="association test only; orientation is not fitted as a committor model")


def reflection_statistics(rows:list[dict],nboot:int=5000)->dict:
    neg=np.asarray([r["committor_resolved"] for r in rows if np.isclose(r["q"],-.2) and r["n_resolved"]>0]);pos=np.asarray([r["committor_resolved"] for r in rows if np.isclose(r["q"],.2) and r["n_resolved"]>0])
    point=float(np.mean(neg)+np.mean(pos)-1);rng=np.random.default_rng(202611233)
    draw=rng.choice(neg,(nboot,len(neg)),replace=True).mean(axis=1)+rng.choice(pos,(nboot,len(pos)),replace=True).mean(axis=1)-1
    return dict(definition="mean h(q=-0.2)+mean h(q=+0.2)-1; zero under reflection symmetry",difference=point,ci95=np.percentile(draw,[2.5,97.5]).tolist(),n_negative=len(neg),n_positive=len(pos))


def make_figure(rows:list[dict],variance:list[dict],orientation:dict,path:Path)->None:
    fig,axes=plt.subplots(1,3,figsize=(13,4));colors=plt.cm.viridis(np.linspace(.1,.9,len(LEVELS)))
    for color,level in zip(colors,LEVELS):
        part=sorted([r for r in rows if np.isclose(r["q"],level)],key=lambda r:r["a_left"])
        for i,r in enumerate(part):
            axes[0].plot([level-.012,level+.012],[r["censoring_lower_bound"]]*2,color=color,alpha=0)
            axes[0].errorbar(level+(i-1.5)*.012,r["committor_resolved"],yerr=None,fmt="o",color=color)
            axes[0].vlines(level+(i-1.5)*.012,r["censoring_lower_bound"],r["censoring_upper_bound"],color=color,lw=2)
    axes[0].plot([-.6,.6],[0,1],"k:",alpha=.4);axes[0].set(xlabel="q",ylabel="resolved-shot committor",title="Full-state estimates\n(lines: censoring bounds)",ylim=(-.05,1.05))
    axes[1].bar([str(r["q"]) for r in variance],[r["corrected_committor_sd"] for r in variance],color=colors)
    axes[1].set(xlabel="exact q",ylabel="corrected between-state SD",title="Fixed-q spread estimator\n(unresolved under censoring)")
    if orientation["tested"]:
        for color,level in zip(colors,LEVELS):
            part=[r for r in rows if np.isclose(r["q"],level)];axes[2].scatter([r["a_left"] for r in part],[r["committor_resolved"] for r in part],color=color,label=f"q={level:g}")
        axes[2].set(xlabel="source-school gap orientation",ylabel="resolved-shot committor",title="Predeclared orientation diagnostic",ylim=(-.05,1.05));axes[2].legend(frameon=False,fontsize=8)
    else:
        axes[2].axis("off");axes[2].text(.5,.55,"Orientation test not opened",ha="center",va="center",fontsize=15,transform=axes[2].transAxes)
        axes[2].text(.5,.42,"Fixed-q committor spread was not resolved",ha="center",va="center",fontsize=11,transform=axes[2].transAxes)
    fig.tight_layout();fig.savefig(path,dpi=220);plt.close(fig)


def analyze(origins:list[dict],output:Path)->dict:
    shots=load_shots(output);rows=state_summaries(origins,shots);levels=level_summaries(rows);variance=variance_components(rows);orientation=orientation_test(rows,variance);reflection=reflection_statistics(rows)
    write_csv(output/"shots.csv",shots);write_csv(output/"state_committors.csv",rows);write_csv(output/"q_level_summary.csv",levels);write_csv(output/"fixed_q_variance.csv",variance)
    (output/"orientation_test.json").write_text(json.dumps(orientation,indent=2)+"\n");(output/"reflection_statistics.json").write_text(json.dumps(reflection,indent=2)+"\n")
    make_figure(rows,variance,orientation,output/"full_state_committor.png")
    summary=dict(n_states=len(rows),n_shots=len(shots),outcomes=dict(A=sum(r["outcome"]=="A" for r in shots),B=sum(r["outcome"]=="B" for r in shots),censored=sum(r["outcome"]=="censored" for r in shots)),
                 q_level_summary=levels,state_committors=rows,fixed_q_variance=variance,orientation_test=orientation,reflection_statistics=reflection)
    (output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n");return summary


def self_test()->None:
    assert A_BOUND<-.4<B_BOUND and A_BOUND<.2<B_BOUND
    assert clopper(0,8)[0]==0 and clopper(8,8)[1]==1
    fake=[dict(q=q,n_resolved=8,committor_resolved=p) for q,p in [(-.2,.2),(-.2,.4),(.2,.6),(.2,.8)]]
    assert abs(reflection_statistics(fake,100)["difference"])<1e-12
    print("PASS: basin, exact-binomial interval, and reflection-statistic tests")


def main()->None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path,default=RAW_DEFAULT);p.add_argument("--cache-q",type=Path,default=CACHE_Q);p.add_argument("--output",type=Path,default=OUT_DEFAULT)
    p.add_argument("--workers",type=int,default=9);p.add_argument("--stage",choices=("audit","pilot","expand","analyze","all"),default="all");p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    args.output.mkdir(parents=True,exist_ok=True);inventory,selected=inventory_and_select(args.cache_q);write_csv(args.output/"coverage_inventory.csv",inventory);write_csv(args.output/"selected_states.csv",selected)
    mod=bc.load_gap20_production(bc.SOURCE);para=bc.parameters(mod);step,adapter_text,adapter_sha=bc.extract_step(mod,bc.SOURCE);origins=load_origins(args.raw_dir,selected,step,para)
    audit_result=audit(args.raw_dir,mod,step,para,origins,adapter_sha);(args.output/"restart_audit.json").write_text(json.dumps(audit_result,indent=2)+"\n")
    if not audit_result["native_replay_passed"] or not audit_result["reflection_one_step_passed"]:raise RuntimeError("Restart/reflection audit failed")
    write_csv(args.output/"state_metadata.csv",[{k:v for k,v in r.items() if k not in ("state","previous","target","dis_last")} for r in origins])
    if args.stage=="audit":return
    all_ids={r["state_id"] for r in origins}
    if args.stage in ("pilot","all"):shoot(origins,all_ids,range(PILOT_SHOTS),adapter_text,args.output,args.workers)
    pilot_rows=state_summaries(origins,load_shots(args.output));pilot_variance=variance_components(pilot_rows);expand_ids=set(choose_expansion(pilot_rows,pilot_variance))
    (args.output/"expansion_selection.json").write_text(json.dumps(dict(selected_state_ids=sorted(expand_ids),rule="one exact-q level with largest corrected pilot variance, requiring >=3 states with >=75% resolution and >=2 intermediate estimates; selected states expanded from 8 to 32 shots"),indent=2)+"\n")
    if args.stage in ("expand","all") and expand_ids:shoot(origins,expand_ids,range(PILOT_SHOTS,EXPANDED_SHOTS),adapter_text,args.output,args.workers)
    analyze(origins,args.output)


if __name__=="__main__":main()
