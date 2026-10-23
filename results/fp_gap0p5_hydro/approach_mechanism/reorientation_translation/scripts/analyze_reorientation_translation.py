#!/usr/bin/env python3
"""Validate collective reorientation followed by inward translation.

No potential is fitted and no simulation is run.  The script uses the archived
coordinate table, compact coordinate caches, and saved state/rdot trajectories.
"""
from __future__ import annotations
import argparse, csv, json, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
COORD_SCRIPTS=HERE.parent.parent/"collective_coordinates"/"scripts"
sys.path.insert(0,str(COORD_SCRIPTS))
from analyze_collective_coordinates import (  # noqa: E402
    ALL_SEEDS, TRAIN_SEEDS, TEST_SEEDS, auc_score, fit_logistic,
    predict_logistic, read_csv, write_csv,
)

LEADS=(-2.0,-1.0)
DT=0.001
FEATURE_DT=0.01
IF=1e-2
DELTA=1e-2


def bootstrap_mean(v,seed=20261023):
    v=np.asarray(v,float); rng=np.random.default_rng(seed)
    b=np.mean(v[rng.integers(0,len(v),(30000,len(v)))],axis=1)
    return float(np.percentile(b,2.5)),float(np.percentile(b,97.5))


def add_sigma(rows):
    for r in rows:
        r["sigma_perp"]=float(r["Lambda"])*float(r["distance"])


def geometric_decomposition(samples):
    lookup={(int(r["seed"]),int(r["pair_id"]),r["condition"],float(r["lead"])):r for r in samples}
    pair_rows=[]; seed_rows=[]; summary=[]
    for lead in LEADS:
      for seed in ALL_SEEDS:
        ids=sorted({int(r["pair_id"]) for r in samples if int(r["seed"])==seed and float(r["lead"])==lead})
        local=[]
        for pair_id in ids:
          e=lookup[(seed,pair_id,"event",lead)]; c=lookup[(seed,pair_id,"control",lead)]
          dl=float(np.log(float(e["Lambda"]))-np.log(float(c["Lambda"])))
          ds=float(np.log(float(e["sigma_perp"]))-np.log(float(c["sigma_perp"])))
          nd=float(-(np.log(float(e["distance"]))-np.log(float(c["distance"]))))
          row=dict(seed=seed,pair_id=pair_id,lead=lead,delta_log_Lambda=dl,
                   delta_log_sigma_perp=ds,minus_delta_log_distance=nd,
                   identity_residual=dl-ds-nd)
          pair_rows.append(row); local.append([dl,ds,nd])
        m=np.mean(local,axis=0)
        seed_rows.append(dict(seed=seed,lead=lead,delta_log_Lambda=m[0],
            delta_log_sigma_perp=m[1],minus_delta_log_distance=m[2]))
      q=[r for r in seed_rows if r["lead"]==lead]
      means={k:float(np.mean([r[k] for r in q])) for k in ("delta_log_Lambda","delta_log_sigma_perp","minus_delta_log_distance")}
      for key,val in means.items():
        lo,hi=bootstrap_mean([r[key] for r in q])
        summary.append(dict(lead=lead,component=key,mean=val,seed_bootstrap_ci_low=lo,
                            seed_bootstrap_ci_high=hi,independent_seeds=len(q)))
      summary.append(dict(lead=lead,component="distance_fraction_of_mean_delta_log_Lambda",
          mean=means["minus_delta_log_distance"]/means["delta_log_Lambda"],
          seed_bootstrap_ci_low="",seed_bootstrap_ci_high="",independent_seeds=len(q)))
    return pair_rows,seed_rows,summary


MODEL_SPECS={
    "psi":("cos_psi","sin_psi"),
    "psi_Lambda":("cos_psi","sin_psi","Lambda"),
    "psi_distance":("cos_psi","sin_psi","distance"),
    "psi_sigma_perp":("cos_psi","sin_psi","sigma_perp"),
    "psi_distance_sigma_perp":("cos_psi","sin_psi","distance","sigma_perp"),
}


def logistic_comparison(samples):
    pooled=[]; byseed=[]; models={}
    for lead in LEADS:
      train=[r for r in samples if float(r["lead"])==lead and int(r["seed"]) in TRAIN_SEEDS]
      test=[r for r in samples if float(r["lead"])==lead and int(r["seed"]) in TEST_SEEDS]
      for name,keys in MODEL_SPECS.items():
        xtr=np.array([[float(r[k]) for k in keys] for r in train]); ytr=np.array([r["condition"]=="event" for r in train])
        xte=np.array([[float(r[k]) for k in keys] for r in test]); yte=np.array([r["condition"]=="event" for r in test])
        model=fit_logistic(xtr,ytr); models[(lead,name)]=model
        prob=predict_logistic(model,xte)
        pooled.append(dict(lead=lead,model=name,features=";".join(keys),n_train=len(train),n_test=len(test),
          auc=auc_score(yte,prob),brier=float(np.mean((prob-yte)**2)),
          nll=float(-np.mean(yte*np.log(prob+1e-12)+(1-yte)*np.log(1-prob+1e-12)))))
        for seed in TEST_SEEDS:
          q=[r for r in test if int(r["seed"])==seed]
          x=np.array([[float(r[k]) for k in keys] for r in q]); y=np.array([r["condition"]=="event" for r in q])
          p=predict_logistic(model,x)
          byseed.append(dict(lead=lead,model=name,seed=seed,n=len(q),auc=auc_score(y,p),
            brier=float(np.mean((p-y)**2)),nll=float(-np.mean(y*np.log(p+1e-12)+(1-y)*np.log(1-p+1e-12)))))
    contrasts=[]
    for lead in LEADS:
      q=[r for r in byseed if r["lead"]==lead]
      comparisons=[(name,"psi") for name in MODEL_SPECS if name!="psi"]
      comparisons += [("psi_distance_sigma_perp","psi_distance"),("psi_Lambda","psi_distance")]
      for name,baseline in comparisons:
        base={int(r["seed"]):float(r["auc"]) for r in q if r["model"]==baseline}
        other={int(r["seed"]):float(r["auc"]) for r in q if r["model"]==name}
        diff=np.array([other[s]-base[s] for s in TEST_SEEDS]); lo,hi=bootstrap_mean(diff)
        contrasts.append(dict(lead=lead,comparison=f"{name}_minus_{baseline}",mean_seed_auc_difference=float(np.mean(diff)),
          seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,heldout_seeds=len(diff),
          seed_differences=";".join(f"{x:.8g}" for x in diff)))
    return pooled,byseed,contrasts


def load_cache(cache_dir,seed):
    p=cache_dir/f"collective_coordinates_velocity_v2_seed{seed:02d}.npz"
    with np.load(p) as z: return {k:z[k] for k in z.files}


def prospective_distance_matched_risk_set(event_pairs,cache_dir,horizon=2.0):
    """Match A-2 event origins to controls crossing-free only over [t0,t0+2]."""
    rows=[]; h=int(round(horizon/FEATURE_DT))
    for seed in ALL_SEEDS:
      z=load_cache(cache_dir,seed)
      events=[r for r in event_pairs if int(r["seed"])==seed and r["condition"]=="event"]
      # Each pair list already contains the thinned, qualified first crossings.
      candidates={1:[], -1:[]}
      for i in range(int(200/FEATURE_DT),len(z["time"])-h,10):
        side=int(z["side"][i]) if np.isfinite(z["side"][i]) else 0
        if side in candidates and np.all(z["side"][i:i+h+1]==side): candidates[side].append(i)
      available={s:set(v) for s,v in candidates.items()}
      # Standardize matching coordinates within seed, preserving circular phase.
      allidx=np.array(candidates[1]+candidates[-1],int)
      mat=np.column_stack((np.sin(z["phase"][allidx]),np.cos(z["phase"][allidx]),z["distance"][allidx]))
      scale=np.std(mat,axis=0); scale[scale<1e-9]=1
      for event in events:
        side=1 if event["origin"]=="left" else -1
        i=int(round((float(event["anchor_time"])-2)/FEATURE_DT))
        cand=np.array(sorted(available[side]),int)
        if not len(cand): continue
        target=np.array([np.sin(z["phase"][i]),np.cos(z["phase"][i]),z["distance"][i]])
        cm=np.column_stack((np.sin(z["phase"][cand]),np.cos(z["phase"][cand]),z["distance"][cand]))
        cost=np.linalg.norm((cm-target)/scale,axis=1); j=int(cand[np.argmin(cost)])
        available[side]={k for k in available[side] if abs(k-j)>50}
        for condition,idx,md in (("event",i,0.0),("control",j,float(np.min(cost)))):
          rows.append(dict(seed=seed,pair_id=event["pair_id"],condition=condition,origin=event["origin"],
            risk_time=float(z["time"][idx]),future_horizon=horizon,match_distance=md,
            phase=float(z["phase"][idx]),distance=float(z["distance"][idx]),psi=float(z["psi"][idx]),
            cos_psi=float(np.cos(z["psi"][idx])),sin_psi=float(np.sin(z["psi"][idx])),
            outcome_first_crossing_within_horizon=int(condition=="event")))
    return rows


def risk_set_tests(rows):
    train=[r for r in rows if int(r["seed"]) in TRAIN_SEEDS]
    test=[r for r in rows if int(r["seed"]) in TEST_SEEDS]
    keys=("cos_psi","sin_psi")
    model=fit_logistic(np.array([[r[k] for k in keys] for r in train]),np.array([r["condition"]=="event" for r in train]))
    out=[]
    for label,q in [("pooled",test)]+[(f"seed_{s}",[r for r in test if int(r["seed"])==s]) for s in TEST_SEEDS]:
      x=np.array([[r[k] for k in keys] for r in q]); y=np.array([r["condition"]=="event" for r in q]); p=predict_logistic(model,x)
      out.append(dict(sample=label,n=len(q),auc=auc_score(y,p),brier=float(np.mean((p-y)**2)),
        nll=float(-np.mean(y*np.log(p+1e-12)+(1-y)*np.log(1-p+1e-12)))))
    paired=[]
    for metric in ("cos_psi","distance"):
      seedvals=[]
      for seed in ALL_SEEDS:
        q=[r for r in rows if int(r["seed"])==seed]; by={(r["pair_id"],r["condition"]):r for r in q}
        ids={r["pair_id"] for r in q}; diff=[float(by[(i,"event")][metric])-float(by[(i,"control")][metric]) for i in ids if (i,"event") in by and (i,"control") in by]
        if diff: seedvals.append(float(np.mean(diff)))
      lo,hi=bootstrap_mean(seedvals)
      paired.append(dict(metric=metric,event_minus_control=float(np.mean(seedvals)),seed_bootstrap_ci_low=lo,
                         seed_bootstrap_ci_high=hi,independent_seeds=len(seedvals)))
    aucs=np.array([r["auc"] for r in out if r["sample"].startswith("seed_")]); lo,hi=bootstrap_mean(aucs)
    return out,paired,dict(mean_seed_auc=float(np.mean(aucs)),seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,heldout_seeds=len(aucs))


def direct_velocity_x_batch(state,batch=64):
    """Direct dipolar x velocity for (3,N,T) states, batched over T."""
    ntime=state.shape[2]; out=np.empty((state.shape[1],ntime)); eye=np.eye(state.shape[1],dtype=bool)
    for a in range(0,ntime,batch):
      b=min(ntime,a+batch)
      z=(state[0,:,a:b]+1j*state[1,:,a:b]).T  # B,N
      th=state[2,:,a:b].T
      dz=z[:,:,None]-z[:,None,:]
      den=dz*dz+DELTA*DELTA
      term=np.exp(1j*th)[:,None,:]/den
      term[:,eye]=0
      u=np.sum(term,axis=2)*IF/np.pi
      out[:,a:b]=u.real.T
    return out


def integrate_segments(event_pairs,raw_dir):
    rows=[]; audit=[]
    for seed in ALL_SEEDS:
      print(f"kinematics seed {seed}",flush=True)
      pairs=[r for r in event_pairs if int(r["seed"])==seed]
      raw=raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
      import scipy.io as sio
      data=sio.loadmat(raw,variable_names=["state","rdot"],squeeze_me=True)
      state=np.asarray(data["state"]); rdot=np.asarray(data["rdot"])
      intervals=[]; union=[]
      for r in pairs:
        i0=int(round((float(r["anchor_time"])-2)/DT)); i1=int(round((float(r["anchor_time"])-1)/DT))
        origins=np.arange(i0,i1,dtype=int); intervals.append((r,i0,i1,origins)); union.extend(origins.tolist())
      unique=np.array(sorted(set(union)),int); pos={k:i for i,k in enumerate(unique)}
      st=state[:,:,unique]; dest=rdot[:,:,unique+1]
      direct=direct_velocity_x_batch(st)
      ex=np.cos(st[2]); ey=np.sin(st[2])
      frozen=np.hypot(dest[0],dest[1])<1e-13
      self_only=(np.hypot(dest[0]-ex,dest[1]-ey)<1e-11) & ~frozen
      normal=~frozen & ~self_only
      direct_eff=np.where(normal,direct,0.0)
      self_eff=np.where(frozen,0.0,ex)
      boundary_eff=dest[0]-self_eff-direct_eff
      max_identity=float(np.max(np.abs(dest[0]-self_eff-direct_eff-boundary_eff)))
      audit.append(dict(seed=seed,unique_native_updates=len(unique),frozen_fraction=float(np.mean(frozen)),
        self_only_override_fraction=float(np.mean(self_only)),max_velocity_identity_error=max_identity))
      for r,i0,i1,origins in intervals:
        jj=np.array([pos[k] for k in origins]); side=1.0 if r["origin"]=="left" else -1.0
        cself=-side*float(np.sum(np.mean(self_eff[:,jj],axis=0))*DT)
        cdirect=-side*float(np.sum(np.mean(direct_eff[:,jj],axis=0))*DT)
        cboundary=-side*float(np.sum(np.mean(boundary_eff[:,jj],axis=0))*DT)
        d0=-side*float(np.mean(state[0,:,i0])); d1=-side*float(np.mean(state[0,:,i1]))
        s0=float(np.std(state[0,:,i0])); s1=float(np.std(state[0,:,i1]))
        actual=d1-d0; residual=actual-cself-cdirect-cboundary
        rows.append(dict(seed=seed,pair_id=r["pair_id"],condition=r["condition"],origin=r["origin"],
          start_time=float(r["anchor_time"])-2,end_time=float(r["anchor_time"])-1,
          distance_start=d0,distance_end=d1,actual_delta_distance=actual,
          self_propulsion_delta_distance=cself,direct_hydro_delta_distance=cdirect,
          boundary_hydro_effective_delta_distance=cboundary,integration_residual=residual,
          sigma_perp_start=s0,sigma_perp_end=s1,actual_delta_sigma_perp=s1-s0,
          frozen_update_fraction=float(np.mean(frozen[:,jj])),self_only_override_fraction=float(np.mean(self_only[:,jj]))))
      del state,rdot,data,st,dest,direct
    return rows,audit


def segment_summaries(rows):
    metrics=("actual_delta_distance","self_propulsion_delta_distance","direct_hydro_delta_distance",
             "boundary_hydro_effective_delta_distance","integration_residual","actual_delta_sigma_perp")
    group=[]; contrast=[]
    for cond in ("event","control"):
      for key in metrics:
        seedmeans=[]
        for seed in ALL_SEEDS:
          q=[float(r[key]) for r in rows if r["condition"]==cond and int(r["seed"])==seed]
          if q: seedmeans.append(np.mean(q))
        lo,hi=bootstrap_mean(seedmeans)
        group.append(dict(condition=cond,metric=key,mean_of_seed_means=float(np.mean(seedmeans)),
          seed_bootstrap_ci_low=lo,seed_bootstrap_ci_high=hi,independent_seeds=len(seedmeans)))
    for key in metrics:
      vals=[]
      for seed in ALL_SEEDS:
        e=[float(r[key]) for r in rows if r["condition"]=="event" and int(r["seed"])==seed]
        c=[float(r[key]) for r in rows if r["condition"]=="control" and int(r["seed"])==seed]
        if e and c: vals.append(np.mean(e)-np.mean(c))
      lo,hi=bootstrap_mean(vals)
      contrast.append(dict(metric=key,event_minus_control=float(np.mean(vals)),seed_bootstrap_ci_low=lo,
                           seed_bootstrap_ci_high=hi,independent_seeds=len(vals)))
    return group,contrast


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--analysis-dir",type=Path,required=True)
    ap.add_argument("--raw-dir",type=Path,required=True); ap.add_argument("--cache-dir",type=Path,required=True)
    ap.add_argument("--skip-kinematics",action="store_true")
    args=ap.parse_args(); out=args.analysis_dir.resolve(); (out/"tables").mkdir(parents=True,exist_ok=True)
    prior=out.parent/"collective_coordinates"/"tables"
    samples=read_csv(prior/"coordinate_samples.csv"); add_sigma(samples)
    pair_geom,seed_geom,geom_summary=geometric_decomposition(samples)
    pooled,byseed,auc_contrasts=logistic_comparison(samples)
    event_pairs=read_csv(prior/"matched_event_controls.csv")
    risk=prospective_distance_matched_risk_set(event_pairs,args.cache_dir)
    risk_metrics,risk_paired,risk_seed_summary=risk_set_tests(risk)
    if args.skip_kinematics:
      kin=read_csv(out/"tables"/"segment_kinematics.csv"); audit=read_csv(out/"tables"/"kinematic_identity_audit.csv")
    else:
      kin,audit=integrate_segments(event_pairs,args.raw_dir)
    kin_group,kin_contrast=segment_summaries(kin)
    for name,rows in (("geometric_pair_decomposition.csv",pair_geom),("geometric_seed_decomposition.csv",seed_geom),
      ("geometric_decomposition_summary.csv",geom_summary),("coordinate_model_pooled.csv",pooled),
      ("coordinate_model_by_seed.csv",byseed),("coordinate_model_auc_contrasts.csv",auc_contrasts),
      ("distance_matched_risk_set.csv",risk),("distance_matched_risk_metrics.csv",risk_metrics),
      ("distance_matched_paired_contrasts.csv",risk_paired),("segment_kinematics.csv",kin),
      ("kinematic_identity_audit.csv",audit),("segment_kinematics_group.csv",kin_group),
      ("segment_kinematics_contrasts.csv",kin_contrast)):
        write_csv(out/"tables"/name,rows)
    write_csv(out/"tables"/"distance_matched_seed_summary.csv",[risk_seed_summary])
    config=dict(parent_commit="2496cec",seeds=list(ALL_SEEDS),training_seeds=list(TRAIN_SEEDS),heldout_seeds=list(TEST_SEEDS),
      leads=list(LEADS),risk_horizon=2.0,n_event_control_pairs=len(event_pairs)//2,
      kinematics="exact native saved update accounting over A-2 to A-1; direct hydro reconstructed from source equation; effective BEM is saved velocity remainder after explicit simulator overrides",
      no_new_simulations=True,no_potential_fit=True)
    (out/"configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    print(json.dumps(dict(samples=len(samples),risk_rows=len(risk),kinematic_segments=len(kin)),indent=2))

if __name__=="__main__": main()
