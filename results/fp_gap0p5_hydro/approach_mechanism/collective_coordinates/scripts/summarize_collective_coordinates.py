#!/usr/bin/env python3
"""Robust summaries and seed-level validation for collective-coordinate test."""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_collective_coordinates import (ALL_SEEDS, LEADS, MODELS, TRAIN_SEEDS,
    TEST_SEEDS, auc_score, fit_logistic, predict_logistic, read_csv, write_csv)


ANGULAR=("wall","alignment","attraction","direct_hydro","boundary_hydro")


def bootstrap_ci(values,seed=20261022):
    v=np.asarray(values,float); rng=np.random.default_rng(seed)
    b=np.mean(v[rng.integers(0,len(v),(20000,len(v)))],axis=1)
    return np.percentile(b,[2.5,97.5])


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("analysis_dir",type=Path); args=ap.parse_args()
    out=args.analysis_dir.resolve(); tab=out/"tables"
    samples=read_csv(tab/"coordinate_samples.csv")
    replay=read_csv(tab/"deterministic_coordinate_dynamics.csv")
    lookup={(r["seed"],r["pair_id"],r["condition"],r["lead"]):r for r in samples}
    for r in replay:
        q=lookup[(r["seed"],r["pair_id"],r["condition"],r["lead"])]
        r["phase"]=q["phase"]; r["polarization"]=q["polarization"]
        angular=sum(float(r[f"psi_{k}"]) for k in ANGULAR)
        r["psi_angular_total"]=angular
        r["opening_closing_angular_total"]=-np.sign(float(r["psi"]))*angular
        r["wall_following_closing_angular_total"]=-np.sign(float(r["wall_following_error"]))*angular
    write_csv(tab/"deterministic_coordinate_dynamics.csv",replay)

    metrics=[f"psi_{k}" for k in ANGULAR]+["psi_angular_total","psi_deterministic_total",
        "opening_closing_angular_total","wall_following_closing_angular_total",
        "opening_closing_wall","wall_following_closing_wall",
        "Lambda_self","Lambda_direct_hydro","Lambda_boundary_hydro","Lambda_deterministic_total"]
    robust=[]
    for cond in ("event","control"):
      for lead in LEADS:
       q=[r for r in replay if r["condition"]==cond and float(r["lead"])==lead]
       for key in metrics:
        v=np.array([float(r[key]) for r in q]); robust.append(dict(condition=cond,lead=lead,metric=key,n=len(v),
            mean=float(np.mean(v)),median=float(np.median(v)),q25=float(np.percentile(v,25)),q75=float(np.percentile(v,75)),minimum=float(np.min(v)),maximum=float(np.max(v))))
    write_csv(tab/"robust_dynamics_summary.csv",robust)

    contrasts=[]
    for lead in LEADS:
      for key in metrics:
       diff=[]
       for seed in ALL_SEEDS:
        e=[r for r in replay if int(r["seed"])==seed and r["condition"]=="event" and float(r["lead"])==lead]
        c=[r for r in replay if int(r["seed"])==seed and r["condition"]=="control" and float(r["lead"])==lead]
        if e and c: diff.append(float(e[0][key])-float(c[0][key]))
       lo,hi=bootstrap_ci(diff)
       contrasts.append(dict(lead=lead,metric=key,event_minus_control_mean=float(np.mean(diff)),
          event_minus_control_median=float(np.median(diff)),seed_bootstrap_ci_low=float(lo),seed_bootstrap_ci_high=float(hi),independent_seeds=len(diff)))
    write_csv(tab/"paired_dynamics_contrasts.csv",contrasts)

    dependence=[]
    for cond in ("event","control"):
      for lead in LEADS:
       q=[r for r in replay if r["condition"]==cond and float(r["lead"])==lead]
       lam=np.array([float(r["Lambda"]) for r in q]); phase=np.array([float(r["phase"]) for r in q])
       x=np.column_stack(((lam-lam.mean())/lam.std(),np.sin(phase),np.cos(phase)))
       X=np.column_stack((np.ones(len(x)),x))
       for key in ("opening_closing_angular_total","wall_following_closing_angular_total","opening_closing_wall","Lambda_deterministic_total"):
        y=np.array([float(r[key]) for r in q]); b=np.linalg.lstsq(X,y,rcond=None)[0]; pred=X@b
        r2=1-float(np.sum((y-pred)**2)/np.sum((y-y.mean())**2)) if np.std(y)>0 else 0
        dependence.append(dict(condition=cond,lead=lead,metric=key,n=len(y),r2_Lambda_and_phase=r2,
            coefficient_standardized_Lambda=float(b[1]),coefficient_sin_phase=float(b[2]),coefficient_cos_phase=float(b[3])))
    write_csv(tab/"dynamics_dependence.csv",dependence)

    byseed=[]; aggregate=[]
    for lead in LEADS:
      train=[r for r in samples if float(r["lead"])==lead and int(r["seed"]) in TRAIN_SEEDS]
      for name,keys in MODELS.items():
       Xtr=np.array([[float(r[k]) for k in keys] for r in train]); ytr=np.array([r["condition"]=="event" for r in train])
       model=fit_logistic(Xtr,ytr); aucs=[]
       for seed in TEST_SEEDS:
        te=[r for r in samples if float(r["lead"])==lead and int(r["seed"])==seed]
        X=np.array([[float(r[k]) for k in keys] for r in te]); y=np.array([r["condition"]=="event" for r in te])
        prob=predict_logistic(model,X); auc=auc_score(y,prob); aucs.append(auc)
        byseed.append(dict(lead=lead,model=name,seed=seed,n=len(te),auc=auc,
            brier=float(np.mean((prob-y)**2)),nll=float(-np.mean(y*np.log(prob+1e-12)+(1-y)*np.log(1-prob+1e-12)))))
       lo,hi=bootstrap_ci(aucs)
       aggregate.append(dict(lead=lead,model=name,mean_seed_auc=float(np.mean(aucs)),
          seed_bootstrap_ci_low=float(lo),seed_bootstrap_ci_high=float(hi),heldout_seeds=len(aucs)))
    write_csv(tab/"heldout_prediction_by_seed.csv",byseed); write_csv(tab/"heldout_prediction_seed_summary.csv",aggregate)

    ablations=[]
    for lead in LEADS:
      train=[r for r in samples if float(r["lead"])==lead and int(r["seed"]) in TRAIN_SEEDS]
      test=[r for r in samples if float(r["lead"])==lead and int(r["seed"]) in TEST_SEEDS]
      for name,keys in (("psi_only",("cos_psi","sin_psi")),("Lambda_only",("Lambda",))):
       Xtr=np.array([[float(r[k]) for k in keys] for r in train]); ytr=np.array([r["condition"]=="event" for r in train])
       Xte=np.array([[float(r[k]) for k in keys] for r in test]); yte=np.array([r["condition"]=="event" for r in test])
       prob=predict_logistic(fit_logistic(Xtr,ytr),Xte)
       ablations.append(dict(lead=lead,model=name,features=";".join(keys),n_train=len(train),n_test=len(test),
           auc=auc_score(yte,prob),brier=float(np.mean((prob-yte)**2)),
           nll=float(-np.mean(yte*np.log(prob+1e-12)+(1-yte)*np.log(1-prob+1e-12)))))
    write_csv(tab/"heldout_coordinate_ablation.csv",ablations)

    phase_match=np.array([float(r["phase_match_distance"]) for r in read_csv(tab/"matched_event_controls.csv") if r["condition"]=="control"])
    pred=read_csv(tab/"heldout_prediction.csv"); paired=read_csv(tab/"paired_coordinate_statistics.csv"); memory=read_csv(tab/"memory_validation.csv")
    summary=[]
    def add(metric,value,detail=""): summary.append(dict(metric=metric,value=value,detail=detail))
    add("n_nucleation_control_pairs",len(phase_match),"20 independent seeds")
    add("phase_match_median_chord_distance",float(np.median(phase_match)))
    add("phase_match_p99_chord_distance",float(np.percentile(phase_match,99)))
    for r in pred: add(f"heldout_auc_{r['model']}_lead{abs(float(r['lead'])):g}",r["auc"],"seeds 15-19; trained on 0-14")
    for r in paired:
        if r["metric"] in ("cos_psi","Lambda"):
            add(f"event_minus_control_{r['metric']}_lead{abs(float(r['lead'])):g}",r["event_minus_control"],f"95% seed bootstrap {r['seed_bootstrap_ci_low']} to {r['seed_bootstrap_ci_high']}")
    for target in ("dpsi_next","dLambda_next"):
        a=next(r for r in memory if r["model"]=="state_only" and r["target"]==target)
        h=next(r for r in memory if r["model"]=="state_plus_history" and r["target"]==target)
        ph=next(r for r in memory if r["model"]=="state_phase_history" and r["target"]==target)
        add(f"history_rmse_reduction_{target}",1-float(h["normalized_rmse"])/float(a["normalized_rmse"]),"held-out seeds")
        add(f"phase_history_rmse_reduction_{target}",1-float(ph["normalized_rmse"])/float(a["normalized_rmse"]),"held-out seeds")
    val=read_csv(tab/"replay_validation.csv")
    add("replay_snapshots",len(val)); add("replay_noise_z_mean",float(np.mean([float(r["noise_z_mean"]) for r in val])))
    add("replay_noise_z_sd",float(np.mean([float(r["noise_z_sd"]) for r in val]))); add("replay_overrides",sum(int(r["overridden"]) for r in val))
    write_csv(out/"summary.csv",summary)
    print(f"wrote robust summaries for {len(samples)//2} event-control lead samples")

if __name__=="__main__": main()
