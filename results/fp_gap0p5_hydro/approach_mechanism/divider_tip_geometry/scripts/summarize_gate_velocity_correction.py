#!/usr/bin/env python3
"""Compare observed-velocity and counterfactual-velocity tangent drifts.

Only class frequencies and the prespecified angular/spatial prospective
associations are recomputed.  Original tables are read-only inputs.
"""
from __future__ import annotations

import argparse,csv
from collections import Counter
from pathlib import Path

import numpy as np

from summarize_tip_mechanism import (association,class_frequencies,
    event_control_matching,geometry_matched_class_association,read,write)


def key(r):
    return (r["seed"],r["pair_id"],r["condition"],r["observation_offset"],r["tip"],r["epsilon"])


def transition_table(old,new):
    a={key(r):r for r in old};b={key(r):r for r in new};counts=Counter()
    for k in sorted(a.keys()&b.keys()):counts[(a[k]["classification"],b[k]["classification"],k[-1])]+=1
    return [dict(epsilon=e,old_class=o,corrected_class=n,count=c)
            for (o,n,e),c in sorted(counts.items())]


def drift_change(old,new):
    a={key(r):r for r in old};rows=[]
    for r in new:
      if r["classification"]=="invalid":continue
      o=a[key(r)]; tip=r["tip"];observed=float(r["alpha_plus_dot"] if tip=="upper" else r["alpha_minus_dot"])
      for side in ("out","in"):
       trial=float(r[f"{side}_counterfactual_alpha_dot"])
       rows.append(dict(seed=r["seed"],pair_id=r["pair_id"],condition=r["condition"],
         observation_offset=r["observation_offset"],tip=tip,epsilon=r["epsilon"],side=side,
         observed_alpha_dot=observed,counterfactual_alpha_dot=trial,
         alpha_dot_change=trial-observed,old_mu=float(o[f"mu_{side}"]),
         corrected_mu=float(r[f"mu_{side}"]),mu_change=float(r[f"mu_{side}"])-float(o[f"mu_{side}"])))
    return rows


def metric_map(rows):
    out={}
    for r in rows:
      k=(r.get("analysis"),str(r.get("offset")),str(r.get("horizon")),r.get("outcome"))
      out[k]=r
    return out


def comparison(old_assoc,new_assoc):
    a=metric_map(old_assoc);b=metric_map(new_assoc);rows=[]
    for k in sorted(a.keys()&b.keys(),key=str):
      x,y=a[k],b[k]
      rows.append(dict(analysis=k[0],offset=k[1],horizon=k[2],outcome=k[3],
        old_seed_mean_risk_difference=x["seed_mean_risk_difference"],
        corrected_seed_mean_risk_difference=y["seed_mean_risk_difference"],
        change=float(y["seed_mean_risk_difference"])-float(x["seed_mean_risk_difference"]),
        corrected_ci_low=y["seed_bootstrap_ci_low"],corrected_ci_high=y["seed_bootstrap_ci_high"],
        corrected_n_pairs=y.get("n_pairs",""),corrected_match_distance=y.get("mean_standardized_match_distance","")))
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--analysis-dir",type=Path,required=True);args=ap.parse_args()
    root=args.analysis_dir.resolve();tab=root/"tables"
    old=read(tab/"tangent_snapshot_replay.csv");new=read(tab/"tangent_snapshot_replay_corrected.csv")
    oldp=[r for r in old if r["epsilon"]=="0.001"];newp=[r for r in new if r["epsilon"]=="0.001"]
    freq=class_frequencies(newp);unadj=association(newp);pairs,matched=geometry_matched_class_association(newp)
    ecpairs,ecsummary=event_control_matching(newp)
    angular=[r for r in unadj+matched if r["outcome"]=="angular_entry"]
    spatial=[r for r in unadj+matched if r["outcome"] in ("focal_crossing","any_crossing")]
    write(tab/"corrected_class_frequencies.csv",freq)
    write(tab/"corrected_angular_entry_association.csv",angular)
    write(tab/"corrected_spatial_crossing_association.csv",spatial)
    write(tab/"corrected_geometry_matched_class_pairs.csv",pairs)
    write(tab/"corrected_event_control_geometry_matches.csv",ecpairs)
    write(tab/"corrected_event_control_matched_summary.csv",ecsummary)
    trans=transition_table(old,new);write(tab/"gate_velocity_class_transitions.csv",trans)
    changes=drift_change(old,new)
    oldassoc=read(tab/"prospective_association.csv")
    write(tab/"gate_velocity_association_comparison.csv",comparison(oldassoc,unadj+matched))
    overrides=sum(int(r.get("out_override",0) or 0)+int(r.get("in_override",0) or 0)
                  for r in new if r["classification"]!="invalid")
    valid=sum(r["classification"]!="invalid" for r in new)
    eps1=[r for r in trans if r["epsilon"]=="0.001"]
    changed=sum(int(r["count"]) for r in eps1 if r["old_class"]!=r["corrected_class"])
    pin_old=sum(r["classification"]=="pinning" for r in oldp if r["classification"]!="invalid")
    pin_new=sum(r["classification"]=="pinning" for r in newp if r["classification"]!="invalid")
    dc=np.array([float(r["alpha_dot_change"]) for r in changes if r["epsilon"]=="0.001"])
    audit=[dict(metric="valid_trial_sides_all_eps",value=valid*2),
      dict(metric="collision_or_speed_override_activations",value=overrides),
      dict(metric="production_override_translation_rule",value="rdot=e when override is active"),
      dict(metric="production_override_angular_rule",value="remove direct and BEM hydro torque"),
      dict(metric="epsilon_1e3_class_changes",value=changed),
      dict(metric="epsilon_1e3_valid_evaluations",value=sum(int(r["count"]) for r in eps1 if r["old_class"]!="invalid")),
      dict(metric="old_pinning_count",value=pin_old),dict(metric="corrected_pinning_count",value=pin_new),
      dict(metric="alpha_dot_change_median",value=np.median(dc)),
      dict(metric="alpha_dot_change_abs_q95",value=np.quantile(np.abs(dc),.95)),
      dict(metric="alpha_dot_change_abs_max",value=np.max(np.abs(dc)))]
    write(tab/"counterfactual_gate_velocity_audit.csv",audit)

if __name__=="__main__":main()
