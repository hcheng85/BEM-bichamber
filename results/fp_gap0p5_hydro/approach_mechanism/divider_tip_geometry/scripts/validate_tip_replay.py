#!/usr/bin/env python3
"""Direct validation of focal affine replay against full simulator functions."""
from __future__ import annotations

import argparse, csv
from pathlib import Path
import sys

import numpy as np
import scipy.io as sio

HERE=Path(__file__).resolve().parent
REPLAY=HERE.parents[1]/"pre_crossing_turning"/"scripts"
sys.path.insert(0,str(REPLAY))
from replay_dynamics import (boundary_hydrodynamics, direct_hydrodynamics,
                             social_terms, wall_avoidance)  # noqa: E402
from exact_tip_dynamics import (Parameters, analytic_geometry, chamber_side,
    FocalBEM, focal_replay, local_to_global, tangent_geometry, visible_neighbors)


def read(path):
    with open(path,newline="") as f:return list(csv.DictReader(f))


def write(path,rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--analysis-dir",type=Path,required=True)
    ap.add_argument("--raw-dir",type=Path,required=True);ap.add_argument("--geometry",type=Path,required=True)
    ap.add_argument("--lu",type=Path,required=True);ap.add_argument("--piv",type=Path,required=True);args=ap.parse_args()
    root=args.analysis_dir.resolve();cand=read(root/"tables"/"geometry_candidates.csv")
    p=Parameters();geom=analytic_geometry(p);g=sio.loadmat(args.geometry);normal=np.asarray(g["normal_vec"])
    boundary=dict(x=np.asarray(g["Sx"]).reshape(-1),y=np.asarray(g["Sy"]).reshape(-1),
        weights=np.asarray(g["weights"]).reshape(-1),nx=normal[:,0],ny=normal[:,1])
    lu=np.load(args.lu,mmap_mode="r");piv=np.load(args.piv);rows=[]
    for seed in (0,4,8,12,16):
        choices=[r for r in cand if int(r["seed"])==seed and r["condition"]=="event" and r["offset"]=="1.8" and r["valid"]=="1"]
        r=choices[0];idx=int(r["index"]);focal=int(r["fish"]);origin=r["origin"];side=chamber_side(origin)
        raw=args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
        st=state[:,:,idx];vel=(state[:2,:,idx+1]-state[:2,:,idx])/p.dt
        gate=tangent_geometry(st[:2,focal],vel[:,focal],origin,p)
        nextgate=tangent_geometry(state[:2,focal,idx+1],np.zeros(2),origin,p)
        visible=visible_neighbors(st[:2],geom,p);bem=FocalBEM.build(st,focal,boundary,lu,piv,p)
        for tip,alpha,adot in (("upper",gate["alpha_plus"],gate["alpha_plus_dot"]),
                               ("lower",gate["alpha_minus"],gate["alpha_minus_dot"])):
          fd=((nextgate["alpha_plus"]-alpha)/p.dt if tip=="upper" else (nextgate["alpha_minus"]-alpha)/p.dt)
          for direction,sgn in (("out",1 if tip=="upper" else -1),("in",-1 if tip=="upper" else 1)):
            theta=local_to_global(alpha+sgn*1e-3,side);trial=st.copy();trial[2,focal]=theta
            rep=focal_replay(st,focal,theta,visible,bem,geom,p)
            ub,wb=boundary_hydrodynamics(trial,boundary,lu,piv,p)
            ud,wd=direct_hydrodynamics(trial,p);al,at,_=social_terms(trial,geom,p);wa,_=wall_avoidance(trial,geom,p)
            rawtotal=wa[focal]+al[focal]+at[focal]+wd[focal]+wb[focal]
            fulltotal=float(np.clip(rawtotal,-np.pi/p.dt,np.pi/p.dt))
            rdot=np.array([np.cos(theta),np.sin(theta)])+ud[:,focal]+ub[:,focal]
            # The chosen validation snapshots are away from safeguard activation;
            # retain an explicit flag instead of silently comparing another rule.
            override=bool(rep["override"])
            if override:fulltotal=float(np.clip(wa[focal]+al[focal]+at[focal],-np.pi/p.dt,np.pi/p.dt))
            rows.append(dict(seed=seed,index=idx,focal_fish=focal,tip=tip,direction=direction,
              boundary_velocity_max_abs_error=np.max(np.abs(rep["rdot"]-(np.array([np.cos(theta),np.sin(theta)])+ud[:,focal]+ub[:,focal]))) if not override else np.nan,
              boundary_torque_abs_error=abs(rep["terms"]["boundary_hydro"]-wb[focal]),
              direct_torque_abs_error=abs(rep["terms"]["direct_hydro"]-wd[focal]),
              alignment_abs_error=abs(rep["terms"]["alignment"]-al[focal]),
              attraction_abs_error=abs(rep["terms"]["attraction"]-at[focal]),
              wall_abs_error=abs(rep["terms"]["wall"]-wa[focal]),
              total_abs_error=abs(rep["total"]-fulltotal),override=int(override),
              alpha_dot_analytic=adot,alpha_dot_forward_difference=fd,
              alpha_dot_abs_error=abs(adot-fd)))
    write(root/"tables"/"exact_replay_validation.csv",rows)

if __name__=="__main__":main()
