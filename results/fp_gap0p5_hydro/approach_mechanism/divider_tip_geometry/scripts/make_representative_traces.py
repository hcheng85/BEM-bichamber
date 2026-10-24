#!/usr/bin/env python3
"""Exact pre-crossing traces for a bounded set of representative events."""
from __future__ import annotations

import argparse,csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

from exact_tip_dynamics import (Parameters,analytic_geometry,chamber_side,classify,exact_wall,
    FocalBEM,focal_replay,global_to_local,local_to_global,ray_wall_branch,tangent_geometry,visible_neighbors)

DT=.001

def read(path):
 with open(path,newline="") as f:return list(csv.DictReader(f))

def write(path,rows):
 fields=[]
 for r in rows:
  for k in r:
   if k not in fields:fields.append(k)
 with open(path,"w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def choose(rows):
 q=[r for r in rows if r["epsilon"]=="0.001" and r["condition"]=="event" and
    r["observation_offset"]=="1.5" and r["initially_outside"]=="1" and
    r["focal_crossing_0p5"]=="1" and r["classification"]!="invalid"]
 chosen=[];used=set()
 targets=[("pinning","upper"),("pinning","lower"),("transmission","upper"),
          ("transmission","lower"),("transmission","upper"),("deflection","upper"),
          ("deflection","lower"),("deflection","lower")]
 for cls,tip in targets:
  candidates=[r for r in q if r["classification"]==cls and r["tip"]==tip and int(r["seed"]) not in used]
  if not candidates:candidates=[r for r in q if r["classification"]==cls and r["tip"]==tip]
  if not candidates:continue
  # Prefer modest initial error, avoiding a cherry-picked extreme.
  candidates.sort(key=lambda r:abs(float(r["initial_q"])-np.median([float(x["initial_q"]) for x in candidates])))
  pick=candidates[0];chosen.append(pick);used.add(int(pick["seed"]))
 return chosen

def main():
 ap=argparse.ArgumentParser();ap.add_argument("--analysis-dir",type=Path,required=True)
 ap.add_argument("--raw-dir",type=Path,required=True);ap.add_argument("--geometry",type=Path,required=True)
 ap.add_argument("--lu",type=Path,required=True);ap.add_argument("--piv",type=Path,required=True);args=ap.parse_args()
 root=args.analysis_dir.resolve();allrows=read(root/"tables"/"tangent_snapshot_replay_corrected.csv");chosen=choose(allrows)
 p=Parameters();geom=analytic_geometry(p);g=sio.loadmat(args.geometry);normal=np.asarray(g["normal_vec"])
 boundary=dict(x=np.asarray(g["Sx"]).reshape(-1),y=np.asarray(g["Sy"]).reshape(-1),weights=np.asarray(g["weights"]).reshape(-1),nx=normal[:,0],ny=normal[:,1])
 lu=np.load(args.lu,mmap_mode="r");piv=np.load(args.piv);trace=[];crossings=[]
 for seed in sorted({int(r["seed"]) for r in chosen}):
  raw=args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat";state=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
  for sel in [r for r in chosen if int(r["seed"])==seed]:
   focal=int(sel["focal_fish"]);side=chamber_side(sel["origin"]);A=float(sel["risk_time"])+2.;tip=sel["tip"]
   i0=int(round((A-.6)/DT));i1=int(round((A+.02)/DT));inds=np.arange(i0,i1+1,10)
   # Exact crossing index and location for the fixed focal fish.
   fine=np.arange(i0,i1+1);hit=fine[side*state[0,focal,fine]>=0]
   crossidx=int(hit[0]) if len(hit) else None
   crossings.append(dict(seed=seed,pair_id=sel["pair_id"],focal_fish=focal,tip=tip,
     reference_class=sel["classification"],crossing_time=(crossidx*DT-A if crossidx else ""),
     crossing_x=(state[0,focal,crossidx] if crossidx else ""),crossing_y=(state[1,focal,crossidx] if crossidx else "")))
   for idx in inds:
    st=state[:,:,idx];vel=(state[:2,:,idx+1]-state[:2,:,idx])/DT;gate=tangent_geometry(st[:2,focal],vel[:,focal],sel["origin"],p)
    theta_local=float(global_to_local(st[2,focal],side));alpha=gate["alpha_plus"] if tip=="upper" else gate["alpha_minus"]
    q=(theta_local-alpha) if tip=="upper" else (alpha-theta_local);sgn=1 if tip=="upper" else -1
    visible=visible_neighbors(st[:2],geom,p);bem=FocalBEM.build(st,focal,boundary,lu,piv,p);rep=focal_replay(st,focal,st[2,focal],visible,bem,geom,p)
    wall_q=sgn*side*rep["terms"]["wall"]
    nonwall=rep["terms"]["alignment"]+rep["terms"]["attraction"]+rep["terms"]["direct_hydro_effective"]+rep["terms"]["boundary_hydro_effective"]
    nonwall_q=sgn*side*nonwall
    branch=ray_wall_branch(st[:,focal],geom,p)["branch"]
    if gate["valid"]:
     if tip=="upper":
      tout,tin=alpha+1e-3,alpha-1e-3
      rout=focal_replay(st,focal,local_to_global(tout,side),visible,bem,geom,p)
      rin=focal_replay(st,focal,local_to_global(tin,side),visible,bem,geom,p)
      gout=tangent_geometry(st[:2,focal],rout["rdot"],sel["origin"],p)
      gin=tangent_geometry(st[:2,focal],rin["rdot"],sel["origin"],p)
      muout=side*rout["total"]-gout["alpha_plus_dot"];muin=side*rin["total"]-gin["alpha_plus_dot"]
     else:
      tout,tin=alpha-1e-3,alpha+1e-3
      rout=focal_replay(st,focal,local_to_global(tout,side),visible,bem,geom,p)
      rin=focal_replay(st,focal,local_to_global(tin,side),visible,bem,geom,p)
      gout=tangent_geometry(st[:2,focal],rout["rdot"],sel["origin"],p)
      gin=tangent_geometry(st[:2,focal],rin["rdot"],sel["origin"],p)
      muout=gout["alpha_minus_dot"]-side*rout["total"];muin=gin["alpha_minus_dot"]-side*rin["total"]
     instclass=classify(muout,muin)
    else:muout=muin=np.nan;instclass="invalid"
    trace.append(dict(seed=seed,pair_id=sel["pair_id"],focal_fish=focal,tip=tip,
      reference_class=sel["classification"],time_from_crossing=idx*DT-A,q=q,
      corridor_width=gate["width"],corridor_valid=int(gate["valid"]),selected_ray_branch=branch,
      instantaneous_tangent_class=instclass,mu_out=muout,mu_in=muin,
      q_directed_wall_torque=wall_q,q_directed_nonwall_torque=nonwall_q,
      q_directed_gate_velocity=(-gate["alpha_plus_dot"] if tip=="upper" else gate["alpha_minus_dot"]),
      divider_distance=gate["d"],x=st[0,focal],y=st[1,focal],crossed=int(side*st[0,focal]>=0)))
 write(root/"tables"/"representative_trace_selection.csv",chosen);write(root/"tables"/"representative_traces.csv",trace);write(root/"tables"/"representative_crossings.csv",crossings)
 figdir=root/"figures";figdir.mkdir(exist_ok=True)
 fig,axs=plt.subplots(4,2,figsize=(12,12),sharex=True);axs=axs.ravel()
 for ax,sel in zip(axs,chosen):
  q=[r for r in trace if r["seed"]==int(sel["seed"]) and r["pair_id"]==sel["pair_id"] and r["tip"]==sel["tip"]]
  t=[r["time_from_crossing"] for r in q];ax.plot(t,[r["q"] for r in q],color="k",lw=2,label="q to tangent");ax.axhline(0,color="k",ls=":",lw=.8);ax.axvline(0,color="#d62728",lw=1)
  ax2=ax.twinx();ax2.plot(t,[r["q_directed_wall_torque"] for r in q],color="#1b9e77",alpha=.8,label="wall");ax2.plot(t,[r["q_directed_nonwall_torque"] for r in q],color="#7570b3",alpha=.75,label="non-wall")
  ax.set_title(f"seed {sel['seed']} {sel['tip']} {sel['classification']}");ax.set_ylabel("q");ax2.set_ylabel("q-directed torque")
 axs[-2].set_xlabel("time from first crossing");axs[-1].set_xlabel("time from first crossing");fig.tight_layout();fig.savefig(figdir/"04_representative_q_torque_traces.png",dpi=220);fig.savefig(figdir/"04_representative_q_torque_traces.pdf");plt.close(fig)
 branches=sorted({r["selected_ray_branch"] for r in trace});bc={b:i for i,b in enumerate(branches)}
 fig,axs=plt.subplots(4,2,figsize=(12,12),sharex=True);axs=axs.ravel()
 for ax,sel in zip(axs,chosen):
  q=[r for r in trace if r["seed"]==int(sel["seed"]) and r["pair_id"]==sel["pair_id"] and r["tip"]==sel["tip"]];t=[r["time_from_crossing"] for r in q]
  ax.plot(t,[r["corridor_width"] for r in q],label="corridor width",color="#1f77b4");ax.plot(t,[r["divider_distance"] for r in q],label="divider distance",color="#ff7f0e");ax.axvline(0,color="#d62728",lw=1)
  ax2=ax.twinx();ax2.step(t,[bc[r["selected_ray_branch"]] for r in q],where="mid",color=".45",alpha=.5);ax2.set_yticks(range(len(branches)),branches,fontsize=5)
  ax.set_title(f"seed {sel['seed']} y_cross={next((c['crossing_y'] for c in crossings if c['seed']==int(sel['seed']) and c['pair_id']==sel['pair_id'] and c['tip']==sel['tip']), ''):.3g}");ax.legend(fontsize=6)
 axs[-2].set_xlabel("time from first crossing");axs[-1].set_xlabel("time from first crossing");fig.tight_layout();fig.savefig(figdir/"05_representative_geometry_branch_traces.png",dpi=220);fig.savefig(figdir/"05_representative_geometry_branch_traces.pdf");plt.close(fig)

if __name__=="__main__":main()
