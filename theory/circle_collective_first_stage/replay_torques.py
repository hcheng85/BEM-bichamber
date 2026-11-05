#!/usr/bin/env python3
"""Replay exact microscopic torques at a bounded, state-stratified snapshot set."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import os
from pathlib import Path
import sys
import types

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree

from analyze_collective import classify
from matv5_memmap import open_state_rdot
from process_collective import FILE_PATTERN, seed_from_name


def load_production(path: Path):
    # tqdm is imported by the simulator but is only used by optional geometry
    # construction; provide a no-op shim when replaying an existing LU factor.
    if 'tqdm' not in sys.modules:
        shim=types.ModuleType('tqdm'); shim.tqdm=lambda iterable,*a,**k: iterable
        sys.modules['tqdm']=shim
    old=os.getcwd(); os.chdir(path.parent)
    try:
        spec=importlib.util.spec_from_file_location("circle_production",path)
        mod=importlib.util.module_from_spec(spec); assert spec.loader
        spec.loader.exec_module(mod)
    finally: os.chdir(old)
    return mod


def social_split(mod, vn: np.ndarray, state: np.ndarray, Ip: float) -> tuple[np.ndarray,np.ndarray]:
    x,y,a=state; n=len(a); align=np.zeros(n); attraction=np.zeros(n)
    for i in range(n):
        js=np.flatnonzero(vn[i]); dx=x[js]-x[i]; dy=y[js]-y[i]
        rho=np.hypot(dx,dy); bearing=mod.getangle(np.full(len(js),a[i]),dx,dy,rho)
        visual=1+np.cos(bearing); denom=np.sum(visual)
        if denom>0:
            align[i]=np.sum(Ip*np.sin(a[js]-a[i])*visual)/denom
            attraction[i]=np.sum(rho*np.sin(bearing)*visual)/denom
    return align,attraction


def stratified_indices(label: np.ndarray, time: np.ndarray, n_each: int=10) -> list[int]:
    chosen=[]
    for value in (-1,0,1):
        idx=np.flatnonzero((label==value)&(time>=20)&(time<=999))
        if not len(idx): continue
        # Equal-quantile indices span the full trajectory and avoid cherry-picking.
        positions=np.linspace(0,len(idx)-1,min(n_each,len(idx))).round().astype(int)
        chosen.extend(idx[positions].tolist())
    return sorted(set(chosen))


def persistent_label(label: np.ndarray, dt: float, window_s: float=10.0) -> np.ndarray:
    """Retain school/mill labels only when locally persistent; dips are transitional."""
    width=max(3,int(round(window_s/dt)))
    stable=np.zeros(len(label),dtype=np.int8)
    for value in (-1,1):
        fraction=uniform_filter1d((label==value).astype(float),width,mode='nearest')
        stable[fraction>=.9]=value
    return stable


def write_csv(path: Path, rows: list[dict]):
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--outputs',type=Path,default=Path(__file__).resolve().parent/'outputs')
    args=p.parse_args(); out=args.outputs; out.mkdir(parents=True,exist_ok=True)
    mod=load_production(args.source.resolve())
    para=dict(Ip=4.0,In=.5,R=mod.R,light=1.0,dt=.001,Iw=.94,If=.01,
              num=100,delta=.01,total_time=1000,gap=mod.gap,
              radius_ratio=50,gap_ratio=50)
    rows=[]; fish_rows=[]
    files=sorted(args.data_dir.glob(FILE_PATTERN),key=seed_from_name)
    print(f'found {len(files)} trajectory files in {args.data_dir}',flush=True)
    if not files: raise FileNotFoundError(args.data_dir/FILE_PATTERN)
    for file in files:
        seed=seed_from_name(file); z=np.load(out/f'metrics_seed{seed}.npz')
        dt_metric=float(z['time'][1]-z['time'][0]); raw_label,_,_=classify(z['P'],z['Pi'],dt_metric)
        label=persistent_label(raw_label,dt_metric)
        selected=stratified_indices(label,z['time'])
        state_mm,rdot_mm,_=open_state_rdot(file)
        print(f'seed {seed}: {len(selected)} snapshots',flush=True)
        for mi in selected:
            step=int(z['step'][mi]); state=np.asarray(state_mm[:,:,step],dtype=float)
            x,y,theta=state; c=np.cos(theta); sn=np.sin(theta); e=np.vstack((c,sn))
            direct_u,direct_w=mod.GetHydro(state,para)
            bem_u,bem_w=mod.GetHydroBEM(state,para)
            vn=mod.GetVoronoi(state[:2],para)
            align,attraction=social_split(mod,vn,state,para['Ip'])
            social_direct=np.asarray(mod.GetVision(vn,state,para)).ravel()
            wall=np.asarray(mod.GetWallAvoid(state,para)).ravel()
            direct_w=np.asarray(direct_w).ravel(); bem_w=np.asarray(bem_w).ravel()
            trial_rdot=e+direct_u+bem_u
            jump=np.asarray(mod.Check_jumping(state,trial_rdot,para),dtype=bool)
            extreme=np.any(np.abs(trial_rdot)>500,axis=0)
            override=jump|extreme
            effective_rdot=trial_rdot.copy(); effective_rdot[:,override]=e[:,override]
            d_w=direct_w.copy(); b_w=bem_w.copy(); d_w[override]=0; b_w[override]=0
            components={'alignment':align,'attraction':attraction,'wall':wall,
                        'direct_hydro':d_w,'bem_hydro':b_w}
            raw_total=sum(components.values()); effective=np.clip(raw_total,-np.pi/para['dt'],np.pi/para['dt'])
            clip_residual=effective-raw_total; components['clip_residual']=clip_residual

            # Reproduce the later positional safety override.  rdot at index
            # step+1 is the velocity that advances state[step] to state[step+1].
            next_r=state[:2]+effective_rdot*para['dt']
            dis_new=np.hypot(next_r[0],next_r[1])-mod.R
            dis_old=np.hypot(x,y)-mod.R
            sticky=(dis_new-mod.safedis>=0)&(dis_old<dis_new)
            effective_rdot[:,sticky]=0
            stored=np.asarray(rdot_mm[:,:,step+1],dtype=float)
            velocity_rmse=float(np.sqrt(np.mean((effective_rdot-stored)**2)))
            velocity_max_error=float(np.max(np.abs(effective_rdot-stored)))

            px,py=np.mean(c),np.mean(sn); P=np.hypot(px,py); Q=P*P
            projection_factor=(2/para['num'])*(px*(-sn)+py*c)
            bq={name:float(np.sum(projection_factor*w)) for name,w in components.items()}
            bq_noise=-para['In']**2*(Q-1/para['num'])
            bq_det=float(np.sum(projection_factor*effective)); bq_total=bq_det+bq_noise
            rr=np.hypot(x,y); cov=np.cov(state[:2]); eig=np.linalg.eigvalsh(cov)
            aspect=float(np.sqrt(max(eig[-1],1e-15)/max(eig[0],1e-15)))
            nn=float(np.mean(cKDTree(state[:2].T).query(state[:2].T,k=2)[0][:,1]))
            edges=np.transpose(np.nonzero(np.triu(vn,1)))
            neighbor_cos=float(np.mean(np.cos(theta[edges[:,0]]-theta[edges[:,1]])))
            m_heading=float(np.mean((x*sn-y*c)/rr))
            recent=max(0,mi-int(round(1/dt_metric)))
            row=dict(seed=seed,step=step,time=float(z['time'][mi]),
                macrostate={-1:'milling',0:'intermediate',1:'schooling'}[int(label[mi])],
                P=P,Q=Q,Pi=float(z['Pi'][mi]),M=m_heading,
                recent_dQ_1=float(Q-float(z['Q'][recent])),radial_mean=float(np.mean(rr/mod.R)),
                radial_std=float(np.std(rr/mod.R)),com_radius=float(np.hypot(np.mean(x),np.mean(y))/mod.R),
                aspect_ratio=aspect,mean_nearest_neighbor=nn,neighbor_heading_cos=neighbor_cos,
                bQ_alignment=bq['alignment'],bQ_attraction=bq['attraction'],bQ_wall=bq['wall'],
                bQ_direct_hydro=bq['direct_hydro'],bQ_bem_hydro=bq['bem_hydro'],
                bQ_clip_residual=bq['clip_residual'],bQ_deterministic=bq_det,
                bQ_noise_ito=bq_noise,bQ_total=bq_total,
                social_split_max_error=float(np.max(np.abs(align+attraction-social_direct))),
                n_collision_override=int(np.count_nonzero(jump)),n_extreme_override=int(np.count_nonzero(extreme)),
                n_sticky_override=int(np.count_nonzero(sticky)),velocity_replay_rmse=velocity_rmse,
                velocity_replay_max_error=velocity_max_error)
            rows.append(row)
            for i in range(para['num']):
                fish_rows.append(dict(seed=seed,step=step,fish=i,theta=theta[i],radius=rr[i],
                    omega_alignment=align[i],omega_attraction=attraction[i],omega_wall=wall[i],
                    omega_direct_hydro=d_w[i],omega_bem_hydro=b_w[i],omega_clip_residual=clip_residual[i],
                    omega_effective=effective[i],collision_override=int(override[i]),sticky_override=int(sticky[i])))
    write_csv(out/'torque_snapshot_summary.csv',rows)
    write_csv(out/'torque_fish_components.csv',fish_rows)
    print(f'saved {len(rows)} snapshots and {len(fish_rows)} fish rows')


if __name__=='__main__': main()
