#!/usr/bin/env python3
"""Test the proposed linear polarization response on milling ensembles."""
from __future__ import annotations

import argparse, csv, json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from validate_drift_consistency import ExactReplay, write_csv


def fit_equivariant(rows):
    p=np.array([[r['p_x'],r['p_y']] for r in rows]); f=np.array([[r['polarization_source_x'],r['polarization_source_y']] for r in rows])
    jp=np.column_stack([-p[:,1],p[:,0]])
    X=np.column_stack([p.ravel(),jp.ravel()]).reshape(-1,2)  # replaced below by explicit stacking
    # Each vector observation contributes x and y scalar equations.
    X=np.vstack([np.column_stack([p[:,0],-p[:,1]]),np.column_stack([p[:,1],p[:,0]])])
    y=np.r_[f[:,0],f[:,1]]
    a,b=np.linalg.lstsq(X,y,rcond=None)[0]
    return float(a),float(b)


def predict(rows,a,b):
    p=np.array([[r['p_x'],r['p_y']] for r in rows]); f=np.array([[r['polarization_source_x'],r['polarization_source_y']] for r in rows])
    pred=a*p+b*np.column_stack([-p[:,1],p[:,0]])
    sse=float(np.sum((f-pred)**2)); zero=float(np.sum(f**2)); centered=float(np.sum((f-f.mean(0))**2))
    return dict(n=len(rows),rmse=float(np.sqrt(np.mean((f-pred)**2))),
                r2_against_zero=1-sse/zero,r2_centered=1-sse/centered)


def main():
    p=argparse.ArgumentParser(); p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--source',type=Path,default=HERE.parents[1].parent/'BEM_circular_parallel_instance.py')
    p.add_argument('--output-dir',type=Path,default=HERE/'outputs'/'drift_consistency')
    args=p.parse_args(); args.output_dir.mkdir(parents=True,exist_ok=True)
    old=list(csv.DictReader((args.output_dir/'uniform_replay_snapshots.csv').open()))
    wanted=[r for r in old if r['interval'] in ('milling','milling_post_transition')]
    replay=ExactReplay(args.data_dir,args.source.resolve(),HERE/'outputs'); rows=[]
    for j,r in enumerate(wanted,1):
        row=replay.evaluate(int(r['seed']),int(r['step'])); rows.append(dict(interval=r['interval'],**row))
        if j%200==0: print(f'{j}/{len(wanted)} milling response states',flush=True)
    write_csv(args.output_dir/'polarization_response_snapshots.csv',rows)
    by_seed={s:[r for r in rows if r['seed']==s] for s in (7,8,9)}
    tests=[]
    for train,test in [(7,8),(8,7)]:
        a,b=fit_equivariant(by_seed[train]); tests.append(dict(train=f'seed{train}',test=f'seed{test}',a=a,b=b,**predict(by_seed[test],a,b)))
    train=by_seed[7]+by_seed[8]; a,b=fit_equivariant(train)
    tests.append(dict(train='seed7+seed8',test='seed7+seed8',a=a,b=b,**predict(train,a,b)))
    tests.append(dict(train='seed7+seed8',test='seed9_post',a=a,b=b,**predict(by_seed[9],a,b)))

    # Contiguous trajectory-block bootstrap for the real growth coefficient a.
    blocks=[]
    for seed in (7,8):
        rr=by_seed[seed]; cuts=np.array_split(np.arange(len(rr)),10); blocks.extend([[rr[i] for i in cut] for cut in cuts])
    rng=np.random.default_rng(29104); boot=[]
    for _ in range(4000):
        sampled=[]
        for idx in rng.integers(0,len(blocks),len(blocks)): sampled.extend(blocks[idx])
        boot.append(fit_equivariant(sampled)[0])
    ci=[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))]
    # Residual correlations flag hidden geometry without expanding the model.
    pred=np.array([[a*r['p_x']-b*r['p_y'],b*r['p_x']+a*r['p_y']] for r in train])
    observed=np.array([[r['polarization_source_x'],r['polarization_source_y']] for r in train]); resid=np.linalg.norm(observed-pred,axis=1)
    correlations={k:float(np.corrcoef(resid,[r[k] for r in train])[0,1]) for k in ['radial_mean','radial_std','com_radius','M','neighbor_heading_cos']}
    result=dict(model='F(p) = a p + b J p, zero intercept, rotationally equivariant',
        training_seeds=[7,8],a=a,b=b,a_block_bootstrap_95pct=ci,noise_decay_threshold=0.5**2/2,
        a_minus_noise_threshold=a-0.5**2/2,eigenvalues_real_part=a,
        heldout_tests=tests,residual_magnitude_correlations=correlations,
        verdict='not_closed' if min(t['r2_against_zero'] for t in tests if t['test'] in ('seed7','seed8','seed9_post'))<0.5 else 'provisionally_predictive')
    (args.output_dir/'polarization_response_matrix.json').write_text(json.dumps(result,indent=2)+'\n')

    fig,axs=plt.subplots(1,2,figsize=(10,4.2)); pnorm=np.sqrt(np.array([r['Q'] for r in train]));
    pvec=np.array([[r['p_x'],r['p_y']] for r in train]); fvec=observed; q=np.maximum(pnorm**2,1e-12)
    radial=np.sum(pvec*fvec,axis=1)/np.sqrt(q); tangential=(pvec[:,0]*fvec[:,1]-pvec[:,1]*fvec[:,0])/np.sqrt(q)
    axs[0].scatter(pnorm,radial,s=12,alpha=.45); axs[0].plot(np.sort(pnorm),a*np.sort(pnorm),color='black',label='linear response')
    axs[1].scatter(pnorm,tangential,s=12,alpha=.45); axs[1].plot(np.sort(pnorm),b*np.sort(pnorm),color='black',label='linear response')
    axs[0].set(xlabel='|p|',ylabel='radial deterministic source',title='growth component')
    axs[1].set(xlabel='|p|',ylabel='tangential deterministic source',title='rotation component')
    for ax in axs: ax.legend(); ax.grid(alpha=.15)
    fig.tight_layout(); fig.savefig(args.output_dir/'polarization_response_matrix.png',dpi=220); plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
