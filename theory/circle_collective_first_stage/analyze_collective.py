#!/usr/bin/env python3
"""Classify collective states and make first-stage descriptive figures."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.ndimage import uniform_filter1d

from matv5_memmap import open_state_rdot
from process_collective import FILE_PATTERN, seed_from_name, wrap


SCHOOL_P, SCHOOL_PI = 0.70, 0.75
MILL_P, MILL_PI = 0.50, 0.75
SMOOTH_S = 1.0


def classify(P: np.ndarray, Pi: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    width = max(1, int(round(SMOOTH_S/dt)))
    Ps = uniform_filter1d(P.astype(float), width, mode="nearest")
    Pis = uniform_filter1d(np.nan_to_num(Pi, nan=0.0).astype(float), width, mode="nearest")
    label = np.zeros(len(Ps), dtype=np.int8)
    label[(Ps >= SCHOOL_P) & (Pis <= SCHOOL_PI)] = 1
    label[(Ps <= MILL_P) & (Pis >= MILL_PI)] = -1
    return label, Ps, Pis


def runs(label: np.ndarray) -> list[tuple[int, int, int]]:
    edges = np.r_[0, np.flatnonzero(label[1:] != label[:-1])+1, len(label)]
    return [(int(a), int(b), int(label[a])) for a,b in zip(edges[:-1],edges[1:])]


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows: return
    with path.open("w", newline="") as handle:
        writer=csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--outputs", type=Path, default=Path(__file__).resolve().parent/"outputs")
    args=p.parse_args(); out=args.outputs; out.mkdir(parents=True,exist_ok=True)
    files=sorted(args.data_dir.glob(FILE_PATTERN),key=seed_from_name)
    zs=[np.load(out/f"metrics_seed{seed_from_name(f)}.npz") for f in files]
    labels=[]; occupancy=[]; intervals=[]; transitions=[]
    fig,axes=plt.subplots(5,2,figsize=(12,14),sharex=True,sharey=True)
    for ax,path,z in zip(axes.ravel(),files,zs):
        seed=seed_from_name(path); dt=float(z['time'][1]-z['time'][0])
        lab,Ps,Pis=classify(z['P'],z['Pi'],dt); labels.append(lab)
        ax.plot(z['time'],Ps,lw=.65,label='P')
        ax.plot(z['time'],Pis,lw=.65,label=r'$\Pi$')
        ax.fill_between(z['time'],0,1,where=lab==1,color='#4c78a8',alpha=.08)
        ax.fill_between(z['time'],0,1,where=lab==-1,color='#e45756',alpha=.08)
        ax.set_title(f"seed {seed}"); ax.set_ylim(0,1); ax.grid(alpha=.15)
        rr=runs(lab)
        for a,b,kind in rr:
            intervals.append(dict(seed=seed,start_time=float(z['time'][a]),
                end_time=float(z['time'][b-1]),duration=(b-a)*dt,
                state={1:'schooling',-1:'milling',0:'intermediate'}[kind]))
        valid=z['time']>=20
        occupancy.append(dict(seed=seed,
            schooling_fraction=float(np.mean(lab[valid]==1)),
            milling_fraction=float(np.mean(lab[valid]==-1)),
            intermediate_fraction=float(np.mean(lab[valid]==0)),
            mean_P=float(np.mean(z['P'][valid])),mean_Pi=float(np.nanmean(z['Pi'][valid])),
            mean_Q=float(np.mean(z['Q'][valid])),mean_M=float(np.mean(z['M_heading'][valid]))))

        # A transition requires 20 s of persistent old/new state (>=90% each),
        # preventing brief threshold excursions from being called transitions.
        w=int(round(20/dt))
        for i in range(w,len(lab)-w):
            old=int(np.sign(np.mean(lab[i-w:i]))); new=int(np.sign(np.mean(lab[i:i+w])))
            if old*new==-1 and np.mean(lab[i-w:i]==old)>=.9 and np.mean(lab[i:i+w]==new)>=.9:
                if not transitions or transitions[-1]['seed']!=seed or i-transitions[-1]['metric_index']>w:
                    transitions.append(dict(seed=seed,metric_index=i,time=float(z['time'][i]),
                        from_state='schooling' if old==1 else 'milling',
                        to_state='schooling' if new==1 else 'milling'))
    axes[0,0].legend(ncol=2,fontsize=8); fig.supxlabel('time'); fig.supylabel('order parameter')
    fig.tight_layout(); fig.savefig(out/'time_series_all_seeds.png',dpi=220); plt.close(fig)
    write_csv(out/'state_occupancy.csv',occupancy); write_csv(out/'state_intervals.csv',intervals)
    write_csv(out/'transition_events.csv',transitions)

    fig,axes=plt.subplots(5,2,figsize=(12,14),sharex=True,sharey=True)
    for ax,path,z in zip(axes.ravel(),files,zs):
        seed=seed_from_name(path)
        ax.plot(z['time'],z['Q'],lw=.55,label=r'$Q=P^2$')
        ax.plot(z['time'],z['M_heading'],lw=.55,label=r'$M$ (signed, circle center)')
        ax.set_title(f"seed {seed}"); ax.set_ylim(-1,1); ax.grid(alpha=.15)
    axes[0,0].legend(ncol=2,fontsize=8); fig.supxlabel('time'); fig.supylabel('projected coordinate')
    fig.tight_layout(); fig.savefig(out/'Q_M_time_series.png',dpi=220); plt.close(fig)
    for event in transitions:
        seed=event['seed']; z=zs[seed]; center=event['time']; m=(z['time']>=center-30)&(z['time']<=center+30)
        fig,ax=plt.subplots(figsize=(8,4)); ax.plot(z['time'][m],z['P'][m],label='P'); ax.plot(z['time'][m],z['Pi'][m],label=r'$\Pi$')
        ax.plot(z['time'][m],z['Q'][m],label='Q',alpha=.8); ax.plot(z['time'][m],z['M_heading'][m],label='M',alpha=.75)
        ax.axvline(center,color='black',ls='--',lw=1,label='persistent-state boundary')
        ax.set(xlabel='time',ylabel='order parameter',title=f"seed {seed}: {event['from_state']} to {event['to_state']}")
        ax.legend(ncol=3); ax.grid(alpha=.15); fig.tight_layout(); fig.savefig(out/f'transition_seed{seed}_{center:.2f}.png',dpi=220); plt.close(fig)

    # Joint distribution after formation transient.
    Pall=np.concatenate([z['P'][z['time']>=20] for z in zs])
    Piall=np.concatenate([z['Pi'][z['time']>=20] for z in zs])
    fig,ax=plt.subplots(figsize=(6,5))
    h=ax.hist2d(Pall,Piall,bins=90,range=[[0,1],[0,1]],norm=LogNorm(),cmap='magma')
    ax.axvline(SCHOOL_P,color='#4c78a8',ls='--',lw=1); ax.axvline(MILL_P,color='#e45756',ls='--',lw=1)
    ax.axhline(MILL_PI,color='white',ls=':',lw=1)
    ax.set(xlabel='polarization P',ylabel=r'original milling $\Pi$',title='All seeds, t >= 20')
    fig.colorbar(h[3],ax=ax,label='samples'); fig.tight_layout(); fig.savefig(out/'joint_P_Pi.png',dpi=220); plt.close(fig)

    # State-conditioned particle-level radial and relative-heading histograms.
    radial_edges=np.linspace(0,1.05,107); angle_edges=np.linspace(-np.pi,np.pi,145)
    hist_r=np.zeros((3,10,len(radial_edges)-1),dtype=np.int64)
    hist_a=np.zeros((3,10,len(angle_edges)-1),dtype=np.int64)
    label_to_row={-1:0,0:1,1:2}
    for path,z,lab in zip(files,zs,labels):
        seed=seed_from_name(path); state,_,_=open_state_rdot(path)
        # The stored descriptive distribution cadence is 0.1 time units.
        steps=np.arange(0,state.shape[2],100,dtype=np.int64)
        for start in range(0,len(steps),1000):
            st=steps[start:start+1000]; s=np.asarray(state[:,:,st])
            # Advanced indexing moves time first in some NumPy releases.
            if s.shape[0]!=3: s=np.moveaxis(s,-1,0)
            x,y,th=s[0],s[1],s[2]
            metric_i=st//10; l=lab[metric_i]
            for state_value,row in label_to_row.items():
                choose=l==state_value
                if not np.any(choose): continue
                rr=np.hypot(x[:,choose],y[:,choose])/7.158910531638177
                aa=wrap(th[:,choose]-np.arctan2(y[:,choose],x[:,choose]))
                hist_r[row,seed]+=np.histogram(rr,radial_edges)[0]
                hist_a[row,seed]+=np.histogram(aa,angle_edges)[0]
    np.savez_compressed(out/'state_conditioned_distributions.npz',radial_edges=radial_edges,
        heading_edges=angle_edges,radial_hist_by_state_seed=hist_r,heading_hist_by_state_seed=hist_a,
        state_names=np.array(['milling','intermediate','schooling']))
    fig,axs=plt.subplots(1,2,figsize=(11,4.2)); colors=['#e45756','#999999','#4c78a8']
    for row,(name,color) in enumerate(zip(['milling','intermediate','schooling'],colors)):
        rden=hist_r[row].sum(0); rden=rden/(rden.sum()*np.diff(radial_edges))
        aden=hist_a[row].sum(0); aden=aden/(aden.sum()*np.diff(angle_edges))
        axs[0].plot((radial_edges[:-1]+radial_edges[1:])/2,rden,label=name,color=color)
        axs[1].plot((angle_edges[:-1]+angle_edges[1:])/2,aden,label=name,color=color)
    axs[0].set(xlabel='r/R',ylabel='density',title='radial distribution')
    axs[1].set(xlabel=r'heading minus radial angle',ylabel='density',title='relative-heading distribution')
    axs[1].legend(); [a.grid(alpha=.15) for a in axs]; fig.tight_layout(); fig.savefig(out/'radial_heading_distributions.png',dpi=220); plt.close(fig)

    config=dict(metric_dt=float(zs[0]['time'][1]-zs[0]['time'][0]),smooth_s=SMOOTH_S,
        schooling_P_min=SCHOOL_P,schooling_Pi_max=SCHOOL_PI,milling_P_max=MILL_P,
        milling_Pi_min=MILL_PI,persistence_s=20,formation_transient_s=20)
    import json; (out/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    print(f"saved {out}; persistent transitions={len(transitions)}")


if __name__=='__main__': main()
