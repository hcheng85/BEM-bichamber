#!/usr/bin/env python3
"""Summarize projected Q drift and screen (Q,M) closure without a potential fit."""
from __future__ import annotations

import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def read(path):
    rows=list(csv.DictReader(path.open())); keys=list(rows[0]); data={}
    for k in keys:
        try: data[k]=np.array([float(r[k]) for r in rows])
        except ValueError: data[k]=np.array([r[k] for r in rows])
    return data


def seed_boot_mean(x,seed,n=4000,rng=None):
    rng=np.random.default_rng(87231) if rng is None else rng; us=np.unique(seed)
    means=np.array([np.mean(x[seed==s]) for s in us]); draws=rng.choice(means,(n,len(means)),replace=True).mean(1)
    return float(np.mean(means)),float(np.quantile(draws,.025)),float(np.quantile(draws,.975))


def cluster_boot_corr(x,y,seed,n=4000):
    rng=np.random.default_rng(9102); us=np.unique(seed); vals=[]
    for _ in range(n):
        take=rng.choice(us,len(us),replace=True); xx=[]; yy=[]
        for s in take: xx.extend(x[seed==s]); yy.extend(y[seed==s])
        if np.std(xx)>0 and np.std(yy)>0: vals.append(np.corrcoef(xx,yy)[0,1])
    return float(np.corrcoef(x,y)[0,1]),float(np.quantile(vals,.025)),float(np.quantile(vals,.975))


def main():
    p=argparse.ArgumentParser(); p.add_argument('--outputs',type=Path,default=Path(__file__).resolve().parent/'outputs')
    args=p.parse_args(); out=args.outputs; d=read(out/'torque_snapshot_summary.csv')
    seed=d['seed'].astype(int); states=['schooling','intermediate','milling']
    comps=['bQ_alignment','bQ_attraction','bQ_wall','bQ_direct_hydro','bQ_bem_hydro','bQ_noise_ito']
    labels=['alignment','attraction','wall','direct hydro','BEM hydro','noise (Itô)']
    summary=[]
    for state in states:
        m=d['macrostate']==state
        for k in comps+['bQ_deterministic','bQ_total']:
            mean,lo,hi=seed_boot_mean(d[k][m],seed[m])
            summary.append(dict(section='state_drift',state=state,quantity=k,value=mean,ci_low=lo,ci_high=hi,n=int(m.sum()),n_seeds=len(np.unique(seed[m]))))

    # Fixed, predeclared bins; no data-driven bandwidth or fitted coefficients.
    qedges=np.array([0,.1,.3,.6,.8,1.000001]); medges=np.array([-1,-.5,0,.5,1.000001])
    qbin=np.digitize(d['Q'],qedges)-1; mbin=np.digitize(d['M'],medges)-1; cell=qbin*10+mbin
    eligible=np.zeros(len(seed),bool)
    for c in np.unique(cell):
        m=cell==c
        if m.sum()>=4 and len(np.unique(seed[m]))>=2:
            eligible[m]=True
    # Cell fixed effects plus within-cell linear Q and M terms remove residual
    # gradients caused by finite bin width. This is a conditioning operation,
    # not a fitted dynamical model or potential.
    cells=np.unique(cell[eligible]); onehot=np.column_stack([(cell[eligible]==c).astype(float) for c in cells])
    X=np.column_stack([onehot,d['Q'][eligible],d['M'][eligible]])
    def nuisance_residual(values):
        beta=np.linalg.lstsq(X,values[eligible],rcond=None)[0]
        return values[eligible]-X@beta
    residual_eligible=nuisance_residual(d['bQ_total'])
    residual=np.full(len(seed),np.nan); residual[eligible]=residual_eligible
    fraction=float(np.var(residual[eligible],ddof=1)/np.var(d['bQ_total'][eligible],ddof=1))
    summary.append(dict(section='closure',state='all',quantity='within_QM_unexplained_variance_fraction',
                        value=fraction,ci_low='',ci_high='',n=int(eligible.sum()),n_seeds=len(np.unique(seed[eligible]))))
    features=['recent_dQ_1','radial_std','com_radius','aspect_ratio','mean_nearest_neighbor','neighbor_heading_cos','Pi']
    closure_rows=[]
    for f in features:
        # Apply the identical nuisance conditioning to each covariate.
        fr=nuisance_residual(d[f])
        r,lo,hi=cluster_boot_corr(fr,residual_eligible,seed[eligible])
        closure_rows.append(dict(feature=f,correlation_with_bQ_residual=r,ci_low=lo,ci_high=hi,
                                 n=int(eligible.sum()),n_seeds=len(np.unique(seed[eligible]))))
        summary.append(dict(section='closure',state='all',quantity='residual_corr_'+f,value=r,ci_low=lo,ci_high=hi,
                            n=int(eligible.sum()),n_seeds=len(np.unique(seed[eligible]))))

    with (out/'summary.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(summary[0]),lineterminator='\n'); w.writeheader(); w.writerows(summary)
    with (out/'closure_screen.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(closure_rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(closure_rows)

    fig,ax=plt.subplots(figsize=(10,4.8)); x=np.arange(len(comps)); width=.24
    colors=['#4c78a8','#999999','#e45756']
    for j,(state,color) in enumerate(zip(states,colors)):
        vals=[]; los=[]; his=[]
        for k in comps:
            row=next(r for r in summary if r['state']==state and r['quantity']==k)
            vals.append(row['value']); los.append(row['value']-row['ci_low']); his.append(row['ci_high']-row['value'])
        ax.bar(x+(j-1)*width,vals,width,label=state,color=color,yerr=np.array([los,his]),capsize=2)
    ax.axhline(0,color='black',lw=.7); ax.set_xticks(x,labels,rotation=20,ha='right')
    ax.set_ylabel(r'contribution to $b_Q$'); ax.set_title('Exact projected instantaneous drift (seed-bootstrap 95% CI)')
    ax.legend(); ax.grid(axis='y',alpha=.18); fig.tight_layout(); fig.savefig(out/'bQ_torque_decomposition.png',dpi=220); plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(10.5,4.2)); color_map={'schooling':'#4c78a8','intermediate':'#999999','milling':'#e45756'}
    for state in states:
        m=(d['macrostate']==state)&eligible
        axs[0].scatter(d['recent_dQ_1'][m],residual[m],s=28,alpha=.75,label=state,color=color_map[state])
        axs[1].scatter(d['radial_std'][m],residual[m],s=28,alpha=.75,label=state,color=color_map[state])
    axs[0].set(xlabel=r'recent $\Delta Q$ (1 time unit)',ylabel=r'within-bin residual in $b_Q$')
    axs[1].set(xlabel='radial spread / R',ylabel=r'within-bin residual in $b_Q$')
    axs[1].legend(); [a.axhline(0,color='black',lw=.6) for a in axs]; [a.grid(alpha=.15) for a in axs]
    fig.suptitle(r'Closure screen after conditioning on fixed $(Q,M)$ cells'); fig.tight_layout(); fig.savefig(out/'closure_residuals.png',dpi=220); plt.close(fig)

    audit=dict(n_snapshots=len(seed),n_fish_rows=len(seed)*100,
        maximum_social_split_error=float(np.max(d['social_split_max_error'])),
        maximum_velocity_replay_error=float(np.max(d['velocity_replay_max_error'])),
        snapshots_with_collision_override=int(np.count_nonzero(d['n_collision_override'])),
        snapshots_with_extreme_override=int(np.count_nonzero(d['n_extreme_override'])),
        snapshots_with_sticky_override=int(np.count_nonzero(d['n_sticky_override'])),
        eligible_closure_snapshots=int(eligible.sum()),fixed_Q_edges=qedges.tolist(),fixed_M_edges=medges.tolist())
    (out/'torque_replay_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
