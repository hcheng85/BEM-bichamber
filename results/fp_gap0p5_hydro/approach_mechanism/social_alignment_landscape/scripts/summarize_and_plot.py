#!/usr/bin/env python3
"""Summaries and decisive figures for the social-alignment landscape test."""
from __future__ import annotations
import argparse,csv,sys
from pathlib import Path
import numpy as np,pandas as pd,scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from analyze_alignment_landscape import (Parameters,analytic_geometry,visible_neighbors,coefficients,
    local_torque,GRID,read_csv,wrap)

RNG=np.random.default_rng(20261023)
def ci(v,n=30000):
    v=np.asarray(v,float);v=v[np.isfinite(v)];b=np.mean(v[RNG.integers(0,len(v),(n,len(v)))],1)
    return float(np.mean(v)),float(np.percentile(b,2.5)),float(np.percentile(b,97.5))
def write(path,rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('analysis_dir',type=Path);ap.add_argument('--raw-dir',type=Path,required=True)
    ap.add_argument('--risk-set',type=Path,required=True);args=ap.parse_args();root=args.analysis_dir.resolve();tab=root/'tables';fig=root/'figures';fig.mkdir(exist_ok=True)
    l=pd.read_csv(tab/'landscape_features.csv.gz');st=pd.read_csv(tab/'frozen_environment_summary.csv')
    po=pd.read_csv(tab/'prediction_pooled.csv');bs=pd.read_csv(tab/'prediction_by_seed.csv');pc=pd.read_csv(tab/'prediction_contrasts.csv')
    outside=l[l.initially_admissible==0].copy();entrants=outside[outside.outcome_gate_entry==1]

    overview=[]
    def add(name,value,lo='',hi='',n=9120):overview.append(dict(metric=name,value=value,seed_ci_low=lo,seed_ci_high=hi,n=n))
    add('harmonic_max_abs_error_max',l.harmonic_max_abs_error.max())
    add('has_stable_orientation_fraction',(l.n_stable>0).mean())
    add('stable_orientation_in_gate_fraction',l.stable_in_gate.mean())
    add('outside_gate_genuine_barrier_fraction',outside.barrier_exists.mean(),n=len(outside))
    add('outside_gate_median_barrier_height',outside.barrier_height.median(),n=outside.barrier_height.notna().sum())
    add('median_absolute_cycle_tilt',l.cycle_tilt.abs().median())
    entry_seed=[pd.to_numeric(entrants[entrants.seed==s].entry_time).median() for s in range(20)]
    m,lo,hi=ci(entry_seed);add('median_gate_entry_time_seed_mean',m,lo,hi,len(entrants))
    tau=outside[(outside.alignment_closing_rate>0)].initial_abs_eta/outside[(outside.alignment_closing_rate>0)].alignment_closing_rate
    add('pooled_median_alignment_reorientation_time',tau.median(),n=len(tau))
    write(tab/'landscape_overview.csv',overview)

    features=('initial_abs_eta','initial_distance','alignment_closing_rate','stable_in_gate','stable_distance_to_gate','barrier_exists','barrier_height','cycle_tilt')
    contrasts=[]
    for col in features:
        vals=[]
        for seed in range(20):
            q=l[l.seed==seed];vals.append(q[q.condition=='event'][col].mean()-q[q.condition=='control'][col].mean())
        m,lo,hi=ci(vals);contrasts.append(dict(feature=col,event_minus_control=m,seed_ci_low=lo,seed_ci_high=hi,independent_seeds=20))
    write(tab/'event_control_feature_contrasts.csv',contrasts)

    outcome=[]
    for col in features:
        vals=[]
        for seed in range(20):
            q=outside[outside.seed==seed];vals.append(q[q.outcome_gate_entry==1][col].mean()-q[q.outcome_gate_entry==0][col].mean())
        m,lo,hi=ci(vals);outcome.append(dict(feature=col,entry_minus_no_entry=m,seed_ci_low=lo,seed_ci_high=hi,independent_seeds=20))
    write(tab/'entry_feature_contrasts.csv',outcome)

    summary=[]
    base=po[po.model=='baseline_eta_distance'].iloc[0];full=po[po.model=='full_landscape'].iloc[0]
    summary.extend([dict(metric='baseline_heldout_auc',value=base.auc),dict(metric='full_landscape_heldout_auc',value=full.auc),
                    dict(metric='full_minus_baseline_pooled_auc',value=full.auc-base.auc)])
    aucrow=pc[(pc.model=='full_landscape')&(pc.metric=='delta_auc_versus_baseline')].iloc[0]
    summary.append(dict(metric='mean_seed_auc_gain',value=aucrow.mean_seed_difference,seed_ci_low=aucrow.seed_bootstrap_ci_low,seed_ci_high=aucrow.seed_bootstrap_ci_high))
    for tau0 in (0.2,0.5):
        for metric in ('neighbor_jaccard','landscape_correlation','relative_landscape_rms_change','torque_sign_agreement'):
            r=st[(st.time==tau0)&(st.metric==metric)].iloc[0];summary.append(dict(metric=f'{metric}_t{tau0:g}',value=r['mean'],seed_ci_low=r.seed_ci_low,seed_ci_high=r.seed_ci_high))
    summary.append(dict(metric='scientific_outcome',value='B'))
    write(root/'summary.csv',summary)

    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    # Figure 1: verification and representative exact landscapes.
    reps=[(0,0,'event',3,'barrier, then entry'),(0,0,'event',78,'stable direction in gate'),(0,1,'control',88,'barrier, no entry')]
    risk=read_csv(args.risk_set);lookup={(int(r['seed']),int(r['pair_id']),r['condition']):r for r in risk};p=Parameters();geom=analytic_geometry(p)
    state=np.asarray(sio.loadmat(args.raw_dir/'If1e-02_Ip40_In05_R50_num100_wid50_n0_BEM.mat',variable_names=['state'],squeeze_me=True)['state'])
    fig1,axs=plt.subplots(1,4,figsize=(13,3.3),constrained_layout=True)
    axs[0].hist(np.log10(np.maximum(l.harmonic_max_abs_error,1e-16)),bins=35,color='#4d4d4d')
    axs[0].set(xlabel=r'$\log_{10}$ max grid error',ylabel='fish',title='Harmonic identity verification')
    for ax,(seed,pid,cond,fish,title) in zip(axs[1:],reps):
        rr=lookup[(seed,pid,cond)];idx=int(round(float(rr['risk_time'])/.001));st0=state[:,:,idx];vn=visible_neighbors(st0[:2],geom,p)
        c=coefficients(st0,fish,np.flatnonzero(vn[fish]));side=1 if rr['origin']=='left' else -1
        om=local_torque(GRID,c,side,p);V=np.r_[0,np.cumsum(-.5*(om[:-1]+om[1:])*np.diff(GRID))]
        row=l[(l.seed==seed)&(l.pair_id==pid)&(l.condition==cond)&(l.fish==fish)].iloc[0]
        V-=np.interp(row.initial_heading,GRID,V)
        ax.plot(GRID,V,color='#762a83');ax.axvspan(row.gate_low,row.gate_high,color='#1b9e77',alpha=.22,label='opening corridor')
        ax.axvline(row.initial_heading,color='k',ls='--',lw=1,label='initial heading')
        for a in str(row.stable_angles).split(';'):
            if a and a!='nan':ax.plot(float(a),np.interp(float(a),GRID,V),'o',color='#d95f02',ms=4)
        ax.set(xlim=(-np.pi,np.pi),xlabel=r'local heading $\theta$',ylabel=r'$V(\theta)-V(\theta_0)$',title=title)
    axs[1].legend(frameon=False,fontsize=7)
    for ext in ('png','pdf'):fig1.savefig(fig/f'01_exact_landscape_derivation.{ext}',dpi=220)
    plt.close(fig1);del state

    # Figure 2: incremental held-out prediction and barrier information.
    fig2,axs=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    order=['baseline_eta_distance','baseline_plus_drift','baseline_plus_stability','full_landscape','landscape_without_heading']
    names=['eta + distance','+ drift','+ stable orientation','full landscape','landscape, no eta']
    for s in range(15,20):
        q=bs[bs['sample']==f'seed_{s}'].set_index('model').loc[order]
        axs[0].plot(range(len(order)),q.auc,'o-',alpha=.45,lw=1)
    pp=po.set_index('model').loc[order];axs[0].plot(range(len(order)),pp.auc,'ko-',lw=2,label='pooled')
    axs[0].set_xticks(range(len(order)),names,rotation=30,ha='right');axs[0].set(ylabel='held-out AUC',title='Prospective angular entry');axs[0].legend(frameon=False)
    for outcome0,color,label in ((1,'#b2182b','entered'),(0,'#2166ac','did not enter')):
        vals=outside[outside.outcome_gate_entry==outcome0].barrier_height.dropna().clip(upper=4).sort_values().to_numpy()
        axs[1].plot(vals,np.arange(1,len(vals)+1)/len(vals),color=color,label=label)
    axs[1].set(xlabel='instantaneous barrier height (clipped at 4)',ylabel='CDF',title='First-passage barrier');axs[1].legend(frameon=False)
    groups=[]
    for stable in (0,1):
        vals=[]
        for seed in range(20):
            q=outside[(outside.seed==seed)&(outside.stable_in_gate==stable)];vals.append(q.outcome_gate_entry.mean())
        groups.append(ci(vals))
    x=np.arange(2);means=[g[0] for g in groups];lo=[g[0]-g[1] for g in groups];hi=[g[2]-g[0] for g in groups]
    axs[2].bar(x,means,color=['#999999','#1b9e77']);axs[2].errorbar(x,means,yerr=[lo,hi],fmt='none',color='k',capsize=3)
    axs[2].set_xticks(x,['stable direction\noutside gate','stable direction\ninside gate']);axs[2].set(ylabel='probability of angular entry',title='Stable orientation geometry')
    for ext in ('png','pdf'):fig2.savefig(fig/f'02_predictive_information.{ext}',dpi=220)
    plt.close(fig2)

    # Figure 3: frozen-environment closure.
    fig3,axs=plt.subplots(1,2,figsize=(8.5,3.4),constrained_layout=True)
    colors={'neighbor_jaccard':'#1b9e77','landscape_correlation':'#7570b3','torque_sign_agreement':'#d95f02'}
    labels={'neighbor_jaccard':'neighbor retention','landscape_correlation':'landscape correlation','torque_sign_agreement':'torque-sign agreement'}
    for metric in colors:
        q=st[st.metric==metric];axs[0].plot(q.time,q['mean'],marker='o',color=colors[metric],label=labels[metric]);axs[0].fill_between(q.time,q.seed_ci_low,q.seed_ci_high,color=colors[metric],alpha=.15)
    med=float(pd.to_numeric(entrants.entry_time).median());axs[0].axvline(med,color='k',ls='--',label=f'median entry {med:.2f}')
    axs[0].set(xlabel='time after A-2',ylabel='similarity to frozen environment',ylim=(0,1.04),title='Frozen-neighbor validity');axs[0].legend(frameon=False,fontsize=8)
    q=st[st.metric=='relative_landscape_rms_change'];axs[1].plot(q.time,q['mean'],marker='o',color='#e7298a');axs[1].fill_between(q.time,q.seed_ci_low,q.seed_ci_high,color='#e7298a',alpha=.18);axs[1].axvline(med,color='k',ls='--')
    axs[1].set(xlabel='time after A-2',ylabel='relative RMS landscape change',title='Landscape evolution')
    for ext in ('png','pdf'):fig3.savefig(fig/f'03_frozen_environment_test.{ext}',dpi=220)
    plt.close(fig3)

if __name__=='__main__':main()
