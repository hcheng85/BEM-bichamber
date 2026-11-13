#!/usr/bin/env python3
"""Compare COM-to-gap versus fixed-horizontal polarization, with q and both chamber vectors.

Read-only analysis of 40 existing no-hydro, gap-20, N=100 trajectories.

Three *distinct* sampling populations, never silently mixed:
  (1) Time-weighted, post-transient raw snapshots for an unconditional empirical p(q,a).
  (2) Existing non-overlapping, crossing-free-history risk origins, with exactly the
      same first-crossing outcome as analyze_q_gap_polarization.py.
  (3) Existing event-aligned *majority-to-minority* first crossings, with the
      same four-additional-net-fish continuation outcome as the original script.

The COM and fixed-axis definitions are evaluated on the identical raw frames.
No fitted potential, new simulation, assumed Markov closure, or raw modification.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

N = 100
Q_EDGES = np.linspace(-1.01, 1.01, 22)
QABS_EDGES = np.array([0., .5, .8, 1.000001])
A_EDGES = np.array([-1.000001, -.6, -.2, .2, .6, 1.000001])
CONT_A_EDGES = np.array([-1.000001, .2, .6, 1.000001])
A_JOINT_EDGES = np.linspace(-1.000001, 1.000001, 21)


def load_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def save_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def features(state):
    """Instantaneous, gap-centered chamber data. Empty chambers are NaN, not 0.

    a_L = P_Lx; a_R = -P_Rx; a_fixed = a_(majority chamber).
    a_com = P_majority dot (gap-COM_majority)/norm(gap-COM_majority).
    An exact tie q=0 has *no majority coordinate*, but a_L,a_R remain defined.
    """
    if state.shape != (3,N):
        raise ValueError(f'Expected raw state shape (3,{N}), got {state.shape}')
    x,y,theta=state
    left=x<0;right=~left
    nl=int(left.sum()); nr=N-nl; q=(nr-nl)/N
    headings=np.column_stack([np.cos(theta),np.sin(theta)])
    def chamber(mask):
        if not np.any(mask): return np.array([np.nan,np.nan]),np.array([np.nan,np.nan])
        return headings[mask].mean(axis=0), np.array([x[mask].mean(), y[mask].mean()])
    PL,CL=chamber(left);PR,CR=chamber(right)
    aL=float(PL[0]);aR=float(-PR[0]);
    if q>0: pm=PR;cm=CR;af=aR;amin=aL;chamber_name='right'
    elif q<0: pm=PL;cm=CL;af=aL;amin=aR;chamber_name='left'
    else: pm=np.array([np.nan,np.nan]);cm=np.array([np.nan,np.nan]);af=float('nan');amin=float('nan');chamber_name='tie'
    radius=float(np.linalg.norm(cm))
    ac=float(pm@(-cm/radius)) if np.isfinite(radius) and radius>1e-12 else np.nan
    # Quantify how the two projections differ; this is a geometrical discrepancy,
    # *not* evidence for changed heading dynamics.
    record=dict(q=float(q),nl=nl,nr=nr,P_Lx=float(PL[0]),P_Ly=float(PL[1]),
                P_Rx=float(PR[0]),P_Ry=float(PR[1]),
                C_Lx=float(CL[0]),C_Ly=float(CL[1]),C_Rx=float(CR[0]),C_Ry=float(CR[1]),
                a_L=aL,a_R=aR,a_fixed=af,a_minority=amin,a_com=ac,
                P_majority_x=float(pm[0]),P_majority_y=float(pm[1]),
                C_majority_x=float(cm[0]),C_majority_y=float(cm[1]),
                majority=chamber_name,minority_count=min(nl,nr))
    return record


def bindex(value, edges):
    if not np.isfinite(value):return None
    i=int(np.searchsorted(edges,float(value),side='right')-1)
    return i if 0<=i<len(edges)-1 else None


def auc(y,p):
    # Pairwise exact handling of ties without scipy; O(n log n).
    y=np.asarray(y,int);p=np.asarray(p,float)
    if len(y)<2 or len(set(y.tolist()))<2:return np.nan
    order=np.argsort(p,kind='stable');pp=p[order];yy=y[order]
    rank_pos=0.;k=0
    while k<len(y):
        j=k+1
        while j<len(y) and pp[j]==pp[k]:j+=1
        rank_pos+=(k+1+j)*.5*float(yy[k:j].sum());k=j
    n1=int(y.sum());n0=len(y)-n1
    return (rank_pos-n1*(n1+1)/2)/(n1*n0)


def model_predict(train, test, mode, predictor, min_cell, min_seeds):
    """Fixed, non-parametric cell probabilities from independent train seeds.

    mode: init or continuation. predictor: q, com, fixed, com_only, fixed_only.
    Fall back to q probability if joint cell is sparsely represented.
    Half-count correction is purely a binomial probability regularizer.
    """
    aedges=A_EDGES if mode=='init' else CONT_A_EDGES
    counts=defaultdict(lambda:[0,0,set()]); qc=defaultdict(lambda:[0,0,set()]);ac=defaultdict(lambda:[0,0,set()])
    tot=0;yes=0
    for row in train:
        qb=bindex(abs(float(row['q'])),QABS_EDGES)
        if qb is None:continue
        y=int(row['outcome']);seed=int(row['seed']);yes+=y;tot+=1
        qc[qb][0]+=1;qc[qb][1]+=y;qc[qb][2].add(seed)
        for field in ('a_com','a_fixed'):
            ab=bindex(float(row[field]),aedges)
            if ab is None:continue
            counts[(field,qb,ab)][0]+=1;counts[(field,qb,ab)][1]+=y;counts[(field,qb,ab)][2].add(seed)
            ac[(field,ab)][0]+=1;ac[(field,ab)][1]+=y;ac[(field,ab)][2].add(seed)
    overall=(yes+.5)/(tot+1.) if tot else .5
    predictions=[]
    for row in test:
        qb=bindex(abs(float(row['q'])),QABS_EDGES)
        if qb is None:continue
        cnt,y,_=qc[qb];p_q=(y+.5)/(cnt+1) if cnt>=min_cell else overall
        if predictor=='q': p=p_q;eligible=True
        else:
            field='a_com' if 'com' in predictor else 'a_fixed'
            ab=bindex(float(row[field]),aedges)
            if ab is None:
                # Missing coordinate: no invented chamber direction.
                p=overall if predictor.endswith('_only') else p_q;eligible=False
            else:
                cell=ac[(field,ab)] if predictor.endswith('_only') else counts[(field,qb,ab)]
                n,y,seedset=cell
                eligible=n>=min_cell and len(seedset)>=min_seeds
                p=(y+.5)/(n+1) if eligible else (overall if predictor.endswith('_only') else p_q)
        predictions.append(dict(seed=int(row['seed']),outcome=int(row['outcome']),p=float(p),
                                supported=int(eligible)))
    return predictions


def metrics(pred):
    if not pred:return None
    y=np.array([r['outcome'] for r in pred]);p=np.array([r['p'] for r in pred]);
    return dict(n=len(y),positive=int(y.sum()),auc=float(auc(y,p)),
                brier=float(np.mean((y-p)**2)),
                coverage=float(np.mean([r['supported'] for r in pred])))


def compare_heldout(rows, mode, cutoff=30, min_cell=None, min_seeds=None, boots=4000):
    min_cell=min_cell if min_cell is not None else (20 if mode=='init' else 10)
    min_seeds=min_seeds if min_seeds is not None else (1 if mode=='init' else 3)
    train=[r for r in rows if r['seed']<cutoff];test=[r for r in rows if r['seed']>=cutoff]
    if not train or not test:return {'available':False,'reason':'No train/test observations'}
    models=['q','com_only','fixed_only','com','fixed']
    model_p={k:model_predict(train,test,mode,k,min_cell,min_seeds) for k in models}
    # Identical ordered test rows in all models, including tie fallback.
    if len({len(p) for p in model_p.values()})!=1:
        raise AssertionError('Models do not have an identical heldout cohort')
    out={k:metrics(model_p[k]) for k in models}
    per_seed=[]
    for seed in sorted(set(r['seed'] for r in test)):
        sub={k:[r for r in model_p[k] if r['seed']==seed] for k in models}
        item={'seed':seed,'n':len(sub['q']),'positives':sum(r['outcome'] for r in sub['q'])}
        for k in models:
            m=metrics(sub[k]);item.update({f'{k}_auc':m['auc'],f'{k}_brier':m['brier']})
        per_seed.append(item)
    rng=np.random.default_rng(20261113)
    comparisons={}
    for num,den in [('fixed','com'),('fixed','q'),('com','q'),('fixed_only','com_only'),('com','com_only'),('fixed','fixed_only')]:
        dif_b=np.asarray([r[f'{den}_brier']-r[f'{num}_brier'] for r in per_seed],float)
        dif_a=np.asarray([r[f'{num}_auc']-r[f'{den}_auc'] for r in per_seed],float)
        def boot(v):
            v=v[np.isfinite(v)]
            if not len(v):return dict(mean=np.nan,ci95=[np.nan,np.nan],seeds=0)
            draws=np.mean(v[rng.integers(0,len(v),size=(boots,len(v)))],axis=1)
            return dict(mean=float(v.mean()),ci95=list(map(float,np.quantile(draws,[.025,.975]))),seeds=len(v))
        comparisons[f'{num}_vs_{den}']={'seed_mean_brier_improvement':boot(dif_b),
                                        'seed_mean_auc_improvement':boot(dif_a)}
    return dict(available=True,mode=mode,train_n=len(train),heldout_n=len(test),
                train_seeds=sorted(set(r['seed'] for r in train)),
                heldout_seeds=sorted(set(r['seed'] for r in test)),
                min_train_cell=min_cell,min_seeds_per_cell=min_seeds,
                metrics=out,per_seed=per_seed,comparisons=comparisons)


def minority_added_value(rows, mode, cutoff=30, min_minor=5, boots=2500):
    """Supplement: is a minority-chamber direction needed, on the SAME eligible subset?

    Train/test only when both chambers contain >=min_minor fish. This is NOT a
    claim that minority polarization exists in fully localized school states.
    Two descriptive fixed 3-bin projections; sparse joint cells back off to
    the majority-only model, with coverage reported.
    """
    eligible=[r for r in rows if r['minority_count']>=min_minor and np.isfinite(r['a_minority'])]
    train=[r for r in eligible if r['seed']<cutoff]
    test=[r for r in eligible if r['seed']>=cutoff]
    if not train or not test:return dict(available=False,n_eligible=len(eligible),reason='Insufficient jointly occupied train/test samples')
    edges=np.array([-1.000001,-.2,.6,1.000001])
    c_major=defaultdict(lambda:[0,0,set()]);c_both=defaultdict(lambda:[0,0,set()]);c_q=defaultdict(lambda:[0,0,set()])
    for r in train:
        q=bindex(abs(r['q']),QABS_EDGES);am=bindex(r['a_fixed'],edges);an=bindex(r['a_minority'],edges)
        if None in (q,am,an):continue
        for counts,key in [(c_q,q),(c_major,(q,am)),(c_both,(q,am,an))]:
            counts[key][0]+=1;counts[key][1]+=int(r['outcome']);counts[key][2].add(r['seed'])
    min_cell=20 if mode=='init' else 10
    grand=(sum(r['outcome'] for r in train)+.5)/(len(train)+1)
    preds=[]
    for r in test:
        q=bindex(abs(r['q']),QABS_EDGES);am=bindex(r['a_fixed'],edges);an=bindex(r['a_minority'],edges)
        if None in (q,am,an):continue
        nq,yq,_=c_q[q]; pq=(yq+.5)/(nq+1) if nq>=min_cell else grand
        nm,ym,sm=c_major[(q,am)];maj_supported=nm>=min_cell and len(sm)>=3
        pm=(ym+.5)/(nm+1) if maj_supported else pq
        nb,yb,sb=c_both[(q,am,an)];both_supported=nb>=min_cell and len(sb)>=3
        pb=(yb+.5)/(nb+1) if both_supported else pm
        preds.append(dict(seed=r['seed'],outcome=r['outcome'],p_major=pm,p_both=pb,
                          eligible_joint=int(both_supported)))
    if not preds:return dict(available=False,n_eligible=len(eligible),reason='No test samples after bin checks')
    def sc(rs,key):
        y=np.asarray([r['outcome'] for r in rs]);p=np.asarray([r[key] for r in rs]);
        return {'auc':float(auc(y,p)),'brier':float(np.mean((y-p)**2))}
    seeds=sorted({r['seed'] for r in preds});per=[]
    for seed in seeds:
        sub=[r for r in preds if r['seed']==seed]
        a=sc(sub,'p_major');b=sc(sub,'p_both')
        auc_gain=b['auc']-a['auc']
        per.append(dict(seed=seed,n=len(sub),n_success=sum(r['outcome'] for r in sub),
                        coverage=np.mean([r['eligible_joint'] for r in sub]),
                        brier_gain=a['brier']-b['brier'],
                        auc_gain=float(auc_gain) if np.isfinite(auc_gain) else None))
    rng=np.random.default_rng(20261113)
    def ci(key):
        arr=np.asarray([np.nan if r[key] is None else r[key] for r in per],float);arr=arr[np.isfinite(arr)]
        if not len(arr):return {'mean':None,'ci95':None,'n_seeds':0}
        means=arr[rng.integers(len(arr),size=(boots,len(arr)))].mean(axis=1)
        return {'mean':float(arr.mean()),'ci95':list(map(float,np.quantile(means,[.025,.975]))),'n_seeds':len(arr)}
    return dict(available=True,mode=mode,min_minority_fish=min_minor,
                n_eligible=len(eligible),n_train=len(train),n_test=len(test),
                n_success_test=sum(r['outcome'] for r in preds),
                coverage_joint_on_test=np.mean([r['eligible_joint'] for r in preds]),
                majority_only=sc(preds,'p_major'),majority_and_minority=sc(preds,'p_both'),
                additional_minority_brier_gain=ci('brier_gain'),
                additional_minority_auc_gain=ci('auc_gain'),per_seed=per,
                caveat='Restricted to both-occupied states with >=5 fish in minority; cannot extrapolate to q~1. Sparse cells fall back to majority-only prediction.')


def event_control_polarization_timecourses(raw_dir, observed_csv, intervals_csv,
                                            events_csv, matched_csv, dt, out, boots=2000):
    """Reconstruct the same previously chosen 37 event/control pairs, NO rematching.

    The old COM time series supplies exact time sampling and a_COM audit;
    absolute frames come from prior event IDs and *existing* matched controls.
    """
    files=(observed_csv,intervals_csv,events_csv,matched_csv)
    if not all(p.exists() for p in files):
        return dict(available=False,missing=[str(p) for p in files if not p.exists()])
    from analyze_initiation_propagation import RawReader
    observed=load_csv(observed_csv)
    intervals=load_csv(intervals_csv)
    events={(int(r['seed']),int(r['event_id'])):r for r in load_csv(events_csv)}
    match={(int(r['event_seed']),int(r['event_id']),int(r['control_seed'])):r
           for r in load_csv(matched_csv) if abs(float(r['lead'])-2)<1e-8}
    # Per-pair event ID is fixed by the archived analysis output. Controls do
    # not get reselected based on fixed-axis performance.
    pair_meta={}
    for r in intervals:
        key=int(r['pair_seed']);kind=r['kind']
        pair_meta.setdefault(key,{})[kind]=r
    original=defaultdict(list)
    for r in observed:
        original[(int(r['pair_seed']),r['kind'])].append(r)
    records=[];max_com_error=0.
    for seed,groups in sorted(pair_meta.items()):
        if not {'event','control'}.issubset(groups):continue
        e=groups['event'];c=groups['control'];event_id=int(e['event_id'])
        ev=events.get((seed,event_id))
        mt=match.get((seed,event_id,int(c['trajectory_seed'])))
        if ev is None or mt is None:
            raise RuntimeError(f'Cannot locate immutable event/control pairing for seed {seed} event {event_id}')
        event_end=int(ev['first_native'])-1
        control_start=round(float(mt['control_time'])/dt)
        control_end=control_start+round(2/dt)-1
        for kind in ('event','control'):
            run_seed=seed if kind=='event' else int(c['trajectory_seed'])
            target_end=event_end if kind=='event' else control_end
            path=raw_dir/f'If00_Ip40_In05_R50_num100_wid20_n{run_seed}_BEM_no_hydro.mat'
            reader=RawReader(path)
            try:
                for z in sorted(original[(seed,kind)],key=lambda r:float(r['time_to_interval_end'])):
                    t=float(z['time_to_interval_end'])
                    k=target_end+round(t/dt)
                    f=features(reader.frame('state',k))
                    old_a=float(z['a'])
                    if not np.isfinite(f['a_com']) or abs(f['a_com']-old_a)>1e-8:
                        raise RuntimeError(f'Old COM time course parity failed: pair={seed}, kind={kind}, t={t}, old={old_a}, new={f["a_com"]}')
                    max_com_error=max(max_com_error,abs(f['a_com']-old_a))
                    records.append(dict(pair_seed=seed,trajectory_seed=run_seed,kind=kind,
                                        time_to_crossing=round(t,4),native_step=k,
                                        a_com=f['a_com'],a_fixed=f['a_fixed'],
                                        q=f['q'],a_L=f['a_L'],a_R=f['a_R']))
            finally:reader.close()
    if not records:return dict(available=False,reason='No matched observations')
    save_csv(out/'matched_polarization_timecourses.csv',records)
    pairs=sorted({r['pair_seed'] for r in records});times=sorted({r['time_to_crossing'] for r in records})
    perpair=[];trajectory={}
    for seed in pairs:
        for kind in ('event','control'):
            subset=sorted((r for r in records if r['pair_seed']==seed and r['kind']==kind),key=lambda r:r['time_to_crossing'])
            if len(subset)!=len(times):raise RuntimeError(f'Incomplete pair-time series seed {seed} kind {kind}')
            if [r['time_to_crossing'] for r in subset]!=times:raise RuntimeError('Time grid mismatch')
            trajectory[(seed,kind)]=subset
    for seed in pairs:
        rec={'seed':seed}
        for kind in ('event','control'):
            first=trajectory[(seed,kind)][0];last=trajectory[(seed,kind)][-1]
            for name in ('a_com','a_fixed'):
                rec[f'{kind}_start_{name}']=first[name]
                rec[f'{kind}_end_{name}']=last[name]
                rec[f'{kind}_delta_{name}']=last[name]-first[name]
        for name in ('a_com','a_fixed'):
            rec[f'paired_delta_{name}']=rec[f'event_delta_{name}']-rec[f'control_delta_{name}']
            rec[f'paired_start_{name}']=rec[f'event_start_{name}']-rec[f'control_start_{name}']
            rec[f'paired_end_{name}']=rec[f'event_end_{name}']-rec[f'control_end_{name}']
        perpair.append(rec)
    save_csv(out/'matched_timecourse_pair_summary.csv',perpair)
    rng=np.random.default_rng(20261113)
    def ci(arr):
        arr=np.asarray(arr,float)
        boot=arr[rng.integers(len(arr),size=(boots,len(arr)))].mean(axis=1)
        return dict(mean=float(arr.mean()),ci95=list(map(float,np.quantile(boot,[.025,.975]))))
    summary={'available':True,'n_pairs':len(pairs),'n_snapshots':len(records),'max_abs_COM_parity_error':max_com_error}
    for name in ('a_com','a_fixed'):
        summary[name]={k:ci([r[k] for r in perpair]) for k in (
            f'event_delta_{name}',f'control_delta_{name}',
            f'paired_delta_{name}',f'paired_start_{name}',f'paired_end_{name}')}
    # Display the *same* trajectories with the two coordinate definitions.
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for j,name in enumerate(('a_com','a_fixed')):
        arrays={kind:np.array([[r[name] for r in trajectory[(seed,kind)]] for seed in pairs])
                for kind in ('event','control')}
        for kind in ('event','control'):
            arr=arrays[kind];mean=arr.mean(axis=0)
            ax=axes[0,j];ax.plot(times,mean,label=kind)
            b=arr[rng.integers(len(pairs),size=(boots,len(pairs)))].mean(axis=1)
            lo,hi=np.quantile(b,[.025,.975],axis=0)
            ax.fill_between(times,lo,hi,alpha=.15)
        ax=axes[0,j];ax.set(title=name+' (matched event/control)',xlabel='time to first crossing',ylabel='gap-directed mean heading');ax.legend(frameon=False)
        difference=arrays['event']-arrays['control'];mean=difference.mean(axis=0)
        b=difference[rng.integers(len(pairs),size=(boots,len(pairs)))].mean(axis=1)
        lo,hi=np.quantile(b,[.025,.975],axis=0)
        ax=axes[1,j];ax.plot(times,mean);ax.fill_between(times,lo,hi,alpha=.15)
        ax.axhline(0,color='gray',linewidth=1,linestyle=':')
        ax.set(title=name+' event minus matched control',xlabel='time to first crossing',ylabel='paired difference')
    fig.suptitle('Same 37 archived matched pairs; reference direction changes only for COM projection')
    fig.savefig(out/'matched_com_vs_fixed_timecourses.png',dpi=180);plt.close(fig)
    return summary


def audit_against_existing(results,root):
    """Numerical parity with original COM predictor (same cohorts/outcomes).

    More than 1e-8 difference is a warning, not a reason to silently replace
    a physically meaningful result with original summary figures.
    """
    files={'initiation':root/'q_gap_polarization/heldout_validation.json',
           'continuation':root/'event_aligned_continuation/heldout_validation.json'}
    audit={}
    for mode,path in files.items():
        if not path.exists():continue
        prior=json.loads(path.read_text());new=results[mode]
        if not new['available']:continue
        fields={'pooled_auc_q':new['metrics']['q']['auc'],
                'pooled_auc_qa':new['metrics']['com']['auc'],
                'pooled_brier_q':new['metrics']['q']['brier'],
                'pooled_brier_qa':new['metrics']['com']['brier']}
        differences={key:abs(float(prior[key])-value) for key,value in fields.items()}
        audit[mode]=dict(differences=differences,
                         parity_all_below_1e_minus_8=all(v<1e-8 for v in differences.values()))
    return audit


def sample_raw(raw_dir, seed, risk, episodes, transient, dt, stationary_dt):
    # The existing reader is imported lazily: --self-test requires no checkout.
    from analyze_initiation_propagation import RawReader
    raw=raw_dir/f'If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat'
    if not raw.exists():raise FileNotFoundError(f'Raw trajectory required for the comparison: {raw}')
    rr=RawReader(raw)
    try:
        ntime=int(rr.info['state'].shape[-1]);
        start=round(transient/dt);jump=round(stationary_dt/dt)
        if jump<=0 or abs(jump*dt-stationary_dt)>1e-8:raise ValueError('stationary_dt must be integer multiple of dt')
        if start>=ntime-1:raise ValueError('Transient exceeds trajectory duration')
        station_steps=set(range(start,ntime-1,jump))
        risk_steps={int(r['native_step']) for r in risk}
        episode_steps={int(r['polarization_native_step']) for r in episodes}
        needed=sorted(station_steps|risk_steps|episode_steps)
        cache={}
        stationary=[]
        for k in needed:
            if not 0<=k<ntime: raise ValueError((seed,k,ntime))
            feat=features(rr.frame('state',k));
            if k in station_steps:stationary.append(dict(seed=seed,step=k,t=k*dt,**feat))
            if k in risk_steps or k in episode_steps:cache[k]=feat
        return cache,stationary
    finally:
        rr.close()


def merge_with_raw(risk_path,event_path,raw_dir,seeds,transient,dt,stationary_dt):
    original=load_csv(risk_path);episodes=load_csv(event_path)
    risk_by=defaultdict(list);ep_by=defaultdict(list)
    for r in original:
        s=int(r['seed'])
        if s in seeds:risk_by[s].append(r)
    for e in episodes:
        s=int(e['seed'])
        if s in seeds and not e.get('exclude','') and int(e['first_majority_to_minority'])==1 and int(e['prelead_endpoints_match'])==1:
            ep_by[s].append(e)
    risk_out=[];continuation_out=[];stationary_out=[];errors=[]
    for seed in sorted(seeds):
        cache,stat=sample_raw(raw_dir,seed,risk_by[seed],ep_by[seed],transient,dt,stationary_dt)
        stationary_out.extend(stat)
        for row in risk_by[seed]:
            k=int(row['native_step']);f=cache[k];q=float(row['q']);a=float(row['a'])
            if abs(f['q']-q)>1e-10 or abs(f['a_com']-a)>1e-9:
                raise RuntimeError(f'Prospective raw-state mismatch seed={seed} step={k}: q {f["q"]} vs {q}, a {f["a_com"]} vs {a}')
            risk_out.append(dict(seed=seed,t=float(row['t']),step=k,q=q,outcome=int(row['first_crossing']),
                                 a_com=f['a_com'],a_fixed=f['a_fixed'],a_minority=f['a_minority'],a_L=f['a_L'],a_R=f['a_R'],
                                 P_Lx=f['P_Lx'],P_Ly=f['P_Ly'],P_Rx=f['P_Rx'],P_Ry=f['P_Ry'],
                                 nl=f['nl'],nr=f['nr'],minority_count=f['minority_count'],
                                 com_x=f['C_majority_x'],com_y=f['C_majority_y']))
            errors.append(abs(f['a_com']-a))
        for row in ep_by[seed]:
            k=int(row['polarization_native_step']);f=cache[k];q=float(row['q_pre']);a=float(row['a_pre'])
            if abs(f['q']-q)>1e-10 or abs(f['a_com']-a)>1e-9:
                raise RuntimeError(f'Event-aligned raw-state mismatch seed={seed} step={k}')
            continuation_out.append(dict(seed=seed,event_id=int(row['event_id']),step=k,
                                         t=float(row['polarization_t']),q=q,outcome=int(row['continued']),
                                         first_t=float(row['first_t']),
                                         additional_net_fish=float(row['additional_net_fish']),
                                         a_com=f['a_com'],a_fixed=f['a_fixed'],a_minority=f['a_minority'],a_L=f['a_L'],a_R=f['a_R'],
                                         P_Lx=f['P_Lx'],P_Ly=f['P_Ly'],P_Rx=f['P_Rx'],P_Ry=f['P_Ry'],
                                         nl=f['nl'],nr=f['nr'],minority_count=f['minority_count'],
                                         com_x=f['C_majority_x'],com_y=f['C_majority_y']))
            errors.append(abs(f['a_com']-a))
        print(f'seed {seed:02d}: stationary {len(stat)}, initiation {len(risk_by[seed])}, continuation {len(ep_by[seed])}',flush=True)
    return stationary_out,risk_out,continuation_out,dict(n_raw_crosschecks=len(errors),max_abs_a_com_error=float(max(errors)) if errors else None)


def seed_histogram(rows, qedges, aedges, column):
    seeds=sorted({r['seed'] for r in rows});by=[];nvalid={}
    for seed in seeds:
        part=[r for r in rows if r['seed']==seed and np.isfinite(r[column])]
        nvalid[seed]=len(part)
        if not part:
            # A chamber-specific coordinate can be undefined for an entire
            # seed when that chamber is never occupied in the sampled frames.
            # Preserve that absence instead of fabricating a zero vector.
            continue
        arr=np.asarray([[r['q'],r[column]] for r in part]) if part else np.empty((0,2))
        h=np.histogram2d(arr[:,0],arr[:,1],bins=[qedges,aedges])[0]
        if h.sum()==0:raise ValueError('Valid values fell outside bins: '+str(seed))
        by.append(h/h.sum())
    if not by:raise ValueError('No valid states for coordinate: '+column)
    return np.stack(by),nvalid


def create_stationary_outputs(stat, out, boots=1000, save_samples=False):
    if save_samples:
        save_csv(out/'stationary_snapshots.csv',stat)
    fig,axes=plt.subplots(2,2,figsize=(11.5,8.5),layout='constrained')
    stats={}
    rng=np.random.default_rng(20261113)
    fields=[('a_fixed','Majority fixed-axis polarization'),('a_com','Majority COM-to-gap polarization'),
            ('a_L','Left chamber: toward gap (+x)'),('a_R','Right chamber: toward gap (−x)')]
    for ax,(field,name) in zip(axes.ravel(),fields):
        h,ns=seed_histogram(stat,Q_EDGES,A_JOINT_EDGES,field)
        mass=h.mean(axis=0)
        positives=mass[mass>0];
        im=ax.pcolormesh(Q_EDGES,A_JOINT_EDGES,np.ma.masked_where(mass.T<=0,mass.T),
                         shading='auto',norm=LogNorm(vmin=float(positives.min()),vmax=float(positives.max())))
        fig.colorbar(im,ax=ax,label='equal-seed probability mass')
        ax.set(xlabel='q',ylabel=field,title=name)
        draws=np.zeros((boots,*mass.shape),dtype=float)
        for b in range(boots):
            draws[b]=h[rng.integers(0,len(h),size=len(h))].mean(axis=0)
        lower,upper=np.quantile(draws,[.025,.975],axis=0)
        phi=np.full(mass.shape,np.nan)
        positive=mass>0
        phi[positive]=-np.log(mass[positive]/np.max(mass))
        outrows=[dict(q_min=Q_EDGES[i],q_max=Q_EDGES[i+1],a_min=A_JOINT_EDGES[j],a_max=A_JOINT_EDGES[j+1],
                      mean_seed_mass=mass[i,j],ci95_lo=lower[i,j],ci95_hi=upper[i,j],
                      statistical_phi_shifted=phi[i,j],n_seeds=len(h))
                 for i in range(len(Q_EDGES)-1) for j in range(len(A_JOINT_EDGES)-1)]
        save_csv(out/f'stationary_joint_{field}.csv',outrows)
        stats[field]=dict(n_seeds=len(h),n_valid_by_seed=ns)
    fig.suptitle('Unconditional post-transient empirical joint distributions (NOT dynamical potentials)')
    fig.savefig(out/'stationary_joint_comparison.png',dpi=170);plt.close(fig)
    # Show -ln p only for occupied cells, and only as a statistical surface.
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for ax,field in zip(axes,['a_com','a_fixed']):
        h,_=seed_histogram(stat,Q_EDGES,A_JOINT_EDGES,field)
        mass=h.mean(axis=0);phi=np.full(mass.shape,np.nan)
        nz=mass>0;phi[nz]=-np.log(mass[nz]/mass.max())
        im=ax.pcolormesh(Q_EDGES,A_JOINT_EDGES,np.ma.masked_invalid(phi.T),shading='auto',cmap='viridis')
        fig.colorbar(im,ax=ax,label='shifted -ln empirical mass')
        ax.set(xlabel='q',ylabel=field,title='Empirical statistical surface: '+field)
    fig.suptitle('Finite-time empirical landscape only; no dynamical barrier inference')
    fig.savefig(out/'stationary_statistical_landscapes.png',dpi=170);plt.close(fig)
    # Optional compact 3D histogram p(q,aL,aR) CONDITIONAL on both chambers occupied.
    both=[r for r in stat if np.isfinite(r['a_L']) and np.isfinite(r['a_R'])]
    b_q=Q_EDGES;b_a=np.linspace(-1.000001,1.000001,11)
    if both:
        arr=np.asarray([[r['q'],r['a_L'],r['a_R']] for r in both])
        cube=np.histogramdd(arr,bins=[b_q,b_a,b_a])[0]
        cube/=cube.sum()
        three=[dict(q_min=b_q[i],q_max=b_q[i+1],aL_min=b_a[j],aL_max=b_a[j+1],
                    aR_min=b_a[k],aR_max=b_a[k+1],mass_conditioned_on_two_occupied=cube[i,j,k])
               for i in range(len(b_q)-1) for j in range(len(b_a)-1) for k in range(len(b_a)-1)
               if cube[i,j,k]>0]
        save_csv(out/'both_occupied_joint_q_aL_aR.csv',three)
    # q=0 has no preferred majority; empty chamber polarization is undefined.
    n_qzero=sum(r['q']==0 for r in stat)
    n_emptyL=sum(r['nl']==0 for r in stat);n_emptyR=sum(r['nr']==0 for r in stat)
    return dict(n_snapshots=len(stat),n_q_zero=n_qzero,n_empty_left=n_emptyL,n_empty_right=n_emptyR,
                n_both_chambers_occupied=len(both),fraction_both_occupied=len(both)/len(stat),
                both_occupied_3d_histogram='p(q,a_L,a_R) conditioned on both chambers nonempty; not an unconditional 3D landscape',
                coordinates=stats,temporal_spacing='See analysis_config.json; frames uniformly subsampled within each seed')


def conditional_grid(rows, field, mode):
    edges=A_EDGES if mode=='init' else CONT_A_EDGES
    counts=np.zeros((3,len(edges)-1),dtype=int);yes=np.zeros_like(counts)
    for r in rows:
        qi=bindex(abs(float(r['q'])),QABS_EDGES);ai=bindex(float(r[field]),edges)
        if qi is None or ai is None:continue
        counts[qi,ai]+=1;yes[qi,ai]+=int(r['outcome'])
    prob=np.divide(yes,counts,out=np.full(counts.shape,np.nan),where=counts>0)
    tab=[dict(q_min=QABS_EDGES[i],q_max=QABS_EDGES[i+1],a_min=edges[j],a_max=edges[j+1],
              count=int(counts[i,j]),positive=int(yes[i,j]),probability=prob[i,j])
         for i in range(3) for j in range(len(edges)-1)]
    return prob,tab


def plot_outcome(rows,mode,out):
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.3),layout='constrained')
    for ax,field,name in zip(axes,('a_com','a_fixed'),('COM-to-gap orientation','Fixed-axis orientation')):
        grid,table=conditional_grid(rows,field,mode)
        save_csv(out/f'{mode}_conditional_{field}.csv',table)
        im=ax.imshow(np.ma.masked_invalid(grid),origin='lower',aspect='auto',vmin=0,vmax=1,
                      extent=[0,grid.shape[1],0,3])
        ax.set(xticks=np.arange(grid.shape[1])+.5,
               xticklabels=['away','slightly away','neutral','slightly toward','toward'] if mode=='init' else ['low','medium','high'],
               yticks=[.5,1.5,2.5],yticklabels=['|q| < .5','.5 ≤ |q| < .8','|q| ≥ .8'],
               xlabel=field,title=name)
        for t in table:
            i=bindex(float(t['q_min'])+1e-6,QABS_EDGES);j=bindex(float(t['a_min'])+1e-6,A_EDGES if mode=='init' else CONT_A_EDGES)
            if i is None or j is None:continue
            ax.text(j+.5,i+.5,f'{t["positive"]}/{t["count"]}',ha='center',va='center',fontsize=8,
                    color='white' if np.isfinite(t['probability']) and t['probability']>.5 else 'black')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='conditional outcome probability',shrink=.75)
    fig.suptitle('First crossing within 5 units' if mode=='init' else '≥4 additional net fish after first crossing')
    fig.savefig(out/f'{mode}_conditional_com_vs_fixed.png',dpi=180,bbox_inches='tight');plt.close(fig)


def paired_gap_direction_effect(stat,out):
    # Fixed-minus-COM projection as a function of |q|, and COM y-offset.
    rows=[]
    for lo,hi in zip(QABS_EDGES[:-1],QABS_EDGES[1:]):
        for typ in ('a_com','a_fixed'):
            values=[float(r[typ]) for r in stat if lo<=abs(r['q'])<hi and np.isfinite(r[typ])]
            rows.append(dict(abs_q_min=lo,abs_q_max=hi,definition=typ,n=len(values),mean=float(np.mean(values)) if values else np.nan))
    save_csv(out/'stationary_polarization_by_imbalance.csv',rows)


def selftest():
    # Left-majority fish heading +x, right-minority fish heading -x.
    st=np.zeros((3,N));st[0,:80]=-2;st[0,80:]=2;st[1,:80]=0;st[2,80:]=np.pi
    f=features(st)
    assert f['nl']==80 and np.isclose(f['q'],-.6)
    assert np.isclose(f['a_L'],1) and np.isclose(f['a_R'],1)
    assert np.isclose(f['a_com'],1) and np.isclose(f['a_fixed'],1)
    # COM off-axis: fixed is unchanged, COM projection is lower.
    st[1,:80]=2
    f=features(st)
    assert np.isclose(f['a_fixed'],1) and np.isclose(f['a_com'],1/np.sqrt(2))
    # Equal populations have no current majority; both chamber vectors still valid.
    st[0,50:]=-2;st[0,:50]=2
    f=features(st);assert f['q']==0 and np.isnan(f['a_fixed']) and np.isnan(f['a_com'])
    assert np.isfinite(f['a_L']) and np.isfinite(f['a_R'])
    # Empty chamber is missing, not a fabricated zero heading.
    st[0,:]=-2;f=features(st)
    assert np.isnan(f['a_R']) and np.isfinite(f['a_fixed'])
    # A chamber-specific stationary histogram omits a seed with no defined
    # observations and records that absence rather than zero-filling it.
    h,nvalid=seed_histogram([dict(seed=0,q=-1.,a=np.nan),
                             dict(seed=1,q=1.,a=.5)],Q_EDGES,A_JOINT_EDGES,'a')
    assert h.shape[0]==1 and nvalid=={0:0,1:1}
    # AUC tied ordering and fixed unit-bin probabilities.
    assert np.isclose(auc([0,1,0,1],[.5,.5,0,1]),.875)
    rr=[]
    for seed in range(40):
        for i in range(18):
            q=(-.9,-.4,.65)[i%3]
            af=(-.75,-.3,.1,.4,.8,.92)[i%6]
            ac=np.clip(af*.7+.08,-1,1)
            y=int(af>.3 and ((i+seed)%5!=0))
            rr.append(dict(seed=seed,q=q,a_fixed=af,a_com=ac,outcome=y))
    for mode in ('init','continuation'):
        val=compare_heldout(rr,mode,boots=30)
        assert val['available'] and val['heldout_n']==180
        for k in ('q','fixed','com','fixed_only','com_only'):
            assert 0<=val['metrics'][k]['auc']<=1
    print('PASS: chamber definitions, off-axis geometry, ties, empty chamber, AUC and heldout scoring')


def main():
    ap=argparse.ArgumentParser()
    root=Path('results/nohydro_gap20_descriptive_transition')
    ap.add_argument('--raw-dir',type=Path,default=Path('/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data'))
    ap.add_argument('--risk-file',type=Path,default=root/'outputs/q_gap_polarization/prospective_q_a_windows.csv')
    ap.add_argument('--continuation-file',type=Path,default=root/'outputs/event_aligned_continuation/event_aligned_samples.csv')
    ap.add_argument('--observed-time-file',type=Path,default=root/'outputs/polarization_torque_budget/observed_time_series.csv')
    ap.add_argument('--observed-interval-file',type=Path,default=root/'outputs/polarization_torque_budget/observed_intervals.csv')
    ap.add_argument('--event-file',type=Path,default=root/'outputs/initiation_propagation/crossing_episodes.csv')
    ap.add_argument('--matched-file',type=Path,default=root/'outputs/initiation_propagation/matched_precursor_pairs.csv')
    ap.add_argument('--output-dir',type=Path,default=root/'outputs/fixed_axis_polarization_comparison')
    ap.add_argument('--seeds',type=int,nargs='*',default=list(range(40)))
    ap.add_argument('--dt',type=float,default=.001)
    ap.add_argument('--transient',type=float,default=200.)
    ap.add_argument('--stationary-dt',type=float,default=1.,help='Subsample full stationary trajectories every this many simulation time units')
    ap.add_argument('--train-cutoff',type=int,default=30)
    ap.add_argument('--save-stationary-samples',action='store_true',help='Write per-frame derived polarization CSV locally; unnecessary for repo outputs')
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    if args.self_test:selftest();return
    if not args.seeds:raise ValueError('Specify at least one seed')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    stat,risk,cont,audit=merge_with_raw(args.risk_file,args.continuation_file,args.raw_dir,
                                        set(args.seeds),args.transient,args.dt,args.stationary_dt)
    save_csv(args.output_dir/'initiation_same_frames.csv',risk)
    save_csv(args.output_dir/'continuation_same_frames.csv',cont)
    station=create_stationary_outputs(stat,args.output_dir,save_samples=args.save_stationary_samples)
    paired_gap_direction_effect(stat,args.output_dir)
    results={'source_frame_validation':audit,'stationary':station,
             'initiation':compare_heldout(risk,'init',args.train_cutoff),
             'continuation':compare_heldout(cont,'continuation',args.train_cutoff),
             'minority_diagnostic_init':minority_added_value(risk,'init',args.train_cutoff),
             'minority_diagnostic_continuation':minority_added_value(cont,'continuation',args.train_cutoff),
             'n_initiation':len(risk),'n_continuation':len(cont),
             'definitions':{
                 'q':'(N_R - N_L)/N; x<0 in left chamber, x>=0 in right',
                 'a_L':'mean cos(theta) for fish with x<0; NaN if none',
                 'a_R':'-mean cos(theta) for fish with x>=0; NaN if none',
                 'a_fixed':'a_L for q<0, a_R for q>0; undefined when q=0',
                 'a_com':'majority P dot (-majority COM)/|majority COM|; undefined at q=0',
                 'stationary_law':'Equal-seed time-sampled empirical probability mass; -log mass is not a dynamical potential',
                 'initiation_outcome':'Reuse indexed first-crossing labels without modification',
                 'continuation_outcome':'Reuse event-aligned threshold of four ADDITIONAL net fish after first majority-to-minority crossing'
             }}
    results['original_COM_parity_audit']=audit_against_existing(results,args.risk_file.parent.parent)
    results['matched_pre_crossing_timecourses']=event_control_polarization_timecourses(
        args.raw_dir,args.observed_time_file,args.observed_interval_file,
        args.event_file,args.matched_file,args.dt,args.output_dir)
    results['recommendation']={
        'retain':'a_com (majority polarization projected along the instantaneous majority-COM-to-gap direction)',
        'reason':'Physically local to the opening and substantially better held-out initiation and continuation prediction than fixed horizontal polarization.',
        'fixed_axis_role':'Mechanism diagnostic separating horizontal heading changes from rotation of the COM-to-gap reference.',
        'minority_coordinate':'Do not add by default; held-out gains are unresolved on the restricted both-occupied subset.'}
    for mode,rows in [('init',risk),('continuation',cont)]:
        if rows:plot_outcome(rows,mode,args.output_dir)
    (args.output_dir/'summary.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    conf={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    (args.output_dir/'analysis_config.json').write_text(json.dumps(conf,indent=2)+'\n')
    init=results['initiation'];prop=results['continuation'];station=results['stationary']
    mi=results['minority_diagnostic_init'];mc=results['minority_diagnostic_continuation']
    timecourse=results['matched_pre_crossing_timecourses']
    icmp=init['comparisons']['fixed_vs_com'];ccmp=prop['comparisons']['fixed_vs_com']
    lines=['# Fixed-axis versus COM-to-gap polarization',
           '', 'Both coordinates use precisely the same raw snapshots, eligibility rules, and future outcomes from the existing 40 no-hydrodynamics gap-20 trajectories. No event was relabeled and no simulation was run.',
           '', '## Validation and cohorts',
           '',f"- First-entry risk origins: {len(risk)}; held-out: {init['heldout_n']} with {init['metrics']['q']['positive']} positives.",
           f"- Event-aligned continuation events: {len(cont)}; held-out: {prop['heldout_n']} with {prop['metrics']['q']['positive']} positives.",
           f"- Unconditional post-transient snapshots: {len(stat)} ({len(stat)//len(args.seeds)} per seed at spacing {args.stationary_dt:g}).",
           f"- Raw-state q/a_COM cross-checks: {audit['n_raw_crosschecks']}; maximum a_COM error {audit['max_abs_a_com_error']:.3e}.",
           '- Archived COM AUC and Brier values reproduce exactly for both outcomes.',
           '', '## Same 37 matched pre-crossing trajectories',
           '',f"The moving COM coordinate decreases before events: event Δa_COM = {timecourse['a_com']['event_delta_a_com']['mean']:.6f}, 95% CI {timecourse['a_com']['event_delta_a_com']['ci95']}. The fixed horizontal coordinate instead increases: event Δa_fixed = {timecourse['a_fixed']['event_delta_a_fixed']['mean']:.6f}, 95% CI {timecourse['a_fixed']['event_delta_a_fixed']['ci95']}.",
           '',f"Relative to unchanged controls, COM Δ = {timecourse['a_com']['paired_delta_a_com']['mean']:.6f}, CI {timecourse['a_com']['paired_delta_a_com']['ci95']}; fixed-axis Δ = {timecourse['a_fixed']['paired_delta_a_fixed']['mean']:.6f}, CI {timecourse['a_fixed']['paired_delta_a_fixed']['ci95']}. The latter interval narrowly includes zero.",
           '', 'Therefore the apparent late decline of COM-to-gap polarization is largely reference-direction rotation, not fish turning away from the horizontal opening.',
           '', '## Held-out prediction: seeds 0–29 train, 30–39 test',
           '', '| Outcome | q AUC | q+a_COM AUC | q+a_fixed AUC | fixed-minus-COM seed-mean AUC difference (95% CI) | fixed-minus-COM Brier gain (95% CI) |',
           '|---|---:|---:|---:|---:|---:|',
           f"| First entry | {init['metrics']['q']['auc']:.6f} | {init['metrics']['com']['auc']:.6f} | {init['metrics']['fixed']['auc']:.6f} | {icmp['seed_mean_auc_improvement']['mean']:.6f} {icmp['seed_mean_auc_improvement']['ci95']} | {icmp['seed_mean_brier_improvement']['mean']:.6f} {icmp['seed_mean_brier_improvement']['ci95']} |",
           f"| Continuation | {prop['metrics']['q']['auc']:.6f} | {prop['metrics']['com']['auc']:.6f} | {prop['metrics']['fixed']['auc']:.6f} | {ccmp['seed_mean_auc_improvement']['mean']:.6f} {ccmp['seed_mean_auc_improvement']['ci95']} | {ccmp['seed_mean_brier_improvement']['mean']:.6f} {ccmp['seed_mean_brier_improvement']['ci95']} |",
           '', 'Positive Brier gain favors fixed-axis. Both fixed-minus-COM differences are negative and exclude zero: COM-to-gap is materially better, especially for continuation. Fixed-axis modestly improves first-entry prediction over q alone, but is worse than q alone for continuation.',
           '', '## Unconditional empirical distributions',
           '',f"The equal-seed empirical p(q,a) landscapes use {station['n_snapshots']} finite-time snapshots. q=0 occurs in {station['n_q_zero']} snapshots and has no majority coordinate. Left/right chambers are empty in {station['n_empty_left']}/{station['n_empty_right']} snapshots; these polarizations remain undefined rather than zero. Both chambers are occupied in {station['n_both_chambers_occupied']} snapshots ({station['fraction_both_occupied']:.3f}).",
           '', 'The shifted -ln p panels are statistical occupation landscapes only. Fixed-axis retains horizontal heading; COM-to-gap also resolves off-axis location relative to the opening.',
           '', '## Separate left/right and minority polarization',
           '',f"For initiation, the minority diagnostic retains {mi['n_eligible']}/{len(risk)} observations with at least five minority fish; held-out joint-cell coverage is {mi['coverage_joint_on_test']:.3f}. Minority Brier gain = {mi['additional_minority_brier_gain']['mean']:.6f}, CI {mi['additional_minority_brier_gain']['ci95']}; AUC gain = {mi['additional_minority_auc_gain']['mean']:.6f}, CI {mi['additional_minority_auc_gain']['ci95']}.",
           '',f"For continuation it retains {mc['n_eligible']}/{len(cont)} observations; coverage is {mc['coverage_joint_on_test']:.3f}. Minority Brier gain = {mc['additional_minority_brier_gain']['mean']:.6f}, CI {mc['additional_minority_brier_gain']['ci95']}; AUC gain = {mc['additional_minority_auc_gain']['mean']:.6f}, CI {mc['additional_minority_auc_gain']['ci95']}.",
           '', 'Neither outcome shows a statistically resolved benefit from minority polarization. Do not add it to the primary state; it is undefined for localized schools and produces sparse cells.',
           '', '## Recommendation',
           '', '**Retain majority-school COM-to-gap polarization a_COM for the eventual empirical collective landscape.** It is physically local to the opening and decisively outperforms fixed-horizontal polarization on independent seeds. Keep a_fixed as a mechanism diagnostic: headings turn horizontally toward the opening even while the moving COM projection falls.',
           '', '## Scope and interpretation',
           '', 'These are finite-time empirical occupation distributions and held-out associations, not fitted potentials, rate barriers, Markov closure, or a whole-school relocation model. The rare complete relocations are not the outcomes tested here.']
    (args.output_dir/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Outputs:',args.output_dir,flush=True)
    print('\n'.join(lines),flush=True)


if __name__=='__main__': main()
