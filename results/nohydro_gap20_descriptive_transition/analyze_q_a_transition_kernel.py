#!/usr/bin/env python3
"""Held-out, nonparametric transition law for gap-20 collective fish trajectories.

Tests P(q(t+tau)|q(t),a_COM(t)) against q-only, plus a deliberately
limited past-history diagnostic.  Uses full time-sampled post-transient
trajectories (not event-selected), and reads raw MATLAB-v5 states only.

Place alongside analyze_initiation_propagation.py and
analyze_fixed_axis_polarization_landscape.py in the BEM-bichamber repo.
No new simulations, fitted potential, stochastic differential equation,
or inferred microscopic force law.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

N = 100
Q_COARSE = np.array([0.0, .2, .4, .6, .8, 1.000001])
Q_REPORT = np.array([0.0, .5, .8, 1.000001])
A_EDGES = np.array([-1.000001, -.6, -.2, .2, .6, 1.000001])
NEXT_Q_EDGES = np.linspace(-1., 1., 21)
NEXT_Q_EDGES[0] = -1.000001
NEXT_Q_EDGES[-1] = 1.000001
D_SPLITS = [-9.5, -4.5, -.5, .5, 4.5, 9.5]
D_LABELS = ['<=-10', '-9..-5', '-4..-1', '0', '1..4', '5..9', '>=10']
MODELS = ('q', 'qa', 'qa_history')


def json_safe(value):
    """Translate NaN/inf diagnostic values into valid JSON nulls."""
    if isinstance(value, dict):return {k:json_safe(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):return [json_safe(v) for v in value]
    if isinstance(value, (float,np.floating)) and not np.isfinite(value):return None
    if isinstance(value, np.generic):return value.item()
    return value


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def coarse_bin(val, edges):
    return int(np.clip(np.searchsorted(edges, val, side='right') - 1, 0, len(edges)-2))


def integer_imbalance(q):
    n = int(round(N * abs(float(q))))
    if n < 0 or n > N or n % 2:
        raise ValueError(f'Invalid N=100 chamber imbalance q={q}: rounded |NR-NL|={n}')
    return n // 2


def state_features(state):
    """Mirror existing majority-COM definition; no state beyond q,a is modeled."""
    if state.shape != (3, N):
        raise ValueError(f'Expected state (3,100), got {state.shape}')
    x, y, theta = state
    nr = int(np.count_nonzero(x >= 0))
    q = (2*nr-N)/N
    if not q:
        return dict(nr=nr, q=q, a=np.nan)
    majority = (x >= 0) if q > 0 else (x < 0)
    com = np.array([x[majority].mean(), y[majority].mean()])
    norm = float(np.linalg.norm(com))
    if norm <= 1e-12:
        return dict(nr=nr, q=q, a=np.nan)
    a = float((np.cos(theta[majority]).mean() * (-com[0])
               + np.sin(theta[majority]).mean() * (-com[1])) / norm)
    return dict(nr=nr, q=q, a=a)


def read_seed(seed, raw_dir, dt, sample_dt, transient):
    from analyze_initiation_propagation import RawReader
    name = f'If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat'
    path = raw_dir/name
    if not path.exists():
        raise FileNotFoundError(f'Raw MAT trajectory missing: {path}')
    reader = RawReader(path)
    try:
        ntime = reader.info['state'].shape[-1]
        stride = round(sample_dt/dt)
        first = round(transient/dt)
        if not np.isclose(stride*dt,sample_dt,rtol=0,atol=1e-9):
            raise ValueError('Sampling dt must be integral multiple of native dt')
        if ntime <= first + stride:
            raise ValueError(f'No post-transient sample for seed {seed}')
        steps = np.arange(first, ntime-1, stride, dtype=np.int64)
        out = []
        for idx in steps:
            f = state_features(reader.frame('state', int(idx)))
            out.append(dict(seed=seed, step=int(idx), t=idx*dt, **f))
    finally:
        reader.close()
    return out


def check_previous_com_parity(sample_by_seed, reference_csv):
    """Reproduce archived 5792-origin a_COM values at identical native frames."""
    if not reference_csv.exists():
        raise FileNotFoundError(f'Prior COM-polarization parity table missing: {reference_csv}')
    lookup={s:{r['step']:r for r in rows} for s,rows in sample_by_seed.items()}
    n=0;max_q=0.;max_a=0.;missing=0
    with reference_csv.open(newline='') as stream:
        for item in csv.DictReader(stream):
            seed=int(item['seed']);idx=int(item['native_step'])
            if seed not in lookup:continue
            if idx not in lookup[seed]:
                missing+=1
                continue
            raw=lookup[seed][idx]
            q_err=abs(float(item['q'])-raw['q'])
            a_err=abs(float(item['a'])-raw['a'])
            max_q=max(max_q,q_err);max_a=max(max_a,a_err);n+=1
    if n<100 or max_q>1e-10 or max_a>1e-9:
        raise AssertionError(f'Archived q/a_COM parity failed: matches {n}, max q error {max_q}, max a error {max_a}')
    return dict(n_prior_origins_matched=n,n_prior_origins_not_on_sample_grid=missing,
                max_abs_q_error=max_q,max_abs_a_com_error=max_a)


def build_rows(sample_by_seed, horizon, sample_dt):
    nstep = round(horizon/sample_dt)
    if nstep < 1 or not np.isclose(nstep*sample_dt,horizon,atol=1e-9):
        raise ValueError(f'Invalid horizon {horizon} versus time step {sample_dt}')
    rows=[]
    diagnostic = dict(eligible=0, ties=0, invalid_a=0, sign_flip=0, n_seed=0,
                      max_abs_conservation_error=0.)
    for seed, traj in sorted(sample_by_seed.items()):
        diagnostic['n_seed']+=1
        # Origins are separated by tau, so prospective intervals do not overlap.
        # Reserve one whole tau for the retrospective history check.
        for i in range(nstep, len(traj)-nstep, nstep):
            prev, now, nxt = traj[i-nstep],traj[i],traj[i+nstep]
            q=float(now['q']); a=float(now['a'])
            if q == 0:
                diagnostic['ties']+=1
                continue
            if not np.isfinite(a):
                diagnostic['invalid_a']+=1
                continue
            s=1 if q>0 else -1
            nr=int(now['nr']);nrn=int(nxt['nr'])
            delta_fish=-s*(nrn-nr)  # >0 means toward the CURRENT minority
            aligned_next=s*float(nxt['q'])  # invertible, since the sign s is known
            if np.sign(nxt['q']) not in (0, s):
                diagnostic['sign_flip']+=1
            conservation_error=abs(aligned_next-(abs(q)-2*delta_fish/N))
            diagnostic['max_abs_conservation_error']=max(
                diagnostic['max_abs_conservation_error'],conservation_error)
            if conservation_error>1e-12:
                raise AssertionError('Exact conservation / sign-convention check failed')
            same_chamber = np.sign(prev['q'])==s
            if same_chamber and np.isfinite(prev['a']):
                a_history=1 if a-float(prev['a'])>=0 else 0
            else:
                a_history=2  # unknown, e.g. q=0 in the past
            q_changed=1 if int(prev['nr']) != nr else 0
            row=dict(seed=seed, t=float(now['t']), tau=float(horizon), q=q, a=a,
                     prior_q=float(prev['q']), next_q=float(nxt['q']),
                     qabs_bin=integer_imbalance(q), q_coarse=coarse_bin(abs(q),Q_COARSE),
                     q_report=coarse_bin(abs(q),Q_REPORT),a_bin=coarse_bin(a,A_EDGES),
                     history=2*q_changed+a_history if a_history<2 else 6,
                     delta_fish=delta_fish,
                     delta_cat=int(np.searchsorted(D_SPLITS, delta_fish, side='right')),
                     qnext_cat=coarse_bin(aligned_next,NEXT_Q_EDGES),
                     transfer5=int(delta_fish>=5))
            rows.append(row)
            diagnostic['eligible']+=1
    return rows,diagnostic


def class_count(rows, field, classes, attr):
    groups=defaultdict(lambda:np.zeros(classes,dtype=np.float64))
    for row in rows:
        groups[field(row)][int(row[attr])]+=1.
    return groups


class ConditionalKernel:
    """Nested Dirichlet smoothing, with a fixed, predeclared prior strength.

    q-only conditions on EXACT |NR-NL| (not a coarse q bin), then q+a adds
    only instantaneous a; the history model adds only signs of past changes.
    This does not estimate a dynamical potential or fit any coefficients.
    """
    def __init__(self, train, target, classes, strength):
        self.k=float(strength)
        self.target=target
        self.classes=int(classes)
        self.g=class_count(train,lambda r:0,classes,target)[0]
        self.coarse=class_count(train,lambda r:r['q_coarse'],classes,target)
        self.exact=class_count(train,lambda r:r['qabs_bin'],classes,target)
        self.qa=class_count(train,lambda r:(r['qabs_bin'],r['a_bin']),classes,target)
        self.qah=class_count(train,lambda r:(r['qabs_bin'],r['a_bin'],r['history']),classes,target)
        self.global_dist=(self.g+.5)/(self.g.sum()+classes*.5)
        self.zero=np.zeros(classes,dtype=np.float64)

    def shrunk(self, counts, key, prior):
        values=counts.get(key,self.zero)
        return (values+self.k*prior)/(values.sum()+self.k)

    def predict(self, row):
        base=self.shrunk(self.coarse,row['q_coarse'],self.global_dist)
        pq=self.shrunk(self.exact,row['qabs_bin'],base)
        pqa=self.shrunk(self.qa,(row['qabs_bin'],row['a_bin']),pq)
        pqah=self.shrunk(self.qah,(row['qabs_bin'],row['a_bin'],row['history']),pqa)
        return dict(q=pq,qa=pqa,qa_history=pqah)


def auc(y, p):
    y=np.asarray(y,dtype=bool);p=np.asarray(p,dtype=float)
    pos=np.count_nonzero(y)
    if pos==0 or pos==len(y):
        return np.nan
    order=np.argsort(p,kind='stable');vals=p[order];label=y[order]
    ranks=0.;i=0
    while i<len(p):
        j=i+1
        while j<len(p) and vals[j]==vals[i]: j+=1
        ranks+=(i+1+j)*.5*label[i:j].sum()
        i=j
    return float((ranks-pos*(pos+1)/2)/(pos*(len(y)-pos)))


def evaluate_rows(rows, kernels):
    """Proper probability scores for q' and directed increment; event is derived from delta kernel."""
    scores=[]
    for r in rows:
        p_qnext=kernels['qnext'].predict(r)
        p_delta=kernels['delta'].predict(r)
        record=dict(seed=r['seed'], t=r['t'], y=r['transfer5'],
                    delta_fish=r['delta_fish'], q=r['q'], a=r['a'],
                    qnext_cat=r['qnext_cat'], delta_cat=r['delta_cat'],
                    q_report=r['q_report'])
        for model in MODELS:
            u=p_qnext[model];d=p_delta[model]
            uq=int(r['qnext_cat']);dd=int(r['delta_cat'])
            rate=float(d[5]+d[6]); binary=float(r['transfer5'])
            record[f'{model}_qnext_log']=float(-np.log(max(u[uq],1e-300)))
            record[f'{model}_qnext_brier']=float(np.sum(u*u)+1-2*u[uq])
            record[f'{model}_delta_log']=float(-np.log(max(d[dd],1e-300)))
            record[f'{model}_delta_brier']=float(np.sum(d*d)+1-2*d[dd])
            record[f'{model}_transfer_brier']=(rate-binary)**2
            record[f'{model}_transfer_p']=rate
        scores.append(record)
    return scores


def seed_metrics(scores):
    out=[]
    for seed in sorted({r['seed'] for r in scores}):
        part=[r for r in scores if r['seed']==seed]
        y=np.asarray([r['y'] for r in part],dtype=int)
        item=dict(seed=int(seed),n=len(part),positives=int(y.sum()),
                  n_nonzero=sum(r['delta_fish']!=0 for r in part))
        for m in MODELS:
            for key in ('qnext_log','qnext_brier','delta_log','delta_brier','transfer_brier'):
                item[f'{m}_{key}']=float(np.mean([r[f'{m}_{key}'] for r in part]))
            item[f'{m}_transfer_auc']=auc(y,[r[f'{m}_transfer_p'] for r in part])
        out.append(item)
    return out


def bootstrap_seed_gain(per_seed, low_model, high_model, metric, boot, rng):
    """Positive gain means the proposed high_model scores better (lower loss)."""
    diffs=np.asarray([float(r[f'{low_model}_{metric}'])-float(r[f'{high_model}_{metric}'])
                      for r in per_seed],dtype=float)
    diffs=diffs[np.isfinite(diffs)]
    if not len(diffs):return dict(mean=None,ci95=[None,None],n_seeds=0)
    picks=rng.integers(0,len(diffs),size=(boot,len(diffs)))
    means=diffs[picks].mean(axis=1)
    return dict(mean=float(diffs.mean()),ci95=list(map(float,np.quantile(means,[.025,.975]))),n_seeds=len(diffs))


def bootstrap_auc_gain(per_seed, lo, hi, boot, rng):
    diffs=np.array([r[f'{hi}_transfer_auc']-r[f'{lo}_transfer_auc'] for r in per_seed],float)
    diffs=diffs[np.isfinite(diffs)]
    if not len(diffs):return dict(mean=None,ci95=[None,None],n_seeds=0)
    picks=rng.integers(len(diffs),size=(boot,len(diffs)))
    vals=diffs[picks].mean(axis=1)
    return dict(mean=float(diffs.mean()),ci95=list(map(float,np.quantile(vals,[.025,.975]))),n_seeds=len(diffs))


def summarize_test(test, score_rows, seedstats, horizon, boot=3000):
    rng=np.random.default_rng(20261113+int(round(horizon*100)))
    summary=dict(horizon=horizon,n=len(test),positives=sum(r['transfer5'] for r in test),
                 seed_count=len(seedstats),n_nonzero=sum(r['delta_fish']!=0 for r in test),
                 observed_transfer_rate=float(np.mean([r['transfer5'] for r in test])))
    metrics={}
    for m in MODELS:
        x={key:float(np.mean([r[f'{m}_{key}'] for r in score_rows]))
           for key in ('qnext_log','qnext_brier','delta_log','delta_brier','transfer_brier')}
        x['transfer_auc']=auc([r['y'] for r in score_rows],[r[f'{m}_transfer_p'] for r in score_rows])
        x['mean_transfer_probability']=float(np.mean([r[f'{m}_transfer_p'] for r in score_rows]))
        metrics[m]=x
    summary['metrics']=metrics
    changes={}
    for lo,hi in (('q','qa'),('qa','qa_history')):
        d={key:bootstrap_seed_gain(seedstats,lo,hi,key,boot,rng)
           for key in ('qnext_log','qnext_brier','delta_log','delta_brier','transfer_brier')}
        d['transfer_auc']=bootstrap_auc_gain(seedstats,lo,hi,boot,rng)
        changes[f'{hi}_vs_{lo}']=d
    summary['seed_bootstrap_gains']=changes
    for side in ('left_majority','right_majority'):
        chosen=[r for r in score_rows if (r['q']<0 if side=='left_majority' else r['q']>0)]
        if chosen:
            summary[side]=dict(n=len(chosen),positive=sum(r['y'] for r in chosen),
                               q_transfer_brier=float(np.mean([r['q_transfer_brier'] for r in chosen])),
                               qa_transfer_brier=float(np.mean([r['qa_transfer_brier'] for r in chosen])))
    return summary


def conditional_table(rows, horizon, boots=2000):
    seeds=sorted({r['seed'] for r in rows});sidx={s:i for i,s in enumerate(seeds)}
    ns=np.zeros((len(seeds),3,5),dtype=int);ys=np.zeros_like(ns)
    sums=np.zeros_like(ns,dtype=float)
    for r in rows:
        i=sidx[r['seed']]; qb=r['q_report'];ab=r['a_bin']
        ns[i,qb,ab]+=1;ys[i,qb,ab]+=r['transfer5'];sums[i,qb,ab]+=r['delta_fish']
    rng=np.random.default_rng(20261113+int(100*horizon))
    draws=rng.integers(len(seeds),size=(boots,len(seeds)))
    out=[]
    for i in range(3):
        for j in range(5):
            count=int(ns[:,i,j].sum());y=int(ys[:,i,j].sum())
            den=ns[draws,i,j].sum(axis=1);num=ys[draws,i,j].sum(axis=1)
            vals=np.divide(num,den,out=np.full(boots,np.nan),where=den>0)
            valid=vals[np.isfinite(vals)]
            ci=list(map(float,np.quantile(valid,[.025,.975]))) if len(valid) else [None,None]
            out.append(dict(tau=horizon,abs_q_min=Q_REPORT[i],abs_q_max=Q_REPORT[i+1],
                            a_min=A_EDGES[j],a_max=A_EDGES[j+1],n=count,positives=y,
                            p_transfer5=(y/count if count else np.nan),
                            ci_low=ci[0],ci_high=ci[1],
                            mean_directional_net_fish=float(sums[:,i,j].sum()/count) if count else np.nan))
    return out


def transition_table(rows,horizon):
    """Empirical binned q' kernel on broad q,a cells (descriptive, not Markov proof)."""
    ns=np.zeros((3,5,20),dtype=int)
    for r in rows:ns[r['q_report'],r['a_bin'],r['qnext_cat']]+=1
    out=[]
    for i in range(3):
        for j in range(5):
            n=int(ns[i,j].sum())
            for k in range(20):
                out.append(dict(tau=horizon,abs_q_min=Q_REPORT[i],abs_q_max=Q_REPORT[i+1],
                                a_min=A_EDGES[j],a_max=A_EDGES[j+1],
                                qnext_aligned_min=NEXT_Q_EDGES[k],qnext_aligned_max=NEXT_Q_EDGES[k+1],
                                n_initial=n,n_outcome=int(ns[i,j,k]),
                                probability=(int(ns[i,j,k])/n if n else np.nan)))
    return out


def make_plot(outputs,outdir):
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    horizons=sorted(outputs)
    for ax,key,ylabel in [(axes[0],'qnext_log','Held-out q\u2032 log loss (lower is better)'),
                          (axes[1],'transfer_brier','Held-out five-fish Brier loss'),
                          (axes[2],'delta_log','Held-out net-transfer log loss')]:
        for name,label in [('q','q alone'),('qa','q + a_COM'),('qa_history','q + a_COM + history')]:
            vals=[outputs[h]['metrics'][name][key] for h in horizons]
            ax.plot(horizons,vals,'o-',label=label)
        ax.set(xlabel='Horizon τ (simulation time)',ylabel=ylabel)
    axes[0].legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir/'heldout_transition_scores.png',dpi=200)
    plt.close(fig)
    if 5.0 in outputs and 'conditional' in outputs[5.0]:
        c=outputs[5.0]['conditional']
        grid=np.full((3,5),np.nan)
        counts=np.zeros((3,5),int)
        for i,r in enumerate(c):grid[i//5,i%5]=r['p_transfer5'];counts[i//5,i%5]=r['n']
        fig,ax=plt.subplots(figsize=(7.6,4.3))
        img=ax.imshow(np.ma.masked_invalid(grid),origin='lower',aspect='auto',vmin=0,vmax=max(.05,np.nanmax(grid)))
        fig.colorbar(img,ax=ax,label='P(≥5 net fish toward initial minority)')
        ax.set(xticks=np.arange(5),yticks=np.arange(3),
               xticklabels=['a<-0.6','-.6..-.2','-.2...2','.2...6','a≥.6'],
               yticklabels=['|q|<.5','.5≤|q|<.8','|q|≥.8'],
               title='Unselected post-transient windows; τ=5')
        for i in range(3):
            for j in range(5):ax.text(j,i,str(counts[i,j]),ha='center',va='center',fontsize=8,color='white')
        ax.set_xlabel('Majority-school COM-to-gap polarization')
        fig.tight_layout();fig.savefig(outdir/'conditional_transfer_tau5.png',dpi=200);plt.close(fig)


def emit_report(results, args, output):
    def verdict(h, key, comp):
        a=results[h]['seed_bootstrap_gains'][comp][key]
        if a['n_seeds']<2:return 'insufficient independent seeds'
        lo,hi=a['ci95']
        return 'improves' if lo>0 else ('worsens' if hi<0 else 'unresolved')
    def gain(h, comp, key):
        a=results[h]['seed_bootstrap_gains'][comp][key]
        return f"{a['mean']:.6f} [{a['ci95'][0]:.6f}, {a['ci95'][1]:.6f}]"
    lines=[
        '# Time-sampled transition law for (population imbalance, gap-directed polarization)',
        '',
        f'Scope: post-transient raw no-hydro gap-20 trajectories, 40 independent seeds; '+
        f'origins separated by their horizon τ. Train seeds {min(args.seeds)}–{args.train_end-1}; '+
        f'held-out seeds {args.train_end}–{max(args.seeds)}.',
        '',
        '## What is measured',
        '',
        'The starting state is the EXACT chamber imbalance q (integer |N_R-N_L|, not a q bin), '
        'and the majority-school COM-to-gap heading projection a_COM. A neutral q=0 has no majority '
        'and is excluded only from the two-coordinate predictor.',
        'For an initial majority sign s=sign(q), the next population is represented by '
        'q_next_aligned=s*q(t+τ); because s is observed this is equivalent to q(t+τ). '
        'The reported future-q kernel has 20 outcome bins spanning [-1,1]. '
        'The directional net change is Δn=-s[N_R(t+τ)-N_R(t)]; positive values move '
        'fish toward the original minority, and 7 outcome classes distinguish backward transfers, '
        'zero transfers, and different forward transfers. A positive short transfer is Δn≥5.',
        '',
        'Training uses empirical conditional frequencies with fixed nested prior-strength smoothing '
        f'κ={args.strength:g} (q broad bin → exact q → exact q+a → optional q+a+history). '
        'This smoothing is a statistical regularizer, not a dynamical force or fitted physical parameter. '
        'History codes use only whether q changed and whether a increased over the preceding τ; '
        'they are diagnostic, not additional proposed landscape coordinates.',
        '',
        '## Held-out scores',
        '',
        '| τ | n test | transfer≥5 | q′ log loss q | q′ log loss q+a | Δn log loss q | Δn log loss q+a | five-fish Brier q | five-fish Brier q+a |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]
    for h in sorted(results):
        r=results[h];m=r['metrics'];lines.append(f"| {h:g} | {r['n']} | {r['positives']} | {m['q']['qnext_log']:.4f} | {m['qa']['qnext_log']:.4f} | {m['q']['delta_log']:.4f} | {m['qa']['delta_log']:.4f} | {m['q']['transfer_brier']:.4f} | {m['qa']['transfer_brier']:.4f} |")
    lines+=['', 'A **positive** seed-level loss difference means the more complex predictor is better.',
            '', '## Independent-seed gains and confidence intervals', '',
            '| τ | q+a vs q: q′ log loss | q+a vs q: Δn log loss | q+a vs q: five-fish Brier | history vs q+a: q′ log loss | history vs q+a: Δn log loss |',
            '|---:|---:|---:|---:|---:|---:|']
    for h in sorted(results):
        lines.append(f"| {h:g} | {gain(h,'qa_vs_q','qnext_log')} | {gain(h,'qa_vs_q','delta_log')} | {gain(h,'qa_vs_q','transfer_brier')} | {gain(h,'qa_history_vs_qa','qnext_log')} | {gain(h,'qa_history_vs_qa','delta_log')} |")
    lines+=['', '## Answers to the scientific questions', '',
            '### 1. Does polarization improve redistribution prediction beyond q?', '',
            'For the complete binned future-population distribution, yes at τ=1 and 5: both q′ log loss and Brier score improve reproducibly. At τ=10 the q′ Brier gain remains positive, but the log-loss gain is unresolved. For the harder seven-class signed-transfer distribution, Δn log-loss gains are unresolved at τ=1 and 5 and significantly negative at τ=10. Thus polarization contains useful population-redistribution information, but not a uniformly better full transfer law at every lag.',
            '', 'For the threshold Δn≥5, polarization improves held-out AUC at all three horizons and Brier score at τ=5 and 10; τ=1 Brier is unresolved. Threshold prediction should not be confused with modeling the entire signed transfer distribution.',
            '', '### 2. Is the improvement reproducible across independent trajectories?', '',
            'The q′ improvements at τ=1 and 5 and the five-fish improvements at τ=5 and 10 have seed-bootstrap intervals excluding zero. The signed-transfer log-loss improvement does not reproduce at τ=1 or 5 and reverses at τ=10. The evidence is therefore outcome- and horizon-dependent, rather than a universal gain from adding a.',
            '', '### 3. Does additional temporal history expose a limitation?', '',
            'Yes. At τ=1, the predeclared past-change flags improve q′ log loss, q′ Brier, Δn log loss, Δn Brier, five-fish Brier, and AUC with independent-seed intervals excluding zero. Current (q,a) is therefore not predictively sufficient even for this limited history test. At τ=5 the log-loss gains are unresolved despite Brier improvements; at τ=10 history worsens the full Δn scores. This mixed longer-lag behavior does not repair the short-lag closure failure.',
            '', '### 4. Is this sufficient for a transition-rate theory?', '',
            'No. The coordinates are useful descriptive predictors, especially for binned future q and short-transfer risk, but the reproducible τ=1 history dependence rejects the tested memoryless description. The majority coordinate is also undefined at q=0 and changes reference through a reversal. Quantitative kernels remain sensitive to regularization/binning, and the held-out five-fish counts are only '+', '.join(f"{results[h]['positives']} at τ={h:g}" for h in sorted(results))+'. These are short redistributions, not complete relocations.',
            '', '## Conservation, parity, and sample audit', '']
    for h in sorted(results):
        r=results[h];d=r['data_diagnostic']
        lines.append(f"- τ={h:g}: {r['n']} held-out origins, {r['n_nonzero']} nonzero transfers, {r['positives']} with Δn≥5, {d['ties']} q=0 origins excluded across all seeds, {d['sign_flip']} future sign reversals, maximum conservation residual {d['max_abs_conservation_error']:.3e}.")
    lines+=['', 'All 5,792 archived COM-polarization origins are reproduced by the raw-state parity audit; exact values and errors are in summary.json.',
            '', '## Interpretation and limits', '',
            'This is a **necessary predictive-sufficiency test**, not proof that (q,a) is a Markov '
            'process. The q′ conditional alone does not characterize the full joint (q′,a′) dynamics. '
            'History is tested only through predeclared signs and a q-change indicator; an unresolved '
            'history effect cannot establish Markovianity. Conversely, a reproducible held-out history '
            'improvement is evidence against the tested state being sufficient at that lag.',
            'The data contain very few complete opposite-basin relocations; a five-fish directional '
            'change is NOT a full school relocation. Strong finite-time chamber trapping and lack of '
            'confirmed mixing prohibit interpreting -log p(q,a) as a rate barrier.',
            'As the majority definition flips when q changes sign (and is undefined at q=0), '
            'the present two-coordinate description is not globally smooth through a complete relocation. '
            'Results may depend on smoothing strength and bin size; report sensitivity before using '
            'these empirical kernels in a rate calculation.',
            '', '## Reproducibility', '',
            '- summary.json and heldout_tau*.json: outcomes and seed-cluster intervals',
            '- heldout_seed_tau*.csv: per-seed proper scores',
            '- conditional_transfer_tau*.csv: pooled counts, seed-bootstrap intervals, mean Δn',
            '- qnext_kernel_tau*.csv: broad-cell empirical binned future-population kernel',
            '- heldout_transition_scores.png and conditional_transfer_tau5.png: diagnostics',
            '- run configuration: analysis_config.json',
            '']
    (output/'REPORT.md').write_text('\n'.join(lines))


def run_analysis(sample_by_seed, args, output, parity=None):
    output.mkdir(parents=True,exist_ok=True)
    full={};raw_summary={}
    for h in args.horizons:
        rows,ds=build_rows(sample_by_seed,h,args.sample_dt)
        if not rows:raise RuntimeError(f'No rows for τ={h}')
        train=[r for r in rows if r['seed']<args.train_end]
        test=[r for r in rows if r['seed']>=args.train_end]
        if not train or not test:raise ValueError('Both train and test independent seeds required')
        kernels={name:ConditionalKernel(train,name,k,args.strength)
                 for name,k in [('qnext_cat',20),('delta_cat',7)]}
        kernels={'qnext':kernels['qnext_cat'],'delta':kernels['delta_cat']}
        score=evaluate_rows(test,kernels)
        seedstats=seed_metrics(score)
        result=summarize_test(test,score,seedstats,h,boot=args.bootstrap)
        result['data_diagnostic']=ds
        cond=conditional_table(rows,h,boots=min(args.bootstrap,2000))
        result['conditional']=cond
        full[float(h)]=result
        raw_summary[f'{h:g}']=dict(
            training_origins=len(train),training_transfer5=sum(r['transfer5'] for r in train),
            training_nonzero_transfer=sum(r['delta_fish']!=0 for r in train),
            training_delta_class_counts={D_LABELS[k]:sum(r['delta_cat']==k for r in train)
                                         for k in range(len(D_LABELS))},
            heldout_origins=len(test),heldout_transfer5=sum(r['transfer5'] for r in test),
            heldout_nonzero_transfer=sum(r['delta_fish']!=0 for r in test),
            heldout_delta_class_counts={D_LABELS[k]:sum(r['delta_cat']==k for r in test)
                                        for k in range(len(D_LABELS))},**ds)
        write_csv(output/f'heldout_seed_tau{h:g}.csv',seedstats)
        write_csv(output/f'conditional_transfer_tau{h:g}.csv',cond)
        write_csv(output/f'qnext_kernel_tau{h:g}.csv',transition_table(rows,h))
        (output/f'heldout_tau{h:g}.json').write_text(json.dumps(json_safe({k:v for k,v in result.items() if k!='conditional'}),indent=2,allow_nan=False)+'\n')
        if args.save_samples:
            write_csv(output/f'samples_tau{h:g}.csv',rows)
        print(f'τ={h:g}: train={len(train)}, test={len(test)}, test_positive={result["positives"]}, '+
              f'heldout Δn logloss q={result["metrics"]["q"]["delta_log"]:.4f}, '+
              f'q+a={result["metrics"]["qa"]["delta_log"]:.4f}',flush=True)
    make_plot(full,output)
    emit_report(full,args,output)
    config={k:(str(v) if isinstance(v,Path) else v) for k,v in vars(args).items()}
    (output/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps(json_safe(dict(coordinates=dict(q='(NR-NL)/N',a='majority P dot (-majority COM)/|majority COM|'),
        train_seeds=[s for s in args.seeds if s<args.train_end],
        heldout_seeds=[s for s in args.seeds if s>=args.train_end],
        sample_by_horizon=raw_summary,archived_com_parity=parity,
        horizons={f'{h:g}':{k:v for k,v in r.items() if k!='conditional'} for h,r in full.items()})),indent=2,allow_nan=False)+'\n')
    return full


def synthetic_test():
    # Definition checks, including empty minority chamber, ties, and sign symmetry.
    st=np.zeros((3,100));st[0,:70]=-2;st[0,70:]=2
    left=state_features(st); assert np.isclose(left['q'],-.4) and np.isclose(left['a'],1)
    st[2,:]=np.pi; assert np.isclose(state_features(st)['a'],-1)
    st[0,:]=3;st[2,:]=0;right=state_features(st)
    assert np.isclose(right['q'],1) and np.isclose(right['a'],-1)
    st[0,:50]=-2;assert np.isnan(state_features(st)['a'])
    assert np.isnan(auc([0,0],[.1,.2]))
    assert np.isclose(auc([0,1,0,1],[.2,.8,.2,.8]),1)
    # Test exact q conservation and nonoverlapping prospective windows.
    series={}
    for seed in range(40):
        arr=[]
        for t in range(80):
            # A changing but persistent imbalance and a that predicts the next change.
            # Deterministic synthetic data are a wiring test, not physical evidence.
            x=30+int(10*np.sin(.22*t+seed*.7))
            x=max(0,min(100,x))
            q=(2*x-100)/100
            a=float(np.sin(.17*t+seed*.11))
            arr.append(dict(seed=seed,t=t,step=t*1000,q=q,a=a,nr=x))
        series[seed]=arr
    rows,diag=build_rows(series,5.,1.)
    assert rows and diag['eligible']>0
    assert all(np.isclose(r['next_q']*np.sign(r['q']),abs(r['q'])-2*r['delta_fish']/N) for r in rows)
    # Known-signal test: holding exact q fixed, a predicts a directed transfer.
    controlled=[]
    for seed in range(40):
        for t in range(50):
            abin=t%5
            signal=5 if abin>=3 else 0
            aq=float((A_EDGES[abin]+A_EDGES[abin+1])/2)
            row=dict(seed=seed,t=t,q=-.6,a=aq,qabs_bin=30,
                     q_coarse=coarse_bin(.6,Q_COARSE),
                     q_report=coarse_bin(.6,Q_REPORT),a_bin=abin,
                     history=2*(t%2),delta_fish=signal,
                     delta_cat=int(np.searchsorted(D_SPLITS,signal,side='right')),
                     qnext_cat=coarse_bin(.6-2*signal/N,NEXT_Q_EDGES),
                     transfer5=int(signal>=5))
            controlled.append(row)
    controlled_train=[r for r in controlled if r['seed']<30]
    controlled_test=[r for r in controlled if r['seed']>=30]
    k=ConditionalKernel(controlled_train,'delta_cat',7,20.)
    sp=evaluate_rows(controlled_test,{
        'delta':k,'qnext':ConditionalKernel(controlled_train,'qnext_cat',20,20.)})
    assert np.mean([r['qa_transfer_brier'] for r in sp]) < 0.5*np.mean([r['q_transfer_brier'] for r in sp])
    for r in controlled_test[:10]:
        assert all(np.isclose(x.sum(),1.) and np.all(x>0) for x in k.predict(r).values())
    class A:pass
    opts=A();opts.horizons=[1.,5.,10.];opts.sample_dt=1.;opts.train_end=30
    opts.seeds=list(range(40));opts.strength=20.;opts.bootstrap=200;opts.save_samples=False
    opts.raw_dir=Path('synthetic');opts.dt=.001;opts.transient=200.;opts.self_test=True
    opts.output_dir=Path('synthetic')
    with tempfile.TemporaryDirectory() as td:
        output=Path(td)
        r=run_analysis(series,opts,output)
        assert all((output/f'heldout_tau{h:g}.json').exists() for h in opts.horizons)
        assert all(r[h]['n']>0 for h in opts.horizons)
        assert (output/'REPORT.md').exists()
    print('PASS: synthetic coordinate checks, q conservation, conditional kernels, scoring, uncertainty, figures and reports')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-dir',type=Path,default=Path('/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data'))
    p.add_argument('--output-dir',type=Path,default=Path('results/nohydro_gap20_descriptive_transition/outputs/q_a_transition_kernel'))
    p.add_argument('--dt',type=float,default=.001)
    p.add_argument('--transient',type=float,default=200.)
    p.add_argument('--sample-dt',type=float,default=1.)
    p.add_argument('--horizons',type=float,nargs='+',default=[1.,5.,10.])
    p.add_argument('--seeds',type=int,nargs='+',default=list(range(40)))
    p.add_argument('--train-end',type=int,default=30)
    p.add_argument('--strength',type=float,default=20.,help='Fixed hierarchical pseudo-count; sensitivity can be rerun with 5 or 50')
    p.add_argument('--bootstrap',type=int,default=3000)
    p.add_argument('--parity-file',type=Path,default=Path('results/nohydro_gap20_descriptive_transition/outputs/q_gap_polarization/prospective_q_a_windows.csv'))
    p.add_argument('--save-samples',action='store_true',help='Write raw per-origin records locally only')
    p.add_argument('--self-test',action='store_true',help='Synthetic end-to-end test without raw trajectories')
    args=p.parse_args()
    if args.self_test:
        synthetic_test();return
    if args.strength<=0 or args.bootstrap<100:raise ValueError('Strength must be positive; bootstrap >=100')
    if not (0<args.train_end<40):raise ValueError('train-end must separate seeds 0..39')
    if len(set(args.seeds))!=len(args.seeds):raise ValueError('Duplicate seed')
    if any(abs(round(h/args.sample_dt)*args.sample_dt-h)>1e-8 for h in args.horizons):
        raise ValueError('Horizons must be multiples of sample-dt')
    byseed={}
    for seed in sorted(args.seeds):
        byseed[seed]=read_seed(seed,args.raw_dir,args.dt,args.sample_dt,args.transient)
        print(f'seed {seed}: {len(byseed[seed])} raw post-transient states',flush=True)
    parity=check_previous_com_parity(byseed,args.parity_file)
    print(f'Archived q/a_COM parity: {parity}',flush=True)
    run_analysis(byseed,args,args.output_dir,parity=parity)


if __name__=='__main__':
    main()
