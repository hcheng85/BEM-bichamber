#!/usr/bin/env python3
"""Source-grounded collective (q,a,b) pathway and diagnostic of stochastic closure.

Uses only archived no-hydrodynamics gap-20 raw MAT trajectories and previously
indexed crossing episodes. No simulator runs, effective potential, or fitted
phenomenological drift. Output is descriptive and an independently held-out
conditional-probability *diagnostic*, not a learned physical theory.

Place alongside analyze_initiation_propagation.py in
results/nohydro_gap20_descriptive_transition/ of BEM-bichamber.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

N = 100
DT = 0.001
GAP = 2.8635642126552705
ELL = GAP
Q_EDGES = np.array([0., 0.5, 0.8, 1.0000001])
A_EDGES = np.array([-1.0000001, -0.3, 0.3, 1.0000001])
LEADS = (-5., -4., -3., -2., -1., -0.5, -0.2, 0., 0.2, 0.5, 1., 2., 3., 5.)
DEFAULT_RAW = Path('/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data')
HERE = Path(__file__).resolve().parent
EVENT_FILE = HERE / 'outputs/initiation_propagation/crossing_episodes.csv'
CONTROL_FILE = HERE / 'outputs/initiation_propagation/matched_precursor_pairs.csv'


def load_csv(path: Path) -> list[dict]:
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def save_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cols, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def raw_file(raw_dir: Path, seed: int) -> Path:
    return raw_dir / f'If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat'


def nr_count(st: np.ndarray) -> int:
    return int(np.count_nonzero(st[0] >= 0))


def features(st: np.ndarray, source_sign: int, ell: float = ELL, width: float = GAP) -> dict:
    """Source-fixed geometrical coordinates; sign(q)>0 denotes right source.

    x toward opening is xi = sign(q)*x_lab in the majority source, so
    0 < xi <= ell defines its approach layer and the progress to the opening
    is (ell-xi)/ell.  This is W/(N ell) from the previous conservation audit.
    COM-to-gap polarization uses the same physical definition as the earlier
    q+a analysis.  If no source fish remain, a is undefined (not zero).
    """
    st = np.asarray(st)
    if st.shape != (3, N):
        raise ValueError(f'Expected (3,{N}) state, got {st.shape}')
    if source_sign not in (-1, 1):
        raise ValueError('source_sign must be -1 or +1')
    x, y, theta = st
    q = 2. * nr_count(st) / N - 1.
    u = source_sign * x  # source is on side u>=0, geometric distance to x=0
    source = (x >= 0) if source_sign > 0 else (x < 0)
    approach = source & (u <= ell) & (np.abs(y) <= width / 2.)
    prog = np.clip((ell - u) / ell, 0., 1.)
    b = float(np.sum(prog[approach]) / N)
    n_approach = int(np.count_nonzero(approach))
    a = float('nan')
    circulation = float('nan')
    if np.any(source):
        com = st[:2, source].mean(axis=1)
        radius = float(np.linalg.norm(com))
        if radius > 1e-12:
            gap_direction = -com / radius
            h = np.array([np.mean(np.cos(theta[source])), np.mean(np.sin(theta[source]))])
            a = float(h @ gap_direction)
            # Angular component relative to the source school COM; a possible
            # omitted variable for diagnosing memory / circulation.
            circulation = float(gap_direction[0] * h[1] - gap_direction[1] * h[0])
    j_heading = float(-source_sign * np.sum(np.cos(theta[approach])) / N)
    # -source_sign is the x-direction from source chamber toward central gap.
    out = dict(q=q, q_source=q * source_sign, a=a, b=b,
               n_source=int(np.sum(source)), n_approach=n_approach,
               j_heading=j_heading, tangent_p=circulation)
    if not (0. <= b <= 1. + 1e-12):
        raise AssertionError('Weighted storage out of range')
    return out


def boot_seed_means(rows: list[dict], value: str, nboot: int = 1500, seed: int = 20261118) -> dict:
    by = defaultdict(list)
    for r in rows:
        v = float(r[value])
        if np.isfinite(v):
            by[int(r['seed'])].append(v)
    arr = np.array([np.mean(v) for v in by.values()])
    if not len(arr):
        return dict(mean=None, ci95=[None, None], seeds=0)
    rng = np.random.default_rng(seed)
    means = arr[rng.integers(0, len(arr), (nboot, len(arr)))].mean(axis=1)
    return dict(mean=float(arr.mean()), ci95=np.percentile(means, [2.5, 97.5]).tolist(), seeds=len(arr))


def source_events(events: list[dict]) -> list[dict]:
    """First crossing originates in the pre-event *majority*, not future label."""
    primary = []
    for r in events:
        q = float(r['q_pre'])
        if abs(q) < 1e-12:
            continue
        source_sign = 1 if q > 0 else -1
        if int(r['direction']) != -source_sign:
            continue
        primary.append(dict(seed=int(r['seed']), event_id=int(r['event_id']),
                            step=int(r['first_native']), source_sign=source_sign,
                            q_pre=q, pre_quiet=float(r.get('pre_quiet_duration', 'inf'))))
    return primary


def extract_event_paths(readers: dict, events: list[dict], horizon_steps: int,
                        leads: tuple[float, ...]) -> tuple[list[dict], list[dict]]:
    paths, events_out = [], []
    for ev in events:
        rd = readers[ev['seed']]
        k, sign = ev['step'], ev['source_sign']
        total = rd.info['state'].shape[-1]
        if k < 1 or k + horizon_steps >= total:
            continue
        pre = rd.frame('state', k-1)
        at = rd.frame('state', k)
        end = rd.frame('state', k+horizon_steps)
        if abs(features(pre, sign)['q'] - ev['q_pre']) > 1e-10:
            raise RuntimeError(f'Indexed q mismatch: seed={ev["seed"]}, event={ev["event_id"]}')
        if -sign * (nr_count(at) - nr_count(pre)) <= 0:
            raise RuntimeError(f'Indexed first crossing has wrong direction: {ev}')
        additional = -sign * (nr_count(end) - nr_count(at))
        success = int(additional >= 4)
        metadata = dict(seed=ev['seed'], event_id=ev['event_id'], source_sign=sign,
                        first_native=k, additional_fish=additional, success4=success,
                        initial_fish=-sign*(nr_count(at)-nr_count(pre)),
                        abs_q_pre=abs(ev['q_pre']))
        events_out.append(metadata)
        for tau in leads:
            idx = k + int(round(tau / DT))
            if not 0 <= idx < total:
                continue
            row = dict(metadata, tau=tau, t=idx*DT)
            row.update(features(rd.frame('state', idx), sign))
            paths.append(row)
    return paths, events_out


def matched_controls(readers: dict, pairs: list[dict], valid_events: set,
                     event_steps: dict, leads: tuple[float, ...] = (2., 1.)) -> list[dict]:
    rows = []
    seen = set()
    for pair in pairs:
        seed, eid = int(pair['event_seed']), int(pair['event_id'])
        lead = float(pair['lead'])
        if (seed, eid) not in valid_events or lead not in leads:
            continue
        key = (seed, eid, lead)
        if key in seen:
            raise RuntimeError(f'Duplicate matched pair: {key}')
        seen.add(key)
        q = float(pair['q'])
        sign = 1 if q > 0 else -1
        ev_step = event_steps[(seed, eid)] - int(round(lead / DT))
        ctrl_seed = int(pair['control_seed'])
        ctrl_step = int(round(float(pair['control_time'])/DT))
        ev = features(readers[seed].frame('state', ev_step), sign)
        ct = features(readers[ctrl_seed].frame('state', ctrl_step), sign)
        if abs(ev['q']-q)>1e-10 or abs(ct['q']-q)>1e-10:
            raise RuntimeError(f'Exact-q event/control match fails at {key}: {ev["q"]}, {ct["q"]}, {q}')
        rows.append(dict(seed=seed, event_id=eid, lead=lead, control_seed=ctrl_seed,
                         q=q, a_event=ev['a'], a_control=ct['a'],
                         b_event=ev['b'], b_control=ct['b'],
                         db=ev['b']-ct['b'], da=ev['a']-ct['a']))
    return rows



def prospective_grid(readers: dict, tmin: float, tmax: float,
                     grid_dt: float, future: float) -> list[dict]:
    """Fixed pre-specified times; all labels use *future* raw population only."""
    if min(tmin, grid_dt, future) <= 0:
        raise ValueError('Require positive stationary start, grid spacing, horizon')
    steps = np.arange(round(tmin/DT), round(tmax/DT)+1, round(grid_dt/DT), dtype=int)
    future_step = int(round(future / DT))
    previous_step = future_step
    rows = []
    for seed, rd in sorted(readers.items()):
        limit = rd.info['state'].shape[-1]
        for k in steps:
            if k < previous_step or k+future_step >= limit:
                continue
            curr = rd.frame('state', int(k))
            q = 2*nr_count(curr)/N - 1
            if abs(q)<1e-12:
                continue  # explicitly exclude ties; do not invent a direction
            sign = 1 if q > 0 else -1
            f = features(curr, sign)
            before = features(rd.frame('state', int(k-previous_step)), sign)
            after = rd.frame('state', int(k+future_step))
            net = -sign * (nr_count(after) - nr_count(curr))
            rows.append(dict(seed=seed, t=float(k*DT), step=int(k), source_sign=sign,
                             future_net=int(net), future4=int(net >= 4),
                             b_prev=before['b'], a_prev=before['a'],
                             **f))
    return rows


def bin_id(x: float, edges: np.ndarray) -> int | None:
    if not np.isfinite(x):
        return None
    idx = int(np.searchsorted(edges, x, side='right')-1)
    return idx if 0 <= idx < len(edges)-1 else None


def qa_key(row: dict) -> tuple[int, int] | None:
    q = bin_id(abs(float(row['q'])), Q_EDGES)
    a = bin_id(float(row['a']), A_EDGES)
    return (q,a) if q is not None and a is not None else None


def fit_diagnostic(train: list[dict], min_n: int = 24,
                   min_seeds: int = 4) -> tuple[dict, list[dict]]:
    groups = defaultdict(list)
    for r in train:
        k=qa_key(r)
        if k is not None:
            groups[k].append(r)
    global_y = np.mean([r['future4'] for r in train])
    global_n = np.mean([r['future_net'] for r in train])
    model, support = {}, []
    for key,rs in groups.items():
        n=len(rs); ys=np.asarray([r['future4'] for r in rs],float)
        nets=np.asarray([r['future_net'] for r in rs],float)
        base_p=float((ys.sum()+.5)/(n+1))
        base_net=float(nets.mean())
        bs=np.asarray([r['b'] for r in rs],float)
        cut=float(np.median(bs))  # training-only unsupervised median
        low=[r for r in rs if r['b']<=cut]
        high=[r for r in rs if r['b']>cut]
        def okay(xs):
            return len(xs)>=min_n and len({int(r['seed']) for r in xs})>=min_seeds
        supported=bool(okay(low) and okay(high) and high and np.max(bs)>np.min(bs))
        splits = {}
        for label,sub in ((0,low),(1,high)):
            if supported:
                ny=np.asarray([r['future4'] for r in sub],float)
                nn=np.asarray([r['future_net'] for r in sub],float)
                # Fixed small shrinkage avoids extreme probabilities on thin bins;
                # this is a benchmark, NOT a phenomenological flux law.
                splits[label]=dict(p=float((ny.sum()+1.*base_p)/(len(sub)+1.)),
                                   net=float((nn.sum()+5.*base_net)/(len(sub)+5.)))
        model[key]=dict(count=n, p=base_p, net=base_net, cut=cut,
                        supported=supported, split=splits)
        support.append(dict(iq=key[0], ia=key[1], count=n,
                            n_seeds=len({r['seed'] for r in rs}),
                            b_median=cut, n_low=len(low), n_high=len(high),
                            n_low_seeds=len({r['seed'] for r in low}),
                            n_high_seeds=len({r['seed'] for r in high}),
                            supported=int(supported),
                            mean_net_low=float(np.mean([r['future_net'] for r in low])) if low else None,
                            mean_net_high=float(np.mean([r['future_net'] for r in high])) if high else None))
    return dict(cells=model, global_p=float(global_y), global_net=float(global_n)), support


def eval_diagnostic(test: list[dict], model: dict) -> tuple[list[dict],dict]:
    preds=[]
    for row in test:
        key=qa_key(row)
        if key is None: continue
        fitted=model['cells'].get(key)
        base_p=model['global_p'] if fitted is None else fitted['p']
        base_n=model['global_net'] if fitted is None else fitted['net']
        supported=bool(fitted and fitted['supported'])
        split=(int(row['b']>fitted['cut']) if supported else -1)
        aug=fitted['split'][split] if supported else dict(p=base_p,net=base_n)
        preds.append(dict(seed=int(row['seed']),t=row['t'], q=row['q'],a=row['a'],b=row['b'],
                          b_prev=row['b_prev'],a_prev=row['a_prev'],
                          future_net=row['future_net'],future4=row['future4'],
                          iq=key[0],ia=key[1], b_bin=split,supported=int(supported),
                          p_qa=base_p,p_qab=aug['p'],net_qa=base_n,net_qab=aug['net']))
    by=defaultdict(list)
    for r in preds: by[r['seed']].append(r)
    seed_scores=[]
    for seed,rs in sorted(by.items()):
        yt=np.array([r['future4'] for r in rs],float)
        yn=np.array([r['future_net'] for r in rs],float)
        p0=np.array([r['p_qa'] for r in rs]);p1=np.array([r['p_qab'] for r in rs])
        n0=np.array([r['net_qa'] for r in rs]);n1=np.array([r['net_qab'] for r in rs])
        seed_scores.append(dict(seed=seed,n=len(rs),n_success=int(yt.sum()),
                                n_supported=int(sum(r['supported'] for r in rs)),
                                brier_qa=float(np.mean((yt-p0)**2)),
                                brier_qab=float(np.mean((yt-p1)**2)),
                                mse_net_qa=float(np.mean((yn-n0)**2)),
                                mse_net_qab=float(np.mean((yn-n1)**2)),
                                gain_brier=float(np.mean((yt-p0)**2-(yt-p1)**2)),
                                gain_mse_net=float(np.mean((yn-n0)**2-(yn-n1)**2))))
    s=dict(heldout_n=len(preds),n_seeds=len(seed_scores),
           heldout_supported=int(sum(r['supported'] for r in preds)),
           brier_improvement=boot_seed_means(seed_scores,'gain_brier'),
           net_mse_improvement=boot_seed_means(seed_scores,'gain_mse_net'),
           seed_scores=seed_scores)
    if preds:
        y=np.array([r['future4'] for r in preds],float)
        n=np.array([r['future_net'] for r in preds],float)
        s.update(brier_qa=float(np.mean([(r['p_qa']-r['future4'])**2 for r in preds])),
                 brier_qab=float(np.mean([(r['p_qab']-r['future4'])**2 for r in preds])),
                 mse_net_qa=float(np.mean([(r['net_qa']-r['future_net'])**2 for r in preds])),
                 mse_net_qab=float(np.mean([(r['net_qab']-r['future_net'])**2 for r in preds])),
                 successes=int(y.sum()))
    return preds,s


def memory_audit(preds: list[dict], boot=1200) -> dict:
    """Lagged state vs held-out residual: a necessary, not sufficient closure check."""
    rows=[r for r in preds if r['supported'] and np.isfinite(r['a_prev'])]
    if len(rows)<30:
        return dict(available=False,reason='Insufficient supported heldout states')
    def correlation(rr, history):
        x=np.array([r[history] for r in rr],float)
        res=np.array([r['future_net']-r['net_qab'] for r in rr],float)
        if np.std(x)<1e-12 or np.std(res)<1e-12: return float('nan')
        return float(spearmanr(x,res).statistic)
    seeds=sorted({r['seed'] for r in rows})
    out=dict(available=True,n=len(rows),n_seeds=len(seeds))
    rng=np.random.default_rng(20261118)
    for name in ['b_prev','a_prev']:
        obs=correlation(rows,name)
        draws=[]
        by={seed:[r for r in rows if r['seed']==seed] for seed in seeds}
        for _ in range(boot):
            rr=[]
            for seed in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[seed])
            corr=correlation(rr,name)
            if np.isfinite(corr): draws.append(corr)
        out[name]=dict(spearman_heldout_residual=obs if np.isfinite(obs) else None,
                       seed_bootstrap_ci95=np.percentile(draws,[2.5,97.5]).tolist() if draws else None)
    out['caveat']='Residual association can also arise from coarse binning; a null is not proof of Markov closure.'
    return out


def draw_path(paths: list[dict], output: Path):
    if not paths:return
    fig,ax=plt.subplots(3,1,figsize=(8.2,9),sharex=True)
    for i,(name,title) in enumerate([('a','Gap-directed source polarization a'),
                                     ('b','Weighted approach storage b'),
                                     ('q_source','Imbalance relative to original source')]):
        for success,color,label in [(0,'#C07950','<4 additional fish'),(1,'#357A99','>=4 additional fish')]:
            time=sorted({r['tau'] for r in paths})
            m=[];low=[];high=[]
            for tau in time:
                group=[r for r in paths if r['success4']==success and r['tau']==tau]
                d=boot_seed_means(group,name,nboot=600,seed=20261118+int(10*tau)+100)
                m.append(d['mean'] if d['mean'] is not None else np.nan)
                lo,hi=d['ci95'];low.append(lo if lo is not None else np.nan);high.append(hi if hi is not None else np.nan)
            ax[i].plot(time,m,'o-',label=label,color=color,markersize=3)
            ax[i].fill_between(time,low,high,color=color,alpha=.18,linewidth=0)
        ax[i].axvline(0,color='0.35',lw=1,ls='--')
        ax[i].set_ylabel(title);ax[i].grid(alpha=.2)
    ax[-1].set_xlabel('Time relative to first crossing (native units)')
    ax[0].legend(frameon=False)
    fig.suptitle('Event-aligned trajectories: source chamber fixed at t=0-')
    fig.tight_layout();fig.savefig(output,dpi=180);plt.close(fig)


def draw_stratified_path(paths: list[dict], output: Path):
    """Prevent success/failure differences being driven only by initial loading."""
    fig, axes = plt.subplots(2, 3, figsize=(12.8, 6.8), sharex=True)
    bands = ((0.,.5), (.5,.8), (.8,1.00001))
    for col, (lo,hi) in enumerate(bands):
        subset = [r for r in paths if lo<=r['abs_q_pre']<hi]
        for row, (var,ylabel) in enumerate((('a','source polarization a'),('b','weighted storage b'))):
            ax=axes[row,col]
            for status,color in ((0,'#C07950'),(1,'#357A99')):
                part=[r for r in subset if r['success4']==status]
                seeds={r['seed'] for r in part}
                events={(r['seed'],r['event_id']) for r in part}
                if len(events)<5 or len(seeds)<3: continue
                time=sorted({r['tau'] for r in part})
                means=[boot_seed_means([r for r in part if r['tau']==t],var,nboot=350)['mean'] for t in time]
                ax.plot(time,means,'o-',markersize=2.5,color=color,
                        label=f'{"continue" if status else "failed"}: {len(events)} events')
            ax.set_title(f'{lo:g} <= |q(0-)| < {hi:g}')
            ax.axvline(0,color='0.45',linestyle='--',lw=1)
            ax.grid(alpha=.18)
            if col==0:ax.set_ylabel(ylabel)
            if row==1:ax.set_xlabel('Relative time')
            ax.legend(frameon=False,fontsize=7)
    fig.suptitle('Event-aligned pathways stratified by pre-entry imbalance')
    fig.tight_layout();fig.savefig(output,dpi=180);plt.close(fig)

def draw_binned(support: list[dict], out: Path):
    valid=[r for r in support if r['supported']]
    fig,ax=plt.subplots(figsize=(7,4))
    x=np.arange(len(valid))
    if valid:
        ax.plot(x,[r['mean_net_low'] for r in valid],'o',label='Below median b')
        ax.plot(x,[r['mean_net_high'] for r in valid],'s',label='Above median b')
        ax.set_xticks(x,[f'q{r["iq"]}/a{r["ia"]}' for r in valid],rotation=45)
        ax.legend(frameon=False)
    ax.set_ylabel('TRAINING mean future net fish (+ toward balance)')
    ax.set_title('Conditional storage contrast, descriptive only')
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out,dpi=180);plt.close(fig)


def report(summary: dict) -> str:
    ep=summary['events'];ct=summary['matched_controls'];d=summary['heldout'];mem=summary['memory']
    g=d['brier_improvement'];f=d['net_mse_improvement']
    def fmt(n):return 'NA' if n is None else f'{n:.5g}'
    lines=[
        '# Collective pathway and (q,a,b) closure diagnostic', '',
        '## Scope',
        'Read-only analysis of existing N=100, gap-20, no-hydrodynamics trajectories. No new microscopic simulation, fitted potential, fitted drift, or causal intervention.',
        f"{ep['primary_events']} fixed-common-cohort majority-to-minority first crossings; {ep['continued4']} transferred at least four ADDITIONAL fish during the next five units.",
        f"Independent original seeds represented: {ep['seeds']}. Event-label groups are retrospective, not an online classifier.", '',
        '## Physical coordinates',
        '- `q=(N_R-N_L)/N` is signed population imbalance. Source side is fixed from q immediately before the first crossing (or from current q on prospective time grids); ties are excluded.',
        '- `a=P_source . (-COM_source/|COM_source|)` is the gap-directed source polarization. When no source fish remain or its COM is undefined, a is missing.',
        '- `b=(1/(N*ell)) sum_{source approach} (ell-u_i)`, with `u_i=source_sign*x_i` and `0<=u_i<=ell`, `|y_i|<=gap/2`. This equals the weighted-storage observable W/(N*ell) from the prior conservation identity.',
        '- Event figures freeze the ORIGINAL source chamber so q, a, and b remain interpretable through changes of majority.', '',
        '## Event-aligned dynamics',
        'The figure compares continued and failed crossing episodes at fixed offsets on one common cohort: no prior crossing for at least five units and post-transient throughout. Uncertainty is clustered by original seed.',
        'Chronology must be read from the curves: the code deliberately does not assert that a changes before b.', '',
        '## Exact-q matched controls',
        f"Matched measurements: {ct['n']}, at leads 2 and 1. Controls were selected using stable future no-flux conditions: retrospective case-control comparisons, NOT an unbiased crossing hazard.",
    ]
    for k in ('1.0','2.0'):
        x=ct['by_lead'].get(k)
        if x:
            v=x['db'];lines.append(f"- Lead {k}: b(event)-b(control)={fmt(v['mean'])} (seed CI [{fmt(v['ci95'][0])}, {fmt(v['ci95'][1])}]); {v['seeds']} event seeds.")
    lines += ['', '## Prospective held-out diagnostic',
             'Fixed time grid, no use of future event direction; seeds 0-29 are training and 30-39 are held out. Within each predeclared q/a cell, training b is split at its median; cells need adequate sample and independent-seed support. Unsupported cells fall back to q+a.',
             f"Held-out samples {d['heldout_n']}, supported q/a/b cells {d['heldout_supported']}; supported coverage {fmt(d['heldout_supported']/max(d['heldout_n'],1))}.",
             f"P(>=4 future net fish) Brier: q+a {fmt(d.get('brier_qa'))}, q+a+b {fmt(d.get('brier_qab'))}.",
             f"Seed-mean improvement (positive favors b): {fmt(g['mean'])}, CI [{fmt(g['ci95'][0])}, {fmt(g['ci95'][1])}].",
             f"Future-net-fish MSE: q+a {fmt(d.get('mse_net_qa'))}, q+a+b {fmt(d.get('mse_net_qab'))}.",
             f"Seed-mean net-MSE improvement: {fmt(f['mean'])}, CI [{fmt(f['ci95'][0])}, {fmt(f['ci95'][1])}].",
             'These are nonparametric diagnostic conditional means, NOT a physically derived closure.', '',
             '## Residual memory',
             f"Supported held-out states: {mem.get('n','NA')}. Association between lagged a/b and future-net residual is a necessary closure check, not a sufficiency proof.",
    ]
    for k in ('a_prev','b_prev'):
        if mem.get(k):
            z=mem[k]; lines.append(f"- {k}: Spearman {fmt(z['spearman_heldout_residual'])}, seed CI {z['seed_bootstrap_ci95']}.")
    lines += ['', '## Interpretation rule',
              'Only describe an observed reorientation-accumulation-release sequence if the event-aligned ordering, control contrast and held-out storage test support it jointly. Even then, chronology alone is not causal evidence. If storage adds no held-out information, reject (q,a,b) as an adequate predictive state and revisit spatial/circulation coordinates. Full opposite-chamber relocation is not measured here.', '',
              '## Reproducibility',
              'Outputs include `event_paths.csv`, `event_outcomes.csv`, `matched_b_controls.csv`, `prospective_qab.csv`, `heldout_predictions.csv`, `qa_b_support.csv`, `summary.json`, `figures/` (including initial-imbalance-stratified pathways). Original MAT data and source files remain unchanged.']
    return '\n'.join(lines)+'\n'


def run(args):
    from analyze_initiation_propagation import RawReader
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    ev=source_events(load_csv(args.events))
    if not ev:raise RuntimeError('No primary events selected')
    pair=load_csv(args.controls)
    eligible={(e['seed'],e['event_id']) for e in ev
              if e['seed'] in args.seeds and e['pre_quiet']>=abs(min(LEADS))
              and e['step']>=round((args.transient+abs(min(LEADS)))/DT)}
    needed=set(args.seeds)
    if not args.skip_controls:
        for p in pair:
            if (int(p['event_seed']),int(p['event_id'])) in eligible and float(p['lead']) in (1.,2.):
                needed.add(int(p['control_seed']))
    readers={}
    try:
        for seed in sorted(needed):
            path=raw_file(args.raw_dir,seed)
            if not path.is_file():raise FileNotFoundError(f'Missing archived trajectory {path}')
            readers[seed]=RawReader(path)
        # Fixed common cohort with >=5 units free of prior crossings; prevents
        # predecessor episodes from masquerading as pre-entry reorientation.
        primary=[e for e in ev if e['seed'] in args.seeds
                 and e['pre_quiet'] >= abs(min(LEADS))
                 and e['step'] >= round((args.transient+abs(min(LEADS)))/DT)]
        paths, outcomes=extract_event_paths(readers,primary,int(round(args.horizon/DT)),LEADS)
        if not outcomes:
            raise RuntimeError('No common-cohort events; check lead and quiet criteria')
        valid={(r['seed'],r['event_id']) for r in outcomes}
        steps={(e['seed'],e['event_id']):e['step'] for e in primary}
        ctrl=matched_controls(readers,pair,valid,steps) if not args.skip_controls else []
        grid=prospective_grid({seed:readers[seed] for seed in args.seeds},args.transient,
                              args.last_grid,args.grid_step,args.horizon)
    finally:
        for rd in readers.values():rd.close()
    save_csv(out/'event_paths.csv',paths)
    save_csv(out/'event_outcomes.csv',outcomes)
    save_csv(out/'matched_b_controls.csv',ctrl)
    save_csv(out/'prospective_qab.csv',grid)
    model,support=fit_diagnostic([r for r in grid if r['seed']<args.train_end],
                                 min_n=args.min_cell_n,min_seeds=args.min_cell_seeds)
    preds,evalsum=eval_diagnostic([r for r in grid if r['seed']>=args.train_end],model)
    save_csv(out/'qa_b_support.csv',support)
    save_csv(out/'heldout_predictions.csv',preds)
    fig=out/'figures';fig.mkdir(exist_ok=True)
    draw_path(paths,fig/'event_pathway.png')
    draw_stratified_path(paths,fig/'event_pathway_by_imbalance.png')
    draw_binned(support,fig/'train_storage_contrast.png')
    by_lead={}
    for lead in sorted({r['lead'] for r in ctrl}):
        c=[r for r in ctrl if r['lead']==lead]
        by_lead[str(lead)]={'db':boot_seed_means(c,'db'), 'da':boot_seed_means(c,'da')}
    result=dict(configuration=dict(raw_dir=str(args.raw_dir),original_N=N,gap=GAP,ell=ELL,
                                   seeds=args.seeds,dt=DT,event_horizon=args.horizon,
                                   leads=list(LEADS),grid_step=args.grid_step,
                                   prospective_range=[args.transient,args.last_grid],
                                   train_end=args.train_end,min_cell_n=args.min_cell_n,
                                   min_cell_seeds=args.min_cell_seeds),
                events=dict(primary_events=len(outcomes),continued4=sum(r['success4'] for r in outcomes),
                            seeds=len({r['seed'] for r in outcomes})),
                matched_controls=dict(n=len(ctrl),by_lead=by_lead),
                prospective_n=len(grid), heldout=evalsum,memory=memory_audit(preds))
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (out/'REPORT.md').write_text(report(result))
    print(json.dumps({'events':result['events'],'matched_pairs':len(ctrl),
                      'prospective_n':len(grid),'heldout':{k:v for k,v in evalsum.items() if k!='seed_scores'},
                      'memory':result['memory']},indent=2))


def self_test():
    rng=np.random.default_rng(9)
    st=np.zeros((3,N));st[0,:70]=5.;st[0,70:]=-5.;st[1,:]=0.;st[2,:]=np.pi
    f=features(st,1)
    assert np.isclose(f['q'],.4) and np.isclose(f['a'],1) and f['b']==0
    st[0,0:10]=.5
    f2=features(st,1)
    assert np.isclose(f2['b'],10*(ELL-.5)/ELL/N) and f2['n_approach']==10
    st2=st.copy();st2[0,:]=-st[0];st2[2,:]=0.
    g=features(st2,-1)
    assert np.isclose(g['b'],f2['b']) and np.isclose(g['q'],-f2['q'])
    # Prospective checks: gain is not guaranteed, but positive-storage labels
    # should be extractable, with no future labels used for choosing b splits.
    rows=[]
    for seed in range(40):
        for t in range(50):
            q=.2 if t%2 else .65
            a=.1 if t%3 else .6
            b=(.005 if t%2 else .15)+(.001*(seed%3))
            y=int(rng.random()<(.15+.5*(b>.05)))
            rows.append(dict(seed=seed,t=t,q=q,a=a,b=b,b_prev=b*.8,a_prev=a,
                             future4=y,future_net=int(y*4)))
    m,s=fit_diagnostic([r for r in rows if r['seed']<30],min_n=5,min_seeds=4)
    preds,res=eval_diagnostic([r for r in rows if r['seed']>=30],m)
    assert len(preds)==500 and res['n_seeds']==10
    assert any(r['supported'] for r in s)
    assert all(0<=p['p_qab']<=1 for p in preds)
    assert source_events([{'seed':'3','event_id':'4','first_native':'7000','q_pre':'-0.6','direction':'1'},
                          {'seed':'3','event_id':'5','first_native':'8000','q_pre':'0.6','direction':'1'}])==[
                            dict(seed=3,event_id=4,step=7000,source_sign=-1,q_pre=-.6,pre_quiet=float('inf'))]
    print('PASS: weighted storage; left/right symmetry; event direction; train-only b bins; heldout diagnostics')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,default=DEFAULT_RAW)
    parser.add_argument('--events',type=Path,default=EVENT_FILE)
    parser.add_argument('--controls',type=Path,default=CONTROL_FILE)
    parser.add_argument('--output',type=Path,default=HERE/'outputs/collective_qab_pathway')
    parser.add_argument('--seeds',type=int,nargs='*',default=list(range(40)))
    parser.add_argument('--horizon',type=float,default=5.)
    parser.add_argument('--transient',type=float,default=200.)
    parser.add_argument('--last-grid',type=float,default=990.)
    parser.add_argument('--grid-step',type=float,default=5.)
    parser.add_argument('--train-end',type=int,default=30)
    parser.add_argument('--min-cell-n',type=int,default=24)
    parser.add_argument('--min-cell-seeds',type=int,default=4)
    parser.add_argument('--skip-controls',action='store_true')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if abs(args.horizon/DT-round(args.horizon/DT))>1e-9:
            raise ValueError('horizon must be multiple of native dt')
        if not args.seeds:raise ValueError('No seeds selected')
        run(args)


if __name__=='__main__':main()
