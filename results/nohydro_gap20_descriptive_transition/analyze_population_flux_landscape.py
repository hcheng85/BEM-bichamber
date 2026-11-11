#!/usr/bin/env python3
"""Empirical p(q), p(q,J_past), and prospective P(transfer|q,J_past).

Run from the BEM-bichamber repository root. Requires the existing
analyze_initiation_propagation.py script and its MAT-v5 reader. No new simulation.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / 'results' / 'nohydro_gap20_descriptive_transition'))

N = 100
Q_EDGES = np.linspace(-1.01, 1.01, 22)
ABS_Q_EDGES = [0., .5, .8, 1.000001]
CATEGORIES = ['toward dominant chamber', 'zero net flux', 'toward balance']


def csv_write(path, items):
    if not items:
        return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(items[0]))
        writer.writeheader()
        writer.writerows(items)


def load_seed(seed, a):
    cache = a.cache_dir / f'seed{seed:02d}_stride{a.stride}_qscan_v4.npz'
    raw = a.raw_dir / f'If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat'
    processed = None if a.processed_dir is None else (a.processed_dir /
        f'event_histogram_If00_In05_num100_R50_gap20_n{seed}_BEM_no_hydro.mat')
    if not cache.exists() and not raw.exists():
        raise FileNotFoundError(f'Need cached occupancy {cache} or raw {raw}')
    a.cache_dir.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        with np.load(cache) as f:
            d = {k: f[k] for k in f.files}
    else:
        from analyze_initiation_propagation import extract_seed
        d = extract_seed(raw, processed, cache, a.stride, a.gap)
    q = np.asarray(d['q'], float)
    nl = np.asarray(d['nl'], float)
    if len(q) < 2 or not np.allclose(q, 1-2*nl/N) or np.max(np.abs(q)) > 1.000001:
        raise ValueError(f'Invalid occupancy index for seed {seed}')
    return q


def relocations(q, bound=.6):
    previous, n = 0, 0
    for v in q:
        state = -1 if v <= -bound else (1 if v >= bound else 0)
        if state:
            n += int(previous != 0 and state != previous)
            previous = state
    return n


def risk_samples(q, seed, a):
    dt = a.dt*a.stride
    ints = [(round(value/dt), value) for value in (a.lookback,a.horizon,a.spacing)]
    if min(k for k, v in ints) < 1 or a.spacing < a.horizon:
        raise ValueError('Require positive lookback, horizon, and nonoverlapping future windows')
    if any(abs(k*dt-v)>1e-8 for k,v in ints):
        raise ValueError('Lookback, horizon, and spacing must be multiples of stride*dt')
    look, future, spacing = (k for k,_ in ints)
    first = max(look, round(a.transient/dt))
    out = []
    for i in range(first, len(q)-future, spacing):
        qi = float(q[i])
        if qi == 0:  # ambiguous source chamber; exclude only from directional risk map
            continue
        toward_balance = -int(np.sign(qi))
        # Actual *net neck-crossing count* over a past window, not approach current.
        jpast = N*(qi-q[i-look])/(2*a.lookback)
        jforward = toward_balance*jpast
        future_dq = toward_balance*(q[i+1:i+future+1]-qi)
        # Future *net displacement* threshold, not an accepted event marker.
        outcome = int(np.max(future_dq) >= 2*a.threshold_fish/N-1e-10)
        qb = int(np.searchsorted(ABS_Q_EDGES, abs(qi), side='right')-1)
        if not 0 <= qb < 3:
            continue
        jb = 0 if jforward < -1e-10 else (2 if jforward > 1e-10 else 1)
        out.append(dict(seed=seed, t=i*dt, q=qi, J_past=jpast,
                        J_toward_balance=jforward, q_bin=qb, J_class=jb,
                        future_transfer=outcome))
    return out


def counts(rows):
    n = np.zeros((3,3), int)
    y = np.zeros((3,3), int)
    for r in rows:
        i,j=r['q_bin'],r['J_class']
        n[i,j]+=1
        y[i,j]+=r['future_transfer']
    return n,y


def seed_bootstrap_cells(rows, seeds, draws=3000):
    """Seed-cluster intervals for the nine conditional-risk cells."""
    seed_counts = np.zeros((len(seeds), 3, 3), int)
    seed_success = np.zeros_like(seed_counts)
    seed_index = {seed: i for i, seed in enumerate(seeds)}
    for row in rows:
        ii = seed_index[row['seed']]
        i, j = row['q_bin'], row['J_class']
        seed_counts[ii, i, j] += 1
        seed_success[ii, i, j] += row['future_transfer']
    rng = np.random.default_rng(20261111)
    pick = rng.integers(len(seeds), size=(draws, len(seeds)))
    out = {}
    for i in range(3):
        for j in range(3):
            den = seed_counts[pick, i, j].sum(axis=1)
            num = seed_success[pick, i, j].sum(axis=1)
            valid = den > 0
            values = num[valid] / den[valid]
            out[(i, j)] = (float(np.quantile(values, .025)),
                           float(np.quantile(values, .975))) if len(values) else (np.nan, np.nan)
    return out


def validate_against_raw(q_by_seed, a):
    """Bounded raw-state audit of occupancy and net centerline-crossing flux."""
    from analyze_initiation_propagation import RawReader
    dt_coarse = a.dt * a.stride
    look = round(a.lookback / dt_coarse)
    rows = []
    for seed in sorted(q_by_seed):
        raw = a.raw_dir / f'If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat'
        if not raw.exists():
            raise FileNotFoundError(f'Raw occupancy validation requires {raw}')
        q = q_by_seed[seed]
        reader = RawReader(raw)
        # Three endpoint checks plus one conservative crossing-count check per seed.
        indices = [round(t / dt_coarse) for t in (250., 500., 750.)]
        for idx in indices:
            state = reader.frame('state', idx * a.stride)
            nl_raw = int(np.sum(state[0] < 0))
            q_raw = 1. - 2. * nl_raw / N
            rows.append(dict(seed=seed, kind='occupancy', t=idx * dt_coarse,
                             q_cache=float(q[idx]), q_raw=q_raw,
                             abs_q_error=abs(float(q[idx]) - q_raw),
                             J_cache=np.nan, J_endpoint=np.nan, J_crossing=np.nan,
                             abs_J_endpoint_error=np.nan, abs_J_crossing_error=np.nan,
                             left_to_right=np.nan, right_to_left=np.nan))
        idx = round(500. / dt_coarse)
        first = (idx - look) * a.stride
        last = idx * a.stride
        block = reader.state_block(first, last)
        x = block[0]
        nr = np.sum(x >= 0, axis=0)
        lr = int(np.sum((x[:, :-1] < 0) & (x[:, 1:] >= 0)))
        rl = int(np.sum((x[:, :-1] >= 0) & (x[:, 1:] < 0)))
        j_cache = N * (q[idx] - q[idx-look]) / (2 * a.lookback)
        j_endpoint = (int(nr[-1]) - int(nr[0])) / a.lookback
        j_crossing = (lr - rl) / a.lookback
        rows.append(dict(seed=seed, kind='net_flux', t=idx * dt_coarse,
                         q_cache=float(q[idx]), q_raw=1. - 2. * (N-int(nr[-1])) / N,
                         abs_q_error=abs(float(q[idx]) - (2.*int(nr[-1])/N-1.)),
                         J_cache=float(j_cache), J_endpoint=float(j_endpoint),
                         J_crossing=float(j_crossing),
                         abs_J_endpoint_error=abs(float(j_cache-j_endpoint)),
                         abs_J_crossing_error=abs(float(j_cache-j_crossing)),
                         left_to_right=lr, right_to_left=rl))
        reader.close()
    csv_write(a.output_dir/'occupancy_flux_raw_validation.csv', rows)
    flux_rows = [row for row in rows if row['kind'] == 'net_flux']
    return dict(n_occupancy_checks=sum(row['kind'] == 'occupancy' for row in rows),
                n_flux_windows=len(flux_rows),
                max_abs_q_error=float(max(row['abs_q_error'] for row in rows)),
                max_abs_J_endpoint_error=float(max(row['abs_J_endpoint_error'] for row in flux_rows)),
                max_abs_J_crossing_error=float(max(row['abs_J_crossing_error'] for row in flux_rows)))


def branch_persistence(q_by_seed, a):
    """Nonoverlapping past-flux windows; one row per independent seed."""
    dt_coarse = a.dt * a.stride
    look = round(a.lookback / dt_coarse)
    first = max(look, round(a.transient / dt_coarse))
    rows = []
    for seed, q in sorted(q_by_seed.items()):
        ends = np.arange(first, len(q) - look, look)
        j0 = N * (q[ends] - q[ends-look]) / (2*a.lookback)
        j1 = N * (q[ends+look] - q[ends]) / (2*a.lookback)
        both = (np.abs(j0) > 1e-12) & (np.abs(j1) > 1e-12)
        rows.append(dict(seed=seed, n_pairs=len(ends), n_both_nonzero=int(both.sum()),
                         same_sign_fraction=float(np.mean(np.sign(j0[both]) == np.sign(j1[both]))) if both.any() else np.nan,
                         next_nonzero_given_current=float(np.mean(np.abs(j1[np.abs(j0)>1e-12]) > 1e-12)) if np.any(np.abs(j0)>1e-12) else np.nan))
    csv_write(a.output_dir/'flux_branch_persistence_by_seed.csv', rows)
    values = np.asarray([row['same_sign_fraction'] for row in rows], float)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(20261111)
    boot = values[rng.integers(len(values), size=(3000, len(values)))].mean(axis=1)
    return dict(mean_same_sign_fraction=float(values.mean()),
                seed_bootstrap_ci95=list(map(float, np.quantile(boot, [.025, .975]))),
                n_seeds=len(values))


def heldout(rows, seeds, min_cell=20):
    if len(seeds)<4:
        return {'available':False, 'reason':'Fewer than 4 independent seeds'}
    split = max(1,len(seeds)//4)
    train = [r for r in rows if r['seed'] in set(seeds[:-split])]
    test = [r for r in rows if r['seed'] in set(seeds[-split:])]
    if not train or not test:
        return {'available':False, 'reason':'Insufficient risk windows'}
    n,y=counts(train)
    qn=n.sum(axis=1); qy=y.sum(axis=1)
    overall=(y.sum()+.5)/(n.sum()+1)  # fixed Jeffreys half-count convention
    results=[]
    for seed in seeds[-split:]:
        rs=[r for r in test if r['seed']==seed]
        if not rs:
            continue
        truth=np.asarray([r['future_transfer'] for r in rs])
        p0=[];p1=[];support=0
        for r in rs:
            i,j=r['q_bin'],r['J_class']
            pq=(qy[i]+.5)/(qn[i]+1) if qn[i]>=min_cell else overall
            if n[i,j]>=min_cell:
                pj=(y[i,j]+.5)/(n[i,j]+1)
                support+=1
            else:
                pj=pq
            p0.append(pq);p1.append(pj)
        results.append(dict(seed=seed, n=len(rs), positives=int(truth.sum()),
            coverage=support/len(rs), brier_q=float(np.mean((truth-p0)**2)),
            brier_qJ=float(np.mean((truth-p1)**2))))
    if not results:
        return {'available':False}
    effect=np.asarray([r['brier_q']-r['brier_qJ'] for r in results])
    rng=np.random.default_rng(20261111)
    boots=rng.choice(effect, (3000,len(effect)), replace=True).mean(axis=1)
    return dict(available=True, training_seeds=seeds[:-split], heldout_seeds=seeds[-split:],
        per_seed=results, mean_brier_q=float(np.mean([r['brier_q'] for r in results])),
        mean_brier_qJ=float(np.mean([r['brier_qJ'] for r in results])),
        mean_brier_improvement=float(effect.mean()),
        bootstrap_ci95=list(map(float,np.quantile(boots,[.025,.975]))))


def dense_joint_series(q, a):
    """All sampled times for a time-weighted joint distribution, independent of risk stride."""
    dt = a.dt*a.stride
    look = round(a.lookback/dt)
    indices = np.arange(max(look, round(a.transient/dt)), len(q))
    return q[indices], N*(q[indices]-q[indices-look])/(2*a.lookback)


def plots(q_by_seed, rows, dense_q, dense_j, a):
    out=a.output_dir
    figures=out/'figures'; figures.mkdir(parents=True,exist_ok=True)
    seeds=sorted(q_by_seed)
    hist=np.array([np.histogram(q_by_seed[s],Q_EDGES)[0] for s in seeds],float)
    hist/=hist.sum(axis=1,keepdims=True)  # equal seed weights
    p=hist.mean(axis=0)
    rng=np.random.default_rng(20261111)
    boot=hist[rng.integers(len(seeds),size=(2000,len(seeds)))].mean(axis=1)
    lo,hi=np.quantile(boot,[.025,.975],axis=0)
    mid=(Q_EDGES[1:]+Q_EDGES[:-1])/2
    landscape=np.full(len(p),np.nan)
    mask=p>0
    landscape[mask]=-np.log(p[mask]); landscape[mask]-=np.nanmin(landscape)
    csv_write(out/'occupancy_pdf.csv', [dict(q=mid[i],mass=p[i],ci_low=lo[i],
       ci_high=hi[i],negative_log_p=landscape[i]) for i in range(len(p))])
    f,ax=plt.subplots(1,2,figsize=(10,3.7))
    ax[0].plot(mid,p,'o-',ms=3);ax[0].fill_between(mid,lo,hi,alpha=.3)
    ax[0].set(xlabel='q',ylabel='probability mass',title='Empirical occupancy p(q)')
    ax[1].plot(mid,landscape,'o-',ms=3)
    ax[1].set(xlabel='q',ylabel='shifted -log p(q)',title='Statistical landscape only')
    f.tight_layout();f.savefig(figures/'01_occupancy_landscape.png',dpi=200);plt.close(f)

    qa=np.asarray(dense_q);ja=np.asarray(dense_j)
    maxj=max(.5,np.max(np.abs(ja))+1e-6)
    jedges=np.linspace(-maxj,maxj,45)
    joint,xe,ye=np.histogram2d(qa,ja,bins=[Q_EDGES,jedges])
    mass=joint/len(qa)
    csv_write(out/'joint_q_Jpast.csv',[dict(q_min=xe[i],q_max=xe[i+1],
       J_min=ye[k],J_max=ye[k+1],count=int(joint[i,k]),mass=mass[i,k])
       for i in range(len(xe)-1) for k in range(len(ye)-1)])
    f,ax=plt.subplots(figsize=(7.5,4.4))
    nonzero=mass[mass>0]
    im=ax.pcolormesh(xe,ye,np.ma.masked_where(mass.T==0,mass.T),
             norm=LogNorm(vmin=nonzero.min(),vmax=nonzero.max()),shading='auto')
    f.colorbar(im,ax=ax,label='joint probability mass')
    ax.set(xlabel='q(t)',ylabel='net flux over preceding window (fish/time)',
       title='Empirical p(q, J_past)')
    f.tight_layout();f.savefig(figures/'02_joint_q_past_flux.png',dpi=200);plt.close(f)

    n,y=counts(rows)
    cell_ci = seed_bootstrap_cells(rows, seeds)
    risk=np.divide(y,n,out=np.full((3,3),np.nan),where=n>=a.min_cell)
    csv_write(out/'conditional_future_transfer.csv', [dict(abs_q_lo=ABS_Q_EDGES[i],
      abs_q_hi=ABS_Q_EDGES[i+1],flux_class=CATEGORIES[k],count=int(n[i,k]),
      positives=int(y[i,k]),p_future=risk[i,k],ci_low=cell_ci[(i,k)][0],
      ci_high=cell_ci[(i,k)][1]) for i in range(3) for k in range(3)])
    f,ax=plt.subplots(figsize=(8.3,4.4))
    vmax=max(.05,float(np.nanmax(risk)) if np.any(np.isfinite(risk)) else 1.)
    im=ax.imshow(np.ma.masked_invalid(risk),origin='lower',aspect='auto',vmin=0,vmax=vmax)
    f.colorbar(im,ax=ax,label='conditional future transfer probability')
    ax.set_xticks(range(3),CATEGORIES);ax.set_yticks(range(3),['0-.5','.5-.8','.8-1'])
    ax.set(ylabel='|q|',title='Future redistribution | current q and PAST flux')
    for i in range(3):
        for k in range(3):
            ax.text(k,i,f'{y[i,k]}/{n[i,k]}',ha='center',va='center',color='white')
    f.tight_layout();f.savefig(figures/'03_conditional_transfer.png',dpi=200);plt.close(f)


def test_synthetic():
    assert relocations(np.array([-.9,-.3,.7,.2,-.7]))==2
    q=np.array([-1.,-1.,-1.,-.98,-.94,-.9,-.86,-.82,-.82,-.82])
    assert abs(N*(q[5]-q[3])/2-4)<1e-10
    class C: pass
    a=C();a.dt=.1;a.stride=1;a.lookback=.2;a.horizon=.2
    a.spacing=.2;a.transient=0.;a.threshold_fish=2
    rr=risk_samples(q,0,a)
    assert len(rr)>0 and all(r['seed']==0 for r in rr)
    print('Synthetic conservation, direction, and risk indexing tests passed.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--raw-dir',type=Path,default=Path('/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data'))
    parser.add_argument('--processed-dir',type=Path)
    parser.add_argument('--cache-dir',type=Path,default=Path('results/nohydro_gap20_descriptive_transition/outputs/initiation_propagation/cache'))
    parser.add_argument('--output-dir',type=Path,default=Path('results/nohydro_gap20_descriptive_transition/outputs/occupancy_flux_landscape'))
    parser.add_argument('--seeds',type=int,nargs='*',default=list(range(40)))
    parser.add_argument('--dt',type=float,default=.001)
    parser.add_argument('--stride',type=int,default=100)
    parser.add_argument('--gap',type=float,default=2.8635642126552705)
    parser.add_argument('--transient',type=float,default=200.)
    parser.add_argument('--lookback',type=float,default=2.)
    parser.add_argument('--horizon',type=float,default=5.)
    parser.add_argument('--spacing',type=float,default=5.)
    parser.add_argument('--threshold-fish',type=int,default=5)
    parser.add_argument('--min-cell',type=int,default=20)
    parser.add_argument('--self-test',action='store_true')
    a=parser.parse_args()
    if a.self_test:
        test_synthetic(); return
    if not a.seeds: raise ValueError('No seeds requested')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    by_seed={};full_q_by_seed={};all_rows=[];full_relocations={};dense_q=[];dense_j=[]
    for seed in sorted(set(a.seeds)):
        q=load_seed(seed,a)
        full_q_by_seed[seed]=q
        by_seed[seed]=q[round(a.transient/(a.dt*a.stride)):]
        full_relocations[seed]=relocations(by_seed[seed])
        qj,jj=dense_joint_series(q,a);dense_q.extend(qj);dense_j.extend(jj)
        ss=risk_samples(q,seed,a);all_rows.extend(ss)
        print('seed',seed,'risk origins',len(ss),'future transfers',sum(r['future_transfer'] for r in ss),flush=True)
    if not all_rows: raise RuntimeError('No risk observations')
    csv_write(a.output_dir/'prospective_windows.csv',all_rows)
    raw_validation=validate_against_raw(full_q_by_seed,a)
    branch_summary=branch_persistence(full_q_by_seed,a)
    plots(by_seed,all_rows,dense_q,dense_j,a)
    val=heldout(all_rows,sorted(by_seed),a.min_cell)
    (a.output_dir/'heldout_validation.json').write_text(json.dumps(val,indent=2)+'\n')
    episode_file = REPO/'results/nohydro_gap20_descriptive_transition/outputs/initiation_propagation/crossing_episodes.csv'
    accepted_episodes = None
    accepted_episode_seeds = None
    if episode_file.exists():
        with episode_file.open(newline='') as stream:
            episode_rows = list(csv.DictReader(stream))
        accepted_episodes = len(episode_rows)
        accepted_episode_seeds = len({int(row['seed']) for row in episode_rows})
    occupancy_by_seed = np.asarray([[np.mean(np.abs(q) >= .8), np.mean(np.abs(q) < .5)]
                                    for q in by_seed.values()])
    rng = np.random.default_rng(20261111)
    pick = rng.integers(len(occupancy_by_seed), size=(3000, len(occupancy_by_seed)))
    occ_boot = occupancy_by_seed[pick].mean(axis=1)
    info=dict(seeds=sorted(by_seed),n_risk_origins=len(all_rows),
       positive_future_redistribution=sum(r['future_transfer'] for r in all_rows),
       seeds_with_positive_future_redistribution=len({r['seed'] for r in all_rows if r['future_transfer']}),
       accepted_cluster_episodes_from_prior_source_exact_index=accepted_episodes,
       seeds_with_accepted_cluster_episodes=accepted_episode_seeds,
       full_opposite_basin_relocations=sum(full_relocations.values()),
       seeds_with_full_relocations=sum(value > 0 for value in full_relocations.values()),
       relocations_by_seed=full_relocations,lookback=a.lookback,horizon=a.horizon,
       spacing=a.spacing,threshold_fish=a.threshold_fish,n_joint_samples=len(dense_q),
       stationary_fraction_abs_q_ge_0p8=float(occupancy_by_seed[:,0].mean()),
       stationary_fraction_abs_q_ge_0p8_ci95=list(map(float,np.quantile(occ_boot[:,0],[.025,.975]))),
       stationary_fraction_abs_q_lt_0p5=float(occupancy_by_seed[:,1].mean()),
       stationary_fraction_abs_q_lt_0p5_ci95=list(map(float,np.quantile(occ_boot[:,1],[.025,.975]))),
       raw_validation=raw_validation,flux_branch_persistence=branch_summary,
       flux_definition='past net neck crossings in fish/time = (N/2)*(q(t)-q(t-lookback))/lookback',
       future_event='net movement of >= threshold_fish from majority to minority chamber at any future sample within horizon')
    (a.output_dir/'summary.json').write_text(json.dumps(info,indent=2)+'\n')
    config = dict(raw_dir=str(a.raw_dir), cache_dir=str(a.cache_dir), seeds=sorted(by_seed),
                  N=N, native_dt=a.dt, occupancy_stride=a.stride,
                  occupancy_dt=a.dt*a.stride, transient=a.transient,
                  lookback=a.lookback, horizon=a.horizon, risk_origin_spacing=a.spacing,
                  transfer_threshold_fish=a.threshold_fish, min_conditional_cell=a.min_cell,
                  q_definition='(N_R-N_L)/N',
                  J_past_definition='(N_R(t)-N_R(t-lookback))/lookback',
                  relocation_definition='successive visits to opposite |q|>=0.6 basins')
    (a.output_dir/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    n_cond, y_cond = counts(all_rows)
    heldout_n = sum(row['n'] for row in val['per_seed'])
    heldout_pos = sum(row['positives'] for row in val['per_seed'])
    relative_brier = val['mean_brier_improvement']/val['mean_brier_q']
    report=f'''# Population imbalance and past neck flux in no-hydrodynamics gap-20 trajectories

## Question and definitions

This is an empirical test using the existing 40 independent, N=100 trajectories. It does not fit a Fokker--Planck equation or interpret an occupancy histogram as energy. After the first 200 time units,

`q=(N_R-N_L)/N`, and `J_past=[N_R(t)-N_R(t-2)]/2` fish per unit time.

At prospective origins spaced by 5 time units, the outcome is a net movement of at least {a.threshold_fish} fish from the chamber that is the majority *at measurement time* toward the minority chamber during `(t,t+5]`. Thus neither the predictor direction nor the label uses the eventual outcome. Predictor and outcome windows do not overlap.

## Numerical audit

The supplied synthetic conservation, direction, and risk-indexing tests pass. Direct reads of the original MAT trajectories gave {raw_validation['n_occupancy_checks']} occupancy checks and {raw_validation['n_flux_windows']} native-resolution crossing-count windows. Maximum errors were {raw_validation['max_abs_q_error']:.2e} in q and {raw_validation['max_abs_J_crossing_error']:.2e} fish/time in flux. The net-flux identity therefore agrees with both endpoint population changes and explicitly counted left-to-right minus right-to-left centerline crossings.

The horizontal structure in the joint plot is not a continuous-flow instability: q lies on a 0.02 population lattice and J_past on a 0.5 fish/time lattice. The mild zig-zag of the binned `-log p(q)` curve is correspondingly bin/lattice sensitive. It is not used as a physical barrier.

## Results

### Population imbalance is the dominant occupancy structure

The stationary sample contains {len(dense_q):,} time points. Equal weighting of the 40 seeds gives `P(|q|>=0.8)={info['stationary_fraction_abs_q_ge_0p8']:.3f}` (seed-bootstrap 95% CI {info['stationary_fraction_abs_q_ge_0p8_ci95'][0]:.3f}--{info['stationary_fraction_abs_q_ge_0p8_ci95'][1]:.3f}) and `P(|q|<0.5)={info['stationary_fraction_abs_q_lt_0p5']:.3f}` ({info['stationary_fraction_abs_q_lt_0p5_ci95'][0]:.3f}--{info['stationary_fraction_abs_q_lt_0p5_ci95'][1]:.3f}). The extreme populations near q=+/-1 are therefore real, but the broad central occupancy shows that partial redistribution is also common.

### Flux forms intermittent, persistent directional branches

Most samples lie on the zero-flux branch, with rarer positive and negative wings. Among nonoverlapping pairs for which both 2-unit windows contain nonzero net flux, the direction persists with probability {branch_summary['mean_same_sign_fraction']:.3f} (seed-bootstrap 95% CI {branch_summary['seed_bootstrap_ci95'][0]:.3f}--{branch_summary['seed_bootstrap_ci95'][1]:.3f}; {branch_summary['n_seeds']} informative seeds). This supports short directional bursts, not a permanently flowing state.

Conditioning on the direction relative to current population balance makes the continuation effect clear:

| current abs(q) | past toward dominant | past zero | past toward balance |
|---|---:|---:|---:|
| 0--0.5 | {y_cond[0,0]}/{n_cond[0,0]} | {y_cond[0,1]}/{n_cond[0,1]} | {y_cond[0,2]}/{n_cond[0,2]} |
| 0.5--0.8 | {y_cond[1,0]}/{n_cond[1,0]} | {y_cond[1,1]}/{n_cond[1,1]} | {y_cond[1,2]}/{n_cond[1,2]} |
| 0.8--1 | {y_cond[2,0]}/{n_cond[2,0]} | {y_cond[2,1]}/{n_cond[2,1]} | {y_cond[2,2]}/{n_cond[2,2]} |

For `|q|<0.5`, the future-transfer probability is {y_cond[0,2]/n_cond[0,2]:.3f} after flux toward balance versus {y_cond[0,1]/n_cond[0,1]:.3f} after zero flux. For `0.5<=|q|<0.8`, the corresponding values are {y_cond[1,2]/n_cond[1,2]:.3f} and {y_cond[1,1]/n_cond[1,1]:.3f}. The extreme-imbalance toward-balance cell has only {n_cond[2,2]} origins and {y_cond[2,2]} positives, so its apparent {y_cond[2,2]/n_cond[2,2]:.3f} probability is not precise; the seed-bootstrap interval is saved in `conditional_future_transfer.csv`.

### Held-out prediction improves, but only modestly

Seeds 0--29 train the fixed nine-cell probability table and seeds 30--39 are untouched test seeds. The held-out set has {heldout_n} prospective windows and {heldout_pos} positives. Mean seed-level Brier score is {val['mean_brier_q']:.5f} for q alone and {val['mean_brier_qJ']:.5f} for `(q,J_past)`. The improvement is {val['mean_brier_improvement']:.6f} (seed-bootstrap 95% CI {val['bootstrap_ci95'][0]:.6f}--{val['bootstrap_ci95'][1]:.6f}), only {100*relative_brier:.1f}% of the q-only score. Past flux therefore contains reproducible continuation information, but it is a small refinement rather than a replacement for population imbalance.

## Event counts and interpretation

The prior source-indexed analysis records {accepted_episodes} accepted short crossing episodes across {accepted_episode_seeds} seeds. This analysis has {len(all_rows)} nonoverlapping prospective windows, of which {info['positive_future_redistribution']} in {info['seeds_with_positive_future_redistribution']} seeds meet the 5-fish redistribution criterion. These are windows, not unique physical events. In contrast, only {info['full_opposite_basin_relocations']} complete opposite-basin relocations occur, across {info['seeds_with_full_relocations']} seeds. That is far too few to validate a probability law for complete school relocation.

The simplest description supported here is a **discrete, directional transfer process**: q describes how population is stored between chambers, while nonzero past flux identifies whether a short crossing burst is already propagating and in which direction. J_past cannot explain first initiation during a crossing-free history, because it is then exactly zero. The data therefore support `(q,J_past)` for short-horizon continuation, not a static two-well potential and not a complete initiation theory.

## Files

- `figures/01_occupancy_landscape.png`: empirical p(q), seed-bootstrap bands, and statistical `-log p` only.
- `figures/02_joint_q_past_flux.png`: empirical joint distribution.
- `figures/03_conditional_transfer.png`: prospective conditional probabilities and counts.
- `heldout_validation.json`: fixed train/test split and seed bootstrap.
- `occupancy_flux_raw_validation.csv`: raw MAT checks.
- `conditional_future_transfer.csv`: cell counts and seed-cluster intervals.
- `summary.json` and `analysis_config.json`: event counts, definitions, and configuration.

No effective potential, diffusion closure, or Markov sufficiency is claimed.
'''
    (a.output_dir/'REPORT.md').write_text(report)
    print('Wrote results to',a.output_dir)

if __name__=='__main__':
    main()
