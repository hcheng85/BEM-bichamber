#!/usr/bin/env python3
"""Exact frozen-neighbor social-alignment landscape and closure test."""
from __future__ import annotations
import argparse, csv, gzip, json, sys, time
from pathlib import Path
import numpy as np
import scipy.io as sio

HERE=Path(__file__).resolve().parent
APPROACH=HERE.parents[1]
sys.path.insert(0,str(APPROACH/"pre_crossing_turning"/"scripts"))
from replay_dynamics import Parameters, analytic_geometry, visible_neighbors  # noqa:E402
sys.path.insert(0,str(APPROACH/"collective_coordinates"/"scripts"))
from analyze_collective_coordinates import auc_score,fit_logistic,predict_logistic  # noqa:E402

DT=0.001
TIMES=(0.0,0.1,0.2,0.5,1.0)
TRAIN=tuple(range(15)); TEST=tuple(range(15,20)); SEEDS=tuple(range(20))
GRID=np.linspace(-np.pi,np.pi,721)

def read_csv(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:return list(csv.DictReader(f))

def write_csv(path,rows):
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True);fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

def wrap(x):return np.arctan2(np.sin(x),np.cos(x))

def angle_to_interval(theta,lo,hi):
    inside=(theta>=lo)&(theta<=hi)
    dlo=wrap(theta-lo);dhi=wrap(theta-hi)
    return np.where(inside,0,np.where(np.abs(dlo)<=np.abs(dhi),dlo,dhi))

def angular_geometry(state,origin,p):
    side=1.0 if origin=='left' else -1.0
    d=-side*state[0];heading=np.arctan2(np.sin(state[2]),side*np.cos(state[2]))
    clo=-p.gap/2-p.tip_r;chi=p.gap/2+p.tip_r
    alo=np.arctan2(clo-state[1],d);ahi=np.arctan2(chi-state[1],d)
    blo=np.arcsin(np.clip(p.tip_r/np.hypot(d,clo-state[1]),0,1))
    bhi=np.arcsin(np.clip(p.tip_r/np.hypot(d,chi-state[1]),0,1))
    lo=alo+blo;hi=ahi-bhi;eta=angle_to_interval(heading,lo,hi)
    return side,d,heading,lo,hi,eta

def coefficients(state,i,neighbors):
    js=np.asarray(neighbors,int);z=state[0]+1j*state[1]
    phi=np.angle(z[js]-z[i]);th=state[2,js]
    return dict(n=len(js),A1=np.sum(np.exp(1j*th)),A2=np.sum(np.exp(1j*(th+phi))),
                B=np.sum(np.exp(1j*phi)),C=.5*np.sum(np.sin(th-phi)),neighbors=tuple(js.tolist()))

def harmonic_global(theta,c,p):
    theta=np.asarray(theta);den=c['n']+np.real(c['B']*np.exp(-1j*theta))
    num=(np.imag(c['A1']*np.exp(-1j*theta))+.5*np.imag(c['A2']*np.exp(-2j*theta))+c['C'])
    return p.Ip*num/den

def direct_global(theta,state,i,neighbors,p):
    theta=np.asarray(theta);js=np.asarray(neighbors,int);z=state[0]+1j*state[1]
    phi=np.angle(z[js]-z[i]);th=state[2,js]
    rel=phi[:,None]-theta[None,:];visual=1+np.cos(rel)
    den=np.sum(visual,axis=0)
    return p.Ip*np.sum(np.sin(th[:,None]-theta[None,:])*visual,axis=0)/den

def local_torque(local_theta,c,side,p):
    global_theta=local_theta if side>0 else np.pi-local_theta
    return side*harmonic_global(global_theta,c,p)

def roots(theta,omega,stable):
    out=[]
    for k in range(len(theta)-1):
        a,b=omega[k],omega[k+1]
        hit=(a>0 and b<=0) if stable else (a<0 and b>=0)
        if hit:
            frac=abs(a)/(abs(a)+abs(b)) if abs(a)+abs(b)>0 else .5
            out.append(float(theta[k]+frac*(theta[k+1]-theta[k])))
    return out

def circ_dist(a,b):return abs(float(wrap(a-b)))

def interval_distance(a,lo,hi):return abs(float(angle_to_interval(np.asarray(a),lo,hi)))

def arc_contains(x,a,b):
    """Whether x lies on the positive/counterclockwise arc a -> b."""
    return float((x-a)%(2*np.pi)) <= float((b-a)%(2*np.pi))+1e-12

def integrate_arc(c,side,p,a,b,direction,n=600):
    length=(b-a)%(2*np.pi) if direction>0 else (a-b)%(2*np.pi)
    ss=np.linspace(0,length,n);th=a+direction*ss
    om=local_torque(wrap(th),c,side,p)
    return float(np.trapezoid(-om*direction,ss))

def path_barrier(c,side,p,start,target,direction,n=600):
    """Maximum uphill V along one directed path from start to target."""
    length=(target-start)%(2*np.pi) if direction>0 else (start-target)%(2*np.pi)
    ss=np.linspace(0,length,n);th=start+direction*ss
    om=local_torque(wrap(th),c,side,p)
    inc=-.5*(om[:-1]+om[1:])*direction*np.diff(ss)
    cumulative=np.r_[0.0,np.cumsum(inc)]
    return float(max(0.0,np.max(cumulative)))

def landscape_features(state,i,origin,c,p):
    side,d,heading,lo,hi,eta=angular_geometry(state[:,i],origin,p)
    omega=local_torque(GRID,c,side,p)
    direct=direct_global(GRID if side>0 else np.pi-GRID,state,i,c['neighbors'],p)*side
    verify=float(np.nanmax(np.abs(omega-direct)))
    V=np.zeros_like(GRID);V[1:]=np.cumsum(-.5*(omega[:-1]+omega[1:])*np.diff(GRID))
    stable=roots(GRID,omega,True);unstable=roots(GRID,omega,False)
    stable_dist=min((interval_distance(a,lo,hi) for a in stable),default=np.nan)
    stable_gate=any(lo<=a<=hi for a in stable)
    actual=float(local_torque(np.asarray([heading]),c,side,p)[0])
    closing=float(-np.sign(eta)*actual) if eta!=0 else 0.0
    barrier_exists=0;barrier=np.nan
    if stable:
        if eta==0:barrier=0.0
        else:
            barrier=min(path_barrier(c,side,p,heading,target,direction)
                        for target in (lo,hi) for direction in (-1,+1))
            barrier_exists=int(barrier>1e-6)
    return dict(initial_distance=float(d),initial_heading=float(heading),gate_low=float(lo),gate_high=float(hi),
        initial_eta=float(eta),initial_abs_eta=abs(float(eta)),alignment_torque=float(actual),
        alignment_closing_rate=closing,n_neighbors=c['n'],harmonic_max_abs_error=verify,
        n_stable=len(stable),n_unstable=len(unstable),stable_in_gate=int(stable_gate),
        stable_distance_to_gate=float(stable_dist) if np.isfinite(stable_dist) else np.nan,
        barrier_exists=barrier_exists,barrier_height=float(barrier) if np.isfinite(barrier) else np.nan,
        basin_intersects_gate=int(stable_gate),
        cycle_tilt=float(V[-1]-V[0]),torque_rms=float(np.sqrt(np.mean(omega**2))),
        stable_angles=';'.join(f'{a:.8g}' for a in stable),unstable_angles=';'.join(f'{a:.8g}' for a in unstable))

def curve_similarity(c0,c1,side,p):
    a=local_torque(GRID,c0,side,p);b=local_torque(GRID,c1,side,p)
    if np.std(a)<1e-12 or np.std(b)<1e-12:corr=np.nan
    else:corr=float(np.corrcoef(a,b)[0,1])
    rel=float(np.sqrt(np.mean((a-b)**2))/(np.sqrt(np.mean(a*a))+1e-12))
    return corr,rel

def extract(risk,fish_table,raw_dir,p,geom):
    fishkey={(int(r['seed']),int(r['pair_id']),r['condition'],int(r['leading_rank'])):r for r in fish_table}
    landscape=[];stability=[];examples=[]
    for seed in SEEDS:
        print(f'seed {seed}: load',flush=True);t=time.time()
        raw=raw_dir/f'If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat'
        state=np.asarray(sio.loadmat(raw,variable_names=['state'],squeeze_me=True)['state'])
        rows=[r for r in risk if int(r['seed'])==seed]
        for rr in rows:
            i0=int(round(float(rr['risk_time'])/DT));origin=rr['origin'];side=1.0 if origin=='left' else -1.0
            st0=state[:,:,i0];vn0=visible_neighbors(st0[:2],geom,p)
            fixed=np.argsort(-side*st0[0])[:10]
            coeff0={int(i):coefficients(st0,int(i),np.flatnonzero(vn0[int(i)])) for i in fixed}
            for rank,i in enumerate(fixed,1):
                base=fishkey[(seed,int(rr['pair_id']),rr['condition'],rank)]
                feat=landscape_features(st0,int(i),origin,coeff0[int(i)],p)
                landscape.append(dict(seed=seed,pair_id=rr['pair_id'],condition=rr['condition'],origin=origin,
                    fish=int(i),leading_rank=rank,outcome_gate_entry=int(base['gate_entered']),
                    initially_admissible=int(base['initially_admissible']),entry_time=base['entry_time_from_A_minus_2'],**feat))
            for tau in TIMES:
                idx=i0+int(round(tau/DT));st=state[:,:,idx];vn=visible_neighbors(st[:2],geom,p)
                for rank,i in enumerate(fixed,1):
                    c0=coeff0[int(i)];c1=coefficients(st,int(i),np.flatnonzero(vn[int(i)]))
                    j0=set(c0['neighbors']);j1=set(c1['neighbors']);union=j0|j1
                    jac=len(j0&j1)/len(union) if union else 1.0
                    corr,rel=curve_similarity(c0,c1,side,p)
                    _,_,heading,_,_,_=angular_geometry(st[:,i],origin,p)
                    frozen=float(local_torque(np.asarray([heading]),c0,side,p)[0])
                    current=float(local_torque(np.asarray([heading]),c1,side,p)[0])
                    stability.append(dict(seed=seed,pair_id=rr['pair_id'],condition=rr['condition'],fish=int(i),leading_rank=rank,
                        time=tau,neighbor_jaccard=jac,landscape_correlation=corr,relative_landscape_rms_change=rel,
                        frozen_torque_at_current_heading=frozen,current_torque_at_current_heading=current,
                        torque_abs_error=abs(frozen-current),torque_sign_agreement=int(np.sign(frozen)==np.sign(current))))
        print(f' seed {seed} done {time.time()-t:.1f}s',flush=True);del state
    return landscape,stability

MODEL_SPECS={
 'baseline_eta_distance':('initial_abs_eta','initial_distance'),
 'baseline_plus_drift':('initial_abs_eta','initial_distance','alignment_closing_rate'),
 'baseline_plus_stability':('initial_abs_eta','initial_distance','alignment_closing_rate','stable_distance_filled','stable_in_gate','has_stable'),
 'full_landscape':('initial_abs_eta','initial_distance','alignment_closing_rate','stable_distance_filled','stable_in_gate','has_stable','barrier_height_filled','barrier_exists','abs_cycle_tilt'),
 'landscape_without_heading':('initial_distance','alignment_closing_rate','stable_distance_filled','stable_in_gate','has_stable','barrier_height_filled','barrier_exists','abs_cycle_tilt')}

def predictive(rows):
    q=[]
    for r in rows:
        if int(r['initially_admissible']):continue
        x=dict(r);x['has_stable']=int(int(r['n_stable'])>0)
        x['stable_distance_filled']=float(r['stable_distance_to_gate']) if np.isfinite(float(r['stable_distance_to_gate'])) else np.pi
        x['barrier_height_filled']=float(r['barrier_height']) if np.isfinite(float(r['barrier_height'])) else 0.0
        x['abs_cycle_tilt']=abs(float(r['cycle_tilt']));q.append(x)
    train=[r for r in q if int(r['seed']) in TRAIN];test=[r for r in q if int(r['seed']) in TEST]
    pooled=[];byseed=[]
    for name,keys in MODEL_SPECS.items():
        xtr=np.array([[float(r[k]) for k in keys] for r in train]);ytr=np.array([int(r['outcome_gate_entry']) for r in train])
        model=fit_logistic(xtr,ytr)
        for label,qq in [('pooled',test)]+[(f'seed_{s}',[r for r in test if int(r['seed'])==s]) for s in TEST]:
            x=np.array([[float(r[k]) for k in keys] for r in qq]);y=np.array([int(r['outcome_gate_entry']) for r in qq]);prob=predict_logistic(model,x)
            row=dict(model=name,sample=label,n=len(qq),positives=int(y.sum()),auc=auc_score(y,prob),
                brier=float(np.mean((prob-y)**2)),nll=float(-np.mean(y*np.log(prob+1e-12)+(1-y)*np.log(1-prob+1e-12))))
            (pooled if label=='pooled' else byseed).append(row)
    contrasts=[];rng=np.random.default_rng(20261023)
    for name in MODEL_SPECS:
        if name=='baseline_eta_distance':continue
        a={r['sample']:r for r in byseed if r['model']==name};b={r['sample']:r for r in byseed if r['model']=='baseline_eta_distance'}
        for metric in ('auc','brier','nll'):
            d=np.array([float(a[f'seed_{s}'][metric])-float(b[f'seed_{s}'][metric]) for s in TEST])
            boot=np.mean(d[rng.integers(0,len(d),(30000,len(d)))],axis=1)
            contrasts.append(dict(model=name,metric=f'delta_{metric}_versus_baseline',mean_seed_difference=float(d.mean()),
                seed_bootstrap_ci_low=float(np.percentile(boot,2.5)),seed_bootstrap_ci_high=float(np.percentile(boot,97.5)),heldout_seeds=5,
                seed_differences=';'.join(f'{v:.8g}' for v in d)))
    return q,pooled,byseed,contrasts

def stability_summary(rows):
    out=[];rng=np.random.default_rng(20261023)
    for tau in TIMES:
        q=[r for r in rows if float(r['time'])==tau]
        for metric in ('neighbor_jaccard','landscape_correlation','relative_landscape_rms_change','torque_abs_error','torque_sign_agreement'):
            vals=[]
            for seed in SEEDS:
                v=np.array([float(r[metric]) for r in q if int(r['seed'])==seed]);vals.append(float(np.nanmean(v)))
            vals=np.asarray(vals);boot=np.nanmean(vals[rng.integers(0,20,(30000,20))],axis=1)
            out.append(dict(time=tau,metric=metric,mean=float(np.nanmean(vals)),seed_ci_low=float(np.nanpercentile(boot,2.5)),
                seed_ci_high=float(np.nanpercentile(boot,97.5)),independent_seeds=20))
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--raw-dir',type=Path,required=True)
    ap.add_argument('--risk-set',type=Path,required=True);ap.add_argument('--fish-table',type=Path,required=True);args=ap.parse_args()
    out=args.out.resolve();(out/'tables').mkdir(parents=True,exist_ok=True)
    p=Parameters();geom=analytic_geometry(p);risk=read_csv(args.risk_set);fish=read_csv(args.fish_table)
    landscape,stability=extract(risk,fish,args.raw_dir,p,geom)
    predrows,pooled,byseed,contrasts=predictive(landscape);stable_summary=stability_summary(stability)
    write_csv(out/'tables'/'landscape_features.csv.gz',landscape);write_csv(out/'tables'/'frozen_environment.csv.gz',stability)
    write_csv(out/'tables'/'prediction_pooled.csv',pooled);write_csv(out/'tables'/'prediction_by_seed.csv',byseed)
    write_csv(out/'tables'/'prediction_contrasts.csv',contrasts);write_csv(out/'tables'/'frozen_environment_summary.csv',stable_summary)
    config=dict(parent_branch='analysis/leading-fish-angular-first-passage',parent_commit='e0e067298096955b5174a6a341fac9f0cce792c9',
        seeds=list(SEEDS),training_seeds=list(TRAIN),heldout_seeds=list(TEST),n_pairs=456,n_fish=len(landscape),
        heading_grid_points=len(GRID),stability_times=list(TIMES),potential='exact integral of frozen social-alignment torque only',
        empirical_potential_fitted=False,new_simulation=False)
    (out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')

if __name__=='__main__':main()
