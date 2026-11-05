#!/usr/bin/env python3
"""Branch-free COM-relative heading harmonics and exact power budgets.

Read-only, bounded replay of the archived circular-confinement trajectories.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

HERE=Path(__file__).resolve().parent
CIRCLE=HERE.parent
sys.path.insert(0,str(CIRCLE))
sys.path.insert(0,str(CIRCLE/'com_relative_turning'))
from matv5_memmap import open_state_rdot  # noqa: E402
from process_collective import FILE_PATTERN,seed_from_name  # noqa: E402
from analyze_com_relative_turning import COMReplay  # noqa: E402

DT=.001; IN=.5; N=100; R=7.158910531638177
SOURCES=('alignment','attraction','wall','direct_hydro','bem_hydro','clip_residual')


def write_csv(path:Path,rows:list[dict[str,Any]]):
    if not rows:return
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def wrap(x):return np.arctan2(np.sin(x),np.cos(x))


def geometry(state:np.ndarray,velocity:np.ndarray|None=None):
    x,y,theta=state;yx=x-x.mean();yy=y-y.mean();rr2=yx*yx+yy*yy
    if np.any(rr2==0):raise FloatingPointError('fish exactly at COM: polar angle undefined')
    phi=np.arctan2(yy,yx);alpha=theta-phi
    result=dict(phi=phi,alpha=alpha,rho=np.sqrt(rr2),min_r_over_R=float(np.sqrt(rr2).min()/R))
    if velocity is not None:
        vr=velocity-velocity.mean(axis=1)[:,None]
        result['omega']=(yx*vr[1]-yy*vr[0])/rr2
    return result


def mode_values(state:np.ndarray,velocity:np.ndarray|None=None):
    g=geometry(state,velocity);phi=g['phi'];u=np.exp(1j*g['alpha']);rho=g['rho'];radial=rho-rho.mean()
    A=[];B=[];H=[];Q=[];V=[]
    mc=float(np.mean(np.sin(g['alpha'])));s=1 if mc>=0 else -1;v=wrap(g['alpha']-s*np.pi/2)
    b0=np.mean(u)
    for m in (1,2):
        f=np.exp(-1j*m*phi);a=np.mean(f);b=np.mean(u*f)
        A.append(a);B.append(b);H.append(b-b0*a);Q.append(np.mean(radial*f));V.append(np.mean(v*f))
    return {**g,'A':np.array(A),'B':np.array(B),'B0':b0,'H':np.array(H),'Q':np.array(Q),'V':np.array(V),
            'P':float(abs(np.mean(np.exp(1j*state[2])))),'M_C':mc}


def stream_seed(path:Path,stride:int=50,t0:float=20,t1:float=999.9):
    state,_,_=open_state_rdot(path);steps=np.arange(round(t0/DT),round(t1/DT)+1,stride,dtype=np.int64)
    n=len(steps);out={k:np.empty(n,np.float32) for k in ('P','M_C','min_r_over_R')}
    for k in ('A','B','H','Q','V'):out[k]=np.empty((n,2),np.complex64)
    for j,step in enumerate(steps):
        mv=mode_values(np.asarray(state[:,:,step],float))
        for k in ('P','M_C','min_r_over_R'):out[k][j]=mv[k]
        for k in ('A','B','H','Q','V'):out[k][j]=mv[k]
    out.update(step=steps,time=steps*DT)
    return out


def source_budget(state:np.ndarray,velocity:np.ndarray,components:dict[str,np.ndarray],m:int):
    mv=mode_values(state,velocity);u=np.exp(1j*mv['alpha']);f=np.exp(-1j*m*mv['phi']);A=np.mean(f);B0=mv['B0'];H=mv['H'][m-1]
    omega=mv['omega'];derivatives={}
    for name in SOURCES:derivatives[name]=1j*np.mean(u*(f-A)*np.asarray(components[name]))
    derivatives['geometry']=(-1j*(m+1)*np.mean(u*f*omega)+1j*A*np.mean(u*omega)+1j*m*B0*np.mean(omega*f))
    derivatives['angular_diffusion']=-.5*IN**2*H
    powers={k:float(2*np.real(np.conj(H)*v)) for k,v in derivatives.items()}
    powers['drift_total']=float(sum(powers.values()))
    powers['finite_N_injection']=float(IN**2/N*(1-abs(A)**2))
    powers['expected_total']=powers['drift_total']+powers['finite_N_injection']
    return mv,H,derivatives,powers


def finite_step_identity(state:np.ndarray,next_state:np.ndarray,omega_eff:np.ndarray,clip:np.ndarray,m:int):
    mv=mode_values(state);H=mv['H'][m-1]
    gp=geometry(next_state);phip=gp['phi'];fp=np.exp(-1j*m*phip);Ap=np.mean(fp)
    theta_star=state[2]+omega_eff*DT
    c=np.exp(1j*(theta_star-phip))*(fp-Ap);decay=np.exp(-.5*IN**2*DT)
    expected=decay*np.mean(c);variance=(1-np.exp(-IN**2*DT))*(1-abs(Ap)**2)/N
    up=np.exp(1j*(next_state[2]-phip));actual=np.mean(up*fp)-np.mean(up)*Ap
    # Exact aggregate effect of clipping on conditional mean-square amplitude.
    theta_raw=state[2]+(omega_eff-clip)*DT
    c_raw=np.exp(1j*(theta_raw-phip))*(fp-Ap)
    exact_clip_power_delta=abs(decay*np.mean(c))**2-abs(decay*np.mean(c_raw))**2
    return dict(H=H,c=c,A_next=Ap,expected=expected,variance=float(variance),actual=actual,
        exact_expected_power=float(abs(expected)**2+variance),actual_power=float(abs(actual)**2),
        exact_clip_power_delta=float(exact_clip_power_delta))


def replay_snapshot(replay:COMReplay,seed:int,step:int):
    base,detail=replay.evaluate(seed,step);state_mm,_=replay.arrays(seed)
    state=np.asarray(state_mm[:,:,step],float);next_state=np.asarray(state_mm[:,:,step+1],float)
    out={'seed':seed,'step':step,'time':step*DT,'P':base['P'],'Pi':base['Pi'],'M_C':base['M_C'],
         'n_angular_clip':base['n_angular_clip'],'n_collision_override':base['n_collision_override'],
         'n_sticky_override':base['n_sticky_override'],'velocity_replay_max_error':base['velocity_replay_max_error'],
         'social_split_max_error':base['social_split_max_error']}
    clip=np.asarray(detail['components']['clip_residual']);out['clip_max_abs_rate']=float(np.max(abs(clip)))
    out['clip_sum_abs_angle']=float(np.sum(abs(clip))*DT)
    for m in (1,2):
        mv,H,derivs,powers=source_budget(state,detail['effective_velocity'],detail['components'],m)
        fs=finite_step_identity(state,next_state,detail['omega_effective'],clip,m)
        current_power=abs(H)**2;ct_delta=powers['expected_total']*DT
        exact_delta=fs['exact_expected_power']-current_power;actual_delta=fs['actual_power']-current_power
        prefix=f'H{m}_';out[prefix+'real']=H.real;out[prefix+'imag']=H.imag;out[prefix+'abs']=abs(H)
        out[f'A{m}_abs']=abs(mv['A'][m-1]);out[f'Q{m}_abs']=abs(mv['Q'][m-1]);out[f'V{m}_abs']=abs(mv['V'][m-1])
        for name in (*SOURCES,'geometry','angular_diffusion','drift_total','finite_N_injection','expected_total'):
            out[prefix+'power_'+name]=powers[name]
        out[prefix+'ct_pred_delta_power']=ct_delta
        out[prefix+'exact_conditional_delta_power']=exact_delta
        out[prefix+'discretization_residual']=exact_delta-ct_delta
        out[prefix+'actual_delta_power']=actual_delta
        out[prefix+'realized_noise_residual']=actual_delta-exact_delta
        out[prefix+'conditional_variance']=fs['variance'];out[prefix+'exact_clip_power_delta']=fs['exact_clip_power_delta']
        out[prefix+'actual_next_real']=fs['actual'].real;out[prefix+'actual_next_imag']=fs['actual'].imag
        out[prefix+'expected_next_real']=fs['expected'].real;out[prefix+'expected_next_imag']=fs['expected'].imag
    out['min_r_com_over_R']=mv['min_r_over_R']
    return out


def replay_intervals(replay:COMReplay,outdir:Path):
    specs={7:[(100,160),(470,530),(800,860)],8:[(100,160),(470,530),(800,860)],
           9:[(250,360)],0:[(100,160)]};rows=[]
    for seed,intervals in specs.items():
        for interval,(t0,t1) in enumerate(intervals):
            for step in np.arange(round(t0/DT),round(t1/DT)+1,200,dtype=np.int64):
                r=replay_snapshot(replay,seed,int(step));r.update(interval=interval,t0=t0,t1=t1);rows.append(r)
        print(f'replayed seed {seed}: {sum(r["seed"]==seed for r in rows)} snapshots',flush=True)
    write_csv(outdir/'branch_free_snapshot_budget.csv',rows);return rows


def redraw_validation(replay:COMReplay,outdir:Path,n_draw:int=20000):
    selected=[(7,120000),(8,120000),(9,280000),(9,309000),(9,312000),(9,320000)]
    rng=np.random.default_rng(20261105);rows=[]
    for seed,step in selected:
        _,detail=replay.evaluate(seed,step);state_mm,_=replay.arrays(seed)
        state=np.asarray(state_mm[:,:,step],float);next_state=np.asarray(state_mm[:,:,step+1],float)
        for m in (1,2):
            fs=finite_step_identity(state,next_state,detail['omega_effective'],detail['components']['clip_residual'],m)
            z=rng.standard_normal((n_draw,N));draw=np.mean(fs['c'][None,:]*np.exp(1j*IN*np.sqrt(DT)*z),axis=1)
            sample_mean=np.mean(draw);sample_var=float(np.mean(abs(draw-sample_mean)**2));se=np.sqrt(fs['variance']/n_draw)
            rows.append(dict(seed=seed,step=step,time=step*DT,m=m,n_draw=n_draw,
                analytic_mean_real=fs['expected'].real,analytic_mean_imag=fs['expected'].imag,
                sample_mean_real=sample_mean.real,sample_mean_imag=sample_mean.imag,
                mean_abs_error=abs(sample_mean-fs['expected']),mean_error_over_se=abs(sample_mean-fs['expected'])/max(se,1e-30),
                analytic_variance=fs['variance'],sample_variance=sample_var,
                variance_relative_error=(sample_var-fs['variance'])/fs['variance']))
    write_csv(outdir/'conditional_redraw_validation.csv',rows);return rows


def interval_summary(rows:list[dict]):
    out=[]
    for seed,interval in sorted(set((r['seed'],r['interval']) for r in rows)):
        use=[r for r in rows if r['seed']==seed and r['interval']==interval]
        rec=dict(seed=seed,interval=interval,t0=use[0]['t0'],t1=use[0]['t1'],n=len(use),
            mean_P=float(np.mean([r['P'] for r in use])),mean_M_C=float(np.mean([r['M_C'] for r in use])),
            mean_abs_H1=float(np.mean([r['H1_abs'] for r in use])),mean_abs_H2=float(np.mean([r['H2_abs'] for r in use])),
            mean_abs_Q2=float(np.mean([r['Q2_abs'] for r in use])),mean_abs_V2=float(np.mean([r['V2_abs'] for r in use])),
            clip_fish_count=int(sum(r['n_angular_clip'] for r in use)),clip_snapshot_count=int(sum(r['n_angular_clip']>0 for r in use)),
            max_clip_abs_rate=max(r['clip_max_abs_rate'] for r in use),sum_clip_abs_angle=float(sum(r['clip_sum_abs_angle'] for r in use)))
        for m in (1,2):
            for name in (*SOURCES,'geometry','angular_diffusion','drift_total','finite_N_injection','expected_total'):
                rec[f'mean_H{m}_power_{name}']=float(np.mean([r[f'H{m}_power_{name}'] for r in use]))
            for name in ('ct_pred_delta_power','exact_conditional_delta_power','discretization_residual','actual_delta_power','realized_noise_residual','exact_clip_power_delta'):
                rec[f'mean_H{m}_{name}']=float(np.mean([r[f'H{m}_{name}'] for r in use]))
        out.append(rec)
    return out


def seed9_windows(rows:list[dict],series:dict):
    width=6.4;target_t0=306.;ti=int(np.argmin(abs(series['time']-target_t0)))
    stable=(series['time']>=250)&(series['time']<300);sp=max(float(np.std(series['P'][stable])),1e-8);sm=max(float(np.std(series['M_C'][stable])),1e-8)
    candidates=np.arange(250.,300.-width+1e-9,width);dist=[]
    for t0 in candidates:
        j=int(np.argmin(abs(series['time']-t0)))
        dist.append(((float(series['P'][j]-series['P'][ti])/sp)**2+(float(series['M_C'][j]-series['M_C'][ti])/sm)**2)**.5)
    control_t0=float(candidates[int(np.argmin(dist))]);match_distance=float(min(dist))
    windows=[('school',250,300),('matched_school_control',control_t0,control_t0+width),('late_school',300,306),
             ('precollapse',target_t0,312.4),('collapse',312.4,315),('early_mill',315,330),('mill',330,360)]
    out=[]
    for label,t0,t1 in windows:
        rs=[r for r in rows if r['seed']==9 and t0<=r['time']<t1];u=(series['time']>=t0)&(series['time']<t1)
        j0=int(np.argmin(abs(series['time']-t0)))
        rec=dict(window=label,t0=t0,t1=t1,n_snapshot=len(rs),initial_P=float(series['P'][j0]),initial_M_C=float(series['M_C'][j0]),
            state_match_distance=match_distance if label=='matched_school_control' else 0.0 if label=='precollapse' else np.nan,
            mean_P=float(np.mean(series['P'][u])),
            mean_M_C=float(np.mean(series['M_C'][u])),mean_abs_H1=float(np.mean(abs(series['H'][u,0]))),
            mean_abs_H2=float(np.mean(abs(series['H'][u,1]))),mean_abs_Q2=float(np.mean(abs(series['Q'][u,1]))),
            mean_abs_V2=float(np.mean(abs(series['V'][u,1]))),clip_fish_count=int(sum(r['n_angular_clip'] for r in rs)))
        for name in (*SOURCES,'geometry','angular_diffusion','drift_total','finite_N_injection','expected_total'):
            rec['mean_H2_power_'+name]=float(np.mean([r['H2_power_'+name] for r in rs]))
        out.append(rec)
    return out


def source_state_summary(rows:list[dict]):
    groups={'school_controls':[r for r in rows if (r['seed']==0) or (r['seed']==9 and 250<=r['time']<300)],
            'stable_mills':[r for r in rows if r['seed'] in (7,8)]}
    out=[]
    for state,blocks in groups.items():
        for name in (*SOURCES,'geometry','angular_diffusion','drift_total','finite_N_injection','expected_total'):
            # Independent seeds, not the three within-seed blocks, are uncertainty units.
            seedmeans=[]
            for seed in sorted(set(r['seed'] for r in blocks)):
                seedmeans.append(np.mean([r['H2_power_'+name] for r in blocks if r['seed']==seed]))
            mean=float(np.mean(seedmeans));lo=float(min(seedmeans));hi=float(max(seedmeans))
            out.append(dict(state=state,source=name,n_independent_seeds=len(seedmeans),mean=mean,seed_mean_min=lo,seed_mean_max=hi,
                same_sign_all_seeds=int(all(x>0 for x in seedmeans) or all(x<0 for x in seedmeans))))
    return out


def plots(series:dict[int,dict],rows:list[dict],intervals:list[dict],redraw:list[dict],outdir:Path):
    d=series[9];u=(d['time']>=250)&(d['time']<=360)
    fig,axs=plt.subplots(3,1,figsize=(12,9),sharex=True)
    axs[0].plot(d['time'][u],d['P'][u],label='P');axs[0].plot(d['time'][u],d['M_C'][u],label='M_C');axs[0].legend(ncol=2)
    axs[1].plot(d['time'][u],abs(d['H'][u,1]),label='|H2| branch-free');axs[1].plot(d['time'][u],abs(d['H'][u,0]),label='|H1| control');axs[1].legend(ncol=2)
    axs[2].plot(d['time'][u],abs(d['Q'][u,1]),label='|Q2|');axs[2].plot(d['time'][u],abs(d['V'][u,1]),label='|V2| wrapped',alpha=.75);axs[2].legend(ncol=2)
    for ax in axs:ax.axvline(312.41,color='black',ls='--',lw=.8);ax.grid(alpha=.2)
    axs[-1].set_xlabel('time');fig.tight_layout();fig.savefig(outdir/'seed9_branch_free_timeline.png',dpi=220);plt.close(fig)

    stable=[r for r in source_state_summary(rows) if r['state']=='stable_mills'];names=[r['source'] for r in stable]
    means=np.array([r['mean'] for r in stable]);lo=np.array([r['seed_mean_min'] for r in stable]);hi=np.array([r['seed_mean_max'] for r in stable])
    fig,ax=plt.subplots(figsize=(11,5));x=np.arange(len(names));ax.bar(x,means,color='#4c78a8')
    ax.errorbar(x,means,yerr=np.vstack((means-lo,hi-means)),fmt='none',ecolor='black',capsize=4)
    ax.axhline(0,color='black',lw=.8);ax.set_xticks(x,[n.replace('_','\n') for n in names]);ax.set_ylabel('H2 power production')
    ax.set_title('Stable mills: mean and range of independent seed means (n=2)');ax.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(outdir/'stable_mill_H2_source_budget.png',dpi=220);plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax,m in zip(axs,(1,2)):
        rr=[r for r in redraw if r['m']==m];ax.scatter([r['analytic_variance'] for r in rr],[r['sample_variance'] for r in rr])
        vals=[r['analytic_variance'] for r in rr]+[r['sample_variance'] for r in rr];a=min(vals);b=max(vals);ax.plot([a,b],[a,b],'k--',lw=.8)
        ax.set(xlabel='analytic variance',ylabel='redraw variance',title=f'H{m}: 20,000 redraws');ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(outdir/'conditional_identity_validation.png',dpi=220);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--source',type=Path,default=CIRCLE.parents[1].parent/'BEM_circular_parallel_instance.py')
    p.add_argument('--outdir',type=Path,default=HERE/'outputs');args=p.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    files=sorted(args.data_dir.glob(FILE_PATTERN),key=seed_from_name);paths={seed_from_name(f):f for f in files}
    if len(paths)!=10:raise FileNotFoundError(f'Expected 10 raw trajectories, found {len(paths)} in {args.data_dir}')
    series={}
    for seed in (0,7,8,9):
        series[seed]=stream_seed(paths[seed]);np.savez_compressed(args.outdir/f'branch_free_series_seed{seed}.npz',**series[seed]);print(f'streamed seed {seed}',flush=True)
    source=args.source.resolve();expected=json.loads((CIRCLE/'simulation_parameters.json').read_text())['source_sha256'];sha=hashlib.sha256(source.read_bytes()).hexdigest()
    if sha!=expected:raise RuntimeError(f'production source hash mismatch: {sha}')
    replay=COMReplay(args.data_dir,source,CIRCLE/'outputs');rows=replay_intervals(replay,args.outdir);redraw=redraw_validation(replay,args.outdir)
    intervals=interval_summary(rows);write_csv(args.outdir/'interval_source_budget.csv',intervals)
    write_csv(args.outdir/'state_source_budget.csv',source_state_summary(rows))
    write_csv(args.outdir/'seed9_transition_windows.csv',seed9_windows(rows,series[9]));plots(series,rows,intervals,redraw,args.outdir)
    disc=np.array([r['H2_discretization_residual'] for r in rows]);exact=np.array([r['H2_exact_conditional_delta_power'] for r in rows]);noise=np.array([r['H2_realized_noise_residual'] for r in rows])
    closure=np.array([r['H2_actual_delta_power']-r['H2_ct_pred_delta_power']-r['H2_discretization_residual']-r['H2_realized_noise_residual'] for r in rows])
    phi_test=np.concatenate((np.linspace(-2.7,-1.8,70),np.linspace(.2,2.4,30)));u_test=np.full(N,1j)
    tangent_residual=max(abs(np.mean(u_test*np.exp(-1j*m*phi_test))-np.mean(u_test)*np.mean(np.exp(-1j*m*phi_test))) for m in (1,2))
    audit=dict(source=str(source),source_sha256=sha,update_order='noise RNG called before final update, but position noise is identically zero and noise never enters rDot/collision/retry; x[n+1] is conditionally deterministic',
        n_replay_snapshots=len(rows),n_redraw_snapshots=6,redraws_per_snapshot=20000,
        max_velocity_replay_error=max(r['velocity_replay_max_error'] for r in rows),max_social_split_error=max(r['social_split_max_error'] for r in rows),
        max_H2_conditional_discretization_residual=float(np.max(abs(disc))),rms_H2_conditional_discretization_residual=float(np.sqrt(np.mean(disc**2))),
        median_abs_H2_conditional_discretization_residual=float(np.median(abs(disc))),p95_abs_H2_conditional_discretization_residual=float(np.quantile(abs(disc),.95)),
        mean_abs_discretization_over_mean_abs_exact_delta=float(np.mean(abs(disc))/np.mean(abs(exact))),
        max_H2_algebraic_closure_error=float(np.max(abs(closure))),mean_H2_realized_noise_residual=float(np.mean(noise)),
        se_H2_realized_noise_residual=float(stats.sem(noise)),
        max_redraw_mean_error_over_se=max(r['mean_error_over_se'] for r in redraw),max_redraw_variance_relative_error=max(abs(r['variance_relative_error']) for r in redraw),
        max_perfect_tangent_H1_H2_residual=float(tangent_residual),
        total_clipped_fish=sum(r['n_angular_clip'] for r in rows),clip_snapshots=sum(r['n_angular_clip']>0 for r in rows),
        minimum_r_com_over_R=min(r['min_r_com_over_R'] for r in rows))
    (args.outdir/'validation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (args.outdir/'analysis_config.json').write_text(json.dumps(dict(stage_A_dt=.05,stage_B_dt=.2,dt=DT,In=IN,N=N,
        intervals={'7':[[100,160],[470,530],[800,860]],'8':[[100,160],[470,530],[800,860]],'9':[[250,360]],'0':[[100,160]]}),indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)


if __name__=='__main__':main()
