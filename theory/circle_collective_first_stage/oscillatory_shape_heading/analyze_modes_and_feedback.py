#!/usr/bin/env python3
"""Targeted test of a co-rotating shape--heading feedback in circular schools."""
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
from scipy import signal, stats

HERE=Path(__file__).resolve().parent
CIRCLE=HERE.parent
sys.path.insert(0,str(CIRCLE))
sys.path.insert(0,str(CIRCLE/'com_relative_turning'))

from matv5_memmap import open_state_rdot  # noqa: E402
from process_collective import FILE_PATTERN,seed_from_name  # noqa: E402
from analyze_com_relative_turning import COMReplay  # noqa: E402

R=7.158910531638177; DT=.001; N=100; CUTOFF=.005*R
REFERENCE_OMEGA=.55169; REFERENCE_PERIOD=11.4


def write_csv(path:Path,rows:list[dict[str,Any]]):
    if not rows:return
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def wrap(x):return np.arctan2(np.sin(x),np.cos(x))


def frame_modes(state:np.ndarray,velocity:np.ndarray,cutoff:float=CUTOFF):
    x,y,th=state; cx,cy=np.mean(x),np.mean(y); yx=x-cx; yy=y-cy
    vx,vy=velocity; vx=vx-np.mean(vx); vy=vy-np.mean(vy)
    rr=np.hypot(yx,yy); valid=rr>cutoff; phi=np.arctan2(yy,yx); alpha=wrap(th-phi)
    mc=float(np.mean(np.sin(alpha[valid]))); s=-1 if mc<0 else 1
    rho=rr[valid]; ph=phi[valid]; al=alpha[valid]; radial=rho-np.mean(rho)
    omega=(yx[valid]*vy[valid]-yy[valid]*vx[valid])/(rr[valid]**2)
    vhead=wrap(al-s*np.pi/2)
    Q=np.array([np.mean(radial*np.exp(-1j*m*ph)) for m in range(1,5)])
    A=np.array([np.mean(np.exp(-1j*m*ph)) for m in range(1,5)])
    V=np.array([np.mean(vhead*np.exp(-1j*m*ph)) for m in range(1,5)])
    ur=(yx[valid]*vx[valid]+yy[valid]*vy[valid])/rr[valid]; urbar=np.mean(ur)
    Qdot=np.array([np.mean(((ur-urbar)-1j*m*radial*omega)*np.exp(-1j*m*ph)) for m in range(1,5)])
    p=abs(np.mean(np.exp(1j*th)))
    return dict(P=p,M_C=mc,s=s,n_valid=int(valid.sum()),n_excluded=int((~valid).sum()),
        min_r_com_over_R=float(rr.min()/R),mean_r_com=float(np.mean(rho)),mean_omega=float(np.mean(omega)),
        Q=Q,A=A,V=V,Qdot=Qdot,valid=valid,phi=phi,alpha=alpha,omega=omega)


def stream_seed(path:Path,stride:int,t0:float=20,t1:float=999.9):
    state,rdot,_=open_state_rdot(path); steps=np.arange(round(t0/DT),round(t1/DT)+1,stride,dtype=np.int64)
    n=len(steps); scalar={k:np.empty(n,np.float32) for k in ['P','M_C','mean_r_com','mean_omega','min_r_com_over_R']}
    excluded=np.empty(n,np.uint8); sgn=np.empty(n,np.int8)
    complex_fields={k:np.empty((n,4),np.complex64) for k in ['Q','A','V','Qdot']}
    for j,step in enumerate(steps):
        fm=frame_modes(np.asarray(state[:,:,step],float),np.asarray(rdot[:,:,step+1],float))
        for k in scalar:scalar[k][j]=fm[k]
        excluded[j]=fm['n_excluded'];sgn[j]=fm['s']
        for k in complex_fields:complex_fields[k][j]=fm[k]
    time=steps*DT; psi=np.zeros(n,float)
    if n>1:psi[1:]=np.cumsum(.5*(scalar['mean_omega'][1:]+scalar['mean_omega'][:-1])*np.diff(time))
    co={k:complex_fields[k]*np.exp(1j*psi[:,None]*np.arange(1,5)[None,:]) for k in ['Q','A','V']}
    return dict(step=steps,time=time,psi=psi.astype(np.float32),s=sgn,n_excluded=excluded,**scalar,**complex_fields,
                Q_co=co['Q'].astype(np.complex64),A_co=co['A'].astype(np.complex64),V_co=co['V'].astype(np.complex64))


def save_modes(path:Path,data:dict):np.savez_compressed(path,**data)


def spectral_summary(seed:int,data:dict,intervals:list[tuple[float,float]],dt:float=.05):
    rows=[]; fs=1/dt
    for idx,(t0,t1) in enumerate(intervals):
        use=(data['time']>=t0)&(data['time']<=t1)
        for mode in range(4):
            z=data['Q_co'][use,mode].astype(complex); z=z-np.mean(z)
            # A full-window periodogram retains the 1/(t1-t0) frequency
            # resolution needed to distinguish the proposed 11.4 time-unit
            # period from nearby oscillations. Welch's default segments were
            # too short for this bounded diagnostic.
            f,S=signal.periodogram(z,fs=fs,detrend='linear',return_onesided=False,scaling='density')
            positive=(f>.02)&(f<.5); ip=np.flatnonzero(positive)[np.argmax(S[positive])]
            target=np.argmin(np.abs(f-REFERENCE_OMEGA/(2*np.pi)))
            rows.append(dict(seed=seed,interval=idx,t0=t0,t1=t1,mode=mode+1,
                mean_abs_Q=float(np.mean(np.abs(data['Q'][use,mode]))),mean_abs_A=float(np.mean(np.abs(data['A'][use,mode]))),
                mean_abs_V=float(np.nanmean(np.abs(data['V'][use,mode]))),peak_frequency=float(f[ip]),
                peak_period=float(1/f[ip]),peak_psd=float(S[ip].real),reference_frequency=REFERENCE_OMEGA/(2*np.pi),
                psd_at_reference=float(S[target].real),reference_to_peak_psd=float(S[target].real/max(S[ip].real,1e-30))))
    return rows


def mode_at_state(state,velocity,s_fixed:int|None=None,cutoff=CUTOFF):
    fm=frame_modes(state,velocity,cutoff); s=fm['s'] if s_fixed is None else s_fixed
    valid=fm['valid']; ph=fm['phi'][valid]; al=fm['alpha'][valid]; v=wrap(al-s*np.pi/2)
    fm['s']=s;fm['V2']=np.mean(v*np.exp(-2j*ph));return fm


def feedback_snapshot(replay:COMReplay,seed:int,step:int,psi:float):
    row,detail=replay.evaluate(seed,step); state_mm,rdot_mm=replay.arrays(seed)
    state=np.asarray(state_mm[:,:,step],float); velocity=detail['effective_velocity']
    fm=mode_at_state(state,velocity); valid=fm['valid']; ph=fm['phi'][valid]; om=fm['omega']; al=fm['alpha'][valid]
    s=fm['s']; vh=wrap(al-s*np.pi/2); basis=np.exp(-2j*ph); phase=np.exp(2j*psi)
    components={name:np.mean(np.asarray(value)[valid]*basis) for name,value in detail['components'].items()}
    geometry=np.mean((-om-2j*vh*om)*basis)
    total=sum(components.values())+geometry
    Q2=fm['Q'][1];V2=fm['V2'];Qdot=fm['Qdot'][1]

    # Native-step identity audit including inferred angular noise. Use a fixed
    # fish set and exclude the rare wrap-branch crossings.
    next_state=np.asarray(state_mm[:,:,step+1],float);next_vel=np.asarray(rdot_mm[:,:,step+2],float)
    x,y=next_state[:2];cx,cy=np.mean(x),np.mean(y);rrn=np.hypot(x-cx,y-cy);common=valid&(rrn>CUTOFF)
    phi0=fm['phi'][common];alpha0=fm['alpha'][common];v0=wrap(alpha0-s*np.pi/2);b0=np.exp(-2j*phi0)
    phin=np.arctan2(y-cy,x-cx);vn=wrap(wrap(next_state[2]-phin)-s*np.pi/2)
    V0=np.mean(v0*b0);Vn=np.mean(vn[common]*np.exp(-2j*phin[common]))
    omega_eff=detail['omega_effective'];noise=wrap(next_state[2]-state[2]-omega_eff*DT)
    ccomp={name:np.mean(np.asarray(value)[common]*b0) for name,value in detail['components'].items()}
    # Direct indexing is clearer and robust when the common mask changes.
    yx=state[0]-np.mean(state[0]);yy=state[1]-np.mean(state[1]);vrel=velocity-np.mean(velocity,axis=1)[:,None]
    omfull=np.full(N,np.nan);omfull[common]=(yx[common]*vrel[1,common]-yy[common]*vrel[0,common])/(yx[common]**2+yy[common]**2)
    geom0=np.mean((-omfull[common]-2j*v0*omfull[common])*b0)
    drift0=sum(ccomp.values())+geom0;noiseproj=np.mean(noise[common]*b0)
    audit_residual=(Vn-V0)-(drift0*DT+noiseproj)
    # Flag a wrapped-coordinate branch whenever the represented increment
    # differs by ~2pi from the locally predicted unwrapped increment. Merely
    # testing |vn-v0|>pi misses a nearly-pi clipped turn that crosses the cut.
    predicted_v_increment=(omega_eff[common]-omfull[common])*DT+noise[common]
    branch=int(np.count_nonzero(np.abs((vn[common]-v0)-predicted_v_increment)>np.pi))
    out={"seed":seed,"step":step,"time":step*DT,"P":row['P'],"Pi":row['Pi'],"M_C":fm['M_C'],"s":s,
        "n_valid":fm['n_valid'],"n_excluded":fm['n_excluded'],"min_r_com_over_R":fm['min_r_com_over_R'],
        "psi":psi,"Q2_real":Q2.real,"Q2_imag":Q2.imag,"V2_real":V2.real,"V2_imag":V2.imag,
        "Q2_co_real":(Q2*phase).real,"Q2_co_imag":(Q2*phase).imag,"V2_co_real":(V2*phase).real,"V2_co_imag":(V2*phase).imag,
        "Q2dot_real":Qdot.real,"Q2dot_imag":Qdot.imag,"V2dot_total_real":total.real,"V2dot_total_imag":total.imag,
        "V2_native_audit_abs_residual":abs(audit_residual),
        "V2_native_audit_branch_safe_abs_residual":abs(audit_residual) if branch==0 else np.nan,
        "V2_native_branch_crossings":branch,
        "velocity_replay_max_error":row['velocity_replay_max_error'],"social_split_max_error":row['social_split_max_error']}
    for name,value in components.items():
        out[f'T2_{name}_real']=value.real;out[f'T2_{name}_imag']=value.imag
        out[f'V2_power_{name}']=float(np.real(np.conj(V2)*value))
    out['T2_geometry_real']=geometry.real;out['T2_geometry_imag']=geometry.imag
    out['V2_power_geometry']=float(np.real(np.conj(V2)*geometry));out['V2_power_total']=float(np.real(np.conj(V2)*total))
    out['Q2_amplitude_rate']=float(np.real(np.conj(Q2)*Qdot)/max(abs(Q2),1e-30))
    return out


def replay_intervals(replay:COMReplay,modes:dict[int,dict],outdir:Path):
    specs={7:[(100,160),(470,530),(800,860)],8:[(100,160),(470,530),(800,860)],9:[(250,360)],0:[(100,160)]}
    rows=[]
    for seed,intervals in specs.items():
        data=modes[seed]
        for interval,(t0,t1) in enumerate(intervals):
            for step in np.arange(round(t0/DT),round(t1/DT)+1,200,dtype=np.int64):
                time=step*DT;psi=float(np.interp(time,data['time'],data['psi']))
                r=feedback_snapshot(replay,seed,int(step),psi);r['interval']=interval;r['t0']=t0;r['t1']=t1;rows.append(r)
        print(f"source replay seed {seed}: {sum(r['seed']==seed for r in rows)} snapshots",flush=True)
    write_csv(outdir/'torque_mode_snapshots.csv',rows)
    return rows


def interval_feedback_summary(rows:list[dict]):
    fields=['alignment','attraction','wall','direct_hydro','bem_hydro','clip_residual','geometry','total']
    out=[]
    for seed in sorted(set(r['seed'] for r in rows)):
        for interval in sorted(set(r['interval'] for r in rows if r['seed']==seed)):
            use=[r for r in rows if r['seed']==seed and r['interval']==interval]
            rec=dict(seed=seed,interval=interval,t0=use[0]['t0'],t1=use[0]['t1'],n=len(use),
                mean_abs_Q2=float(np.mean([abs(complex(r['Q2_real'],r['Q2_imag'])) for r in use])),
                mean_abs_V2=float(np.mean([abs(complex(r['V2_real'],r['V2_imag'])) for r in use])),
                mean_Q2_amplitude_rate=float(np.mean([r['Q2_amplitude_rate'] for r in use])),
                max_velocity_error=max(r['velocity_replay_max_error'] for r in use),
                max_social_error=max(r['social_split_max_error'] for r in use),
                max_native_V2_residual_branch_safe=float(np.nanmax([r['V2_native_audit_branch_safe_abs_residual'] for r in use])),
                total_branch_crossings=sum(r['V2_native_branch_crossings'] for r in use))
            for f in fields:rec['mean_V2_power_'+f]=float(np.mean([r['V2_power_'+f] for r in use]))
            out.append(rec)
    return out


def seed9_window_summary(rows:list[dict],modes:dict):
    """Prespecified chronology around the archived seed-9 transition."""
    windows=[('school',250,300),('late_school',300,306),('circulation_loss',306,311.5),
             ('polarization_collapse',311.5,315),('early_mill',315,330),('mill',330,360)]
    out=[]
    for label,t0,t1 in windows:
        rs=[r for r in rows if r['seed']==9 and t0<=r['time']<t1]
        use=(modes['time']>=t0)&(modes['time']<t1)
        rec=dict(window=label,t0=t0,t1=t1,n_torque=len(rs),
            mean_P=float(np.mean(modes['P'][use])),mean_M_C=float(np.mean(modes['M_C'][use])),
            mean_abs_Q1=float(np.mean(np.abs(modes['Q'][use,0]))),
            mean_abs_Q2=float(np.mean(np.abs(modes['Q'][use,1]))),
            mean_abs_Q3=float(np.mean(np.abs(modes['Q'][use,2]))),
            mean_abs_Q4=float(np.mean(np.abs(modes['Q'][use,3]))),
            mean_abs_V2=float(np.mean(np.abs(modes['V'][use,1]))))
        for f in ('alignment','attraction','wall','direct_hydro','bem_hydro','geometry','total'):
            rec['mean_V2_power_'+f]=float(np.mean([r['V2_power_'+f] for r in rs]))
        rec['mean_Q2_amplitude_rate']=float(np.mean([r['Q2_amplitude_rate'] for r in rs]))
        out.append(rec)
    return out


def cross_phase_summary(rows:list[dict]):
    """Descriptive complex Fourier phase at the proposed and observed Q2 frequencies."""
    out=[]
    keys=('alignment','attraction','wall','direct_hydro','bem_hydro','geometry')
    groups=sorted(set((r['seed'],r['interval']) for r in rows))
    for seed,interval in groups:
        rs=[r for r in rows if r['seed']==seed and r['interval']==interval]
        t=np.array([r['time'] for r in rs]);psi=np.array([r['psi'] for r in rs]);phase=np.exp(2j*psi)
        q=np.array([complex(r['Q2_real'],r['Q2_imag']) for r in rs])*phase
        v=np.array([complex(r['V2_real'],r['V2_imag']) for r in rs])*phase
        sources={k:np.array([complex(r['T2_'+k+'_real'],r['T2_'+k+'_imag']) for r in rs])*phase for k in keys}
        sources['geometry']=np.array([complex(r['T2_geometry_real'],r['T2_geometry_imag']) for r in rs])*phase
        dt=float(np.median(np.diff(t)));f,S=signal.periodogram(q,fs=1/dt,detrend='linear',return_onesided=False)
        ok=(f>.02)&(f<.5);fpeak=float(f[np.flatnonzero(ok)[np.argmax(S[ok])]])
        for which,freq in (('observed_Q2_peak',fpeak),('proposed',REFERENCE_OMEGA/(2*np.pi))):
            carrier=np.exp(-2j*np.pi*freq*(t-t[0]))
            def coef(z):return np.sum(signal.detrend(z,type='linear')*carrier)/len(z)
            fq=coef(q);fv=coef(v)
            rec=dict(seed=seed,interval=interval,t0=rs[0]['t0'],t1=rs[0]['t1'],frequency_kind=which,
                frequency=freq,Q2_fourier_amplitude=float(abs(fq)),V2_fourier_amplitude=float(abs(fv)),
                V2_phase_minus_Q2=float(np.angle(fv*np.conj(fq))))
            for k,z in sources.items():
                fk=coef(z);rec[k+'_fourier_amplitude']=float(abs(fk));rec[k+'_phase_minus_Q2']=float(np.angle(fk*np.conj(fq)))
            out.append(rec)
    return out


def stable_mill_summary(interval_rows:list[dict]):
    """Treat six independent time blocks as the uncertainty units."""
    use=[r for r in interval_rows if r['seed'] in (7,8)]
    fields=['mean_Q2_amplitude_rate']+[f'mean_V2_power_{k}' for k in
        ('alignment','attraction','wall','direct_hydro','bem_hydro','geometry','total')]
    out=[]
    for field in fields:
        x=np.array([r[field] for r in use],float);se=stats.sem(x);ci=stats.t.interval(.95,len(x)-1,loc=np.mean(x),scale=se)
        out.append(dict(metric=field,n_intervals=len(x),mean=float(np.mean(x)),standard_error=float(se),
            ci95_low=float(ci[0]),ci95_high=float(ci[1]),n_positive=int(np.sum(x>0)),n_negative=int(np.sum(x<0))))
    return out


def plots(modes:dict[int,dict],rows:list[dict],outdir:Path):
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=False)
    for ax,seed,window in ((axs[0],7,(100,200)),(axs[1],9,(250,360))):
        d=modes[seed];u=(d['time']>=window[0])&(d['time']<=window[1])
        for m in range(4):ax.plot(d['time'][u],np.abs(d['Q'][u,m])/d['mean_r_com'][u],label=f'm={m+1}',lw=.9)
        ax.set(ylabel=r'$|Q_m|/\bar\rho$',title=f'seed {seed}');ax.grid(alpha=.2);ax.legend(ncol=4)
    axs[1].axvline(312.41,color='black',ls='--',lw=.8);axs[1].set_xlabel('time')
    fig.tight_layout();fig.savefig(outdir/'shape_modes_stable_and_transition.png',dpi=220);plt.close(fig)

    fig,axs=plt.subplots(2,2,figsize=(11,9))
    for row,(seed,window) in zip(axs,[(7,(100,160)),(8,(100,160)),(9,(250,330)),(0,(100,160))]):
        d=modes[seed];u=(d['time']>=window[0])&(d['time']<=window[1])
        q=d['Q_co'][u,1]/np.mean(d['mean_r_com'][u]);v=d['V_co'][u,1]
        row[0].plot(q.real,q.imag,lw=.7);row[0].set_title(f'seed {seed}: co-rotating Q2');row[0].set_aspect('equal',adjustable='datalim')
        row[1].plot(v.real,v.imag,lw=.7);row[1].set_title(f'seed {seed}: co-rotating V2');row[1].set_aspect('equal',adjustable='datalim')
        for ax in row:ax.grid(alpha=.2);ax.set_xlabel('real');ax.set_ylabel('imaginary')
    fig.tight_layout();fig.savefig(outdir/'corotating_Q2_V2.png',dpi=220);plt.close(fig)

    tr=[r for r in rows if r['seed']==9];t=np.array([r['time'] for r in tr])
    fig,axs=plt.subplots(4,1,figsize=(12,12),sharex=True)
    d=modes[9];u=(d['time']>=250)&(d['time']<=360)
    axs[0].plot(d['time'][u],d['P'][u],label='P');axs[0].plot(d['time'][u],d['M_C'][u],label='M_C');axs[0].legend(ncol=2)
    axs[1].plot(d['time'][u],np.abs(d['Q'][u,1])/d['mean_r_com'][u],label='|Q2|/mean rho');axs[1].plot(d['time'][u],np.abs(d['V'][u,1]),label='|V2|');axs[1].legend(ncol=2)
    for f,l in [('alignment','alignment'),('direct_hydro','direct hydro'),('wall','wall'),('geometry','geometry')]:
        axs[2].plot(t,[r['V2_power_'+f] for r in tr],label=l,lw=.8)
    axs[2].legend(ncol=4);axs[2].set_ylabel('V2 modal power')
    axs[3].plot(t,[r['Q2_amplitude_rate'] for r in tr],label='Q2 amplitude rate');axs[3].plot(t,[r['V2_power_total'] for r in tr],label='total V2 power');axs[3].legend(ncol=2)
    for ax in axs:ax.axvline(312.41,color='black',ls='--',lw=.8);ax.grid(alpha=.2)
    axs[3].set_xlabel('time');fig.tight_layout();fig.savefig(outdir/'seed9_mode_feedback_timeline.png',dpi=220);plt.close(fig)

    sm=stable_mill_summary(interval_feedback_summary(rows))
    keep=[r for r in sm if r['metric'].startswith('mean_V2_power_')]
    labels=[r['metric'].removeprefix('mean_V2_power_').replace('_','\n') for r in keep]
    means=np.array([r['mean'] for r in keep]);lo=np.array([r['ci95_low'] for r in keep]);hi=np.array([r['ci95_high'] for r in keep])
    fig,ax=plt.subplots(figsize=(10,5));x=np.arange(len(keep))
    ax.bar(x,means,color=['#4c78a8','#f58518','#54a24b','#e45756','#72b7b2','#b279a2','#777777'])
    ax.errorbar(x,means,yerr=np.vstack((means-lo,hi-means)),fmt='none',ecolor='black',capsize=3,lw=1)
    ax.axhline(0,color='black',lw=.8);ax.set_xticks(x,labels);ax.set_ylabel(r'$\langle\mathrm{Re}(V_2^*D_k)\rangle$')
    ax.set_title('Stable mills: source-resolved heading-mode power\nmean and 95% t interval over six time blocks')
    ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(outdir/'stable_mill_modal_power.png',dpi=220);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--source',type=Path,default=CIRCLE.parents[1].parent/'BEM_circular_parallel_instance.py')
    p.add_argument('--outdir',type=Path,default=HERE/'outputs');args=p.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    files=sorted(args.data_dir.glob(FILE_PATTERN),key=seed_from_name);paths={seed_from_name(f):f for f in files}
    if not all(s in paths for s in (0,7,8,9)):raise FileNotFoundError(args.data_dir)
    modes={}
    for seed in (0,7,8,9):
        modes[seed]=stream_seed(paths[seed],50);save_modes(args.outdir/f'modes_seed{seed}.npz',modes[seed]);print(f'streamed modes seed {seed}',flush=True)
    spectral=[]
    spectral+=spectral_summary(7,modes[7],[(100,160),(470,530),(800,860)],.05)
    spectral+=spectral_summary(8,modes[8],[(100,160),(470,530),(800,860)],.05)
    spectral+=spectral_summary(9,modes[9],[(250,310),(315,360)],.05)
    spectral+=spectral_summary(0,modes[0],[(100,160)],.05);write_csv(args.outdir/'mode_spectral_summary.csv',spectral)

    # Sampling convergence at dt=.025 on one stable mill and the transition.
    conv=[]
    for seed,t0,t1 in ((7,100,200),(9,250,360)):
        fine=stream_seed(paths[seed],25,t0,t1);coarse=modes[seed];use=(coarse['time']>=t0)&(coarse['time']<=t1)
        for m in range(4):
            c=np.abs(coarse['Q'][use,m]);f=np.abs(fine['Q'][::2,m]);n=min(len(c),len(f))
            conv.append(dict(seed=seed,mode=m+1,t0=t0,t1=t1,n_common=n,max_common_frame_difference=float(np.max(np.abs(c[:n]-f[:n]))),
                mean_coarse=float(np.mean(c)),mean_fine=float(np.mean(np.abs(fine['Q'][:,m]))),relative_mean_difference=float((np.mean(np.abs(fine['Q'][:,m]))-np.mean(c))/max(np.mean(c),1e-30))))
    write_csv(args.outdir/'sampling_convergence.csv',conv)

    source=args.source.resolve();params=json.loads((CIRCLE/'simulation_parameters.json').read_text());sha=hashlib.sha256(source.read_bytes()).hexdigest()
    if sha!=params['source_sha256']:raise RuntimeError(f'source hash mismatch {sha}')
    replay=COMReplay(args.data_dir,source,CIRCLE/'outputs');rows=replay_intervals(replay,modes,args.outdir)
    summary=interval_feedback_summary(rows);write_csv(args.outdir/'feedback_interval_summary.csv',summary)
    write_csv(args.outdir/'stable_mill_feedback_summary.csv',stable_mill_summary(summary))
    write_csv(args.outdir/'seed9_feedback_windows.csv',seed9_window_summary(rows,modes[9]))
    write_csv(args.outdir/'cross_phase_summary.csv',cross_phase_summary(rows))
    plots(modes,rows,args.outdir)
    audit=dict(source=str(source),source_sha256=sha,n_torque_snapshots=len(rows),cutoff_over_R=CUTOFF/R,
        excluded_fish_snapshots=sum(r['n_excluded'] for r in rows),minimum_r_com_over_R=min(r['min_r_com_over_R'] for r in rows),
        max_velocity_replay_error=max(r['velocity_replay_max_error'] for r in rows),max_social_split_error=max(r['social_split_max_error'] for r in rows),
        max_native_V2_identity_residual_branch_safe=float(np.nanmax([r['V2_native_audit_branch_safe_abs_residual'] for r in rows])),
        native_V2_branch_crossings=sum(r['V2_native_branch_crossings'] for r in rows),
        co_rotation='psi(t)=integral mean_i omega_i^C dt; tilde Z_m=Z_m exp(+i m psi)',
        exact_V2_drift='mean[(Omega-omega_C)exp(-2i phi_C) - 2i v omega_C exp(-2i phi_C)]')
    (args.outdir/'validation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    config=dict(stage_A_dt=.05,stage_A_fine_dt=.025,stage_B_dt=.2,stage_B_intervals={"7":[[100,160],[470,530],[800,860]],"8":[[100,160],[470,530],[800,860]],"9":[[250,360]],"0":[[100,160]]},
        cutoff_over_R=CUTOFF/R,reference_real_lambda=.00701,reference_imag_lambda=REFERENCE_OMEGA,reference_period=REFERENCE_PERIOD)
    (args.outdir/'analysis_config.json').write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps({'outdir':str(args.outdir),'torque_snapshots':len(rows)},indent=2),flush=True)

if __name__=='__main__':main()
