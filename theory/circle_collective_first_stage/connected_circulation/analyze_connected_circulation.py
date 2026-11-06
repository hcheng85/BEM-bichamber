#!/usr/bin/env python3
"""Connected circulation: streaming geometry, source replay, and balance audit.

The production simulator is imported read-only through the established replay.
No new trajectory is generated. Native-step replays are deliberately bounded.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
sys.path[:0]=[str(BASE),str(BASE/'com_relative_turning')]
from matv5_memmap import open_state_rdot
from process_collective import FILE_PATTERN,seed_from_name
from analyze_com_relative_turning import COMReplay

DT=.001; SIGMA=.5; N=100; R=7.158910531638177
SOURCES=('alignment','attraction','wall','direct_hydro','bem_hydro','clip_residual')

def write_csv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def read_csv(path):
    out=[]
    for row in csv.DictReader(Path(path).open()):
        out.append({k:float(v) if k not in ('group','source','series','criterion','kind') else v for k,v in row.items()})
    return out

def stream(path,stride=50,t0=20.,t1=999.9):
    sm,_,_=open_state_rdot(path);steps=np.arange(round(t0/DT),round(t1/DT)+1,stride)
    keys=('P','MC','G','K','p_real','p_imag','A1_real','A1_imag','C_real','C_imag',
          'com_x','com_y','com_radius_R','mean_rho','com_heading_angle','abs_A1','H2_power','min_rho_R','bound_C')
    out={k:np.empty(len(steps)) for k in keys};err=0.
    for a in range(0,len(steps),512):
        ss=steps[a:a+512];st=np.asarray(sm[:,:,ss],float);x,y,th=st
        cx=x.mean(axis=0);cy=y.mean(axis=0);xx=x-cx;yy=y-cy;rho=np.hypot(xx,yy)
        if np.any(rho==0):raise RuntimeError('Undefined polar angle at COM')
        phi=np.arctan2(yy,xx);e=np.exp(1j*th);f=np.exp(-1j*phi)
        p=e.mean(axis=0);A=f.mean(axis=0);B=(e*f).mean(axis=0);C=B-p*A
        MC=np.sin(th-phi).mean(axis=0);G=(p*A).imag;K=C.imag
        err=max(err,float(np.max(abs(MC-G-K))))
        A2=(f*f).mean(axis=0);H2=(e*f**3).mean(axis=0)-B*A2
        vals=dict(P=abs(p),MC=MC,G=G,K=K,p_real=p.real,p_imag=p.imag,A1_real=A.real,A1_imag=A.imag,
            C_real=C.real,C_imag=C.imag,com_x=cx,com_y=cy,com_radius_R=np.hypot(cx,cy)/R,
            mean_rho=rho.mean(axis=0),com_heading_angle=np.angle(p*np.exp(-1j*np.arctan2(cy,cx))),
            abs_A1=abs(A),H2_power=abs(H2)**2,min_rho_R=rho.min(axis=0)/R,
            bound_C=np.sqrt(np.maximum(0,(1-abs(p)**2)*(1-abs(A)**2))))
        for k in keys:out[k][a:a+len(ss)]=vals[k]
    return dict(step=steps,time=steps*DT,**out),err

def phasor_moments(a):
    """Exact complex sum moments for independent exp(i Normal(0,sigma²dt))."""
    r=np.exp(-.5*SIGMA**2*DT);v1=-np.expm1(-SIGMA**2*DT)
    mu=r*a.sum();v=v1*np.sum(abs(a)**2);w=-r*r*v1*np.sum(a*a)
    third=-r*v1*v1*np.sum(a*abs(a)**2)
    vp=2*abs(mu)**2*v+2*np.real(np.conj(mu)**2*w)+4*np.real(np.conj(mu)*third)+v*v+abs(w)**2-v1**4*np.sum(abs(a)**4)
    return mu,float(v),float(max(0,vp)),float(max(0,.5*(v-w.real)))

def validate_identities(paths,out):
    """Independent finite differences and noise redraws; no BEM solves."""
    rng=np.random.default_rng(56029);rows=[]
    for seed,step in [(7,120000),(8,112400),(9,310000),(9,313000)]:
        sm,rm,_=open_state_rdot(paths[seed]);st=np.asarray(sm[:,:,step],float);vel=np.asarray(rm[:,:,step+1],float)
        th=st[2];pos=st[:2];xx,yy=pos-pos.mean(axis=1)[:,None];phi=np.arctan2(yy,xx)
        vr=vel-vel.mean(axis=1)[:,None];omega=(xx*vr[1]-yy*vr[0])/(xx*xx+yy*yy)
        e=np.exp(1j*th);f=np.exp(-1j*phi);p=e.mean();A=f.mean()
        # Arbitrary deterministic turning probes the analytical derivative;
        # source-specific coefficients subsequently reuse the same identity.
        turn=rng.normal(size=N);analytic=1j*np.mean(e*(f-A)*turn)-1j*np.mean((e-p)*f*omega)
        h=1e-6;cs=[]
        for sign in (-1,1):
            xp=pos+sign*h*vel;xp-=xp.mean(axis=1)[:,None];fp=np.exp(-1j*np.arctan2(xp[1],xp[0]));ep=np.exp(1j*(th+sign*h*turn))
            cs.append((ep*fp).mean()-ep.mean()*fp.mean())
        derivative=(cs[1]-cs[0])/(2*h)
        ac=e*(f-A)/N;ah=e*f*(f*f-(f*f).mean())/N
        eps=SIGMA*np.sqrt(DT)*rng.normal(size=(30000,N));z=np.exp(1j*eps)
        _,_,_,vc=phasor_moments(ac);_,_,vh,_=phasor_moments(ah)
        samplec=(z@ac).imag;sampleh=abs(z@ah)**2
        rows.append(dict(seed=seed,step=step,derivative_error=abs(derivative-analytic),
            C_im_variance=vc,C_im_redraw_variance=float(samplec.var()),H2_power_variance=vh,H2_power_redraw_variance=float(sampleh.var()),redraws=len(z)))
    write_csv(out/'independent_identity_tests.csv',rows)

def finite_geometry_audit(replay,out):
    """A handful of exact split steps to check polar-geometry linearization."""
    rows=[]
    for seed,t in [(9,t) for t in (282.,286.,300.,306.,308.,309.,310.,311.,312.,313.,315.,320.)]+[(7,101.),(8,112.4)]:
        step=round(t/DT);_,d=replay.evaluate(seed,step);sm,_=replay.arrays(seed)
        st=np.asarray(sm[:,:,step],float);nxt=np.asarray(sm[:,:,step+1],float)
        e=np.exp(1j*st[2]);ph=np.arctan2(st[1]-st[1].mean(),st[0]-st[0].mean());f=np.exp(-1j*ph)
        pn=np.arctan2(nxt[1]-nxt[1].mean(),nxt[0]-nxt[0].mean());fn=np.exp(-1j*pn)
        current=np.mean(e*(f-f.mean()));after_position=np.mean(e*(fn-fn.mean()))
        star=np.mean(e*np.exp(1j*DT*d['omega_effective'])*(fn-fn.mean()));decay=np.exp(-.5*SIGMA**2*DT)
        rr=snapshot(replay,seed,step)
        rows.append(dict(seed=seed,time=t,continuous_geometry=rr['bC_geometry'],finite_geometry=(after_position-current).imag/DT,
            continuous_angular=sum(rr['bC_'+k] for k in SOURCES),finite_angular=(star-after_position).imag/DT,
            finite_diffusion=(decay-1)*star.imag/DT,total_finite=(decay*star-current).imag/DT))
    write_csv(out/'finite_geometry_audit.csv',rows)

def snapshot(replay,seed,step):
    base,d=replay.evaluate(seed,step);sm,rdot=replay.arrays(seed)
    st=np.asarray(sm[:,:,step],float);nxt=np.asarray(sm[:,:,step+1],float)
    x,y,theta=st;xx=x-x.mean();yy=y-y.mean();rr2=xx*xx+yy*yy
    phi=np.arctan2(yy,xx);e=np.exp(1j*theta);f=np.exp(-1j*phi)
    p=e.mean();A=f.mean();B=(e*f).mean();C=B-p*A;u=e*f;f2=f*f;A2=f2.mean();H=(u*f2).mean()-B*A2
    vr=d['effective_velocity']-d['effective_velocity'].mean(axis=1)[:,None]
    omega=(xx*vr[1]-yy*vr[0])/rr2
    dc={k:1j*np.mean(e*(f-A)*d['components'][k]) for k in SOURCES}
    dc['geometry']=-1j*np.mean((e-p)*f*omega);dc['diffusion']=-.5*SIGMA**2*C
    hp={k:2*np.real(np.conj(H)*1j*np.mean(u*(f2-A2)*d['components'][k])) for k in SOURCES}
    dhg=-3j*np.mean(u*f2*omega)+1j*A2*np.mean(u*omega)+2j*B*np.mean(omega*f2)
    hp['geometry']=float(2*np.real(np.conj(H)*dhg));hp['diffusion']=-SIGMA**2*abs(H)**2
    hp['injection']=SIGMA**2/N*(1-abs(A2)**2)
    phin=np.arctan2(nxt[1]-nxt[1].mean(),nxt[0]-nxt[0].mean());fn=np.exp(-1j*phin)
    en=np.exp(1j*nxt[2]);estar=np.exp(1j*(theta+d['omega_effective']*DT))
    cn=(en*fn).mean()-en.mean()*fn.mean();hn=(en*fn**3).mean()-(en*fn).mean()*(fn*fn).mean()
    ac=estar*(fn-fn.mean())/N;ah=estar*fn*(fn*fn-(fn*fn).mean())/N
    muc,vc,_,vki=phasor_moments(ac);muh,vh,vhp,_=phasor_moments(ah)
    next_h_power=abs(muh)**2+vh;expected_h=(next_h_power-abs(H)**2)/DT
    r=dict(seed=seed,step=step,time=step*DT,P=abs(p),MC=B.imag,G=(p*A).imag,K=C.imag,
        H2_power=abs(H)**2,min_rho_R=np.sqrt(rr2.min())/R,
        C_finite_drift=(muc.imag-C.imag)/DT,C_observed_delta=cn.imag-C.imag,
        C_noise_delta=cn.imag-muc.imag,C_noise_variance=vki,
        H2_finite_drift=expected_h,H2_observed_delta=abs(hn)**2-abs(H)**2,
        H2_noise_delta=abs(hn)**2-next_h_power,H2_noise_variance=vhp,
        velocity_error=base['velocity_replay_max_error'],social_error=base['social_split_max_error'],
        position_error=float(np.max(abs(nxt[:2]-st[:2]-DT*d['effective_velocity']))),
        clip_fish=base['n_angular_clip'],clip_abs_angle=float(np.sum(abs(d['components']['clip_residual']))*DT),
        sticky_fish=base['n_sticky_override'],collision_fish=base['n_collision_override'])
    for k,z in dc.items():r['bC_'+k]=float(z.imag)
    r['bC_total']=float(sum(dc.values()).imag)
    for k,z in hp.items():r['bH2_'+k]=float(z)
    r['bH2_total']=float(sum(hp.values()))
    r['C_discretization_drift']=r['C_finite_drift']-r['bC_total'];r['H2_discretization_drift']=r['H2_finite_drift']-r['bH2_total']
    # Identity between decomposed full circulation and existing replay.
    r['MC_drift_replay_error']=abs(base['bMC_total']-(np.mean(1j*e*f*(d['omega_effective']-omega))-.5*SIGMA**2*B).imag)
    replay.cache.clear() # avoid retaining every 100-fish replay in memory
    return r

def replay_grid(replay,out):
    # Same stable blocks as the old analysis at a 20x finer grid. A bounded
    # native-grid portion in each seed supplies a true telescoping balance.
    specs=[('mill7',7,100.,160.,10),('mill8',8,100.,160.,10),
           ('transition',9,300.,330.,20),('school0',0,100.,130.,20),('school9',9,282.,288.4,20),
           ('native7',7,100.,102.,1),('native8',8,100.,102.,1)]
    rows=[];seen=set();path=out/'replay.csv';start=time.monotonic()
    if path.exists():
        rows=read_csv(path);seen={(int(r['seed']),int(r['step'])) for r in rows}
    for label,seed,t0,t1,stride in specs:
        for step in range(round(t0/DT),round(t1/DT),stride):
            if (seed,step) in seen:continue
            row=snapshot(replay,seed,step);rows.append(row);seen.add((seed,step))
            if len(rows)%300==0:
                write_csv(path,rows);print(f'{label}: {len(rows)} snapshots; elapsed {time.monotonic()-start:.0f}s',flush=True)
        write_csv(path,rows);print(f'completed {label}',flush=True)
    return rows,specs

def geometry_analysis(series,fine,out):
    # Equal-duration, nonoverlapping candidate controls. Match only the initial
    # state; stability over the window defines a reference, not a prospective
    # classifier cohort. Reflect chirality using the initial MC sign only.
    target={k:float(np.interp(306.,fine['time'],fine[k])) for k in ('P','MC','com_radius_R','mean_rho')}
    refs={k:np.concatenate([series[s][k][(series[s]['time']>=100)&(series[s]['time']<900)] for s in range(7)]) for k in target}
    scales={k:max(float(np.std(v)),1e-8) for k,v in refs.items()}
    controls=[]
    for seed in range(7):
        d=series[seed];candidates=[]
        for t in np.arange(100.,880.,10.):
            idx=(d['time']>=t)&(d['time']<=t+6.4)
            if np.min(d['P'][idx])<.8:continue # define stable controls, report conditioning
            initial={k:float(np.interp(t,d['time'],d[k])) for k in target}
            mirror=1 if initial['MC']*target['MC']>=0 else -1
            initial['MC']*=mirror
            dist=sum(((initial[k]-target[k])/scales[k])**2 for k in target)**.5
            candidates.append((dist,t,initial,mirror))
        if not candidates:continue
        dist,t,initial,mirror=min(candidates,key=lambda z:z[0]);r=dict(seed=seed,t0=t,t1=t+6.4,mirror=mirror,match_distance=dist)
        for k in ('P','MC','G','K','com_radius_R','mean_rho'):
            mult=mirror if k in ('MC','G','K') else 1
            r[k+'_initial']=mult*float(np.interp(t,d['time'],d[k]));r[k+'_delta']=mult*float(np.interp(t+6.4,d['time'],d[k]))-r[k+'_initial']
        controls.append(r)
    write_csv(out/'matched_controls.csv',controls)
    windows=[(250,300),(300,306),(306,310),(310,312),(312,312.4),(312.4,315),(315,330)]
    wins=[]
    for a,b in windows:
        r=dict(t0=a,t1=b)
        for k in ('P','MC','G','K','com_radius_R','mean_rho'):
            r[k+'_initial']=float(np.interp(a,fine['time'],fine[k]));r[k+'_final']=float(np.interp(b,fine['time'],fine[k]));r[k+'_delta']=r[k+'_final']-r[k+'_initial']
        wins.append(r)
    write_csv(out/'seed9_changes.csv',wins)
    # Threshold chronology with persistence 0.2. No smoothing or fitted changepoint.
    times=[]
    for field,cut,sign in [('P',.8,-1),('P',.7,-1),('P',.5,-1),('K',-.05,-1),('K',-.1,-1),('K',-.2,-1),('MC',0.,-1),('MC',-.2,-1),('G',.05,-1),('G',0.,-1)]:
        mask=(fine[field]<cut) if sign<0 else fine[field]>cut
        valid=np.convolve(mask.astype(int),np.ones(20,dtype=int),mode='valid')==20
        ix=np.flatnonzero(valid&(fine['time'][:len(valid)]>=300)&(fine['time'][:len(valid)]<=330))
        times.append(dict(series=field,threshold=cut,persistence=.2,time=float(fine['time'][ix[0]]) if len(ix) else np.nan))
    write_csv(out/'chronology.csv',times)
    onset=[]
    for field in ('K','G','P'):
        base=float(np.interp(306.,fine['time'],fine[field]))
        for drop in (.025,.05,.1):
            mask=fine[field]<=base-drop
            valid=np.convolve(mask.astype(int),np.ones(20,dtype=int),mode='valid')==20
            ix=np.flatnonzero(valid&(fine['time'][:len(valid)]>=306)&(fine['time'][:len(valid)]<=330))
            onset.append(dict(series=field,baseline_time=306.,drop=drop,persistence=.2,time=float(fine['time'][ix[0]]) if len(ix) else np.nan))
    write_csv(out/'equal_drop_onsets.csv',onset)
    # Seven independent control seeds; bootstrap only these seed-level deltas.
    rng=np.random.default_rng(68031);unc=[]
    for k in ('P','MC','G','K'):
        vals=np.array([r[k+'_delta'] for r in controls]);boots=rng.choice(vals,(10000,len(vals)),replace=True).mean(axis=1)
        event=float(np.interp(312.4,fine['time'],fine[k])-np.interp(306,fine['time'],fine[k]))
        unc.append(dict(series=k,event_delta=event,control_mean=float(vals.mean()),control_min=float(vals.min()),control_max=float(vals.max()),
            control_mean_ci_low=float(np.quantile(boots,.025)),control_mean_ci_high=float(np.quantile(boots,.975)),
            event_minus_control_ci_low=float(event-np.quantile(boots,.975)),event_minus_control_ci_high=float(event-np.quantile(boots,.025)),n_control_seeds=len(vals),n_event_seeds=1))
    write_csv(out/'control_uncertainty.csv',unc)
    # Additional descriptive specificity check over all disjoint reference
    # windows, not just the best match. Thresholds are the chronology levels.
    screen=[]
    for seed in range(7):
        d=series[seed]
        for cutoff in (-.05,-.1,-.2):
            eligible=hits=0
            for a in np.arange(100.,880.,10.):
                u=(d['time']>=a)&(d['time']<=a+6.4)
                if np.min(d['P'][u])<.8:continue
                k=d['K'][u];sign=1 if k[0]>=0 else -1;k=k*sign
                if k[0]<.1:continue
                eligible+=1;hits+=int(np.any(np.convolve((k<cutoff).astype(int),np.ones(4,dtype=int),mode='valid')==4))
            screen.append(dict(seed=seed,threshold=cutoff,n_stable_windows=eligible,n_reversals=hits))
    write_csv(out/'stable_school_reversal_screen.csv',screen)
    return controls

def summarize_replay(rows,out):
    budgets=[]
    for label,seed,a,b in [('mill7',7,100,160),('mill8',8,100,160),('school0',0,100,130),('school9',9,282,288.4),
                         ('pre1',9,300,306),('pre2',9,306,310),('pre3',9,310,312.4),('collapse',9,312.4,315),('post',9,315,330)]:
        grid=10 if seed in (7,8) else 20
        use=[r for r in rows if int(r['seed'])==seed and a<=r['time']<b and int(r['step'])%grid==0]
        rec=dict(group=label,seed=seed,t0=a,t1=b,n=len(use))
        for k in [*['bC_'+s for s in (*SOURCES,'geometry','diffusion','total')],'C_finite_drift','C_discretization_drift']:
            rec[k]=float(np.mean([r[k] for r in use]))
        budgets.append(rec)
    write_csv(out/'connected_source_budget.csv',budgets)
    # Same-grid offset sensitivity is a sampling diagnostic, not independent
    # seed uncertainty. Native 2-unit blocks also permit exact endpoint balance.
    balance=[]
    for seed in (7,8):
        for a,b,finegrid in [(100,160,10),(100,102,1)]:
            full=[r for r in rows if int(r['seed'])==seed and a<=r['time']<b and int(r['step'])%finegrid==0]
            for stride in ([200,50,10] if finegrid==10 else [200,50,10,1]):
                offset_means=[]
                for offset in range(0,stride,finegrid):
                    sub=[r for r in full if int(r['step'])%stride==offset]
                    if sub:offset_means.append(np.mean([r['H2_finite_drift'] for r in sub]))
                chosen=[r for r in full if int(r['step'])%stride==0]
                rec=dict(seed=seed,t0=a,t1=b,sample_dt=stride*DT,n=len(chosen),
                    mean_H2_continuous=float(np.mean([r['bH2_total'] for r in chosen])),
                    mean_H2_exact=float(np.mean([r['H2_finite_drift'] for r in chosen])),
                    mean_C_exact=float(np.mean([r['C_finite_drift'] for r in chosen])),
                    offset_min=float(np.min(offset_means)),offset_max=float(np.max(offset_means)),
                    actual_endpoint_H2_rate=np.nan,noise_H2_rate=np.nan,noise_H2_rate_se=np.nan,
                    actual_endpoint_C_rate=np.nan,noise_C_rate=np.nan,noise_C_rate_se=np.nan)
                if stride==1:
                    T=b-a;rec['actual_endpoint_H2_rate']=sum(r['H2_observed_delta'] for r in chosen)/T
                    rec['noise_H2_rate']=sum(r['H2_noise_delta'] for r in chosen)/T
                    rec['noise_H2_rate_se']=np.sqrt(sum(r['H2_noise_variance'] for r in chosen))/T
                    rec['actual_endpoint_C_rate']=sum(r['C_observed_delta'] for r in chosen)/T
                    rec['noise_C_rate']=sum(r['C_noise_delta'] for r in chosen)/T
                    rec['noise_C_rate_se']=np.sqrt(sum(r['C_noise_variance'] for r in chosen))/T
                balance.append(rec)
    write_csv(out/'sampling_balance.csv',balance)
    return budgets,balance

def stationary_geometry(series,rows,out):
    records=[]
    for seed in (7,8):
        d=series[seed]
        for a,b in [(100,160),(100,500),(500,900)]:
            use=(d['time']>=a)&(d['time']<b)
            rec=dict(seed=seed,t0=a,t1=b,mean_power=float(d['H2_power'][use].mean()),
                endpoint_rate=float((np.interp(b,d['time'],d['H2_power'])-np.interp(a,d['time'],d['H2_power']))/(b-a)),
                mean_P=float(d['P'][use].mean()),mean_MC=float(d['MC'][use].mean()))
            records.append(rec)
    write_csv(out/'stationary_occupation.csv',records)
    blocks=[]
    for seed in (7,8):
        for a in np.arange(100,160,10):
            use=[r for r in rows if int(r['seed'])==seed and a<=r['time']<a+10 and int(r['step'])%10==0]
            blocks.append(dict(seed=seed,t0=a,t1=a+10,n=len(use),mean_exact_H2=float(np.mean([r['H2_finite_drift'] for r in use])),
                mean_continuous_H2=float(np.mean([r['bH2_total'] for r in use])),mean_exact_C=float(np.mean([r['C_finite_drift'] for r in use]))))
    write_csv(out/'stationary_blocks.csv',blocks)

def make_plots(series,fine,controls,budgets,balance,out):
    fig,axes=plt.subplots(2,1,figsize=(11,5),sharex=True);d=series[9]
    for k in ('MC','G','K'):axes[0].plot(d['time'],d[k],label=k,lw=.6)
    axes[0].legend();axes[0].set_ylabel('signed circulation')
    axes[1].plot(d['time'],d['P'],label='P',lw=.6);axes[1].plot(d['time'],d['com_radius_R'],label='COM radius / R',lw=.6)
    axes[1].legend();axes[1].set_xlabel('time')
    for ax in axes:ax.axvspan(300,330,color='gray',alpha=.2);ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'seed9_full_record.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(4,1,figsize=(11,10),sharex=True)
    for k in ('MC','G','K'):axes[0].plot(fine['time'],fine[k],label={'MC':'total circulation','G':'polarization × shape','K':'connected circulation'}[k])
    axes[0].legend();axes[0].set_ylabel('signed circulation')
    axes[1].plot(fine['time'],fine['P'],label='P');axes[1].plot(fine['time'],fine['abs_A1'],label='|A1|');axes[1].legend()
    axes[2].plot(fine['time'],fine['com_radius_R'],label='COM radius / R');axes[2].plot(fine['time'],fine['mean_rho']/R,label='mean school radius / R');axes[2].legend()
    axes[3].plot(fine['time'],np.cos(fine['com_heading_angle']));axes[3].set_ylabel('cos(heading − COM bearing)');axes[3].set_xlabel('time')
    for ax in axes:ax.axvline(312.41,color='k',ls='--',lw=.7);ax.grid(alpha=.2);ax.set_xlim(300,330)
    fig.tight_layout();fig.savefig(out/'seed9_circulation_geometry.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4));t=np.arange(0,6.401,.05)
    for ax,k in zip(axes,('P','G','K')):
        for r in controls:
            d=series[int(r['seed'])];vals=np.interp(r['t0']+t,d['time'],d[k])*(r['mirror'] if k in ('G','K') else 1);ax.plot(t,vals-vals[0],color='gray',alpha=.5,lw=.8)
        vals=np.interp(306+t,fine['time'],fine[k]);ax.plot(t,vals-vals[0],color='red',lw=2,label='seed 9 pre-collapse')
        ax.set(title=k,xlabel='time after matched start',ylabel='change from initial state');ax.grid(alpha=.2)
    axes[0].legend();fig.tight_layout();fig.savefig(out/'matched_control_changes.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    src=(*SOURCES,'geometry','diffusion');labels=[s.replace('_','\n') for s in src]
    for ax,group in zip(axes,('pre2','pre3')):
        row=next(r for r in budgets if r['group']==group);ax.bar(np.arange(len(src)),[row['bC_'+s] for s in src]);ax.set_xticks(np.arange(len(src)),labels,rotation=35)
        ax.axhline(0,color='black',lw=.6);ax.set(title=f"seed 9: {row['t0']}–{row['t1']}",ylabel='drift of connected circulation')
    fig.tight_layout();fig.savefig(out/'connected_drift_sources.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,seed in zip(axes,(7,8)):
        subset=sorted([r for r in balance if r['seed']==seed and r['t1']==160],key=lambda r:r['sample_dt'])
        xx=np.array([r['sample_dt'] for r in subset]);yy=np.array([r['mean_H2_exact'] for r in subset])
        ax.plot(xx,yy,'o-',label='exact conditional drift, offset 0')
        ax.fill_between(xx,[r['offset_min'] for r in subset],[r['offset_max'] for r in subset],alpha=.2,label='range across grid offsets')
        ax.axhline(0,color='black',lw=.6);ax.set(xscale='log',xlabel='sampling interval',ylabel='mean H2 power drift',title=f'seed {seed}: t=100–160')
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.tight_layout();fig.savefig(out/'stationary_sampling_audit.png',dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);ap.add_argument('--stage',choices=('geometry','replay','summary','all'),default='all');args=ap.parse_args()
    out=HERE/'outputs';out.mkdir(parents=True,exist_ok=True)
    paths={seed_from_name(p):p for p in args.data_dir.glob(FILE_PATTERN)}
    if set(paths)!=set(range(10)):raise FileNotFoundError(f'Missing raw seeds: {set(range(10))-set(paths)} in {args.data_dir}')
    series={};audit={}
    if args.stage in ('geometry','all'):
        maxerr=0.
        for seed,path in sorted(paths.items()):
            d,err=stream(path);series[seed]=d;maxerr=max(maxerr,err);np.savez_compressed(out/f'series_seed{seed}.npz',**d);print(f'streamed seed {seed}',flush=True)
        fine,err=stream(paths[9],10,250,360);np.savez_compressed(out/'seed9_fine.npz',**fine)
        audit=dict(max_circulation_identity_error=max(maxerr,err),min_rho_R=min(d['min_rho_R'].min() for d in series.values()))
        (out/'geometry_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    else:
        series={s:dict(np.load(out/f'series_seed{s}.npz')) for s in range(10)};fine=dict(np.load(out/'seed9_fine.npz'))
    controls=geometry_analysis(series,fine,out)
    if args.stage in ('geometry','all'):validate_identities(paths,out)
    if args.stage=='geometry':return
    if args.stage in ('replay','all'):
        source=BASE.parents[1].parent/'BEM_circular_parallel_instance.py';sha=hashlib.sha256(source.read_bytes()).hexdigest()
        assert sha==json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256']
        replay=COMReplay(args.data_dir,source,BASE/'outputs');rows,specs=replay_grid(replay,out);finite_geometry_audit(replay,out)
    else:rows=read_csv(out/'replay.csv')
    budgets,balance=summarize_replay(rows,out);stationary_geometry(series,rows,out);make_plots(series,fine,controls,budgets,balance,out)
    (out/'replay_audit.json').write_text(json.dumps(dict(n_snapshots=len(rows),max_velocity_error=max(r['velocity_error'] for r in rows),
        max_position_error=max(r['position_error'] for r in rows),max_social_error=max(r['social_error'] for r in rows),
        max_MC_drift_replay_error=max(r['MC_drift_replay_error'] for r in rows),
        min_replay_rho_R=min(r['min_rho_R'] for r in rows),
        clipped_fish_snapshots=sum(r['clip_fish'] for r in rows),
        clipping_snapshots=sum(r['clip_fish']>0 for r in rows),
        sticky_fish_snapshots=sum(r['sticky_fish'] for r in rows),
        collision_fish_snapshots=sum(r['collision_fish'] for r in rows),
        max_C_step_identity_residual=max(abs(r['C_observed_delta']-DT*r['C_finite_drift']-r['C_noise_delta']) for r in rows),
        max_H2_step_identity_residual=max(abs(r['H2_observed_delta']-DT*r['H2_finite_drift']-r['H2_noise_delta']) for r in rows)),indent=2)+'\n')
    config=dict(N=N,R=R,In=SIGMA,dt=DT,raw_file_pattern=FILE_PATTERN,
        raw_directory=str(args.data_dir),geometry_stride=50,seed9_fine_stride=10,
        bootstrap_seed=68031,bootstrap_draws=10000,noise_redraw_seed=56029,
        production_source_sha256=json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256'],
        replay_source='theory/circle_collective_first_stage/com_relative_turning/analyze_com_relative_turning.py',
        snapshot_schedule='See replay_grid: fixed 18920 distinct archived states; no simulated trajectories',
        controls='Seven independent seed-level stable-window contrasts; one transition; no population CI for transitions')
    (out/'analysis_configuration.json').write_text(json.dumps(config,indent=2)+'\n')

if __name__=='__main__':main()
