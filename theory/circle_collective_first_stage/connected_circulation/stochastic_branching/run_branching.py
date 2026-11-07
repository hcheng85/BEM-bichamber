#!/usr/bin/env python3
"""Bounded source-exact stochastic branching from archived circle states."""
from __future__ import annotations
from pathlib import Path
import argparse,csv,hashlib,json,multiprocessing as mp,sys,time
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
sys.path[:0]=[str(BASE),str(BASE/'com_relative_turning'),str(HERE.parent/'same_noise_counterfactual')]
from analyze_com_relative_turning import COMReplay
from replay_torques import load_production
from replay_baseline import source_step

DT=.001;IN=.5;N=100;R=7.158910531638177;STEPS=1000
TIMES=(309.,310.);BRANCHES=('baseline','baseline_duplicate','plusK','minusK','control_plus','control_minus')
METRICS=('K','P','MC','com_x','com_y','com_radius_R','mean_radius','radial_std','sigma_major','sigma_minor','aspect')
_MOD=_PARA=_STEP=None

def save_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def coordinate(st):
    xy=st[:2]-st[:2].mean(axis=1)[:,None];phi=np.arctan2(xy[1],xy[0]);e=np.exp(1j*st[2]);f=np.exp(-1j*phi)
    p=e.mean();K=np.mean(e*(f-f.mean())).imag
    return p,float(K),f

def metric(st):
    xy=st[:2]-st[:2].mean(axis=1)[:,None];rho=np.linalg.norm(xy,axis=0);phi=np.arctan2(xy[1],xy[0]);e=np.exp(1j*st[2]);f=np.exp(-1j*phi)
    p=e.mean();B=np.mean(e*f);K=(B-p*f.mean()).imag;cov=xy@xy.T/N;ev=np.linalg.eigvalsh(cov)
    return np.array([K,abs(p),B.imag,*st[:2].mean(axis=1),np.linalg.norm(st[:2].mean(axis=1))/R,rho.mean(),rho.std(),np.sqrt(ev[1]),np.sqrt(ev[0]),np.sqrt(ev[1]/max(ev[0],1e-30))])

def exact_constraints(delta,theta,f,p0,Ktarget):
    e=np.exp(1j*(theta+delta));p=e.mean();K=np.mean(e*(f-f.mean())).imag
    return np.array([p.real-p0.real,p.imag-p0.imag,K-Ktarget])

def exact_jac(delta,theta,f):
    e=np.exp(1j*(theta+delta));return np.vstack([-e.imag/N,e.real/N,np.real(e*(f-f.mean()))/N])

def min_k_perturbation(st,target):
    theta=st[2];p0,K0,f=coordinate(st);J=exact_jac(np.zeros(N),theta,f)
    A=J[:2];_,s,vh=np.linalg.svd(A,full_matrices=False);rank=np.sum(s>s[0]*1e-12)
    g=J[2];h=g-vh[:rank].T@(vh[:rank]@g);x0=target*h/(g@h)
    fun=lambda x:.5*x@x;jac=lambda x:x
    con={'type':'eq','fun':lambda x:exact_constraints(x,theta,f,p0,K0+target),'jac':lambda x:exact_jac(x,theta,f)}
    z=minimize(fun,x0,jac=jac,constraints=con,method='SLSQP',options=dict(ftol=1e-13,maxiter=2000,disp=False))
    if not z.success or np.max(abs(con['fun'](z.x)))>1e-11:raise RuntimeError(f'K solve failed: {z.message}, {con["fun"](z.x)}')
    return z.x,dict(iterations=z.nit,objective=z.fun,constraint_error=float(np.max(abs(con['fun'](z.x)))))

def exact_control(st,rms,raw):
    theta=st[2];p0,K0,f=coordinate(st);J=exact_jac(np.zeros(N),theta,f)
    _,s,vh=np.linalg.svd(J,full_matrices=False);rank=np.sum(s>s[0]*1e-12)
    q=raw-vh[:rank].T@(vh[:rank]@raw);q*=rms/np.sqrt(np.mean(q*q));r2=N*rms*rms
    def c(x):return np.r_[exact_constraints(x,theta,f,p0,K0),x@x-r2]
    def cj(x):return np.vstack([exact_jac(x,theta,f),2*x])
    z=minimize(lambda x:.5*np.sum((x-q)**2),q,jac=lambda x:x-q,constraints={'type':'eq','fun':c,'jac':cj},
        method='SLSQP',options=dict(ftol=1e-13,maxiter=3000,disp=False))
    if not z.success or np.max(abs(c(z.x)))>1e-10:raise RuntimeError(f'control solve failed: {z.message}, {c(z.x)}')
    return z.x,dict(iterations=z.nit,objective=z.fun,constraint_error=float(np.max(abs(c(z.x)))))

def init_worker(source,extracted):
    global _MOD,_PARA,_STEP
    _MOD=load_production(Path(source));_PARA=dict(Ip=4.,In=IN,R=_MOD.R,light=1.,dt=DT,Iw=.94,If=.01,num=N,delta=.01,total_time=1000,gap=_MOD.gap,radius_ratio=50,gap_ratio=50)
    env=vars(_MOD);exec(Path(extracted).read_text(),env);_STEP=env['archived_step']

def task_worker(task):
    t0,rep,noise_seed,states,dis0=task;rng=np.random.default_rng(noise_seed);noise=IN*np.sqrt(DT)*rng.normal(size=(STEPS,N))
    current={k:v.copy() for k,v in states.items()};memory={k:dis0.copy() for k in BRANCHES}
    data=np.empty((len(BRANCHES),STEPS+1,len(METRICS)));flags={k:dict(wall_clip=0,angular_clip=0,collision=0,safety_candidates=0) for k in BRANCHES}
    for j,k in enumerate(BRANCHES):data[j,0]=metric(current[k])
    duplicate_max=0.;start=time.monotonic()
    for n in range(STEPS):
        inc=np.zeros((3,N));inc[2]=noise[n]
        for j,k in enumerate(BRANCHES):
            current[k],_,_,memory[k],fl=_STEP(current[k],_PARA,inc,memory[k]);data[j,n+1]=metric(current[k])
            for x in flags[k]:flags[k][x]+=fl[x]
        duplicate_max=max(duplicate_max,float(np.max(abs(current['baseline']-current['baseline_duplicate']))))
    return dict(t0=t0,replicate=rep,noise_seed=noise_seed,data=data,duplicate_max=duplicate_max,flags=flags,seconds=time.monotonic()-start)

def prepare(replay,step_func,out):
    sm,_=replay.arrays(9);rng=np.random.default_rng(20261106);rows=[];states_by_time={};dis_by_time={};arrays={}
    for t0 in TIMES:
        step=round(t0/DT);st=np.asarray(sm[:,:,step],float);zero=np.zeros((3,N))
        _,_,_,dis0,_=step_func(np.asarray(sm[:,:,step-1],float),replay.para,zero,None)
        plus,pa=min_k_perturbation(st,.02);minus,ma=min_k_perturbation(st,-.02)
        cp,cpa=exact_control(st,np.sqrt(np.mean(plus*plus)),rng.normal(size=N));cm,cma=exact_control(st,np.sqrt(np.mean(minus*minus)),rng.normal(size=N))
        deltas=dict(baseline=np.zeros(N),baseline_duplicate=np.zeros(N),plusK=plus,minusK=minus,control_plus=cp,control_minus=cm)
        p0,K0,_=coordinate(st);branches={}
        for name,d in deltas.items():
            x=st.copy();x[2]+=d;branches[name]=x;p,K,_=coordinate(x)
            rows.append(dict(time=t0,branch=name,deltaK=K-K0,p_real_error=p.real-p0.real,p_imag_error=p.imag-p0.imag,
                rms_heading=float(np.sqrt(np.mean(d*d))),max_abs_heading=float(np.max(abs(d))),
                first_order_dK=float(exact_jac(np.zeros(N),st[2],coordinate(st)[2])[2]@d),
                iterations=np.nan,objective=np.nan,constraint_error=np.nan))
            arrays[f't{int(t0)}_{name}']=d
        rows[-4].update(**pa);rows[-3].update(**ma);rows[-2].update(**cpa);rows[-1].update(**cma)
        states_by_time[t0]=branches;dis_by_time[t0]=dis0
    save_csv(out/'perturbations.csv',rows);np.savez_compressed(out/'perturbation_vectors.npz',**arrays)
    return states_by_time,dis_by_time,rows

def summarize(out,perturb_rows):
    files=sorted(out.glob('trajectory_t*_r*.npz'));records=[];trajectories={}
    for file in files:
        z=np.load(file,allow_pickle=False);t0=float(z['t0']);rep=int(z['replicate']);data=z['data'];trajectories[(t0,rep)]=data
        for j,b in enumerate(BRANCHES):
            records.append(dict(time=t0,replicate=rep,branch=b,initial_K=data[j,0,0],final_K=data[j,-1,0],initial_P=data[j,0,1],final_P=data[j,-1,1],
                initial_MC=data[j,0,2],final_MC=data[j,-1,2],initial_com_radius_R=data[j,0,5],final_com_radius_R=data[j,-1,5],
                initial_mean_radius=data[j,0,6],final_mean_radius=data[j,-1,6],initial_aspect=data[j,0,10],final_aspect=data[j,-1,10]))
    save_csv(out/'endpoint_metrics.csv',records)
    paired=[]
    for t0 in TIMES:
        for rep in range(4):
            d=trajectories[(t0,rep)];base=d[0]
            for j,b in enumerate(BRANCHES[2:],2):
                r=dict(time=t0,replicate=rep,branch=b)
                for q,k in enumerate(METRICS):r['final_delta_'+k]=d[j,-1,q]-base[-1,q]
                r['K_gain']=(d[j,-1,0]-base[-1,0])/(d[j,0,0]-base[0,0]) if abs(d[j,0,0]-base[0,0])>1e-12 else np.nan
                r['K_auc']=float(np.trapezoid(d[j,:,0]-base[:,0],dx=DT));paired.append(r)
    save_csv(out/'paired_responses.csv',paired)
    # Equal-RMS controls distinguish a selective K displacement from the effect
    # of applying an arbitrary polarization-preserving heading perturbation.
    adjusted=[]
    for t0 in TIMES:
        for rep in range(4):
            d=trajectories[(t0,rep)]
            for target,control in (('plusK','control_plus'),('minusK','control_minus')):
                jt=BRANCHES.index(target);jc=BRANCHES.index(control)
                r=dict(time=t0,replicate=rep,contrast=f'{target}_minus_{control}')
                for q,k in enumerate(METRICS):r['final_delta_'+k]=d[jt,-1,q]-d[jc,-1,q]
                r['K_auc']=float(np.trapezoid(d[jt,:,0]-d[jc,:,0],dx=DT));adjusted.append(r)
    save_csv(out/'control_adjusted_responses.csv',adjusted)

    # Fixed prospective times make persistence/decay visible without choosing a
    # response time after inspecting each realization.
    timepoints=[]
    for t0 in TIMES:
        for rep in range(4):
            d=trajectories[(t0,rep)];base=d[0]
            for tau in (0.,.1,.25,.5,1.):
                n=int(round(tau/DT))
                for j,b in enumerate(BRANCHES[2:],2):
                    timepoints.append(dict(time=t0,replicate=rep,tau=tau,branch=b,
                        delta_K=d[j,n,0]-base[n,0],delta_P=d[j,n,1]-base[n,1],delta_MC=d[j,n,2]-base[n,2]))
    save_csv(out/'timepoint_responses.csv',timepoints)
    sums=[]
    for t0 in TIMES:
        for b in BRANCHES[2:]:
            use=[x for x in paired if x['time']==t0 and x['branch']==b]
            r=dict(time=t0,branch=b,n=len(use))
            for k in ('final_delta_K','final_delta_P','final_delta_MC','final_delta_com_radius_R','final_delta_mean_radius','final_delta_aspect','K_gain','K_auc'):
                v=np.array([x[k] for x in use]);finite=v[np.isfinite(v)]
                r[k+'_mean']=float(np.mean(finite)) if len(finite) else np.nan
                r[k+'_sd']=float(np.std(finite,ddof=1)) if len(finite)>1 else np.nan
            sums.append(r)
    save_csv(out/'summary.csv',sums)
    adj_sums=[]
    for t0 in TIMES:
        for contrast in ('plusK_minus_control_plus','minusK_minus_control_minus'):
            use=[x for x in adjusted if x['time']==t0 and x['contrast']==contrast]
            r=dict(time=t0,contrast=contrast,n=len(use))
            for k in ('final_delta_K','final_delta_P','final_delta_MC','final_delta_com_radius_R','final_delta_mean_radius','final_delta_aspect','K_auc'):
                v=np.array([x[k] for x in use]);r[k+'_mean']=float(v.mean());r[k+'_sd']=float(v.std(ddof=1))
            adj_sums.append(r)
    save_csv(out/'control_adjusted_summary.csv',adj_sums)

    tp_sums=[]
    for t0 in TIMES:
        for tau in (0.,.1,.25,.5,1.):
            for b in BRANCHES[2:]:
                use=[x for x in timepoints if x['time']==t0 and x['tau']==tau and x['branch']==b]
                r=dict(time=t0,tau=tau,branch=b,n=len(use))
                for k in ('delta_K','delta_P','delta_MC'):
                    v=np.array([x[k] for x in use]);r[k+'_mean']=float(v.mean());r[k+'_sd']=float(v.std(ddof=1))
                tp_sums.append(r)
    save_csv(out/'timepoint_summary.csv',tp_sums)
    # Small decisive figures: paired time traces and endpoint distributions.
    tt=np.arange(STEPS+1)*DT;fig,ax=plt.subplots(2,3,figsize=(13,7),sharex=True)
    for row,t0 in enumerate(TIMES):
        for col,q in enumerate((0,1,2)):
            for rep in range(4):
                d=trajectories[(t0,rep)];
                for j,b in enumerate(('plusK','minusK','control_plus','control_minus'),2):
                    ax[row,col].plot(t0+tt,d[j,:,q]-d[0,:,q],color={'plusK':'C3','minusK':'C0','control_plus':'C1','control_minus':'C2'}[b],alpha=.35,lw=.8)
            ax[row,col].axhline(0,color='k',lw=.5);ax[row,col].set(title=f't0={t0:g}, '+METRICS[q],ylabel='branch − same-noise baseline');ax[row,col].grid(alpha=.2)
    for a in ax[-1]:a.set_xlabel('time')
    fig.tight_layout();fig.savefig(out/'paired_time_evolution.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for a,t0 in zip(ax,TIMES):
        for i,b in enumerate(BRANCHES[2:]):
            vals=[x['final_delta_K'] for x in paired if x['time']==t0 and x['branch']==b];a.scatter(np.full(4,i),vals);a.plot([i-.15,i+.15],[np.mean(vals)]*2,color='k')
        a.set_xticks(range(4),BRANCHES[2:],rotation=25);a.axhline(0,color='k',lw=.5);a.set(title=f't0={t0:g}',ylabel='final K − baseline K');a.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(out/'endpoint_K_responses.png',dpi=180);plt.close(fig)
    labels=('K','P','MC','COM radius / R','mean radius','aspect')
    keys=('final_delta_K','final_delta_P','final_delta_MC','final_delta_com_radius_R','final_delta_mean_radius','final_delta_aspect')
    fig,ax=plt.subplots(2,3,figsize=(12,7))
    for a,label,key in zip(ax.flat,labels,keys):
        x=0
        for t0,color in zip(TIMES,('C4','C5')):
            for contrast in ('plusK_minus_control_plus','minusK_minus_control_minus'):
                vals=np.array([r[key] for r in adjusted if r['time']==t0 and r['contrast']==contrast])
                a.scatter(np.full(len(vals),x),vals,color=color,alpha=.8)
                a.plot([x-.13,x+.13],[vals.mean()]*2,color='k',lw=1.2);x+=1
        a.axhline(0,color='k',lw=.5);a.set_title(label);a.grid(axis='y',alpha=.2)
        a.set_xticks(range(4),('309 +K','309 -K','310 +K','310 -K'),rotation=25)
    fig.suptitle('Targeted perturbation minus its equal-RMS p,K-preserving control')
    fig.tight_layout();fig.savefig(out/'control_adjusted_endpoints.png',dpi=180);plt.close(fig)
    return sums,paired

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--stage',choices=('all','setup','simulate','summary'),default='all');args=ap.parse_args()
    out=HERE/'outputs';out.mkdir(parents=True,exist_ok=True);source=BASE.parents[1].parent/'BEM_circular_parallel_instance.py'
    sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha==json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256']
    replay=COMReplay(args.data_dir,source,BASE/'outputs');step,step_sha=source_step(replay.mod,source,out)
    states,mem,pert=prepare(replay,step,out)
    if args.stage=='setup':return
    tasks=[]
    for t0 in TIMES:
        for rep in range(4):tasks.append((t0,rep,910000+int(t0)*100+rep,states[t0],mem[t0]))
    if args.stage in ('all','simulate'):
        todo=[x for x in tasks if not (out/f'trajectory_t{int(x[0])}_r{x[1]}.npz').exists()]
        if todo:
            ctx=mp.get_context('spawn');start=time.monotonic();aud=[]
            with ctx.Pool(min(args.workers,len(todo)),initializer=init_worker,initargs=(str(source),str(out/'extracted_production_step.py'))) as pool:
                for result in pool.imap_unordered(task_worker,todo):
                    file=out/f'trajectory_t{int(result["t0"])}_r{result["replicate"]}.npz'
                    np.savez_compressed(file,t0=result['t0'],replicate=result['replicate'],noise_seed=result['noise_seed'],data=result['data'],metrics=np.array(METRICS),branches=np.array(BRANCHES))
                    aud.append({k:v for k,v in result.items() if k not in ('data','flags')}|{'flags':result['flags']})
                    print(f'finished t={result["t0"]}, replicate={result["replicate"]}, duplicate={result["duplicate_max"]:.3g}',flush=True)
            old=[] if not (out/'simulation_audit.json').exists() else json.loads((out/'simulation_audit.json').read_text())['runs']
            allruns=old+aud;(out/'simulation_audit.json').write_text(json.dumps(dict(source_sha256=sha,step_sha256=step_sha,steps_per_branch=STEPS,
                branches=list(BRANCHES),workers=args.workers,wall_seconds=time.monotonic()-start,runs=allruns),indent=2)+'\n')
        else:
            allruns=json.loads((out/'simulation_audit.json').read_text())['runs']
        if any(r['duplicate_max']>1e-13 for r in allruns):raise RuntimeError('Duplicate baseline validation failed')
    if args.stage=='simulate':return
    sums,paired=summarize(out,pert)
    (out/'analysis_configuration.json').write_text(json.dumps(dict(source_sha256=sha,step_sha256=step_sha,dt=DT,steps=STEPS,replicates=4,
        initial_times=TIMES,branches=BRANCHES,noise='fresh PCG64 Gaussian increments; shared across branches within replicate',
        interpretation='bounded one-time-unit branching; not archived-transition replay or transition probability'),indent=2)+'\n')

if __name__=='__main__':main()
