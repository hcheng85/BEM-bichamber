#!/usr/bin/env python3
"""Bounded source-exact test of COM transport and K-coordinate sufficiency."""
from __future__ import annotations
from pathlib import Path
import argparse,ast,copy,csv,hashlib,json,multiprocessing as mp,sys,time
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1];PILOT=HERE.parent/'stochastic_branching'
sys.path[:0]=[str(BASE),str(BASE/'com_relative_turning'),str(HERE.parent/'same_noise_counterfactual'),str(PILOT)]
from analyze_com_relative_turning import COMReplay
from replay_torques import load_production
from run_branching import coordinate,metric,exact_constraints,exact_jac,min_k_perturbation,save_csv

DT=.001;IN=.5;N=100;R=7.158910531638177;STEPS=1000;TIMES=(309.,310.);REPS=4;NDIR=3;RMS_FACTOR=1.5
BRANCHES=('baseline','baseline_duplicate')+tuple(f'{s}_d{j}' for s in ('plus','minus') for j in range(NDIR))
METRICS=('K','P','MC','com_x','com_y','com_radius_R','mean_radius','radial_std','sigma_major','sigma_minor','aspect')
SOURCES=('swimming','direct_hydro','bem_hydro','override_residual','actual')
_MOD=_PARA=_STEP=None

def source_step_components(mod,source,out):
    """Extract the production step, adding only diagnostic return values."""
    tree=ast.parse(source.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='step_Computation')
    loop=next(x for x in fn.body if isinstance(x,ast.For));first=next(i for i,x in enumerate(loop.body) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='CurState')
    body=copy.deepcopy(loop.body[first:]);replacements=0
    for node in body:
        if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='wNoise':
            node.value=ast.Name(id='noise_increment',ctx=ast.Load());replacements+=1
    assert replacements==1
    prefix=ast.parse('''def archived_step_components(current, para, noise_increment, dis_last):
    state = np.empty((3, para["num"], 2))
    state[:,:,0] = current
    rDotList = np.empty((2, para["num"], 2))
    count = 0
''').body[0]
    suffix=ast.parse('''
swim_mean = np.mean(e, axis=1)
direct_mean = np.mean(U, axis=1)
bem_mean = np.mean(U_BEM, axis=1)
actual_mean = np.mean(rDot, axis=1)
override_mean = actual_mean - swim_mean - direct_mean - bem_mean
return state[:,:,count].copy(), rDot.copy(), thetaDot.copy(), dis_last.copy(), dict(wall_clip=int(np.sum(abs(wWall)>=np.pi/para["dt"])), angular_clip=int(np.sum(abs(thetaDot)>=np.pi/para["dt"])), collision=int(if_jump.sum()), safety_candidates=int(np.sum(dis>=safedis))), np.vstack([swim_mean,direct_mean,bem_mean,override_mean,actual_mean])
''').body
    prefix.body+=body+suffix;module=ast.fix_missing_locations(ast.Module(body=[prefix],type_ignores=[]));text=ast.unparse(module)+'\n'
    (out/'extracted_production_step_components.py').write_text(text);env=dict(vars(mod));exec(compile(module,str(source)+'::archived_step_components','exec'),env)
    return env['archived_step_components'],hashlib.sha256(text.encode()).hexdigest()

def fixed_rms_k_perturbation(st,target,rms,raw):
    """Exact p,K,RMS perturbation nearest an outcome-independent raw direction."""
    theta=st[2];p0,K0,f=coordinate(st);dmin,_=min_k_perturbation(st,target);J=exact_jac(dmin,theta,f)
    B=np.vstack([J,dmin]);_,s,vh=np.linalg.svd(B,full_matrices=False);rank=np.sum(s>s[0]*1e-12)
    h=raw-vh[:rank].T@(vh[:rank]@raw);h/=np.linalg.norm(h)
    extra=np.sqrt(max(N*rms*rms-dmin@dmin,0.));q=dmin+extra*h;r2=N*rms*rms
    def c(x):return np.r_[exact_constraints(x,theta,f,p0,K0+target),x@x-r2]
    def cj(x):return np.vstack([exact_jac(x,theta,f),2*x])
    z=minimize(lambda x:.5*np.sum((x-q)**2),q,jac=lambda x:x-q,constraints={'type':'eq','fun':c,'jac':cj},method='SLSQP',options=dict(ftol=1e-13,maxiter=4000,disp=False))
    err=float(np.max(abs(c(z.x))))
    if not z.success or err>1e-10:raise RuntimeError(f'fixed-RMS solve failed: {z.message}, error={err}')
    return z.x,dict(iterations=z.nit,objective=z.fun,constraint_error=err)

def prepare(replay,step,out):
    sm,_=replay.arrays(9);rows=[];distances=[];states={};memory={};vectors={}
    for t0 in TIMES:
        n=round(t0/DT);st=np.asarray(sm[:,:,n],float);zero=np.zeros((3,N));_,_,_,dis0,_,_=step(np.asarray(sm[:,:,n-1],float),replay.para,zero,None)
        mins={sgn:min_k_perturbation(st,.02 if sgn=='plus' else -.02)[0] for sgn in ('plus','minus')}
        rms=RMS_FACTOR*max(np.sqrt(np.mean(x*x)) for x in mins.values());deltas={'baseline':np.zeros(N),'baseline_duplicate':np.zeros(N)}
        # Fixed seeds and raw directions are independent of every trajectory outcome.
        for j in range(NDIR):
            raw=np.random.default_rng(2026110700+int(t0)*10+j).normal(size=N)
            for sgn,target in (('plus',.02),('minus',-.02)):
                d,a=fixed_rms_k_perturbation(st,target,rms,raw);name=f'{sgn}_d{j}';deltas[name]=d
                vectors[f't{int(t0)}_{name}']=d
                p0,K0,_=coordinate(st);x=st.copy();x[2]+=d;p,K,_=coordinate(x)
                rows.append(dict(time=t0,branch=name,direction=j,target_deltaK=target,deltaK=K-K0,p_real_error=p.real-p0.real,p_imag_error=p.imag-p0.imag,rms_heading=np.sqrt(np.mean(d*d)),max_abs_heading=np.max(abs(d)),**a))
        branches={}
        for name,d in deltas.items():x=st.copy();x[2]+=d;branches[name]=x
        for sign in ('plus','minus'):
            for a in range(NDIR):
                for b in range(a+1,NDIR):
                    da=deltas[f'{sign}_d{a}'];db=deltas[f'{sign}_d{b}'];distances.append(dict(time=t0,sign=sign,direction_a=a,direction_b=b,rms_separation=np.sqrt(np.mean((da-db)**2)),cosine_similarity=(da@db)/np.sqrt((da@da)*(db@db))))
        states[t0]=branches;memory[t0]=dis0
    save_csv(out/'perturbations.csv',rows);save_csv(out/'perturbation_distances.csv',distances);np.savez_compressed(out/'perturbation_vectors.npz',**vectors)
    return states,memory,rows

def init_worker(source,extracted):
    global _MOD,_PARA,_STEP
    _MOD=load_production(Path(source));_PARA=dict(Ip=4.,In=IN,R=_MOD.R,light=1.,dt=DT,Iw=.94,If=.01,num=N,delta=.01,total_time=1000,gap=_MOD.gap,radius_ratio=50,gap_ratio=50)
    env=vars(_MOD);exec(Path(extracted).read_text(),env);_STEP=env['archived_step_components']

def task_worker(task):
    t0,rep,noise_seed,states,dis0=task;rng=np.random.default_rng(noise_seed);noise=IN*np.sqrt(DT)*rng.normal(size=(STEPS,N));start=time.monotonic()
    current={k:v.copy() for k,v in states.items()};memory={k:dis0.copy() for k in BRANCHES}
    data=np.empty((len(BRANCHES),STEPS+1,len(METRICS)));vel=np.empty((len(BRANCHES),STEPS,len(SOURCES),2));flags={k:dict(wall_clip=0,angular_clip=0,collision=0,safety_candidates=0) for k in BRANCHES}
    for j,k in enumerate(BRANCHES):data[j,0]=metric(current[k])
    duplicate_max=0.
    for n in range(STEPS):
        inc=np.zeros((3,N));inc[2]=noise[n]
        for j,k in enumerate(BRANCHES):
            current[k],_,_,memory[k],fl,vel[j,n]=_STEP(current[k],_PARA,inc,memory[k]);data[j,n+1]=metric(current[k])
            for q in flags[k]:flags[k][q]+=fl[q]
        duplicate_max=max(duplicate_max,float(np.max(abs(current['baseline']-current['baseline_duplicate']))))
    return dict(t0=t0,replicate=rep,noise_seed=noise_seed,data=data,velocities=vel,duplicate_max=duplicate_max,flags=flags,seconds=time.monotonic()-start)

def instantaneous_rows(replay,step,out,new_states,new_mem):
    sm,_=replay.arrays(9);zero=np.zeros((3,N));rows=[]
    old=np.load(PILOT/'outputs/perturbation_vectors.npz')
    for t0 in TIMES:
        st=np.asarray(sm[:,:,round(t0/DT)],float);p0,K0,_=coordinate(st)
        old_names=('baseline','plusK','minusK','control_plus','control_minus')
        families=[('existing',[(name,np.zeros(N) if name=='baseline' else old[f't{int(t0)}_{name}']) for name in old_names]),('sameK',[(name,new_states[t0][name][2]-st[2]) for name in BRANCHES if name!='baseline_duplicate'])]
        for family,items in families:
            vals={}
            for name,d in items:
                x=st.copy();x[2]+=d;_,_,_,_,fl,v=step(x,replay.para,zero,new_mem[t0]);p,K,_=coordinate(x);vals[name]=(p,K,v,fl)
            base=vals['baseline']
            for name,(p,K,v,fl) in vals.items():
                c=st[:2].mean(axis=1);er=c/max(np.linalg.norm(c),1e-30);dv=v-base[2]
                row=dict(family=family,time=t0,branch=name,deltaK=K-K0,p_real_error=p.real-p0.real,p_imag_error=p.imag-p0.imag)
                for si,s in enumerate(SOURCES):
                    row[f'delta_{s}_vx']=dv[si,0];row[f'delta_{s}_vy']=dv[si,1];row[f'delta_{s}_vr']=dv[si]@er
                row.update({f'flag_{k}':z for k,z in fl.items()});rows.append(row)
    save_csv(out/'instantaneous_velocity.csv',rows)

def decompose(data,vel):
    """Exact vector and scalar-radius displacement decomposition."""
    c=data[:,[3,4]];increments=DT*vel[:,:4,:];vector=increments.sum(axis=0)
    r0=np.linalg.norm(c[:-1],axis=1);r1=np.linalg.norm(c[1:],axis=1);g=(c[:-1]+c[1:])/np.maximum((r0+r1)[:,None],1e-30)
    radial=np.einsum('tsi,ti->s',increments,g)
    measured=c[-1]-c[0];measured_radial=np.linalg.norm(c[-1])-np.linalg.norm(c[0])
    return vector,radial,measured,measured_radial

def radial_step_components(data,vel):
    """Per-step exact radial displacement assigned to vector components."""
    c=data[:,[3,4]];r0=np.linalg.norm(c[:-1],axis=1);r1=np.linalg.norm(c[1:],axis=1);g=(c[:-1]+c[1:])/np.maximum((r0+r1)[:,None],1e-30)
    return np.einsum('tsi,ti->ts',DT*vel[:,:4,:],g)

def summarize(out):
    files=sorted(out.glob('trajectory_t*_r*.npz'));runs={};rows=[];instant=[]
    for f in files:
        z=np.load(f,allow_pickle=False);t0=float(z['t0']);rep=int(z['replicate']);data=z['data'];vel=z['velocities'];runs[(t0,rep)]=(data,vel)
        for j,b in enumerate(BRANCHES):
            vec,rad,meas,measrad=decompose(data[j],vel[j]);rec=vec.sum(axis=0)
            r=dict(time=t0,replicate=rep,branch=b,initial_K=data[j,0,0],final_K=data[j,-1,0],initial_P=data[j,0,1],final_P=data[j,-1,1],K_auc=float(np.trapezoid(data[j,:,0],dx=DT)),measured_dx=meas[0],measured_dy=meas[1],measured_dr=measrad,reconstruction_x=rec[0],reconstruction_y=rec[1],reconstruction_dr=rad.sum(),residual_norm=np.linalg.norm(meas-rec),residual_radial=measrad-rad.sum())
            for q,s in enumerate(SOURCES[:4]):r[f'{s}_dx']=vec[q,0];r[f'{s}_dy']=vec[q,1];r[f'{s}_dr']=rad[q]
            rows.append(r)
    save_csv(out/'transport_decomposition.csv',rows)
    paired=[]
    for t0 in TIMES:
        for rep in range(REPS):
            data,vel=runs[(t0,rep)];base=next(x for x in rows if x['time']==t0 and x['replicate']==rep and x['branch']=='baseline')
            for j,b in enumerate(BRANCHES[2:],2):
                x=next(x for x in rows if x['time']==t0 and x['replicate']==rep and x['branch']==b);sgn=b.split('_')[0];r=dict(time=t0,replicate=rep,branch=b,sign=sgn,direction=int(b[-1]),final_delta_K=data[j,-1,0]-data[0,-1,0],final_delta_P=data[j,-1,1]-data[0,-1,1],K_auc_delta=float(np.trapezoid(data[j,:,0]-data[0,:,0],dx=DT)))
                for k in ('measured_dx','measured_dy','measured_dr')+tuple(f'{s}_{z}' for s in SOURCES[:4] for z in ('dx','dy','dr')):r['delta_'+k]=x[k]-base[k]
                paired.append(r)
    save_csv(out/'paired_responses.csv',paired)
    ds=[]
    for t0 in TIMES:
        for sign in ('plus','minus'):
            for direction in range(NDIR):
                use=[x for x in paired if x['time']==t0 and x['sign']==sign and x['direction']==direction];r=dict(time=t0,sign=sign,direction=direction,n=len(use))
                for k in ('final_delta_K','final_delta_P','K_auc_delta','delta_measured_dx','delta_measured_dy','delta_measured_dr','delta_swimming_dr','delta_direct_hydro_dr','delta_bem_hydro_dr','delta_override_residual_dr'):
                    v=np.array([x[k] for x in use]);r[k+'_mean']=v.mean();r[k+'_sd']=v.std(ddof=1)
                ds.append(r)
    save_csv(out/'direction_summary.csv',ds)
    suff=[]
    for t0 in TIMES:
        for sign in ('plus','minus'):
            use=[x for x in ds if x['time']==t0 and x['sign']==sign];r=dict(time=t0,sign=sign,n_directions=NDIR,n_noise=REPS)
            for k in ('final_delta_K','final_delta_P','K_auc_delta','delta_measured_dr','delta_swimming_dr','delta_direct_hydro_dr','delta_bem_hydro_dr'):
                means=np.array([x[k+'_mean'] for x in use]);sds=np.array([x[k+'_sd'] for x in use]);r[k+'_direction_sd']=means.std(ddof=1);r[k+'_mean_noise_sd']=sds.mean();r[k+'_direction_range']=means.max()-means.min();r[k+'_grand_mean']=means.mean()
            suff.append(r)
    save_csv(out/'sufficiency_summary.csv',suff)
    cumulative={};time_rows=[];reversal=[]
    for t0 in TIMES:
        for sign in ('plus','minus'):
            by_rep=[]
            for rep in range(REPS):
                data,vel=runs[(t0,rep)];base_steps=radial_step_components(data[0],vel[0]);direction_curves=[]
                for j,b in enumerate(BRANCHES[2:],2):
                    if not b.startswith(sign):continue
                    direction_curves.append(np.cumsum(radial_step_components(data[j],vel[j])-base_steps,axis=0))
                by_rep.append(np.mean(direction_curves,axis=0))
            arr=np.asarray(by_rep);cumulative[(t0,sign)]=arr
            total=arr[:,:,:4].sum(axis=2).mean(axis=0);final_sign=np.sign(total[-1]);opposite=np.where(total*final_sign<0)[0];first_final=np.where(total*final_sign>0)[0]
            reversal.append(dict(time=t0,sign=sign,first_step_delta_radius=total[0],final_delta_radius=total[-1],first_final_sign_time=(first_final[0]+1)*DT if len(first_final) else np.nan,last_opposite_sign_time=(opposite[-1]+1)*DT if len(opposite) else np.nan))
            for tau in (.1,.25,.5,1.):
                n=int(round(tau/DT))-1
                for q,s in enumerate(SOURCES[:4]):time_rows.append(dict(time=t0,sign=sign,tau=tau,source=s,mean=arr[:,n,q].mean(),sd=arr[:,n,q].std(ddof=1)))
    save_csv(out/'cumulative_transport_summary.csv',time_rows)
    save_csv(out/'radial_reversal_summary.csv',reversal)
    audit=json.loads((out/'simulation_audit.json').read_text());flag_totals={k:0 for k in ('wall_clip','angular_clip','collision','safety_candidates')}
    for run in audit['runs']:
        for b,f in run['flags'].items():
            if b=='baseline_duplicate':continue
            for k in flag_totals:flag_totals[k]+=f[k]
    (out/'reconstruction_audit.json').write_text(json.dumps(dict(max_vector_residual=max(x['residual_norm'] for x in rows),max_radial_residual=max(abs(x['residual_radial']) for x in rows),duplicate_max=max(x['duplicate_max'] for x in audit['runs']),nonduplicate_flag_totals=flag_totals),indent=2)+'\n')
    # Figures
    tt=np.arange(STEPS+1)*DT;fig,ax=plt.subplots(2,2,figsize=(11,7),sharex='col')
    for row,t0 in enumerate(TIMES):
        vals=[]
        for rep in range(REPS):
            data,_=runs[(t0,rep)]
            for j,b in enumerate(BRANCHES[2:],2):ax[row,0].plot(tt,data[j,:,0]-data[0,:,0],color='C3' if b.startswith('plus') else 'C0',alpha=.22,lw=.8)
        ax[row,0].axhline(0,color='k',lw=.5);ax[row,0].set(title=f't0={t0:g}: K - baseline',ylabel='Delta K');ax[row,0].grid(alpha=.2)
        for sign,color in (('plus','C3'),('minus','C0')):
            for direction in range(NDIR):
                u=[x for x in paired if x['time']==t0 and x['sign']==sign and x['direction']==direction]
                ax[row,1].scatter([x['K_auc_delta'] for x in u],[x['delta_measured_dr'] for x in u],color=color,marker=('o','s','^')[direction],alpha=.8,label=f'{sign} d{direction}' if row==0 else None)
        ax[row,1].axhline(0,color='k',lw=.5);ax[row,1].axvline(0,color='k',lw=.5);ax[row,1].set(title=f't0={t0:g}: integrated response',xlabel='integral Delta K dt',ylabel='Delta COM radius');ax[row,1].grid(alpha=.2)
    ax[1,0].set_xlabel('elapsed time');ax[0,1].legend(ncol=2,fontsize=8);fig.tight_layout();fig.savefig(out/'sameK_responses.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(11,8),sharey='row')
    for col,t0 in enumerate(TIMES):
        for row,sign in enumerate(('plus','minus')):
            u=[x for x in ds if x['time']==t0 and x['sign']==sign];x=np.arange(NDIR);bottom=np.zeros(NDIR)
            for key,label,color in (('delta_swimming_dr_mean','swimming','C2'),('delta_direct_hydro_dr_mean','direct hydro','C1'),('delta_bem_hydro_dr_mean','BEM hydro','C4'),('delta_override_residual_dr_mean','override residual','C7')):
                v=np.array([z[key] for z in u]);ax[row,col].bar(x,v,bottom=bottom,label=label if col==0 and row==0 else None,color=color);bottom+=v
            meas=np.array([z['delta_measured_dr_mean'] for z in u]);ax[row,col].scatter(x,meas,color='k',zorder=4,label='measured' if col==0 and row==0 else None);ax[row,col].axhline(0,color='k',lw=.5);ax[row,col].set(title=f't0={t0:g}, {sign}K',xticks=x,ylabel='Delta COM radius');ax[row,col].grid(axis='y',alpha=.2)
    for a in ax[1]:a.set_xlabel('prespecified direction')
    ax[0,0].legend(fontsize=8);fig.subplots_adjust(hspace=.32,wspace=.16);fig.savefig(out/'radial_transport_decomposition.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(11,7),sharex=True)
    tt=np.arange(1,STEPS+1)*DT
    for col,t0 in enumerate(TIMES):
        for row,sign in enumerate(('plus','minus')):
            arr=cumulative[(t0,sign)]
            for q,(label,color) in enumerate((('swimming','C2'),('direct hydro','C1'),('BEM hydro','C4'))):
                mean=arr[:,:,q].mean(axis=0);ax[row,col].plot(tt,mean,label=label if row==0 and col==0 else None,color=color)
            total=arr[:,:,:4].sum(axis=2).mean(axis=0);ax[row,col].plot(tt,total,color='k',lw=1.3,label='total' if row==0 and col==0 else None)
            ax[row,col].axhline(0,color='k',lw=.5);ax[row,col].set(title=f't0={t0:g}, {sign}K',ylabel='cumulative Delta COM radius');ax[row,col].grid(alpha=.2)
    for a in ax[1]:a.set_xlabel('elapsed time')
    ax[0,0].legend(fontsize=8);fig.tight_layout();fig.savefig(out/'cumulative_radial_transport.png',dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--stage',choices=('all','setup','simulate','summary'),default='all');args=ap.parse_args()
    out=HERE/'outputs';out.mkdir(parents=True,exist_ok=True);source=BASE.parents[1].parent/'BEM_circular_parallel_instance.py';sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha==json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256']
    replay=COMReplay(args.data_dir,source,BASE/'outputs');step,step_sha=source_step_components(replay.mod,source,out);states,mem,pert=prepare(replay,step,out);instantaneous_rows(replay,step,out,states,mem)
    env=dict(vars(replay.mod));exec((PILOT/'outputs/extracted_production_step.py').read_text(),env);old_step=env['archived_step'];equiv=[];zero=np.zeros((3,N))
    for t0 in TIMES:
        old=old_step(states[t0]['baseline'],replay.para,zero,mem[t0]);new=step(states[t0]['baseline'],replay.para,zero,mem[t0])
        equiv.append(dict(time=t0,state_max=np.max(abs(old[0]-new[0])),velocity_max=np.max(abs(old[1]-new[1])),turning_max=np.max(abs(old[2]-new[2])),memory_max=np.max(abs(old[3]-new[3])),flags_equal=old[4]==new[4]))
    (out/'source_equivalence.json').write_text(json.dumps(equiv,indent=2)+'\n')
    if any(max(x['state_max'],x['velocity_max'],x['turning_max'],x['memory_max'])>0 or not x['flags_equal'] for x in equiv):raise RuntimeError('diagnostic production-step extraction changed the update')
    (out/'analysis_configuration.json').write_text(json.dumps(dict(source_sha256=sha,step_sha256=step_sha,dt=DT,steps=STEPS,replicates=REPS,directions_per_sign=NDIR,rms_factor=RMS_FACTOR,initial_times=TIMES,branches=BRANCHES,noise='fresh PCG64 Gaussian; identical across branches within each block',scope='one-time-unit same-K sufficiency and exact COM transport decomposition'),indent=2)+'\n')
    if args.stage=='setup':return
    tasks=[(t0,rep,960000+int(t0)*100+rep,states[t0],mem[t0]) for t0 in TIMES for rep in range(REPS)]
    if args.stage in ('all','simulate'):
        todo=[x for x in tasks if not (out/f'trajectory_t{int(x[0])}_r{x[1]}.npz').exists()]
        if todo:
            start=time.monotonic();aud=[];ctx=mp.get_context('spawn')
            with ctx.Pool(min(args.workers,len(todo)),initializer=init_worker,initargs=(str(source),str(out/'extracted_production_step_components.py'))) as pool:
                for z in pool.imap_unordered(task_worker,todo):
                    np.savez_compressed(out/f'trajectory_t{int(z["t0"])}_r{z["replicate"]}.npz',t0=z['t0'],replicate=z['replicate'],noise_seed=z['noise_seed'],data=z['data'],velocities=z['velocities'],metrics=np.array(METRICS),sources=np.array(SOURCES),branches=np.array(BRANCHES));aud.append({k:v for k,v in z.items() if k not in ('data','velocities','flags')}|{'flags':z['flags']});print(f'finished t={z["t0"]}, rep={z["replicate"]}, duplicate={z["duplicate_max"]:.3g}',flush=True)
            old=[] if not (out/'simulation_audit.json').exists() else json.loads((out/'simulation_audit.json').read_text())['runs'];runs=old+aud;(out/'simulation_audit.json').write_text(json.dumps(dict(source_sha256=sha,step_sha256=step_sha,workers=args.workers,wall_seconds=time.monotonic()-start,runs=runs),indent=2)+'\n')
        else:runs=json.loads((out/'simulation_audit.json').read_text())['runs']
        if any(z['duplicate_max']>1e-13 for z in runs):raise RuntimeError('duplicate validation failed')
    if args.stage=='simulate':return
    summarize(out)

if __name__=='__main__':main()
