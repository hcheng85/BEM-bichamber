#!/usr/bin/env python3
"""Bounded archived-state heading perturbations; production code is read-only."""
from pathlib import Path
import argparse,csv,hashlib,json,sys
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
sys.path[:0]=[str(BASE),str(BASE/'com_relative_turning')]
from analyze_com_relative_turning import COMReplay
from replay_torques import social_split

DT=.001; IN=.5; IP=4.; N=100
SRC=('alignment','attraction','wall','direct_hydro','bem_hydro','clip_residual','geometry','diffusion')

def save(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def quantities(st):
    theta=st[2];pos=st[:2]-st[:2].mean(axis=1)[:,None]
    phi=np.arctan2(pos[1],pos[0]);e=np.exp(1j*theta);f=np.exp(-1j*phi)
    p=e.mean();A=f.mean();C=np.mean(e*(f-A));Theta=np.angle(p)
    delta=np.angle(e*np.exp(-1j*Theta));z=np.cos(Theta-phi)
    approx=np.mean((delta-delta.mean())*(z-z.mean()))
    grad=np.real(e*(f-A))/N;long=np.array([np.cos(Theta),np.sin(Theta)])@pos
    return dict(e=e,f=f,p=p,A=A,K=C.imag,grad=grad,Theta=Theta,delta=delta,
                approx=approx,long=long,pos=pos)

def weights(st,vn):
    # Same angular weights as getangle/GetVision; no hard field-of-view cutoff.
    dx=st[0][None,:]-st[0][:,None];dy=st[1][None,:]-st[1][:,None]
    beta=np.arctan2(dy,dx)-st[2][:,None]
    w=vn*(1+np.cos(beta));den=w.sum(axis=1)
    if np.any(den<=0):raise ValueError('Zero visual denominator')
    W=w/den[:,None]
    dw=vn*np.sin(beta);dW=(dw-W*dw.sum(axis=1)[:,None])/den[:,None]
    diff=st[2][None,:]-st[2][:,None]
    align=IP*np.sum(W*np.sin(diff),axis=1)
    J=IP*W*np.cos(diff)
    J[np.diag_indices(N)]=-J.sum(axis=1)+IP*np.sum(dW*np.sin(diff),axis=1)
    a=np.stack([np.sum(W*dx,axis=1),np.sum(W*dy,axis=1)])
    return W,align,J,a

def evaluate(replay,st):
    mod,para=replay.mod,replay.para;q=quantities(st)
    vn=mod.GetVoronoi(st[:2],para);W,align,J,a=weights(st,vn)
    alignment,attraction=social_split(mod,vn,st,IP)
    visual=mod.GetVision(vn,st,para).ravel()
    du,dw=mod.GetHydro(st,para);bu,bw=mod.GetHydroBEM(st,para)
    propulsion=np.vstack((np.cos(st[2]),np.sin(st[2])))
    vel=propulsion+du+bu
    collision=np.asarray(mod.Check_jumping(st,vel,para),bool)
    extreme=np.any(abs(vel)>500,axis=0);override=collision|extreme
    vel[:,override]=propulsion[:,override]
    dw=np.asarray(dw).ravel().copy();bw=np.asarray(bw).ravel().copy();dw[override]=0;bw[override]=0
    comp=dict(alignment=alignment,attraction=attraction,wall=mod.GetWallAvoid(st,para).ravel(),direct_hydro=dw,bem_hydro=bw)
    raw=sum(comp.values());eff=np.clip(raw,-np.pi/DT,np.pi/DT);comp['clip_residual']=eff-raw
    new=st[:2]+DT*vel;oldr=np.hypot(*st[:2])-mod.R;newr=np.hypot(*new)-mod.R
    sticky=(newr-mod.safedis>=0)&(oldr<newr);vel[:,sticky]=0
    vr=vel-vel.mean(axis=1)[:,None];xx,yy=q['pos'];r2=xx*xx+yy*yy
    if np.min(r2)<=0:raise ValueError('Undefined COM polar angle')
    omega=(xx*vr[1]-yy*vr[0])/r2
    b={s:float(q['grad']@v) for s,v in comp.items()}
    b['geometry']=float((-1j*np.mean((q['e']-q['p'])*q['f']*omega)).imag)
    b['diffusion']=-.5*IN**2*q['K'];b['total']=sum(b.values())
    nextpos=st[:2]+DT*vel;nextpos-=nextpos.mean(axis=1)[:,None]
    fn=np.exp(-1j*np.arctan2(nextpos[1],nextpos[0]))
    b['finite_total']=(np.exp(-.5*IN**2*DT)*np.mean(np.exp(1j*(st[2]+DT*eff))*(fn-fn.mean())).imag-q['K'])/DT
    return dict(q=q,W=W,J=J,a=a,vn=vn,vel=vel,comp=comp,b=b,
        clip=eff!=raw,override=override,sticky=sticky,
        visual_error=float(np.max(abs(alignment+attraction-visual))),weight_error=float(np.max(abs(align-alignment))))

def spectral_audit(st,vn):
    W,_,J,a=weights(st,vn);linear=IP*(W-np.eye(N))
    sync=st.copy();sync[2]=np.angle(np.mean(np.exp(1j*st[2])))
    Ws,_,Js,_=weights(sync,vn)
    # Independent centered numerical Jacobians of the nonlinear alignment law.
    errors=[]
    for z,jac in ((st,J),(sync,Js)):
        fd=np.empty((N,N));h=1e-6
        for i in range(N):
            plus=z.copy();minus=z.copy();plus[2,i]+=h;minus[2,i]-=h
            fd[:,i]=(weights(plus,vn)[1]-weights(minus,vn)[1])/(2*h)
        errors.append(float(np.max(abs(fd-jac))))
    return dict(row_sum_error=float(np.max(abs(W.sum(axis=1)-1))),
        frozen_W_max_real=float(np.max(np.linalg.eigvals(linear).real)),
        sync_max_real=float(np.max(np.linalg.eigvals(Js).real)),
        actual_J_max_real=float(np.max(np.linalg.eigvals(J).real)),
        actual_J_vs_WI=float(np.max(abs(J-linear))),
        sync_identity_error=float(np.max(abs(Js-IP*(Ws-np.eye(N))))),
        actual_J_fd_error=errors[0],sync_J_fd_error=errors[1])

def transport(replay,controls,out):
    rows=[];maps={};times=np.arange(0,5.001,.02);edges=np.linspace(-2.5,2.5,16)
    for seed,start,mirror in [(9,306.,1)]+[(int(c['seed']),float(c['t0']),int(c['mirror'])) for c in controls]:
        sm,_=replay.arrays(seed);initial=np.asarray(sm[:,:,round(start/DT)],float);qi=quantities(initial)
        front=qi['long']>=0;heat=[]
        for tau in times:
            st=np.asarray(sm[:,:,round((start+tau)/DT)],float);q=quantities(st)
            delta=mirror*q['delta'];long=q['long'];bins=np.digitize(long/np.std(long),edges)-1
            heat.append([float(np.mean(delta[bins==k])) if np.any(bins==k) else np.nan for k in range(len(edges)-1)])
            rows.append(dict(seed=seed,time=start+tau,tau=tau,mirror=mirror,K=mirror*q['K'],K_small=mirror*q['approx'],P=abs(q['p']),
                fixed_front_delta=float(delta[front].mean()),fixed_back_delta=float(delta[~front].mean()),
                rolling_front_delta=float(delta[long>=0].mean()),rolling_back_delta=float(delta[long<0].mean())))
        maps[f'seed{seed}']=np.array(heat)
    save(out/'transport_traces.csv',rows);np.savez_compressed(out/'transport_maps.npz',tau=times,edges=edges,**maps)
    # Fixed predeclared lags; correlation is descriptive, not a fitted speed.
    corr=[]
    for seed in [9]+list(range(7)):
        z=[r for r in rows if r['seed']==seed];f=np.array([r['fixed_front_delta'] for r in z]);b=np.array([r['fixed_back_delta'] for r in z])
        for lag in (-.5,-.2,0,.2,.5):
            k=round(lag/.02)
            ff,bb=(f[:-k],b[k:]) if k>0 else ((f[-k:],b[:k]) if k<0 else (f,b))
            corr.append(dict(seed=seed,lag=lag,front_to_later_back_corr=float(np.corrcoef(ff,bb)[0,1])))
    save(out/'transport_lag_checks.csv',corr)
    return rows,maps,times,edges

def summarize(rows,states,out):
    contrasts=[];rng=np.random.default_rng(11062026)
    for tau in (0,2,3,4,5):
        for source in (*SRC,'total','finite_total'):
            event=next(r['lambda_'+source] for r in rows if r['seed']==9 and r['tau']==tau and r['epsilon']==.001)
            vals=np.array([r['lambda_'+source] for r in rows if r['seed']!=9 and r['tau']==tau and r['epsilon']==.001])
            boots=rng.choice(vals,(10000,len(vals)),replace=True).mean(axis=1)
            contrasts.append(dict(tau=tau,event_time=306+tau,source=source,event=event,control_mean=vals.mean(),control_min=vals.min(),control_max=vals.max(),
                event_minus_control=event-vals.mean(),difference_ci_low=event-np.quantile(boots,.975),difference_ci_high=event-np.quantile(boots,.025)))
    save(out/'seed_control_contrasts.csv',contrasts)
    averages=[]
    for source in ('total','finite_total'):
        means={seed:np.mean([r['lambda_'+source] for r in rows if r['seed']==seed and r['epsilon']==.001]) for seed in [9]+list(range(7))}
        vals=np.array([means[s] for s in range(7)]);boots=rng.choice(vals,(10000,7),replace=True).mean(axis=1)
        averages.append(dict(source=source,event_five_snapshot_mean=means[9],control_mean=vals.mean(),
            difference=means[9]-vals.mean(),difference_ci_low=means[9]-np.quantile(boots,.975),difference_ci_high=means[9]-np.quantile(boots,.025)))
    save(out/'seed_averaged_contrasts.csv',averages)
    checks=[]
    for st in states:
        r=[x for x in rows if x['seed']==st['seed'] and x['tau']==st['tau']]
        ref=next(x for x in r if x['epsilon']==.0005)
        checks.append(dict(seed=st['seed'],tau=st['tau'],max_lambda_total_difference=max(abs(x['lambda_total']-ref['lambda_total']) for x in r),
            max_lambda_source_difference=max(abs(x['lambda_'+s]-ref['lambda_'+s]) for x in r for s in SRC),
            max_p_change=max(x['p_even_change'] for x in r),max_graph_changes=max(x['graph_changes'] for x in r),
            max_clip_changes=max(x['clip_changes'] for x in r),max_override_changes=max(x['override_changes'] for x in r)))
    save(out/'convergence.csv',checks)
    return contrasts

def plot(rows,states,contrasts,traces,maps,times,edges,out):
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for seed in [*range(7),9]:
        a=[r for r in states if r['seed']==seed];color='crimson' if seed==9 else 'gray'
        ax[0].plot([r['time']-r['tau']+r['tau'] if seed==9 else 306+r['tau'] for r in a],[r['mirror']*r['K'] for r in a],color=color,alpha=1 if seed==9 else .35)
        ax[0].scatter([306+r['tau'] for r in a],[r['mirror']*r['K_small'] for r in a],color=color,marker='x')
    ax[0].set(xlabel='event-equivalent time',ylabel='K (line), small-heading approximation (×)')
    for s in ('total','finite_total'):
        a=[r for r in contrasts if r['source']==s];x=np.array([r['event_time'] for r in a])
        ax[1].plot(x,[r['event'] for r in a],'o-',label='event '+s)
    a=[r for r in contrasts if r['source']=='total'];ax[1].plot(x,[r['control_mean'] for r in a],color='gray',label='control mean')
    ax[1].fill_between(x,[r['control_min'] for r in a],[r['control_max'] for r in a],color='gray',alpha=.2,label='control seed range')
    ax[1].axhline(0,color='k',lw=.6);ax[1].set(xlabel='event time',ylabel='directional response per change in K');ax[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'local_response.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4));a=[r for r in rows if r['seed']==9 and r['epsilon']==.001]
    for s in ('alignment','direct_hydro','wall','attraction','geometry','diffusion'):
        ax[0].plot([306+r['tau'] for r in a],[r['lambda_'+s] for r in a],'o-',label=s)
    ax[0].axhline(0,color='k',lw=.6);ax[0].set(xlabel='event time',ylabel='source directional response');ax[0].legend(fontsize=8)
    for seed in [*range(7),9]:
        z=[r for r in states if r['seed']==seed];ax[1].plot([306+r['tau'] for r in z],[r['mean_a_forward'] for r in z],color='crimson' if seed==9 else 'gray',alpha=1 if seed==9 else .35)
    ax[1].axhline(0,color='k',lw=.6);ax[1].set(xlabel='event-equivalent time',ylabel='mean neighbor displacement along heading')
    fig.tight_layout();fig.savefig(out/'sources_and_transport_bias.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(11,7))
    for axis,seed in zip(ax[0],(9,6)):
        im=axis.pcolormesh(306+times,(edges[1:]+edges[:-1])/2,maps[f'seed{seed}'].T,cmap='coolwarm',vmin=-.8,vmax=.8,shading='nearest')
        axis.set(title=f'seed {seed}',ylabel='longitudinal position / std; front > 0');fig.colorbar(im,ax=axis,label='heading deviation [rad]')
    for axis,seed in zip(ax[1],(9,6)):
        a=[r for r in traces if r['seed']==seed]
        for k in ('fixed_front_delta','fixed_back_delta'):axis.plot([306+r['tau'] for r in a],[r[k] for r in a],label=k)
        axis.legend(fontsize=8);axis.set(xlabel='event-equivalent time',ylabel='mean heading deviation [rad]')
    fig.tight_layout();fig.savefig(out/'fore_aft_transport.png',dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);args=ap.parse_args()
    out=HERE/'outputs';out.mkdir(parents=True,exist_ok=True)
    source=BASE.parents[1].parent/'BEM_circular_parallel_instance.py'
    sha=hashlib.sha256(source.read_bytes()).hexdigest();assert sha==json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256']
    replay=COMReplay(args.data_dir,source,BASE/'outputs')
    controls=list(csv.DictReader((HERE.parent/'outputs/matched_controls.csv').open()))
    cases=[(9,306.,1)]+[(int(c['seed']),float(c['t0']),int(c['mirror'])) for c in controls]
    rows=[];states=[];saved={};maxvel=maxcomp=0.
    for seed,start,mirror in cases:
        sm,rm=replay.arrays(seed)
        for tau in (0,2,3,4,5):
            step=round((start+tau)/DT);st=np.asarray(sm[:,:,step],float);base=evaluate(replay,st);q=base['q']
            original,detail=replay.evaluate(seed,step)
            maxvel=max(maxvel,float(np.max(abs(base['vel']-rm[:,:,step+1]))))
            maxcomp=max(maxcomp,max(float(np.max(abs(base['comp'][k]-detail['components'][k]))) for k in detail['components']))
            constraints=np.vstack((-np.sin(st[2])/N,np.cos(st[2])/N))
            u,sv,vh=np.linalg.svd(constraints,full_matrices=False);rank=int(np.sum(sv>sv[0]*1e-12))
            h=q['grad']-vh[:rank].T@(vh[:rank]@q['grad']);norm=np.sqrt(np.mean(h*h))
            if norm<1e-12:raise ValueError('Degenerate constrained direction')
            h/=norm;projection=q['grad']@h
            spec=spectral_audit(st,base['vn'])
            forward=np.array([np.cos(q['Theta']),np.sin(q['Theta'])])@base['a']
            rec=dict(seed=seed,time=start+tau,tau=tau,mirror=mirror,K=q['K'],K_small=q['approx'],P=abs(q['p']),
                max_abs_delta=float(np.max(abs(q['delta']))),rms_delta=float(np.sqrt(np.mean(q['delta']**2))),
                constraint_error=float(np.max(abs(constraints@h))),dK_dh=projection,rank=rank,
                mean_a_forward=float(forward.mean()),fraction_a_forward_positive=float(np.mean(forward>0)),
                characteristic_forward=-IP*float(forward.mean()),base_clip=int(base['clip'].sum()),base_override=int(base['override'].sum()),
                visual_error=base['visual_error'],weight_error=base['weight_error'],**spec)
            states.append(rec)
            name=f's{seed}_t{step}';saved[name+'_W']=base['W'];saved[name+'_a']=base['a'];saved[name+'_h']=h;saved[name+'_delta']=q['delta'];saved[name+'_long']=q['long']
            for eps in (.0005,.001,.002):
                trial=[]
                for sign in (-1,1):
                    stp=st.copy();stp[2]+=sign*eps*h;trial.append(evaluate(replay,stp))
                minus,plus=trial;dk=plus['q']['K']-minus['q']['K']
                r=dict(seed=seed,time=start+tau,tau=tau,epsilon=eps,K_difference=dk,K_direction_derivative=dk/(2*eps),
                    p_odd_derivative=abs(plus['q']['p']-minus['q']['p'])/(2*eps),
                    p_even_change=max(abs(x['q']['p']-q['p']) for x in trial),
                    graph_changes=sum(int(np.count_nonzero(x['vn']!=base['vn'])) for x in trial),
                    clip_changes=sum(int(np.count_nonzero(x['clip']!=base['clip'])) for x in trial),
                    override_changes=sum(int(np.count_nonzero(x['override']!=base['override'])+np.count_nonzero(x['sticky']!=base['sticky'])) for x in trial),
                    max_clip_fish=max(int(x['clip'].sum()) for x in trial),
                    weight_change=max(float(np.max(abs(x['W']-base['W']))) for x in trial))
                for s in (*SRC,'total','finite_total'):r['lambda_'+s]=(plus['b'][s]-minus['b'][s])/dk
                for s in base['comp']:
                    r['frozen_projection_'+s]=float(q['grad']@(plus['comp'][s]-minus['comp'][s]))/(2*eps*projection)
                rows.append(r)
            print(f'seed {seed}, t={start+tau}: response {rows[-2]["lambda_total"]:.5f}',flush=True)
    save(out/'states.csv',states);save(out/'responses.csv',rows);np.savez_compressed(out/'weights_directions.npz',**saved)
    traces,maps,times,edges=transport(replay,controls,out)
    contrasts=summarize(rows,states,out);plot(rows,states,contrasts,traces,maps,times,edges,out)
    (out/'audit.json').write_text(json.dumps(dict(production_sha256=sha,max_velocity_error=maxvel,max_component_error=maxcomp,
        snapshots=len(states),trial_rhs_calls=6*len(states),epsilons=[.0005,.001,.002],controls='Inherited seven seed-level matched windows; no rematching',
        data_dir=str(args.data_dir),dt=DT,In=IN,Ip=IP,bootstrap_seed=11062026,bootstrap_draws=10000),indent=2)+'\n')

if __name__=='__main__':main()
