#!/usr/bin/env python3
"""Falsifiable closure test for z=(m,P_L,P_R), with whole-seed validation."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DT = .01
TRAIN = tuple(range(15))
TEST = tuple(range(15,20))
LAGS = (.1,.5,1.0)
MIXED_EDGES = np.array([-1.000001,-.75,-.4,-.2,0,.2,.4,.75,1.000001])
DOM = .4
NBOOT = 1000
MIN_TRAIN_ROW = 500
MIN_TEST_ROW = 200
RNG = np.random.default_rng(20261029)


def write_csv(path, rows):
    if not rows: return
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def boot_mean(values, nboot=NBOOT):
    v=np.asarray(values,float);v=v[np.isfinite(v)]
    if not len(v): return (np.nan,np.nan,np.nan,0)
    b=np.mean(v[RNG.integers(0,len(v),size=(nboot,len(v)))],axis=1)
    return float(np.mean(v)),float(np.quantile(b,.025)),float(np.quantile(b,.975)),len(v)


def row_norm(c):
    rows=c.sum(1);p=np.divide(c,rows[:,None],out=np.zeros_like(c),where=rows[:,None]>0)
    p[rows==0,rows==0]=1
    return p,rows


def js_div(p,q):
    p=np.asarray(p,float);q=np.asarray(q,float)
    p=p/p.sum();q=q/q.sum();mid=.5*(p+q)
    def kl(a,b):
        z=a>0;return np.sum(a[z]*np.log(a[z]/b[z]))
    return float(.5*kl(p,mid)+.5*kl(q,mid))


class Discretizer:
    def __init__(self,data):
        m=np.concatenate([data["m"][s] for s in TRAIN])
        pl=np.concatenate([data["P_L"][s] for s in TRAIN])
        pr=np.concatenate([data["P_R"][s] for s in TRAIN])
        self.bound_left=np.quantile(pl[m==-1],[1/3,2/3])
        self.bound_right=np.quantile(pr[m==1],[1/3,2/3])
        self.pl_med=np.zeros(8);self.pr_med=np.zeros(8)
        mixed=(m>-1)&(m<1);mb=np.clip(np.searchsorted(MIXED_EDGES,m,side="right")-1,0,7)
        for k in range(8):
            z=mixed&(mb==k)
            self.pl_med[k]=np.median(pl[z]);self.pr_med[k]=np.median(pr[z])
        # 3 left-boundary + 8*4 mixed + 3 right-boundary = 38 states.
        self.nstate=38;self.ngroup=10
        self.state_group=np.empty(38,dtype=int)
        self.state_group[:3]=0
        for k in range(8):self.state_group[3+4*k:3+4*k+4]=k+1
        self.state_group[35:]=9
        self.group_center=np.r_[-1,.5*(MIXED_EDGES[:-1]+MIXED_EDGES[1:]),1]

    def state(self,m,pl,pr):
        m=np.asarray(m);out=np.empty(m.shape,dtype=np.int16)
        left=m==-1;right=m==1;mixed=~(left|right)
        out[left]=np.searchsorted(self.bound_left,pl[left],side="right")
        out[right]=35+np.searchsorted(self.bound_right,pr[right],side="right")
        mb=np.clip(np.searchsorted(MIXED_EDGES,m[mixed],side="right")-1,0,7)
        ph=(pl[mixed]>self.pl_med[mb]).astype(int);rh=(pr[mixed]>self.pr_med[mb]).astype(int)
        out[mixed]=3+4*mb+2*ph+rh
        return out

    def mstate(self,m):
        m=np.asarray(m);out=np.empty(m.shape,dtype=np.int16)
        left=m==-1;right=m==1;mixed=~(left|right)
        out[left]=0;out[right]=9
        out[mixed]=1+np.clip(np.searchsorted(MIXED_EDGES,m[mixed],side="right")-1,0,7)
        return out

    def json(self):
        return {"mixed_m_edges":MIXED_EDGES.tolist(),"boundary_left_P_tertiles":self.bound_left.tolist(),
                "boundary_right_P_tertiles":self.bound_right.tolist(),"mixed_PL_medians":self.pl_med.tolist(),
                "mixed_PR_medians":self.pr_med.tolist(),"n_z_states":38,"n_m_states":10}


def load_data(path):
    with np.load(path) as z:
        data={k:np.asarray(z[k]) for k in z.files}
    if data["m"].shape!=(20,80001):raise ValueError(data["m"].shape)
    return data


def state_arrays(data,disc,model):
    result={}
    for s in range(20):
        result[s]=disc.state(data["m"][s],data["P_L"][s],data["P_R"][s]) if model=="mPLPR" else disc.mstate(data["m"][s])
    return result


def counts_by_seed(states,seeds,lag,nstate):
    out={}
    for s in seeds:
        x=states[s]
        out[s]=np.bincount(x[:-lag]*nstate+x[lag:],minlength=nstate*nstate).reshape(nstate,nstate).astype(float)
    return out


def collapse(matrix,mapping,ngroup):
    out=np.zeros((matrix.shape[0],ngroup))
    for k in range(ngroup):out[:,k]=matrix[:,mapping==k].sum(1)
    return out


def ck_metric(c1,c2,mapping,centers,scope,minrow=MIN_TRAIN_ROW):
    p1,r1=row_norm(c1);p2,r2=row_norm(c2);composed=p1@p1
    proj2=collapse(p2,mapping,len(centers));projc=collapse(composed,mapping,len(centers))
    tvp=.5*np.abs(proj2-projc).sum(1);tvf=.5*np.abs(p2-composed).sum(1)
    valid=(r1>=minrow)&(r2>=minrow)
    if scope=="corridor":valid &= np.abs(centers[mapping])<DOM
    w=r2[valid];w=w/w.sum()
    return float(w@tvp[valid]),float(w@tvf[valid]),int(valid.sum()),int(r2[valid].sum())


def kernel_floor(train,test,mapping,centers,scope):
    pt,rt=row_norm(train);pv,rv=row_norm(test)
    a=collapse(pt,mapping,len(centers));b=collapse(pv,mapping,len(centers))
    tvp=.5*np.abs(a-b).sum(1);tvf=.5*np.abs(pt-pv).sum(1)
    valid=(rt>=MIN_TRAIN_ROW)&(rv>=MIN_TEST_ROW)
    if scope=="corridor":valid &= np.abs(centers[mapping])<DOM
    w=rv[valid];w=w/w.sum()
    return float(w@tvp[valid]),float(w@tvf[valid]),int(valid.sum())


def ck_analysis(states_by_model,disc):
    rows=[];kernels={}
    for tau in LAGS:
        lag=int(round(tau/DT))
        for model,states in states_by_model.items():
            n=38 if model=="mPLPR" else 10
            mapping=disc.state_group if model=="mPLPR" else np.arange(10)
            one=counts_by_seed(states,TRAIN,lag,n);two=counts_by_seed(states,TRAIN,2*lag,n)
            test=counts_by_seed(states,TEST,2*lag,n);kernels[(model,tau)]=row_norm(sum(one.values()))[0]
            for scope in ("all","corridor"):
                point=ck_metric(sum(one.values()),sum(two.values()),mapping,disc.group_center,scope)
                floor=kernel_floor(sum(two.values()),sum(test.values()),mapping,disc.group_center,scope)
                bp=[];bf=[]
                seeds=np.asarray(TRAIN)
                for _ in range(300):
                    draw=RNG.choice(seeds,len(seeds),replace=True)
                    val=ck_metric(sum(one[int(s)] for s in draw),sum(two[int(s)] for s in draw),mapping,disc.group_center,scope)
                    bp.append(val[0]);bf.append(val[1])
                rows.append(dict(model=model,tau=tau,scope=scope,projected_m_ck_tv=point[0],
                    projected_ci_low=np.quantile(bp,.025),projected_ci_high=np.quantile(bp,.975),
                    full_state_ck_tv=point[1],full_ci_low=np.quantile(bf,.025),full_ci_high=np.quantile(bf,.975),
                    supported_states=point[2],supported_samples=point[3],train_test_projected_tv=floor[0],
                    train_test_full_tv=floor[1],train_test_supported_states=floor[2]))
    return rows,kernels


def corridor_samples(data,states):
    rows=[];censored=0;durations=[]
    for seed in range(20):
        m=data["m"][seed];start=None;buf=[]
        for i in range(1,len(m)):
            if start is None and abs(m[i])<=.2 and abs(m[i-1])>.2:
                start=i;buf=[]
            if start is not None:
                if (i-start)%10==0 and abs(m[i])<=.2:
                    buf.append(i)
                if m[i]>=DOM or m[i]<=-DOM:
                    right=int(m[i]>=DOM);durations.append(dict(seed=seed,duration=(i-start)*DT,right=right))
                    for j in buf:
                        rows.append(dict(seed=seed,index=j,state=int(states[seed][j]),right=right,
                            dm=float(data["dm_past"][seed,j]),heading=float(data["net_heading_right"][seed,j])))
                    start=None;buf=[]
        if start is not None:censored+=1
    return rows,censored,durations


def committor(kernel,mapping,centers):
    left=centers[mapping]<=-DOM;right=centers[mapping]>=DOM;inner=~(left|right)
    q=np.zeros(len(mapping));q[right]=1
    A=np.eye(inner.sum())-kernel[np.ix_(inner,inner)];b=kernel[np.ix_(inner,right)].sum(1)
    try:q[inner]=np.linalg.solve(A,b)
    except np.linalg.LinAlgError:q[inner]=np.linalg.lstsq(A,b,rcond=None)[0]
    return np.clip(q,0,1)


def seed_scores(rows,prediction,name,common=None):
    out=[]
    for seed in TEST:
        idx=[i for i,r in enumerate(rows) if r["seed"]==seed and (common is None or common[i])]
        if not idx:continue
        y=np.array([rows[i]["right"] for i in idx]);p=np.clip(prediction[idx],1e-6,1-1e-6)
        out.append(dict(model=name,seed=seed,n=len(idx),brier=float(np.mean((y-p)**2)),
                        log_loss=float(-np.mean(y*np.log(p)+(1-y)*np.log(1-p)))))
    return out


def predictive_analysis(data,disc,states_by_model,kernels):
    all_samples={};duration_rows=None;censored=None
    for model in ("m","mPLPR"):
        rows,cens,durs=corridor_samples(data,states_by_model[model]);all_samples[model]=rows
        if model=="mPLPR":duration_rows=durs;censored=cens
    score_rows=[]
    # Operator committor uses exactly the finite-state transition kernel.
    for tau in LAGS:
        for model in ("m","mPLPR"):
            mapping=disc.state_group if model=="mPLPR" else np.arange(10)
            q=committor(kernels[(model,tau)],mapping,disc.group_center)
            rows=all_samples[model];pred=np.array([q[r["state"]] for r in rows])
            score_rows += [dict(kind="operator_committor",tau=tau,**r) for r in seed_scores(rows,pred,model)]

    # Empirical state-matched diagnostics: same held-out observations and no future predictor leakage.
    train=[r for r in all_samples["mPLPR"] if r["seed"] in TRAIN]
    test=[r for r in all_samples["mPLPR"] if r["seed"] in TEST]
    n=38
    def fit_index(rows,kind):
        st=np.array([r["state"] for r in rows])
        if kind=="state":return st,n
        if kind=="past":
            cat=np.array([0 if r["dm"]<0 else 2 if r["dm"]>0 else 1 for r in rows]);return 3*st+cat,3*n
        cat=np.array([int(r["heading"]>=0) for r in rows]);return 2*st+cat,2*n
    fits={};test_index={};supports=[]
    ytr=np.array([r["right"] for r in train]);yte=np.array([r["right"] for r in test])
    for kind in ("state","past","heading"):
        ix,nx=fit_index(train,kind);cnt=np.bincount(ix,minlength=nx);succ=np.bincount(ix,weights=ytr,minlength=nx)
        fits[kind]=(succ+1)/(cnt+2);ti,_=fit_index(test,kind);test_index[kind]=ti;supports.append(cnt[ti]>=50)
    common=np.logical_and.reduce(supports)
    for kind,label in (("state","state_empirical"),("past","state_plus_past_dm"),("heading","state_plus_signed_heading")):
        pred=fits[kind][test_index[kind]]
        score_rows += [dict(kind="memory_diagnostic",tau="",**r) for r in seed_scores(test,pred,label,common)]

    # Supported within-state probability splits, retained as direct evidence of hidden direction/history.
    contrast=[]
    for state in range(38):
        for variable in ("dm","heading"):
            neg=[r for r in train if r["state"]==state and r[variable]<0]
            pos=[r for r in train if r["state"]==state and r[variable]>0]
            if min(len(neg),len(pos))<100:continue
            contrast.append(dict(state=state,variable=variable,n_negative=len(neg),n_positive=len(pos),
                p_right_negative=np.mean([r["right"] for r in neg]),p_right_positive=np.mean([r["right"] for r in pos]),
                positive_minus_negative=np.mean([r["right"] for r in pos])-np.mean([r["right"] for r in neg])))
    return score_rows,contrast,all_samples["mPLPR"],duration_rows,censored


def km(episodes,times):
    d=np.array([r["duration"] for r in episodes]);e=np.array([r["event"] for r in episodes])
    surv=np.ones(len(times));value=1.
    for k,t in enumerate(times):
        lo=times[k-1] if k else -1
        for et in np.unique(d[(d>lo)&(d<=t)&(e==1)]):
            risk=np.sum(d>=et);death=np.sum((d==et)&(e==1));value*=1-death/risk
        surv[k]=value
    return surv


def dominance_episodes(m,states,seed):
    regime=0;start=None;out=[]
    for i,v in enumerate(m):
        new=1 if v>=DOM else -1 if v<=-DOM else regime
        if regime==0 and new:regime,start=new,i
        elif regime and new!=regime:
            out.append(dict(seed=seed,direction="LR" if regime==-1 else "RL",duration=(i-start)*DT,event=1,state=int(states[start])))
            regime,start=new,i
    if regime:
        out.append(dict(seed=seed,direction="LR" if regime==-1 else "RL",duration=(len(m)-1-start)*DT,event=0,state=int(states[start])))
    return out


def model_survival(kernel,mapping,centers,starts,direction,times,tau):
    target=centers[mapping]>=DOM if direction=="LR" else centers[mapping]<=-DOM
    transient=np.flatnonzero(~target);lookup={s:i for i,s in enumerate(transient)}
    ss=[lookup[s] for s in starts if s in lookup];sub=kernel[np.ix_(transient,transient)]
    v=np.bincount(ss,minlength=len(transient)).astype(float);v/=v.sum();out=[];step=0
    for t in times:
        want=int(round(t/tau))
        while step<want:v=v@sub;step+=1
        out.append(v.sum())
    return np.asarray(out)


def fpt_analysis(data,disc,states_by_model,kernels):
    all_ep={}
    for model in ("m","mPLPR"):
        all_ep[model]=sum((dominance_episodes(data["m"][s],states_by_model[model][s],s) for s in range(20)),[])
    rows=[];curves=[]
    for tau in (.5,1.0):
        times=np.arange(0,180.0001,tau)
        for direction in ("LR","RL"):
            obs=[r for r in all_ep["m"] if r["seed"] in TEST and r["direction"]==direction]
            observed=km(obs,times);obs_rmst=float(np.trapezoid(observed,times))
            curves += [dict(tau=tau,direction=direction,time=t,series="observed",survival=v) for t,v in zip(times,observed)]
            for model in ("m","mPLPR"):
                mapping=disc.state_group if model=="mPLPR" else np.arange(10)
                ep=[r for r in all_ep[model] if r["seed"] in TEST and r["direction"]==direction]
                pred=model_survival(kernels[(model,tau)],mapping,disc.group_center,[r["state"] for r in ep],direction,times,tau)
                pred_rmst=float(np.trapezoid(pred,times));curves += [dict(tau=tau,direction=direction,time=t,series=model,survival=v) for t,v in zip(times,pred)]
                seed_diff=[];seed_max=[]
                for s in TEST:
                    oe=[r for r in all_ep["m"] if r["seed"]==s and r["direction"]==direction]
                    me=[r for r in all_ep[model] if r["seed"]==s and r["direction"]==direction]
                    if not oe or not me:continue
                    os=km(oe,times);ps=model_survival(kernels[(model,tau)],mapping,disc.group_center,[r["state"] for r in me],direction,times,tau)
                    seed_diff.append(np.trapezoid(ps-os,times));seed_max.append(np.max(np.abs(ps-os)))
                md=boot_mean(seed_diff);mm=boot_mean(seed_max)
                rows.append(dict(tau=tau,direction=direction,model=model,heldout_episodes=len(ep),events=sum(r["event"] for r in ep),
                    censored=sum(1-r["event"] for r in ep),observed_rmst=obs_rmst,predicted_rmst=pred_rmst,
                    pooled_survival_max_error=np.max(np.abs(pred-observed)),seed_mean_rmst_error=md[0],rmst_error_ci_low=md[1],
                    rmst_error_ci_high=md[2],seed_mean_max_error=mm[0],max_error_ci_low=mm[1],max_error_ci_high=mm[2]))
    return rows,curves


def nucleation_analysis(data,disc,states):
    horizon=50;stride=10;records=[]
    for seed in range(20):
        m=data["m"][seed]
        for i in range(0,len(m)-horizon,stride):
            if abs(m[i])!=1:continue
            left=m[i]==-1;h=float(data["heading_to_gap_L"][seed,i] if left else data["heading_to_gap_R"][seed,i])
            outcome=int(np.any(m[i+1:i+horizon+1]!=m[i]))
            records.append(dict(seed=seed,state=int(states[seed][i]),heading=h,outcome=outcome,right=outcome))
    train=[r for r in records if r["seed"] in TRAIN];test=[r for r in records if r["seed"] in TEST]
    y=np.array([r["outcome"] for r in train]);st=np.array([r["state"] for r in train]);hs=np.array([int(r["heading"]>=0) for r in train])
    def fit(ix,n):
        c=np.bincount(ix,minlength=n);s=np.bincount(ix,weights=y,minlength=n);return (s+1)/(c+2),c
    p0,c0=fit(st,38);p1,c1=fit(2*st+hs,76)
    tst=np.array([r["state"] for r in test]);ths=np.array([int(r["heading"]>=0) for r in test]);common=(c0[tst]>=100)&(c1[2*tst+ths]>=100)
    rows=[]
    for name,pred in (("occupied_polarization_state",p0[tst]),("plus_heading_to_gap_sign",p1[2*tst+ths])):
        rows += seed_scores(test,pred,name,common)
    return rows,records


def nucleation_rates(records,data):
    rows=[]
    for split,seeds in (("train",TRAIN),("heldout",TEST)):
        q=[r for r in records if r["seed"] in seeds]
        for heading_class,keep in (("all",lambda r:True),("toward_gap",lambda r:r["heading"]>=0),("away_from_gap",lambda r:r["heading"]<0)):
            z=[r for r in q if keep(r)]
            rows.append(dict(split=split,heading_class=heading_class,n_frames=len(z),
                exit_within_0p5=np.mean([r["outcome"] for r in z])))
    return rows


def occupancy_support(data,out):
    rows=[];slices=((-0.75,-.2),(-.2,.2),(.2,.75))
    for lo,hi in slices:
        for bins in (8,12,16):
            hist=[]
            for seeds in (TRAIN,TEST):
                m=np.concatenate([data["m"][s] for s in seeds]);pl=np.concatenate([data["P_L"][s] for s in seeds]);pr=np.concatenate([data["P_R"][s] for s in seeds])
                z=(m>=lo)&(m<hi);hist.append(np.histogram2d(pl[z],pr[z],bins=bins,range=((0,1),(0,1)))[0])
            rows.append(dict(m_low=lo,m_high=hi,bins=bins,train_samples=int(hist[0].sum()),test_samples=int(hist[1].sum()),train_test_js=js_div(hist[0].ravel()+.5,hist[1].ravel()+.5)))
    return rows


def figures(data,ck,score,fpt,curves,nucleation,out):
    out.mkdir(parents=True,exist_ok=True)
    slices=((-0.75,-.2),(-.2,.2),(.2,.75));fig,axs=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
    for row,seeds in enumerate((TRAIN,TEST)):
        m=np.concatenate([data["m"][s] for s in seeds]);pl=np.concatenate([data["P_L"][s] for s in seeds]);pr=np.concatenate([data["P_R"][s] for s in seeds])
        for col,(lo,hi) in enumerate(slices):
            z=(m>=lo)&(m<hi);h=np.histogram2d(pl[z],pr[z],bins=12,range=((0,1),(0,1)))[0];threshold=100 if row==0 else 30
            u=np.full(h.shape,np.nan);ok=h>=threshold;u[ok]=-np.log(h[ok]/h.sum());u[ok]-=np.nanmin(u)
            im=axs[row,col].imshow(u.T,origin="lower",extent=(0,1,0,1),aspect="auto",cmap="viridis_r")
            axs[row,col].set(title=f"{['train','held out'][row]}: {lo}<m<{hi}\nn={int(h.sum())}",xlabel="$P_L$",ylabel="$P_R$")
            fig.colorbar(im,ax=axs[row,col],label="$U_{occ}$")
        for side,color in ((-1,"#0072B2"),(1,"#D55E00")):
            z=m==side;p=pl[z] if side==-1 else pr[z]
            axs[row,3].hist(p,bins=np.linspace(0,1,31),density=True,histtype="step",lw=2,color=color,label=f"m={side:+d}, n={len(p)}")
        axs[row,3].set(title=f"{['train','held out'][row]} boundary strata",xlabel="occupied-group polarization",ylabel="density");axs[row,3].legend(fontsize=8)
    fig.savefig(out/"01_occupancy_landscape.png",dpi=210);fig.savefig(out/"01_occupancy_landscape.pdf");plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for model,color in (("m","#777777"),("mPLPR","#0072B2")):
        z=[r for r in ck if r["model"]==model and r["scope"]=="corridor"]
        axs[0].plot([r["tau"] for r in z],[r["projected_m_ck_tv"] for r in z],"o-",label=model,color=color)
        axs[0].plot([r["tau"] for r in z],[r["train_test_projected_tv"] for r in z],"--",color=color,alpha=.6)
    axs[0].set(xlabel="lag",ylabel="projected-m CK total variation",title="Solid: CK error; dashed: train-test floor");axs[0].legend();axs[0].grid(alpha=.2)
    mem=[r for r in score if r["kind"]=="memory_diagnostic"]
    names=["state_empirical","state_plus_past_dm","state_plus_signed_heading"]
    means=[];lo=[];hi=[]
    for name in names:
        vals=[r["brier"] for r in mem if r["model"]==name];b=boot_mean(vals);means.append(b[0]);lo.append(b[0]-b[1]);hi.append(b[2]-b[0])
    axs[1].bar(range(3),means,yerr=np.array([lo,hi]),color=["#777777","#E69F00","#009E73"],capsize=4)
    axs[1].set(xticks=range(3),xticklabels=["state",r"+ past $\Delta m$","+ signed heading"],ylabel="held-out Brier",title="Same supported held-out samples")
    fig.savefig(out/"02_closure_and_memory.png",dpi=210);fig.savefig(out/"02_closure_and_memory.pdf");plt.close(fig)

    fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for col,direction in enumerate(("LR","RL")):
        z=[r for r in curves if r["tau"]==.5 and r["direction"]==direction]
        for series,color,ls in (("observed","k","-"),("m","#777777","--"),("mPLPR","#0072B2","-")):
            q=[r for r in z if r["series"]==series];axs[0,col].plot([r["time"] for r in q],[r["survival"] for r in q],label=series,color=color,ls=ls)
        axs[0,col].set(xlim=(0,120),ylim=(0,1.02),xlabel="time",ylabel="survival",title=f"held-out {direction}, kernel lag 0.5");axs[0,col].legend();axs[0,col].grid(alpha=.2)
    nrows=nucleation;names=["occupied_polarization_state","plus_heading_to_gap_sign"]
    vals=[];err=[]
    for name in names:
        q=[r["brier"] for r in nrows if r["model"]==name];b=boot_mean(q);vals.append(b[0]);err.append([b[0]-b[1],b[2]-b[0]])
    axs[1,0].bar(range(2),vals,yerr=np.array(err).T,color=["#777777","#009E73"],capsize=4);axs[1,0].set(xticks=range(2),xticklabels=["occupied $P$","+ heading sign"],ylabel="held-out Brier",title="Boundary nucleation within 0.5")
    op=[r for r in score if r["kind"]=="operator_committor" and r["tau"]==.5]
    vals=[];err=[]
    for name in ("m","mPLPR"):
        q=[r["brier"] for r in op if r["model"]==name];b=boot_mean(q);vals.append(b[0]);err.append([b[0]-b[1],b[2]-b[0]])
    axs[1,1].bar(range(2),vals,yerr=np.array(err).T,color=["#777777","#0072B2"],capsize=4);axs[1,1].set(xticks=range(2),xticklabels=["m","$(m,P_L,P_R)$"],ylabel="held-out Brier",title="Operator committor, lag 0.5")
    fig.savefig(out/"03_heldout_passage.png",dpi=210);fig.savefig(out/"03_heldout_passage.pdf");plt.close(fig)


def summarize(score,ck,fpt,nucleation,contrast):
    rows=[]
    for kind in ("operator_committor","memory_diagnostic"):
        for model in sorted({r["model"] for r in score if r["kind"]==kind}):
            taus=sorted({r["tau"] for r in score if r["kind"]==kind and r["model"]==model},key=lambda x:str(x))
            for tau in taus:
                q=[r for r in score if r["kind"]==kind and r["model"]==model and r["tau"]==tau]
                for metric in ("brier","log_loss"):
                    b=boot_mean([r[metric] for r in q]);rows.append(dict(section=kind,model=model,tau=tau,metric=metric,mean=b[0],ci_low=b[1],ci_high=b[2],n_seeds=b[3],n_samples=sum(r["n"] for r in q)))
    for model in sorted({r["model"] for r in nucleation}):
        q=[r for r in nucleation if r["model"]==model]
        for metric in ("brier","log_loss"):
            b=boot_mean([r[metric] for r in q]);rows.append(dict(section="nucleation",model=model,tau=.5,metric=metric,mean=b[0],ci_low=b[1],ci_high=b[2],n_seeds=b[3],n_samples=sum(r["n"] for r in q)))
    return rows


def paired_differences(score,nucleation):
    rows=[]
    comparisons=[]
    for tau in LAGS:
        comparisons.append(("operator_committor",tau,"mPLPR","m"))
    comparisons += [("memory_diagnostic","","state_plus_past_dm","state_empirical"),
                    ("memory_diagnostic","","state_plus_signed_heading","state_empirical")]
    for kind,tau,a,b in comparisons:
        for metric in ("brier","log_loss"):
            aa={r["seed"]:r[metric] for r in score if r["kind"]==kind and r["tau"]==tau and r["model"]==a}
            bb={r["seed"]:r[metric] for r in score if r["kind"]==kind and r["tau"]==tau and r["model"]==b}
            seeds=sorted(set(aa)&set(bb));v=[aa[s]-bb[s] for s in seeds];q=boot_mean(v)
            rows.append(dict(section=kind,tau=tau,model_a=a,model_b=b,metric=metric,
                difference_a_minus_b=q[0],ci_low=q[1],ci_high=q[2],n_seeds=q[3]))
    for metric in ("brier","log_loss"):
        aa={r["seed"]:r[metric] for r in nucleation if r["model"]=="plus_heading_to_gap_sign"}
        bb={r["seed"]:r[metric] for r in nucleation if r["model"]=="occupied_polarization_state"}
        seeds=sorted(set(aa)&set(bb));q=boot_mean([aa[s]-bb[s] for s in seeds])
        rows.append(dict(section="nucleation",tau=.5,model_a="plus_heading_to_gap_sign",
            model_b="occupied_polarization_state",metric=metric,difference_a_minus_b=q[0],ci_low=q[1],ci_high=q[2],n_seeds=q[3]))
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument("--features",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--qj-results",type=Path,required=True)
    args=p.parse_args();out=args.output_dir;tables=out/"tables";figdir=out/"figures";tables.mkdir(parents=True,exist_ok=True)
    data=load_data(args.features);disc=Discretizer(data)
    states={model:state_arrays(data,disc,model) for model in ("m","mPLPR")}
    ck,kernels=ck_analysis(states,disc)
    score,contrast,samples,durations,censored=predictive_analysis(data,disc,states,kernels)
    fpt,curves=fpt_analysis(data,disc,states,kernels)
    nucleation,nuc_records=nucleation_analysis(data,disc,states["mPLPR"])
    support=occupancy_support(data,out)
    summary=summarize(score,ck,fpt,nucleation,contrast)
    write_csv(tables/"markov_consistency.csv",ck);write_csv(tables/"predictive_scores_by_seed.csv",score)
    write_csv(tables/"direction_history_contrasts.csv",contrast);write_csv(tables/"first_passage_validation.csv",fpt)
    write_csv(tables/"survival_curves.csv",curves);write_csv(tables/"nucleation_scores_by_seed.csv",nucleation)
    write_csv(tables/"nucleation_rates.csv",nucleation_rates(nuc_records,data))
    write_csv(tables/"occupancy_support.csv",support);write_csv(tables/"summary.csv",summary)
    write_csv(tables/"paired_differences.csv",paired_differences(score,nucleation))
    # Existing validated benchmark copied verbatim into a compact comparison table.
    benchmark=[]
    with (args.qj_results/"tables"/"splitting_scores.csv").open() as f:
        for r in csv.DictReader(f):benchmark.append(dict(source="existing_closure_qJ",test="splitting",**r))
    with (args.qj_results/"tables"/"first_passage_validation.csv").open() as f:
        for r in csv.DictReader(f):benchmark.append(dict(source="existing_closure_qJ",test="first_passage",**r))
    write_csv(tables/"qJ_benchmark.csv",benchmark)
    figures(data,ck,score,fpt,curves,nucleation,figdir)
    config={"train_seeds":list(TRAIN),"test_seeds":list(TEST),"burn_time":200,"dt":DT,"lags":list(LAGS),
        "dominance_threshold":DOM,"past_dm_lag":.1,"state_discretization":disc.json(),
        "empty_chamber_polarization":"undefined NaN; boundary strata use occupied group only",
        "decision_scope":"nucleation separate from recruitment; no fitted potential", "censored_corridor_excursions":censored,
        "completed_corridor_excursions":len(durations),"raw_parameter_points":1}
    (out/"analysis_config.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
