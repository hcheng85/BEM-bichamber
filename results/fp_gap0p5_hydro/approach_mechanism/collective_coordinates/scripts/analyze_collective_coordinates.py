#!/usr/bin/env python3
"""Test collective angular mismatch and wall-normal extension before nucleation.

The script is intentionally data-driven only at the level of prespecified
coordinates.  It does not infer a potential.  Controls are prospective and
matched within seed on chamber and vertical oscillation phase only.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parent
TURNING_SCRIPTS = HERE.parent.parent / "pre_crossing_turning" / "scripts"
sys.path.insert(0, str(TURNING_SCRIPTS))
from replay_dynamics import (  # noqa: E402
    Parameters, analytic_geometry, boundary_hydrodynamics,
    direct_hydrodynamics, jumping_mask, social_terms, wall_avoidance,
)

LEADS = (-2.0, -1.0)
DT_FEATURE = 0.01
NATIVE_DT = 0.001
TRAIN_SEEDS = tuple(range(15))
TEST_SEEDS = tuple(range(15, 20))
ALL_SEEDS = tuple(range(20))


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def coordinate_series(state, rdot, stride=10):
    """Return prespecified collective coordinates at a fixed native stride."""
    st = state[:, :, ::stride]
    vel = rdot[:, :, ::stride]
    x, y, th = st
    left = np.all(x < 0, axis=0)
    right = np.all(x >= 0, axis=0)
    exact = left | right
    side = np.where(left, 1.0, np.where(right, -1.0, np.nan))

    xc = np.mean(x, axis=0); yc = np.mean(y, axis=0)
    hx = np.mean(np.cos(th), axis=0); hy = np.mean(np.sin(th), axis=0)
    local_hx = side * hx
    local_xc = side * xc
    distance = -local_xc
    target = np.arctan2(-yc, distance)
    heading = np.arctan2(hy, local_hx)
    psi = wrap(heading - target)

    dx = x - xc[None, :]; dy = y - yc[None, :]
    cxx = np.mean(dx * dx, axis=0)
    cyy = np.mean(dy * dy, axis=0)
    cxy = np.mean(dx * dy, axis=0)
    disc = np.sqrt(np.maximum((cxx-cyy)**2 + 4*cxy*cxy, 0))
    eig_hi = np.maximum((cxx+cyy+disc)/2, 0)
    eig_lo = np.maximum((cxx+cyy-disc)/2, 1e-12)
    aspect = np.sqrt(eig_hi/eig_lo)
    wall_extent = np.sqrt(np.maximum(cxx, 0))
    lam = wall_extent/distance

    radius = np.sqrt(dx*dx+dy*dy)
    speed = np.sqrt(vel[0]*vel[0]+vel[1]*vel[1])
    cross = dx*vel[1]-dy*vel[0]
    denom = radius*speed
    ell_lab = np.mean(np.divide(cross, denom, out=np.zeros_like(cross),
                                where=denom > 1e-12), axis=0)
    ell = side*ell_lab
    pol = np.hypot(hx, hy)

    # Vertical phase uses the actual center position and velocity.  Scaling is
    # seed-local, and the pair (sin phase, cos phase) is used for matching.
    vy = np.gradient(yc, DT_FEATURE)
    sy = max(float(np.nanstd(yc[exact])), 1e-9)
    sv = max(float(np.nanstd(vy[exact])), 1e-9)
    phase = np.arctan2(vy/sv, yc/sy)
    for arr in (psi, lam, aspect, ell, pol, distance, phase, xc, yc, vy):
        arr[~exact] = np.nan
    return dict(time=np.arange(st.shape[2])*DT_FEATURE, exact=exact,
                side=side, psi=psi, Lambda=lam, ell=ell, M=np.abs(ell),
                aspect=aspect, polarization=pol, distance=distance,
                phase=phase, xc=xc, yc=yc, vy=vy)


def load_or_make_cache(seed, raw_dir, cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir/f"collective_coordinates_velocity_v2_seed{seed:02d}.npz"
    if path.exists():
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    raw = raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
    data=sio.loadmat(raw, variable_names=["state","rdot"],squeeze_me=True)
    state=np.asarray(data["state"]); rdot=np.asarray(data["rdot"])
    z = coordinate_series(state,rdot)
    np.savez_compressed(path, **z)
    return z


def valid_control_mask(z, event_times, pre=2.2, post=3.0):
    exact = z["exact"].astype(bool); n = len(exact)
    a = int(round(pre/DT_FEATURE)); b = int(round(post/DT_FEATURE))
    bad = np.flatnonzero(~exact)
    diff = np.zeros(n+1, int)
    for j in bad:
        lo=max(0, j-b); hi=min(n, j+a+1); diff[lo]+=1; diff[hi]-=1
    valid = np.cumsum(diff[:-1]) == 0
    valid[:a] = False; valid[-b:] = False
    valid &= z["time"] >= 202.2
    for t in event_times:
        valid &= np.abs(z["time"]-t) > 5.0
    valid &= (np.arange(n) % 10) == 0
    return valid


def phase_vector(z, idx):
    return np.array([np.sin(z["phase"][idx]), np.cos(z["phase"][idx])])


def select_event_control_pairs(episodes, caches):
    """One prospective phase-matched control per qualified first crossing."""
    rows=[]
    for seed in ALL_SEEDS:
        z=caches[seed]
        ev=[r for r in episodes if int(r["seed"])==seed
            and float(r["full_dwell_before"]) >= 2.2
            and 202.2 <= float(r["first_cross_time"]) <= 997.0]
        ev.sort(key=lambda r: float(r["first_cross_time"]))
        # Remove strongly overlapping nucleations without using their outcome.
        kept=[]
        for r in ev:
            if not kept or float(r["first_cross_time"])-float(kept[-1]["first_cross_time"]) >= 5:
                kept.append(r)
        event_times=[float(r["first_cross_time"]) for r in episodes if int(r["seed"])==seed]
        valid=valid_control_mask(z,event_times)
        available=set(np.flatnonzero(valid).tolist())
        for pair_id,e in enumerate(kept):
            t=float(e["first_cross_time"]); idx=int(round((t-2)/DT_FEATURE))
            origin=e["origin"]; side=1 if origin=="left" else -1
            cand=np.array([j for j in available if int(z["side"][j])==side],int)
            if not len(cand):
                continue
            target=phase_vector(z,idx)
            dist=np.linalg.norm(np.column_stack((np.sin(z["phase"][cand]),
                                                  np.cos(z["phase"][cand])))-target,axis=1)
            j=int(cand[np.argmin(dist)])
            # Controls are unique and separated by 0.5 time units.  A wider
            # exclusion depleted one chamber's phase support in a few seeds.
            available={k for k in available if abs(k-j)>50}
            for condition,anchor,match_distance in (("event",t,0.0),
                                                     ("control",float(z["time"][j])+2,float(np.min(dist)))):
                rows.append(dict(seed=seed,pair_id=pair_id,condition=condition,
                                 origin=origin,anchor_time=anchor,
                                 phase_match_distance=match_distance,
                                 event_success=e["success"] if condition=="event" else ""))
    return rows


def sample_rows(pairs, caches):
    rows=[]
    for p in pairs:
        z=caches[int(p["seed"])]
        for lead in LEADS:
            idx=int(round((float(p["anchor_time"])+lead)/DT_FEATURE))
            row=dict(p,lead=lead,index=idx)
            for key in ("psi","Lambda","ell","M","aspect","polarization",
                        "distance","phase","yc","vy"):
                row[key]=float(z[key][idx])
            row["cos_psi"]=float(np.cos(row["psi"])); row["sin_psi"]=float(np.sin(row["psi"]))
            rows.append(row)
    return rows


def auc_score(y, score):
    y=np.asarray(y,int); score=np.asarray(score,float)
    order=np.argsort(score,kind="mergesort"); ranks=np.empty(len(score),float)
    ranks[order]=np.arange(1,len(score)+1)
    # Mean ranks for ties.
    vals,inv=np.unique(score,return_inverse=True)
    for k in range(len(vals)): ranks[inv==k]=np.mean(ranks[inv==k])
    n1=np.sum(y==1); n0=np.sum(y==0)
    return float((np.sum(ranks[y==1])-n1*(n1+1)/2)/(n1*n0))


def fit_logistic(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    mu=x.mean(0); sd=x.std(0); sd[sd<1e-9]=1
    z=(x-mu)/sd; X=np.column_stack((np.ones(len(z)),z))
    def obj(b):
        a=X@b
        return float(np.sum(np.logaddexp(0,a)-y*a)+0.5e-4*np.sum(b[1:]**2))
    res=minimize(obj,np.zeros(X.shape[1]),method="BFGS")
    return res.x,mu,sd


def predict_logistic(model,x):
    b,mu,sd=model; a=np.column_stack((np.ones(len(x)),(x-mu)/sd))@b
    return 1/(1+np.exp(-np.clip(a,-30,30)))


MODELS={
    "psi_Lambda": ("cos_psi","sin_psi","Lambda"),
    "angular_momentum": ("ell","M"),
    "aspect_ratio": ("aspect",),
}


def predictive_tests(samples):
    rows=[]
    for lead in LEADS:
        train=[r for r in samples if r["lead"]==lead and int(r["seed"]) in TRAIN_SEEDS]
        test=[r for r in samples if r["lead"]==lead and int(r["seed"]) in TEST_SEEDS]
        for name,keys in MODELS.items():
            Xtr=np.array([[r[k] for k in keys] for r in train]); ytr=np.array([r["condition"]=="event" for r in train])
            Xte=np.array([[r[k] for k in keys] for r in test]); yte=np.array([r["condition"]=="event" for r in test])
            model=fit_logistic(Xtr,ytr); prob=predict_logistic(model,Xte)
            rows.append(dict(lead=lead,model=name,features=";".join(keys),
                train_seeds=len(TRAIN_SEEDS),test_seeds=len(TEST_SEEDS),n_train=len(train),n_test=len(test),
                auc=auc_score(yte,prob),brier=float(np.mean((prob-yte)**2)),
                nll=float(-np.mean(yte*np.log(prob+1e-12)+(1-yte)*np.log(1-prob+1e-12)))) )
    return rows


def paired_statistics(samples):
    rows=[]
    for lead in LEADS:
        for metric in ("cos_psi","Lambda","ell","M","aspect"):
            by={}
            for r in samples:
                if r["lead"]==lead: by[(int(r["seed"]),int(r["pair_id"]),r["condition"])]=float(r[metric])
            diffs_by_seed=[]
            for seed in ALL_SEEDS:
                d=[]
                ids=sorted({k[1] for k in by if k[0]==seed})
                for i in ids:
                    if (seed,i,"event") in by and (seed,i,"control") in by:
                        d.append(by[(seed,i,"event")]-by[(seed,i,"control")])
                if d: diffs_by_seed.append((seed,float(np.mean(d))))
            vals=np.array([v for _,v in diffs_by_seed]); rng=np.random.default_rng(20261022)
            boot=np.mean(vals[rng.integers(0,len(vals),(20000,len(vals)))],axis=1)
            rows.append(dict(lead=lead,metric=metric,event_minus_control=float(np.mean(vals)),
                seed_bootstrap_ci_low=float(np.percentile(boot,2.5)),
                seed_bootstrap_ci_high=float(np.percentile(boot,97.5)),independent_seeds=len(vals)))
    return rows


def ridge_fit(x,y,alpha=1e-4):
    mu=x.mean(0); sd=x.std(0); sd[sd<1e-9]=1
    X=np.column_stack((np.ones(len(x)),(x-mu)/sd)); reg=np.eye(X.shape[1]); reg[0,0]=0
    b=np.linalg.solve(X.T@X+alpha*reg,X.T@y)
    return b,mu,sd


def ridge_predict(model,x):
    b,mu,sd=model
    return np.column_stack((np.ones(len(x)),(x-mu)/sd))@b


def memory_rows(caches, tau=0.1):
    k=int(round(tau/DT_FEATURE)); rows=[]
    for seed,z in caches.items():
        n=len(z["time"])
        # Origins sampled every 0.2 units; require a continuous exact-side window.
        for i in range(k,n-k,20):
            if not np.all(z["exact"][i-k:i+k+1]): continue
            if not (z["side"][i-k]==z["side"][i]==z["side"][i+k]): continue
            dpsi_prev=float(wrap(z["psi"][i]-z["psi"][i-k]))
            dpsi_next=float(wrap(z["psi"][i+k]-z["psi"][i]))
            rows.append(dict(seed=seed,psi=float(z["psi"][i]),Lambda=float(z["Lambda"][i]),
                sinpsi=float(np.sin(z["psi"][i])),cospsi=float(np.cos(z["psi"][i])),
                phase=float(z["phase"][i]),sinphase=float(np.sin(z["phase"][i])),
                cosphase=float(np.cos(z["phase"][i])),dpsi_prev=dpsi_prev,
                dLambda_prev=float(z["Lambda"][i]-z["Lambda"][i-k]),
                dpsi_next=dpsi_next,dLambda_next=float(z["Lambda"][i+k]-z["Lambda"][i])))
    return rows


def memory_tests(rows):
    train=[r for r in rows if r["seed"] in TRAIN_SEEDS]; test=[r for r in rows if r["seed"] in TEST_SEEDS]
    specs={
        "state_only":("sinpsi","cospsi","Lambda"),
        "state_plus_phase":("sinpsi","cospsi","Lambda","sinphase","cosphase"),
        "state_plus_history":("sinpsi","cospsi","Lambda","dpsi_prev","dLambda_prev"),
        "state_phase_history":("sinpsi","cospsi","Lambda","sinphase","cosphase","dpsi_prev","dLambda_prev"),
    }
    out=[]
    ykeys=("dpsi_next","dLambda_next")
    for name,keys in specs.items():
        xtr=np.array([[r[k] for k in keys] for r in train]); xte=np.array([[r[k] for k in keys] for r in test])
        ytr=np.array([[r[k] for k in ykeys] for r in train]); yte=np.array([[r[k] for k in ykeys] for r in test])
        pred=ridge_predict(ridge_fit(xtr,ytr),xte)
        for j,target in enumerate(ykeys):
            rmse=float(np.sqrt(np.mean((pred[:,j]-yte[:,j])**2)))
            scale=float(np.std(yte[:,j])); r2=1-float(np.sum((pred[:,j]-yte[:,j])**2)/np.sum((yte[:,j]-np.mean(yte[:,j]))**2))
            out.append(dict(model=name,target=target,n_train=len(train),n_test=len(test),rmse=rmse,
                            normalized_rmse=rmse/scale,r2=r2))
    # Fixed-phase fits: current state versus current state + history in four bins.
    for q in range(4):
        lo=-np.pi+q*np.pi/2; hi=lo+np.pi/2
        tr=[r for r in train if lo <= r["phase"] < hi]; te=[r for r in test if lo <= r["phase"] < hi]
        for name,keys in (("fixed_phase_state",("sinpsi","cospsi","Lambda")),
                          ("fixed_phase_history",("sinpsi","cospsi","Lambda","dpsi_prev","dLambda_prev"))):
            xtr=np.array([[r[k] for k in keys] for r in tr]); xte=np.array([[r[k] for k in keys] for r in te])
            ytr=np.array([[r[k] for k in ykeys] for r in tr]); yte=np.array([[r[k] for k in ykeys] for r in te])
            pred=ridge_predict(ridge_fit(xtr,ytr),xte)
            for j,target in enumerate(ykeys):
                rmse=float(np.sqrt(np.mean((pred[:,j]-yte[:,j])**2))); scale=float(np.std(yte[:,j]))
                out.append(dict(model=name,target=target,phase_bin=q,n_train=len(tr),n_test=len(te),rmse=rmse,
                                normalized_rmse=rmse/scale,r2=1-float(np.sum((pred[:,j]-yte[:,j])**2)/np.sum((yte[:,j]-np.mean(yte[:,j]))**2))))
    return out


def collective_heading_rate(theta, omega):
    z=np.mean(np.exp(1j*theta)); dz=np.mean(1j*np.exp(1j*theta)*omega)
    return float(np.imag(np.conj(z)*dz)/(abs(z)**2))


def target_rate(state, velocity, side):
    xc=float(np.mean(state[0])); yc=float(np.mean(state[1]))
    ux=side*float(np.mean(velocity[0])); uy=float(np.mean(velocity[1]))
    d=-side*xc
    return float((-d*uy-yc*ux)/(d*d+yc*yc))


def lambda_rate(state, velocity, side):
    xx=side*state[0]; uu=side*velocity[0]
    xc=float(np.mean(xx)); uc=float(np.mean(uu)); dev=xx-xc
    sigma=float(np.sqrt(np.mean(dev*dev))); d=-xc
    sigmadot=float(np.mean(dev*(uu-uc))/sigma) if sigma>1e-12 else 0.0
    ddot=-uc
    return float(sigmadot/d-sigma*ddot/(d*d))


def replay_selected(pairs, raw_dir, boundary, lu, piv, p):
    """Replay one representative event/control pair per seed at both leads."""
    geom=analytic_geometry(p); selected=[]
    for seed in ALL_SEEDS:
        ep=[r for r in pairs if int(r["seed"])==seed and r["condition"]=="event"]
        if not ep: continue
        # Deterministic and reproducible: temporally central qualified event.
        e=ep[len(ep)//2]
        selected.extend([r for r in pairs if int(r["seed"])==seed and int(r["pair_id"])==int(e["pair_id"])])
    rows=[]; validation=[]
    for seed in ALL_SEEDS:
        chosen=[r for r in selected if int(r["seed"])==seed]
        if not chosen: continue
        print(f"replay seed {seed}", flush=True)
        raw=raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state_all=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
        for a in chosen:
            side=1.0 if a["origin"]=="left" else -1.0
            for lead in LEADS:
                idx=int(round((float(a["anchor_time"])+lead)/p.dt))
                state=state_all[:,:,idx]; nxt=state_all[:,:,idx+1]
                align,attract,_=social_terms(state,geom,p)
                ud,wd=direct_hydrodynamics(state,p); ub,wb=boundary_hydrodynamics(state,boundary,lu,piv,p)
                ww,wall_distance=wall_avoidance(state,geom,p)
                terms={"wall":ww,"alignment":align,"attraction":attract,
                       "direct_hydro":wd,"boundary_hydro":wb}
                total=np.clip(sum(terms.values()),-np.pi/p.dt,np.pi/p.dt)
                selfv=np.vstack((np.cos(state[2]),np.sin(state[2])))
                rdot=selfv+ud+ub
                overridden=jumping_mask(state,rdot,geom,p) | np.any(np.abs(rdot)>500,axis=0)
                total[overridden]=np.clip(ww+align+attract,-np.pi/p.dt,np.pi/p.dt)[overridden]
                wd_eff=wd.copy(); wb_eff=wb.copy(); wd_eff[overridden]=0; wb_eff[overridden]=0
                terms["direct_hydro"]=wd_eff; terms["boundary_hydro"]=wb_eff
                obs=nxt[2]-state[2]; noise=obs-total*p.dt
                hx=np.mean(np.cos(state[2])); hy=np.mean(np.sin(state[2]))
                theta_local=np.arctan2(hy,side*hx)
                xc=np.mean(state[0]); yc=np.mean(state[1]); target=np.arctan2(-yc,-side*xc)
                psi=float(wrap(theta_local-target))
                # Nearest wall-following tangent to the outer circle at COM.
                radial=np.arctan2(yc,side*xc)
                tangents=np.array([wrap(radial+np.pi/2),wrap(radial-np.pi/2)])
                terr=wrap(theta_local-tangents); wall_err=float(terr[np.argmin(np.abs(terr))])
                row=dict(seed=seed,pair_id=a["pair_id"],condition=a["condition"],origin=a["origin"],
                         lead=lead,psi=psi,Lambda=float(np.std(state[0])/abs(xc)),
                         wall_following_error=wall_err,mean_wall_distance=float(np.mean(wall_distance)),
                         overridden=int(np.sum(overridden)))
                for key,w in terms.items():
                    rate=side*collective_heading_rate(state[2],w)
                    row[f"psi_{key}"]=rate
                    row[f"opening_closing_{key}"]=-np.sign(psi)*rate if psi != 0 else 0.0
                    row[f"wall_following_closing_{key}"]=-np.sign(wall_err)*rate if wall_err != 0 else 0.0
                for key,u in (("self",selfv),("direct_hydro",ud),("boundary_hydro",ub)):
                    ar=target_rate(state,u,side)
                    row[f"psi_target_{key}"]=-ar
                    row[f"Lambda_{key}"]=lambda_rate(state,u,side)
                row["psi_deterministic_total"]=sum(row[f"psi_{k}"] for k in terms)+sum(row[f"psi_target_{k}"] for k in ("self","direct_hydro","boundary_hydro"))
                row["Lambda_deterministic_total"]=sum(row[f"Lambda_{k}"] for k in ("self","direct_hydro","boundary_hydro"))
                rows.append(row)
                validation.append(dict(seed=seed,condition=a["condition"],lead=lead,
                    noise_z_mean=float(np.mean(noise/(p.In*np.sqrt(p.dt)))),
                    noise_z_sd=float(np.std(noise/(p.In*np.sqrt(p.dt)),ddof=1)),
                    rms_update_residual=float(np.sqrt(np.mean((obs-total*p.dt-noise)**2))),overridden=int(np.sum(overridden))))
        del state_all
    return selected,rows,validation


def make_figures(out,samples,pred,paired,replay,memory):
    colors={"event":"#0072B2","control":"#999999"}
    fig,ax=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
    metrics=[("cos_psi",r"heading toward opening, $\cos\psi$"),("Lambda",r"relative wall-normal extent, $\Lambda$"),
             ("M",r"angular momentum magnitude, $M$"),("aspect",r"aspect ratio, $A$")]
    for a,(key,label) in zip(ax.flat,metrics):
        for cond in colors:
            vals=[]
            for lead in LEADS:
                seedmeans=[]
                for seed in ALL_SEEDS:
                    q=[float(r[key]) for r in samples if r["condition"]==cond and r["lead"]==lead and int(r["seed"])==seed]
                    if q: seedmeans.append(np.mean(q))
                vals.append(np.mean(seedmeans))
            a.plot([-x for x in LEADS],vals,"o-",color=colors[cond],label=cond)
        a.set_xlabel("time before first crossing / pseudo-anchor"); a.set_ylabel(label); a.grid(alpha=.2)
    ax[0,0].legend(frameon=False)
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"01_coordinate_precursors.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    names=list(MODELS)
    for j,lead in enumerate(LEADS):
        q=[r for r in pred if r["lead"]==lead]
        axes[j].bar(np.arange(len(q)),[r["auc"] for r in q],color=["#0072B2","#D55E00","#999999"])
        axes[j].set_xticks(np.arange(len(q)),[r["model"].replace("_","\n") for r in q]); axes[j].set_ylim(.35,1)
        axes[j].axhline(.5,color="k",lw=.8); axes[j].set_ylabel("held-out-seed AUC"); axes[j].set_title(f"lead {abs(lead):g}")
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"02_heldout_prediction.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    keys=["psi_wall","psi_alignment","psi_attraction","psi_direct_hydro","psi_boundary_hydro",
          "psi_target_self","psi_target_direct_hydro","psi_target_boundary_hydro"]
    labels=["wall","alignment","attraction","direct hydro","BEM hydro","target: self","target: direct","target: BEM"]
    for c,cond in enumerate(("event","control")):
        rr=[r for r in replay if r["condition"]==cond]
        width=.1; x=np.arange(len(LEADS))
        for j,(key,label) in enumerate(zip(keys,labels)):
            axes[c].bar(x+(j-3.5)*width,[np.median([float(r[key]) for r in rr if float(r["lead"])==lead]) for lead in LEADS],width,label=label)
        axes[c].axhline(0,color="k",lw=.8); axes[c].set_xticks(x,[str(v) for v in LEADS]); axes[c].set_title(cond)
        axes[c].set_ylabel(r"median contribution to $\dot\psi$"); axes[c].set_xlabel("lead")
    axes[1].legend(frameon=False,fontsize=7,ncol=2)
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"03_psi_dynamics_decomposition.{ext}",dpi=220)
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    base={r["target"]:r for r in memory if r["model"]=="state_only"}
    phase={r["target"]:r for r in memory if r["model"]=="state_plus_phase"}
    hist={r["target"]:r for r in memory if r["model"]=="state_plus_history"}
    both={r["target"]:r for r in memory if r["model"]=="state_phase_history"}
    for a,target in zip(axes,("dpsi_next","dLambda_next")):
        a.bar(np.arange(4),[base[target]["normalized_rmse"],phase[target]["normalized_rmse"],hist[target]["normalized_rmse"],both[target]["normalized_rmse"]],color=["#999999","#56B4E9","#E69F00","#009E73"])
        a.set_xticks(np.arange(4),["state","+phase","+history","+both"],rotation=20); a.set_ylabel("held-out normalized RMSE"); a.set_title(target); a.set_ylim(0,1.1)
    for ext in ("png","pdf"): fig.savefig(out/"figures"/f"04_short_time_memory.{ext}",dpi=220)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--analysis-dir",type=Path,required=True)
    ap.add_argument("--raw-dir",type=Path,required=True)
    ap.add_argument("--cache-dir",type=Path,required=True)
    ap.add_argument("--geometry-file",type=Path,required=True)
    ap.add_argument("--lu-file",type=Path,required=True)
    ap.add_argument("--piv-file",type=Path,required=True)
    ap.add_argument("--skip-replay",action="store_true",help="reuse existing deterministic snapshot tables")
    args=ap.parse_args(); out=args.analysis_dir.resolve(); (out/"tables").mkdir(parents=True,exist_ok=True); (out/"figures").mkdir(exist_ok=True)
    caches={}
    for seed in ALL_SEEDS:
        print(f"coordinates seed {seed}",flush=True)
        caches[seed]=load_or_make_cache(seed,args.raw_dir,args.cache_dir)
    episodes=read_csv(out.parent/"tables"/"nucleation_episodes.csv")
    pairs=select_event_control_pairs(episodes,caches); samples=sample_rows(pairs,caches)
    pred=predictive_tests(samples); paired=paired_statistics(samples)
    memrows=memory_rows(caches); mem=memory_tests(memrows)

    if args.skip_replay:
        selected=read_csv(out/"tables"/"selected_replay_pairs.csv")
        replay=read_csv(out/"tables"/"deterministic_coordinate_dynamics.csv")
        validation=read_csv(out/"tables"/"replay_validation.csv")
    else:
        g=sio.loadmat(args.geometry_file); normal=np.asarray(g["normal_vec"])
        boundary=dict(x=np.asarray(g["Sx"]).reshape(-1),y=np.asarray(g["Sy"]).reshape(-1),
                      weights=np.asarray(g["weights"]).reshape(-1),nx=normal[:,0],ny=normal[:,1])
        lu=np.load(args.lu_file,mmap_mode="r"); piv=np.load(args.piv_file)
        selected,replay,validation=replay_selected(pairs,args.raw_dir,boundary,lu,piv,Parameters())

    write_csv(out/"tables"/"matched_event_controls.csv",pairs)
    write_csv(out/"tables"/"coordinate_samples.csv",samples)
    write_csv(out/"tables"/"paired_coordinate_statistics.csv",paired)
    write_csv(out/"tables"/"heldout_prediction.csv",pred)
    write_csv(out/"tables"/"selected_replay_pairs.csv",selected)
    write_csv(out/"tables"/"deterministic_coordinate_dynamics.csv",replay)
    write_csv(out/"tables"/"replay_validation.csv",validation)
    write_csv(out/"tables"/"memory_validation.csv",mem)
    make_figures(out,samples,pred,paired,replay,mem)
    config=dict(dataset="20 existing gap-50 full-BEM hydrodynamic trajectories",new_simulations=False,
        leads=list(LEADS),training_seeds=list(TRAIN_SEEDS),heldout_seeds=list(TEST_SEEDS),
        control_matching="within seed and origin chamber; nearest vertical phase at A-2; no heading, extension, angular momentum, aspect ratio, or outcome matching",
        coordinate_definitions={"psi":"wrapped reflected collective heading minus COM-to-opening direction",
          "Lambda":"sqrt(C_xx)/abs(x_center)","ell":"Filella signed normalized angular momentum using saved full rdot after chamber reflection",
          "M":"abs(ell)","A":"sqrt(largest covariance eigenvalue/smallest covariance eigenvalue)"},
        replay="one temporally central event-control pair per seed at leads 2 and 1; 80 snapshots total",
        parameters=Parameters().__dict__)
    (out/"configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    print(json.dumps(dict(pairs=len(pairs),samples=len(samples),memory_samples=len(memrows),replay_snapshots=len(replay)),indent=2))


if __name__=="__main__": main()
