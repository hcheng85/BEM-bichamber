#!/usr/bin/env python3
"""First-principles kinetic first exit for one source-side approach attempt.

The PDE discretization is a monotone continuous-time lattice generator with
outgoing kinetic boundary conditions.  Archived trajectories are read only for
validation; no microscopic simulation or fitted coefficient is used.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.special import ndtr

from analyze_aperture_capture_current import HALF_WIDTH, IN_NOISE, capture_geometry
from analyze_approach_current_balance import DT, ELL, N, load_gap20_production, raw_path, replay_torques
from analyze_initiation_propagation import RawReader


OMEGA_NODES=np.asarray([-3.,-2.,-1.,0.,1.,2.,3.])
FINAL_GRID=(31,19,36)
CONVERGENCE_GRIDS=((15,9,16),(23,13,24),(31,19,36))
NOISE_LEVELS=(0.5,0.25,0.125,0.0625,0.0)
SHORT_FLIGHT_CUTOFFS=(0.25,0.5,1.0,2.0,5.0)


def write_csv(path: Path,rows: list[dict]) -> None:
    if not rows:return
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n");writer.writeheader();writer.writerows(rows)


def grid_coordinates(shape: tuple[int,int,int]) -> tuple[np.ndarray,np.ndarray,np.ndarray,float,float,float]:
    nz,ny,np_=shape;dz=ELL/nz;dy=2*HALF_WIDTH/ny;dp=2*np.pi/np_
    z=-ELL+(np.arange(nz)+.5)*dz;y=-HALF_WIDTH+(np.arange(ny)+.5)*dy;phi=-np.pi+(np.arange(np_)+.5)*dp
    return z,y,phi,dz,dy,dp


def solve_backward(shape: tuple[int,int,int],noise: float,omega: float=0.) -> dict:
    """Solve committor, capture-weighted time, and exit time on a CTMC grid."""
    nz,ny,np_=shape;z,y,phi,dz,dy,dp=grid_coordinates(shape);n=nz*ny*np_
    rows=[];cols=[];data=[];rhs=np.zeros(n)
    def idx(i,j,k):return (i*ny+j)*np_+k
    diffusion=.5*noise*noise/(dp*dp)
    for i in range(nz):
        for j in range(ny):
            for k,p in enumerate(phi):
                row=idx(i,j,k);diag=0.;u=math.cos(p);v=math.sin(p)
                rz=abs(u)/dz
                if rz:
                    diag+=rz;ii=i+(1 if u>0 else -1)
                    if 0<=ii<nz:rows.append(row);cols.append(idx(ii,j,k));data.append(-rz)
                    elif u>0:rhs[row]+=rz  # physical gap boundary
                ry=abs(v)/dy
                if ry:
                    diag+=ry;jj=j+(1 if v>0 else -1)
                    if 0<=jj<ny:rows.append(row);cols.append(idx(i,jj,k));data.append(-ry)
                    # lateral exits are failure, so no RHS contribution
                if diffusion:
                    for kk in ((k-1)%np_,(k+1)%np_):
                        diag+=diffusion;rows.append(row);cols.append(idx(i,j,kk));data.append(-diffusion)
                ra=abs(omega)/dp
                if ra:
                    kk=(k+(1 if omega>0 else -1))%np_;diag+=ra;rows.append(row);cols.append(idx(i,j,kk));data.append(-ra)
                rows.append(row);cols.append(row);data.append(diag)
    matrix=coo_matrix((data,(rows,cols)),shape=(n,n)).tocsc();factor=splu(matrix)
    q=np.clip(factor.solve(rhs),0,1);m=np.maximum(factor.solve(q),0);texit=np.maximum(factor.solve(np.ones(n)),0)
    return dict(shape=shape,z=z,y=y,phi=phi,q=q.reshape(shape),m=m.reshape(shape),texit=texit.reshape(shape),
                noise=noise,omega=omega,matrix_nnz=int(matrix.nnz))


def ballistic_indicator(z: np.ndarray,y: np.ndarray,phi: np.ndarray) -> np.ndarray:
    u=np.cos(phi);hit=y-z*np.sin(phi)/u
    return ((u>0)&(np.abs(hit)<=HALF_WIDTH)).astype(float)


def grid_ballistic_error(solution: dict) -> dict:
    zz,yy,pp=np.meshgrid(solution["z"],solution["y"],solution["phi"],indexing="ij")
    exact=ballistic_indicator(zz,yy,pp);q=solution["q"]
    return dict(l1=float(np.mean(np.abs(q-exact))),classification_error=float(np.mean((q>=.5)!=(exact>.5))),
                interior_l1=float(np.mean(np.abs(q-exact)[np.abs(np.abs(yy-zz*np.tan(pp))-HALF_WIDTH)>.15])))


def make_interpolator(solutions: list[dict],field: str):
    omega=np.asarray([s["omega"] for s in solutions]);base=solutions[0];phi=base["phi"]
    values=np.stack([s[field] for s in solutions]);values=np.concatenate([values,values[:,:,:,0:1]],axis=3)
    phi_ext=np.r_[phi,phi[0]+2*np.pi]
    interp=RegularGridInterpolator((omega,base["z"],base["y"],phi_ext),values,bounds_error=False,fill_value=None)
    def evaluate(om,z,y,p):
        om=np.clip(np.asarray(om),omega[0],omega[-1]);z=np.clip(np.asarray(z),base["z"][0],base["z"][-1]);y=np.clip(np.asarray(y),base["y"][0],base["y"][-1])
        p=(np.asarray(p)-phi[0])%(2*np.pi)+phi[0]
        return interp(np.column_stack([om,z,y,p]))
    return evaluate


def leading_straight_asymptotic(z,y,phi,noise):
    u=np.cos(phi);out=np.zeros_like(np.asarray(z,dtype=float));valid=u>0
    tau=-np.asarray(z)[valid]/u[valid];hit=np.asarray(y)[valid]-np.asarray(z)[valid]*np.tan(np.asarray(phi)[valid])
    sigma=noise*np.sqrt(tau**3/(3*u[valid]**2))
    deterministic=sigma<1e-14;prob=np.empty_like(tau);prob[deterministic]=(np.abs(hit[deterministic])<=HALF_WIDTH)
    prob[~deterministic]=ndtr((HALF_WIDTH-hit[~deterministic])/sigma[~deterministic])-ndtr((-HALF_WIDTH-hit[~deterministic])/sigma[~deterministic])
    out[valid]=prob
    return out


def convergence_audit(output: Path) -> tuple[list[dict],list[dict],dict]:
    rows=[];noise_rows=[];solutions={}
    probes=np.asarray([[-1.0,0.,0.],[-1.0,1.0,.4],[-1.0,1.1,.4],[-2.0,0.,.5],[-.4,-1.1,-.35]])
    for shape in CONVERGENCE_GRIDS:
        sol=solve_backward(shape,IN_NOISE,0.);solutions[(shape,IN_NOISE,0.)]=sol
        interp=make_interpolator([sol],"q");values=interp(np.zeros(len(probes)),probes[:,0],probes[:,1],probes[:,2])
        err=grid_ballistic_error(sol)
        rows.append(dict(nz=shape[0],ny=shape[1],nphi=shape[2],noise=IN_NOISE,omega=0.,nnz=sol["matrix_nnz"],**err,
                         **{f"probe_{i}":float(v) for i,v in enumerate(values)}))
    for i,row in enumerate(rows):
        row["max_probe_change_from_previous"]=None if i==0 else float(max(abs(row[f"probe_{j}"]-rows[i-1][f"probe_{j}"]) for j in range(len(probes))))
    shape=FINAL_GRID
    for noise in NOISE_LEVELS:
        key=(shape,noise,0.);sol=solutions.get(key) or solve_backward(shape,noise,0.);solutions[key]=sol
        err=grid_ballistic_error(sol);zz,yy,pp=np.meshgrid(sol["z"],sol["y"],sol["phi"],indexing="ij")
        tau=np.where(np.cos(pp)>0,-zz/np.cos(pp),np.inf);short=(tau<=1)&(np.cos(pp)>.2)
        approx=leading_straight_asymptotic(zz.ravel(),yy.ravel(),pp.ravel(),noise).reshape(shape)
        exact=ballistic_indicator(zz,yy,pp)
        noise_rows.append(dict(noise=noise,rotational_diffusion=.5*noise*noise,**err,
                               short_flight_asymptotic_rmse=float(np.sqrt(np.mean((sol["q"][short]-approx[short])**2))),
                               short_flight_asymptotic_to_ray_l1=float(np.mean(np.abs(approx[short]-exact[short]))),
                               short_flight_n=int(short.sum())))
    return rows,noise_rows,solutions[(FINAL_GRID,IN_NOISE,0.)]


def segment_exit(z0,y0,z1,y1):
    candidates=[];dz=z1-z0;dy=y1-y0
    if dz>0:
        t=-z0/dz
        if 0<=t<=1 and abs(y0+t*dy)<=HALF_WIDTH+1e-10:candidates.append((t,"gap"))
    if dz<0:
        t=(-ELL-z0)/dz
        if 0<=t<=1:candidates.append((t,"upstream"))
    if dy>0:
        t=(HALF_WIDTH-y0)/dy
        if 0<=t<=1:candidates.append((t,"lateral"))
    if dy<0:
        t=(-HALF_WIDTH-y0)/dy
        if 0<=t<=1:candidates.append((t,"lateral"))
    return min(candidates) if candidates else (1.,"unclassified")


def straight_exit(z: float,y: float,phi: float) -> tuple[float,str]:
    """Exact first boundary reached by a straight inward characteristic."""
    u=math.cos(phi);v=math.sin(phi)
    if u<=0:return float("inf"),"upstream"
    candidates=[(-z/u,"gap")]
    if v>0:candidates.append(((HALF_WIDTH-y)/v,"lateral"))
    elif v<0:candidates.append(((-HALF_WIDTH-y)/v,"lateral"))
    return min(candidates)


def deterministic_arc_exit(z: np.ndarray,y: np.ndarray,phi: np.ndarray,omega: float,dt: float=.002) -> tuple[np.ndarray,np.ndarray]:
    """First exit along exact constant-turn circular characteristics.

    One full angular period suffices: a characteristic still inside after a
    period is a closed orbit and is trapped in the zero-diffusion problem.
    """
    z=np.asarray(z,float);y=np.asarray(y,float);phi=np.asarray(phi,float)
    n=len(z);duration=np.full(n,np.nan);kind=np.full(n,"trapped",dtype=object);alive=np.ones(n,bool)
    period=2*np.pi/abs(omega);previous_z=z.copy();previous_y=y.copy()
    for t in np.arange(dt,period+1.5*dt,dt):
        current_z=z+(np.sin(phi+omega*t)-np.sin(phi))/omega
        current_y=y+(np.cos(phi)-np.cos(phi+omega*t))/omega
        outside=(current_z>=0)|(current_z<-ELL)|(np.abs(current_y)>HALF_WIDTH)
        exiting=np.flatnonzero(alive&outside)
        for i in exiting:
            fraction,label=segment_exit(previous_z[i],previous_y[i],current_z[i],current_y[i])
            duration[i]=t-dt+fraction*dt;kind[i]=label
        alive[exiting]=False
        previous_z[alive]=current_z[alive];previous_y[alive]=current_y[alive]
        if not alive.any():break
    return duration,kind


def constant_turn_audit() -> tuple[list[dict],dict]:
    """Compare zero-noise PDE classifications with exact circular arcs."""
    rng=np.random.default_rng(20261119);n=4000
    z=rng.uniform(-ELL,0,n);y=rng.uniform(-HALF_WIDTH,HALF_WIDTH,n);phi=rng.uniform(-np.pi,np.pi,n)
    rows=[];solutions={}
    for omega in (-1.,1.):
        print(f"zero-noise circular-characteristic audit omega={omega:g}",flush=True)
        sol=solve_backward(FINAL_GRID,0.,omega);solutions[omega]=sol
        duration,kind=deterministic_arc_exit(z,y,phi,omega);exiting=kind!="trapped";capture=kind=="gap"
        q=make_interpolator([sol],"q")(np.full(n,omega),z,y,phi)
        rows.append(dict(omega=omega,n=n,n_exiting=int(exiting.sum()),n_trapped=int((~exiting).sum()),
                         trapped_fraction=float(np.mean(~exiting)),capture_fraction_among_exiting=float(np.mean(capture[exiting])),
                         pde_classification_error_excluding_trapped=float(np.mean((q[exiting]>=.5)!=capture[exiting])),
                         mean_pde_probability_on_trapped=float(np.mean(q[~exiting])) if np.any(~exiting) else float("nan"),
                         arc_median_exit_time=float(np.median(duration[exiting]))))
    return rows,solutions


def extract_attempts(raw_dir: Path,risk_file: Path,source_file: Path,horizon: float=5.) -> list[dict]:
    risk=[r for r in csv.DictReader(risk_file.open()) if int(r["seed"])>=30];by=defaultdict(list)
    for r in risk:by[int(r["seed"])].append(r)
    with contextlib.redirect_stdout(io.StringIO()):mod=load_gap20_production(source_file)
    para=dict(Ip=4,In=.5,R=mod.R,light=1,dt=DT,Iw=.94,If=0,num=N,delta=.01,radius_ratio=50,gap_ratio=20)
    max_steps=int(round(horizon/DT));out=[]
    for count,seed in enumerate(sorted(by),1):
        print(f"archived attempts {count}/{len(by)} seed={seed}",flush=True);reader=RawReader(raw_path(raw_dir,seed))
        try:
            for origin in by[seed]:
                step=int(origin["native_step"]);state=reader.frame("state",step);source_sign=1 if float(origin["q"])>0 else -1;s=-source_sign
                cap=capture_geometry(state,source_sign);indices=np.flatnonzero(cap["inside"]&(cap["u"]>0))
                if not len(indices):continue
                _,_,_,omega_lab=replay_torques(mod,state,para);omega_local=s*omega_lab
                block=reader.state_block(step,step+max_steps);z=s*block[0,indices,:];y=block[1,indices,:]
                inside=(z<0)&(z>=-ELL)&(np.abs(y)<=HALF_WIDTH)
                for row,fish in enumerate(indices):
                    exits=np.flatnonzero(~inside[row,1:])
                    kind="censored";duration=horizon
                    if len(exits):
                        k=int(exits[0]+1);fraction,kind=segment_exit(z[row,k-1],y[row,k-1],z[row,k],y[row,k]);duration=(k-1+fraction)*DT
                    theta=state[2,fish];phi=math.atan2(math.sin(theta),s*math.cos(theta))
                    ballistic_time,ballistic_kind=straight_exit(float(s*state[0,fish]),float(state[1,fish]),phi)
                    out.append(dict(seed=seed,origin_step=step,fish=int(fish),z=float(s*state[0,fish]),y=float(state[1,fish]),phi=phi,
                                    omega_local=float(omega_local[fish]),ballistic_capture=int(cap["mask"][fish]),
                                    ballistic_exit_time=ballistic_time,ballistic_exit_kind=ballistic_kind,
                                    actual_capture=int(kind=="gap"),exit_kind=kind,actual_exit_time=duration,censored=int(kind=="censored")))
        finally:reader.close()
    return out


def attach_predictions(attempts: list[dict],zero: dict,omega_solutions: list[dict]) -> None:
    om=np.asarray([r["omega_local"] for r in attempts]);z=np.asarray([r["z"] for r in attempts]);y=np.asarray([r["y"] for r in attempts]);p=np.asarray([r["phi"] for r in attempts]);zeros=np.zeros(len(attempts))
    fields={}
    for name in ("q","m","texit"):
        fields[("zero",name)]=make_interpolator([zero],name)(zeros,z,y,p)
        fields[("frozen",name)]=make_interpolator(omega_solutions,name)(om,z,y,p)
    for i,row in enumerate(attempts):
        for model in ("zero","frozen"):
            q=float(np.clip(fields[(model,"q")][i],0,1));m=float(max(fields[(model,"m")][i],0));te=float(max(fields[(model,"texit")][i],0))
            row[f"{model}_capture_probability"]=q;row[f"{model}_mean_exit_time"]=te;row[f"{model}_conditional_capture_time"]=m/q if q>1e-8 else float("nan")
        row["omega_clipped"]=int(abs(row["omega_local"])>OMEGA_NODES[-1])


def seed_bootstrap_validation(rows: list[dict],nboot=3000) -> tuple[list[dict],dict]:
    uncensored=[r for r in rows if not r["censored"]];seeds=sorted({r["seed"] for r in uncensored});by={s:[r for r in uncensored if r["seed"]==s] for s in seeds}
    def stats(rr,model):
        y=np.asarray([r["actual_capture"] for r in rr]);q=np.asarray([r[f"{model}_capture_probability"] for r in rr]);t=np.asarray([r["actual_exit_time"] for r in rr]);tp=np.asarray([r[f"{model}_mean_exit_time"] for r in rr])
        captured=[r for r in rr if r["actual_capture"] and np.isfinite(r[f"{model}_conditional_capture_time"])]
        time_bias=float(np.mean(tp-t));time_mae=float(np.mean(np.abs(tp-t)))
        cap_time_bias=float(np.mean([r[f"{model}_conditional_capture_time"]-r["actual_exit_time"] for r in captured])) if captured else float("nan")
        return dict(mean_probability=float(q.mean()),observed_capture=float(y.mean()),calibration_bias=float(q.mean()-y.mean()),brier=float(np.mean((q-y)**2)),exit_time_bias=time_bias,exit_time_mae=time_mae,capture_time_bias=cap_time_bias)
    observed={m:stats(uncensored,m) for m in ("zero","frozen")};rng=np.random.default_rng(20261118);draws={m:defaultdict(list) for m in ("zero","frozen")};diff=defaultdict(list)
    for _ in range(nboot):
        rr=[]
        for s in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[int(s)])
        ss={m:stats(rr,m) for m in ("zero","frozen")}
        for m in ss:
            for k,v in ss[m].items():draws[m][k].append(v)
        for k in ss["zero"]:diff[k].append(ss["zero"][k]-ss["frozen"][k])
    table=[]
    for model in ("zero","frozen"):
        for metric,value in observed[model].items():table.append(dict(model=model,metric=metric,estimate=value,ci95=np.percentile(draws[model][metric],[2.5,97.5]).tolist()))
    comparison={k:dict(zero_minus_frozen=observed["zero"][k]-observed["frozen"][k],ci95=np.percentile(v,[2.5,97.5]).tolist()) for k,v in diff.items()}
    return table,comparison


def short_flight_validation(rows: list[dict],nboot=2000) -> list[dict]:
    """Seed-bootstrap validation stratified only by initial ballistic duration."""
    rng=np.random.default_rng(20261120);out=[]
    for cutoff in SHORT_FLIGHT_CUTOFFS:
        part=[r for r in rows if not r["censored"] and r["ballistic_exit_time"]<=cutoff]
        seeds=sorted({r["seed"] for r in part});by={s:[r for r in part if r["seed"]==s] for s in seeds}
        def stats(rr,model):
            yy=np.asarray([r["actual_capture"] for r in rr]);qq=np.asarray([r[f"{model}_capture_probability"] for r in rr]);tt=np.asarray([r["actual_exit_time"] for r in rr]);tp=np.asarray([r[f"{model}_mean_exit_time"] for r in rr])
            return dict(observed_capture=float(yy.mean()),mean_probability=float(qq.mean()),calibration_bias=float(np.mean(qq-yy)),brier=float(np.mean((qq-yy)**2)),exit_time_mae=float(np.mean(np.abs(tp-tt))))
        observed={m:stats(part,m) for m in ("zero","frozen")};draws={m:defaultdict(list) for m in observed}
        for _ in range(nboot):
            rr=[]
            for s in rng.choice(seeds,len(seeds),replace=True):rr.extend(by[int(s)])
            for model in observed:
                for metric,value in stats(rr,model).items():draws[model][metric].append(value)
        for model in observed:
            for metric,value in observed[model].items():
                out.append(dict(ballistic_time_max=cutoff,model=model,n=len(part),n_seeds=len(seeds),metric=metric,estimate=value,
                                ci95=json.dumps(np.percentile(draws[model][metric],[2.5,97.5]).tolist())))
    return out


def calibration_rows(rows: list[dict],model: str,bins=5) -> list[dict]:
    part=[r for r in rows if not r["censored"]];order=sorted(part,key=lambda r:r[f"{model}_capture_probability"]);groups=np.array_split(np.asarray(order,dtype=object),bins);out=[]
    for i,g in enumerate(groups):
        out.append(dict(model=model,bin=i,n=len(g),n_seeds=len({r["seed"] for r in g}),mean_probability=float(np.mean([r[f"{model}_capture_probability"] for r in g])),observed_capture=float(np.mean([r["actual_capture"] for r in g])),mean_actual_exit_time=float(np.mean([r["actual_exit_time"] for r in g])),mean_predicted_exit_time=float(np.mean([r[f"{model}_mean_exit_time"] for r in g]))))
    return out


def make_figure(zero: dict,noise_rows: list[dict],calibration: list[dict],validation: list[dict],short_validation: list[dict],output: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10,7.6))
    k=int(np.argmin(abs(zero["phi"]-.5)));im=axes[0,0].pcolormesh(zero["z"],zero["y"],zero["q"][:,:,k].T,shading="auto",vmin=0,vmax=1,cmap="viridis");fig.colorbar(im,ax=axes[0,0],label="capture probability");axes[0,0].set(xlabel="z",ylabel="y",title=f"zero torque, In=0.5, phi={zero['phi'][k]:.2f}")
    axes[0,1].plot([r["noise"] for r in noise_rows],[r["classification_error"] for r in noise_rows],"o-",label="PDE classification vs ray");axes[0,1].plot([r["noise"] for r in noise_rows],[r["short_flight_asymptotic_to_ray_l1"] for r in noise_rows],"s--",label="asymptotic probability vs ray");axes[0,1].set(xlabel="angular noise In",ylabel="error",title="zero-noise ray limit");axes[0,1].legend(frameon=False)
    ax=axes[1,0]
    for model,color in (("zero","tab:orange"),("frozen","tab:blue")):
        part=[r for r in calibration if r["model"]==model];ax.plot([r["mean_probability"] for r in part],[r["observed_capture"] for r in part],"o-",label=model,color=color)
    ax.plot([0,1],[0,1],"k:");ax.set(xlabel="predicted one-attempt capture",ylabel="observed",title="archived held-out approaches");ax.legend(frameon=False)
    ax=axes[1,1]
    for model,color in (("zero","tab:orange"),("frozen","tab:blue")):
        part=[r for r in short_validation if r["model"]==model and r["metric"]=="calibration_bias"]
        ax.plot([r["ballistic_time_max"] for r in part],[r["estimate"] for r in part],"o-",color=color,label=model)
        lo=[json.loads(r["ci95"])[0] for r in part];hi=[json.loads(r["ci95"])[1] for r in part];ax.fill_between([r["ballistic_time_max"] for r in part],lo,hi,color=color,alpha=.15)
    ax.axhline(0,color="k",ls=":");ax.set(xlabel="initial ballistic exit time cutoff",ylabel="predicted - observed",title="regime of one-attempt validity");ax.legend(frameon=False)
    fig.tight_layout();output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    sol=solve_backward((9,7,12),0.,0.);err=grid_ballistic_error(sol)
    assert err["classification_error"]<.2 and np.all((sol["q"]>=0)&(sol["q"]<=1))
    t,k=segment_exit(-.1,0,.1,0);assert k=="gap" and abs(t-.5)<1e-12
    z=np.asarray([-.5]);y=np.asarray([.1]);p=np.asarray([.2]);duration,kind=deterministic_arc_exit(z,y,p,1.)
    assert np.isfinite(duration[0]) and kind[0] in ("gap","lateral","upstream")
    print("PASS: monotone kinetic solve, ballistic and circular characteristics, and segment exit")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--risk-file",type=Path);p.add_argument("--production-source",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.risk_file,args.production_source,args.output)):p.error("all paths required")
    args.output.mkdir(parents=True,exist_ok=True)
    convergence,noise_rows,zero=convergence_audit(args.output);write_csv(args.output/"grid_convergence.csv",convergence);write_csv(args.output/"noise_asymptotics.csv",noise_rows)
    circular_rows,_=constant_turn_audit();write_csv(args.output/"constant_turn_validation.csv",circular_rows)
    omega_solutions=[]
    for omega in OMEGA_NODES:
        print(f"constant-turn solve omega={omega:g}",flush=True);omega_solutions.append(zero if omega==0 else solve_backward(FINAL_GRID,IN_NOISE,float(omega)))
    attempts=extract_attempts(args.raw_dir,args.risk_file,args.production_source);attach_predictions(attempts,zero,omega_solutions)
    validation,comparison=seed_bootstrap_validation(attempts);calibration=calibration_rows(attempts,"zero")+calibration_rows(attempts,"frozen")
    short_validation=short_flight_validation(attempts)
    write_csv(args.output/"trajectory_validation.csv",validation);write_csv(args.output/"calibration_bins.csv",calibration);write_csv(args.output/"short_flight_validation.csv",short_validation)
    make_figure(zero,noise_rows,calibration,validation,short_validation,args.output/"first_exit_solution.png")
    summary=dict(grid=dict(final_shape=FINAL_GRID,omega_nodes=OMEGA_NODES.tolist(),noise_levels=list(NOISE_LEVELS)),
                 attempts=dict(n=len(attempts),n_seeds=len({r['seed'] for r in attempts}),uncensored=sum(not r['censored'] for r in attempts),captured=sum(r['actual_capture'] for r in attempts),censored=sum(r['censored'] for r in attempts),omega_clip_fraction=float(np.mean([r['omega_clipped'] for r in attempts]))),
                 validation=validation,zero_minus_frozen=comparison,noise_asymptotics=noise_rows,constant_turn_validation=circular_rows,short_flight_validation=short_validation,
                 interpretation="one approach attempt; lateral/upstream exit is failure for this attempt but permits later renewal")
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),risk_file=str(args.risk_file),production_source=str(args.production_source),heldout_seeds="30--39",dt=DT,ell=ELL,half_width=HALF_WIDTH,angular_noise=IN_NOISE,attempt_horizon=5.,fitted_coefficients=0,new_microscopic_simulations=0)
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
