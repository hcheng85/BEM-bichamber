#!/usr/bin/env python3
"""Bounded validation of angular persistence and prospective phase opportunities.

No BEM solve, fitted coefficient, or new simulation is used.  Native saved
states supply exact positions/headings and the effective translational residual.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

DT = 1e-3
STRIDE = 10
SDT = DT*STRIDE
BURN = 200.0
GAP = 7.158910531638177
RADIUS = GAP
TIP_R = .5
REFERENCE_DISTANCES = (.5, 1.0, 2.0)
PRIMARY_D = 1.0
HYSTERESIS_Y = -TIP_R


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="") as f:
        w=csv.DictWriter(f, fields, lineterminator="\n"); w.writeheader(); w.writerows(rows)


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def corridor(x, y, theta, sign):
    """Exact rounded-tip corridor in chamber-reflected angular coordinates."""
    d = -sign*np.asarray(x)
    y = np.asarray(y); theta = np.asarray(theta)
    cm = -GAP/2-TIP_R; cp = GAP/2+TIP_R
    vm = cm-y; vp = cp-y
    rm = np.hypot(d, vm); rp = np.hypot(d, vp)
    domain = (d > 0) & (rm > TIP_R) & (rp > TIP_R)
    lo = np.full(np.broadcast(d, y).shape, np.nan, dtype=float)
    hi = np.full_like(lo, np.nan)
    lo[domain] = np.arctan2(vm[domain], d[domain])+np.arcsin(TIP_R/rm[domain])
    hi[domain] = np.arctan2(vp[domain], d[domain])-np.arcsin(TIP_R/rp[domain])
    width = hi-lo
    valid = domain & np.isfinite(width) & (width > 0)
    center = .5*(lo+hi); half=.5*width
    theta_local = np.arctan2(np.sin(theta), sign*np.cos(theta))
    delta = wrap(theta_local-center)
    margin = half-np.abs(delta)
    inside = valid & (margin >= 0)
    return dict(d=d, lo=lo, hi=hi, width=width, center=center, half=half,
                theta_local=theta_local, margin=margin, valid=valid, inside=inside)


def first_crossing_native(state, approximate_time, sign):
    center=int(round(approximate_time/DT)); lo=max(0,center-12); hi=min(state.shape[2],center+2)
    crossed=sign*state[0,:,lo:hi] >= 0
    hits=np.argwhere(crossed)
    if not len(hits):
        return None, None, 0
    rel=int(hits[:,1].min()); fish=np.flatnonzero(crossed[:,rel])
    # One representative event identity: most inward among simultaneous fish.
    localx=sign*state[0,fish,lo+rel]
    focal=int(fish[np.argmax(localx)])
    return lo+rel, focal, len(fish)


def residence_metrics(inside):
    inside=np.asarray(inside,bool)
    n=len(inside)
    if n == 0:
        return dict(residence_fraction=np.nan,total_residence=0.,initial_residence=0.,
                    final_residence=0.,exit_count=0,persistent=0)
    exits=int(np.sum(inside[:-1] & ~inside[1:]))
    first_false=np.flatnonzero(~inside)
    initial=(first_false[0] if inside[0] and len(first_false) else n if inside[0] else 0)*DT
    trailing=0
    for value in inside[::-1]:
        if not value: break
        trailing += 1
    return dict(residence_fraction=float(np.mean(inside)),total_residence=float(np.sum(inside)*DT),
                initial_residence=float(initial),final_residence=float(trailing*DT),
                exit_count=exits,persistent=int(np.all(inside)))


def crossing_rows(seed, state, episodes):
    rows=[]
    for event_id, ep in enumerate(episodes):
        origin=ep["origin"]; sign=1. if origin=="left" else -1.
        cross, focal, tied=first_crossing_native(state,float(ep["first_cross_time"]),sign)
        if cross is None:
            for ref in REFERENCE_DISTANCES:
                rows.append(dict(seed=seed,event_id=event_id,origin=origin,reference_distance=ref,
                                 status="crossing_not_resolved",first_cross_time=ep["first_cross_time"]))
            continue
        full_start=max(0,int(round((float(ep["first_cross_time"])-float(ep["full_dwell_before"]))/DT)))
        d=-sign*state[0,focal,full_start:cross+1]
        for ref in REFERENCE_DISTANCES:
            if d[0] <= ref:
                status="reference_left_censored"; refidx=None
            else:
                enter=np.flatnonzero((d[1:]<=ref)&(d[:-1]>ref))
                refidx=full_start+1+int(enter[0]) if len(enter) else None
                status="ok" if refidx is not None else "reference_not_reached"
            base=dict(seed=seed,event_id=event_id,origin=origin,first_cross_time=cross*DT,
                      focal_fish=focal,simultaneous_first_crossers=tied,
                      reference_distance=ref,status=status)
            if refidx is None:
                rows.append(base); continue
            sl=slice(refidx,cross)  # exclude d<=0 crossing state where the cone is undefined
            x=state[0,focal,sl]; y=state[1,focal,sl]; th=state[2,focal,sl]
            g=corridor(x,y,th,sign)
            if not bool(g["valid"][0]):
                base["status"]="invalid_reference_corridor"; rows.append(base); continue
            frozen_center=float(g["center"][0]); frozen_half=float(g["half"][0])
            frozen_inside=np.abs(wrap(g["theta_local"]-frozen_center)) <= frozen_half
            instant_inside=g["inside"]
            frozen=residence_metrics(frozen_inside)
            instant=residence_metrics(instant_inside)
            steps=cross-refidx
            self_disp=DT*np.sum(sign*np.cos(state[2,focal,refidx:cross]))
            total_disp=sign*(state[0,focal,cross]-state[0,focal,refidx])
            max_abs=max(abs(float(g["lo"][0])),abs(float(g["hi"][0])))
            safe_time=(float(g["d"][0])/np.cos(max_abs)) if max_abs<np.pi/2 else np.inf
            base.update(reference_time=refidx*DT,flight_time=steps*DT,
                        reference_d=float(g["d"][0]),reference_y=float(y[0]),
                        reference_width=float(g["width"][0]),
                        reference_margin=float(g["margin"][0]),
                        valid_corridor_fraction=float(np.mean(g["valid"])),
                        self_propulsion_displacement=float(self_disp),
                        translational_hydrodynamics=float(total_disp-self_disp),
                        total_inward_displacement=float(total_disp),
                        isolated_safe_flight_bound=float(safe_time),
                        **{f"frozen_{k}":v for k,v in frozen.items()},
                        **{f"instant_{k}":v for k,v in instant.items()})
            rows.append(base)
    return rows


def section_cycles(state):
    idx=np.arange(int(BURN/DT),state.shape[2]-1,STRIDE)
    x=state[0][:,idx]; y=state[1][:,idx]
    left=np.all(x<0,axis=0); right=np.all(x>=0,axis=0)
    full=np.where(left,1,np.where(right,-1,0))  # sign points inward
    yc=np.mean(y,axis=0)
    up=np.flatnonzero((full[1:]!=0)&(full[1:]==full[:-1])&(yc[:-1]<0)&(yc[1:]>=0))+1
    accepted=[]
    for q in up:
        if not accepted:
            accepted.append(int(q)); continue
        prev=accepted[-1]; prevsign=full[prev]
        reset=np.any((full[prev:q+1]==prevsign)&(yc[prev:q+1]<=HYSTERESIS_Y))
        if full[q]!=prevsign or reset:
            accepted.append(int(q))
    cycles=[]
    for ci,(a,b) in enumerate(zip(accepted[:-1],accepted[1:])):
        sign=float(full[a]); origin="left" if sign==1 else "right"
        local=sign*x[:,a:b+1]
        hit=np.argwhere(local>=0)
        if len(hit):
            rel=int(hit[:,1].min()); end=a+rel; outcome=1
            first_fish=np.flatnonzero(local[:,rel]>=0)
            first=int(first_fish[np.argmax(local[first_fish,rel])])
        else:
            end=b; outcome=0; first=-1
        cycles.append(dict(cycle_id=ci,start_sample=a,end_sample=end,return_sample=b,
                           start_index=int(idx[a]),end_index=int(idx[end]),return_index=int(idx[b]),
                           origin=origin,sign=sign,outcome_crossing=outcome,first_crosser=first,
                           section_y=float(yc[a]),cycle_duration=(idx[end]-idx[a])*DT,
                           complete_return=int(not outcome)))
    return cycles, idx


def opportunity_rows(seed,state,cycles):
    fish_rows=[]; cycle_rows=[]
    for cyc in cycles:
        sign=float(cyc["sign"]); start=int(cyc["start_index"]); end=int(cyc["end_index"])
        inds=np.arange(start,end+1,STRIDE)
        x=state[0][:,inds]; y=state[1][:,inds]; th=state[2][:,inds]
        d=-sign*x
        for ref in REFERENCE_DISTANCES:
            rows=[]
            for fish in range(state.shape[1]):
                g=corridor(x[fish],y[fish],th[fish],sign)
                eligible=(d[fish]<=ref)&(d[fish]>0)&g["valid"]
                hit=np.flatnonzero(eligible)
                if not len(hit): continue
                j=int(hit[0]); crossed=bool(np.any(sign*state[0,fish,inds[j]:end+1]>=0))
                row=dict(seed=seed,cycle_id=cyc["cycle_id"],origin=cyc["origin"],
                         reference_distance=ref,fish=fish,opportunity_time=inds[j]*DT,
                         time_from_section=(inds[j]-start)*DT,censored_inside_shell=int(d[fish,0]<=ref),
                         distance=float(d[fish,j]),height=float(y[fish,j]),
                         corridor_width=float(g["width"][j]),angular_margin=float(g["margin"][j]),
                         inside_corridor=int(g["inside"][j]),candidate_crossed=int(crossed),
                         cycle_crossed=cyc["outcome_crossing"],
                         candidate_was_first_crosser=int(cyc["outcome_crossing"] and fish==cyc["first_crosser"]))
                rows.append(row); fish_rows.append(row)
            eligible=[r for r in rows]
            if eligible:
                earliest=min(float(r["opportunity_time"]) for r in eligible)
                tied=[r for r in eligible if float(r["opportunity_time"])==earliest]
                lead=min(tied,key=lambda r:float(r["distance"]))
                cycle_rows.append(dict(**lead,n_fish_opportunities=len(eligible),
                                       n_distinct_opportunity_fish=len({int(r["fish"]) for r in eligible})))
            else:
                cycle_rows.append(dict(seed=seed,cycle_id=cyc["cycle_id"],origin=cyc["origin"],
                    reference_distance=ref,fish="",opportunity_time="",time_from_section="",
                    censored_inside_shell="",distance="",height="",corridor_width="",angular_margin="",
                    inside_corridor="",candidate_crossed="",cycle_crossed=cyc["outcome_crossing"],
                    candidate_was_first_crosser="",n_fish_opportunities=0,n_distinct_opportunity_fish=0))
    return fish_rows,cycle_rows


def auc(y,score):
    y=np.asarray(y,int); score=np.asarray(score,float)
    pos=score[y==1]; neg=score[y==0]
    if not len(pos) or not len(neg): return np.nan
    return float(np.mean(pos[:,None]>neg[None,:])+.5*np.mean(pos[:,None]==neg[None,:]))


def seed_boot(values,seed_ids,nboot=20000,rng_seed=20261024):
    values=np.asarray(values,float); seed_ids=np.asarray(seed_ids,int)
    seeds=np.unique(seed_ids); means=[]
    for s in seeds:
        q=values[seed_ids==s]
        if np.any(np.isfinite(q)): means.append(np.nanmean(q))
    means=np.asarray(means)
    if not len(means): return np.nan,np.nan,np.nan,0
    rng=np.random.default_rng(rng_seed)
    boot=np.mean(means[rng.integers(0,len(means),size=(nboot,len(means)))],axis=1)
    return float(np.mean(means)),float(np.quantile(boot,.025)),float(np.quantile(boot,.975)),len(means)


def summarize_crossings(rows):
    out=[]
    for ref in REFERENCE_DISTANCES:
        q=[r for r in rows if float(r["reference_distance"])==ref]
        ok=[r for r in q if r["status"]=="ok"]
        for metric in ("frozen_persistent","instant_persistent","frozen_residence_fraction",
                       "instant_residence_fraction","frozen_exit_count","instant_exit_count",
                       "flight_time","self_propulsion_displacement","translational_hydrodynamics",
                       "isolated_safe_flight_bound"):
            vals=[float(r[metric]) for r in ok]; seeds=[int(r["seed"]) for r in ok]
            mean,lo,hi,nseed=seed_boot(vals,seeds)
            out.append(dict(reference_distance=ref,metric=metric,n_events=len(ok),n_seeds=nseed,
                            seed_mean=mean,ci_low=lo,ci_high=hi))
        status={s:sum(r["status"]==s for r in q) for s in sorted({r["status"] for r in q})}
        for s,n in status.items(): out.append(dict(reference_distance=ref,metric=f"status_{s}",n_events=n,n_seeds="",seed_mean=n/len(q),ci_low="",ci_high=""))
    return out


def summarize_opportunities(rows):
    out=[]
    for ref in REFERENCE_DISTANCES:
        q=[r for r in rows if float(r["reference_distance"])==ref and int(r["n_fish_opportunities"])>0]
        seeds=sorted({int(r["seed"]) for r in q})
        seed_auc=[]; seed_diff=[]; seed_margin=[]; seed_width=[]
        auc_seeds=[]; diff_seeds=[]; contrast_seeds=[]
        for s in seeds:
            z=[r for r in q if int(r["seed"])==s]
            y=np.array([int(r["cycle_crossed"]) for r in z]); score=np.array([float(r["angular_margin"]) for r in z])
            a=auc(y,score)
            if np.isfinite(a): seed_auc.append(a); auc_seeds.append(s)
            inside=np.array([int(r["inside_corridor"]) for r in z])
            if np.any(inside==1) and np.any(inside==0):
                seed_diff.append(np.mean(y[inside==1])-np.mean(y[inside==0])); diff_seeds.append(s)
            if np.any(y==1) and np.any(y==0):
                seed_margin.append(np.mean(score[y==1])-np.mean(score[y==0]))
                widths=np.array([float(r["corridor_width"]) for r in z])
                seed_width.append(np.mean(widths[y==1])-np.mean(widths[y==0]))
                contrast_seeds.append(s)
        for metric,vals,sids in (("seed_auc_angular_margin",seed_auc,auc_seeds),
                                 ("crossing_rate_inside_minus_outside",seed_diff,diff_seeds),
                                 ("margin_cross_minus_noncross",seed_margin,contrast_seeds),
                                 ("width_cross_minus_noncross",seed_width,contrast_seeds)):
            mean,lo,hi,nseed=seed_boot(vals,sids) if vals else (np.nan,np.nan,np.nan,0)
            out.append(dict(reference_distance=ref,metric=metric,n_cycles=len(q),n_seeds=nseed,
                            seed_mean=mean,ci_low=lo,ci_high=hi))
        out.append(dict(reference_distance=ref,metric="cycles_with_opportunity",n_cycles=len(q),n_seeds=len(seeds),
                        seed_mean=len(q),ci_low="",ci_high=""))
        descriptive=(
            ("crossing_fraction",lambda r: int(r["cycle_crossed"]),None),
            ("inside_corridor_fraction",lambda r: int(r["inside_corridor"]),None),
            ("earliest_candidate_crossed_fraction",lambda r: int(r["candidate_crossed"]),None),
            ("mean_number_opportunity_fish",lambda r: int(r["n_fish_opportunities"]),None),
            ("earliest_candidate_is_first_crosser_given_crossing",
             lambda r: int(r["candidate_was_first_crosser"]),lambda r: int(r["cycle_crossed"])==1),
        )
        for metric,extract,keep in descriptive:
            seed_values=[]; value_seeds=[]
            for s in seeds:
                z=[r for r in q if int(r["seed"])==s and (keep is None or keep(r))]
                if z:
                    seed_values.append(np.mean([extract(r) for r in z]));value_seeds.append(s)
            mean,lo,hi,nseed=seed_boot(seed_values,value_seeds)
            out.append(dict(reference_distance=ref,metric=metric,n_cycles=len(q),n_seeds=nseed,
                            seed_mean=mean,ci_low=lo,ci_high=hi))
    return out


def seed_level_summary(crossings, opportunities):
    """Compact independent-seed table underlying all interval estimates."""
    out=[]
    for ref in REFERENCE_DISTANCES:
        for seed in range(20):
            c=[r for r in crossings if r["status"]=="ok" and int(r["seed"])==seed
               and float(r["reference_distance"])==ref]
            p=[r for r in opportunities if int(r["seed"])==seed
               and float(r["reference_distance"])==ref and int(r["n_fish_opportunities"])>0]
            row=dict(seed=seed,reference_distance=ref,n_analyzable_crossings=len(c),n_prospective_cycles=len(p))
            for metric in ("frozen_persistent","instant_persistent","frozen_residence_fraction",
                           "instant_residence_fraction","frozen_exit_count","instant_exit_count",
                           "flight_time","self_propulsion_displacement","translational_hydrodynamics",
                           "isolated_safe_flight_bound"):
                row[metric]=np.mean([float(r[metric]) for r in c]) if c else np.nan
            if p:
                y=np.array([int(r["cycle_crossed"]) for r in p]);m=np.array([float(r["angular_margin"]) for r in p])
                inside=np.array([int(r["inside_corridor"]) for r in p])
                row["auc_angular_margin"]=auc(y,m)
                row["crossing_fraction"]=np.mean(y)
                row["inside_corridor_fraction"]=np.mean(inside)
                row["crossing_rate_inside_minus_outside"]=(np.mean(y[inside==1])-np.mean(y[inside==0])
                    if np.any(inside==1) and np.any(inside==0) else np.nan)
                row["mean_number_opportunity_fish"]=np.mean([int(r["n_fish_opportunities"]) for r in p])
            out.append(row)
    return out


def figures(crossings,cycles,summary_cross,summary_opp,out):
    out.mkdir(parents=True,exist_ok=True)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    for j,(metric,title) in enumerate((("instant_residence_fraction","instantaneous corridor"),("frozen_residence_fraction","frozen corridor"))):
        data=[]; labels=[]
        for ref in REFERENCE_DISTANCES:
            vals=[float(r[metric]) for r in crossings if r["status"]=="ok" and float(r["reference_distance"])==ref]
            data.append(vals);labels.append({.5:"$d=r$",1.:"$d=2r$",2.:"$d=4r$"}[ref])
        axs[j].violinplot(data,showmedians=True);axs[j].set_xticks(range(1,4),labels)
        axs[j].set_ylim(-.03,1.03);axs[j].set_ylabel("fraction of reference-to-crossing time inside")
        axs[j].set_title(title);axs[j].grid(alpha=.2)
    fig.suptitle("Actual first crossers: angular residence from first shell entry")
    fig.tight_layout();fig.savefig(out/"01_actual_crossing_persistence.png",dpi=220);fig.savefig(out/"01_actual_crossing_persistence.pdf");plt.close(fig)

    q=[r for r in cycles if float(r["reference_distance"])==PRIMARY_D and int(r["n_fish_opportunities"])>0]
    fig,axs=plt.subplots(1,2,figsize=(11,4.5))
    for outcome,color,label in ((0,"#777777","no crossing"),(1,"#D55E00","crossing")):
        z=[r for r in q if int(r["cycle_crossed"])==outcome]
        axs[0].scatter([float(r["angular_margin"]) for r in z],[float(r["corridor_width"]) for r in z],s=18,alpha=.55,color=color,label=label)
    axs[0].axvline(0,color="k",lw=1);axs[0].set(xlabel="initial angular margin (positive = inside)",ylabel="corridor width",title="Earliest opportunity per cycle");axs[0].legend(frameon=False)
    for inside,color,label in ((0,"#999999","outside"),(1,"#009E73","inside")):
        seed_rates=[]
        for s in sorted({int(r["seed"]) for r in q}):
            z=[r for r in q if int(r["seed"])==s and int(r["inside_corridor"])==inside]
            if z: seed_rates.append(np.mean([int(r["cycle_crossed"]) for r in z]))
        axs[1].scatter(np.full(len(seed_rates),inside)+np.linspace(-.06,.06,max(len(seed_rates),1)),seed_rates,color=color,alpha=.7)
        if seed_rates: axs[1].plot([inside-.15,inside+.15],[np.mean(seed_rates)]*2,color="k",lw=2)
    axs[1].set(xticks=[0,1],xticklabels=["outside","inside"],ylabel="per-seed crossing fraction",title="Prospective spatial outcome",ylim=(-.05,1.05));axs[1].grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/"02_prospective_opportunities.png",dpi=220);fig.savefig(out/"02_prospective_opportunities.pdf");plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--raw-dir",type=Path,required=True)
    ap.add_argument("--episodes",type=Path,required=True);ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--seeds",type=int,nargs="*",default=list(range(20)));ap.add_argument("--force",action="store_true")
    args=ap.parse_args();out=args.out.resolve();cache=out/"tables"/"by_seed";cache.mkdir(parents=True,exist_ok=True)
    episodes=read_csv(args.episodes);all_cross=[];all_fish=[];all_cycles=[];all_cycle_defs=[]
    for seed in args.seeds:
        cf=cache/f"seed_{seed:02d}_crossings.csv";ff=cache/f"seed_{seed:02d}_fish_opportunities.csv";yf=cache/f"seed_{seed:02d}_cycle_opportunities.csv";df=cache/f"seed_{seed:02d}_cycles.csv"
        if cf.exists() and ff.exists() and yf.exists() and df.exists() and not args.force:
            all_cross+=read_csv(cf);all_fish+=read_csv(ff);all_cycles+=read_csv(yf);all_cycle_defs+=read_csv(df);print(f"seed {seed}: checkpoint",flush=True);continue
        raw=args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
        cross=crossing_rows(seed,state,[r for r in episodes if int(r["seed"])==seed])
        cycles,_=section_cycles(state);fish,cycle=opportunity_rows(seed,state,cycles)
        for r in cycles:r["seed"]=seed
        write_csv(cf,cross);write_csv(ff,fish);write_csv(yf,cycle);write_csv(df,cycles)
        all_cross+=cross;all_fish+=fish;all_cycles+=cycle;all_cycle_defs+=cycles
        del state;gc.collect();print(f"seed {seed}: done; crossings {len(cross)//3}, cycles {len(cycles)}",flush=True)
    write_csv(out/"tables"/"actual_crossing_persistence.csv",all_cross)
    write_csv(out/"tables"/"prospective_fish_opportunities.csv",all_fish)
    write_csv(out/"tables"/"prospective_cycle_opportunities.csv",all_cycles)
    write_csv(out/"tables"/"prospective_cycles.csv",all_cycle_defs)
    sc=summarize_crossings(all_cross);so=summarize_opportunities(all_cycles)
    write_csv(out/"tables"/"crossing_persistence_summary.csv",sc);write_csv(out/"tables"/"prospective_summary.csv",so)
    write_csv(out/"tables"/"seed_level_summary.csv",seed_level_summary(all_cross,all_cycles))
    figures(all_cross,all_cycles,sc,so,out/"figures")
    config=dict(native_dt=DT,prospective_stride=STRIDE,section="full-occupancy upward Yc=0",
        hysteresis=f"same-origin Yc <= {HYSTERESIS_Y} before return",reference_distances=REFERENCE_DISTANCES,
        primary_reference=PRIMARY_D,tip_radius=TIP_R,gap=GAP,R=RADIUS,gap_over_R=1,gap_over_2R=.5,
        new_simulations=False,fitted_coefficients=False,seeds=args.seeds)
    (out/"analysis_config.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__":main()
