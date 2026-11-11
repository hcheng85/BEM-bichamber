#!/usr/bin/env python3
"""Conservative longitudinal-current correction for no-hydro gap-20 data.

All transport terms are measured from straight native-step particle segments.
The production simulator and raw trajectories are read-only inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from analyze_approach_current_balance import DT, ELL, N, RawReader, indicator, raw_path, read_csv, write_csv

HALF_WIDTH = ELL / 2
WINDOW_STEPS = 200
Q_EDGES = np.asarray([0.0, 0.25, 0.50, 0.75, 0.90, 1.000001])
JA_EDGES = np.asarray([-0.15, -0.05, 0.0, 0.025, 0.05, 0.10, 0.20, 1.0])


def online_direction(state: np.ndarray) -> tuple[int, float, str]:
    """Direction from information at this frame; never uses a future event."""
    q = float(1 - 2 * np.sum(state[0] < 0) / N)
    if q < 0:
        return 1, q, "majority_left"
    if q > 0:
        return -1, q, "majority_right"
    px = float(np.mean(np.cos(state[2])))
    return (1 if px >= 0 else -1), q, "population_tie_heading"


def segment_budget(states: np.ndarray, velocities: np.ndarray, direction: int) -> dict[str, float]:
    """Exact piecewise-linear finite-window conservation budget.

    Local coordinate z=s*x lies in [-ELL,0] in the source-side rectangle.
    Positive boundary counts are outward from that rectangle.
    """
    s = float(direction)
    z = s * states[0]
    y = states[1]
    z0, z1 = z[:, :-1], z[:, 1:]
    y0, y1 = y[:, :-1], y[:, 1:]
    dz, dy = z1 - z0, y1 - y0
    nsteps = z0.shape[1]
    duration = nsteps * DT

    # Liang-Barsky interval of every native straight segment inside rectangle.
    enter = np.zeros_like(z0)
    leave = np.ones_like(z0)
    valid = np.ones_like(z0, dtype=bool)
    for coord0, delta, lower, upper in ((z0, dz, -ELL, 0.0), (y0, dy, -HALF_WIDTH, HALF_WIDTH)):
        moving = np.abs(delta) > 0
        lo = np.zeros_like(coord0); hi = np.ones_like(coord0)
        a = np.zeros_like(coord0); b = np.zeros_like(coord0)
        a[moving] = (lower - coord0[moving]) / delta[moving]
        b[moving] = (upper - coord0[moving]) / delta[moving]
        lo[moving] = np.minimum(a[moving], b[moving])
        hi[moving] = np.maximum(a[moving], b[moving])
        valid &= moving | ((coord0 >= lower) & (coord0 <= upper))
        enter = np.maximum(enter, lo)
        leave = np.minimum(leave, hi)
    valid &= leave > enter
    path_current_integral = float(np.sum(dz[valid] * (leave[valid] - enter[valid])))

    # Neck crossing at z=0, signed outward.  The y check prevents divider-wall
    # contacts outside the opening from being counted as neck transport.
    neck_mask = ((z0 < 0) & (z1 >= 0)) | ((z0 >= 0) & (z1 < 0))
    t_neck = np.zeros_like(z0)
    nonzero_z = np.abs(dz) > 0
    t_neck[nonzero_z] = -z0[nonzero_z] / dz[nonzero_z]
    y_neck = y0 + t_neck * dy
    neck_mask &= nonzero_z & (t_neck >= 0) & (t_neck <= 1) & (np.abs(y_neck) <= HALF_WIDTH)
    neck_count = float(np.sum(np.sign(dz[neck_mask])))

    lateral_weighted = 0.0
    lateral_count = 0.0
    # Top: upward is outward.  Bottom: downward is outward.
    for boundary, outward_sign in ((HALF_WIDTH, 1.0), (-HALF_WIDTH, -1.0)):
        if boundary > 0:
            cross = ((y0 <= boundary) & (y1 > boundary)) | ((y0 > boundary) & (y1 <= boundary))
        else:
            cross = ((y0 >= boundary) & (y1 < boundary)) | ((y0 < boundary) & (y1 >= boundary))
        nonzero_y = np.abs(dy) > 0
        t = np.zeros_like(y0); t[nonzero_y] = (boundary - y0[nonzero_y]) / dy[nonzero_y]
        z_cross = z0 + t * dz
        cross &= nonzero_y & (t >= 0) & (t <= 1) & (z_cross >= -ELL) & (z_cross < 0)
        signs = outward_sign * np.sign(dy[cross])
        lateral_count += float(np.sum(signs))
        lateral_weighted += float(np.sum(signs * (z_cross[cross] + ELL)))

    inside0 = (z[:, 0] < 0) & (z[:, 0] >= -ELL) & (np.abs(y[:, 0]) <= HALF_WIDTH)
    inside1 = (z[:, -1] < 0) & (z[:, -1] >= -ELL) & (np.abs(y[:, -1]) <= HALF_WIDTH)
    moment0 = float(np.sum((z[:, 0] + ELL) * inside0))
    moment1 = float(np.sum((z[:, -1] + ELL) * inside1))
    moment_change = moment1 - moment0

    # Instantaneous stored-current approximation used in the proposed plug law.
    ja = []
    for k in range(nsteps):
        inside = (z[:, k] < 0) & (z[:, k] >= -ELL) & (np.abs(y[:, k]) <= HALF_WIDTH)
        ja.append(float(np.sum(s * velocities[0, :, k] * inside) / N))
    mean_ja = float(np.mean(ja))
    plug_flux = N * mean_ja / ELL
    path_flux = path_current_integral / (ELL * duration)
    accumulation = moment_change / (ELL * duration)
    lateral = lateral_weighted / (ELL * duration)
    path_correction = path_flux - plug_flux
    neck_flux = neck_count / duration
    reconstructed = plug_flux + path_correction - accumulation - lateral
    closure = neck_flux - reconstructed
    total_correction = accumulation + lateral - path_correction
    return dict(duration=duration, q_start=float(1 - 2 * np.sum(states[0, :, 0] < 0) / N),
                Ja=mean_ja, n_approach_start=int(np.sum(inside0)), n_approach_end=int(np.sum(inside1)),
                plug_flux=plug_flux, path_flux=path_flux, neck_flux=neck_flux,
                weighted_accumulation=accumulation, weighted_lateral_flux=lateral,
                path_discretization_correction=path_correction, total_correction=total_correction,
                reconstructed_neck_flux=reconstructed, closure_residual=closure,
                net_neck_count=neck_count, net_lateral_count=lateral_count,
                weighted_lateral_count=lateral_weighted, weighted_moment_change=moment_change)


def choose_first_by_seed(rows: list[dict[str, str]], predicate=lambda r: True) -> list[dict[str, str]]:
    out=[]
    for seed in sorted({int(r["seed"]) for r in rows if predicate(r)}):
        options=[r for r in rows if int(r["seed"])==seed and predicate(r)]
        out.append(min(options,key=lambda r:(float(r["time"]),int(r["event_id"]))))
    return out


def build_intervals(events: list[dict[str, str]], matched: list[dict[str, str]]) -> list[dict[str, Any]]:
    intervals=[]
    for ev in choose_first_by_seed(events):
        first=int(ev["first_native"])
        intervals.append(dict(phase="pre_crossing", trajectory_seed=int(ev["seed"]), anchor_seed=int(ev["seed"]),
                              event_id=int(ev["event_id"]), start=first-2001, stop=first-1,
                              eventual_direction=int(ev["direction"])))
    for ev in choose_first_by_seed(events,lambda r:int(r["continuation"])==1 and float(r["duration"])>=1.0):
        first=int(ev["first_native"])
        intervals.append(dict(phase="sustained_crossing", trajectory_seed=int(ev["seed"]),anchor_seed=int(ev["seed"]),
                              event_id=int(ev["event_id"]),start=first+200,stop=first+1000,
                              eventual_direction=int(ev["direction"])))
    for ev in choose_first_by_seed(events,lambda r:int(r["reversal"])==1 or float(r["duration"])>=0.5):
        end=int(round((float(ev["time"])+float(ev["duration"]))/DT))
        intervals.append(dict(phase="termination_or_reversal",trajectory_seed=int(ev["seed"]),anchor_seed=int(ev["seed"]),
                              event_id=int(ev["event_id"]),start=end-500,stop=end+500,
                              eventual_direction=int(ev["direction"])))

    # One lead-2 control per initiating seed, matching the previous analysis.
    event_lookup={(int(r["seed"]),int(r["event_id"])):r for r in events}
    lead2=[r for r in matched if math.isclose(float(r["lead"]),2.0)]
    for anchor in sorted({int(r["event_seed"]) for r in lead2}):
        options=[r for r in lead2 if int(r["event_seed"])==anchor]
        chosen=min(options,key=lambda r:(float(event_lookup[(anchor,int(r["event_id"]))]["time"]),int(r["event_id"])))
        start=int(round(float(chosen["control_time"])/DT))
        intervals.append(dict(phase="matched_non_crossing",trajectory_seed=int(chosen["control_seed"]),anchor_seed=anchor,
                              event_id=int(chosen["event_id"]),start=start,stop=start+2000,
                              eventual_direction=int(chosen["direction"])))
    return intervals


def extract_windows(raw_dir: Path, intervals: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows=[]
    for idx,meta in enumerate(intervals):
        if idx % 25 == 0: print(f"conservative intervals {idx}/{len(intervals)}",flush=True)
        reader=RawReader(raw_path(raw_dir,meta["trajectory_seed"]))
        states=reader.block("state",meta["start"],meta["stop"])
        velocities=reader.block("rdot",meta["start"]+1,meta["stop"]+1);reader.close()
        for offset in range(0,meta["stop"]-meta["start"],WINDOW_STEPS):
            if offset+WINDOW_STEPS>=states.shape[2]: break
            st=states[:,:,offset:offset+WINDOW_STEPS+1]
            vv=velocities[:,:,offset:offset+WINDOW_STEPS]
            direction,q,rule=online_direction(st[:,:,0])
            budget=segment_budget(st,vv,direction)
            rows.append(dict(phase=meta["phase"],trajectory_seed=meta["trajectory_seed"],anchor_seed=meta["anchor_seed"],
                             event_id=meta["event_id"],time=(meta["start"]+offset)*DT,direction=direction,
                             direction_rule=rule,eventual_direction=meta["eventual_direction"],
                             agrees_with_event_direction=int(direction==meta["eventual_direction"]),**budget))
    return rows


def bootstrap_seed_mean(rows: list[dict[str, Any]], key: str, nboot: int=5000) -> tuple[float,float,float]:
    seed_values=[]
    for seed in sorted({int(r["trajectory_seed"]) for r in rows}):
        seed_values.append(float(np.mean([float(r[key]) for r in rows if int(r["trajectory_seed"])==seed])))
    values=np.asarray(seed_values);rng=np.random.default_rng(20261112)
    draws=rng.choice(values,(nboot,len(values)),replace=True).mean(axis=1)
    return float(values.mean()),*map(float,np.quantile(draws,[.025,.975]))


def phase_summary(rows: list[dict[str, Any]], out: Path) -> tuple[list[dict[str, Any]],dict[str,Any]]:
    table=[];summary={}
    keys=("neck_flux","plug_flux","total_correction","weighted_accumulation","weighted_lateral_flux","path_discretization_correction")
    for phase in ("pre_crossing","sustained_crossing","termination_or_reversal","matched_non_crossing"):
        part=[r for r in rows if r["phase"]==phase]
        info={"n_windows":len(part),"n_trajectory_seeds":len({r['trajectory_seed'] for r in part})}
        for key in keys:
            est,lo,hi=bootstrap_seed_mean(part,key);info[key]={"estimate":est,"ci95":[lo,hi]}
            table.append(dict(phase=phase,quantity=key,estimate=est,ci95_low=lo,ci95_high=hi,
                              n_windows=len(part),n_seeds=info["n_trajectory_seeds"]))
        for key in ("weighted_accumulation","weighted_lateral_flux"):
            seed_means=[np.mean([r[key] for r in part if int(r["trajectory_seed"])==seed]) for seed in sorted({int(r["trajectory_seed"]) for r in part})]
            info[key]["seed_fraction_positive"]=float(np.mean(np.asarray(seed_means)>0))
            info[key]["seed_mean_absolute"]=float(np.mean(np.abs(seed_means)))
        seed_order=sorted({int(r["trajectory_seed"]) for r in part})
        seed_acc=np.asarray([np.mean([r["weighted_accumulation"] for r in part if int(r["trajectory_seed"])==seed]) for seed in seed_order])
        seed_lat=np.asarray([np.mean([r["weighted_lateral_flux"] for r in part if int(r["trajectory_seed"])==seed]) for seed in seed_order])
        info["seed_correlation_accumulation_lateral"]=float(np.corrcoef(seed_acc,seed_lat)[0,1])
        info["online_event_direction_agreement"]=float(np.mean([r["agrees_with_event_direction"] for r in part]))
        summary[phase]=info
    write_csv(out/"phase_contribution_summary.csv",table)
    return table,summary


def conditional_closure(rows: list[dict[str,Any]], out: Path) -> tuple[list[dict[str,Any]],dict[str,Any]]:
    # Symmetry reduction uses |q| only; no information beyond q is introduced.
    for r in rows:
        r["q_source"]=abs(float(r["q_start"]))
        r["q_bin"]=int(np.clip(np.searchsorted(Q_EDGES,r["q_source"],side="right")-1,0,len(Q_EDGES)-2))
        r["Ja_bin"]=int(np.clip(np.searchsorted(JA_EDGES,float(r["Ja"]),side="right")-1,0,len(JA_EDGES)-2))
    train=[r for r in rows if int(r["trajectory_seed"])<30]
    test=[r for r in rows if int(r["trajectory_seed"])>=30]
    cells=[];lookup={};component_keys=("total_correction","weighted_accumulation","weighted_lateral_flux","path_discretization_correction")
    for qb in range(len(Q_EDGES)-1):
        for jb in range(len(JA_EDGES)-1):
            tr=[r for r in train if r["q_bin"]==qb and r["Ja_bin"]==jb]
            te=[r for r in test if r["q_bin"]==qb and r["Ja_bin"]==jb]
            nseed=len({r["trajectory_seed"] for r in tr})
            means={key:float(np.mean([r[key] for r in tr])) if tr else float("nan") for key in component_keys}
            supported=nseed>=3 and bool(te)
            if supported: lookup[(qb,jb)]=means
            row=dict(q_bin=qb,Ja_bin=jb,q_low=Q_EDGES[qb],q_high=Q_EDGES[qb+1],
                              Ja_low=JA_EDGES[jb],Ja_high=JA_EDGES[jb+1],train_n=len(tr),train_n_seeds=nseed,
                              test_n=len(te),test_n_seeds=len({r['trajectory_seed'] for r in te}),supported=int(supported))
            for key in component_keys:
                row[f"train_mean_{key}"]=means[key]
                row[f"test_mean_{key}"]=float(np.mean([r[key] for r in te])) if te else float("nan")
            cells.append(row)
    write_csv(out/"conditional_correction_cells.csv",cells)
    evaluated=[]
    for r in test:
        key=(r["q_bin"],r["Ja_bin"])
        if key not in lookup: continue
        predicted=lookup[key];predicted_correction=predicted["total_correction"]
        evaluated.append({**r,"predicted_correction":predicted_correction,
                          **{f"predicted_{name}":value for name,value in predicted.items()},
                          "conditional_flux_prediction":float(r["plug_flux"])-predicted_correction,
                          "plug_abs_error":abs(float(r["plug_flux"])-float(r["neck_flux"])),
                          "conditional_abs_error":abs(float(r["plug_flux"])-predicted_correction-float(r["neck_flux"]))})
    write_csv(out/"heldout_conditional_closure.csv",evaluated)
    seed_rows=[]
    for seed in sorted({int(r["trajectory_seed"]) for r in evaluated}):
        p=[r for r in evaluated if int(r["trajectory_seed"])==seed]
        actual=np.asarray([r["neck_flux"] for r in p]);base=np.asarray([r["plug_flux"] for r in p]);corr=np.asarray([r["conditional_flux_prediction"] for r in p])
        seed_rows.append(dict(seed=seed,n=len(p),plug_mae=float(np.mean(np.abs(base-actual))),
                              conditional_mae=float(np.mean(np.abs(corr-actual))),
                              mae_improvement=float(np.mean(np.abs(base-actual)-np.abs(corr-actual))),
                              plug_bias=float(np.mean(base-actual)),conditional_bias=float(np.mean(corr-actual)),
                              plug_spearman=float(spearmanr(base,actual).statistic) if np.std(actual)>0 else float("nan"),
                              conditional_spearman=float(spearmanr(corr,actual).statistic) if np.std(actual)>0 else float("nan")))
    write_csv(out/"heldout_conditional_seed_summary.csv",seed_rows)
    imp=bootstrap_seed_mean([{**r,"trajectory_seed":r["seed"]} for r in seed_rows],"mae_improvement")
    supported_cells=[r for r in cells if r["supported"]]
    component_transfer={}
    for key in component_keys:
        train_means=np.asarray([r[f"train_mean_{key}"] for r in supported_cells]);test_means=np.asarray([r[f"test_mean_{key}"] for r in supported_cells])
        residual=np.asarray([r[key]-r[f"predicted_{key}"] for r in evaluated])
        component_transfer[key]=dict(cell_mean_spearman=float(spearmanr(train_means,test_means).statistic) if len(supported_cells)>1 else None,
                                     heldout_residual_rms=float(np.sqrt(np.mean(residual**2))),heldout_residual_mean=float(np.mean(residual)))
    summary=dict(q_edges=Q_EDGES.tolist(),Ja_edges=JA_EDGES.tolist(),minimum_train_seeds_per_cell=3,
                 n_train=len(train),n_test=len(test),n_evaluated=len(evaluated),heldout_coverage=len(evaluated)/len(test),
                 n_supported_cells=len(supported_cells),cell_mean_spearman=component_transfer["total_correction"]["cell_mean_spearman"],
                 component_transfer=component_transfer,
                 seed_mean_mae_improvement=imp[0],seed_mean_mae_improvement_ci95=[imp[1],imp[2]],
                 plug_seed_mean_mae=float(np.mean([r["plug_mae"] for r in seed_rows])),
                 conditional_seed_mean_mae=float(np.mean([r["conditional_mae"] for r in seed_rows])),
                 plug_seed_mean_bias=float(np.mean([r["plug_bias"] for r in seed_rows])),
                 conditional_seed_mean_bias=float(np.mean([r["conditional_bias"] for r in seed_rows])),
                 plug_seed_median_spearman=float(np.nanmedian([r["plug_spearman"] for r in seed_rows])),
                 conditional_seed_median_spearman=float(np.nanmedian([r["conditional_spearman"] for r in seed_rows])))
    (out/"conditional_closure_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return evaluated,summary


def make_figure(rows: list[dict[str,Any]], phase_table: list[dict[str,Any]], heldout: list[dict[str,Any]], out: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(11,8.5));axes=axes.ravel()
    residual=np.abs([r["closure_residual"] for r in rows]);axes[0].hist(np.maximum(residual,1e-18),bins=np.logspace(-18,-10,25));axes[0].set(xscale="log",xlabel="absolute conservative closure residual",ylabel="windows")
    phases=("pre_crossing","sustained_crossing","termination_or_reversal","matched_non_crossing")
    labels=("pre","sustained","termination/reversal","control");keys=("weighted_accumulation","weighted_lateral_flux","total_correction")
    x=np.arange(len(phases));colors=("tab:blue","tab:orange","tab:green")
    for j,(key,color) in enumerate(zip(keys,colors)):
        values=np.array([next(r["estimate"] for r in phase_table if r["phase"]==p and r["quantity"]==key) for p in phases])
        axes[1].bar(x+(j-1)*.24,values,width=.23,label=key.replace("weighted_","").replace("_flux",""),color=color)
    axes[1].axhline(0,color="k",lw=.8);axes[1].set_xticks(x,labels,rotation=20);axes[1].set_ylabel("mean contribution to plug discrepancy");axes[1].legend(frameon=False)
    actual=np.asarray([r["neck_flux"] for r in heldout]);base=np.asarray([r["plug_flux"] for r in heldout]);corrected=np.asarray([r["conditional_flux_prediction"] for r in heldout])
    lim=max(np.quantile(np.abs(np.r_[actual,base,corrected]),.98),1)
    axes[2].scatter(base,actual,s=12,alpha=.3,label="plug");axes[2].scatter(corrected,actual,s=12,alpha=.3,label="conditional correction");axes[2].plot([-lim,lim],[-lim,lim],"k:");axes[2].set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel="held-out predicted neck flux",ylabel="observed neck flux");axes[2].legend(frameon=False)
    axes[3].scatter([r["weighted_accumulation"] for r in rows],[r["plug_flux"]-r["neck_flux"] for r in rows],s=10,alpha=.25)
    axes[3].set(xlabel="weighted population accumulation",ylabel="plug-flow discrepancy")
    fig.tight_layout();fig.savefig(out/"figures"/"conservation_correction.png",dpi=220);plt.close(fig)


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path,required=True);p.add_argument("--previous",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);(args.output/"figures").mkdir(exist_ok=True)
    events=read_csv(args.previous/"crossing_episodes.csv");matched=read_csv(args.previous/"matched_precursor_pairs.csv")
    intervals=build_intervals(events,matched);write_csv(args.output/"analysis_intervals.csv",intervals)
    rows=extract_windows(args.raw_dir,intervals);write_csv(args.output/"conservative_window_budgets.csv",rows)
    phase_table,phase_info=phase_summary(rows,args.output)
    heldout,closure_info=conditional_closure(rows,args.output)
    make_figure(rows,phase_table,heldout,args.output)
    summary=dict(n_intervals=len(intervals),n_windows=len(rows),n_raw_seeds=len({r['trajectory_seed'] for r in rows}),
                 direction_rule="majority chamber toward gap; instantaneous mean heading only when q=0",
                 maximum_absolute_conservative_closure=float(np.max(np.abs([r["closure_residual"] for r in rows]))),
                 rms_conservative_closure=float(np.sqrt(np.mean(np.asarray([r["closure_residual"] for r in rows])**2))),
                 phase_summary=phase_info,conditional_closure=closure_info)
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),previous=str(args.previous),N=N,dt=DT,ell=ELL,width=ELL,window_duration=WINDOW_STEPS*DT,
                q_bins=Q_EDGES.tolist(),Ja_bins=JA_EDGES.tolist(),train_seeds="0--29",heldout_seeds="30--39",
                direction_assignment="current majority chamber; current mean heading tie-break; never eventual direction")
    (args.output/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")


if __name__=="__main__": main()
