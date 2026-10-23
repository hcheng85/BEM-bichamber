#!/usr/bin/env python3
"""Fixed-identity leading-fish angular first-passage analysis.

Uses only archived trajectories.  The exact RNG stream supplies the stochastic
angular increments; the saved heading increments then give the exact applied
deterministic impulse, including simulator clipping and overrides.  Individual
source terms are replayed from the simulator equations on a 0.005 time grid.
"""
from __future__ import annotations

import argparse, csv, gzip, json, sys, time
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.spatial import Delaunay

HERE = Path(__file__).resolve().parent
REPLAY_DIR = HERE.parents[1] / "pre_crossing_turning" / "scripts"
sys.path.insert(0, str(REPLAY_DIR))
from replay_dynamics import (  # noqa: E402
    Parameters, analytic_geometry, boundary_hydrodynamics, direct_hydrodynamics,
    hits_circle, intersects, ray_wall_intersection, replay_snapshot,
)

DT = 0.001
REPLAY_STRIDE = 5
N_LEADING = 10
SEEDS = tuple(range(20))
RNG_DISCONTINUITY = {10: 600000, 11: 600000, 13: 600000,
                     15: 600000, 17: 600000, 19: 600000}


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    opener=gzip.open if str(path).endswith(".gz") else open
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields: fields.append(key)
    with opener(path, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def angular_geometry(state, origin, gap, tip_r=0.5):
    """Exact tangent interval to the two rounded divider tips in local coordinates."""
    side = 1.0 if origin == "left" else -1.0
    d = -side * state[0]
    # Reflect the chamber-normal component only: (hx,hy)->(side*hx,hy).
    # Thus right-chamber angular rates acquire a factor -1.
    heading = np.arctan2(np.sin(state[2]), side*np.cos(state[2]))
    cy_lo, cy_hi = -gap/2-tip_r, gap/2+tip_r
    a_lo = np.arctan2(cy_lo-state[1], d)
    a_hi = np.arctan2(cy_hi-state[1], d)
    b_lo = np.arcsin(np.clip(tip_r/np.hypot(d, cy_lo-state[1]), 0, 1))
    b_hi = np.arcsin(np.clip(tip_r/np.hypot(d, cy_hi-state[1]), 0, 1))
    lo, hi = a_lo+b_lo, a_hi-b_hi
    inside = (heading >= lo) & (heading <= hi)
    dlo, dhi = wrap(heading-lo), wrap(heading-hi)
    eta = np.where(inside, 0.0, np.where(np.abs(dlo) <= np.abs(dhi), dlo, dhi))
    return side, d, heading, lo, hi, eta, inside


def social_selected(state, selected, geom, p):
    """Exact GetVision alignment/attraction, evaluated only for selected targets."""
    loc, angle = state[:2], state[2]
    tri = Delaunay(loc.T)
    neigh = [set() for _ in range(loc.shape[1])]
    for a, b, c in tri.simplices:
        neigh[a].update((b, c)); neigh[b].update((a, c)); neigh[c].update((a, b))
    align, attract = {}, {}
    for i in selected:
        visible = []
        for j in neigh[i]:
            a, b = loc[:, i], loc[:, j]
            blocked = any(intersects(a, b, c, d) for c, d in geom["segments"])
            blocked |= hits_circle(a, b, geom["c6"], p.tip_r)
            blocked |= hits_circle(a, b, geom["c5"], p.tip_r)
            if not blocked:
                visible.append(j)
        if not visible:
            align[i] = attract[i] = 0.0
            continue
        js = np.asarray(visible, int)
        rho = loc[:, js]-loc[:, i, None]
        dist = np.linalg.norm(rho, axis=0)
        ex = np.array([np.cos(angle[i]), np.sin(angle[i])])
        dot = np.clip(np.sum(ex[:, None]*rho, axis=0)/dist, -1, 1)
        sign = np.sign(ex[0]*rho[1]-ex[1]*rho[0]); sign[sign == 0] = 1
        theta = sign*np.arccos(dot)
        visual = 1+np.cos(theta); den = visual.sum()
        align[i] = float(np.sum(p.Ip*np.sin(angle[js]-angle[i])*visual)/den)
        attract[i] = float(np.sum(dist*np.sin(theta)*visual)/den)
    return align, attract


def direct_selected(state, selected, p):
    z, angle = state[0]+1j*state[1], state[2]
    out = {}
    for i in selected:
        mask = np.arange(len(angle)) != i
        dz, aj = z[i]-z[mask], angle[mask]
        out[i] = float(np.sum(np.imag(np.exp(1j*(2*angle[i]+aj))*dz /
                            (dz**2+p.delta**2)**2))*2*p.If/np.pi)
    return out


def wall_selected(state, selected, geom, p):
    out = {}
    for i in selected:
        normal, d = ray_wall_intersection(state[:, i], geom, p)
        phi = wrap(state[2, i]-normal)
        val = p.Iw*np.sign(phi)/d
        if abs(d) <= abs(p.safedis):
            val *= 1+1/d**(p.high_order-1)
        out[i] = float(np.clip(val, -np.pi/p.dt, np.pi/p.dt))
    return out


def rng_noise_for_indices(seed, indices, nfish=100, chunk=20000):
    """Replay simulator RNG and retain Gaussian rows at requested native steps."""
    indices = np.asarray(sorted(set(indices)), int)
    rs = np.random.RandomState(seed)
    alpha = rs.uniform(-np.pi, np.pi, nfish)
    rpos = np.sqrt(rs.uniform(0.2, 1, nfish))
    theta = rs.uniform(-1, 1, nfish)
    out = np.empty((len(indices), nfish), np.float64)
    cursor = 0; take = 0; maximum = int(indices[-1])
    while cursor <= maximum:
        stop = min(maximum+1, cursor+chunk)
        z = rs.normal(0, 1, (stop-cursor, nfish))
        while take < len(indices) and indices[take] < stop:
            out[take] = z[indices[take]-cursor]; take += 1
        cursor = stop
    return indices, out, alpha, rpos, theta


def first_crossing(state, i0, i1, origin):
    side = 1.0 if origin == "left" else -1.0
    localx = side*state[0, :, i0:i1+1]
    crossed = localx >= 0
    hits = np.argwhere(crossed)
    if not len(hits):
        return None, [], None
    rel = int(hits[:, 1].min())
    fish = np.flatnonzero(crossed[:, rel]).tolist()
    return i0+rel, fish, rel*DT


def build_segments(rows, state, seed, p):
    segments=[]; fish_rows=[]; history=[]
    for si, r in enumerate(rows):
        i0=int(round(float(r["risk_time"])/DT)); i1=i0+int(round(2/DT))
        origin=r["origin"]; side=1.0 if origin=="left" else -1.0
        d0=-side*state[0,:,i0]; fixed=np.argsort(d0)[:N_LEADING]
        search_end=min(state.shape[2]-1, i1+20)
        cross_idx,cross_ids,cross_rel=first_crossing(state,i0,search_end,origin)
        if r["condition"]=="control" and cross_idx is not None and cross_idx<=i1:
            raise RuntimeError(f"control crosses inside horizon seed={seed} pair={r['pair_id']}")
        # Event A is the first crossing at the end of the prospective horizon.
        event_cross_ids=cross_ids if r["condition"]=="event" else []
        covered=bool(set(event_cross_ids)&set(fixed))
        seg=dict(index=si,seed=seed,pair_id=int(r["pair_id"]),condition=r["condition"],
                 origin=origin,side=side,i0=i0,i1=i1,fixed=fixed,
                 crossing_index=cross_idx,crossing_ids=event_cross_ids,covered=covered)
        segments.append(seg)
        # Turnover and tracked geometry at a compact 0.1 grid.
        for idx in range(i0,i1+1,100):
            cur=np.argsort(-side*state[0,:,idx])[:N_LEADING]
            retention=len(set(fixed)&set(cur))/N_LEADING
            snap=state[:,fixed,idx]
            _,dd,hh,ll,uu,ee,ii=angular_geometry(snap,origin,p.gap,p.tip_r)
            for rank,fish in enumerate(fixed,1):
                history.append(dict(seed=seed,pair_id=r["pair_id"],condition=r["condition"],
                    origin=origin,time_from_A_minus_2=(idx-i0)*DT,fish=int(fish),leading_rank=rank,
                    leading_retention=retention,x=float(snap[0,rank-1]),y=float(snap[1,rank-1]),
                    heading=float(hh[rank-1]),distance_to_divider=float(dd[rank-1]),
                    signed_angular_error=float(ee[rank-1]),admissible=int(ii[rank-1])))
        st=state[:,fixed,i0:i1+1]
        _,dist,heading,lo,hi,eta,inside=angular_geometry(st,origin,p.gap,p.tip_r)
        for rank,fish in enumerate(fixed,1):
            ent=np.flatnonzero(inside[rank-1])
            entry_rel=int(ent[0]) if len(ent) else None
            crossed_rel=np.flatnonzero(side*state[0,fish,i0:i1+1] >= 0)
            reaches=bool(len(crossed_rel))
            eta0=float(eta[rank-1,0]); sign0=float(np.sign(eta0))
            target_lower = eta0 < 0
            end_rel=entry_rel if entry_rel is not None else i1-i0
            target=np.where(target_lower,lo[rank-1,:end_rel+1],hi[rank-1,:end_rel+1])
            alpha_change=float(np.sum(wrap(np.diff(target)))) if len(target)>1 else 0.0
            fish_rows.append(dict(seed=seed,pair_id=r["pair_id"],condition=r["condition"],origin=origin,
                fish=int(fish),leading_rank=rank,is_eventual_first_crosser=int(fish in event_cross_ids),
                event_first_crosser_covered=int(covered) if r["condition"]=="event" else "",
                initial_distance=float(dist[rank-1,0]),initial_heading=float(heading[rank-1,0]),
                initial_eta=eta0,initial_abs_eta=abs(eta0),initially_admissible=int(inside[rank-1,0]),
                gate_entered=int(entry_rel is not None),entry_time_from_A_minus_2=(entry_rel*DT if entry_rel is not None else ""),
                reached_divider_through_opening=int(reaches),reference_alpha_change_to_entry=alpha_change,
                endpoint_rel=end_rel,initial_eta_sign=sign0,rng_valid="",
                observed_local_heading_change_to_entry="",deterministic_local_impulse_to_entry="",
                stochastic_local_increment_to_entry="",eta_budget_residual="",
                wall_impulse_to_entry=0.0,alignment_impulse_to_entry=0.0,attraction_impulse_to_entry=0.0,
                direct_hydro_impulse_to_entry=0.0,boundary_hydro_residual_impulse_to_entry="",
                wall_impulse_full=0.0,alignment_impulse_full=0.0,attraction_impulse_full=0.0,
                direct_hydro_impulse_full=0.0,boundary_hydro_residual_impulse_full=""))
    return segments,fish_rows,history


def analyze_seed(seed, rows, raw, p, geom, validation_context):
    print(f"seed {seed}: loading",flush=True); t0=time.time()
    state=np.asarray(sio.loadmat(raw,variable_names=["state"],squeeze_me=True)["state"])
    segments,fish_rows,history=build_segments(rows,state,seed,p)
    requested=[]
    for s in segments: requested.extend(range(s["i0"],s["i1"]))
    noise_idx,z,alpha0,rpos0,theta0=rng_noise_for_indices(seed,requested)
    pos={int(k):j for j,k in enumerate(noise_idx)}
    # Exact initialization audit (the two radial arrays are reconstructed below).
    rr=p.R/1.5*rpos0; tt=np.pi*.4*theta0
    initx=-rr*np.cos(tt); inity=rr*np.sin(tt)
    init_err=max(float(np.max(np.abs(wrap(state[2,:,0]-alpha0)))),
                 float(np.max(np.abs(state[0,:,0]-initx))),float(np.max(np.abs(state[1,:,0]-inity))))

    # Exact observed/noise/deterministic budgets for each fixed fish.
    row_lookup={(int(r["pair_id"]),r["condition"],int(r["fish"])):r for r in fish_rows}
    for s in segments:
        inds=np.arange(s["i0"],s["i1"],dtype=int)
        zz=z[[pos[int(k)] for k in inds]]
        for fish in s["fixed"]:
            row=row_lookup[(s["pair_id"],s["condition"],int(fish))]
            rng_valid=not (seed in RNG_DISCONTINUITY and s["i1"]>RNG_DISCONTINUITY[seed])
            row["rng_valid"]=int(rng_valid)
            end=int(row["endpoint_rel"])
            obs=np.diff(state[2,fish,s["i0"]:s["i1"]+1])
            noi=p.In*np.sqrt(p.dt)*zz[:,fish]
            det=obs-noi
            e=max(0,end)
            if not rng_valid:
                row["observed_local_heading_change_to_entry"]=float(s["side"]*np.sum(obs[:e]))
                row["stochastic_local_increment_to_entry"]=np.nan
                row["deterministic_local_impulse_to_entry"]=np.nan
                row["eta_budget_residual"]=np.nan
                row["deterministic_local_impulse_full"]=np.nan
                row["stochastic_local_increment_full"]=np.nan
                row["observed_local_heading_change_full"]=float(s["side"]*np.sum(obs))
                continue
            row["observed_local_heading_change_to_entry"]=float(s["side"]*np.sum(obs[:e]))
            row["stochastic_local_increment_to_entry"]=float(s["side"]*np.sum(noi[:e]))
            row["deterministic_local_impulse_to_entry"]=float(s["side"]*np.sum(det[:e]))
            predicted=(row["deterministic_local_impulse_to_entry"]+
                       row["stochastic_local_increment_to_entry"]-
                       float(row["reference_alpha_change_to_entry"]))
            actual=row["observed_local_heading_change_to_entry"]-float(row["reference_alpha_change_to_entry"])
            row["eta_budget_residual"]=float(actual-predicted)
            row["deterministic_local_impulse_full"]=float(s["side"]*np.sum(det))
            row["stochastic_local_increment_full"]=float(s["side"]*np.sum(noi))
            row["observed_local_heading_change_full"]=float(s["side"]*np.sum(obs))

    # Replay non-BEM source equations on a 0.005 grid. Each frame is computed once.
    active=defaultdict(list)
    for s in segments:
        for idx in range(s["i0"],s["i1"],REPLAY_STRIDE):
            for fish in s["fixed"]:
                active[idx].append((s, int(fish)))
    print(f"seed {seed}: replaying {len(active)} unique frames",flush=True)
    for count,idx in enumerate(sorted(active)):
        refs=active[idx]; selected=sorted({fish for _,fish in refs})
        st=state[:,:,idx]
        al,at=social_selected(st,selected,geom,p); di=direct_selected(st,selected,p); wa=wall_selected(st,selected,geom,p)
        for s,fish in refs:
            row=row_lookup[(s["pair_id"],s["condition"],fish)]
            factor=s["side"]*REPLAY_STRIDE*DT
            rel=idx-s["i0"]
            for key,val in (("wall",wa[fish]),("alignment",al[fish]),("attraction",at[fish]),("direct_hydro",di[fish])):
                row[f"{key}_impulse_full"] += factor*val
                if rel < int(row["endpoint_rel"]): row[f"{key}_impulse_to_entry"] += factor*val
        if count and count%2000==0: print(f"  {count}/{len(active)}",flush=True)

    for row in fish_rows:
        known_entry=sum(float(row[f"{k}_impulse_to_entry"]) for k in ("wall","alignment","attraction","direct_hydro"))
        known_full=sum(float(row[f"{k}_impulse_full"]) for k in ("wall","alignment","attraction","direct_hydro"))
        row["boundary_hydro_residual_impulse_to_entry"]=float(row["deterministic_local_impulse_to_entry"])-known_entry
        row["boundary_hydro_residual_impulse_full"]=float(row["deterministic_local_impulse_full"])-known_full
        sg=float(row["initial_eta_sign"])
        for suffix in ("to_entry","full"):
            row[f"deterministic_closing_{suffix}"]=-sg*float(row[f"deterministic_local_impulse_{suffix}"])
            row[f"stochastic_closing_{suffix}"]=-sg*float(row[f"stochastic_local_increment_{suffix}"])
        row["geometry_closing_to_entry"]=sg*float(row["reference_alpha_change_to_entry"])

    # Two exact full-BEM validation snapshots per seed: compare source replay with RNG-derived residual.
    validation=[]
    for s in (segments[0],segments[len(segments)//2]):
        idx=s["i0"]+1000; rep=replay_snapshot(state[:,:,idx],state[:,:,idx+1],**validation_context,p=p,geom=geom)
        zi=z[pos[idx]]; expected=p.In*np.sqrt(p.dt)*zi
        inferred=(state[2,:,idx+1]-state[2,:,idx]-expected)/p.dt
        known=rep["terms"]["wall"]+rep["terms"]["alignment"]+rep["terms"]["attraction"]+rep["terms"]["direct_hydro_effective"]
        bem_inferred=inferred-known
        raw=(rep["terms"]["wall"]+rep["terms"]["alignment"]+rep["terms"]["attraction"]+
             rep["terms"]["direct_hydro"]+rep["terms"]["boundary_hydro"])
        clipped=np.abs(raw)>np.pi/p.dt
        valid=~clipped & ~rep["overridden"]
        chosen=np.argsort(-s["side"]*state[0,:,idx])[:N_LEADING]
        al,at=social_selected(state[:,:,idx],chosen,geom,p); di=direct_selected(state[:,:,idx],chosen,p); wa=wall_selected(state[:,:,idx],chosen,geom,p)
        rng_valid=not (seed in RNG_DISCONTINUITY and idx>=RNG_DISCONTINUITY[seed])
        validation.append(dict(seed=seed,index=idx,time=idx*DT,rng_valid=int(rng_valid),initialization_max_abs_error=init_err,
            rng_noise_max_abs_error=float(np.max(np.abs(rep["noise_increment"]-expected))),
            deterministic_total_max_abs_error=float(np.max(np.abs(inferred-rep["deterministic_total"]))),
            boundary_torque_rmse_unclipped=float(np.sqrt(np.mean((bem_inferred[valid]-rep["terms"]["boundary_hydro_effective"][valid])**2))),
            clipped_count=int(np.sum(clipped)),override_count=int(np.sum(rep["overridden"])),
            selected_wall_max_abs_error=float(max(abs(wa[i]-rep["terms"]["wall"][i]) for i in chosen)),
            selected_alignment_max_abs_error=float(max(abs(al[i]-rep["terms"]["alignment"][i]) for i in chosen)),
            selected_attraction_max_abs_error=float(max(abs(at[i]-rep["terms"]["attraction"][i]) for i in chosen)),
            selected_direct_max_abs_error=float(max(abs(di[i]-rep["terms"]["direct_hydro"][i]) for i in chosen))))
    print(f"seed {seed}: done in {time.time()-t0:.1f}s",flush=True)
    return fish_rows,history,validation


def seed_statistics(fish_rows, history, segments):
    seed_rows=[]; pair_rows=[]
    metrics=("initial_abs_eta","gate_entered","reached_divider_through_opening",
             "deterministic_closing_full","stochastic_closing_full",
             "wall_impulse_full","alignment_impulse_full","attraction_impulse_full",
             "direct_hydro_impulse_full","boundary_hydro_residual_impulse_full")
    for seed in SEEDS:
        for cond in ("event","control"):
            q=[r for r in fish_rows if int(r["seed"])==seed and r["condition"]==cond]
            row=dict(seed=seed,condition=cond,n_pairs=len({int(r["pair_id"]) for r in q}),n_fixed_fish=len(q))
            for m in metrics: row[m]=float(np.nanmean([float(r[m]) for r in q]))
            seed_rows.append(row)
    for key in sorted({(int(r["seed"]),int(r["pair_id"])) for r in fish_rows}):
        for cond in ("event","control"):
            q=[r for r in fish_rows if (int(r["seed"]),int(r["pair_id"]))==key and r["condition"]==cond]
            row=dict(seed=key[0],pair_id=key[1],condition=cond)
            for m in metrics: row[m]=float(np.mean([float(r[m]) for r in q]))
            pair_rows.append(row)
    contrasts=[]; rng=np.random.default_rng(20261022)
    for m in metrics:
        vals=[]
        for seed in SEEDS:
            e=next(r for r in seed_rows if r["seed"]==seed and r["condition"]=="event")
            c=next(r for r in seed_rows if r["seed"]==seed and r["condition"]=="control")
            vals.append(e[m]-c[m])
        vals=np.asarray(vals); boot=np.mean(vals[rng.integers(0,len(vals),(30000,len(vals)))],axis=1)
        contrasts.append(dict(metric=m,event_minus_control_mean=float(vals.mean()),
            seed_bootstrap_ci_low=float(np.percentile(boot,2.5)),seed_bootstrap_ci_high=float(np.percentile(boot,97.5)),
            independent_seeds=len(vals),seed_differences=";".join(f"{x:.8g}" for x in vals)))
    # Coverage and turnover are event/pair-level but uncertainty remains seed-level.
    coverage=[]
    for seed in SEEDS:
        ev=[s for s in segments if s["seed"]==seed and s["condition"]=="event"]
        coverage.append(np.mean([s["covered"] for s in ev]))
    boot=np.mean(np.asarray(coverage)[rng.integers(0,20,(30000,20))],axis=1)
    coverage_summary=dict(metric="eventual_first_crosser_in_fixed_ten",mean=float(np.mean(coverage)),
        seed_bootstrap_ci_low=float(np.percentile(boot,2.5)),seed_bootstrap_ci_high=float(np.percentile(boot,97.5)),independent_seeds=20)
    return pair_rows,seed_rows,contrasts,coverage_summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True); ap.add_argument("--raw-dir",type=Path,required=True)
    ap.add_argument("--risk-set",type=Path,required=True); ap.add_argument("--geometry",type=Path,required=True)
    ap.add_argument("--lu",type=Path,required=True); ap.add_argument("--piv",type=Path,required=True)
    ap.add_argument("--seeds",type=int,nargs="*",default=list(SEEDS))
    args=ap.parse_args(); out=args.out.resolve(); (out/"tables").mkdir(parents=True,exist_ok=True)
    p=Parameters(); geom=analytic_geometry(p)
    risk=read_csv(args.risk_set)
    g=sio.loadmat(args.geometry); normal=np.asarray(g["normal_vec"])
    boundary=dict(x=np.asarray(g["Sx"]).reshape(-1),y=np.asarray(g["Sy"]).reshape(-1),
        weights=np.asarray(g["weights"]).reshape(-1),nx=normal[:,0],ny=normal[:,1])
    context=dict(boundary=boundary,lu=np.load(args.lu,mmap_mode="r"),piv=np.load(args.piv))
    all_fish=[]; all_history=[]; validations=[]; all_segments=[]
    chosen=tuple(args.seeds)
    for seed in chosen:
        rows=[r for r in risk if int(r["seed"])==seed]
        raw=args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        fish,hist,val=analyze_seed(seed,rows,raw,p,geom,context)
        all_fish.extend(fish); all_history.extend(hist); validations.extend(val)
        # Compact segment reconstruction from fish table for coverage accounting.
        for cond in ("event","control"):
            for pid in sorted({int(r["pair_id"]) for r in fish if r["condition"]==cond}):
                q=[r for r in fish if int(r["pair_id"])==pid and r["condition"]==cond]
                all_segments.append(dict(seed=seed,pair_id=pid,condition=cond,
                    covered=bool(int(q[0]["event_first_crosser_covered"])) if cond=="event" else False))
    if chosen != SEEDS:
        write_csv(out/"tables"/"fixed_leading_fish.csv",all_fish)
        write_csv(out/"tables"/"tracked_fish_history.csv.gz",all_history)
        write_csv(out/"tables"/"replay_validation.csv",validations)
        return
    pair_rows,seed_rows,contrasts,coverage=seed_statistics(all_fish,all_history,all_segments)
    write_csv(out/"tables"/"fixed_leading_fish.csv",all_fish)
    write_csv(out/"tables"/"tracked_fish_history.csv.gz",all_history)
    write_csv(out/"tables"/"pair_summary.csv",pair_rows); write_csv(out/"tables"/"seed_summary.csv",seed_rows)
    write_csv(out/"tables"/"seed_paired_contrasts.csv",contrasts); write_csv(out/"tables"/"replay_validation.csv",validations)
    write_csv(out/"summary.csv",[coverage]+contrasts)
    config=dict(source_commit="4259ba7",n_seeds=20,n_pairs=456,n_leading=10,dt=DT,replay_stride=REPLAY_STRIDE,
        replay_dt=DT*REPLAY_STRIDE,noise_reconstruction="RandomState(seed), 300 initialization uniforms, then 100 normals/step",
        rng_discontinuity_steps=RNG_DISCONTINUITY,
        risk_set=str(args.risk_set),raw_dir=str(args.raw_dir),potential_fitted=False,new_simulation=False)
    (out/"configuration.json").write_text(json.dumps(config,indent=2)+"\n")

if __name__=="__main__": main()
