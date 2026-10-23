#!/usr/bin/env python3
"""Prospective rounded-tip one-sided-drift test on archived BEM trajectories."""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np
import scipy.io as sio

from exact_tip_dynamics import (
    Parameters, analytic_geometry, chamber_side, classify, exact_wall,
    FocalBEM, focal_replay, global_to_local, local_to_global, ray_wall_branch,
    tangent_geometry, visible_neighbors, wrap,
)

DT = 1e-3
OFFSETS = (0.0, 1.5, 1.8)
EPSILONS = (1e-3, 5e-3)


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def corridor_series(state, fish, indices, origin, p):
    side = chamber_side(origin)
    x = state[0, fish, indices]
    y = state[1, fish, indices]
    heading = global_to_local(state[2, fish, indices], side)
    d = -side*x
    vlo = -p.gap/2-p.tip_r-y
    vhi = p.gap/2+p.tip_r-y
    rlo = np.hypot(d, vlo)
    rhi = np.hypot(d, vhi)
    valid = (d > 0) & (rlo > p.tip_r) & (rhi > p.tip_r)
    lo = np.full_like(d, np.nan)
    hi = np.full_like(d, np.nan)
    lo[valid] = np.arctan2(vlo[valid], d[valid])+np.arcsin(p.tip_r/rlo[valid])
    hi[valid] = np.arctan2(vhi[valid], d[valid])-np.arcsin(p.tip_r/rhi[valid])
    valid &= hi > lo
    inside = valid & (heading >= lo) & (heading <= hi)
    return valid, inside, lo, hi


def outcomes(state, idx, fish, origin, p, horizon):
    steps = int(round(horizon/DT))
    inds = np.arange(idx+1, min(idx+steps+1, state.shape[2]))
    _, inside, _, _ = corridor_series(state, fish, inds, origin, p)
    side = chamber_side(origin)
    focal_cross = bool(np.any(side*state[0, fish, inds] >= 0))
    any_cross = bool(np.any(side*state[0, :, inds] >= 0))
    return int(np.any(inside)), int(focal_cross), int(any_cross)


def evaluate_snapshot(state, idx, focal, meta, boundary, lu, piv, p, geom):
    origin = meta["origin"]
    side = chamber_side(origin)
    velocity = (state[:2, :, idx+1]-state[:2, :, idx])/DT
    focal_velocity = velocity[:, focal]
    gate = tangent_geometry(state[:2, focal, idx], focal_velocity, origin, p)
    theta_local = float(global_to_local(state[2, focal, idx], side))
    current_fish = state[:, focal, idx]
    current_hit = ray_wall_branch(current_fish, geom, p)
    phase_direction = "up" if np.mean(velocity[1]) >= 0 else "down"
    base = dict(seed=int(meta["seed"]), pair_id=int(meta["pair_id"]),
                condition=meta["condition"], origin=origin,
                risk_time=float(meta["risk_time"]), observation_offset=float(meta["offset"]),
                observation_time=idx*DT, focal_fish=int(focal), x=float(state[0, focal, idx]),
                y=float(state[1, focal, idx]), theta_global=float(state[2, focal, idx]),
                theta_local=theta_local, divider_distance=float(gate["d"]),
                collective_vertical_phase=phase_direction,
                collective_vertical_velocity=float(np.mean(velocity[1])),
                alpha_minus=gate["alpha_minus"], alpha_plus=gate["alpha_plus"],
                alpha_minus_dot=gate["alpha_minus_dot"], alpha_plus_dot=gate["alpha_plus_dot"],
                corridor_width=gate["width"], corridor_valid=int(gate["valid"]),
                current_ray_branch=current_hit["branch"],
                current_ray_distance=current_hit["distance"], Iw=p.Iw)
    for horizon in (0.2, 0.5):
        ent, focal_cross, any_cross = outcomes(state, idx, focal, origin, p, horizon)
        tag = str(horizon).replace(".", "p")
        base[f"angular_entry_{tag}"] = ent
        base[f"focal_crossing_{tag}"] = focal_cross
        base[f"any_crossing_{tag}"] = any_cross
    if not gate["valid"]:
        rows = []
        for tip in ("upper", "lower"):
            for eps in EPSILONS:
                rows.append(dict(base, tip=tip, epsilon=eps, initial_q=np.nan,
                                 initially_outside="", line_of_sight=0,
                                 mu_out=np.nan, mu_in=np.nan, classification="invalid"))
        return rows

    # Mid-corridor must pass between the rounded tips rather than select a near solid wall.
    midpoint = .5*(gate["alpha_minus"]+gate["alpha_plus"])
    mid_global = local_to_global(midpoint, side)
    mid_hit = ray_wall_branch(np.r_[state[:2, focal, idx], mid_global], geom, p)
    los = mid_hit["valid"] and mid_hit["branch"] not in {
        "upper_tip", "lower_tip", "left_vertical", "right_vertical"}
    base["midcorridor_ray_branch"] = mid_hit["branch"]
    base["line_of_sight"] = int(los)
    visible = visible_neighbors(state[:2, :, idx], geom, p)
    bem = FocalBEM.build(state[:, :, idx], focal, boundary, lu, piv, p)
    rows = []
    for tip in ("upper", "lower"):
        if tip == "upper":
            alpha, alpha_dot = gate["alpha_plus"], gate["alpha_plus_dot"]
            initial_q = theta_local-alpha
            direction = ((+1, "out"), (-1, "in"))
        else:
            alpha, alpha_dot = gate["alpha_minus"], gate["alpha_minus_dot"]
            initial_q = alpha-theta_local
            direction = ((-1, "out"), (+1, "in"))
        for eps in EPSILONS:
            evaluated = {}
            for sign, label in direction:
                trial_local = alpha+sign*eps
                trial_global = local_to_global(trial_local, side)
                rep = focal_replay(state[:, :, idx], focal, trial_global,
                                   visible, bem, geom, p)
                local_rate = side*rep["total"]
                mu = local_rate-alpha_dot if tip == "upper" else alpha_dot-local_rate
                evaluated[label] = (mu, rep, trial_local, trial_global)
            mu_out, out_rep = evaluated["out"][0], evaluated["out"][1]
            mu_in, in_rep = evaluated["in"][0], evaluated["in"][1]
            row = dict(base, tip=tip, epsilon=eps, initial_q=float(initial_q),
                       initial_abs_angular_error=float(max(initial_q, 0.0)),
                       initially_outside=int(initial_q > 0), line_of_sight=int(los),
                       mu_out=float(mu_out), mu_in=float(mu_in),
                       classification=(classify(mu_out, mu_in) if los else "invalid"),
                       out_ray_branch=out_rep["hit"]["branch"],
                       in_ray_branch=in_rep["hit"]["branch"],
                       out_ray_distance=out_rep["hit"]["distance"],
                       in_ray_distance=in_rep["hit"]["distance"],
                       out_override=int(out_rep["override"]), in_override=int(in_rep["override"]),
                       out_total_raw=out_rep["total_raw"], in_total_raw=in_rep["total_raw"])
            for side_name, rep in (("out", out_rep), ("in", in_rep)):
                for term, value in rep["terms"].items():
                    row[f"{side_name}_{term}"] = value
                row[f"{side_name}_total"] = rep["total"]
            rows.append(row)
    return rows


def analyze_seed(seed, risk_rows, raw, boundary, lu, piv, p, geom, offsets):
    started = time.time()
    print(f"seed {seed}: loading {raw.name}", flush=True)
    state = np.asarray(sio.loadmat(raw, variable_names=["state"], squeeze_me=True)["state"])
    rows = []
    geometry_candidates = []
    for count, rr in enumerate(risk_rows):
        side = chamber_side(rr["origin"])
        for offset in offsets:
            idx = int(round((float(rr["risk_time"])+offset)/DT))
            # Prospective leading fish: closest fish in the origin chamber at this snapshot.
            d = -side*state[0, :, idx]
            eligible = np.flatnonzero(d > 0)
            focal = int(eligible[np.argmin(d[eligible])])
            meta = dict(rr, seed=seed, offset=offset)
            snaprows = evaluate_snapshot(state, idx, focal, meta, boundary, lu, piv, p, geom)
            rows.extend(snaprows)
            velocity = (state[:2, focal, idx+1]-state[:2, focal, idx])/DT
            gate = tangent_geometry(state[:2, focal, idx], velocity, rr["origin"], p)
            geometry_candidates.append(dict(seed=seed, pair_id=rr["pair_id"],
                condition=rr["condition"], origin=rr["origin"], offset=offset, index=idx,
                fish=focal, x=state[0, focal, idx], y=state[1, focal, idx],
                width=gate["width"], valid=int(gate["valid"])))
        if (count+1) % 10 == 0:
            print(f"  {count+1}/{len(risk_rows)} risk rows", flush=True)
    print(f"seed {seed}: {len(rows)} tangent rows in {time.time()-started:.1f}s", flush=True)
    return rows, geometry_candidates


def planar_identity(p, geom):
    # Left chamber, upper flat divider wall x=-r; wall-only Euler integration.
    dt = 1e-5
    fish = np.array([-2.5, 5.2, 0.18], float)
    d0 = -p.tip_r-fish[0]
    theta0 = fish[2]
    rows = []
    for step in range(500000):
        omega, hit = exact_wall(fish, geom, p)
        d = -p.tip_r-fish[0]
        predicted = theta0+p.Iw*np.log(d0/d)
        if step % 100 == 0:
            rows.append(dict(step=step, time=step*dt, x=fish[0], y=fish[1],
                             theta=fish[2], perpendicular_distance=d,
                             ray_distance=hit["distance"], branch=hit["branch"],
                             amplified=int(hit.get("amplified", False)),
                             clipped=int(hit.get("wall_clipped", False)),
                             identity_prediction=predicted,
                             identity_error=wrap(fish[2]-predicted)))
        if hit["branch"] != "left_vertical" or hit["distance"] <= .11 or d <= .15:
            break
        fish[:2] += dt*np.array([np.cos(fish[2]), np.sin(fish[2])])
        fish[2] += dt*omega
    return rows


def tangent_sweeps(candidates, geom, p):
    valid = [r for r in candidates if int(r["valid"])]
    picks = []
    upper = [r for r in valid if float(r["y"]) > 0]
    lower = [r for r in valid if float(r["y"]) < 0]
    if upper: picks.append(("upper_archived", max(upper, key=lambda r: float(r["y"]))))
    if lower: picks.append(("lower_archived", min(lower, key=lambda r: float(r["y"]))))
    closed = [r for r in candidates if not int(r["valid"])]
    if closed: picks.append(("closed_archived", min(closed, key=lambda r: float(r["width"]) if np.isfinite(float(r["width"])) else -999)))
    rows = []
    for name, r in picks:
        origin = r["origin"]
        side = chamber_side(origin)
        pos = np.array([float(r["x"]), float(r["y"])])
        gate = tangent_geometry(pos, np.zeros(2), origin, p)
        boundaries = (("upper", gate["alpha_plus"]), ("lower", gate["alpha_minus"]))
        for tip, alpha in boundaries:
            if not np.isfinite(alpha):
                continue
            for delta in np.linspace(-.03, .03, 121):
                theta_local = alpha+delta
                fish = np.r_[pos, local_to_global(theta_local, side)]
                wall, hit = exact_wall(fish, geom, p)
                rows.append(dict(case=name, source_seed=r["seed"], source_pair=r["pair_id"],
                    source_condition=r["condition"], origin=origin, x=pos[0], y=pos[1],
                    corridor_valid=int(gate["valid"]), corridor_width=gate["width"],
                    tip=tip, tangent_angle=alpha, delta=delta, theta_local=theta_local,
                    ray_branch=hit["branch"], ray_distance=hit["distance"], wall_torque=wall,
                    amplified=int(hit.get("amplified", False)), clipped=int(hit.get("wall_clipped", False))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--risk-set", type=Path, required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--lu", type=Path, required=True)
    ap.add_argument("--piv", type=Path, required=True)
    ap.add_argument("--seeds", type=int, nargs="*", default=list(range(20)))
    ap.add_argument("--offsets", type=float, nargs="*", default=list(OFFSETS))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    out = args.out.resolve()
    (out/"tables"/"by_seed").mkdir(parents=True, exist_ok=True)
    p = Parameters()
    geom = analytic_geometry(p)
    g = sio.loadmat(args.geometry)
    normal = np.asarray(g["normal_vec"])
    boundary = dict(x=np.asarray(g["Sx"]).reshape(-1), y=np.asarray(g["Sy"]).reshape(-1),
                    weights=np.asarray(g["weights"]).reshape(-1),
                    nx=normal[:, 0], ny=normal[:, 1])
    lu = np.load(args.lu, mmap_mode="r")
    piv = np.load(args.piv)
    risk = read_csv(args.risk_set)
    all_rows, all_candidates = [], []
    for seed in args.seeds:
        target = out/"tables"/"by_seed"/f"seed_{seed:02d}.csv"
        ctarget = out/"tables"/"by_seed"/f"seed_{seed:02d}_geometry.csv"
        if target.exists() and ctarget.exists() and not args.force:
            rows = read_csv(target)
            candidates = read_csv(ctarget)
            print(f"seed {seed}: using checkpoint", flush=True)
        else:
            seedrisk = [r for r in risk if int(r["seed"]) == seed]
            raw = args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
            rows, candidates = analyze_seed(seed, seedrisk, raw, boundary, lu, piv, p, geom,
                                            tuple(args.offsets))
            write_csv(target, rows)
            write_csv(ctarget, candidates)
        all_rows.extend(rows)
        all_candidates.extend(candidates)
    write_csv(out/"tables"/"tangent_snapshot_replay.csv", all_rows)
    write_csv(out/"tables"/"geometry_candidates.csv", all_candidates)
    ident = planar_identity(p, geom)
    write_csv(out/"tables"/"planar_wall_identity.csv", ident)
    sweeps = tangent_sweeps(all_candidates, geom, p)
    write_csv(out/"tables"/"ray_branch_sweeps.csv", sweeps)
    config = dict(source_branch="analysis/social-alignment-angular-landscape",
                  source_commit="e8fef2091dc949e72dcb7bd5bf453f925b32d8ff",
                  seeds=args.seeds, offsets=args.offsets, epsilons=list(EPSILONS),
                  focal_selection="closest fish still in origin chamber at observation time",
                  no_fitted_coefficients=True, new_simulations=False,
                  parameters=p.__dict__, raw_dir=str(args.raw_dir), risk_set=str(args.risk_set),
                  production_source="../BEM_bichamber_parallel_instance.py",
                  exact_replay_source="pre_crossing_turning/scripts/replay_dynamics.py")
    (out/"configuration.json").write_text(json.dumps(config, indent=2)+"\n")


if __name__ == "__main__":
    main()
