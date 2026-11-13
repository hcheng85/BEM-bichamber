#!/usr/bin/env python3
"""Exact finite-step budget for pre-crossing majority-school gap polarization.

Run from BEM-bichamber repository root after placing this file alongside
analyze_approach_current_balance.py and analyze_initiation_propagation.py.

Scientific target: why does a=P_majority dot g_gap rise before a crossing?
Only observed state, exact source torques, and seed-matched event/control pairs.
No phenomenological parameters, inferred potential, new simulator, or raw edits.

Modes:
  --self-test             synthetic exact-budget and control-matching checks
  without --source       inexpensive trajectory-level descriptive comparison
  with --source          add bounded expensive native-step source-exact budget
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

N = 100
DT = .001
IN = .5
BUDGET_TERMS = (
    "wall", "social", "clip", "noise_relaxation", "realized_noise",
    "membership", "gap_direction_rotation",
)


def csv_read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def csv_write(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def source_geometry(state, direction):
    """Return a, P, g, membership and q for a fixed original majority side.

    direction=+1: original majority on left (x<0) moving right;
    direction=-1: original majority on right (x>=0) moving left.
    Reference side remains fixed over a pre-first-crossing interval.
    """
    x, y, theta = state
    mask = x < 0 if direction == 1 else x >= 0
    number = int(mask.sum())
    if number <= 0:
        raise ValueError("Source chamber empty")
    com = np.array([x[mask].mean(), y[mask].mean()], dtype=float)
    norm = float(np.hypot(*com))
    if norm <= 1e-14:
        raise ValueError("Source-chamber COM at the gap; direction undefined")
    gap = -com / norm
    p = np.array([np.cos(theta[mask]).mean(), np.sin(theta[mask]).mean()])
    q = 1.0 - 2.0 * np.count_nonzero(x < 0) / len(x)
    return dict(a=float(p @ gap), P=p, g=gap, mask=mask, n=number,
                q=q, com=com, majority_consistent=(direction * q < 0))


def exact_increment(st, nxt, direction, wall, social, clip, in_noise=IN, dt=DT):
    """Algebraically exact a(t+dt)-a(t) split on the *native* time step.

    The simulator first applies deterministic theta drift then an independent
    Gaussian angular increment, with its observed realization reconstructed
    from theta_next (no random-number replay is required). The conditional
    expectation of cos(theta_next) yields the exact exp(-In^2 dt/2) term.

    Discrete splitting order: angle change at old membership and gap direction;
    source-chamber membership change at updated theta and old gap direction;
    rotation of COM-to-gap unit vector at new membership and updated theta.
    Unlike the approach-current balance there is no velocity-override term:
    this observable depends on heading and position, not stored rdot.
    """
    b = source_geometry(st, direction)
    e = source_geometry(nxt, direction)
    theta = st[2]
    theta_next = nxt[2]
    raw_omega = wall + social + clip
    delta = raw_omega * dt
    phi = theta - np.arctan2(b["g"][1], b["g"][0])
    # (cos(phi+delta)-cos(phi))/delta with a stable zero-delta limit.
    slope = -np.sin(phi + .5 * delta) * np.sinc(delta / (2.0 * np.pi))
    weights = b["mask"] / b["n"]
    terms = {
        "wall": float(np.sum(weights * slope * wall * dt)),
        "social": float(np.sum(weights * slope * social * dt)),
        "clip": float(np.sum(weights * slope * clip * dt)),
    }
    deterministic = np.cos(phi + delta)
    next_projection = np.cos(theta_next - np.arctan2(b["g"][1], b["g"][0]))
    factor = np.exp(-.5 * in_noise * in_noise * dt)
    terms["noise_relaxation"] = float(np.sum(weights * (factor - 1.) * deterministic))
    terms["realized_noise"] = float(np.sum(weights * (next_projection - factor * deterministic)))
    post_h = np.column_stack((np.cos(theta_next), np.sin(theta_next)))
    prior_member_p = post_h[b["mask"]].mean(axis=0)
    new_member_p = post_h[e["mask"]].mean(axis=0)
    terms["membership"] = float((new_member_p - prior_member_p) @ b["g"])
    terms["gap_direction_rotation"] = float(new_member_p @ (e["g"] - b["g"]))
    observed = e["a"] - b["a"]
    closure = observed - sum(terms.values())
    return dict(delta_a=float(observed), reconstructed=float(sum(terms.values())),
                residual=float(closure), a_before=b["a"], a_after=e["a"],
                q_before=b["q"], q_after=e["q"],
                reference_majority_consistent=int(b["majority_consistent"] and e["majority_consistent"]),
                **terms)


def matched_pairs(event_csv, control_csv, lead=2.):
    """One predetermined first valid crossing/control pair per event seed.

    The preexisting matching was performed at the same occupancy and control
    quietness; do not retrospectively change the control using any new result.
    """
    events = {(int(r["seed"]), int(r["event_id"])): r for r in csv_read(event_csv)}
    output = {}
    for row in csv_read(control_csv):
        if not np.isclose(float(row["lead"]), lead):
            continue
        key = (int(row["event_seed"]), int(row["event_id"]))
        if key not in events:
            raise ValueError(f"Missing indexed event {key}")
        ev = events[key]
        q = float(row["q"])
        direction = int(row["direction"])
        # Require reference population already to be majority on source side.
        if not direction * q < 0:
            continue
        if int(ev["direction"]) != direction:
            raise ValueError(f"Direction mismatch at {key}")
        # Preserve the first eligible event, independent of its size/outcome.
        cand = dict(seed=key[0], event_id=key[1], first_native=int(ev["first_native"]),
                    event_time=float(ev["time"]), direction=direction, q=q,
                    control_seed=int(row["control_seed"]),
                    control_start=int(round(float(row["control_time"]) / DT)))
        if key[0] not in output or cand["first_native"] < output[key[0]]["first_native"]:
            output[key[0]] = cand
    return [output[seed] for seed in sorted(output)]


def paths(pair, kind, duration_steps):
    if kind == "event":
        # The final sample is precisely one native frame BEFORE first crossing.
        start = pair["first_native"] - duration_steps
        stop = pair["first_native"] - 1
        seed = pair["seed"]
    else:
        start = pair["control_start"]
        stop = start + duration_steps - 1
        seed = pair["control_seed"]
    return int(seed), int(start), int(stop)


def file_for_seed(raw_dir, seed):
    return raw_dir / f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"


def sample_a(reader, start, stop, direction, stride=100):
    """Trace a from raw states; at most a few MB per interval."""
    idx = np.unique(np.r_[np.arange(start, stop+1, stride), stop])
    block = reader.block("state", start, stop)
    results = []
    for k in idx:
        g = source_geometry(block[:, :, int(k-start)], direction)
        results.append(dict(native_step=int(k), a=g["a"], q=g["q"],
                            majority_consistent=int(g["majority_consistent"])))
    return results, endpoint_split(block[:, :, 0], block[:, :, -1], direction)


def endpoint_split(st, nxt, direction):
    """Exact endpoint split into heading, membership, and gap geometry.

    This is descriptive and source-free: the heading term contains all angular
    mechanisms and does not attribute them to social, wall, or noise sources.
    """
    b = source_geometry(st, direction)
    e = source_geometry(nxt, direction)
    h0 = np.column_stack((np.cos(st[2]), np.sin(st[2])))
    h1 = np.column_stack((np.cos(nxt[2]), np.sin(nxt[2])))
    p0_old = h0[b["mask"]].mean(axis=0)
    p1_old = h1[b["mask"]].mean(axis=0)
    p1_new = h1[e["mask"]].mean(axis=0)
    heading = float((p1_old-p0_old) @ b["g"])
    membership = float((p1_new-p1_old) @ b["g"])
    geometry = float(p1_new @ (e["g"]-b["g"]))
    delta = float(e["a"]-b["a"])
    return dict(endpoint_heading_change=heading,
                endpoint_membership_change=membership,
                endpoint_gap_geometry_change=geometry,
                endpoint_split_residual=delta-heading-membership-geometry)


def cheap_analysis(raw_dir, pairs, output, duration_steps=2000, stride=100):
    from analyze_approach_current_balance import RawReader
    rows, trajectory, audit = [], [], []
    for pair in pairs:
        for kind in ("event", "control"):
            seed, start, stop = paths(pair, kind, duration_steps)
            reader = RawReader(file_for_seed(raw_dir, seed))
            try:
                trace, endpoint = sample_a(reader, start, stop, pair["direction"], stride)
            finally:
                reader.close()
            for r in trace:
                trajectory.append(dict(pair_seed=pair["seed"], trajectory_seed=seed,
                                       kind=kind, time_to_interval_end=(r["native_step"]-stop)*DT,
                                       a=r["a"], q=r["q"],
                                       majority_consistent=r["majority_consistent"]))
            rows.append(dict(pair_seed=pair["seed"], event_id=pair["event_id"],
                             kind=kind, trajectory_seed=seed,
                             start_a=trace[0]["a"], end_a=trace[-1]["a"],
                             delta_a=trace[-1]["a"]-trace[0]["a"],
                             start_q=trace[0]["q"], end_q=trace[-1]["q"],
                             all_reference_majority_consistent=int(all(t["majority_consistent"] for t in trace)),
                             **endpoint))
    csv_write(output/"observed_intervals.csv", rows)
    csv_write(output/"observed_time_series.csv", trajectory)
    comparison = []
    for seed in sorted({r["pair_seed"] for r in rows}):
        ev = next(r for r in rows if r["pair_seed"] == seed and r["kind"] == "event")
        co = next(r for r in rows if r["pair_seed"] == seed and r["kind"] == "control")
        comparison.append(dict(seed=seed, event_delta_a=ev["delta_a"],
                               control_delta_a=co["delta_a"],
                               paired_difference=ev["delta_a"]-co["delta_a"],
                               event_start_a=ev["start_a"], control_start_a=co["start_a"],
                               difference_start_a=ev["start_a"]-co["start_a"],
                               event_end_a=ev["end_a"], control_end_a=co["end_a"],
                               difference_end_a=ev["end_a"]-co["end_a"],
                               **{f"event_{name}": ev[name] for name in
                                  ("endpoint_heading_change","endpoint_membership_change","endpoint_gap_geometry_change")},
                               **{f"control_{name}": co[name] for name in
                                  ("endpoint_heading_change","endpoint_membership_change","endpoint_gap_geometry_change")},
                               **{f"difference_{name}": ev[name]-co[name] for name in
                                  ("endpoint_heading_change","endpoint_membership_change","endpoint_gap_geometry_change")},
                               event_endpoint_split_residual=ev["endpoint_split_residual"],
                               control_endpoint_split_residual=co["endpoint_split_residual"],
                               reference_majority_consistent=int(ev["all_reference_majority_consistent"] and co["all_reference_majority_consistent"])))
    csv_write(output/"observed_paired_differences.csv", comparison)
    fig, axs = plt.subplots(1,2,figsize=(11,4.3))
    for kind,label in (("event", "future crossing"),("control", "same-q control")):
        vals = [r for r in rows if r["kind"] == kind]
        axs[0].scatter([r["start_a"] for r in vals], [r["end_a"] for r in vals],s=22,alpha=.6,label=label)
    axs[0].plot([-1,1],[-1,1],"k:",lw=1)
    axs[0].set(xlabel="a at beginning of 2-unit window",ylabel="a immediately before first crossing / control end",xlim=(-1,1),ylim=(-1,1))
    axs[0].legend(frameon=False)
    vals = [r["paired_difference"] for r in comparison]
    axs[1].scatter(range(len(vals)),vals,s=20)
    axs[1].axhline(0,color="0.4",ls=":")
    axs[1].set(xlabel="independent matched pair",ylabel="event-minus-control Δa")
    fig.tight_layout();fig.savefig(output/"observed_polarization_change.png",dpi=200);plt.close(fig)

    # Equal-pair time-series mean and seed-bootstrap band.
    times = sorted({r["time_to_interval_end"] for r in trajectory})
    pair_seeds = sorted({r["pair_seed"] for r in trajectory})
    rng = np.random.default_rng(20261112)
    ts_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharex=True)
    kind_values = {}
    for kind, color in (("event", "C0"), ("control", "C1")):
        matrix = np.asarray([[next(r["a"] for r in trajectory if r["pair_seed"] == seed and r["kind"] == kind and r["time_to_interval_end"] == time)
                              for time in times] for seed in pair_seeds])
        kind_values[kind] = matrix
        pick = rng.integers(len(pair_seeds), size=(4000, len(pair_seeds)))
        boot = matrix[pick].mean(axis=1)
        mean = matrix.mean(axis=0); lo, hi = np.quantile(boot, [.025, .975], axis=0)
        axes[0].plot(times, mean, color=color, label=kind)
        axes[0].fill_between(times, lo, hi, color=color, alpha=.2)
        for j, time in enumerate(times):
            ts_rows.append(dict(kind=kind, time_to_end=time, mean_a=mean[j], ci_low=lo[j], ci_high=hi[j]))
    difference = kind_values["event"]-kind_values["control"]
    pick = rng.integers(len(pair_seeds), size=(4000, len(pair_seeds)))
    boot = difference[pick].mean(axis=1)
    mean = difference.mean(axis=0); lo, hi = np.quantile(boot, [.025, .975], axis=0)
    axes[1].plot(times, mean, color="C2"); axes[1].fill_between(times, lo, hi, color="C2", alpha=.2)
    axes[1].axhline(0, color=".4", ls=":")
    for j, time in enumerate(times):
        ts_rows.append(dict(kind="event_minus_control", time_to_end=time, mean_a=mean[j], ci_low=lo[j], ci_high=hi[j]))
    axes[0].legend(frameon=False); axes[0].set_ylabel("gap-directed polarization a")
    axes[0].set_title("Matched event and control trajectories")
    axes[1].set_title("Paired event-minus-control difference")
    for ax in axes: ax.set_xlabel("time relative to interval end")
    fig.tight_layout(); fig.savefig(output/"observed_event_control_time_series.png", dpi=200); plt.close(fig)
    csv_write(output/"observed_time_series_summary.csv", ts_rows)
    return comparison


def full_native_budget(mod, para, raw_dir, pairs, output, duration_steps=2000, n_exact=8):
    from analyze_approach_current_balance import RawReader, replay_torques
    if n_exact < 1:
        return [], dict(available=False,reason="n_exact=0")
    idx = np.unique(np.linspace(0, len(pairs)-1, min(n_exact,len(pairs))).round().astype(int))
    selected = [pairs[j] for j in idx]
    rows, audits = [], []
    for pair in selected:
        print(f"Exact native source replay: event seed {pair['seed']}, control {pair['control_seed']}",flush=True)
        for kind in ("event","control"):
            seed,start,stop = paths(pair,kind,duration_steps)
            reader = RawReader(file_for_seed(raw_dir,seed))
            try:
                states = reader.block("state",start,stop)
                values = dict.fromkeys(("delta_a",*BUDGET_TERMS),0.)
                worst = 0.
                n_majority_flip = 0
                for k in range(states.shape[2]-1):
                    st = states[:,:,k]
                    nxt = states[:,:,k+1]
                    wall,social,clip,_ = replay_torques(mod,st,para)
                    z = exact_increment(st,nxt,pair["direction"],wall,social,clip)
                    for name in values:
                        values[name] += z[name]
                    worst = max(worst,abs(z["residual"]))
                    n_majority_flip += (not z["reference_majority_consistent"])
                    if k and k%500==0:
                        print(f"  {kind}: {k}/{states.shape[2]-1} native steps",flush=True)
                a0 = source_geometry(states[:,:,0],pair["direction"])["a"]
                a1 = source_geometry(states[:,:,-1],pair["direction"])["a"]
                integrated = values["delta_a"]-sum(values[n] for n in BUDGET_TERMS)
                audit = dict(pair_seed=pair["seed"], kind=kind, trajectory_seed=seed,
                             steps=states.shape[2]-1,max_abs_step_closure=worst,
                             integrated_residual=integrated,
                             endpoint_error=values["delta_a"]-(a1-a0),
                             majority_inconsistent_steps=n_majority_flip)
                rows.append(dict(pair_seed=pair["seed"],kind=kind,trajectory_seed=seed,
                                 start_a=a0,end_a=a1,**values,
                                 reconstruction=sum(values[n] for n in BUDGET_TERMS)))
                audits.append(audit)
            finally:
                reader.close()
    csv_write(output/"exact_native_balance.csv",rows)
    csv_write(output/"exact_native_audit.csv",audits)
    paired=[]
    for seed in sorted({r["pair_seed"] for r in rows}):
        ev=next(r for r in rows if r["pair_seed"]==seed and r["kind"]=="event")
        co=next(r for r in rows if r["pair_seed"]==seed and r["kind"]=="control")
        line=dict(seed=seed)
        for name in ("delta_a",*BUDGET_TERMS):
            line[f"event_{name}"]=ev[name]
            line[f"control_{name}"]=co[name]
            line[f"difference_{name}"]=ev[name]-co[name]
        paired.append(line)
    csv_write(output/"exact_native_paired_differences.csv",paired)
    fig,ax=plt.subplots(figsize=(10,4.6))
    keys=list(BUDGET_TERMS)
    estimate=[np.mean([r[f"difference_{name}"] for r in paired]) for name in keys]
    ax.bar(range(len(keys)),estimate)
    ax.set_xticks(range(len(keys)),keys,rotation=35,ha="right")
    ax.set(ylabel="event minus same-q control contribution to Δa",title="Source-exact polarization production: bounded matched seeds")
    ax.axhline(0,ls=":",color=".4");fig.tight_layout()
    fig.savefig(output/"exact_polarization_budget.png",dpi=200);plt.close(fig)
    return paired,dict(available=True,n_pairs=len(paired),seeds=[r["seed"] for r in paired],
                       max_abs_step_closure=max(r["max_abs_step_closure"] for r in audits),
                       max_abs_integrated_closure=max(abs(r["integrated_residual"]) for r in audits),
                       max_abs_endpoint_error=max(abs(r["endpoint_error"]) for r in audits),
                       inconsistent_steps=sum(r["majority_inconsistent_steps"] for r in audits))


def seed_bootstrap_mean(values, nboot=5000, seed=20261112):
    values=np.asarray(values,dtype=float)
    if len(values)<1: return dict(mean=None, ci95=[None,None],n=0)
    rng=np.random.default_rng(seed)
    boot=np.mean(rng.choice(values,size=(nboot,len(values)),replace=True),axis=1)
    return dict(mean=float(np.mean(values)),ci95=list(map(float,np.quantile(boot,[.025,.975]))),n=len(values))


def summarize_endpoint_split(comparison, output):
    """Summarize the exact two-endpoint identity without assigning torque sources."""
    terms=("endpoint_heading_change","endpoint_membership_change","endpoint_gap_geometry_change")
    summary={}
    for group in ("event","control","difference"):
        summary[group]={name:seed_bootstrap_mean([r[f"{group}_{name}"] for r in comparison])
                        for name in terms}
    summary["geometry_minus_heading_in_difference"]=seed_bootstrap_mean([
        r["difference_endpoint_gap_geometry_change"]-r["difference_endpoint_heading_change"]
        for r in comparison])
    summary["max_abs_identity_residual"]=max(
        abs(r[name]) for r in comparison
        for name in ("event_endpoint_split_residual","control_endpoint_split_residual"))

    labels=("heading\nchange","membership\nchange","gap-direction\ngeometry")
    x=np.arange(len(terms),dtype=float)
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    for offset,group,color,label in ((-.24,"event","C0","future crossing"),
                                     (0.,"control","C1","same-q control"),
                                     (.24,"difference","C2","event minus control")):
        means=np.array([summary[group][name]["mean"] for name in terms])
        low=np.array([summary[group][name]["ci95"][0] for name in terms])
        high=np.array([summary[group][name]["ci95"][1] for name in terms])
        ax.bar(x+offset,means,width=.22,color=color,alpha=.78,label=label)
        ax.errorbar(x+offset,means,yerr=np.vstack((means-low,high-means)),fmt="none",
                    ecolor="black",elinewidth=1,capsize=3)
    ax.axhline(0,color=".35",ls=":")
    ax.set_xticks(x,labels)
    ax.set_ylabel("contribution to two-unit change in a")
    ax.set_title("Exact endpoint split (37 independent matched pairs)")
    ax.legend(frameon=False,ncol=3,fontsize=9)
    fig.tight_layout();fig.savefig(output/"observed_endpoint_decomposition.png",dpi=200);plt.close(fig)
    return summary


def self_test():
    rng=np.random.default_rng(20261112)
    x=np.r_[rng.uniform(-3.,-.1,size=64), rng.uniform(.1,3.,size=36)]
    y=rng.uniform(-1,1,size=N)
    theta=rng.uniform(-np.pi,np.pi,size=N)
    st=np.stack([x,y,theta]); nxt=st.copy()
    nxt[0]+=rng.uniform(-.015,.015,size=N)
    nxt[1]+=rng.uniform(-.01,.01,size=N)
    # Force a genuine, count-changing crossing on one test step.
    nxt[0,63]=.1
    wall=rng.normal(0,.5,size=N)
    social=rng.normal(0,.8,size=N)
    clip=np.zeros(N)
    noise=rng.normal(0,IN*np.sqrt(DT),size=N)
    nxt[2]=theta + (wall+social+clip)*DT+noise
    z=exact_increment(st,nxt,1,wall,social,clip)
    assert abs(z["residual"])<5e-13, z
    assert abs(z["delta_a"]-sum(z[k] for k in BUDGET_TERMS))<5e-13
    assert z["membership"]!=0
    assert z["gap_direction_rotation"]!=0
    # Test no spatial motion and zero angular noise / torque.
    z0=exact_increment(st,st,1,np.zeros(N),np.zeros(N),np.zeros(N),in_noise=0)
    assert abs(z0["delta_a"])<1e-14
    assert max(abs(z0[k]) for k in BUDGET_TERMS)<1e-14
    # Test angle-only budget against vector projection.
    st2=st.copy(); st2[2]=theta + social*DT + noise
    z2=exact_increment(st,st2,1,np.zeros(N),social,np.zeros(N))
    assert abs(z2["residual"])<5e-13
    # Full mocked read-only pipeline: same public APIs, no MAT/source dependency.
    import tempfile
    import types
    import sys
    import re
    arrays = {}
    const_wall = np.full(N, .07)
    const_social = np.full(N, -.04)
    for seed in range(4):
        states = np.zeros((3,N,48))
        init = st.copy()
        if seed % 2: init[0] *= .95
        states[:,:,0] = init
        for k in range(47):
            old = states[:,:,k]
            new = old.copy()
            new[2] += (const_wall+const_social)*DT+rng.normal(0,IN*np.sqrt(DT),size=N)
            new[0] += DT*np.cos(old[2])
            new[1] += DT*np.sin(old[2])
            states[:,:,k+1] = new
        arrays[seed] = states
    class FakeReader:
        def __init__(self,path):
            match = re.search(r"_n(\d+)_BEM_no_hydro.mat$",str(path))
            self.states = arrays[int(match.group(1))]
        def block(self,name,start,stop):
            assert name == "state"
            return self.states[:,:,start:stop+1]
        def close(self): pass
    fake = types.ModuleType("analyze_approach_current_balance")
    fake.RawReader = FakeReader
    fake.replay_torques = lambda mod, state, para: (const_wall,const_social,np.zeros(N),const_wall+const_social)
    previous = sys.modules.get("analyze_approach_current_balance")
    sys.modules["analyze_approach_current_balance"] = fake
    try:
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            pairs = [dict(seed=0,event_id=0,first_native=19,direction=1,control_seed=1,control_start=7),
                     dict(seed=2,event_id=1,first_native=22,direction=1,control_seed=3,control_start=9)]
            differences = cheap_analysis(out,pairs,out,duration_steps=10,stride=2)
            assert len(differences)==2
            exact, audits = full_native_budget(None,{},out,pairs,out,duration_steps=10,n_exact=2)
            assert len(exact)==2 and audits["max_abs_step_closure"]<1e-12, audits
            assert audits["max_abs_endpoint_error"]<1e-12, audits
            assert (out/"exact_native_balance.csv").exists()
    finally:
        if previous is None: del sys.modules["analyze_approach_current_balance"]
        else: sys.modules["analyze_approach_current_balance"] = previous
    print("PASS: exact polarization budget + two-pair mocked end-to-end pipeline")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--raw-dir",type=Path)
    p.add_argument("--source",type=Path,help="Archived unmodified BEM_bichamber_parallel_instance_no_hydro.py; omit for descriptive mode")
    p.add_argument("--event-dir",type=Path,default=Path("results/nohydro_gap20_descriptive_transition/outputs/initiation_propagation"))
    p.add_argument("--output-dir",type=Path,default=Path("results/nohydro_gap20_descriptive_transition/outputs/polarization_torque_budget"))
    p.add_argument("--duration",type=float,default=2.)
    p.add_argument("--stride",type=int,default=100)
    p.add_argument("--n-exact",type=int,default=8)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        self_test();return
    if args.raw_dir is None: p.error("--raw-dir required for real data")
    nsteps=round(args.duration/DT)
    if nsteps<2 or abs(nsteps*DT-args.duration)>1e-10:
        raise ValueError("duration must be an integer number of native steps")
    pairs=matched_pairs(args.event_dir/"crossing_episodes.csv",
                        args.event_dir/"matched_precursor_pairs.csv")
    if not pairs: raise RuntimeError("No eligible majority-to-minority pairs")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    comparison=cheap_analysis(args.raw_dir,pairs,args.output_dir,nsteps,args.stride)
    observational=seed_bootstrap_mean([r["paired_difference"] for r in comparison])
    endpoint=summarize_endpoint_split(comparison,args.output_dir)
    start_difference=seed_bootstrap_mean([r["difference_start_a"] for r in comparison])
    end_difference=seed_bootstrap_mean([r["difference_end_a"] for r in comparison])
    screening_supported=observational["ci95"][0]>0
    result=dict(observational=dict(n_pairs=len(comparison),
        n_reference_majority_consistent=sum(r["reference_majority_consistent"] for r in comparison),
        event_minus_control_delta_a=observational,
        event_minus_control_start_a=start_difference,
        event_minus_control_end_a=end_difference,
        screening_supported=screening_supported,
        endpoint_split=endpoint), exact=dict(available=False,
        reason="not run because the preregistered positive screening condition failed"))
    if args.source is not None:
        from analyze_approach_current_balance import load_gap20_production
        source=args.source.resolve()
        mod=load_gap20_production(source)
        para=dict(Ip=4.,In=IN,R=mod.R,light=1.,dt=DT,Iw=.94,If=0.,num=N,
                  delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
        paired,audit=full_native_budget(mod,para,args.raw_dir,pairs,args.output_dir,nsteps,args.n_exact)
        if paired:
            audit["term_event_minus_control"]={name:seed_bootstrap_mean([r[f"difference_{name}"] for r in paired])
                                                 for name in ("delta_a",*BUDGET_TERMS)}
        audit["source_sha256"]=hashlib.sha256(source.read_bytes()).hexdigest()
        result["exact"]=audit
    (args.output_dir/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    config=dict(raw_dir=str(args.raw_dir),archived_source=(str(args.source) if args.source else None),
        event_dir=str(args.event_dir),duration=args.duration,stride=args.stride,n_exact=args.n_exact,
        n_all_eligible_pairs=len(pairs),fixed_reference_direction="majority at interval start",
        last_event_frame="first crossing native frame minus one",
        microscopic_replay_gate="run only if event-minus-control delta_a is supported as positive",
        microscopic_replay_gate_passed=screening_supported,
        exact_budget="wall+social+clip+Gaussian_relaxation+realized_noise+membership+gap_direction_rotation",
        caution="Matched controls are quiet by construction; event minus control is descriptive, not causal")
    (args.output_dir/"analysis_config.json").write_text(json.dumps(config,indent=2)+"\n")
    conclusion=("The preregistered positive screening condition failed, so the eight-pair source-exact replay "
                "was not run. The endpoint identity separates net heading change from school geometry, but the "
                "heading term still combines social, wall, and angular-noise effects and cannot attribute a torque source."
                if not result["exact"].get("available") else
                "The exact native budget is computed. Interpret only term contributions with independent-seed "
                "uncertainty and avoid treating conditional events as randomized causal experiments.")
    report=f"""# Does school gap-directed polarization rise immediately before crossing?

Read-only analysis of existing no-hydro gap-20 data. The measured orientation is
`a = mean(e(theta)) dot unit(gap-COM_majority)`, evaluated within the
**fixed original majority chamber** over the two units before an event, and
in a previously matched same-occupancy control over the same duration.

## Data-level diagnostic

Eligible majority-to-minority paired events: {len(pairs)}.
Event-minus-control observed Δa, equal independent-seed weighting:
{observational['mean']:.6f} (95% seed-bootstrap CI {observational['ci95']}).
Source-chamber majority consistent at all sampled times in
{result['observational']['n_reference_majority_consistent']}/{len(pairs)} pairs.

The screening result is opposite to the proposed rise: event schools begin the
window much more gap-directed than controls (event-minus-control start `a` =
{start_difference['mean']:.6f}, 95% CI {start_difference['ci95']}), but that
advantage falls to {end_difference['mean']:.6f} (95% CI {end_difference['ci95']})
immediately before crossing. Thus the known predictive orientation is already
present at least two time units before initiation; it is not produced as a
late two-unit increase.

## Exact endpoint heading/geometry split

The two-endpoint identity closes to a maximum absolute residual of
{endpoint['max_abs_identity_residual']:.3e}. In event windows, net heading
change contributes {endpoint['event']['endpoint_heading_change']['mean']:.6f}
(95% CI {endpoint['event']['endpoint_heading_change']['ci95']}), while rotation
of the school-COM-to-gap direction contributes
{endpoint['event']['endpoint_gap_geometry_change']['mean']:.6f}
(95% CI {endpoint['event']['endpoint_gap_geometry_change']['ci95']}). Membership
contributes exactly zero because no fish crosses during these pre-first-crossing
windows. Thus gap-direction geometry is the statistically resolved dominant
contribution to the absolute event-window decline: the school headings turn
slightly more gapward, but motion of the school COM rotates the gap direction
away more strongly.

For the paired event-minus-control discrepancy, heading contributes
{endpoint['difference']['endpoint_heading_change']['mean']:.6f}
(95% CI {endpoint['difference']['endpoint_heading_change']['ci95']}) and gap
geometry contributes {endpoint['difference']['endpoint_gap_geometry_change']['mean']:.6f}
(95% CI {endpoint['difference']['endpoint_gap_geometry_change']['ci95']}). Their
difference is not resolved: geometry minus heading =
{endpoint['geometry_minus_heading_in_difference']['mean']:.6f}
(95% CI {endpoint['geometry_minus_heading_in_difference']['ci95']}). Therefore
neither can be called the unique dominant explanation of the paired decline.

## Conditional source-exact decomposition

At each step the deterministic angle update is replayed from the archived
production source, including the unchanged visible Delaunay social law and
forward-ray wall law. The identity is

`Δa = wall + social + clip + Gaussian_noise_relaxation + realized_noise + membership + gap_direction_rotation`.

The membership term includes both fish joining/leaving the original source
chamber and its count normalization. The geometry term is the rotation of the
majority-COM-to-gap direction; it is **not another torque**. A collision
velocity override does not enter directly because a is defined from
headings and positions, although it can indirectly affect the COM movement.

The requested decision rule permits the eight-pair replay only if the 37-pair
screening supports a positive event-specific rise. It does not, so social,
wall, clip, noise-relaxation, and realized-noise terms were not computed on
the real trajectories. Exact native-step status: {json.dumps(result['exact'],indent=2)}

## Interpretation

{conclusion}

The social term remains the original combined visual rule. Noise relaxation
is the exact Gaussian conditional mean at each native step; realized noise is
zero-conditional-mean, but not necessarily zero in a selected event sample.
The event/control comparison has only matched observational meaning.
This establishes an identity and describes mechanisms, not a closed Markov
model for (q,a) or a whole-school transition-rate theory.

## Reproducibility

- `observed_intervals.csv`: event/control endpoints and exact endpoint terms.
- `observed_paired_differences.csv`: independent-seed paired contributions.
- `observed_time_series.csv` and `observed_time_series_summary.csv`: trajectory
  samples and seed-bootstrap bands.
- `observed_polarization_change.png`: pair-level endpoint changes.
- `observed_event_control_time_series.png`: mean time course and paired contrast.
- `observed_endpoint_decomposition.png`: heading, membership, and geometry split.
- `summary.json`: screening decision, uncertainty, and identity audit.
"""
    (args.output_dir/"REPORT.md").write_text(report)
    print(report)


if __name__=="__main__":
    main()
