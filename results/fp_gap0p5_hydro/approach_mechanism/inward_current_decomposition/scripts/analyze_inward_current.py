#!/usr/bin/env python3
"""Source-resolved audit of pre-crossing inward-current generation.

This is deliberately a bounded replay: one matched event/control pair per
selected seed and one full microscopic replay every ``replay_stride`` native
steps.  Native saved states are used for the exact finite-step Gaussian and
translational identities.  No coefficient is fitted.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

HERE = Path(__file__).resolve()
MECH = HERE.parents[2]
REPLAY_DIR = MECH / "pre_crossing_turning" / "scripts"
sys.path.insert(0, str(REPLAY_DIR))
from replay_dynamics import Parameters, analytic_geometry, replay_snapshot  # noqa: E402

DT = 1e-3
DISCONTINUITY = {10: 600000, 11: 600000, 13: 600000,
                 15: 600000, 17: 600000, 19: 600000}
SOURCE_NAMES = ("wall", "alignment", "attraction",
                "direct_hydro_effective", "boundary_hydro_effective")


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fields, lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def rng_noise_for_indices(seed, indices, nfish=100, chunk=20000):
    """Replay the production RandomState stream without retaining unused draws."""
    want = np.asarray(sorted(set(int(i) for i in indices)), dtype=int)
    if len(want) == 0:
        return {}
    rng = np.random.RandomState(seed)
    # Production initialization consumes three uniform vectors before the
    # per-step Gaussian stream (alpha, radial position, initial heading).
    rng.uniform(-np.pi, np.pi, nfish)
    rng.uniform(0.2, 1, nfish)
    rng.uniform(-1, 1, nfish)
    result = {}
    cursor = 0
    target_pos = 0
    last = int(want[-1])
    while cursor <= last:
        count = min(chunk, last-cursor+1)
        block = rng.normal(size=(count, nfish))
        while target_pos < len(want) and want[target_pos] < cursor+count:
            idx = int(want[target_pos])
            result[idx] = block[idx-cursor].copy()
            target_pos += 1
        cursor += count
    return result


def select_pairs(rows, seeds):
    """Select one high-quality matched pair per seed without outcome features."""
    selected = []
    for seed in seeds:
        sr = [r for r in rows if int(r["seed"]) == seed]
        by_pair = {}
        for r in sr:
            by_pair.setdefault(int(r["pair_id"]), []).append(r)
        candidates = []
        for pair, rr in by_pair.items():
            if {r["condition"] for r in rr} != {"event", "control"}:
                continue
            # Exact RNG replay is only valid before the known production restart.
            cutoff = DISCONTINUITY.get(seed, np.inf)
            if any(int(round((float(r["risk_time"])+1.0)/DT)) >= cutoff for r in rr):
                continue
            control = next(r for r in rr if r["condition"] == "control")
            candidates.append((float(control["match_distance"]), pair, rr))
        if not candidates:
            raise RuntimeError(f"no RNG-valid matched pair for seed {seed}")
        candidates.sort(key=lambda x: (x[0], x[1]))
        selected.extend(candidates[0][2])
    return selected


def load_boundary(path):
    g = sio.loadmat(path)
    normal = np.asarray(g["normal_vec"])
    return dict(x=np.asarray(g["Sx"]).reshape(-1),
                y=np.asarray(g["Sy"]).reshape(-1),
                weights=np.asarray(g["weights"]).reshape(-1),
                nx=normal[:, 0], ny=normal[:, 1])


def rolling_leader_decomposition(state, i0, i1, origin_sign, stride):
    """Exactly split rolling closest-ten current changes into turn/replacement."""
    inds = list(range(i0, i1+1, stride))
    if inds[-1] != i1:
        inds.append(i1)
    ids = []
    js = []
    for idx in inds:
        # Positive d means inside the origin chamber; smaller is closer divider.
        d = -origin_sign*state[0, :, idx]
        origin = np.flatnonzero(d >= 0)
        pool = origin if len(origin) >= 10 else np.arange(state.shape[1])
        lead = pool[np.argsort(d[pool])[:10]]
        ids.append(set(int(x) for x in lead))
        js.append(origin_sign*np.mean(np.cos(state[2, lead, idx])))
    turning = replacement = 0.0
    turnover = []
    for a, b, ia, ib in zip(ids[:-1], ids[1:], inds[:-1], inds[1:]):
        common = a & b; enter = b-a; leave = a-b
        turning += origin_sign/10*sum(np.cos(state[2, j, ib])-np.cos(state[2, j, ia])
                                      for j in common)
        replacement += origin_sign/10*(sum(np.cos(state[2, j, ib]) for j in enter)
                                        - sum(np.cos(state[2, j, ia]) for j in leave))
        turnover.append(len(enter)/10)
    return dict(leader_current_initial=js[0], leader_current_final=js[-1],
                leader_current_change=js[-1]-js[0],
                leader_turning_contribution=turning,
                leader_replacement_contribution=replacement,
                leader_split_residual=js[-1]-js[0]-turning-replacement,
                mean_leader_turnover=float(np.mean(turnover)),
                unique_leader_identities=len(set().union(*ids)))


def interp_integral(sample_t, sample_y, native_t, weight=None):
    y = np.interp(native_t, sample_t, sample_y)
    if weight is not None:
        y = y*weight
    return DT*np.sum(y)


def analyze_segment(row, state, noise, boundary, lu, piv, p, geom, stride):
    seed = int(row["seed"]); condition = row["condition"]
    origin = row["origin"]; sign = 1.0 if origin == "left" else -1.0
    i0 = int(round(float(row["risk_time"])/DT)); steps = int(round(1.0/DT)); i1 = i0+steps
    theta = state[2, :, i0:i1+1]
    js = sign*np.mean(np.cos(theta), axis=0)
    z = np.stack([noise[i] for i in range(i0, i1)], axis=1)
    dtheta = np.diff(theta, axis=1)
    omega = (dtheta-p.In*np.sqrt(DT)*z)/DT
    decay = np.exp(-0.5*p.In**2*DT)
    deterministic_heading = theta[:, :-1]+omega*DT
    conditional_next = sign*np.mean(decay*np.cos(deterministic_heading), axis=0)
    deterministic_steering = sign*np.mean(np.cos(deterministic_heading)-np.cos(theta[:, :-1]), axis=0)
    noise_decay = sign*np.mean((decay-1.0)*np.cos(deterministic_heading), axis=0)
    stochastic = js[1:]-conditional_next
    step_residual = np.diff(js)-deterministic_steering-noise_decay-stochastic

    # Exact effective translational current inferred from the actual Euler position update.
    vx = np.diff(state[0, :, i0:i1+1], axis=1)/DT
    jh = sign*np.mean(vx-np.cos(theta[:, :-1]), axis=0)
    x_in = sign*(np.mean(state[0, :, i1])-np.mean(state[0, :, i0]))
    self_displacement = DT*np.sum(js[:-1])
    hydro_displacement = DT*np.sum(jh)

    replay_inds = list(range(i0, i1, stride))
    source_rows = []
    for idx in replay_inds:
        rep = replay_snapshot(state[:, :, idx], state[:, :, idx+1], boundary, lu, piv, p, geom)
        th = state[2, :, idx]
        terms = rep["terms"]
        rec = dict(index=idx, time=(idx-i0)*DT)
        omega_sum = np.zeros(state.shape[1])
        for name in SOURCE_NAMES:
            om = np.asarray(terms[name])
            omega_sum += om
            rec[name] = -sign*np.mean(np.sin(th)*om)
        omega_clip = rep["deterministic_total"]-omega_sum
        rec["clipping_override"] = -sign*np.mean(np.sin(th)*omega_clip)
        rec["total_replay"] = -sign*np.mean(np.sin(th)*rep["deterministic_total"])
        rec["omega_replay_rms_error"] = float(np.sqrt(np.mean((rep["deterministic_total"]-omega[:, idx-i0])**2)))
        rec["direct_translation"] = sign*np.mean(rep["direct_velocity_effective"][0])
        rec["boundary_translation"] = sign*np.mean(rep["boundary_velocity_effective"][0])
        rec["kinematic_hydro"] = jh[idx-i0]
        source_rows.append(rec)

    native_t = np.arange(steps)*DT
    # Discrete summation-by-parts weight: the last angular increment cannot
    # affect a velocity already used within this interval.
    weight = (steps-1-np.arange(steps))*DT
    st = np.array([r["time"] for r in source_rows])
    contributions = {}
    for name in SOURCE_NAMES+("clipping_override",):
        sy = np.array([r[name] for r in source_rows])
        contributions[name] = interp_integral(st, sy, native_t, weight)
    exact_steering = float(np.sum(weight*deterministic_steering))
    exact_noise_decay = float(np.sum(weight*noise_decay))
    exact_stochastic = float(np.sum(weight*stochastic))
    inherited = float(steps*DT*js[0])
    ito_continuous = float(-0.5*p.In**2*DT*np.sum(weight*js[:-1]))
    source_sum = sum(contributions.values())
    direct_trans = interp_integral(st, np.array([r["direct_translation"] for r in source_rows]), native_t)
    boundary_trans = interp_integral(st, np.array([r["boundary_translation"] for r in source_rows]), native_t)

    result = dict(seed=seed, pair_id=int(row["pair_id"]), condition=condition,
                  origin=origin, risk_time=float(row["risk_time"]), replay_stride=stride,
                  actual_inward_displacement=x_in,
                  self_propulsion_displacement=self_displacement,
                  translational_hydrodynamics=hydro_displacement,
                  kinematic_residual=x_in-self_displacement-hydro_displacement,
                  inherited_current_contribution=inherited,
                  deterministic_steering_exact=exact_steering,
                  angular_noise_decay_exact=exact_noise_decay,
                  ito_noise_decay_continuous=ito_continuous,
                  realized_stochastic_residual=exact_stochastic,
                  exact_displacement_reconstruction_residual=(self_displacement-inherited-exact_steering-exact_noise_decay-exact_stochastic),
                  source_sum_sampled=source_sum,
                  finite_step_and_quadrature_remainder=exact_steering-source_sum,
                  direct_translation_sampled=direct_trans,
                  boundary_translation_sampled=boundary_trans,
                  translation_source_remainder=hydro_displacement-direct_trans-boundary_trans,
                  max_finite_step_identity_residual=float(np.max(np.abs(step_residual))),
                  rms_replayed_omega_error=float(np.sqrt(np.mean([r["omega_replay_rms_error"]**2 for r in source_rows]))),
                  initial_Js=js[0], final_Js=js[-1], delta_Js=js[-1]-js[0])
    for name, value in contributions.items():
        result[f"source_{name}_displacement"] = value
    # Unweighted source impulses directly answer whether production precedes approach.
    for name in SOURCE_NAMES+("clipping_override",):
        sy = np.array([r[name] for r in source_rows])
        result[f"impulse_{name}"] = interp_integral(st, sy, native_t)
        early = native_t < .5
        result[f"early_impulse_{name}"] = interp_integral(st, sy, native_t[early])
    result.update(rolling_leader_decomposition(state, i0, i1, sign, stride))
    for r in source_rows:
        r.update(seed=seed, pair_id=int(row["pair_id"]), condition=condition, origin=origin)
    return result, source_rows


def paired_summary(rows, rng_seed=20261023, boots=20000):
    metrics = [k for k in rows[0] if k not in {"seed", "pair_id", "condition", "origin", "risk_time", "replay_stride"}]
    by = {}
    for r in rows:
        by.setdefault(int(r["seed"]), {})[r["condition"]] = r
    seeds = sorted(s for s, x in by.items() if {"event", "control"} <= set(x))
    rng = np.random.default_rng(rng_seed)
    out = []
    for metric in metrics:
        try:
            ev = np.array([float(by[s]["event"][metric]) for s in seeds])
            co = np.array([float(by[s]["control"][metric]) for s in seeds])
        except (ValueError, TypeError):
            continue
        diff = ev-co
        bd = np.mean(diff[rng.integers(0, len(diff), size=(boots, len(diff)))], axis=1)
        out.append(dict(metric=metric, n_seeds=len(seeds), event_mean=np.mean(ev),
                        control_mean=np.mean(co), paired_difference=np.mean(diff),
                        ci_low=np.quantile(bd, .025), ci_high=np.quantile(bd, .975)))
    return out


def seed_paired_differences(rows):
    by = {}
    for r in rows:
        by.setdefault(int(r["seed"]), {})[r["condition"]] = r
    omit = {"seed", "pair_id", "condition", "origin", "risk_time", "replay_stride"}
    metrics = [k for k in rows[0] if k not in omit]
    out = []
    for seed, rr in sorted(by.items()):
        if {"event", "control"} > set(rr):
            continue
        row = dict(seed=seed, event_pair_id=rr["event"]["pair_id"],
                   control_pair_id=rr["control"]["pair_id"])
        for metric in metrics:
            try:
                row[metric] = float(rr["event"][metric])-float(rr["control"][metric])
            except (ValueError, TypeError):
                pass
        out.append(row)
    return out


def plots(rows, paired, out):
    out.mkdir(parents=True, exist_ok=True)
    d = {r["metric"]: r for r in paired}
    panels = [
      (["actual_inward_displacement", "self_propulsion_displacement",
        "inherited_current_contribution", "deterministic_steering_exact",
        "angular_noise_decay_exact", "realized_stochastic_residual",
        "translational_hydrodynamics"],
       ["actual", "self propulsion", "inherited", "exact steering",
        "noise decay", "realized noise", "hydro translation"], "Exact additive audit"),
      (["source_wall_displacement", "source_alignment_displacement",
        "source_attraction_displacement", "source_direct_hydro_effective_displacement",
        "source_boundary_hydro_effective_displacement", "source_clipping_override_displacement",
        "finite_step_and_quadrature_remainder"],
       ["wall", "alignment", "attraction", "direct angular", "BEM angular",
        "clipping / override", "finite-step correction"], "Steering-source resolution")]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True)
    for ax, (metrics, labels, title) in zip(axs, panels):
        y = np.arange(len(metrics)); val = np.array([float(d[m]["paired_difference"]) for m in metrics])
        lo = np.array([float(d[m]["ci_low"]) for m in metrics]); hi = np.array([float(d[m]["ci_high"]) for m in metrics])
        ax.axvline(0, color=".25", lw=1)
        ax.errorbar(val, y, xerr=np.vstack((val-lo, hi-val)), fmt="o", capsize=3, color="#1f77b4")
        ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_title(title)
        ax.set_xlabel("event - matched control displacement")
    fig.suptitle("Pre-crossing inward-current decomposition, A-2 to A-1")
    fig.tight_layout(); fig.savefig(out/"01_displacement_contributions.png", dpi=220); fig.savefig(out/"01_displacement_contributions.pdf"); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    for cond, color in (("event", "#d62728"), ("control", "#1f77b4")):
        rr = [r for r in rows if r["condition"] == cond]
        axs[0].scatter([float(r["leader_turning_contribution"]) for r in rr],
                       [float(r["leader_replacement_contribution"]) for r in rr], label=cond, color=color)
        axs[1].scatter([float(r["inherited_current_contribution"]) for r in rr],
                       [float(r["actual_inward_displacement"]) for r in rr], label=cond, color=color)
    axs[0].axhline(0,color=".6",lw=.8); axs[0].axvline(0,color=".6",lw=.8)
    axs[0].set(xlabel="turning of retained leaders", ylabel="identity replacement", title="Rolling closest ten")
    axs[1].axhline(0,color=".6",lw=.8); axs[1].axvline(0,color=".6",lw=.8)
    axs[1].set(xlabel="inherited-current contribution", ylabel="actual inward displacement", title="Inherited approach")
    for ax in axs: ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out/"02_inheritance_and_leader_renewal.png", dpi=220); fig.savefig(out/"02_inheritance_and_leader_renewal.pdf"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--risk-set", type=Path, required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--lu", type=Path, required=True)
    ap.add_argument("--piv", type=Path, required=True)
    ap.add_argument("--seeds", type=int, nargs="*", default=list(range(20)))
    ap.add_argument("--replay-stride", type=int, default=10)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(); out = args.out.resolve()
    (out/"tables"/"by_seed").mkdir(parents=True, exist_ok=True)
    p = Parameters(); geom = analytic_geometry(p); boundary = load_boundary(args.geometry)
    lu = np.load(args.lu, mmap_mode="r"); piv = np.load(args.piv)
    selected = select_pairs(read_csv(args.risk_set), args.seeds)
    write_csv(out/"tables"/"selected_risk_set.csv", selected)
    all_summary = []; all_sources = []
    for seed in args.seeds:
        target = out/"tables"/"by_seed"/f"seed_{seed:02d}_summary.csv"
        starget = out/"tables"/"by_seed"/f"seed_{seed:02d}_sources.csv"
        if target.exists() and starget.exists() and not args.force:
            all_summary.extend(read_csv(target)); all_sources.extend(read_csv(starget)); print(f"seed {seed}: checkpoint", flush=True); continue
        raw = args.raw_dir/f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state = np.asarray(sio.loadmat(raw, variable_names=["state"], squeeze_me=True)["state"])
        sr = [r for r in selected if int(r["seed"]) == seed]
        required = []
        for r in sr:
            i0 = int(round(float(r["risk_time"])/DT)); required.extend(range(i0, i0+1000))
        noise = rng_noise_for_indices(seed, required, state.shape[1])
        summaries = []; sources = []
        for r in sr:
            a, b = analyze_segment(r, state, noise, boundary, lu, piv, p, geom, args.replay_stride)
            summaries.append(a); sources.extend(b)
        write_csv(target, summaries); write_csv(starget, sources)
        all_summary.extend(summaries); all_sources.extend(sources); print(f"seed {seed}: done", flush=True)
    write_csv(out/"tables"/"segment_decomposition.csv", all_summary)
    write_csv(out/"tables"/"sampled_source_histories.csv", all_sources)
    paired = paired_summary(all_summary)
    write_csv(out/"tables"/"paired_event_control_summary.csv", paired)
    write_csv(out/"tables"/"seed_paired_differences.csv", seed_paired_differences(all_summary))
    headline = {"actual_inward_displacement", "self_propulsion_displacement",
                "translational_hydrodynamics", "inherited_current_contribution",
                "deterministic_steering_exact", "angular_noise_decay_exact",
                "realized_stochastic_residual", "source_wall_displacement",
                "source_alignment_displacement", "source_attraction_displacement",
                "source_direct_hydro_effective_displacement",
                "source_boundary_hydro_effective_displacement",
                "source_clipping_override_displacement",
                "leader_turning_contribution", "leader_replacement_contribution"}
    write_csv(out/"summary.csv", [r for r in paired if r["metric"] in headline])
    plots(all_summary, paired, out/"figures")
    config = dict(seeds=args.seeds, replay_stride=args.replay_stride, native_dt=DT,
                  interval="risk time A-2 through A-1", pair_selection="minimum distance-match cost per seed, no outcome features",
                  geometry_convention="archived gap=R; gap/R=1; gap/(2R)=0.5",
                  fitted_coefficients=False, new_simulations=False)
    (out/"analysis_config.json").write_text(json.dumps(config, indent=2)+"\n")


if __name__ == "__main__":
    main()
