#!/usr/bin/env python3
"""Source-exact tests of the single-fish circular wall-avoidance theory.

The production update reproduced here is the N=1, wall-only specialization of
``BEM_circular_parallel_instance.py``.  In particular, this retains both the
forward-ray torque regularization and the separate positional sticky safeguard.
No coefficient is fitted to the trajectories.
"""
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
from numba import njit
from scipy.linalg import expm


RADIUS = 1.0
SAFE_DISTANCE = 0.1
HIGH_ORDER = 6
PRODUCTION_DT = 1.0e-3
WALL_CLIP = math.pi / PRODUCTION_DT  # hard-coded in GetWallAvoid


def orbit_parameters(iw: float, radius: float = RADIUS) -> dict[str, float]:
    s = math.sqrt(1.0 + iw * iw)
    return {
        "r_star": radius / s,
        "d_star": radius * iw / s,
        "clearance_star": radius * (1.0 - 1.0 / s),
        "omega_star": s / radius,
        "gamma": s / (radius * iw),
        "coupling": s**4 / (radius**2 * iw**2),
        "decay": s / (2.0 * radius * iw),
        "oscillation": math.sqrt((1.0 + iw * iw) * (4.0 * iw * iw + 3.0))
        / (2.0 * radius * iw),
    }


@njit(cache=True)
def wrap_angle(angle: float) -> float:
    return math.atan2(math.sin(angle), math.cos(angle))


@njit(cache=True)
def source_step(
    x: float,
    y: float,
    theta: float,
    previous_proposed_signed_distance: float,
    iw: float,
    inoise: float,
    dt: float,
    gaussian: float,
):
    """One Euler--Maruyama step matching the circular production safeguards."""
    ct = math.cos(theta)
    st = math.sin(theta)
    rr = x * x + y * y
    projection = x * ct + y * st
    disc = projection * projection + RADIUS * RADIUS - rr
    if disc < 0.0:
        disc = 0.0
    distance = -projection + math.sqrt(disc)
    hit_x = x + distance * ct
    hit_y = y + distance * st
    wall_angle = math.atan2(hit_y, hit_x)
    incidence = wrap_angle(theta - wall_angle)
    raw = 0.0
    if incidence > 0.0:
        raw = iw / distance
    elif incidence < 0.0:
        raw = -iw / distance
    correction = distance <= SAFE_DISTANCE
    if correction:
        raw *= 1.0 + 1.0 / distance ** (HIGH_ORDER - 1)
    wall_clipped = min(max(raw, -WALL_CLIP), WALL_CLIP)
    total_clip = math.pi / dt
    omega = min(max(wall_clipped, -total_clip), total_clip)
    clipped = omega != raw

    proposed_x = x + ct * dt
    proposed_y = y + st * dt
    proposed_signed_distance = math.sqrt(proposed_x**2 + proposed_y**2) - RADIUS
    freeze = False
    if proposed_signed_distance + SAFE_DISTANCE >= 0.0:
        if math.isnan(previous_proposed_signed_distance):
            freeze = True
        elif previous_proposed_signed_distance < proposed_signed_distance:
            freeze = True
    if freeze:
        next_x = x
        next_y = y
    else:
        next_x = proposed_x
        next_y = proposed_y
    next_theta = theta + omega * dt + inoise * math.sqrt(dt) * gaussian
    return (
        next_x,
        next_y,
        next_theta,
        proposed_signed_distance,
        distance,
        omega,
        correction,
        clipped,
        freeze,
    )


@njit(cache=True)
def simulate_paths(
    initial: np.ndarray,
    iw: float,
    inoise: float,
    dt: float,
    n_steps: int,
    stride: int,
    random_numbers: np.ndarray,
):
    n = initial.shape[0]
    n_save = n_steps // stride + 1
    states = np.empty((n_save, n, 3), dtype=np.float64)
    diagnostics = np.zeros((n, 3), dtype=np.int64)  # correction, clipping, freeze
    x = initial[:, 0].copy()
    y = initial[:, 1].copy()
    theta = initial[:, 2].copy()
    previous = np.empty(n)
    previous[:] = np.nan
    states[0, :, :] = initial
    save_index = 1
    for step in range(n_steps):
        for j in range(n):
            out = source_step(
                x[j], y[j], theta[j], previous[j], iw, inoise, dt,
                random_numbers[step, j],
            )
            x[j], y[j], theta[j], previous[j] = out[0], out[1], out[2], out[3]
            diagnostics[j, 0] += int(out[6])
            diagnostics[j, 1] += int(out[7])
            diagnostics[j, 2] += int(out[8])
        if (step + 1) % stride == 0:
            states[save_index, :, 0] = x
            states[save_index, :, 1] = y
            states[save_index, :, 2] = theta
            save_index += 1
    return states, diagnostics


@njit(cache=True)
def stochastic_statistics(
    initial: np.ndarray,
    intended_chirality: np.ndarray,
    iw: float,
    inoise: float,
    dt: float,
    burn_steps: int,
    sample_steps: int,
    seed: int,
):
    """Per-trajectory local statistics; contamination is reported, not hidden."""
    np.random.seed(seed)
    n = initial.shape[0]
    x = initial[:, 0].copy()
    y = initial[:, 1].copy()
    theta = initial[:, 2].copy()
    previous = np.empty(n)
    previous[:] = np.nan
    alive = np.ones(n, dtype=np.bool_)
    count = np.zeros(n, dtype=np.int64)
    sums = np.zeros((n, 5))  # u, v, u2, v2, uv
    flags = np.zeros((n, 4), dtype=np.int64)  # correction, clip, freeze, switch
    pars = orbit_parameters_numba(iw)
    r_star = pars[0]
    total = burn_steps + sample_steps
    for step in range(total):
        for j in range(n):
            z = np.random.normal()
            out = source_step(x[j], y[j], theta[j], previous[j], iw, inoise, dt, z)
            x[j], y[j], theta[j], previous[j] = out[0], out[1], out[2], out[3]
            flags[j, 0] += int(out[6])
            flags[j, 1] += int(out[7])
            flags[j, 2] += int(out[8])
            phi = math.atan2(y[j], x[j])
            alpha = wrap_angle(theta[j] - phi)
            current_chirality = 1.0 if math.sin(alpha) >= 0.0 else -1.0
            if current_chirality != intended_chirality[j]:
                flags[j, 3] += 1
                alive[j] = False
            if out[6] or out[7] or out[8]:
                alive[j] = False
            if step >= burn_steps and alive[j]:
                u = math.sqrt(x[j] * x[j] + y[j] * y[j]) - r_star
                v = wrap_angle(alpha - intended_chirality[j] * math.pi / 2.0)
                sums[j, 0] += u
                sums[j, 1] += v
                sums[j, 2] += u * u
                sums[j, 3] += v * v
                sums[j, 4] += u * v
                count[j] += 1
    return count, sums, flags


@njit(cache=True)
def orbit_parameters_numba(iw: float):
    s = math.sqrt(1.0 + iw * iw)
    return np.array([RADIUS / s, RADIUS * iw / s, s / RADIUS])


def make_initial(r: float, alpha: float, phi: float = 0.0) -> np.ndarray:
    return np.array([r * math.cos(phi), r * math.sin(phi), phi + alpha])


def trajectory_coordinates(states: np.ndarray, chirality: float, iw: float):
    x, y, theta = states[:, 0], states[:, 1], states[:, 2]
    r = np.hypot(x, y)
    phi = np.unwrap(np.arctan2(y, x))
    alpha = np.arctan2(np.sin(theta - np.arctan2(y, x)), np.cos(theta - np.arctan2(y, x)))
    u = r - orbit_parameters(iw)["r_star"]
    v = np.arctan2(
        np.sin(alpha - chirality * np.pi / 2),
        np.cos(alpha - chirality * np.pi / 2),
    )
    return r, phi, alpha, u, v


def deterministic_validation(outdir: Path):
    strengths = (0.25, 0.94, 2.0)
    dts = (0.002, 0.001, 0.0005)
    rows: list[dict[str, object]] = []
    traces = {}
    rng = np.random.default_rng(20261104)
    for iw in strengths:
        p = orbit_parameters(iw)
        initial = []
        metadata = []
        for chirality in (-1, 1):
            for radial_factor, angular_offset in ((0.85, 0.15), (1.05, -0.12), (0.7, 0.45), (0.45, -0.65)):
                rr = min(0.88, radial_factor * p["r_star"])
                initial.append(make_initial(rr, chirality * np.pi / 2 + angular_offset))
                metadata.append((chirality, "broad" if abs(angular_offset) > 0.2 else "local"))
        initial = np.asarray(initial)
        for dt in dts:
            t_end = 80.0
            stride = max(1, int(round(0.01 / dt)))
            n_steps = int(round(t_end / dt))
            random = np.zeros((n_steps, len(initial)))
            states, diag = simulate_paths(initial, iw, 0.0, dt, n_steps, stride, random)
            times = np.arange(states.shape[0]) * stride * dt
            for j, (chirality, scope) in enumerate(metadata):
                r, phi, alpha, _, _ = trajectory_coordinates(states[:, j, :], chirality, iw)
                tail = times >= 60.0
                slope = np.polyfit(times[tail], phi[tail], 1)[0]
                final_chirality = 1 if np.mean(np.sin(alpha[tail])) >= 0 else -1
                rows.append({
                    "Iw": iw, "dt": dt, "initial_chirality": chirality,
                    "initial_scope": scope, "final_chirality": final_chirality,
                    "theory_radius": p["r_star"], "measured_radius": float(np.mean(r[tail])),
                    "radius_error": float(np.mean(r[tail]) - p["r_star"]),
                    "theory_frequency": chirality * p["omega_star"],
                    "measured_frequency": float(slope),
                    "frequency_error": float(slope - chirality * p["omega_star"]),
                    "correction_fraction": diag[j, 0] / n_steps,
                    "clip_fraction": diag[j, 1] / n_steps,
                    "freeze_fraction": diag[j, 2] / n_steps,
                })
            if dt == 0.001:
                traces[iw] = (times, states[:, 1, :], diag[1])
    write_csv(outdir / "tables" / "deterministic_orbits.csv", rows)
    return rows, traces


def stability_validation(outdir: Path):
    rows = []
    traces = {}
    for iw in (0.25, 0.94, 2.0):
        p = orbit_parameters(iw)
        for chirality in (-1, 1):
            sgn = float(chirality)
            A = np.array([[0.0, -sgn], [sgn * p["coupling"], -p["gamma"]]])
            vals, vecs = np.linalg.eig(A)
            idx = int(np.argmax(np.imag(vals)))
            for perturbation, delta in (("radial", np.array([0.01, 0.0])), ("angular", np.array([0.0, 0.01]))):
                dt = 0.0005
                stride = 4
                t_end = 12.0
                initial = make_initial(p["r_star"] + delta[0], chirality * np.pi / 2 + delta[1])[None, :]
                random = np.zeros((int(t_end / dt), 1))
                states, diag = simulate_paths(initial, iw, 0.0, dt, int(t_end / dt), stride, random)
                times = np.arange(states.shape[0]) * stride * dt
                _, _, _, u, v = trajectory_coordinates(states[:, 0, :], chirality, iw)
                observed = np.column_stack([u, v])
                predicted = np.vstack([expm(A * t) @ delta for t in times])
                # Compare only while the imposed perturbation remains well above
                # the O(dt) displacement of the Euler map's discrete orbit.
                local_window = times <= min(4.0, 2.0 / p["decay"])
                scale = np.sqrt(np.mean(predicted[local_window]**2))
                normalized_rmse = float(
                    np.sqrt(np.mean((observed[local_window] - predicted[local_window]) ** 2)) / scale
                )
                inv = np.linalg.inv(vecs)
                modal = (inv @ observed.T)[idx]
                use = (
                    (times >= 0.05)
                    & (times <= min(4.0, 2.0 / p["decay"]))
                    & (np.abs(modal) > 1e-10)
                )
                measured_decay = -float(np.polyfit(times[use], np.log(np.abs(modal[use])), 1)[0])
                measured_osc = abs(float(np.polyfit(times[use], np.unwrap(np.angle(modal[use])), 1)[0]))
                rows.append({
                    "Iw": iw, "chirality": chirality, "perturbation": perturbation,
                    "theory_decay": p["decay"], "measured_decay": measured_decay,
                    "theory_oscillation": p["oscillation"], "measured_oscillation": measured_osc,
                    "normalized_linear_response_rmse": normalized_rmse,
                    "correction_steps": int(diag[0, 0]), "freeze_steps": int(diag[0, 2]),
                })
                if chirality == 1 and perturbation == "radial":
                    traces[iw] = (times, observed, predicted)
    write_csv(outdir / "tables" / "linear_stability.csv", rows)
    return rows, traces


def stochastic_validation(outdir: Path, quick: bool = False):
    strengths = (0.25, 0.94, 2.0)
    noises = (0.03, 0.05, 0.10)
    n_per_branch = 16 if quick else 48
    burn_time = 40.0 if quick else 100.0
    sample_time = 80.0 if quick else 300.0
    rows = []
    for iw in strengths:
        p = orbit_parameters(iw)
        for inoise in noises:
            for dt in (0.001,):
                chirality = np.repeat(np.array([-1.0, 1.0]), n_per_branch)
                initial = np.vstack([make_initial(p["r_star"], c * np.pi / 2) for c in chirality])
                count, sums, flags = stochastic_statistics(
                    initial, chirality, iw, inoise, dt,
                    int(burn_time / dt), int(sample_time / dt),
                    710000 + int(iw * 1000) + int(inoise * 10000),
                )
                for branch in (-1, 1):
                    subset = np.where((chirality == branch) & (count > 1))[0]
                    total_branch = np.where(chirality == branch)[0]
                    if len(subset):
                        means = sums[subset, :2] / count[subset, None]
                        second = sums[subset, 2:4] / count[subset, None]
                        variances = second - means**2
                        cov = sums[subset, 4] / count[subset] - means[:, 0] * means[:, 1]
                        var_u = float(np.sum(sums[subset, 2]) / np.sum(count[subset]) - (np.sum(sums[subset, 0]) / np.sum(count[subset]))**2)
                        var_v = float(np.sum(sums[subset, 3]) / np.sum(count[subset]) - (np.sum(sums[subset, 1]) / np.sum(count[subset]))**2)
                        cov_uv = float(np.sum(sums[subset, 4]) / np.sum(count[subset]) - (np.sum(sums[subset, 0]) * np.sum(sums[subset, 1]) / np.sum(count[subset])**2))
                        seed_se_u = float(np.std(variances[:, 0], ddof=1) / np.sqrt(len(subset))) if len(subset) > 1 else float("nan")
                        seed_se_v = float(np.std(variances[:, 1], ddof=1) / np.sqrt(len(subset))) if len(subset) > 1 else float("nan")
                    else:
                        var_u = var_v = cov_uv = seed_se_u = seed_se_v = float("nan")
                    theory_v = inoise**2 * RADIUS * iw / (2 * math.sqrt(1 + iw**2))
                    theory_u = inoise**2 * RADIUS**3 * iw**3 / (2 * (1 + iw**2) ** 2.5)
                    rows.append({
                        "Iw": iw, "In": inoise, "dt": dt, "branch": branch,
                        "n_trajectories": len(total_branch), "n_with_clean_samples": len(subset),
                        "clean_sample_fraction": float(np.sum(count[subset]) / (len(total_branch) * int(sample_time / dt))) if len(subset) else 0.0,
                        "radial_variance": var_u, "radial_seed_se": seed_se_u,
                        "radial_theory": theory_u, "radial_ratio": var_u / theory_u if theory_u else float("nan"),
                        "angular_variance": var_v, "angular_seed_se": seed_se_v,
                        "angular_theory": theory_v, "angular_ratio": var_v / theory_v if theory_v else float("nan"),
                        "covariance": cov_uv,
                        "trajectories_with_contact": int(np.sum(np.any(flags[total_branch][:, [0, 2]] > 0, axis=1))),
                        "trajectories_with_clipping": int(np.sum(flags[total_branch, 1] > 0)),
                        "trajectories_with_switch": int(np.sum(flags[total_branch, 3] > 0)),
                    })
    # Focused dt check in the valid intermediate case.
    iw, inoise = 0.94, 0.05
    p = orbit_parameters(iw)
    for dt in (0.002, 0.0005):
        chirality = np.repeat(np.array([-1.0, 1.0]), n_per_branch)
        initial = np.vstack([make_initial(p["r_star"], c * np.pi / 2) for c in chirality])
        count, sums, flags = stochastic_statistics(
            initial, chirality, iw, inoise, dt,
            int(burn_time / dt), int(sample_time / dt), 880000 + int(dt * 1e7),
        )
        subset = np.where(count > 1)[0]
        var_u = float(np.sum(sums[subset, 2]) / np.sum(count[subset]) - (np.sum(sums[subset, 0]) / np.sum(count[subset]))**2)
        var_v = float(np.sum(sums[subset, 3]) / np.sum(count[subset]) - (np.sum(sums[subset, 1]) / np.sum(count[subset]))**2)
        rows.append({
            "Iw": iw, "In": inoise, "dt": dt, "branch": 0,
            "n_trajectories": len(chirality), "n_with_clean_samples": len(subset),
            "clean_sample_fraction": float(np.sum(count[subset]) / (len(chirality) * int(sample_time / dt))),
            "radial_variance": var_u, "radial_seed_se": float("nan"),
            "radial_theory": inoise**2 * iw**3 / (2 * (1 + iw**2) ** 2.5),
            "radial_ratio": var_u / (inoise**2 * iw**3 / (2 * (1 + iw**2) ** 2.5)),
            "angular_variance": var_v, "angular_seed_se": float("nan"),
            "angular_theory": inoise**2 * iw / (2 * math.sqrt(1 + iw**2)),
            "angular_ratio": var_v / (inoise**2 * iw / (2 * math.sqrt(1 + iw**2))),
            "covariance": float("nan"),
            "trajectories_with_contact": int(np.sum(np.any(flags[:, [0, 2]] > 0, axis=1))),
            "trajectories_with_clipping": int(np.sum(flags[:, 1] > 0),),
            "trajectories_with_switch": int(np.sum(flags[:, 3] > 0)),
        })
    write_csv(outdir / "tables" / "noise_covariance.csv", rows)
    return rows


def regime_map(outdir: Path, quick: bool = False):
    strengths = (0.08, 0.15, 0.25, 0.50, 0.94, 2.0, 4.0)
    noises = (0.03, 0.10, 0.20, 0.30)
    n_per_branch = 6 if quick else 16
    burn_time = 20.0 if quick else 50.0
    sample_time = 40.0 if quick else 100.0
    rows = []
    for iw in strengths:
        p = orbit_parameters(iw)
        for inoise in noises:
            chirality = np.repeat(np.array([-1.0, 1.0]), n_per_branch)
            initial = np.vstack([make_initial(p["r_star"], c * np.pi / 2) for c in chirality])
            count, sums, flags = stochastic_statistics(
                initial, chirality, iw, inoise, 0.001,
                int(burn_time / 0.001), int(sample_time / 0.001),
                990000 + int(iw * 1000) + int(inoise * 10000),
            )
            subset = np.where(count > 1)[0]
            if len(subset):
                var_v = float(np.sum(sums[subset, 3]) / np.sum(count[subset]) - (np.sum(sums[subset, 1]) / np.sum(count[subset]))**2)
                theory_v = inoise**2 * iw / (2 * math.sqrt(1 + iw**2))
                ratio = var_v / theory_v
            else:
                ratio = float("nan")
            rows.append({
                "Iw": iw, "In": inoise,
                "formal_clearance": p["clearance_star"],
                "formal_ray_distance": p["d_star"],
                "clean_sample_fraction": float(np.sum(count) / (len(chirality) * int(sample_time / 0.001))),
                "angular_variance_ratio": ratio,
                "contact_trajectory_fraction": float(np.mean(np.any(flags[:, [0, 2]] > 0, axis=1))),
                "clipped_trajectory_fraction": float(np.mean(flags[:, 1] > 0)),
                "switch_trajectory_fraction": float(np.mean(flags[:, 3] > 0)),
            })
    write_csv(outdir / "tables" / "regime_map.csv", rows)
    return rows


def source_audit(outdir: Path):
    rows = []
    test_cases = [
        (0.0, 0.0, 0.2), (0.35, -0.1, 1.2), (-0.55, 0.25, -2.4),
        (0.88, 0.0, 0.0), (0.0, -0.89, -0.5),
    ]
    for x, y, theta in test_cases:
        projection = x * math.cos(theta) + y * math.sin(theta)
        roots = sorted([
            -projection - math.sqrt(projection**2 + 1 - x*x - y*y),
            -projection + math.sqrt(projection**2 + 1 - x*x - y*y),
        ])
        positive = [root for root in roots if root >= 0]
        d = min(positive)
        hx, hy = x + d * math.cos(theta), y + d * math.sin(theta)
        beta = math.atan2(hy, hx)
        incidence = math.atan2(math.sin(theta-beta), math.cos(theta-beta))
        rows.append({
            "x": x, "y": y, "theta": theta, "negative_root": roots[0],
            "selected_positive_root": d, "hit_radius_error": math.hypot(hx, hy)-1,
            "wall_normal_angle": beta, "wrapped_incidence": incidence,
            "torque_sign": int(np.sign(incidence)),
        })
    write_csv(outdir / "tables" / "source_audit_cases.csv", rows)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def figures(outdir: Path, det_rows, det_traces, stability_rows, stability_traces, noise_rows, regime_rows):
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

    # A: exact ray geometry and both safety conditions.
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)
    x, y, theta = 0.48, 0.20, 1.05
    projection = x*np.cos(theta)+y*np.sin(theta)
    d = -projection + np.sqrt(projection**2 + 1-x*x-y*y)
    hit = np.array([x+d*np.cos(theta), y+d*np.sin(theta)])
    ang = np.linspace(0, 2*np.pi, 500)
    axes[0].plot(np.cos(ang), np.sin(ang), color="0.3")
    axes[0].arrow(x, y, d*np.cos(theta), d*np.sin(theta), width=.008, length_includes_head=True, color="C0")
    axes[0].arrow(hit[0], hit[1], .18*hit[0], .18*hit[1], width=.006, color="C3")
    axes[0].scatter([x, hit[0]], [y, hit[1]], c=["C0", "C3"])
    axes[0].annotate("first positive ray hit", hit, xytext=(-65, 16), textcoords="offset points")
    axes[0].set(title="A. Source-exact forward-ray geometry", aspect="equal", xlim=(-1.2,1.3), ylim=(-1.2,1.2))
    iws = np.logspace(-2, 1, 400)
    rr = 1/np.sqrt(1+iws**2)
    dd = iws/np.sqrt(1+iws**2)
    axes[1].semilogx(iws, 1-rr, label=r"radial clearance $R-r_*$")
    axes[1].semilogx(iws, dd, label=r"ray distance $d_*$")
    axes[1].axhline(.1, color="k", ls="--", label="production safety threshold")
    for iw in (.25,.94,2): axes[1].axvline(iw, color="0.75", lw=.8)
    axes[1].set(xlabel=r"$I_w$", ylabel="distance", title="Two distinct safeguards", ylim=(0,1.02))
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(figdir / "S1_source_geometry_and_safeguards.png", dpi=220)
    plt.close(fig)

    # A/B: requested orbit radius and angular-frequency tests.
    iw_curve = np.linspace(0.15, 2.2, 500)
    for quantity, ylabel, filename in (
        ("radius", r"orbital radius $r$", "A_orbit_radius_vs_Iw.png"),
        ("frequency", r"$|\dot\varphi|$", "B_orbital_frequency_vs_Iw.png"),
    ):
        fig, ax = plt.subplots(figsize=(5.7, 4.2), constrained_layout=True)
        theory = (
            1 / np.sqrt(1 + iw_curve**2)
            if quantity == "radius"
            else np.sqrt(1 + iw_curve**2)
        )
        ax.plot(iw_curve, theory, color="k", label="smooth theory")
        for iw in (.25, .94, 2):
            subset = [r for r in det_rows if r["Iw"] == iw and r["dt"] == .001]
            values = np.array([
                r["measured_radius"] if quantity == "radius" else abs(r["measured_frequency"])
                for r in subset
            ])
            ax.errorbar(iw, np.mean(values), yerr=np.std(values, ddof=1), fmt="o",
                        color="C3" if iw == .25 else "C0", capsize=3)
        ax.axvspan(0, math.sqrt(1/.9**2-1), color="C3", alpha=.12,
                   label="formal orbit inside sticky layer")
        ax.set(xlabel=r"wall strength $I_w$", ylabel=ylabel,
               title=("A. Radius: theory versus source-exact update" if quantity == "radius"
                      else "B. Turning rate: theory versus source-exact update"))
        ax.legend(frameon=False, fontsize=8)
        fig.savefig(figdir / filename, dpi=220)
        plt.close(fig)

    # Supplement: deterministic production trajectories and convergence.
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    for iw, (times, state, diag) in det_traces.items():
        r = np.hypot(state[:,0], state[:,1])
        axes[0].plot(times, r, label=fr"$I_w={iw}$")
        axes[0].axhline(orbit_parameters(iw)["r_star"], color=axes[0].lines[-1].get_color(), ls="--", lw=1)
    axes[0].axhline(.9, color="k", ls=":", label="sticky boundary")
    axes[0].set(xlabel="time", ylabel=r"$r$", title="Source-exact deterministic radius")
    axes[0].legend(frameon=False, fontsize=8)
    for iw in (.25,.94,2):
        subset=[r for r in det_rows if r["Iw"]==iw and r["initial_scope"]=="local" and r["initial_chirality"]==1]
        axes[1].plot([r["dt"] for r in subset], [abs(r["radius_error"]) for r in subset], "o-", label=fr"radius, $I_w={iw}$")
        axes[1].plot([r["dt"] for r in subset], [abs(r["frequency_error"]) for r in subset], "s--", color=axes[1].lines[-1].get_color())
    axes[1].set(xlabel=r"$Delta t$", ylabel="absolute error", yscale="log", xscale="log", title="Timestep convergence (solid radius; dashed frequency)")
    axes[1].legend(frameon=False, fontsize=7)
    fig.savefig(figdir / "S2_deterministic_orbits_and_convergence.png", dpi=220)
    plt.close(fig)

    # C: local linear response.
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    for ax, iw in zip(axes, (.25,.94,2)):
        t, obs, pred = stability_traces[iw]
        ax.plot(t, obs[:,0], label="nonlinear/source")
        ax.plot(t, pred[:,0], "--", label="linear theory")
        ax.set(title=fr"$I_w={iw}$", xlabel="time", ylabel=r"$r-r_*$")
        if iw == .25: ax.text(.03,.95,"contact safeguard active",transform=ax.transAxes,va="top",color="C3")
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(figdir / "C_linear_stability.png", dpi=220)
    plt.close(fig)

    # D: weak-noise covariance.
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    for iw in (.25,.94,2):
        xvals=[]; yv=[]; yr=[]; ev=[]; er=[]; tv=[]; tr=[]
        for inoise in (.03,.05,.1):
            subset=[r for r in noise_rows if r["Iw"]==iw and r["In"]==inoise and r["dt"]==.001 and r["branch"] in (-1,1)]
            xvals.append(inoise**2)
            angular_values = np.array([r["angular_variance"] for r in subset], dtype=float)
            radial_values = np.array([r["radial_variance"] for r in subset], dtype=float)
            yv.append(np.nan if np.all(np.isnan(angular_values)) else np.nanmean(angular_values))
            yr.append(np.nan if np.all(np.isnan(radial_values)) else np.nanmean(radial_values))
            ev.append(np.sqrt(np.nansum([r["angular_seed_se"]**2 for r in subset]))/2)
            er.append(np.sqrt(np.nansum([r["radial_seed_se"]**2 for r in subset]))/2)
            tv.append(subset[0]["angular_theory"]); tr.append(subset[0]["radial_theory"])
        if np.all(np.isnan(yv)):
            c="C3"
            axes[0].plot(xvals,tv,"--",color=c,label=fr"$I_w={iw}$ theory only (sticky)")
            axes[1].plot(xvals,tr,"--",color=c)
        else:
            line=axes[0].errorbar(xvals,yv,yerr=ev,fmt="o-",capsize=2,label=fr"$I_w={iw}$")
            c=line[0].get_color(); axes[0].plot(xvals,tv,"--",color=c)
            axes[1].errorbar(xvals,yr,yerr=er,fmt="o-",capsize=2,label=fr"$I_w={iw}$",color=c)
            axes[1].plot(xvals,tr,"--",color=c)
    axes[0].set(xlabel=r"$I_n^2$",ylabel=r"Var$(\alpha-\alpha_*)$",title="D. Angular variance (points; OU dashed)")
    axes[1].set(xlabel=r"$I_n^2$",ylabel=r"Var$(r-r_*)$",title="Radial variance (points; OU dashed)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(figdir / "D_noise_covariance.png", dpi=220)
    plt.close(fig)

    # E: empirical limits of the local theory.
    strengths=sorted(set(float(r["Iw"]) for r in regime_rows)); noises=sorted(set(float(r["In"]) for r in regime_rows))
    contact=np.full((len(noises),len(strengths)),np.nan); switches=contact.copy(); ratio=contact.copy()
    for row in regime_rows:
        i=noises.index(float(row["In"])); j=strengths.index(float(row["Iw"])); contact[i,j]=row["contact_trajectory_fraction"]; switches[i,j]=row["switch_trajectory_fraction"]; ratio[i,j]=row["angular_variance_ratio"]
    fig, axes=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    for ax,z,title,cmap,vmin,vmax in zip(axes,(contact,switches,np.abs(np.log10(ratio))),("contact/sticky fraction","chirality-switch fraction",r"$|\log_{10}(\mathrm{Var}/\mathrm{OU})|$"),("magma","magma","viridis"),(0,0,0),(1,1,1)):
        im=ax.imshow(z,origin="lower",aspect="auto",cmap=cmap,vmin=vmin,vmax=vmax)
        ax.set_xticks(range(len(strengths)),strengths); ax.set_yticks(range(len(noises)),noises)
        ax.set(xlabel=r"$I_w$",ylabel=r"$I_n$",title="E. "+title)
        fig.colorbar(im,ax=ax,shrink=.82)
    fig.savefig(figdir / "E_validity_regime_map.png",dpi=220)
    plt.close(fig)


def summarize(outdir: Path, det_rows, stability_rows, noise_rows, regime_rows):
    smooth_det = [r for r in det_rows if r["Iw"] in (.94,2.0) and r["dt"]==.001 and r["initial_scope"]=="local"]
    invalid_det = [r for r in det_rows if r["Iw"]==.25 and r["dt"]==.001]
    valid_noise = [r for r in noise_rows if r["Iw"] in (.94,2.0) and r["dt"]==.001 and r["branch"] in (-1,1)]
    summary = {
        "source_file": "BEM_circular_parallel_instance.py",
        "production_dt": PRODUCTION_DT,
        "smooth_orbit_minimum_Iw_from_positional_safety": math.sqrt(1/0.9**2-1),
        "smooth_local_max_abs_radius_error": max(abs(float(r["radius_error"])) for r in smooth_det),
        "smooth_local_max_abs_frequency_error": max(abs(float(r["frequency_error"])) for r in smooth_det),
        "Iw_0p25_mean_freeze_fraction": float(np.mean([r["freeze_fraction"] for r in invalid_det])),
        "valid_noise_angular_ratio_median": float(np.nanmedian([r["angular_ratio"] for r in valid_noise])),
        "valid_noise_radial_ratio_median": float(np.nanmedian([r["radial_ratio"] for r in valid_noise])),
        "valid_noise_max_abs_covariance": float(np.nanmax(np.abs([r["covariance"] for r in valid_noise]))),
        "stability_valid_median_normalized_rmse": float(np.median([r["normalized_linear_response_rmse"] for r in stability_rows if r["Iw"] in (.94,2.0)])),
        "n_regime_cells": len(regime_rows),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_csv(outdir / "summary.csv", [{"metric": k, "value": v} for k,v in summary.items()])
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument("--quick",action="store_true",help="short smoke-test ensemble")
    args=parser.parse_args(); outdir=args.outdir
    (outdir/"tables").mkdir(parents=True,exist_ok=True)
    source_audit(outdir)
    det_rows,det_traces=deterministic_validation(outdir)
    stability_rows,stability_traces=stability_validation(outdir)
    noise_rows=stochastic_validation(outdir,args.quick)
    regime_rows=regime_map(outdir,args.quick)
    figures(outdir,det_rows,det_traces,stability_rows,stability_traces,noise_rows,regime_rows)
    summary=summarize(outdir,det_rows,stability_rows,noise_rows,regime_rows)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
