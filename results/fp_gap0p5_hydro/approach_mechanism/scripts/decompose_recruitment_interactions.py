#!/usr/bin/env python3
"""Decompose social/direct/BEM terms during held-out successful recruitment episodes."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.linalg import lu_solve
from scipy.spatial import Delaunay


R = 7.158910531638177
GAP = R
TIP_R = 0.5
IP = 4.0
IF = 1e-2
DELTA = 1e-2


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def intersects(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> bool:
    def ccw(p, q, r):
        return (r[1] - p[1]) * (q[0] - p[0]) > (q[1] - p[1]) * (r[0] - p[0])
    return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)


def hits_circle(a: np.ndarray, b: np.ndarray, center: np.ndarray, radius: float) -> bool:
    direction = b - a
    f = a - center
    aa = np.dot(direction, direction)
    if aa == 0:
        return False
    bb = 2 * np.dot(f, direction)
    cc = np.dot(f, f) - radius**2
    disc = bb**2 - 4 * aa * cc
    if disc < 0:
        return False
    root = np.sqrt(disc)
    return any(0 <= t <= 1 for t in ((-bb - root) / (2 * aa), (-bb + root) / (2 * aa)))


def visible_neighbors(location: np.ndarray) -> np.ndarray:
    tri = Delaunay(location.T)
    vn = np.zeros((location.shape[1], location.shape[1]), dtype=bool)
    for simplex in tri.simplices:
        for i, j in ((simplex[0], simplex[1]), (simplex[1], simplex[2]), (simplex[2], simplex[0])):
            vn[i, j] = vn[j, i] = True
    outer_y = np.sqrt(R**2 - TIP_R**2)
    cap = GAP / 2 + TIP_R
    segments = ((np.array([-TIP_R, cap]), np.array([-TIP_R, outer_y])),
                (np.array([TIP_R, outer_y]), np.array([TIP_R, cap])),
                (np.array([-TIP_R, -outer_y]), np.array([-TIP_R, -cap])),
                (np.array([TIP_R, -cap]), np.array([TIP_R, -outer_y])))
    circles = (np.array([0.0, cap]), np.array([0.0, -cap]))
    for i in range(len(vn)):
        for j in np.flatnonzero(vn[i]):
            if j <= i:
                continue
            a, b = location[:, i], location[:, j]
            blocked = any(intersects(a, b, c, d) for c, d in segments)
            blocked |= any(hits_circle(a, b, center, TIP_R) for center in circles)
            if blocked:
                vn[i, j] = vn[j, i] = False
    return vn


def social_parts(state: np.ndarray, cross_source: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    location = state[:2]
    angle = state[2]
    vn = visible_neighbors(location)
    same = np.zeros(len(angle)); cross = np.zeros(len(angle))
    for i in range(len(angle)):
        js = np.flatnonzero(vn[i])
        if not len(js):
            continue
        rho = location[:, js] - location[:, i, None]
        distance = np.linalg.norm(rho, axis=0)
        ex = np.array([np.cos(angle[i]), np.sin(angle[i])])
        dot = np.clip((ex[:, None] * rho).sum(axis=0) / distance, -1, 1)
        sign = np.sign(ex[0] * rho[1] - ex[1] * rho[0]); sign[sign == 0] = 1
        rel = sign * np.arccos(dot)
        visual = 1 + np.cos(rel)
        term = (IP * np.sin(angle[js] - angle[i]) + distance * np.sin(rel)) * visual
        denom = visual.sum()
        cross[i] = term[cross_source[js]].sum() / denom
        same[i] = term[~cross_source[js]].sum() / denom
    return same, cross


def direct_parts(state: np.ndarray, cross_source: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    z = state[0] + 1j * state[1]
    angle = state[2]
    n = len(angle)
    us = np.zeros((2, n)); uc = np.zeros((2, n)); ws = np.zeros(n); wc = np.zeros(n)
    for i in range(n):
        js = np.arange(n) != i
        dz = z[i] - z[js]
        source_angle = angle[js]
        uterm = np.exp(1j * source_angle) / (dz**2 + DELTA**2) * IF / np.pi
        wterm = (np.imag(np.exp(1j * (2 * angle[i] + source_angle)) * dz /
                         (dz**2 + DELTA**2) ** 2) * 2 * IF / np.pi)
        is_cross = cross_source[js]
        for target, mask in ((uc, is_cross), (us, ~is_cross)):
            val = uterm[mask].sum()
            target[:, i] = (val.real, -val.imag)
        wc[i] = wterm[is_cross].sum(); ws[i] = wterm[~is_cross].sum()
    return us, uc, ws, wc


def bem_part(state: np.ndarray, source_mask: np.ndarray, geometry: dict[str, np.ndarray],
             lu: np.ndarray, piv: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xcoord, ycoord = geometry["x"], geometry["y"]
    nx, ny, weights = geometry["nx"], geometry["ny"], geometry["weights"]
    z = state[0] + 1j * state[1]
    angle = state[2]
    dz = xcoord[:, None] + 1j * ycoord[:, None] - z[None, source_mask]
    incident = IF * np.exp(1j * angle[source_mask])[None, :] / (dz**2 + DELTA**2) / np.pi
    f = np.sum(nx[:, None] * incident.real - ny[:, None] * incident.imag, axis=1)
    sigma = lu_solve((lu, piv), -f)
    dx = xcoord[None, :] - state[0, :, None]
    dy = ycoord[None, :] - state[1, :, None]
    dz_eval = dx + 1j * dy
    field = sigma[None, :] * weights[None, :] / (2 * np.pi * dz_eval)
    velocity = np.vstack((field.real.sum(axis=1), -field.imag.sum(axis=1)))
    factor = sigma[None, :] * weights[None, :] / (2 * np.pi)
    radius4 = (dx**2 + dy**2) ** 2
    dudx = np.sum(factor * (dx**2 - dy**2) / radius4, axis=1)
    dudy = np.sum(factor * (2 * dx * dy) / radius4, axis=1)
    dvdx = dudy.copy()
    dvdy = np.sum(factor * (dy**2 - dx**2) / radius4, axis=1)
    omega = -(np.sin(angle) * np.cos(angle) * (dudx - dvdy) +
              np.sin(angle)**2 * dvdx - np.cos(angle)**2 * dudy)
    return velocity, omega


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--geometry-file", type=Path, required=True)
    parser.add_argument("--lu-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.analysis_dir.resolve()
    episodes = list(csv.DictReader((root / "tables" / "nucleation_episodes.csv").open()))
    g = sio.loadmat(args.geometry_file)
    normal = np.asarray(g["normal_vec"])
    geometry = {"x": np.asarray(g["Sx"]).reshape(-1), "y": np.asarray(g["Sy"]).reshape(-1),
                "weights": np.asarray(g["weights"]).reshape(-1), "nx": normal[:, 0], "ny": normal[:, 1]}
    lu = np.load(args.lu_dir / "lu_R_50_wid_50.npy", mmap_mode="r")
    piv = np.load(args.lu_dir / "piv_R_50_wid_50.npy")
    rows = []
    for seed in range(15, 20):
        candidates = [r for r in episodes if int(r["seed"]) == seed and r["success"] == "1"
                      and float(r["full_dwell_before"]) >= 20]
        if not candidates:
            continue
        target = np.median([float(r["d_at_cross"]) for r in candidates])
        event = min(candidates, key=lambda r: abs(float(r["d_at_cross"]) - target))
        origin = -1 if event["origin"] == "left" else 1
        raw = args.raw_dir / f"If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat"
        state_all = np.asarray(sio.loadmat(raw, variable_names=["state"], squeeze_me=True)["state"])
        first = int(round(float(event["first_cross_time"]) / .001))
        stages = (1, 5, 20, 50)
        for target_count in stages:
            idx = first
            while idx < state_all.shape[2]:
                cross_source = state_all[0, :, idx] >= 0 if origin == -1 else state_all[0, :, idx] < 0
                if cross_source.sum() >= target_count:
                    break
                idx += 10
            state = state_all[:, :, idx]
            cross_source = state[0] >= 0 if origin == -1 else state[0] < 0
            focal = ~cross_source
            inward = 1.0 if origin == -1 else -1.0
            vs, vc = social_parts(state, cross_source)
            uds, udc, wds, wdc = direct_parts(state, cross_source)
            ubs, wbs = bem_part(state, ~cross_source, geometry, lu, piv)
            ubc, wbc = bem_part(state, cross_source, geometry, lu, piv)
            visibility = visible_neighbors(state[:2])
            # Contribution to growth of the inward heading component
            # c = n_in . e(theta): dc/dt = -n_in,x sin(theta) omega.
            align_factor = -inward * np.sin(state[2, focal])
            rows.append({
                "seed": seed, "origin": event["origin"], "stage_opposite_target": target_count,
                "opposite_count": int(cross_source.sum()), "time_since_first": idx * .001 - float(event["first_cross_time"]),
                "social_same_turn_mean": float(np.mean(vs[focal])), "social_cross_turn_mean": float(np.mean(vc[focal])),
                "social_same_inward_alignment_mean": float(np.mean(align_factor * vs[focal])),
                "social_cross_inward_alignment_mean": float(np.mean(align_factor * vc[focal])),
                "direct_same_turn_mean": float(np.mean(wds[focal])), "direct_cross_turn_mean": float(np.mean(wdc[focal])),
                "direct_same_inward_alignment_mean": float(np.mean(align_factor * wds[focal])),
                "direct_cross_inward_alignment_mean": float(np.mean(align_factor * wdc[focal])),
                "bem_same_turn_mean": float(np.mean(wbs[focal])), "bem_cross_turn_mean": float(np.mean(wbc[focal])),
                "bem_same_inward_alignment_mean": float(np.mean(align_factor * wbs[focal])),
                "bem_cross_inward_alignment_mean": float(np.mean(align_factor * wbc[focal])),
                "direct_same_inward_mean": float(np.mean(inward * uds[0, focal])),
                "direct_cross_inward_mean": float(np.mean(inward * udc[0, focal])),
                "bem_same_inward_mean": float(np.mean(inward * ubs[0, focal])),
                "bem_cross_inward_mean": float(np.mean(inward * ubc[0, focal])),
                "net_cross_hydro_inward_mean": float(np.mean(inward * (udc[0, focal] + ubc[0, focal]))),
                "net_same_hydro_inward_mean": float(np.mean(inward * (uds[0, focal] + ubs[0, focal]))),
                "cross_visible_edges": int(visibility[focal][:, cross_source].sum()),
            })
        del state_all
    write_csv(root / "tables" / "interaction_decomposition.csv", rows)
    summary_rows = []
    rng = np.random.default_rng(3005)
    for stage in (1, 5, 20, 50):
        sub = [r for r in rows if r["stage_opposite_target"] == stage]
        for key in ("social_cross_inward_alignment_mean", "direct_cross_inward_alignment_mean",
                    "bem_cross_inward_alignment_mean", "direct_cross_inward_mean", "bem_cross_inward_mean",
                    "net_cross_hydro_inward_mean", "net_same_hydro_inward_mean", "cross_visible_edges"):
            values = np.array([r[key] for r in sub], float)
            boot = np.mean(values[rng.integers(0, len(values), (5000, len(values)))], axis=1)
            summary_rows.append({"stage_opposite_target": stage, "metric": key,
                                 "mean": float(np.mean(values)), "median": float(np.median(values)),
                                 "seed_bootstrap_mean_ci_low": float(np.percentile(boot, 2.5)),
                                 "seed_bootstrap_mean_ci_high": float(np.percentile(boot, 97.5)),
                                 "independent_seeds": len(values)})
    write_csv(root / "tables" / "interaction_stage_summary.csv", summary_rows)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    metrics = (("social_cross_inward_alignment_mean", "social increase of inward heading"),
               ("direct_cross_inward_mean", "direct hydro inward velocity"),
               ("bem_cross_inward_mean", "BEM-mediated inward velocity"))
    colors = ("#0072B2", "#D55E00", "#009E73")
    for key, label, color in zip((m[0] for m in metrics), (m[1] for m in metrics), colors):
        for seed in sorted({int(r["seed"]) for r in rows}):
            sub = [r for r in rows if int(r["seed"]) == seed]
            axes[0].plot([r["opposite_count"] for r in sub], [r[key] for r in sub], color=color, alpha=.22)
        for stage in (1, 5, 20, 50):
            vals = [r[key] for r in rows if r["stage_opposite_target"] == stage]
            axes[0].scatter(stage, np.median(vals), color=color, s=35, label=label if stage == 1 else None)
    axes[0].axhline(0, color="k", lw=.8)
    axes[0].set_xlabel("fish already in opposite chamber")
    axes[0].set_ylabel("mean contribution on fish remaining in origin chamber")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=.2)
    for seed in sorted({int(r["seed"]) for r in rows}):
        sub = [r for r in rows if int(r["seed"]) == seed]
        axes[1].plot([r["opposite_count"] for r in sub], [r["cross_visible_edges"] for r in sub], marker="o", alpha=.55)
    axes[1].set_xlabel("fish already in opposite chamber")
    axes[1].set_ylabel("cross-side visible-neighbor edges")
    axes[1].grid(alpha=.2)
    for ext in ("png", "pdf"):
        fig.savefig(root / "figures" / f"06_interaction_decomposition.{ext}", dpi=220)
    plt.close(fig)
    print(f"wrote {len(rows)} stage decompositions")


if __name__ == "__main__":
    main()
