#!/usr/bin/env python3
"""Exact focal-fish replay at rounded-divider tangent headings.

The implementation mirrors ``BEM_bichamber_parallel_instance.py`` but avoids a
full BEM solve for every trial heading.  For one focal fish the boundary right
hand side is exactly affine in cos(theta) and sin(theta); one three-column LU
solve therefore gives the BEM field for any focal heading without fitting.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import lu_solve

HERE = Path(__file__).resolve().parent
REPLAY = HERE.parents[1] / "pre_crossing_turning" / "scripts"
sys.path.insert(0, str(REPLAY))
from replay_dynamics import (  # noqa: E402
    Parameters, analytic_geometry, direct_hydrodynamics, hits_circle,
    intersects, visible_neighbors,
)


def wrap(x):
    return np.arctan2(np.sin(x), np.cos(x))


def chamber_side(origin: str) -> float:
    return 1.0 if origin == "left" else -1.0


def local_to_global(theta_local, side):
    return np.arctan2(np.sin(theta_local), side * np.cos(theta_local))


def global_to_local(theta_global, side):
    return np.arctan2(np.sin(theta_global), side * np.cos(theta_global))


def tangent_geometry(position, velocity, origin, p: Parameters):
    """Tangent corridor and its material boundary velocities.

    Local x is distance toward the divider.  A valid corridor requires the
    focal point to be outside both rounded tips and alpha_minus < alpha_plus.
    """
    side = chamber_side(origin)
    x, y = map(float, position)
    vx, vy = map(float, velocity)
    d = -side * x
    dd = -side * vx
    out = {"side": side, "d": d, "d_dot": dd}
    vals = {}
    for label, cy, sign in (("minus", -p.gap / 2 - p.tip_r, +1.0),
                            ("plus", p.gap / 2 + p.tip_r, -1.0)):
        v = cy - y
        dv = -vy
        rho = np.hypot(d, v)
        if d <= 0 or rho <= p.tip_r:
            vals[label] = (np.nan, np.nan, rho)
            continue
        beta = np.arcsin(p.tip_r / rho)
        alpha = np.arctan2(v, d) + sign * beta
        atan_dot = (d * dv - v * dd) / rho**2
        rho_dot = (d * dd + v * dv) / rho
        beta_dot = -p.tip_r * rho_dot / (rho**2 * np.sqrt(1-(p.tip_r/rho)**2))
        alpha_dot = atan_dot + sign * beta_dot
        vals[label] = (alpha, alpha_dot, rho)
    out.update(alpha_minus=vals["minus"][0], alpha_minus_dot=vals["minus"][1],
               rho_minus=vals["minus"][2], alpha_plus=vals["plus"][0],
               alpha_plus_dot=vals["plus"][1], rho_plus=vals["plus"][2])
    out["width"] = out["alpha_plus"] - out["alpha_minus"]
    out["valid"] = bool(np.isfinite(out["width"]) and out["width"] > 0)
    return out


BRANCHES = ("outer_circle", "left_vertical", "right_vertical", "upper_tip",
            "lower_tip", "upper_right_smooth", "upper_left_smooth",
            "lower_left_smooth", "lower_right_smooth")


def _circle_candidate(loc, heading, center, radius, th_min, th_max, inner):
    x0, y0 = loc - center
    ct, st = np.cos(heading), np.sin(heading)
    disc = (x0*ct+y0*st)**2-(x0*x0+y0*y0-radius*radius)
    if disc < 0:
        return None
    root = np.sqrt(disc)
    t = -(x0*ct+y0*st) + (-root if inner else root)
    if t < 0:
        return None
    normal = np.arctan2(y0+t*st, x0+t*ct)
    if normal > th_max or normal < th_min:
        return None
    if inner:
        normal += np.pi
    return float(t), float(normal)


def ray_wall_branch(fish, geom, p: Parameters):
    """Exact GetInterc candidates, augmented only with a branch label."""
    loc = np.asarray(fish[:2], float)
    heading = float(fish[2])
    x0, y0 = loc
    ct, st = np.cos(heading), np.sin(heading)
    candidates = []
    disc = (x0*ct+y0*st)**2-(x0*x0+y0*y0-p.R*p.R)
    if disc >= 0:
        t = -(x0*ct+y0*st)+np.sqrt(disc)
        normal = np.arctan2(y0+t*st, x0+t*ct)
        if t >= 0 and not (abs(normal) < geom["th1"] and abs(normal) > geom["th8"]):
            candidates.append((float(t), float(normal), BRANCHES[0]))
    for branch, xwall, normal in ((BRANCHES[1], -p.tip_r, 0.0),
                                  (BRANCHES[2], p.tip_r, np.pi)):
        if ct != 0:
            t = (xwall-x0)/ct
            yy = y0+t*st
            if t >= 0 and p.gap/2+p.tip_r <= abs(yy) <= geom["node12"][1]:
                candidates.append((float(t), normal, branch))
    specs = ((BRANCHES[3], geom["c6"], p.tip_r, -np.pi, 0, True),
             (BRANCHES[4], geom["c5"], p.tip_r, 0, np.pi, True),
             (BRANCHES[5], geom["c1"], p.smooth_r, geom["th8"], np.pi, False),
             (BRANCHES[6], geom["c2"], p.smooth_r, 0, geom["th1"], False),
             (BRANCHES[7], geom["c3"], p.smooth_r, geom["th2"], 0, False),
             (BRANCHES[8], geom["c4"], p.smooth_r, -np.pi, geom["th7"], False))
    for branch, *args in specs:
        ans = _circle_candidate(loc, heading, *args)
        if ans is not None:
            candidates.append((*ans, branch))
    if not candidates:
        return dict(valid=False, branch="none", distance=np.nan, normal=np.nan)
    distance, normal, branch = min(candidates, key=lambda z: z[0])
    return dict(valid=True, branch=branch, distance=distance, normal=normal)


def exact_wall(fish, geom, p: Parameters):
    hit = ray_wall_branch(fish, geom, p)
    if not hit["valid"]:
        return np.nan, hit
    phi = wrap(fish[2]-hit["normal"])
    omega = p.Iw*np.sign(phi)/hit["distance"]
    amplified = abs(hit["distance"]) <= abs(p.safedis)
    if amplified:
        omega *= 1+1/hit["distance"]**(p.high_order-1)
    omega = float(np.clip(omega, -np.pi/p.dt, np.pi/p.dt))
    hit.update(phi=float(phi), amplified=amplified, wall_clipped=abs(omega) >= np.pi/p.dt)
    return omega, hit


def social_focal(state, focal, theta, visible, p: Parameters):
    js = np.flatnonzero(visible[focal])
    if not len(js):
        return 0.0, 0.0
    rho = state[:2, js]-state[:2, focal, None]
    dist = np.linalg.norm(rho, axis=0)
    bearing = np.arctan2(rho[1], rho[0])
    rel = wrap(bearing-theta)
    visual = 1+np.cos(rel)
    den = visual.sum()
    if den == 0:
        return 0.0, 0.0
    align = np.sum(p.Ip*np.sin(state[2, js]-theta)*visual)/den
    attract = np.sum(dist*np.sin(rel)*visual)/den
    return float(align), float(attract)


def direct_focal(state, focal, theta, p: Parameters):
    z = state[0]+1j*state[1]
    mask = np.arange(state.shape[1]) != focal
    dz = z[focal]-z[mask]
    aj = state[2, mask]
    u = np.sum(np.exp(1j*aj)/(dz**2+p.delta**2))*p.If/np.pi
    omega = np.sum(np.imag(np.exp(1j*(2*theta+aj))*dz /
                             (dz**2+p.delta**2)**2))*2*p.If/np.pi
    return np.array([u.real, -u.imag]), float(omega)


@dataclass
class FocalBEM:
    state: np.ndarray
    focal: int
    boundary: dict
    sigma_other: np.ndarray
    sigma_cos: np.ndarray
    sigma_sin: np.ndarray
    p: Parameters

    @classmethod
    def build(cls, state, focal, boundary, lu, piv, p):
        x, y = boundary["x"], boundary["y"]
        nx, ny = boundary["nx"], boundary["ny"]
        z = state[0]+1j*state[1]
        dz = x[:, None]+1j*y[:, None]-z[None, :]
        inc = p.If*np.exp(1j*state[2])[None, :]/(dz**2+p.delta**2)/np.pi
        f_all = np.sum(nx[:, None]*inc.real-ny[:, None]*inc.imag, axis=1)
        base = p.If/(dz[:, focal]**2+p.delta**2)/np.pi
        f_cos = nx*base.real-ny*base.imag
        ibase = 1j*base
        f_sin = nx*ibase.real-ny*ibase.imag
        th0 = state[2, focal]
        f_other = f_all-np.cos(th0)*f_cos-np.sin(th0)*f_sin
        rhs = -np.column_stack((f_other, f_cos, f_sin))
        sig = lu_solve((lu, piv), rhs, check_finite=False)
        return cls(np.asarray(state), int(focal), boundary, sig[:, 0], sig[:, 1], sig[:, 2], p)

    def evaluate(self, theta):
        sigma = self.sigma_other+np.cos(theta)*self.sigma_cos+np.sin(theta)*self.sigma_sin
        x, y, weights = self.boundary["x"], self.boundary["y"], self.boundary["weights"]
        dx = x-self.state[0, self.focal]
        dy = y-self.state[1, self.focal]
        dz = dx+1j*dy
        factor = sigma*weights/(2*np.pi)
        field = factor/dz
        velocity = np.array([field.real.sum(), -field.imag.sum()])
        r4 = (dx*dx+dy*dy)**2
        dudx = np.sum(factor*(dx*dx-dy*dy)/r4)
        dudy = np.sum(factor*(2*dx*dy)/r4)
        dvdx = dudy
        dvdy = np.sum(factor*(dy*dy-dx*dx)/r4)
        omega = -(np.sin(theta)*np.cos(theta)*(dudx-dvdy)
                  + np.sin(theta)**2*dvdx-np.cos(theta)**2*dudy)
        return velocity, float(omega)


def focal_replay(state, focal, theta_global, visible, bem: FocalBEM, geom, p):
    fish = np.r_[state[:2, focal], theta_global]
    wall, hit = exact_wall(fish, geom, p)
    align, attract = social_focal(state, focal, theta_global, visible, p)
    u_direct, direct = direct_focal(state, focal, theta_global, p)
    u_bem, boundary = bem.evaluate(theta_global)
    terms = dict(wall=wall, alignment=align, attraction=attract,
                 direct_hydro=direct, boundary_hydro=boundary)
    total_raw = sum(terms.values())
    total = float(np.clip(total_raw, -np.pi/p.dt, np.pi/p.dt))
    rdot = np.array([np.cos(theta_global), np.sin(theta_global)])+u_direct+u_bem
    nxt = state[:2, focal]+rdot*p.dt
    jumping = any(intersects(state[:2, focal], nxt, c, d) for c, d in geom["segments"])
    override = bool(jumping or np.any(np.abs(rdot) > 500))
    if override:
        terms["direct_hydro_effective"] = 0.0
        terms["boundary_hydro_effective"] = 0.0
        total = float(np.clip(wall+align+attract, -np.pi/p.dt, np.pi/p.dt))
        rdot = np.array([np.cos(theta_global), np.sin(theta_global)])
    else:
        terms["direct_hydro_effective"] = direct
        terms["boundary_hydro_effective"] = boundary
    return dict(terms=terms, total=total, total_raw=float(total_raw), rdot=rdot,
                override=override, hit=hit)


def classify(mu_out, mu_in):
    if not np.isfinite(mu_out+mu_in):
        return "invalid"
    if mu_out > 0:
        return "deflection" if mu_in >= 0 else "repelling"
    if mu_in > 0:
        return "pinning"
    return "transmission"


__all__ = ["Parameters", "analytic_geometry", "visible_neighbors", "wrap",
           "chamber_side", "local_to_global", "global_to_local",
           "tangent_geometry", "ray_wall_branch", "exact_wall", "FocalBEM",
           "focal_replay", "classify"]
