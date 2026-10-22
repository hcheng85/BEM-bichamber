#!/usr/bin/env python3
"""Exact deterministic angular terms used by the archived gap-50 BEM simulator.

This module is intentionally self-contained.  Equations and geometry mirror
``BEM_bichamber_parallel_instance.py``; it does not import that script because
the archived top-level parameter sweep has since been edited to other gaps.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import lu_solve
from scipy.spatial import Delaunay


@dataclass(frozen=True)
class Parameters:
    R: float = 7.158910531638177
    gap: float = 7.158910531638177
    tip_r: float = 0.5
    smooth_r: float = 0.5
    Ip: float = 4.0
    In: float = 0.5
    Iw: float = 0.94
    If: float = 1e-2
    delta: float = 1e-2
    dt: float = 1e-3
    safedis: float = -0.1
    high_order: int = 6


def analytic_geometry(p: Parameters) -> dict[str, np.ndarray | float]:
    R, r, delta, gap = p.R, p.tip_r, p.smooth_r, p.gap
    aux1 = np.array([r, np.sqrt(R**2-r**2)])
    aux2 = np.array([-r, np.sqrt(R**2-r**2)])
    aux3 = np.array([-r, -np.sqrt(R**2-r**2)])
    aux4 = np.array([r, -np.sqrt(R**2-r**2)])

    theta = np.arctan2(aux2[1], aux2[0])
    cy = np.sqrt((R-delta)**2-(R*np.cos(theta)-delta)**2)
    c2 = np.array([R*np.cos(theta)-delta, cy])
    node12 = np.array([R*np.cos(theta), cy])
    node1 = c2/np.linalg.norm(c2)*R
    th1 = np.arctan2(node1[1]-c2[1], node1[0]-c2[0])

    theta = np.arctan2(aux3[1], aux3[0])
    cy = -np.sqrt((R-delta)**2-(R*np.cos(theta)-delta)**2)
    c3 = np.array([R*np.cos(theta)-delta, cy])
    node3 = np.array([R*np.cos(theta), cy])
    node2 = c3/np.linalg.norm(c3)*R
    th2 = np.arctan2(node2[1]-c3[1], node2[0]-c3[0])

    theta = np.arctan2(aux4[1], aux4[0])
    cy = -np.sqrt((R-delta)**2-(R*np.cos(theta)+delta)**2)
    c4 = np.array([R*np.cos(theta)+delta, cy])
    node6 = np.array([R*np.cos(theta), cy])
    node7 = c4/np.linalg.norm(c4)*R
    th7 = np.arctan2(node7[1]-c4[1], node7[0]-c4[0])

    theta = np.arctan2(aux1[1], aux1[0])
    cy = np.sqrt((R-delta)**2-(R*np.cos(theta)+delta)**2)
    c1 = np.array([R*np.cos(theta)+delta, cy])
    node9 = np.array([R*np.cos(theta), cy])
    node8 = c1/np.linalg.norm(c1)*R
    th8 = np.arctan2(node8[1]-c1[1], node8[0]-c1[0])

    node11 = np.array([-r, gap/2+r]); node10 = np.array([r, gap/2+r])
    node4 = np.array([-r, -gap/2-r]); node5 = np.array([r, -gap/2-r])
    c6 = np.array([0.0, gap/2+r]); c5 = np.array([0.0, -gap/2-r])
    segments = np.array([[node11, node12], [node9, node10],
                         [node3, node4], [node5, node6]])
    return dict(c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6,
                node12=node12, th1=th1, th2=th2, th7=th7, th8=th8,
                segments=segments)


def _ccw(a, b, c):
    return (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])


def intersects(a, b, c, d):
    return _ccw(a, c, d) != _ccw(b, c, d) and _ccw(a, b, c) != _ccw(a, b, d)


def hits_circle(a, b, center, radius):
    d = b-a; f = a-center; aa = np.dot(d, d)
    if aa == 0:
        return False
    bb = 2*np.dot(f, d); cc = np.dot(f, f)-radius**2
    disc = bb**2-4*aa*cc
    if disc < 0:
        return False
    root = np.sqrt(disc)
    return any(0 <= t <= 1 for t in ((-bb-root)/(2*aa), (-bb+root)/(2*aa)))


def visible_neighbors(location, geom, p):
    tri = Delaunay(location.T)
    vn = np.zeros((location.shape[1], location.shape[1]), dtype=bool)
    for simplex in tri.simplices:
        for i, j in ((simplex[0], simplex[1]), (simplex[1], simplex[2]),
                     (simplex[2], simplex[0])):
            vn[i, j] = vn[j, i] = True
    for i in range(len(vn)):
        for j in np.flatnonzero(vn[i]):
            if j <= i:
                continue
            a, b = location[:, i], location[:, j]
            blocked = any(intersects(a, b, c, d) for c, d in geom["segments"])
            blocked |= hits_circle(a, b, geom["c6"], p.tip_r)
            blocked |= hits_circle(a, b, geom["c5"], p.tip_r)
            if blocked:
                vn[i, j] = vn[j, i] = False
    return vn


def social_terms(state, geom, p):
    """Return alignment and attraction pieces of exact GetVision."""
    loc, angle = state[:2], state[2]
    vn = visible_neighbors(loc, geom, p)
    align = np.zeros(len(angle)); attract = np.zeros(len(angle))
    for i in range(len(angle)):
        js = np.flatnonzero(vn[i])
        if not len(js):
            continue
        rho = loc[:, js]-loc[:, i, None]
        dist = np.linalg.norm(rho, axis=0)
        ex = np.array([np.cos(angle[i]), np.sin(angle[i])])
        dot = np.clip(np.sum(ex[:, None]*rho, axis=0)/dist, -1, 1)
        sign = np.sign(ex[0]*rho[1]-ex[1]*rho[0]); sign[sign == 0] = 1
        theta = sign*np.arccos(dot)
        visual = 1+np.cos(theta); denom = visual.sum()
        align[i] = np.sum(p.Ip*np.sin(angle[js]-angle[i])*visual)/denom
        attract[i] = np.sum(dist*np.sin(theta)*visual)/denom
    return align, attract, vn


def direct_hydrodynamics(state, p):
    z = state[0]+1j*state[1]; angle = state[2]; n = len(angle)
    velocity = np.zeros((2, n)); omega = np.zeros(n)
    for i in range(n):
        js = np.arange(n) != i; dz = z[i]-z[js]; aj = angle[js]
        u = np.sum(np.exp(1j*aj)/(dz**2+p.delta**2))*p.If/np.pi
        velocity[:, i] = (u.real, -u.imag)
        omega[i] = np.sum(np.imag(np.exp(1j*(2*angle[i]+aj))*dz /
                                    (dz**2+p.delta**2)**2))*2*p.If/np.pi
    return velocity, omega


def boundary_hydrodynamics(state, boundary, lu, piv, p):
    x, y = boundary["x"], boundary["y"]
    nx, ny, weights = boundary["nx"], boundary["ny"], boundary["weights"]
    z = state[0]+1j*state[1]; angle = state[2]
    dz = x[:, None]+1j*y[:, None]-z[None, :]
    incident = p.If*np.exp(1j*angle)[None, :]/(dz**2+p.delta**2)/np.pi
    f = np.sum(nx[:, None]*incident.real-ny[:, None]*incident.imag, axis=1)
    sigma = lu_solve((lu, piv), -f, check_finite=False)
    dx = x[None, :]-state[0, :, None]; dy = y[None, :]-state[1, :, None]
    dz_eval = dx+1j*dy
    field = sigma[None, :]*weights[None, :]/(2*np.pi*dz_eval)
    velocity = np.vstack((field.real.sum(axis=1), -field.imag.sum(axis=1)))
    factor = sigma[None, :]*weights[None, :]/(2*np.pi)
    r4 = (dx**2+dy**2)**2
    dudx = np.sum(factor*(dx**2-dy**2)/r4, axis=1)
    dudy = np.sum(factor*(2*dx*dy)/r4, axis=1)
    dvdx = dudy.copy()
    dvdy = np.sum(factor*(dy**2-dx**2)/r4, axis=1)
    omega = -(np.sin(angle)*np.cos(angle)*(dudx-dvdy)
              + np.sin(angle)**2*dvdx-np.cos(angle)**2*dudy)
    return velocity, omega


def _circle_intersection(loc, orientation, center, radius, th_min, th_max, inner):
    x0, y0 = loc-center; ct, st = np.cos(orientation), np.sin(orientation)
    disc = (x0*ct+y0*st)**2-(x0**2+y0**2-radius**2)
    if disc < 0:
        return 1000.0, 0.0, False
    root = np.sqrt(disc)
    t = -(x0*ct+y0*st)+(-root if inner else root)
    if t < 0:
        return 1000.0, 0.0, False
    th = np.arctan2(y0+t*st, x0+t*ct)
    if th > th_max or th < th_min:
        return 1000.0, 0.0, False
    if inner:
        th += np.pi
    return t, th, True


def ray_wall_intersection(fish, geom, p):
    loc, th0 = fish[:2], fish[2]; x0, y0 = loc
    ct, st = np.cos(th0), np.sin(th0)
    ds, ths = [], []
    disc = (x0*ct+y0*st)**2-(x0**2+y0**2-p.R**2)
    if disc >= 0:
        t = -(x0*ct+y0*st)+np.sqrt(disc)
        th = np.arctan2(y0+t*st, x0+t*ct)
        if t >= 0 and not (abs(th) < geom["th1"] and abs(th) > geom["th8"]):
            ds.append(t); ths.append(th)
    for xwall, normal, ymin in ((-p.tip_r, 0.0, p.gap/2+p.tip_r),
                                (p.tip_r, np.pi, p.gap/2+p.tip_r)):
        if ct != 0:
            t = (xwall-x0)/ct; yy = y0+t*st
            if t >= 0 and ymin <= abs(yy) <= geom["node12"][1]:
                ds.append(t); ths.append(normal)
    specs = ((geom["c6"], p.tip_r, -np.pi, 0, True),
             (geom["c5"], p.tip_r, 0, np.pi, True),
             (geom["c1"], p.smooth_r, geom["th8"], np.pi, False),
             (geom["c2"], p.smooth_r, 0, geom["th1"], False),
             (geom["c3"], p.smooth_r, geom["th2"], 0, False),
             (geom["c4"], p.smooth_r, -np.pi, geom["th7"], False))
    for args in specs:
        d, th, ok = _circle_intersection(loc, th0, *args)
        if ok:
            ds.append(d); ths.append(th)
    if not ds:
        raise RuntimeError(f"no forward wall intersection for {fish}")
    k = int(np.argmin(ds))
    return ths[k], ds[k]


def wall_avoidance(state, geom, p):
    omega = np.zeros(state.shape[1]); distance = np.zeros_like(omega)
    for i in range(state.shape[1]):
        normal, d = ray_wall_intersection(state[:, i], geom, p)
        phi = np.arctan2(np.sin(state[2, i]-normal), np.cos(state[2, i]-normal))
        w = p.Iw*np.sign(phi)/d
        if abs(d) <= abs(p.safedis):
            w *= 1+1/d**(p.high_order-1)
        omega[i] = np.clip(w, -np.pi/p.dt, np.pi/p.dt); distance[i] = d
    return omega, distance


def jumping_mask(state, rdot, geom, p):
    nxt = state[:2]+rdot*p.dt
    mask = np.zeros(state.shape[1], bool)
    for i in range(state.shape[1]):
        mask[i] = any(intersects(state[:2, i], nxt[:, i], c, d)
                      for c, d in geom["segments"])
    return mask


def replay_snapshot(state, state_next, boundary, lu, piv, p, geom):
    align, attract, vn = social_terms(state, geom, p)
    u_direct, w_direct = direct_hydrodynamics(state, p)
    u_bem, w_bem = boundary_hydrodynamics(state, boundary, lu, piv, p)
    w_wall, wall_distance = wall_avoidance(state, geom, p)
    terms = dict(wall=w_wall, alignment=align, attraction=attract,
                 direct_hydro=w_direct, boundary_hydro=w_bem)
    total_raw = sum(terms.values())
    total = np.clip(total_raw, -np.pi/p.dt, np.pi/p.dt)
    e = np.vstack((np.cos(state[2]), np.sin(state[2])))
    rdot = e+u_direct+u_bem
    overridden = jumping_mask(state, rdot, geom, p) | np.any(np.abs(rdot) > 500, axis=0)
    total_nohydro = np.clip(w_wall+align+attract, -np.pi/p.dt, np.pi/p.dt)
    total[overridden] = total_nohydro[overridden]
    # Attribute simulator overrides explicitly: hydrodynamic angular terms are removed.
    terms["direct_hydro_effective"] = w_direct.copy()
    terms["boundary_hydro_effective"] = w_bem.copy()
    terms["direct_hydro_effective"][overridden] = 0
    terms["boundary_hydro_effective"][overridden] = 0
    observed_increment = state_next[2]-state[2]
    noise_increment = observed_increment-total*p.dt
    return dict(terms=terms, deterministic_total=total,
                observed_increment=observed_increment, noise_increment=noise_increment,
                noise_z=noise_increment/(p.In*np.sqrt(p.dt)), overridden=overridden,
                wall_distance=wall_distance, visible=vn, rdot=rdot)
