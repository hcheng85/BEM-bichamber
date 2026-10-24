"""Exact series calculation for a deliberately reduced angular-gate model.

This is NOT the full BEM-fish simulator.  The swimmer first undergoes circular
Brownian angular preparation during t_prepare, then must remain inside a
stationary geometric angular cone during t_flight.  The first part is a wrapped
heat kernel; the second is absorbing Brownian survival in an interval.
"""

import argparse
from functools import lru_cache
import math
import numpy as np
from scipy.optimize import minimize_scalar


@lru_cache(maxsize=8)
def _quadrature(nodes):
    return np.polynomial.legendre.leggauss(nodes)


def half_cone(gap: float, tip_radius: float, distance: float, height=0.0):
    """Exact tangent half-angle for symmetric rounded tips when height=0."""
    c_minus = -gap / 2 - tip_radius
    c_plus = gap / 2 + tip_radius
    diffs = (c_minus - height, c_plus - height)
    rho = [math.hypot(distance, dy) for dy in diffs]
    if distance <= 0 or min(rho) <= tip_radius:
        raise ValueError("No admissible tangent cone for this location")
    lo = math.atan2(diffs[0], distance) + math.asin(tip_radius / rho[0])
    hi = math.atan2(diffs[1], distance) - math.asin(tip_radius / rho[1])
    if lo >= hi:
        raise ValueError("Geometric tangent cone is closed")
    return (hi-lo)/2, (hi+lo)/2


def angular_gate_probability(D, theta_initial, a, t_prepare, t_flight,
                             pre_drift=0.0, flight_drift=0.0,
                             nodes=512, n_modes=100):
    """Exact angular-only success probability up to quadrature/series error.

    Preparation occurs on the circle with constant drift pre_drift.  After the
    gate opens, the angle diffuses with constant drift flight_drift, killed on
    exiting (-a,a).  This is stricter than mere spatial crossing.
    """
    if D <= 0 or not 0 < a < math.pi or t_prepare <= 0 or t_flight <= 0:
        raise ValueError("Require D>0, 0<a<pi, positive times")
    xg, wg = _quadrature(nodes)
    x = a * xg
    w = a * wg
    mean = theta_initial + pre_drift * t_prepare
    # Wrapped Gaussian.  For D*t<~30, tails outside +/-5 turns are negligible.
    winding = np.arange(-6, 7)[:, None]
    dx = x[None, :] - mean + 2 * math.pi * winding
    p0 = np.exp(-dx**2 / (4 * D * t_prepare)).sum(axis=0)
    p0 /= math.sqrt(4 * math.pi * D * t_prepare)
    if flight_drift == 0.0:
        k = np.arange(n_modes)
        n = 2*k + 1
        eig = np.exp(-D*(n*math.pi/(2*a))**2 * t_flight)
        coeff = (4/ math.pi) * ((-1.)**k) / n * eig
        survival = np.sum(coeff[:, None] * np.cos(np.outer(n*math.pi/(2*a), x)),axis=0)
    else:
        beta = flight_drift/(2*D)
        n = np.arange(1, n_modes+1)
        nu = n*math.pi/(2*a)
        # Sturm-Liouville eigenfunction solution of b*S_x+D*S_xx=S_t.
        amp = (nu/a)*(np.exp(-beta*a)-(-1.)**n*np.exp(beta*a))/(beta**2+nu**2)
        weights = amp*np.exp(-(D*nu**2+flight_drift**2/(4*D))*t_flight)
        survival = np.exp(-beta*x)*np.sum(
            weights[:,None]*np.sin(np.outer(nu,x+a)),axis=0)
    return float(np.dot(w, p0 * survival))


def calculate(gap=7.158910531638177, tip_radius=.5, distance=2.34,
              height=0., angular_mismatch=.5, t_prepare=2., t_flight=None,
              pre_drift=0., flight_drift=0.):
    a, center = half_cone(gap, tip_radius, distance, height)
    max_abs_heading = max(abs(center-a), abs(center+a))
    if max_abs_heading >= math.pi/2:
        raise ValueError("No uniform positive self-propelled x speed inside cone")
    safe_flight_bound = distance / math.cos(max_abs_heading)
    if t_flight is None:
        t_flight = safe_flight_bound
    theta0 = center + a + angular_mismatch
    obj = lambda logD: -angular_gate_probability(
        math.exp(logD), theta0-center, a, t_prepare, t_flight,
        pre_drift=pre_drift, flight_drift=flight_drift)
    opt = minimize_scalar(obj, bounds=(math.log(.0001), math.log(3.)),
                          method="bounded", options={"xatol":1e-8})
    Dopt = math.exp(opt.x)
    # Matched exponential estimate is meaningful only when zero-noise steering
    # still leaves the heading outside the corridor.
    delta = angular_mismatch
    final_centered_heading = (theta0-center+pre_drift*t_prepare+math.pi) % (2*math.pi)-math.pi
    unmet_gap = max(0., abs(final_centered_heading)-a)
    Dbal = a * unmet_gap / (math.pi * math.sqrt(t_prepare*t_flight))
    print(f"gap={gap:.6g}, d={distance:.4f}, height={height:.2f}, a={a:.5f}, "
          f"theta0-center={a+delta:.5f}, tprep={t_prepare:.2f}, "
          f"tflight={t_flight:.2f}, safe_tflight={safe_flight_bound:.2f}, "
          f"predrift={pre_drift:.3f}, flightdrift={flight_drift:.3f}")
    print(f"Dopt={Dopt:.6g}, Inopt={math.sqrt(2*Dopt):.6g}, "
          f"Pmax={-opt.fun:.7g}, balance_D={Dbal:.6g}, "
          f"balance_In={math.sqrt(2*Dbal):.6g}, unmet_gap={unmet_gap:.5f}")
    for In in [.05,.1,.2,.35,.5,.7,1,1.5,2.]:
        P = angular_gate_probability(.5*In*In, theta0-center, a,
                                     t_prepare, t_flight, pre_drift=pre_drift,
                                     flight_drift=flight_drift)
        print(f"{In:.3f},{P:.10g}")
    # Series/numerical integration convergence at the maximizer.
    for nodes, modes in [(128,30), (256,50), (512,100), (1024,150)]:
        P = angular_gate_probability(Dopt,theta0-center,a,t_prepare,t_flight,
                                     pre_drift=pre_drift,flight_drift=flight_drift,
                                     nodes=nodes,n_modes=modes)
        print(f"convergence {nodes:4d} {modes:3d}: {P:.12g}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--distance", type=float, default=2.34)
    parser.add_argument("--height", type=float, default=0.)
    parser.add_argument("--mismatch", type=float, default=.5)
    parser.add_argument("--prepare", type=float, default=2.)
    parser.add_argument("--flight", type=float, default=None)
    parser.add_argument("--pre-drift", type=float, default=0.)
    parser.add_argument("--flight-drift", type=float, default=0.)
    args=parser.parse_args()
    calculate(distance=args.distance,height=args.height,
              angular_mismatch=args.mismatch,t_prepare=args.prepare,
              t_flight=args.flight,pre_drift=args.pre_drift,
              flight_drift=args.flight_drift)
