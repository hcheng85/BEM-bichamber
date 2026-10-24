"""Small independent consistency tests of the analytic gate series."""

import math
from finite_window_gate import half_cone, angular_gate_probability as P


def test_tangent_symmetry_and_closed_cone():
    a, center = half_cone(7.158910531638177, .5, 2.34, 0.)
    assert abs(center) < 1e-14
    assert abs(a-.94349) < 1e-4


def test_probability_bounds_and_limits():
    theta0, a, prep, flight = 1.44349, .94349, 2., 3.99
    low = P(.0005, theta0, a, prep, flight)
    mid = P(.09, theta0, a, prep, flight)
    high = P(2., theta0, a, prep, flight)
    assert 0 <= low < mid <= 1
    assert 0 <= high < mid


def test_flight_drift_symmetry():
    positive = P(.12, 1.44349, .94349, 2., 3.99,
                 pre_drift=-.08, flight_drift=.12)
    reflected = P(.12, -1.44349, .94349, 2., 3.99,
                  pre_drift=.08, flight_drift=-.12)
    assert math.isclose(positive, reflected, rel_tol=1e-11)


def test_zero_drift_limit_of_full_eigenseries():
    exact_odd_modes = P(.12, 1.44349, .94349, 2., 3.99)
    full_modes_tiny_drift = P(.12, 1.44349, .94349, 2., 3.99,
                              flight_drift=1e-8)
    assert abs(exact_odd_modes-full_modes_tiny_drift) < 1e-8


def test_quadrature_convergence():
    a, _ = half_cone(7.158910531638177, .5, 2.34)
    left = P(.09, a+.5, a, 2., 3.99, nodes=128, n_modes=30)
    right = P(.09, a+.5, a, 2., 3.99, nodes=512, n_modes=100)
    assert abs(left-right) < 1e-9
