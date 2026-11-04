#!/usr/bin/env python3
"""Validate exact circular orbits against the original Cartesian ray-wall law.

No fit, social coupling, hydrodynamics or custom boundary override is used.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov


def cartesian_ray_wall(t: float, y: np.ndarray, radius: float, strength: float) -> np.ndarray:
    """Wall avoidance from nearest *forward* intersection of a ray and a disk."""
    del t
    x, ycoord, theta = y
    heading = np.array([np.cos(theta), np.sin(theta)])
    pos = np.array([x, ycoord])
    projection = float(np.dot(pos, heading))
    discriminant = projection**2 + radius**2 - float(np.dot(pos, pos))
    if discriminant <= 0:
        raise ValueError("Position exited the smooth interior circle")
    distance = -projection + np.sqrt(discriminant)
    point = pos + distance * heading
    wall_angle = np.arctan2(point[1], point[0])
    incidence = np.arctan2(np.sin(theta - wall_angle), np.cos(theta - wall_angle))
    wall_torque = strength * np.sign(incidence) / distance
    return np.array([heading[0], heading[1], wall_torque])


def orbit_parameters(radius: float, strength: float) -> tuple[float, float, float, float]:
    s = np.sqrt(1 + strength**2)
    return radius / s, radius * strength / s, s / radius, s / (radius * strength)


def linear_covariance(radius: float, strength: float, noise_amplitude: float):
    s = np.sqrt(1 + strength**2)
    gamma = s / (radius * strength)
    coupling = s**4 / (radius**2 * strength**2)
    jacobian = np.array([[0.0, -1.0], [coupling, -gamma]])
    noise_covariance = np.array([[0.0, 0.0], [0.0, noise_amplitude**2]])
    numeric = solve_continuous_lyapunov(jacobian, -noise_covariance)
    expected = np.diag([noise_amplitude**2 / (2 * gamma * coupling),
                        noise_amplitude**2 / (2 * gamma)])
    if not np.allclose(numeric, expected, rtol=2e-11, atol=2e-13):
        raise AssertionError(f"Lyapunov mismatch: {numeric} vs {expected}")
    eigenvalues = np.linalg.eigvals(jacobian)
    if not np.all(np.real(eigenvalues) < 0):
        raise AssertionError("Orbit is not locally attracting")
    return expected


def main() -> None:
    radius = 1.0
    strengths = (0.25, 0.94, 2.0)
    t_end = 50.0
    t_grid = np.linspace(0, t_end, 1800)
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.7), constrained_layout=True)
    for strength in strengths:
        theory_r, theory_d, theory_omega, _ = orbit_parameters(radius, strength)
        if theory_d <= 0.1:
            raise AssertionError("Circular orbit lies inside the simulator's safety layer")
        linear_covariance(radius, strength, noise_amplitude=0.10)
        for chirality in (+1, -1):
            solution = solve_ivp(
                cartesian_ray_wall, (0.0, t_end),
                [0.85 * theory_r, 0.0, chirality * 1.5],
                args=(radius, strength), t_eval=t_grid, method="DOP853",
                max_step=0.025, rtol=1e-10, atol=1e-12,
            )
            if not solution.success:
                raise RuntimeError(solution.message)
            x, y, _ = solution.y
            r = np.hypot(x, y)
            projection = x * np.cos(solution.y[2]) + y * np.sin(solution.y[2])
            distances = -projection + np.sqrt(projection**2 + radius**2 - r**2)
            if np.min(distances) <= 0.1:
                raise AssertionError("Validation trajectory entered the original safety layer")
            phi = np.unwrap(np.arctan2(y, x))
            tail = t_grid > 35
            measured_r = float(np.mean(r[tail]))
            measured_omega = float(np.polyfit(t_grid[tail], phi[tail], 1)[0])
            radius_error = abs(measured_r - theory_r)
            frequency_error = abs(measured_omega - chirality * theory_omega)
            print(f"Iw={strength:.2f}, chirality={chirality:+d}: "
                  f"r*= {theory_r:.6f}, measured= {measured_r:.6f}, "
                  f"omega*= {chirality*theory_omega:.6f}, measured= {measured_omega:.6f}")
            if radius_error > 1e-5 or frequency_error > 1e-5:
                raise AssertionError(f"Unacceptable orbit discrepancy: {radius_error}, {frequency_error}")
            if chirality == 1:
                (line,) = axes[0].plot(t_grid, r, lw=1.5, label=f"$I_w={strength:g}$")
                axes[0].axhline(theory_r, color=line.get_color(), ls="--", lw=1)
                if strength == 0.94:
                    axes[1].plot(x, y, lw=1.5, label="CCW from Cartesian ray law")
                    circle = plt.Circle((0, 0), theory_r, fill=False,
                                        ls="--", ec="black", lw=1.3, label="predicted orbit")
                    axes[1].add_patch(circle)
            elif strength == 0.94:
                axes[1].plot(x, y, lw=1.5, label="CW from Cartesian ray law")
    axes[0].set(xlabel="time / (R/U)", ylabel="$r/R$", title="Wall-induced orbit selection")
    axes[0].legend(frameon=False)
    angles = np.linspace(0, 2*np.pi, 500)
    axes[1].plot(np.cos(angles), np.sin(angles), color="grey", lw=2,
                 label="outer wall")
    axes[1].set(xlabel="$x/R$", ylabel="$y/R$", title="Two attracting chiral orbits")
    axes[1].axis("equal")
    axes[1].legend(loc="upper right", fontsize=8, frameon=False)
    path = Path(__file__).with_name("orbit_validation.png")
    fig.savefig(path, dpi=180)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
