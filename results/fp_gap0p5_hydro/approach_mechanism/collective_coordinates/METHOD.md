# Method

## Coordinates and chamber convention

Right-origin configurations are reflected through the vertical divider before
constructing signed quantities. In the reflected frame, inward is positive x.

- `psi = wrap(theta_group - alpha_gap)`, where `theta_group` is the argument of
  the mean heading and `alpha_gap` points from the school center to `(0,0)`.
- `Lambda = sqrt(C_xx)/abs(x_center)`, with `C` the population covariance.
- `ell` is the signed mean normalized angular momentum using the saved full
  velocity `rdot`; `M = abs(ell)`.
- `A = sqrt(lambda_max(C)/lambda_min(C))`.

## Event/control design

All qualifying first-crossing events in seeds 0--19 are used after the 200-time
unit burn-in. Events must have at least 2.2 time units of uninterrupted
same-chamber occupancy and are thinned by 5 time units. A prospective control
pseudo-anchor must keep every fish in the same chamber from -2.2 through +3.0,
remain over 5 units from every event A, and share the event's seed and origin
chamber. Nearest matching uses only the unit-circle representation of vertical
oscillation phase derived from center height and vertical center velocity.
Heading, extension, angular momentum, aspect ratio, crossing outcome, and any
post-crossing measurement are excluded from matching and predictors.

There are 456 event/control pairs. Statistics use fixed leads A-2 and A-1.
Logistic comparisons train on seeds 0--14 and are evaluated on untouched seeds
15--19. Models contain only linear standardized versions of the prespecified
coordinates; no neural network, flux variable, or surrogate potential is used.

## Equation replay

One temporally central event/control pair per seed is replayed at both leads.
The existing validated replay reproduces Delaunay visibility, alignment,
attraction, direct dipolar hydrodynamics, the saved BEM boundary solve, exact
rounded-tip wall avoidance, clipping, and simulator overrides.

For each microscopic angular contribution `omega_k`, its contribution to the
group heading is calculated analytically from `Z = mean(exp(i theta_j))`:

```
d arg(Z)/dt = Im(conj(Z) mean(i exp(i theta_j) omega_kj)) / abs(Z)^2.
```

The rotation of the COM-to-gap direction is computed separately from
self-propulsion, direct hydrodynamic velocity, and BEM velocity. For each
translational component, `dLambda/dt` follows by differentiating population
variance and center-to-divider distance. Reported torque summaries are robust
medians because direct hydrodynamics can produce rare large angular rates.

## Memory test

Full-occupancy origins are sampled every 0.2 time units. Linear increment
models at lag 0.1 are trained on seeds 0--14 and tested on 15--19. The baseline
uses current `(sin psi, cos psi, Lambda)`; comparisons add phase, preceding
`(Delta psi, Delta Lambda)`, or both. Separate state/history models in four
phase quadrants test whether fixed-phase sampling restores closure.
