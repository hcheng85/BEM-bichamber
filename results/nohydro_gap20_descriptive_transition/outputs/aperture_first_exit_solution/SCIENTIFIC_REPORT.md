# Single-attempt aperture capture from a kinetic first-exit problem

## Conclusion

The parameter-free single-particle problem can be solved and has a useful but
restricted physical domain.  For a fish already moving inward in the
source-side approach rectangle, angular diffusion broadens the deterministic
aperture ray into a capture probability.  This prediction is approximately
calibrated for the shortest archived approaches (initial ballistic exit time
at most 0.25--0.5), but it overpredicts capture on longer approaches.  Freezing
the exact instantaneous microscopic turning rate reduces the overall
probability bias slightly but severely degrades exit-time prediction.  Thus the
failure is not missing aperture geometry: during a long attempt the many-fish,
wall, and visibility-dependent turning environment evolves and cannot be
represented by zero or frozen constant torque.

This is a theory of capture in **one approach attempt**.  A lateral or upstream
exit is absorbing only for this mathematical attempt, not for the physical
chamber.  The same fish may return in a later, statistically dependent renewal.
The result is therefore not an eventual chamber-transition probability.

## Backward problem and numerical solution

Use source-fixed coordinates, with `z=0` at the opening, `z=-ell` at the back
of the approach layer, and `|y|<h`.  Here
`ell=2.8635642126552705`, `h=ell/2`, unit swimming speed, and the production
angular-noise amplitude is `In=0.5`.  For prescribed constant turning rate
`Omega`,

```text
dz = cos(phi) dt,
dy = sin(phi) dt,
dphi = Omega dt + In dW.
```

The one-attempt capture committor solves

```text
L q = 0,
L = cos(phi) d_z + sin(phi) d_y + Omega d_phi
    + (In^2/2) d_phi_phi.
```

On the outgoing part of `z=0`, `q=1`; on the outgoing lateral and upstream
boundaries, `q=0`; angle is periodic.  The mean first-exit time `T` and the
capture-weighted time `m=E[tau 1_capture]` solve `L T=-1` and `L m=-q` with
zero outgoing boundary data.  The conditional capture time is `m/q`.

The implementation uses a monotone continuous-time lattice generator with
upwind spatial/angular drift and centered angular diffusion.  Sparse LU gives
`q`, `T`, and `m`; no coefficient is estimated from trajectories.  Grid sizes
were 15x9x16, 23x13x24, and 31x19x36.  At `In=0.5`, the largest change among
five off-grid committor probes fell from 0.0539 to 0.0151 on the last
refinement.  The final matrix has 21,204 states and 104,220 nonzeros.  This is
a systematic numerical approximation, not a claimed closed-form finite-noise
solution.

## Straight characteristics and leading angular-noise correction

For zero torque and `u=cos(phi)>0`, the deterministic time and transverse
intercept at the opening are

```text
tau = -z/u,
y_hit = y - z tan(phi).
```

The exact zero-noise capture indicator is `1(|y_hit|<=h)`, provided the gap is
reached before a lateral boundary.  Linearizing in angular noise gives

```text
delta y_hit = (In/u) integral_0^tau W_s ds,
Var(delta y_hit) = In^2 tau^3/(3 u^2).
```

Consequently the leading capture probability is the Gaussian aperture mass

```text
Phi((h-y_hit)/sigma_y) - Phi((-h-y_hit)/sigma_y),
sigma_y^2 = In^2 tau^3/(3 u^2).
```

The leading exit-time fluctuation is

```text
delta tau = In tan(phi) integral_0^tau W_s ds,
Var(delta tau) = In^2 tan(phi)^2 tau^3/3.
```

The controlled regime is `In sqrt(tau) << 1`, with `cos(phi)` bounded away
from zero, transverse spread small compared with the distance to competing
boundaries, and the deterministic path away from boundary corners.  Near an
aperture edge the smoothed probability, rather than a regular expansion of the
indicator, is the appropriate leading result.

As `In` decreases from 0.5 to 0.25, 0.125, and 0.0625, the leading Gaussian
probability's mean absolute difference from the ray indicator on the
short-flight grid falls from 0.0476 to 0.0237, 0.0119, and 0.00654, reaching
zero at `In=0`.  The numerical PDE classification error relative to the ray is
0.00264 at zero noise.  Its nonzero probability L1 plateau (0.0415) is the
expected fixed-grid upwind diffusion; it is not physical angular noise.

## Constant prescribed turning

For `Omega != 0`, the deterministic characteristics are circular arcs:

```text
phi(t) = phi_0 + Omega t,
z(t) = z_0 + [sin(phi_0+Omega t)-sin(phi_0)]/Omega,
y(t) = y_0 + [cos(phi_0)-cos(phi_0+Omega t)]/Omega.
```

An exact characteristic audit on 4,000 uniform states for each of
`Omega=-1,+1` found zero-noise PDE classification errors of 3.13% and 2.99%
among characteristics that exit.  Closed arcs trapped inside the rectangle
occupied 9.88% and 8.90% of the samples.  They make the zero-diffusion limit
singular: any positive angular diffusion eventually releases them, whereas a
deterministic closed orbit has no first exit.  Reported PDE values on this
trapped set are numerical leakage and have no deterministic interpretation.

For a deterministic downstream exit at time `tau_0`, write
`phi_d(s)=phi_0+Omega s`.  The controlled first angular-diffusion correction is

```text
delta tau = In/cos(phi_d(tau_0))
            integral_0^tau_0 sin(phi_d(s)) W_s ds,

delta y_exit = In integral_0^tau_0 K(s) W_s ds,
K(s) = cos(phi_d(s)) + tan(phi_d(tau_0)) sin(phi_d(s)).
```

Using Ito's isometry after exchanging integrations,

```text
Var(delta y_exit) = In^2 integral_0^tau_0
                    [integral_r^tau_0 K(s) ds]^2 dr,

Var(delta tau) = In^2/cos(phi_d(tau_0))^2 integral_0^tau_0
                 [integral_r^tau_0 sin(phi_d(s)) ds]^2 dr.
```

Gaussian aperture integration around the arc endpoint is therefore the
minimal finite-persistence correction.  It reduces to the straight formulas
as `Omega` tends to zero.  It is controlled only away from tangent exits,
corners, trapped arcs, and times over which the prescribed torque changes.

## Archived-trajectory validation

The prospective held-out partition was preserved: original seeds 30--39.
From 376 source-fixed prospective configurations, 6,115 fish were initially
inside the approach rectangle and moving inward.  All left the rectangle
within five time units; 427 (6.983%) exited through the opening.  Uncertainty
was bootstrapped over the ten original seeds, never over fish as independent
replicates.

For zero torque, the full-cohort mean predicted capture was 0.1222
(95% seed-bootstrap CI 0.1000--0.1448), compared with observed 0.06983
(0.04524--0.09475).  The calibration bias was +0.05242
(0.04529--0.05940), and exit-time MAE was 0.1931
(0.1799--0.2086).

The initial ballistic duration isolates the approximation's physical regime:

| maximum ballistic time | attempts | observed | zero-torque prediction | bias (95% CI) | exit-time MAE |
|---:|---:|---:|---:|---:|---:|
| 0.25 | 639 | 0.02034 | 0.02191 | +0.00157 (-0.00537, 0.00723) | 0.1058 |
| 0.50 | 1,318 | 0.01669 | 0.02125 | +0.00455 (0.00106, 0.00766) | 0.1130 |
| 1.00 | 2,566 | 0.03079 | 0.04274 | +0.01196 (0.00843, 0.01591) | 0.1408 |
| 2.00 | 4,581 | 0.05435 | 0.08625 | +0.03190 (0.02498, 0.03814) | 0.1826 |
| 5.00 | 6,115 | 0.06983 | 0.12224 | +0.05242 (0.04539, 0.05966) | 0.1931 |

Thus the shortest-flight result is consistent with calibration, while the
failure grows monotonically as the fish has time to turn in its evolving
environment.  Selection by initial ballistic time uses only the instantaneous
prospective state, not future trajectory information.

As a controlled extension, the exact instantaneous production turning torque
was frozen for each fish and solutions at `Omega=-3,-2,...,3` were
interpolated (only 3/6,115 values required clipping).  It reduced the global
probability bias to +0.03683 (0.02713--0.04529), but the Brier-score difference
from zero torque was unresolved: zero minus frozen = -0.000104
(-0.00310, 0.00267).  More decisively, frozen torque increased exit-time MAE
from 0.1931 to 0.5655; the difference was -0.3724
(-0.3950, -0.3507).  A constant instantaneous torque is therefore not a valid
long-flight correction.  The likely source is time evolution of social input,
wall response, and visibility as the fish and school move; the present test
does not separate those causes.

## What is established and what remains

Established analytically and numerically:

- the parameter-free backward first-exit problem and its committor/time
  solution for the production angular noise;
- recovery of the deterministic ray classification as noise vanishes;
- the leading short-flight angular-diffusion correction and its validity
  conditions;
- deterministic circular-arc characteristics, their perturbation, and the
  trapped-orbit singularity;
- quantitative short-flight agreement, followed by a resolved long-flight
  failure in independent archived seeds.

The remaining collective transition rate cannot be obtained by treating this
committor as an eventual switching probability.  A first-principles theory
must additionally derive (i) the source-conditioned rate and state
distribution of entries into the approach layer, (ii) the time-dependent
conditional torque/visibility process during an attempt, and (iii) the
renewal kernel for lateral/upstream exits and correlated re-entry, coupled to
the changing chamber population.  Only then can approach-arrival flux and
single-attempt capture be composed without double-counting repeated fish.
No fitted collective potential or extra empirical coordinate is supported by
this calculation.

## Reproducibility

`solve_aperture_first_exit.py` contains the solver, self-test, convergence and
characteristic audits, archived validation, and figure generation.  The
production source was read without modification (SHA-256
`5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`).
No new microscopic trajectories were generated and no fitted coefficient was
used.  Numerical summaries in this directory are sufficient to regenerate
the conclusions; raw trajectories and intermediate trajectory samples are
excluded from version control.
