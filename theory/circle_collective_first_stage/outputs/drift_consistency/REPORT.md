# Validation of stationary polarization drift and the seed-9 transition

## Decision

All validation phases specified in `NEXT_VALIDATION_TASK.md` were completed in
order.  Phase 1 validates the source-exact angular update and noise convention.
Phase 2 shows that the earlier positive milling drift was a sparse-snapshot
sampling bias: time-uniform replay is consistent with stationary torque balance.
Phase 3 finds loss of coherent signed circulation before the polarization
collapse, but no distinct prospective impulse in the scalar projected drift
`bQ_total`.

The exact polarization-vector generator is valid, but a closed response
`F(p) = A p` is **not supported** by held-out seeds.  Consequently the fitted
eigenvalue of a pooled response matrix must not be called a stability boundary,
and neither a one-dimensional potential nor a `Q`-only transition mechanism is
justified.

## Immutable inputs

- Ten original circular full-BEM trajectories, seeds 0--9.
- Production source `BEM_circular_parallel_instance.py`, SHA-256
  `2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`.
- `N=100`, `Ip=4`, `In=0.5`, `Iw=0.94`, `If=0.01`, `dt=0.001`,
  nondimensional `R=7.158910531638177`.

The production simulator, forward-ray wall interaction, MAT files, LU factors,
and existing reports were not modified.

## Phase 0: compact-data diagnostic

The committed diagnostic reproduces the original inconsistency.  Full-series
stationary deterministic targets are 0.00804 and 0.01074 for milling seeds 7
and 8, whereas the previous ten selected snapshots gave 0.02389 and 0.06503.
This establishes a sampling question; it is not evidence that milling is
unstable.

## Phase 1: angular-update audit -- passed

The source applies

```text
omega = clip(w_direct + w_BEM + w_wall + w_social, -pi/dt, pi/dt)
theta_next = theta + omega*dt + In*sqrt(dt)*Normal(0,1)
```

If a proposed translational step exits the circle, or if an extreme velocity
safeguard is triggered, direct and BEM angular terms are removed and the update
uses `clip(w_wall+w_social)`.  The later positional safety retry changes
translation only, not the angular increment.  Noise is independently Gaussian
for each fish and step with variance `In**2*dt = 0.00025`.

The audit used 1,050 source-exact steps from milling seeds 7, 8 and the
post-transition interval of seed 9, giving 105,000 fish increments:

| diagnostic | observed | expected |
|---|---:|---:|
| residual mean | 3.38e-5 | 0 (SE 4.88e-5) |
| residual variance / dt | 0.24911 | 0.25 |
| standardized skewness | -0.0023 | 0 |
| excess kurtosis | 0.0119 | 0 |
| normal KS p-value | 0.348 | not rejected |
| maximum native-lag 1--10 correlation | 0.0134 | approximately 0 |

Residual variance divided by `dt` ranges from 0.2453 to 0.2528 across forward
wall-distance quartiles.  The maximum social-split error is `1.38e-13` and the
maximum stored-velocity replay error is `8.88e-16`.  No clip, collision,
extreme-velocity, or sticky override occurred in this Phase-1 sample.  Every
prespecified acceptance check passes.

## Phase 2: stationary drift consistency -- passed

Snapshots are uniform in time, not selected by `Q` or by proximity to a
transition.  Seeds 7 and 8 use 400 states each.  Schooling seeds 0--6 use 200
states each; seed 9 uses separate 400-state pre- and post-transition intervals.
Uncertainty is calculated from ten trajectory-time blocks for each interval.

### Milling balance

| seed | stationary target | replayed deterministic source | 95% block interval | exact finite-step total |
|---|---:|---:|---:|---:|
| 7 | 0.00804 | 0.00951 | [-0.00013, 0.01914] | -0.00054 |
| 8 | 0.01074 | 0.01530 | [-0.00010, 0.03071] | 0.00227 |

The replay means converge toward the targets as the number of snapshots grows
from 25 to 400.  For seed 7 the 200- and 400-state estimates are 0.01793 and
0.00951; for seed 8 they are 0.02006 and 0.01530.  Thus the previous ten-state
positive drift is explained by nonuniform sparse sampling.

The milling torque balance remains physically informative.  Alignment creates
polarization (`+0.03030`, `+0.03913` for seeds 7 and 8), while direct
fish--fish hydrodynamics opposes it (`-0.01814`, `-0.02011`).  Attraction is a
smaller negative term, and wall and BEM angular contributions are nearly zero
in the time average.

### Continuous generator versus the simulated finite step

The continuous Itô generator and the exact Euler-step expectation are both
reported.  Given deterministic headings
`theta_det = theta + omega*dt`, the exact conditional expectation is

```text
E[Q_next | state] = 1/N + exp(-In**2*dt) * (Q_det - 1/N).
```

All eleven separated intervals have their observed endpoint rate inside the
95% trajectory-block interval of this exact finite-step drift.  The continuous
stationary target is also inside the deterministic-source interval for the two
stationary milling seeds.

This distinction matters in schooling: rare large forward-wall torques can
produce a non-negligible finite-`dt` correction even without activating the
hard angular clip.  Two clipped fish occur in the 200-state grids of seeds 3
and 6; no collision, extreme-velocity, or sticky override occurs.  The exact
finite-step audit passes for every schooling and milling interval, including
the separated seed-9 intervals.

## Phase 3: seed-9 transition -- exploratory

Source-exact torques were evaluated every 0.1 time unit over `t=300..325`.
Stable-schooling windows from seeds 0, 1 and 4 provide a descriptive reference.
Using the fixed rule “at least 7 of 10 successive samples beyond three baseline
standard deviations,” the chronology is:

| quantity | first sustained change |
|---|---:|
| mean fish radius decreases | 308.1 |
| signed circulation `M` decreases | 308.6 |
| radial spread increases | 311.4 |
| `P` and `Q` decrease | 311.6 |
| original milling `Pi` increases | 312.2 |
| neighbor heading coherence decreases | 313.6 |
| sustained negative `bQ_total` excursion | not detected |

The mean `abs(M)` over `t=308..312` is approximately 0.740, below the 2.5%
quantile 0.828 of the stable-schooling controls.  In contrast, the center/radial
motion and the large instantaneous direct-hydrodynamic and total-`bQ`
fluctuations overlap the control ranges.  Alignment remains strong until the
polarization collapse, and neighbor coherence decreases later.

Direct evidence therefore supports this limited ordering in the single event:

```text
inward orbital geometry change -> loss of coherent signed circulation
-> polarization collapse -> high Pi / reduced neighbor coherence
```

It does not identify a unique torque impulse that initiates the event.  With
only one transition, this ordering is exploratory and cannot establish a rate,
causal law, or universal precursor.

## Test of the proposed polarization response matrix

For the 800 time-uniform milling states from seeds 7 and 8, the most restrictive
rotationally equivariant linear response was tested:

```text
F(p) = E[sum_i e_perp_i Omega_i / N | p] = a p + b J p.
```

The pooled coefficients are `a=0.1315`, `b=0.0793`; the trajectory-block 95%
interval for `a` is `[0.0521, 0.2210]`.  Although the point estimate is close
to the angular-noise threshold `In**2/2 = 0.125`, it is not a usable stability
result:

- pooled in-sample `R2` against a zero-response baseline is only 0.005;
- training seed 7 and testing seed 8 gives `R2=-0.233`;
- training seed 8 and testing seed 7 gives `R2=-0.185`;
- training seeds 7+8 and testing post-transition seed 9 gives `R2=-0.022`;
- the fitted rotation coefficient reverses strongly between seeds
  (`b=-0.563` versus `+0.598`).

The instantaneous deterministic polarization source is therefore not closed by
`p` alone.  The apparent pooled `a` is a weak average of broad, seed-dependent
responses, not a first-principles eigenvalue.

## Defensible analytic direction

The validated exact starting point is

```text
dp = [F - In**2*p/2] dt + (In/N) sum_i e_perp_i dW_i,
F = sum_i e_perp_i Omega_i / N.
```

A plausible next closure must retain at least a dynamic circulation/orbital
state and probably radial geometry, for example `F(p, M, g)` with `g`
representing the school-wall orbital configuration.  It should be tested with
controlled perturbations of steady milling configurations at fixed geometry
and circulation phase.  The existing natural trajectories do not identify a
time-independent matrix `A`, and no potential should be constructed from them.

## Files

- `angular_update_audit.json`: Phase-1 source/noise audit.
- `stationary_balance_by_seed.csv`: continuous and exact finite-step balances.
- `stationary_balance_convergence.csv/png`: milling convergence.
- `phase2_acceptance.json`: explicit acceptance checks.
- `transition_seed9_torque_timeline.csv/png`: Phase-3 chronology.
- `transition_chronology.json`: fixed-rule change times and control ranges.
- `transition_seed9_control_comparison.png`: event/control comparison.
- `polarization_response_matrix.json/png`: held-out linear-response test.
- `uniform_replay_snapshots.csv`, `polarization_response_snapshots.csv`:
  compact snapshot summaries; no per-fish or raw trajectory data.
