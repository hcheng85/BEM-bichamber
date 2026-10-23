# Reorientation, proximity, and inward translation before first-fish crossing

## Decision

The evidence supports **B at the collective-coordinate level**:

> First-fish crossing is primarily proximity-dependent stochastic nucleation;
> once an event trajectory is underway, its inward approach is driven mainly by
> self-propulsion, but a preceding collective reorientation is not independently
> resolved after initial center distance is controlled.

This is not evidence for a static potential barrier. It also does not mean that
the approach is passive: the exact kinematics show a strong active inward
translation. What fails is the first arrow in the proposed collective chain
`collective reorientation -> inward translation -> crossing`. The middle arrow
is established; collective `psi` does not identify the initiating turn at
comparable starting distance. The earlier leading-edge result therefore remains
the more plausible location of the angular gate.

## 1. What Lambda measured

For every archived event/control pair,

```
Delta log Lambda = Delta log sigma_perp - Delta log d.
```

The unit of replication is the simulation seed: pair differences are averaged
within each of 20 seeds before bootstrap resampling.

| lead | Delta log Lambda | Delta log sigma_perp | -Delta log d | proximity fraction |
|---|---:|---:|---:|---:|
| 2 | 0.0514 [0.0037, 0.0992] | -0.0393 [-0.0708, -0.0076] | 0.0906 [0.0498, 0.1327] | 176% |
| 1 | 0.2470 [0.2066, 0.2876] | 0.0383 [0.0089, 0.0693] | 0.2087 [0.1689, 0.2498] | 84.5% |

Intervals are seed-cluster bootstrap 95% intervals. The requested preliminary
values are reproduced. At lead 2, the numerator actually indicates a narrower
wall-normal school; the positive `Lambda` contrast is more than entirely caused
by events already being closer to the divider. At lead 1, genuine extension is
small but resolved, while about 85% of the mean logarithmic contrast is center
proximity.

`Lambda` should therefore not be called a deformation coordinate without
separating its denominator.

## 2. Coordinate comparison

All models are linear logistic regressions trained on seeds 0--14 and evaluated
without refitting on seeds 15--19. Circular heading enters as `(cos psi,
sin psi)`. The samples and split are identical to the preceding analysis.

| model | pooled AUC, lead 2 | pooled AUC, lead 1 |
|---|---:|---:|
| `psi` | 0.721 | 0.742 |
| `(psi, Lambda)` | 0.689 | 0.768 |
| `(psi, d)` | 0.720 | **0.805** |
| `(psi, sigma_perp)` | 0.707 | 0.739 |
| `(psi, d, sigma_perp)` | 0.721 | **0.805** |

At lead 1, the mean held-out-seed AUC gain of `(psi,Lambda)` over `psi` is
+0.0240, but its five-seed bootstrap interval is [-0.0008, 0.0494]; it is not
statistically established. The gain from adding distance is +0.0579
[0.0141, 0.1018]. Adding `sigma_perp` to `(psi,d)` changes mean seed AUC by
-0.0012 [-0.0158, 0.0134]. Thus wall-normal shape adds no independent
discrimination once heading and center distance are known.

At lead 2, adding `Lambda` actually reduces mean seed AUC by -0.0260
[-0.0452, -0.0035]. Neither distance nor spread gives an established gain over
heading there. The complete per-seed AUC table is retained rather than treating
the pooled AUCs as independent replicates.

## 3. Exact pre-crossing kinematics

All 456 event and 456 phase-matched control segments were integrated over
`A-2` to `A-1` at the native `dt=0.001`. For unit swimming speed,

```
mean_j cos(theta_j) = P cos(psi + alpha_gap),
```

so the self-propulsion term is exactly the requested collective-heading
integral. Direct dipolar velocity is reconstructed from the source equation.
The saved effective `rdot` identifies frozen and self-only overrides; the
remaining velocity is the BEM-mediated contribution actually applied by the
simulator. The three contributions reproduce every saved center displacement
with mean residual of order `10^-17`.

Negative changes mean motion toward the divider:

| contribution from A-2 to A-1 | events | controls | event - control |
|---|---:|---:|---:|
| actual change in d | -0.3034 | +0.0204 | -0.3238 [-0.3519, -0.2950] |
| self propulsion | -0.2455 | +0.0157 | -0.2612 [-0.2849, -0.2368] |
| direct hydrodynamics | -0.0593 | +0.0091 | -0.0684 [-0.0741, -0.0627] |
| effective BEM translation | +0.0014 | -0.0043 | +0.0057 [0.0044, 0.0070] |
| actual change in sigma_perp | +0.0309 | -0.0288 | +0.0597 [0.0450, 0.0735] |

Self-propulsion accounts for about 81% of the extra event inward displacement;
direct hydrodynamics adds about 21%, and effective BEM translation weakly
opposes the contrast. Genuine wall-normal widening occurs, but its length-scale
contrast is only about 18% of the center displacement contrast. The dominant
geometric change is translation, not deformation.

This establishes an exact physical statement: a gap-directed collective
heading produces inward center motion through self-propulsion, and this is what
event trajectories do during the measured interval. It does **not** establish
that collective heading was the cause that selected those trajectories before
they became closer.

## 4. Distance-matched prospective risk set

The sensitivity analysis starts exactly at `t0=A-2`. A positive outcome is a
first crossing during `(t0,t0+2]`. A control has all fish in its origin chamber
at `t0` and no crossing only over that same two-unit future horizon; it is not
required to remain event-free beyond the prediction window. Each control is
matched within seed and origin chamber on vertical phase and `d(t0)`, never on
heading.

Matching is close: median phase-angle mismatch is 0.0096 rad; the 99th
percentile is 0.0494 rad. Event distance remains only 0.0090 smaller than
control distance, compared with distances spanning roughly 1--4.6.

At matched distance, events retain a small but seed-consistent heading contrast:

```
event - control cos(psi) = +0.0599 [0.0466, 0.0744].
```

However, its predictive magnitude is weak. A heading-only model trained on
seeds 0--14 has pooled held-out AUC 0.534, and mean seed AUC 0.536
[0.531, 0.539]. Thus collective heading is statistically detectable but barely
discriminates which equally close, phase-matched configurations cross within
two units. Most of the former heading separation was coupled to configurations
that were already closer.

## Exact, observed, and hypothesized statements

**Exact consequence of the microscopic equations**

- `P cos(psi+alpha_gap)` is the inward self-propelled center velocity.
- The saved update decomposes exactly into self-propulsion, direct dipolar
  translation, effective BEM translation, and documented overrides.

**Directly supported by trajectories**

- Event schools undergo a large inward displacement from A-2 to A-1, about 81%
  of whose event/control contrast is self-propelled.
- Center proximity, not wall-normal spread, explains most of `Lambda`.
- Distance substantially improves held-out prediction at lead 1; spread does
  not add independent information.
- At comparable initial distance, collective `psi` retains only weak predictive
  discrimination.

**Remaining hypothesis**

- A small leading-edge subpopulation may turn first and pull the collective
  center inward before the global heading becomes an independently useful risk
  coordinate. This is consistent with the prior leading-edge angular-gate
  analysis, but the present collective test does not prove the initiation
  mechanism.

## Stopping criterion and bounded next test

The broad collective-coordinate investigation should stop here. A static
potential and a new parameter sweep are not justified. If the initiating turn
must be resolved, the single bounded next test is to integrate the already
reconstructed wall/social/direct/BEM angular torques for the **same fixed
leading ten fish** from `A-2` until their first entry into the tangent cone,
against distance- and phase-matched controls. That test uses existing
trajectories and directly asks whether a cumulative deterministic angular
  impulse precedes the local gate; no additional global diagnostics are needed.
