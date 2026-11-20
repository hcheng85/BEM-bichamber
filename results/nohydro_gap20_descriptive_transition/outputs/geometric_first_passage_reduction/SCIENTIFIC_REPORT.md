# Can source-school geometry support a predictive first-passage theory?

## Decision

The source school has an excellent coherent-aggregate **translational**
description, but the proposed state `(d, psi)` does not provide a closed,
predictive first-passage theory.

Two assumptions fail for identifiable physical reasons:

1. The COM is not the correct encounter surface for an extended school.  A
   leading portion of the school reaches the established approach rectangle in
   every crossing episode, whereas the COM enters it in only 84% of episodes.
2. The collective turning rate is not determined by COM distance and relative
   bearing.  It retains the individual positions and visual neighborhoods in
   the production social and wall rules.  This missing social contribution
   remains outcome-separating on independent seeds after controlling `q`, `d`,
   and `psi`.

Thus `(d,psi)` is a useful kinematic projection and a good short-time
conditional description **when the microscopic torque is supplied**, but it
cannot presently predict attempt frequency or commitment probability from its
own state.  No fitted turning law or extra state hierarchy is introduced to
hide this failure.

## Exact crossing-free equations

For the fish currently in the original source chamber, define

```text
R = source-school COM,               d = |R|,
g = -R/d,                            h = mean_i e_i = P m,
theta = angle(m),                    chi = angle(g),
psi = theta - chi.
```

Let `V=mean_i v_i` be the actual stored source velocity, including collision
overrides, and let `g_perp` be `g` rotated counterclockwise by 90 degrees.
Between crossings, membership is fixed and the exact kinematics are

```text
dR/dt   = V,
dd/dt   = -g . V,
dg/dt   = -(I - g g^T) V/d,
dchi/dt = -g_perp . V/d,
dpsi    = dtheta + (g_perp . V/d) dt.
```

The production angular equation is

```text
d theta_i = Omega_i dt + In dW_i,       In = 0.5,
Omega_i = wall_i + social_i + both production clipping corrections.
```

Ito's formula for the direction of the mean heading gives

```text
dtheta = [sum_i cos(theta_i-theta) Omega_i /(N_s P) + B_Ito] dt
         + In/(N_s P) sum_i cos(theta_i-theta) dW_i,

B_Ito = -In^2/(N_s^2 P^2)
         sum_i (m_perp.e_i_perp)(m.e_i_perp).
```

The wall, social, and clipping parts of the collective angular rate are
therefore exact weighted projections of the corresponding production torques.
The angular-noise quadratic variation is

```text
d[theta]_t = In^2/(N_s^2 P^2)
             sum_i cos^2(theta_i-theta) dt.
```

If the coherent-aggregate approximation `V approximately P m` is valid, the
two geometric equations reduce to

```text
dd/dt   approximately -P cos(psi),
dpsi    approximately [omega_collective + P sin(psi)/d] dt
                        + collective angular noise.
```

The `P sin(psi)/d` term is geometric bearing rotation, not a torque.  For a
misaligned straight path it increases `|psi|`; active turning must compensate
it to keep the school directed toward the opening.

At chamber crossings, membership and these collective quantities jump.  The
present first-passage audit terminates before that jump and therefore does not
fold membership transfer into a smooth equation.

## Numerical identity and velocity audits

The analysis reuses the fixed 337-event cohort from 37 original seeds and the
unchanged archived production source.  It evaluates six strictly pre-crossing
leads and performs no new simulation.

- Maximum stored-position update residual: `6.46e-17`.
- Maximum relative-bearing identity residual: `8.89e-13` in rate units.
- Native-step continuous-kinematic RMSE: `1.46e-4` for `dd/dt` and `5.94e-5`
  for bearing rotation; these are finite-step curvature errors at `dt=0.001`.

The aggregate velocity approximation passes unusually strongly.  On the ten
held-out seeds, mean `cos(angle(V,m))` is at least `0.9999986` at leads 2, 1,
and 0.5.  Mean `|V-h|` is below `4.9e-4`.  Source polarization is also tightly
concentrated: its all-event mean remains `0.966--0.969`, with fifth percentiles
above `0.931` at every evaluated lead.

Consequently the parameter-free unit-speed approximation
`dd/dt=-cos(psi)` has held-out RMSE `0.0200`, compared with an actual
`dd/dt` standard deviation of `0.338`.  Over a 0.1-time-unit horizon, frozen
instantaneous translation predicts the observed distance rate with held-out
correlation `0.9993` and RMSE `0.0157`.  Translation is not the obstruction to
the reduced theory.

## Turning and the microscopic interactions

Using the **actual instantaneous microscopic torque**, the frozen-rate
equation predicts the next 0.1-time-unit `psi` rate on held-out seeds with
correlation `0.827`.  Its RMSE is `0.193`, versus an observed standard deviation
of `0.342`.  This is meaningful short-time predictability, but it is not a
reduced prediction: computing that rate required every fish position, heading,
and visual neighborhood.

The physical torque audit explains the limitation.  In a perfectly compact
school one might replace all individual wall rays by a wall torque evaluated
at `(R,theta)`.  That controlled approximation fails here:

```text
held-out compact-wall correlation       0.690
held-out compact-wall RMSE              0.143
actual wall-rate standard deviation     0.075
median school RMS radius                1.05
median COM distance from opening        2.73.
```

The school is not point-like on the relevant geometric scale.  Different
members see different wall intersections.  The social torque is even less
reducible: the production `GetVision` rule depends on visual neighbors,
relative positions, and headings, none of which are determined by `R` and
`m`.  High polarization makes `V approximately Pm`; it does not determine the
weighted sum of social torques.

To test this empirically without fitting a drift, torque contributions were
projected onto the reflection-invariant alignment rate
`-sin(psi) omega`.  Continued and failed events were then compared inside the
same fixed `q`, `d`, and angular-sector strata.  On held-out seeds:

| lead | social alignment-rate difference | wall difference | geometric difference |
|---:|---:|---:|---:|
| 2.0 | +0.043 `[0.002,0.084]` | +0.012 `[0.003,0.031]` | +0.031 `[0.018,0.049]` |
| 1.0 | +0.067 `[0.037,0.102]` | -0.010 `[-0.038,0.013]` | +0.033 `[0.010,0.056]` |
| 0.5 | +0.026 `[-0.013,0.074]` | -0.011 `[-0.049,0.019]` | +0.013 `[-0.011,0.035]` |
| 0.2 | +0.042 `[-0.011,0.082]` | -0.036 `[-0.069,-0.008]` | -0.003 `[-0.044,0.019]` |

The lead-1 social separation is independently reproduced in training seeds
(`+0.060`, `[0.037,0.088]`).  It is also larger than the residual geometric
separation there.  States with comparable coarse `(q,d,psi)` therefore do not
have the same collective turning environment.  This is a structural closure
failure, not merely excess angular noise.

## COM-entry sensitivity

The established approach rectangle is used unchanged.  No new distance or
occupancy threshold is selected.  At 0.1-time-unit resolution during the last
two crossing-free units:

| detection audit | all events | held-out events |
|---|---:|---:|
| at least one source fish enters the rectangle | 100% | 100% |
| source COM enters the rectangle | 83.7% `[79.6,87.5]` | 84.5% `[78.3,89.5]` |
| matched by the previous prospective COM attempt detector | 64.7% `[58.8,71.7]` | 63.6% `[56.7,71.6]` |

The 55 episodes in which the COM never enters still have a median maximum of
27 source fish simultaneously in the approach rectangle (range 1--42).
Therefore this is not generally a single isolated scout.  A finite leading
part can initiate the physical encounter while the centroid remains outside.

The previous COM-attempt frequency is consequently a reproducible rate for
that specific COM event surface, but not the full physical approach frequency.
Changing the COM threshold would only move the bias and is not pursued.

## Independent-pathway evidence

Successful held-out episodes approach faster even after the same coarse
`q,d,psi` controls.  The continued-minus-failed `dd/dt` differences are
`-0.040` at lead 2 `[-0.076,-0.008]`, `-0.092` at lead 1
`[-0.153,-0.035]`, and `-0.108` at lead 0.5 `[-0.173,-0.061]`.
Together with the social-alignment contrasts, this confirms that the observed
curved encounter is predictive on independent trajectories.

It does **not** establish `(d,psi)` closure.  The successful path differs in
both motion through the geometric state and in a microscopic social-turning
input that the state does not specify.

## Consequence for a first-passage or renewal theory

A closed point-aggregate first-passage problem would require a drift and noise
law for

```text
(d,psi) -> (-P cos psi,
            omega(d,psi) + P sin(psi)/d + angular noise)
```

and an absorbing encounter boundary defined by COM position.  The data reject
both the autonomous `omega(d,psi)` assumption and the COM absorbing boundary.
Fitting an empirical `omega` would conceal the unresolved visual-neighborhood
physics and would violate the first-principles objective.

The simplest physically defensible alternative is not another geometric
threshold.  Attempt initiation should be tied to the already established
fish flux/storage in the physical approach region, which detects the leading
part and obeys the exact conservation law.  Its persistence must then be
coupled to source-school orientation through a derived finite-extent social
turning closure.  The existing approach storage `b` is the natural quantity to
test in that role; it is not a newly invented coordinate.  Until that angular
closure is derived or controlled, a renewal identity can organize attempt and
commitment statistics but cannot predict their rates from first principles.

## Status of claims

- **Exact:** the `(d,psi)` kinematic equations and production-torque projection.
- **Controlled and supported:** `V approximately Pm` and nearly unit-speed
  translation during approaches.
- **Conditionally predictive:** frozen microscopic torque over 0.1 time units.
- **Empirical, independent-seed evidence:** faster convergence and additional
  social alignment distinguish continuation.
- **Rejected:** point-COM encounter detection and autonomous `(d,psi)` turning
  closure.
- **Not constructed:** a fitted potential, phenomenological turning law,
  transition matrix, or new microscopic simulation.
