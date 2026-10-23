# Collective angular mismatch and wall-normal extension before nucleation

## Main result

The pair `(psi, Lambda)` is a meaningful **precursor description**, but not a
closed two-dimensional stochastic state. Here `psi` is the reflected
collective heading relative to the direction from the school center to the
opening, and `Lambda = sqrt(C_xx)/abs(x_center)` is wall-normal extent relative
to center-to-divider distance.

Across 456 first-crossing events and 456 prospective non-crossing controls from
all 20 independent full-BEM trajectories, event configurations point more
toward the opening at both leads. Event minus control in `cos(psi)` is +0.290
(seed-bootstrap 95% interval 0.258--0.323) at lead 2 and +0.298
(0.254--0.341) at lead 1. `Lambda` is not resolved at lead 2: +0.0125
(-0.0013--0.0266), but is clearly larger at lead 1: +0.0643
(0.0518--0.0770). Thus angular alignment appears first; wall-normal extension
becomes informative nearer the first crossing.

A simple logistic comparison trained on seeds 0--14 and evaluated only on
seeds 15--19 gives:

| coordinates | AUC, lead 2 | AUC, lead 1 |
|---|---:|---:|
| `(psi, Lambda)` | 0.689 | **0.768** |
| signed Filella angular momentum and magnitude | 0.682 | 0.681 |
| geometric aspect ratio | 0.425 | 0.648 |

The mean AUC across the five held-out seeds for `(psi, Lambda)` is 0.700
(seed-bootstrap interval 0.657--0.742) at lead 2 and 0.775
(0.729--0.824) at lead 1. An ablation shows that `psi` alone is better at lead
2 (AUC 0.721), whereas adding `Lambda` helps at lead 1 (`psi` alone: 0.742;
pair: 0.768). This is a modest but reproducible precursor, not a deterministic
transition criterion.

## Comparators

The angular-momentum definition was source-audited against Huang et al.'s use
of the Filella model: each contribution uses the saved full velocity,
`(r_j-r_c) cross rdot_j / (norm(r_j-r_c) norm(rdot_j))`, not the heading alone.
Its signed school mean does not shift before nucleation: event minus control is
-0.013 at lead 2 and -0.014 at lead 1, with both confidence intervals spanning
zero. Its magnitude does increase (+0.075 and +0.089), indicating stronger
collective rotation without a preferred rotational sign, but it is less
discriminating than `(psi, Lambda)` at lead 1.

Aspect ratio is not an early precursor. It changes only at lead 1 and its sign
means that imminent schools are more isotropic by this conventional eigenvalue
ratio. It does not identify specifically wall-normal protrusion; `Lambda` does.

## Deterministic dynamics

The original microscopic equations were replayed for one event/control pair
per seed at leads 2 and 1 (80 snapshots). The replay separates wall avoidance,
alignment, attraction, direct hydrodynamic torque, BEM-mediated torque, and the
rotation of the COM-to-opening direction. The angular residual divided by
`In*sqrt(dt)` has mean 0.0019 and mean within-snapshot SD 0.996; no replayed
fish activated an override.

There is **no resolved collective restoring torque that holds the school in a
wall-following orientation before nucleation**:

- At lead 2, the median total angular contribution toward the opening is
  +0.121 for events and +0.171 for controls. Their paired difference is not
  resolved (mean +0.085; interval -0.310--0.494).
- At lead 1, the corresponding medians are +0.153 and +0.238; the paired mean
  is dominated by rare direct-hydrodynamic spikes and its interval is extremely
  wide. Robust medians show no event-specific angular drive.
- Wall torque itself generally rotates toward the opening, but is **weaker** in
  events than controls at lead 2 (paired difference -0.087; interval
  -0.131-- -0.047) and indistinguishable at lead 1. It therefore cannot explain
  why the matched event configuration nucleates.
- The total angular dynamics at lead 2 move event schools away from the nearest
  wall-following tangent relative to controls (paired difference -0.480;
  interval -0.955-- -0.043). This contradicts a stable wall-following angular
  well that must be crossed.
- `Lambda` has a resolved positive deterministic growth rate in events relative
  to controls at lead 1: +0.0757 per time unit (interval 0.0388--0.1154).
  Instantaneously this growth is principally self-propelled kinematics; direct
  hydrodynamic translation is a smaller modifier and BEM translation is small.

Direct hydrodynamic angular rates are intermittent. Two representative
snapshots contain very large values, so the report uses medians and retains
means, ranges, and every raw replay value in the tables. These spikes are real
outputs of the simulator equations but do not establish a reproducible
seed-level initiating torque.

## Markov-memory test

At a 0.1-time-unit lag, a linear state-only model of `(psi, Lambda)` has almost
no held-out increment skill. Adding oscillation phase reduces normalized RMSE
for the next `psi` increment by 6.3%; adding the preceding increment reduces it
by 9.1%, and phase plus history by 11.5%. Thus phase matters for angular motion,
but does not remove angular memory.

The failure is decisive for `Lambda`: adding its preceding increment reduces
normalized RMSE from 0.993 to 0.242 (75.6%), while phase alone does nothing.
The same strong history benefit remains when separate models are fit in each
phase quadrant (normalized RMSE 0.22--0.29 with history versus 0.89--0.98
without it). `Lambda` is a geometrical position variable with ballistic
extension/contraction memory, so an overdamped two-coordinate Markov model
cannot represent its dynamics.

## Scientific decision

1. **Do `(psi, Lambda)` predict nucleation?** Yes, genuinely but moderately.
   `psi` carries early information; `Lambda` adds useful information nearer A.
   The result is prospective, phase-matched, and held out by seed.
2. **Does a deterministic angular restoring mechanism exist?** Not at the
   collective-coordinate level. No event-specific restoring drift toward
   wall-following motion is resolved. Wall geometry steers individual leading
   fish in the earlier analysis, but the collective mean does not behave as a
   particle trapped behind a fixed angular barrier.
3. **Is a 2-D potential defensible?** No. The strong `Lambda` history, residual
   angular memory, phase dependence, and intermittent nonconservative
   hydrodynamic torques rule out an autonomous overdamped potential in
   `(psi,Lambda)`. A periodically driven escape model is closer, but it must be
   phase- and rate-dependent (for example `(psi, Lambda, dLambda, phase)`), not
   a static two-dimensional quasipotential inferred from occupancy.
4. **What is the better first-principles reduction?** The evidence favors
   angular first passage of a small leading-edge subpopulation, gated by the
   oscillation/wall geometry, followed by wall-normal recruitment of the rest
   of the school. Collective `(psi,Lambda)` are useful coarse observables of
   that process, not its closed causal state.

## Scope and limitations

- Event A is the first fish crossing; 340 analyzed episodes later committed and
  116 resolved as failures. The classifier predicts imminent A, not eventual
  whole-school commitment.
- Controls remain entirely in one chamber from 2.2 units before through 3 units
  after their pseudo-anchor. They are matched within seed on origin chamber and
  vertical phase only; no tested coordinate is matched away. Median phase-chord
  mismatch is 0.0016 and the 99th percentile is 0.0109.
- The 80 replay snapshots diagnose instantaneous deterministic terms. They are
  not enough to estimate a smooth drift field, and no potential is fitted.
- All uncertainty treats independent seeds, rather than frames or fish, as the
  replication unit.

The normalized angular-momentum convention follows Eq. 2 of Huang, Ling, and
Kanso (2024), which builds on the hydrodynamic fish model of Filella et al.:
[PNAS paper](https://doi.org/10.1073/pnas.2406293121) and
[Filella et al. PRL](https://doi.org/10.1103/PhysRevLett.120.198101).
