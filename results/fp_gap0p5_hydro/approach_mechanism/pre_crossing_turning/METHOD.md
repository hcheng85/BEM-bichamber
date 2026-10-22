# Method

## Sampling design

- Held-out independent seeds: 15--19.
- Two balanced triplets per seed: one successful nucleation, one resolved
  failed nucleation, and one prospective non-crossing anchor.  Total: 10 per
  condition.
- Within each seed, success and failure share the origin chamber and are paired
  using the state at `A-2`. Controls are required to retain all 100 fish in the
  same chamber from 2.2 time units before to 4 time units after their anchor.
- Matching variables are chamber, vertical center, vertical center velocity
  (oscillation phase proxy), front distance, and horizontal/vertical extent.
  Heading and angular admissibility are deliberately excluded.
- Snapshots are at `A-2`, `A-1`, and `A-0.2`, before the first fish crossing.
  End-censored episodes are excluded.

The leading edge is the ten fish with smallest distance to the divider. This
keeps the description coarse and fixed before looking at outcome.

## Angular geometry

After reflecting right-origin episodes, local inward direction is angle zero.
For a fish at distance `d` and height `y`, the limiting rays are the inward
tangents to the lower and upper rounded-tip circles. If the tip radius is `r`
and their centers are at `y_c = +/- (gap/2+r)`, the lower and upper limits are

```
lower = atan2(y_c_lower-y,d) + asin(r/distance_to_lower_center)
upper = atan2(y_c_upper-y,d) - asin(r/distance_to_upper_center).
```

These are the exact collision-free angular boundaries imposed by the rounded
tips for a point swimmer. Angular error is the shortest periodic rotation to this interval, so it
is never greater than pi. A contribution is positive when it reduces that
error. Fish already inside the interval have zero error and are omitted from
the signed closing-rate average; their fraction is reported separately.

## Equation replay

`scripts/replay_dynamics.py` reproduces:

1. Delaunay visual neighbors with visibility blocked by the divider and rounded
   tips;
2. separate visual alignment and attraction terms;
3. direct dipolar hydrodynamic velocity and torque;
4. the boundary integral solve using the saved geometry and LU factors;
5. wall-avoidance torque from the simulator's forward-ray intersection with the
   outer circle, straight divider segments, rounded tips, and smoothing arcs;
6. the simulator's clipping and hydrodynamic-removal rule for wall jumps.

The angular residual is

```
delta_theta - dt*(wall + alignment + attraction + direct_hydro + boundary_hydro).
```

Dividing by `In*sqrt(dt)` gives pooled mean -0.0055 and standard deviation
1.0004 over 9,000 replayed fish updates. This matches the programmed standard
normal angular increment. No replayed update activated the jump/large-velocity
override.

Uncertainty uses the five independent seeds as the resampling unit. The 100
fish or two events within a seed are not treated as independent replicates.
