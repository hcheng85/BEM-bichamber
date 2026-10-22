# What turns the leading edge toward the gap?

## Conclusion

The trajectories do **not** support a static deterministic angular barrier
created by wall avoidance. In successful episodes, the mean wall contribution
toward the collision-free opening directions is -0.134 at `A-2`, +0.491 at
`A-1`, and +0.081 at `A-0.2` (rad/time). A restoring barrier should remain
directed away from the opening; this torque changes sign with approach phase.

The directly observed mechanism is a **phase-dependent angular gate**. The
leading ten fish progressively enter the cone bounded by rays tangent to the
two rounded tips: 25%, 31%, and 67% are admissibly oriented at `A-2`, `A-1`,
and `A-0.2`. Matched prospective controls remain at 5%, 7%, and 3%; resolved
failures reach 14%, 20%, and 34%. At `A-0.2`, success minus control is +0.64
(seed-bootstrap 95% interval 0.43 to 0.85), and success minus failure is +0.33
(0.08 to 0.61). Every measurement precedes the first crossing.

The replay identifies **wall geometry as the clearest deterministic turning
term at the intermediate approach phase**, not as a barrier. At `A-1`, wall
avoidance contributes +0.491 in successful episodes (interval 0.121 to 0.933).
Direct fish--fish hydrodynamics adds +0.196 on average but is highly variable
(-2.42 to 2.43), while alignment opposes the turn (-0.340). The total is +0.455
but imprecise (-2.11 to 3.12). Failures receive a similar wall turn (+0.546),
so wall steering explains why the front turns during this oscillation phase,
but not why every turn commits to a crossing. Relative to controls, the total
success-minus-control turn is positive in four of five seeds (+0.982), but its
interval (-2.05 to 3.52) crosses zero.

## Physical decomposition

Values below are seed-balanced means for successful nucleations. Positive means
rotation toward the shortest geometrically admissible crossing direction. Once
a fish is inside the cone, its angular error is zero and it is omitted from the
signed closing-rate average; admissible fraction is reported separately.

| contribution (rad/time) | A-2 | A-1 | A-0.2 |
|---|---:|---:|---:|
| wall avoidance | -0.134 | +0.491 | +0.081 |
| social alignment | -0.075 | -0.340 | +0.053 |
| social attraction | +0.179 | +0.100 | -0.015 |
| direct hydrodynamics | -0.419 | +0.196 | -0.441 |
| boundary-mediated hydrodynamics | +0.006 | +0.007 | +0.092 |
| total deterministic | -0.444 | +0.455 | -0.230 |

- **Wall avoidance is not a fixed barrier.** It weakly opposes the opening
  earlier, turns the front toward it at the intermediate phase, and is small
  just before A. This follows from the forward wall-intersection geometry
  changing as the school sweeps past a rounded tip.
- **Direct hydrodynamic torque is intermittent and seed-variable.** It can
  amplify individual turns but has no sign-stable seed-level mean here, so it
  is not established as the initiating cause.
- **Boundary-mediated hydrodynamic torque is small.** The BEM wall-reflection
  field is not the dominant angular steering term in this sample.
- **Social alignment does not initiate the turn**; it opposes it at `A-1`.
  Social attraction is weakly positive early (`A-2`: +0.179, interval 0.012 to
  0.409) and may help draw the front toward the tip region.
- **Noise supplies exploration, not a directed mean torque.** The inferred
  residual matches the programmed `N(0,1)` increment (mean -0.0055, SD 1.0004).
  Its one-step angular scale is `In*sqrt(dt)=0.0158` rad, larger than the mean
  deterministic increment `0.00046` rad at `A-1`; a phase-persistent bias can
  nevertheless accumulate across many steps.

At `A-0.2`, most successful leading fish are already inside the tangent cone.
The reported torque there describes only the remaining outside fish and must
not be read as the torque that produced the already-completed turn.

## Barrier test

A genuine angular barrier requires reproducible restoring drift away from the
opening over the relevant phase range, followed by noise-driven escape against
that drift. The first requirement fails: wall torque changes sign and becomes
significantly positive at `A-1`. The system is not purely free angular diffusion
either, because this phase-specific deterministic wall turn is resolved.

The defensible classification is:

> **noise-assisted, phase-dependent deterministic steering**, with wall
> geometry as the resolved intermediate-phase steering term and direct
> hydrodynamics as an unresolved, intermittent modifier.

This is evidence for a time/phase-dependent angular gate, not evidence for a
conservative angular energy landscape.

## Controls and replay validation

The earlier event set contained 353 successes and 159 failures. That imbalance
is not used here. The new design uses exactly two successes, two resolved
failures, and two prospective controls in each of five independent held-out
seeds. Controls are matched before the anchor on chamber, vertical center and
center velocity (oscillation phase), front distance, and school extent, but
never on heading or the angular gate. End-censored failures are excluded.

The replay follows the simulator update: exact rounded-tip and smoothed wall
intersections, Delaunay visibility, separate alignment and attraction, direct
hydrodynamics, the saved BEM solve, clipping, and jump overrides. Across 9,000
fish updates, the residual divided by `In*sqrt(dt)` has mean -0.0055 and SD
1.0004; no analyzed update activated an override. This is the expected
signature of the programmed angular noise after deterministic terms are
removed.

## Recommended theory and next test

Do **not** fit a one-dimensional angular potential from these results. A static
landscape would average phases at which wall torque changes sign. The justified
first-principles target is a phase-dependent angular escape equation,

```
d psi = b(psi, phase, coarse-front-geometry) dt + In dW,
```

with `b` assembled from the known wall, social, direct-hydrodynamic, and BEM
equations rather than inferred from crossing flux. Before solving it, densely
replay these same saved trajectories between `A-2` and A to measure torque
persistence and integrated angular impulse. No new simulation is required.
That test will determine whether intermittent direct hydrodynamics contributes
a reproducible impulse or whether successful entry is mainly noisy selection
during the wall-steered phase.

## Files

- `figures/01_pre_event_angular_dynamics.*`: angular distance, tangent-cone
  admissibility, and total deterministic turning.
- `figures/02_torque_decomposition_and_restoring.*`: physical decomposition and
  wall restoring behavior.
- `tables/leading_fish_torques.csv`: fish-level replay results.
- `tables/seed_level_summary.csv`, `group_summary.csv`, and
  `paired_contrasts.csv`: seed-balanced inference.
- `tables/selected_triplets.csv`: prospective matching audit.
- `tables/replay_validation.csv` and `source_audit.csv`: implementation audit.
- `scripts/`: complete reproducible analysis; raw trajectories are not copied.
