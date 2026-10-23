# What turns a leading fish toward the opening?

## Scientific conclusion

The trajectories support **deterministic angular steering into the admissible
opening corridor, dominated by social alignment and assisted by net
hydrodynamic torque**, rather than escape over a wall-induced angular barrier.
Angular noise contributes to individual entries, but the data do not show that
noise is necessary or that its cumulative increment distinguishes imminent
crossings from matched controls.

This conclusion is deliberately qualified. The deterministic entry mechanism
is generic: fish in no-crossing controls that enter the angular corridor have
essentially the same torque budget. Imminent-event configurations are already
better aligned at A-2. The present window therefore does **not** identify the
earlier cause that selected those configurations.

Of the three prespecified possibilities, the evidence supports **A,
deterministic steering**, for angular-gate entry. It does not support B as a
necessary-noise claim. It also rejects the strong form of C because the mean
deterministic closing impulse is much larger than the mean stochastic closing
increment among entrants. However, A is not a complete event-selection law.

## Data and prospective design

The analysis uses all 456 first-crossing events and their 456 existing controls
from 20 independent full-BEM seeds. Controls were selected prospectively within
seed and origin chamber, matched on center distance and vertical oscillation
phase at A-2, and had no first crossing during the following two time units.
Heading was not matched.

At A-2, the ten fish closest to the divider were selected and their identities
were held fixed. No fish was added retrospectively. For each fish, the angular
error `eta` is measured relative to the nearest collision-free tangent to the
two rounded divider tips. Angular entry and spatial crossing are distinct
outcomes.

## Fixed identities are an important negative result

- The eventual first crosser was in the fixed ten in only **235/456 = 51.5%**
  of events. The mean across seeds was 0.514, with seed-bootstrap 95% CI
  **[0.467, 0.560]**.
- By A, only **36.3%** of the original ten remained in the instantaneous leading
  ten for events (controls: 37.7%).
- Thus the leading edge is not a stable material subpopulation. A model based on
  a permanently designated leader is not supported.

Figure 1 shows the fixed-identity time histories and the initial rank of every
captured first crosser. Imminent-event fish start closer to the angular gate and
remain more likely to occupy it, even though the identities turn over rapidly.

## Prospective angular differences

Across the fixed ten, with each seed treated as one independent replicate:

| quantity at or after A-2 | event | control | event - control, 95% seed CI |
|---|---:|---:|---:|
| initial mean absolute angular error | 0.554 rad | 0.633 rad | **-0.078 [-0.107, -0.052]** |
| entered angular gate by A | 0.818 | 0.748 | **+0.069 [0.056, 0.085]** |
| crossed spatially by A | 0.052 | 0 | **+0.052 [0.047, 0.057]** |

Among captured first crossers, comparison with the fish of the same A-2 rank in
the matched control gives an even larger initial angular advantage:
`Delta |eta| = -0.292 rad [-0.360, -0.225]`. But neither the total deterministic
closing impulse nor the stochastic closing increment over A-2 to A differs
reliably from that same-rank control: their intervals both cross zero. Much of
the separation therefore predates A-2.

## Exact angular budget and source decomposition

In the chamber-reflected coordinates,

```
d eta_i = (Omega_wall + Omega_alignment + Omega_attraction
           + Omega_direct + Omega_BEM - d alpha_i/dt) dt
          + In dW_i .
```

The simulator seeds NumPy once per trajectory, consumes exactly 300 uniform
draws in initialization, and then draws 100 Gaussian angular increments per
native step. Replaying that stream makes the observed-minus-noise deterministic
increment exact, including clipping and collision overrides. At valid audit
snapshots:

- maximum initialization error: `2.7e-15`;
- maximum noise-increment error: `1.25e-14`;
- maximum deterministic-update error: `1.25e-11`;
- boundary-torque RMSE on unclipped, non-overridden fish: `2.94e-12`;
- selected-fish wall/social/direct source terms agree with full replay to
  floating-point precision.

For fish initially outside the gate that later enter it, seed-cluster means are:

| closing contribution up to first entry | event | control |
|---|---:|---:|
| total deterministic torque | **0.338 rad** | **0.354 rad** |
| stochastic angular increment | 0.022 rad | 0.025 rad |
| motion of tangent corridor | -0.020 rad | -0.022 rad |

The event-control differences for all three entry-budget terms include zero.
The deterministic entry mechanism is therefore reproducible but not specific to
successful nucleation.

The deterministic source budget up to entry is:

| source | event | control | interpretation |
|---|---:|---:|---|
| wall avoidance | 0.020 | 0.038 | small; no angular wall barrier demonstrated |
| social alignment | **0.224** | **0.211** | largest stable closing contribution |
| social attraction | -0.046 | -0.047 | opposes angular entry on average |
| direct hydrodynamics | 0.107 | 0.243 | closing, but variable |
| BEM/override residual | 0.034 | -0.091 | uncertain and cancels direct hydrodynamics |

Direct and boundary-mediated hydrodynamic angular terms can become individually
large and opposite. Their sum and the exact total deterministic impulse are
better conditioned than either term alone. Figure 2 displays this cancellation
separately so that it does not hide the smaller, reproducible alignment term.

## Archived RNG discontinuity

Six files (seeds 10, 11, 13, 15, 17, and 19) have an RNG-stream discontinuity
at native step 600,000, consistent with a resumed run whose RNG state was not
preserved. Geometry and headings remain valid. The deterministic/noise split is
flagged unavailable for 149 of 912 segments crossing that boundary; all 912
segments remain in geometry, identity-turnover, and gate-entry analyses. Torque
contrasts use the remaining 763 valid segments and still cluster uncertainty by
all 20 seeds. No missing stochastic increments were imputed.

## First-principles model justified by the result

The simplest defensible reduction is not a static effective potential. It is a
time-dependent angular first-passage model for a **renewing leading-edge
subpopulation**:

```
d eta_i = b_i(X(t), phase(t)) dt + In dW_i - d alpha_i,
```

where `b_i` is retained from the microscopic wall, alignment, attraction,
direct-hydrodynamic, and BEM torque laws; `alpha_i` is the rounded-tip tangent
direction; and entry into the collision-free angular interval is absorbing for
the angular stage. A separate self-propelled spatial-approach stage follows.
The active set should be proximity-defined and allowed to turn over rather than
assigning a permanent leader.

No angular potential is justified: no restoring wall-following drift or
event-specific one-dimensional drift was demonstrated. A future bounded test,
if desired, should move the observation start earlier for the actual prospective
risk set and retain a rolling proximity-defined active set; it should not fit a
potential to occupancy or crossing flux.

## Direct evidence versus interpretation

- **Exact consequence of simulator equations:** the angular budget identity,
  tangent-gate geometry, and source replay on valid RNG intervals.
- **Directly supported by trajectories:** initial angular advantage, higher
  gate-entry probability, rapid leading-set turnover, dominance of deterministic
  over stochastic closing among entrants, and alignment as the largest stable
  source term.
- **Interpretation:** deterministic alignment/hydrodynamic steering supplies a
  generic angular opportunity; proximity and the renewing leading set determine
  which opportunity becomes a first crossing.
- **Still unresolved:** what created the event group's angular advantage before
  A-2, and whether noise is counterfactually necessary for particular entries.

