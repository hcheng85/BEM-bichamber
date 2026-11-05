# COM-relative turning balance in the circular schooling--milling transition

## Main conclusion

The exact COM-relative identity is numerically valid, but the data do **not**
support a simple microscopic locking mechanism in which each fish's angular
turning rate matches its geometrical orbital rotation.

Instead, the equations and stationary trajectories support a source-specific
balance:

- social alignment strongly reinforces whichever COM-relative circulation is
  already established;
- direct fish--fish hydrodynamic turning and angular-noise decorrelation oppose
  established milling circulation;
- attraction and wall turning provide small positive reinforcement in milling;
- boundary-mediated hydrodynamic turning is negligible in this projection;
- the projected COM-geometrical rotation term supports weak signed tangency in
  schooling but averages to approximately zero in sustained milling.

In seed 9, the old positive circulation begins to reverse before polarization
collapses.  Direct hydrodynamics is negative during this interval, but its mean
lies inside the range of stable-schooling controls.  The unusual observation is
the loss of the normally positive geometrical contribution.  The total drift is
still inside the control range.  Therefore the data show a change in balance
that permits chirality reversal, not a unique event-specific torque impulse.
Once the new negative circulation appears, alignment strongly amplifies it.

This result identifies how milling circulation is maintained and how the one
observed reversal unfolds.  It does not establish a causal transition law or a
predictive stability boundary.

## Exact construction

For each stored state, the analysis uses the instantaneous center of mass and
mean translational velocity,

```text
C     = mean_i(x_i)
vbar  = mean_i(v_i)
y_i   = x_i - C
alpha_i^C = theta_i - atan2(y_iy, y_ix)
omega_i^C = cross(y_i, v_i-vbar) / |y_i|^2.
```

With the source angular dynamics

```text
d theta_i = Omega_i dt + In dW_i,
```

the exact relative dynamics and Itô drift are

```text
d alpha_i^C = (Omega_i-omega_i^C)dt + In dW_i,

b_MC = mean_i[cos(alpha_i^C)(Omega_i-omega_i^C)]
       - In^2 M_C/2,
M_C  = mean_i[sin(alpha_i^C)].
```

The angular part is decomposed without fitted coefficients into alignment,
attraction, wall avoidance, direct hydrodynamics, BEM-mediated hydrodynamics,
and the residual of the production angular clip.  The geometrical term
`-mean(cos(alpha_i^C) omega_i^C)` remains separate.

## Immutable inputs and bounded computation

- Ten original N=100 full-BEM trajectories, seeds 0--9.
- `R=7.158910531638177`, `dt=0.001`, `In=0.5`, `Ip=4`, `Iw=0.94`,
  `If=0.01`.
- Production source `BEM_circular_parallel_instance.py`, SHA-256
  `2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`.
- Velocity convention: `rdot[:,:,step+1]` advances `state[:,:,step]`.

The large MATLAB-v5 arrays were opened with read-only memory maps.  A compact
COM series was sampled every 0.1 time unit.  Source-exact replay used 1,450
stationary snapshots, 301 seed-9 transition snapshots, and 1,155 snapshots in
21 prespecified stable-schooling control windows.  No complete trajectory was
loaded into memory.  The simulator, interaction laws, safety rules, and MAT
files were not modified.

## Replay and stochastic validation

The replay uses the original Delaunay visible-neighbor graph, visual weighting,
wall ray intersection, direct hydrodynamics, BEM solve, collision/extreme
fallbacks, positional retry, and total angular clipping.

Across 2,906 replayed states (290,600 fish updates):

| audit | result |
|---|---:|
| maximum stored/replayed velocity error | `5.33e-15` |
| maximum split/direct social torque error | `2.30e-12` |
| angular clips | 0 |
| collision overrides | 0 |
| extreme-velocity overrides | 0 |
| sticky positional retries | 0 |
| fish excluded at `|y_i| <= 1e-8 R` | 0 |
| minimum sampled `|y_i|/R` | `4.02e-4` |

Six fish-snapshots lie within `1e-3 R` and 512 within `1e-2 R`; they are
reported rather than silently discarded.  None reaches the singular exclusion
threshold.

The finite-step stochastic audit uses the actual deterministic next position
and

```text
E[sin(alpha_i,next^C) | state]
 = exp(-In^2 dt/2) sin(theta_i + Omega_i dt - phi_i,next^C).
```

After normalization by its exact conditional variance, the observed-minus-
expected martingale increment has mean `7.46e-5`, variance `0.997`, skewness
`-0.0157`, and excess kurtosis `-0.130`.  Thus angular noise is correctly
separated from deterministic turning.

Ten long contiguous time blocks quantify temporal uncertainty.  The continuous
generator's mean contains zero in 10 of 11 stationary intervals.  The exception
is schooling seed 0, where near-COM geometry and finite-step curvature make the
continuous estimate less robust.  The source-exact finite-step expected drift
contains the observed endpoint rate in all 11 intervals.  Stationarity is
therefore consistent with the microscopic update when finite-step and noise
effects are treated correctly.

## Stationary physical balance

Because circulation has either sign, contributions below are multiplied by the
interval's persistent chirality.  Positive values reinforce the existing
circulation and negative values oppose it.  Uncertainty is across independent
seed/interval units: eight schooling intervals and three milling intervals.

| contribution | schooling mean [95% interval] | milling mean [95% interval] |
|---|---:|---:|
| alignment | `+0.0432 [0.0254, 0.0609]` | `+0.3017 [0.2880, 0.3154]` |
| attraction | `+0.00052 [-0.00074, 0.00178]` | `+0.00576 [0.00494, 0.00659]` |
| wall | `-0.00022 [-0.00115, 0.00071]` | `+0.00248 [0.00128, 0.00368]` |
| direct hydrodynamics | `-0.0130 [-0.0413, 0.0153]` | `-0.1917 [-0.3254, -0.0579]` |
| BEM hydrodynamics | `1.5e-6 [-1.0e-5, 1.3e-5]` | `-3.1e-7 [-1.4e-6, 0.8e-6]` |
| COM geometry | `+0.0191 [0.0115, 0.0267]` | `0.00004 [-0.0108, 0.0109]` |
| Itô noise decay | `-0.0177 [-0.0247, -0.0107]` | `-0.1109 [-0.1128, -0.1091]` |

Mean `|M_C|` is `0.198 [0.189,0.207]` in schooling and
`0.887 [0.875,0.899]` in milling.  The milling balance is consequently clean:
alignment is the principal circulation-producing angular term, while direct
hydrodynamics and noise are the principal opposing terms.  Their near balance
maintains a stationary high-`|M_C|` state.  Attraction and wall avoidance are
small but consistently reinforcing; BEM angular turning is negligible.

This is not microscopic rate locking.  Across stationary intervals, the mean
per-fish RMS values are approximately

| state | `RMS(Omega)` | `RMS(omega_C)` | `RMS(Omega-omega_C)` | correlation |
|---|---:|---:|---:|---:|
| schooling | 7.35 | 0.755 | 7.39 | 0.0012 |
| milling | 7.76 | 1.02 | 7.75 | 0.0148 |

Instantaneous angular and orbital rates therefore do not track each other fish
by fish.  The small milling geometry contribution also has a geometric reason:
when headings are nearly tangential, `cos(alpha_i^C)` is small, so orbital
rotation contributes little to the drift of the already saturated sine order
parameter.

## Seed-9 transition

The persistent schooling-to-milling boundary is at approximately `t=312.41`.
The source-exact window means show:

| window | `M_C` | `P` | alignment | direct hydro | geometry | total `b_MC` |
|---|---:|---:|---:|---:|---:|---:|
| 300--306 baseline | +0.205 | 0.914 | +0.0269 | -0.0203 | +0.0539 | +0.0355 |
| 306--311.5 circulation loss | +0.144 | 0.918 | +0.0489 | -0.1007 | +0.0021 | -0.0663 |
| 311.5--315 polarization collapse | -0.578 | 0.591 | -0.3252 | +0.1002 | +0.0121 | -0.2000 |
| 315--325 early milling | -0.809 | 0.357 | -0.2930 | +0.1141 | -0.0007 | -0.0619 |

The sign reversal of `M_C` begins before the large polarization loss.  During
the pre-collapse interval, direct hydrodynamic turning pushes against the old
positive chirality and the positive geometrical contribution nearly vanishes.
During and after the sign reversal, alignment changes sign and strongly drives
the new negative circulation; direct hydrodynamics and noise then oppose it.
Alignment therefore consolidates the new milling state rather than providing a
demonstrated initiating impulse.

To avoid interpreting one event in isolation, the 306--311.5 interval was
compared with 21 equally long windows from the seven persistent-schooling
seeds.  The chirality-aligned direct-hydrodynamic mean is `-0.1007`, compared
with a control mean `-0.0514` and control-window 2.5--97.5% range
`[-0.1449,0.0223]`: it is not exceptional.  The event's geometrical mean
`+0.0021` is below the control range `[0.0062,0.0560]`.  The total drift
`-0.0663` remains within the control range `[-0.0733,0.0724]`.

The defensible interpretation is therefore limited: loss of geometrical
stabilization is associated with the observed reversal, and ordinary negative
direct-hydrodynamic fluctuations then outweigh alignment.  Because the total
drift is not outside control variability and only one transition exists, this
does not establish causation or a universal trigger.

## Answers to the scientific questions

### 1. What mechanism do the equations and data support?

They support a **circulation-maintenance and reversal mechanism**, not a static
potential or per-fish angular/orbital locking.  Alignment amplifies an existing
COM-relative chirality.  Direct hydrodynamic turning and angular noise oppose
that ordered circulation.  In the observed transition, the normally positive
schooling geometry contribution disappears before polarization collapses;
alignment subsequently amplifies the opposite chirality and establishes
milling.

### 2. Which microscopic terms drive or oppose COM-relative circulation?

- Sustained milling: alignment is the dominant driver; attraction and wall
  avoidance are small reinforcing terms.
- Direct fish--fish hydrodynamic turning and Itô angular-noise decay are the
  dominant opposing terms.
- BEM-mediated angular turning is negligible in this projection.
- COM-relative geometry reinforces weak schooling chirality on average but is
  neutral in stationary milling.
- No clip or numerical override contributes in the sampled states.

### 3. What remains unresolved before a predictive stability boundary?

`M_C` has not been shown to be a closed Markov variable, and one spontaneous
transition cannot identify a response law.  A predictive boundary requires
multiple independent transitions or controlled perturbations about both
stationary branches, with enough data to measure how the source-exact
alignment/direct-hydrodynamic/geometry balance changes under perturbations and
across parameters.  In particular, one must determine whether loss of the
geometrical contribution reproducibly precedes switching and measure the
growth/decay rate of circulation perturbations.  This dataset establishes the
balance at one parameter point, not its stability derivative.

## Reproducibility map

- `analyze_com_relative_turning.py`: complete streaming and source-exact replay.
- `outputs/analysis_config.json`: sampling and uncertainty configuration.
- `outputs/replay_validation.json`: source hash, indexing, override, singularity,
  velocity, social-torque, and martingale audits.
- `outputs/stationary_balance_by_seed.csv`: interval-level decomposition and
  block uncertainty.
- `outputs/stationary_group_comparison.csv`: seed-level stationary comparison.
- `outputs/seed9_transition_com_drift.csv` and
  `outputs/seed9_transition_window_means.csv`: transition timeline and bounded
  window summaries.
- `outputs/schooling_control_window_means.csv` and
  `outputs/precollapse_control_comparison.csv`: prospective control comparison.
- `outputs/com_relative_series.npz`: compact 0.1-time-unit COM diagnostics; no
  raw trajectories or per-fish histories.

