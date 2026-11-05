# Next numerical validation: stationary drift and transition mechanism

## Why this test comes before another phase-diagram sweep

The ten circular full-BEM simulations (`Ip=4`, `In=0.5`, `Iw=0.94`,
`If=0.01`, `N=100`, `R=7.158910531638177`) reveal two long-lived
states. Seeds 7 and 8 remain milling; seed 9 switches from schooling to
milling at approximately `t=312.41`. The current torque replay uses 170
selected snapshots (only 30 from milling) and finds a mean instantaneous
projected total drift `bQ_total = +0.03257` in milling. In contrast,
the full 980-time-unit post-formation time series for seeds 7 and 8 has
nearly zero net drift in `Q`.

For independent angular Brownian noise with amplitude `In`, the exact
Itô identity for `Q = |N^{-1} sum_i (cos(theta_i),sin(theta_i))|^2` is

\[
\frac{d}{dt}\mathbb E Q
=\mathbb E\left[\frac{2}{N}\sum_i
    (\mathbf p\cdot\mathbf e_{\perp i})\Omega_i\right]
-I_n^2\left(\mathbb E Q-\frac1N\right).
\]

For a *stationary, unconditioned* trajectory the mean deterministic
source must equal `In**2 * (mean(Q) - 1/N)`. The observed full-series
means imply approximately `0.00804` for seed 7 and `0.01074` for
seed 8, versus about `0.02389` and `0.06503`, respectively, in the
present ten-snapshot-per-seed replay. A third milling contribution
comes from seed 9 after its state change. The discrepancy is a
**sampling/source-audit question**, not evidence that the milling
branch is unstable. The schooling snapshot mean `0.1980` is close to
the full-series stationary target of about `0.2022`.

The main scientific aim is to establish whether the exact microscopic
torque balance is internally consistent and, only then, determine
which changes in fish configuration affect polarization growth near
the observed state transition. Do **not** fit a reduced potential or
interpret a single transition as a measured transition rate.

## Phase 0: reproducible check using committed compact data

Run:

```bash
python theory/circle_collective_first_stage/check_stationarity_identity.py
```

This prints per-seed full-series values, corresponding stationary
torque-balance targets, and the small replay averages. Seeds are selected
as persistently schooling/milling using the pre-existing occupancy
table; seed 9 is listed separately, not treated as stationary over its
entire 1000-time-unit history. This diagnostic does **not** prove a
violation of Itô calculus: sparse sampling and finite-time martingale
fluctuations must be quantified in Phase 2.

## Phase 1: audit the original angular update before reinterpreting any drift

Use the original production source
`BEM_circular_parallel_instance.py` (source SHA-256 recorded in
`simulation_parameters.json`) and the *same* MAT files. The simulation
and its forward-ray wall interaction are immutable inputs.

1. Identify and document the precise update used for `theta`,
   including the deterministic angular torque, any angular clipping,
   noise variance and distribution, and operation order. Do **not**
   assume that `In` is the coefficient of `dW` until the code confirms
   this. If it differs, update the projected Itô convention in a
   separate analytical note, not the production simulator.
2. Use the original `GetVision`, `GetWallAvoid`, direct hydrodynamics,
   BEM hydrodynamics, and any actual angular update overrides, following
   `replay_torques.py`. First independently reproduce the total angular
   drift used in the source; agreement of *translational* velocity alone
   does not establish agreement of angular drift.
3. At stratified, reproducibly selected snapshots in steady milling
   seeds 7 and 8 (plus stationary intervals in seed 9), compute the
   observed one-step angular increment as the principal wrapped angle
   difference. Compare `delta_theta - omega_effective*dt` against the
   exact stochastic increment implied by the audited source.
   Record residual mean, variance divided by `dt`, distribution shape,
   time correlations, and any conditional dependence on wall distance.
4. Count numerical override/clip events and report them separately.
   Do **not** 'fix' the fallback boundary, substitute a different wall
   law, or treat a rare numerical safeguard as a physical mechanism.
   Discard or flag such snapshots when a decomposition is ambiguous.

**Acceptance:** the replay must match the angular update and the
original source's noise convention to a stated statistical tolerance.
If it fails, stop the physics interpretation and report the exact
disagreement; do not patch the original simulator silently.

## Phase 2: time-uniform stationary drift consistency

1. Prioritize complete post-formation milling intervals in seeds 7 and
   8. Do not select snapshots for high/low `Q` or near state transitions.
   Evaluate the full microscopic torque on a *uniform-in-time* or
   randomized-uniform grid, initially ~100-200 snapshots per seed,
   increasing until the mean is reasonably converged. Preserve exact
   source interaction laws (including Delaunay visibility, flow and BEM).
2. Recompute `bQ_alignment`, `bQ_attraction`, `bQ_wall`,
   `bQ_direct_hydro`, `bQ_bem_hydro`, clipping residual and the Itô
   angular-noise term from each evaluated state. Check per-snapshot
   algebraic additivity. Estimate time averages for each seed
   separately before combining them.
3. Compare the time-averaged deterministic sum against
   `sigma_theta**2 * (mean(Q) - 1/N)` for the *source-verified* noise
   amplitude `sigma_theta`. Also compare the complete drift average
   against `(Q(t1)-Q(t0))/(t1-t0)`, allowing for the stochastic
   martingale and finite-step bias. Use trajectory-block uncertainty,
   convergence across increasing snapshot counts and a time-window
   sensitivity analysis; a three-seed bootstrap alone is inadequate.
4. Repeat the same check for the continuously schooling seeds as a
   positive control. Keep seed 9's pre/post transition intervals
   separate, because its entire trajectory is not stationary.

**Acceptance:** all compared quantities use the same time interval,
noise convention and physical torque. Stationarity identities must be
consistent within a reported uncertainty; otherwise do not use these
projected drifts to claim a stability mechanism.

## Phase 3: focused seed-9 transition analysis (exploratory)

After Phases 1-2 pass, replay around `t≈312.41` with higher cadence
(start with `t=300..325`, every `0.1` time unit; adjust to measured
variation and computation budget). Retain stable reference windows
from schooling and milling seeds. Examine the chronological changes of

- polarization `P`, squared polarization `Q`, original COM- and
  velocity-based milling `Pi`, signed circle-centered heading milling
  `M`, COM radius, radial spread, and mean fish radius;
- each projected torque contribution to `bQ`, *and* their sum;
- neighbor visibility/heading correlations and distance-to-wall
  statistics, without modifying the interaction;
- absolute/relative times of the circulation change, polarization
  collapse, spatial redistribution and any torque-balance change.

Evaluate whether a prospective precursor is distinct from an ordinary
fluctuation in control schooling windows. One switching event **cannot**
establish general temporal ordering, rates, or causality.

## Output contract

Create a new additive subdirectory
`theory/circle_collective_first_stage/outputs/drift_consistency/` with:

- `angular_update_audit.json`: exact implemented noise/update law,
  source hash, sampled step counts, residual statistics, and override
  counts;
- `stationary_balance_by_seed.csv`: seed, interval, number of replayed
  states, mean `Q`, predicted stationary deterministic source,
  replayed deterministic source and its components, total drift,
  uncertainty estimates and convergence diagnostics;
- `stationary_balance_convergence.png`: means and uncertainty versus
  increasing time-uniform snapshot counts for milling seeds 7 and 8;
- `transition_seed9_torque_timeline.csv` and matching plot, if Phase 3
  is completed;
- `REPORT.md`: whether each validation stage passed, what remains
  ambiguous, and which *specific* analytic closure is plausible.

Raw MAT files, source code and intermediate per-fish arrays are not to
be copied into Git. Keep this run read-only relative to production
simulation and archive. Use a new analysis branch and a normal
non-force push; preserve all existing reports and figures.

## The theoretical decision after validation

The exact polarization-vector generator is

\[
d\mathbf p=\left[\frac1N\sum_i\mathbf e_{\perp i}\Omega_i
-\frac{I_n^2}{2}\mathbf p\right]dt
+\frac{I_n}{N}\sum_i\mathbf e_{\perp i}\,dW_i.
\]

Once the microscopic drift is validated, test whether its conditional
response close to the milling ensemble admits
`E[sum(e_perp*Omega)/N | p] ≈ A p`, with a response matrix `A`
whose entries can ultimately be calculated from a justified
geometric/kinetic model. A candidate *polarization-mode* threshold is
`max(Re(eig(A))) = In**2/2`; it is **not yet** the full stochastic
schooling–milling phase boundary, because the collective variables
may retain hidden spatial memory and noise affects the ensemble.
