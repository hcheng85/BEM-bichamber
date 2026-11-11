# Microscopic balance and closure of the gap-20 approach current

## Question and scope

This bounded extension asks what creates and sustains the source-side approach
current $J_a$, whether the parameter-free plug-flow estimate predicts neck
transport, and whether $(q,J_a)$ is an autonomous Markov state.  It uses only
the existing 40 no-hydrodynamics trajectories (`N=100`, `In=0.5`, gap 20).
No simulation, potential, phenomenological coefficient, or parameter scan was
introduced.

For crossing direction $s=+1$ from the left and $s=-1$ from the right, the
approach indicator and stored-velocity current are

\[
 I_i={\bf 1}\{s x_i<0,\;-s x_i\leq\ell,\;|y_i|\leq w/2\},\qquad
 J_a={s\over N}\sum_i I_i v_{x,i},
 \quad \ell=w=2.8635642126552705.
\]

The source file
`BEM_bichamber_parallel_instance_no_hydro.py` was read from the parent archive
and has SHA-256
`5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`.
It was not changed.  Because that archived file currently initializes gap 5,
the analysis executes it in memory after one audited assignment replacement,
`[50,5,20] -> [50,20,20]`, reproducing the geometry of the gap-20 raw files.

## Exact finite-step balance

Away from a collision override, $v_{x,i}=\cos\theta_i$.  Define the ideal
orientational current

\[
 \widetilde J_a={s\over N}\sum_i I_i\cos\theta_i.
\]

The production angular update is

\[
 \theta_i^*=\theta_i+\Omega_i^{\rm eff}\Delta t,
 \qquad
 \theta_i'=\theta_i^*+I_n\sqrt{\Delta t}\,\xi_i,
\]

where

\[
 \Omega_i^{\rm eff}=\Omega_i^{\rm wall}+\Omega_i^{\rm social}
                     +\Omega_i^{\rm clip}.
\]

The social term is the original combined visual alignment/attraction law; it
is not re-split by an invented convention.  For each deterministic angular
increment, the source allocation uses the exact line factor

\[
 L_i={\cos(\theta_i+\Omega_i^{\rm eff}\Delta t)-\cos\theta_i
            \over \Omega_i^{\rm eff}\Delta t},
\]

with limit $-\sin\theta_i$.  The native-step identity is then

\[
 \Delta J_a=B+W+S+C+R_n+\Delta M+\Delta J_{\rm collision}.
\]

Here $B=sN^{-1}\sum_i(I_i'-I_i)\cos\theta_i'$ is boundary
transport; $W,S,C$ are the line-integrated wall, social, and clipping
contributions; $R_n$ is the exact Gaussian relaxation

\[
 {s\over N}\sum_i I_i
 \left[e^{-I_n^2\Delta t/2}-1\right]\cos\theta_i^*;
\]

$\Delta M$ is the realized zero-conditional-mean angular-noise residual; and
the last term is the measured change in $J_a-\widetilde J_a$ caused by the
production collision velocity override.  This construction respects the
actual update order and stored convention `rdot[:,:,k+1]` advances
`state[:,:,k]`.

The source-exact replay used eight approximately evenly spaced, independent
event seeds and their same-$q$ controls over the two units before first
crossing.  The optimized visibility graph agreed with the production graph at
every audit snapshot.  Maximum native-step and integrated closure errors were
$9.4\times10^{-17}$ and $7.7\times10^{-17}$.  Reconstructed standardized
noise had mean -0.00035 and variance 0.99960.  There were no angular clipping
events in these intervals.

## What grows the current?

The table gives the event-minus-control integrated contribution over the two
pre-crossing units.  Intervals are seed-cluster bootstrap intervals over the
eight source-replayed pairs.

| contribution | paired difference | 95% interval |
|---|---:|---:|
| measured change in $J_a$ | 0.0499 | 0.0101, 0.0943 |
| approach-boundary transport | 0.0396 | 0.0049, 0.0806 |
| social turning | 0.0418 | 0.0169, 0.0642 |
| wall turning | -0.0362 | -0.0569, -0.0156 |
| angular-noise relaxation | -0.0130 | -0.0229, -0.0044 |
| realized stochastic increment | 0.0189 | -0.0135, 0.0598 |
| collision override | -0.00122 | -0.00367, 0.00000 |
| clipping | 0 | 0, 0 |

The mean measured $J_a$ change is 0.0350 in event intervals and -0.0149 in
controls.  Its reproducible positive production is explained by two comparable
effects: fish are recruited into the approach region with gap-directed
headings, and the social law turns fish toward the gap.  Wall turning opposes
this event-specific growth rather than causing it.  Gaussian angular diffusion
also relaxes directional current, as required by the exact Itô factor.  The
realized noise difference is uncertain and the collision correction is small.

Thus the trajectories support **socially maintained directional approach plus
spatial recruitment**, resisted by wall turning and angular decorrelation.
They do not support a wall-driven turn or a reproducible favorable noise kick.
Boundary transport is a membership term, so it must not be misread as a second
angular torque.

## Parameter-free plug-flow audit

The proposed estimate

\[
 F_{\rm neck}\simeq {N J_a\over\ell}
\]

was evaluated only on sustained episodes in held-out seeds 30--39.  There are
164 eligible episodes and 1,777 prospective samples.  Net boundary crossings
were counted over 0.2-unit windows.  The exact approach-region continuity
audit,

\[
 \dot N_a+F_{\rm neck}+F_{\rm upstream}+F_{\rm lateral}=0,
\]

closed identically for every window.

The plug expression is a useful qualitative instantaneous ordering but not a
quantitative flux law:

- In the same 0.2 window, pooled Spearman correlation is 0.494 and the median
  within-seed correlation is 0.498.  However, its mean prediction is 3.57
  fish/time versus an observed 6.60, with seed-mean bias -3.10
  (95% interval -3.54, -2.63) and seed-mean MAE 4.83.
- A flux-weighted physical travel time was calculated directly from fish
  distances and positive velocities; its median is 2.03.  Shifting the target
  neck window by that time reduces pooled correlation to 0.135 and median
  within-seed correlation to 0.042.  The seed-mean bias becomes +1.47
  (0.69, 2.31), and 68.3% of delayed windows contain zero net neck flux.
- Plug error correlates with approach accumulation (Spearman 0.350) and, more
  weakly, lateral loss (0.180), but not upstream exchange (0.031).

The failure after a Lagrangian travel delay is physically informative.  During
sustained collective propagation, $J_a$ is an Eulerian current distributed
across a continuously renewing approach population; it is not a tagged packet
whose mean members arrive together two units later.  Strong accumulation,
lateral exchange, nonuniform longitudinal density, and renewal of fish inside
the region violate plug flow.  The exact continuity equation succeeds, whereas
$N J_a/\ell$ does not supply the missing boundary state.

## Is $(q,J_a)$ Markovian?

The conditional drift audit uses one event/control pair from each of 37
event-producing seeds, sampled every 0.1.  Drift over the next 0.2 was compared
within fixed empirical $(q,J_a)$ cells.  This is a closure diagnostic, not a
fitted drift model.

- A coarse rising-versus-falling sign contrast is unresolved: -0.0049
  (-0.0230, 0.0113).
- After removing the $(q,J_a)$-cell means, pooling each matched event and
  control interval gives a mean within-pair Spearman correlation of 0.155
  (0.064, 0.244).  This is not robust evidence for memory in transition paths:
  event intervals alone give 0.048 (-0.046, 0.139).  The control-only estimate
  is 0.266 (0.011, 0.495), but only 14 controls have nondegenerate variation.
- Within the same cells, future-crossing intervals have future drift 0.0304
  larger than controls (0.0074, 0.0523).
- Low- versus high-seed halves agree only moderately on cell drift
  (13 common cells, correlation 0.547; RMSE 0.022 versus pooled future-drift
  standard deviation 0.064).

Therefore $(q,J_a)$ has a moderately reproducible coarse conditional tendency,
but it is not a validated autonomous Markov state at the tested resolution.
Future-outcome dependence remains after conditioning, while direct history
tests are mixed rather than evidence for strong event-trajectory memory.  The
likely missing state is exactly what the microscopic balance exposes:
approach-region population/geometry and its boundary renewal, plus local
heading organization.  This result does not justify adding a fitted coordinate
here; it specifies why the proposed two-variable closure is presently
insufficient without claiming that history alone is the cause.

## Verdict

1. **Microscopic origin:** social turning and recruitment of already directed
   fish build pre-crossing $J_a$.  Wall turning and angular-noise relaxation
   dissipate it.  Realized noise and collision overrides are not primary.
2. **Propagation:** $J_a$ orders contemporaneous neck transport, but the
   parameter-free plug-flow magnitude and tagged travel-time prediction fail.
   Exact accumulation and side-boundary transport cannot be discarded.
3. **Closure:** $(q,J_a)$ is useful and interpretable, but only moderate
   split-seed drift reproducibility and conditional outcome dependence prevent
   treating it as a validated Markov state.  Strong memory in transition paths
   is not established.

The result narrows the physical picture without inventing a potential: a
socially organized, renewing approach layer feeds the gap, while geometry and
wall/noise relaxation determine whether that layer persists.

## Reproducibility and limitations

- `analyze_approach_current_balance.py` contains all replay, balance, flux, and
  closure calculations.  `--reuse-balance-closure` resumes from compact
  source-replay tables without recomputing torques.
- `native_balance_*` contains the exact eight-pair replay and audit.
- `plug_flow_*` contains held-out seed samples and seed summaries.
- `closure_*` contains cell-level and independent-seed diagnostics.
- `summary.json`, `analysis_configuration.json`, and
  `figures/approach_current_mechanism.png` are the compact entry points.

The exact torque decomposition is intentionally limited to eight independent
pairs because source replay is the expensive step.  Its intervals are therefore
wide.  The held-out plug samples overlap within episodes; inference is reported
at seed level rather than pretending 1,777 independent trials.  Matched
controls were selected to have a stable future, so the event/control comparison
answers a conditional case-control question, not an absolute crossing hazard.
No raw MAT file or cache belongs in this output.
