# First-stage analysis of schooling--milling dynamics in circular confinement

## Scope and data audit

This analysis uses the ten files
`If1e-02_Ip40_In05_R50_num100_wid100_n{0..9}_BEM.mat` in
`/Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data`.  The files were read
only; none was modified.  Each is an uncompressed MATLAB-v5 file of about
4.0 GB, so `state` and `rdot` were accessed with read-only Fortran-order memory
maps and 5,000-frame chunks.  The entire trajectories were never loaded into
RAM.

The production source is `BEM_circular_parallel_instance.py` (SHA-256
`2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`).
The filename's `R50` is the radius *ratio*; the nondimensional simulated radius
is `R=7.158910531638177`.  The audited parameters are `N=100`, `Ip=4`,
`In=0.5`, `Iw=0.94`, `If=0.01`, `dt=0.001`, and duration 1000.

The source and direct numerical checks establish that

- `state[0,i,t]`, `state[1,i,t]`, and `state[2,i,t]` are fish $i$'s
  $x$, $y$, and heading angle;
- `rdot[:,i,t]` is the **actual stored translational velocity after collision,
  extreme-velocity, and positional safety overrides**.  For $t\geq1$, it is
  the velocity used to advance `state[...,t-1]` to `state[...,t]`, rather than
  a freshly evaluated velocity at `state[...,t]`.

## Order parameters

At 0.01-time-unit cadence I computed

\[
P=\left|N^{-1}\sum_i(\cos\theta_i,\sin\theta_i)\right|,\qquad Q=P^2,
\]

the simulator's original milling statistic

\[
\Pi=N^{-1}\sum_i\left|\widehat{(r_i-r_{\rm COM})}
             \times\widehat{\dot r_i}\right|,
\]

and signed milling about the fixed circle center

\[
M=N^{-1}\sum_i \hat r_i\times e(\theta_i).
\]

A velocity-based signed fixed-center milling statistic was also saved as a
diagnostic.  Particle-level radial and heading-minus-radial-angle distributions
were accumulated every 0.1 time units.

The archived `GetGlobalParameters` helper contains legacy hard-coded dimensions
of 10 inside `Reshape` and `np.tile`, so it cannot literally process these
`N=100` arrays.  The calculation above is the same mathematical definition,
generalized vectorially to the actual population.  This helper is not called by
the production time step, so the issue does not alter the simulations.  Also,
the helper's signed `Md` is COM-relative and velocity-based; the fixed-circle,
heading-based $M$ requested here is deliberately a separate quantity.  For
fidelity to that helper, $\Pi(t)$ pairs `state[...,t]` with the stored
`rdot[...,t]`, despite the one-index update convention noted above; the offset
is only one native step (0.001).

For event description only, a one-time-unit moving average was labeled
schooling when $P\geq0.70,\Pi\leq0.75$, milling when
$P\leq0.50,\Pi\geq0.75$, and intermediate otherwise.  A reported transition
requires 90% occupancy of the old and new labels in consecutive 20-time-unit
windows.  These thresholds identify the two visibly separated modes; they are
not a fitted dynamical model.

## Collective states and transitions

After discarding the first 20 time units, the seed-time fractions are 0.727
schooling, 0.270 milling, and 0.0024 intermediate.  Seeds 0--6 remain
predominantly schooling; seeds 7 and 8 remain milling; seed 9 changes from
schooling to milling.  The single persistent transition is in seed 9 at the
nominal threshold time 312.41.  A looser but still separated threshold places
the same transition at 311.3.  A very strict threshold loses the formal event
because the post-transition polarization intermittently exceeds its cutoff, so
the exact onset is threshold-dependent even though the state change is clear.

Only one sustained transition occurs in 9,800 post-formation seed-time units.
Consequently these data reveal coexistence/bistability at this parameter point,
but cannot estimate a transition-rate law or locate a stability boundary.

The particle distributions distinguish the modes physically.  Milling has a
compact annulus with radial-density maximum near $r/R=0.233$; schooling is
farther out, peaking near $r/R=0.500$.  Headings in both modes favor tangential
directions near $\pm\pi/2$, but the schooling peaks are sharper.  A high
value of signed $M$ is therefore not, by itself, a milling identifier: a
polarized school can orbit coherently around the circular wall.

## Exact microscopic torque replay

I replayed 170 state-stratified snapshots.  Schooling and milling snapshots
were required to be locally persistent for ten time units; short threshold
excursions were treated as transitional.  The replay uses the original
Delaunay visible-neighbor graph, angular visibility weights, direct dipolar
hydrodynamics, BEM solve, circular-ray wall interaction, clipping, and update
overrides.  The social law was split algebraically into its alignment and
attraction terms without changing their sum.

Validation is numerical rather than assumed:

- maximum error between split and direct social torque: $5.4\times10^{-14}$;
- maximum error between replayed and stored translational velocity:
  $3.6\times10^{-15}$;
- no selected snapshot used collision, extreme-velocity, or sticky overrides,
  so the selected torque decompositions need no ambiguous clipped attribution.

For every component $k$, the exact projected deterministic contribution is

\[
b_Q^{(k)}={2\over N}\sum_i[p\cdot e_{\perp i}]\,\Omega_i^{(k)},
\]

and angular diffusion supplies

\[
b_Q^{\rm noise}=-I_n^2(Q-1/N).
\]

Seed-bootstrap means (95% intervals are in `summary.csv`) show the following
balance:

| state | alignment | attraction | wall | direct hydro | BEM hydro | noise | total $b_Q$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| schooling | +0.580 | -0.018 | +0.012 | -0.376 | +0.00006 | -0.201 | -0.003 |
| intermediate | +0.419 | -0.017 | -0.027 | -0.206 | -0.0010 | -0.144 | +0.023 |
| milling | +0.041 | -0.003 | +0.001 | +0.006 | +0.000002 | -0.011 | +0.033 |

The strongest first-principles result is therefore a torque balance, not a
potential: social alignment creates polarization, while direct fish--fish
hydrodynamics and angular-noise decorrelation oppose it in the schooling state.
Boundary-mediated angular hydrodynamics is negligible in this particular
projection and parameter set.  In milling, all projected terms are much weaker;
the small positive total drift is consistent with restoration away from the
finite-$N$ low-$Q$ floor rather than an instability toward schooling.

## Does $(Q,M)$ close the dynamics?

This is a preliminary rejection/validation screen, not a fitted surrogate.
Snapshots were conditioned in fixed $Q$-$M$ cells, with cell fixed effects
and within-cell linear $Q,M$ gradients removed.  Seed-cluster bootstrap
correlations between the remaining $b_Q$ residual and recent
$\Delta Q$, radial spread, center radius, aspect ratio, nearest-neighbor
spacing, neighbor alignment, and $\Pi$ all had 95% intervals containing zero.
Thus this bounded sample does **not detect** a reproducible conditional-mean
dependence on those individual variables.

However, 94% of instantaneous $b_Q$ variance remains within the rather broad
$(Q,M)$ cells.  Such variance alone does not disprove Markov closure--a
closed stochastic coordinate may have conditional variability--but it means
that the present 170 snapshots cannot establish closure either.  There are
only ten independent seeds and only one persistent transition.  The defensible
conclusion is therefore: **$(Q,M)$ closure is unresolved, not validated.**
The one observed transition also shows signed-circulation reorganization before
the large fall of $Q$, but one event cannot establish it as a general
precursor.

## Wall-safety audit

The production wall law amplifies torque when the forward-ray intersection
distance satisfies $|d|\leq0.1$, and later rejects outward positional steps
inside the 0.1 boundary layer.  Across all one billion fish-time samples:

- the high-order forward-ray condition is triggered once (seed 0), with
  minimum $d=0.0972446$;
- 9,306 fish-steps have the unique near-wall zero-velocity signature of the
  positional safety override (fraction $9.31\times10^{-6}$);
- no saved point lies outside the circle; the largest radius is 7.09607, below
  $R=7.15891$.

The original wall formula and safety condition were not modified.  The
high-order count is evaluated directly from saved position and heading.  The
positional count is inferred from the stored post-override velocity signature;
the rejected trial velocity itself was not saved, so this qualification should
remain attached to that count.

## Data needed to predict a stability boundary

The current data are a single noisy parameter point with random initial
conditions.  A publishable stability boundary requires, at minimum:

1. **More independent trajectories at this point** (roughly 50 or more, or
   enough to observe tens of sustained switches), because one transition cannot
   constrain rates or uncertainty.
2. **Prepared schooling and milling initial conditions on both sides of a
   targeted parameter bracket**, rather than a broad blind sweep.  Small grids
   in the candidate control parameters $I_p/I_n$, followed by $I_f$ and
   $I_w$ sensitivity, would separate deterministic loss of stability,
   stochastic switching, and basin-selection probability.
3. **Deterministic/noise-reduced perturbation runs** from converged rotating
   schooling and milling states.  Their decay/growth rates or Floquet/Jacobian
   spectra are the data needed for a true deterministic stability boundary.
4. **Denser matched samples at the same $(Q,M)$** but with different radial
   spread, aspect ratio, circulation direction, and recent history.  This is
   required to validate or reject the proposed coarse Markov closure before
   writing any two-dimensional Fokker--Planck model.
5. Saving the per-component angular torques, neighbor graph, and override flags
   during future runs would make this test much cheaper and would remove the
   remaining ambiguity in retrospectively identifying positional safety events.

No effective potential was fitted in this analysis.

## Reproducibility map

- `process_collective.py`: memory-mapped chunk reduction and safety scan.
- `analyze_collective.py`: state/event definitions and distributions.
- `replay_torques.py`: exact source replay and projected-drift decomposition.
- `analyze_drift_closure.py`: seed-bootstrap torque summary and closure screen.
- `outputs/metrics_seed*.npz`: compact per-seed order-parameter series.
- `outputs/summary.csv`, `outputs/closure_screen.csv`: primary numerical results.
- `outputs/torque_snapshot_summary.csv`: snapshot-level replay results.
- `outputs/streaming_inventory.csv`: file/offset/safety audit.
- `outputs/*.png`: decisive plots.
