# Initiation and propagation of no-hydrodynamics gap-20 crossings

## Scientific question and scope

This analysis asks two separate questions in the `N=100`, `In=0.5`, nominal
gap-20 bichamber ensemble:

1. What changes **before the first fish crosses** the divider?
2. Once crossing begins, what distinguishes continued collective transport
   from a short excursion or reversal?

The distinction between chamber population, directed transport, and reversal
is motivated by [Li et al., *Soft Matter* (2022)](https://pubs.rsc.org/en/content/articlelanding/2022/sm/d2sm00134a).
Their experiment concerns active particles in a long channel.  No force-balance
law, oscillation model, or fitted coefficient from that work is imposed on this
two-chamber geometry.

The production simulator and raw MAT files were read only.  The primary sample
is the 20 independent seeds (`n0`--`n19`) for which exact processed occupancy
indices and raw trajectories are both available.  Raw `n20`--`n39` files also
exist, but they do not have corresponding processed event indices and were not
silently mixed into this bounded calculation.  The first 200 time units of
each trajectory were excluded.

## Definitions and source audit

For a left-to-right event let `s=+1`; for right-to-left let `s=-1`.  With
`N_L+N_R=N`,

\[
q=(N_R-N_L)/N.
\]

An episode starts at the first native-resolution centerline crossing following
at least two time units with no centerline crossing.  Subsequent occupancy
changes separated by at most two time units are grouped into the episode.  The
first-crossing time is therefore a physical initiation marker, not a coarse
accepted-event marker.  These are crossing episodes and include isolated or
reversing attempts; they are not all whole-school relocations.

Immediately before initiation, source-chamber gap-directed polarization is

\[
a=\frac{1}{N_s}\sum_{i\in\mathrm{source}}s\cos\theta_i.
\]

The direction-aligned source current and localized approach current are

\[
J_s=\frac{1}{N}\sum_{i\in\mathrm{source}}s v_{x,i},\qquad
J_a=\frac{1}{N}\sum_{i\in\mathrm{approach}}s v_{x,i},
\]

where the prespecified approach zone is on the source side, within one opening
width of the divider, and within half an opening width vertically.  The
post-initiation neck current over a window `w` is the direction-aligned finite
change of `q/2` divided by `w`.

The stored-velocity convention was checked against the simulator:
`rdot[:,:,k+1]` advances `state[:,:,k]`.  In the no-hydrodynamics simulation,
`v_i=(cos(theta_i),sin(theta_i))` except where the production collision override
sets a velocity to zero.  Consequently,

\[
J_s=\frac{N_s}{N}a+J_{\rm collision}.
\]

Thus `a` and `J_s` are **not independent physical mechanisms** in this ensemble.
The RMS collision residual is 0.0023 at lead 2 and 0.0074 at lead 0.1.  The
localized `J_a` adds spatial information that `a` lacks.

## Sample and controls

The analysis found 338 post-transient episodes in 19 event-producing seeds:
292 have positive direction-aligned population change at `+5`, while 46 have
zero or negative change and are labeled reversals.  Median episode duration is
2.0 time units, median gross crossing count is 6, and median directional net
change at `+5` is 3 fish.  Seed event counts range from 1 to 40, so fish,
frames, and episodes were not treated as independent replicates.

For each event and lead time, the control is drawn from a different seed at
exactly the same discrete `q`, using the same prospective direction convention,
and has no occupancy change for two units before through five units after the
control time.  There are 247 same-event pairs available at all leads from 5 to
0.1; 304 pairs are available from lead 2 onward.  Confidence intervals resample
the 19 initiating seeds and therefore represent seed-cluster, not frame-level,
uncertainty.

This control answers a specific question: among configurations with the same
population imbalance, what distinguishes a future crossing from a stable
no-flux configuration?  It is intentionally not an estimate of an unconditional
hazard, and the direction of the future event is used to sign the observables.

## Result 1: a collective directional reversal precedes the first crossing

On the fixed 247-event cohort, the same-`q` discrimination is:

| lead before first crossing | polarization `a` AUC (95% CI) | approach current `J_a` AUC (95% CI) |
|---:|---:|---:|
| 5 | 0.429 (0.362, 0.478) | 0.445 (0.417, 0.471) |
| 4 | 0.514 (0.433, 0.590) | 0.494 (0.468, 0.523) |
| 3 | 0.590 (0.540, 0.640) | 0.570 (0.529, 0.609) |
| 2 | 0.647 (0.593, 0.699) | 0.742 (0.704, 0.781) |
| 1 | 0.749 (0.698, 0.796) | 0.916 (0.885, 0.940) |
| 0.5 | 0.811 (0.770, 0.844) | 0.952 (0.927, 0.977) |

At lead 5, future crossing configurations are still polarized **away** from
their eventual crossing direction.  They are statistically neutral near lead
4, show a weak sign reversal around lead 3, and are clearly gap-directed by
lead 2.  At lead 2, the mean event-control difference in `a` is 0.135 (95% CI
0.035, 0.233).  The lead-2 seed-mean difference is positive in 16 of 19 seeds
for `a` and 17 of 19 for `J_a`; at lead 1 both are positive in 17 of 19 seeds.
The result is therefore not produced by only one high-event-count trajectory.

All of these measurements precede the first centerline crossing, so `a` has
genuine prospective information before neck flux exists.  However, the data do
**not** establish `a` as an earlier trigger than source-side translation:
`J_a` changes over the same lead-3 to lead-2 interval and is more discriminating
thereafter.  The defensible precursor is therefore a coupled reorientation and
localized approach, not polarization alone.

The exact matching on `q` also answers the population-imbalance question.  `q`
describes how loaded each chamber is and defines the available source
population, but the same `q` supports both future crossing and stable no-flux
configurations.  It is not sufficient to initiate or predict an episode by
itself.

## Result 2: current persistence separates continuation from reversal

To avoid correlating a current window only with transport already counted in
that same window, the primary outcomes subtract the population change and time
already observed during the window.  Results are:

| current window after initiation | AUC: continuation vs reversal | correlation with remaining net fish | remaining duration | remaining gross crossings |
|---:|---:|---:|---:|---:|
| 0--0.2 | 0.635 (0.565, 0.685) | 0.261 (0.162, 0.364) | 0.239 (0.132, 0.332) | 0.362 (0.264, 0.452) |
| 0--0.5 | 0.766 (0.719, 0.804) | 0.519 (0.432, 0.592) | 0.433 (0.308, 0.529) | 0.601 (0.503, 0.686) |
| 0--1.0 | 0.852 (0.807, 0.886) | 0.589 (0.526, 0.640) | 0.461 (0.318, 0.564) | 0.630 (0.548, 0.692) |

The 0--0.2 neck current already carries modest information about later cluster
growth and continuation.  A current sustained for 0.5--1.0 time units is much
more informative.  These longer-window values are **ongoing-transition
commitment variables**, not pre-initiation precursors.  Their increasing skill
means that persistent transport distinguishes continuation from reversal; it
does not show that current caused the initial turn.

## Physical interpretation and stopping conclusion

The existing trajectories support the following chronology:

1. At fixed chamber imbalance, the school is initially neutral or oriented
   away from the eventual crossing direction.
2. Roughly three to two time units before the first fish crosses, source-side
   heading reverses toward the gap and localized approach begins.
3. The first fish crosses.
4. Weak early neck current sometimes decays or reverses; current that persists
   for 0.5--1.0 predicts a larger, longer episode and continued net transport.

The evidence therefore supports **directional reorientation plus spatial
localization as initiation**, followed by **current persistence as propagation
and commitment**.  It does not support a one-dimensional population-only
mechanism, a fitted potential, or importing a long-channel force-balance law.
Because heading and approach current appear together at the present resolution,
the microscopic cause of the turn (wall, social alignment, attraction, or
noise) remains unresolved by this bounded analysis.

## Limitations

- There are 338 episodes but only 19 independent event-producing seeds; all
  uncertainty is therefore seed-clustered.
- The primary analysis uses the 20 seeds with processed event indices.  The 20
  additional raw seeds should be indexed before claiming a 40-seed result.
- Direction signing uses the subsequent first-crossing direction.  The reported
  AUCs quantify direction-conditioned precursors, not an online two-direction
  predictor.
- Controls are deliberately stable for a five-unit future window; they answer a
  case-control question rather than the absolute crossing hazard.
- `a` and full source current are kinematically linked in the no-hydrodynamics
  simulator.  Their agreement is an identity audit, not two independent pieces
  of evidence.
- Associations with post-initiation current are predictive, not torque-source
  ablations or causal interventions.

## Reproducibility map

- `analyze_initiation_propagation.py`: raw read-only extraction, event/control
  construction, uncertainty, and figures.
- `summarize_common_precursor_cohort.py`: fixed-cohort sensitivity analysis.
- `crossing_episodes.csv`: episode features and future outcomes.
- `matched_precursor_pairs.csv`: exact-`q` event/control pairs.
- `statistical_summary.csv`: primary seed-bootstrap estimates.
- `common_cohort_summary.csv`: the same event cohort at all lead times.
- `sample_and_identity_audit.json`: counts and no-hydrodynamics current identity.
- `analysis_configuration.json`: paths, indexing, and prespecified definitions.
- `figures/`: chronology, uncertainty, propagation, and representative traces.

No raw trajectory or cache file is intended for version control.
