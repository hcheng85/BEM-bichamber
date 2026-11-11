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
contains all 40 independent raw seeds (`n0`--`n39`).  Existing exact processed
occupancy indices were used for `n0`--`n19`; `n20`--`n39` were indexed directly
from raw positions with a bounded, sequential reader.  The first 200 time units
of each trajectory were excluded.  The complete initial 20-seed result is
preserved in `../initiation_propagation_20seed_reference/`.

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
The RMS collision residual is 0.0019 at lead 2 and 0.0068 at lead 0.1.  The
localized `J_a` adds spatial information that `a` lacks.

## Sample and controls

The analysis found 658 post-transient episodes in 37 event-producing seeds:
567 have positive direction-aligned population change at `+5`, while 91 have
zero or negative change and are labeled reversals.  Median episode duration is
2.0 time units, median gross crossing count is 7, and median directional net
change at `+5` is 4 fish.  Seed event counts are strongly unequal, so fish,
frames, and episodes were not treated as independent replicates.

For each event and lead time, the control is drawn from a different seed at
exactly the same discrete `q`, using the same prospective direction convention,
and has no occupancy change for two units before through five units after the
control time.  There are 542 same-event pairs available at all leads from 5 to
0.1; 641 pairs are available from lead 2 onward.  Confidence intervals resample
the 37 initiating seeds and therefore represent seed-cluster, not frame-level,
uncertainty.

This control answers a specific question: among configurations with the same
population imbalance, what distinguishes a future crossing from a stable
no-flux configuration?  It is intentionally not an estimate of an unconditional
hazard, and the direction of the future event is used to sign the observables.

## Result 1: a collective directional reversal precedes the first crossing

On the fixed 542-event cohort, the same-`q` discrimination is:

| lead before first crossing | polarization `a` AUC (95% CI) | approach current `J_a` AUC (95% CI) |
|---:|---:|---:|
| 5 | 0.433 (0.400, 0.467) | 0.459 (0.445, 0.473) |
| 4 | 0.469 (0.432, 0.500) | 0.472 (0.452, 0.491) |
| 3 | 0.571 (0.537, 0.606) | 0.578 (0.554, 0.603) |
| 2 | 0.683 (0.649, 0.716) | 0.758 (0.728, 0.789) |
| 1 | 0.753 (0.718, 0.786) | 0.918 (0.897, 0.936) |
| 0.5 | 0.779 (0.747, 0.809) | 0.947 (0.931, 0.963) |

At lead 5, future crossing configurations are still polarized **away** from
their eventual crossing direction.  They are statistically neutral near lead
4, reverse significantly toward the future crossing direction by lead 3, and
are clearly gap-directed by lead 2.  On the common cohort at lead 2, the mean
event-control difference in `a` is 0.219 (95% CI 0.151, 0.280).  The lead-2
seed-mean difference is positive in 32 of 37 seeds for `a` and 34 of 37 for
`J_a`; at lead 1 both are positive in 35 of 37 seeds.
The result is therefore not produced by only one high-event-count trajectory.
It also survives the prespecified ensemble extension: relative to the preserved
20-seed analysis, the common-cohort lead-2 AUC changes from 0.647 to 0.683 for
`a` and from 0.742 to 0.758 for `J_a`, without changing the inferred timing.

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
| 0--0.2 | 0.650 (0.606, 0.691) | 0.231 (0.148, 0.309) | 0.190 (0.101, 0.268) | 0.330 (0.265, 0.395) |
| 0--0.5 | 0.789 (0.750, 0.826) | 0.470 (0.393, 0.536) | 0.356 (0.259, 0.438) | 0.558 (0.487, 0.620) |
| 0--1.0 | 0.878 (0.847, 0.908) | 0.564 (0.512, 0.613) | 0.431 (0.344, 0.506) | 0.611 (0.551, 0.663) |

The 0--0.2 neck current already carries modest information about later cluster
growth and continuation.  A current sustained for 0.5--1.0 time units is much
more informative.  These longer-window values are **ongoing-transition
commitment variables**, not pre-initiation precursors.  Their increasing skill
means that persistent transport distinguishes continuation from reversal; it
does not show that current caused the initial turn.
The continuation AUCs likewise change only from 0.635/0.766/0.852 in the
20-seed reference to 0.650/0.789/0.878 in the 40-seed analysis for the
0.2/0.5/1.0 windows.

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

- There are 658 episodes but only 37 independent event-producing seeds; all
  uncertainty is therefore seed-clustered.
- The second half of the ensemble was indexed directly from raw positions;
  processed occupancy files were not required or fabricated.
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
- `figures/`: chronology, uncertainty, propagation, and representative traces;
  `common_cohort_predictive_effects.png` corresponds to the fixed-cohort table
  reported above.

No raw trajectory or cache file is intended for version control.
