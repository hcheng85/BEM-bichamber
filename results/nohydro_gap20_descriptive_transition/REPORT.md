# First descriptive analysis: no-hydrodynamics collective bichamber crossings

> **Raw-trajectory extension (10 November 2026).**  The gap-20 raw trajectories
> subsequently became available.  The initiation/propagation analysis is in
> [`outputs/initiation_propagation/REPORT.md`](outputs/initiation_propagation/REPORT.md).
> It supersedes only the data-access statement below; the first-stage processed
> occupancy results are retained unchanged for provenance.

## Scope

This is the bichamber analysis, not the circular schooling--milling analysis.
The physical event of interest is a **collective chamber crossing**.  Schooling
and milling labels from circular confinement are not reused here.

The present calculation uses the locally available processed ensembles:

- `R=50`, `N=100`, `Ip=4`, `In=0.5`, `Iw=0.94`, `If=0`;
- gap 20 (primary larger-gap case) and gap 10 (comparison);
- 20 independent seeds at each gap, duration 1000 per seed;
- the first 200 time units are excluded.

The processed files contain the exact chamber population `N_L(t)`, accepted
event markers, and the archived net event mark.  They do **not** contain fish
positions, headings, velocities, or fish identities.  Consequently this stage
can establish event statistics and distinguish cluster crossings from full
school relocation, but cannot yet identify a geometric or orientational
mechanism.

## Definitions

The chamber imbalance is

\[
q(t)=\frac{N_R-N_L}{N}=1-\frac{2N_L}{N}.
\]

An archived accepted event is the output of `GetTransitEvent` in
`utilites_histogram.py`: nearby changes are merged over 2 time units and a
5-time-unit before/after occupancy check rejects short returns.  Its mark
`Delta N_L` is the net occupancy change between successive accepted markers.
This mark is a crossing-cluster statistic; it is not automatically a full
school relocation.

A strict collective relocation is defined independently as a first passage
from `q <= -0.6` to `q >= 0.6`, or the reverse.  The relocation duration is the
time between the last sample in the origin basin and the first sample in the
destination basin.  This threshold is prespecified and the script exposes it
as an argument.

## Results

### Gap 20 supplies more crossing statistics, but few full relocations

After transient removal, gap 20 has 283 accepted event markers, compared with
78 at gap 10.  The mean seed-level event rates are

| gap | accepted-event rate | 95% seed-bootstrap interval |
|---:|---:|---:|
| 10 | 0.00487 | 0.00312--0.00669 |
| 20 | 0.01769 | 0.01131--0.02437 |

For nonzero marks, gap 20 has 268 clusters with pooled mean size 9.14 fish and
median 6; 101 contain at least 10 fish, 36 at least 20 fish, and the largest is
46 fish.  Gap 10 has 68 nonzero marks with pooled mean 4.82, median 2, and
maximum 32.  The larger opening therefore increases both accepted-event
frequency and typical net cluster size.

This does **not** yield hundreds of whole-school transitions.  With the strict
`|q|=0.6` definition, gap 20 contains only five relocations, arising from three
of the 20 seeds; gap 10 contains none.  Gap-20 relocation durations are 17.7,
80.7, 96.9, 346.5, and 465.3 time units.  The median is 96.9, but five events
are far too few for a reliable duration law or conditional transition model.

### The stationary occupancy distribution is not a transition mechanism

Gap 20 spends a seed-averaged fraction 0.528 of post-transient time at
`|q|>0.8` and 0.170 at `|q|<0.2`.  Gap 10 gives 0.416 and 0.197 respectively,
with wide seed-to-seed intervals.  The bimodality is consistent with a school
often occupying one chamber, but it cannot by itself define a physical
potential or explain why a crossing begins.

The event-aligned occupancy plot is direction-normalized using the sign of the
nonzero mark.  It is a descriptive chronology only: the marker is constructed
from occupancy changes, so behavior at or after time zero cannot be treated as
a predictor.

## Interpretation

The most important result of this first pass is a distinction of sample size:

1. Gap 20 has enough data to study **cluster crossing statistics** (268
   nonzero post-transient marks in the currently local 20 seeds).
2. It does not yet have enough data to estimate a robust law for **complete
   collective relocation** (five strict relocations).
3. Frames are highly correlated and must not be counted as independent
   transition samples; the independent replication unit is the seed.

Thus, the phrase "more samples" is correct for crossing clusters, but not yet
for full school switching.  A one-dimensional occupancy analysis alone cannot
answer whether approach, collective heading, leading-edge nucleation, or local
deformation causes those events.

## Data-access audit and bounded next stage

The extension inventory records completed raw trajectories for gap 10 and gap
20, `In=0.5`, seeds 20--39, which would double each ensemble to 40 seeds.  It
also records 20 raw `In=0.7` trajectories at each gap.  These files reside on
the workstation/external-data path and are not mounted in the present session.
An attempted key-only SSH connection to the recorded workstation was rejected,
so no raw trajectory was silently substituted or invented.  The added seeds
also have not been event-postprocessed into the local `Processed_data` folder.

Once that raw path is mounted or key access is enabled, the next bounded
analysis should process gap 20, `In=0.5`, all 40 seeds first and compute, before
each event:

- global and chamber-local polarization;
- school COM distance and direction relative to the opening;
- near-gap and leading-edge population;
- spatial extent parallel and normal to the divider;
- heading toward the admissible opening and cluster identity;
- event duration, cluster growth, and propagation, with prospective stable
  controls matched by seed, origin chamber, vertical phase, and proximity.

The `In=0.7` ensemble should remain a separate noise condition, not be pooled
to inflate the sample count.  Full relocations and accepted cluster events
should also remain separate outcomes.

## Reproducibility

Run from the repository checkout:

```bash
python results/nohydro_gap20_descriptive_transition/analyze_nohydro_collective_transition.py \
  --processed-dir ../Processed_data \
  --output-dir results/nohydro_gap20_descriptive_transition/outputs \
  --gaps 10 20
```

Compact outputs include `summary.csv`, `per_seed_summary.csv`,
`accepted_cluster_events.csv`, `dominance_relocations.csv`, the analysis
configuration, and three figures.  Original simulations and processed files
are unchanged.
