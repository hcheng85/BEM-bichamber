# Bottleneck control of alignment-induced population localization

## Decision

**Go:** the archived width series establishes a reproducible, nontrivial
geometric effect.  Strong alignment does not produce the same macroscopic
population state independent of the chamber opening.  In the
no-hydrodynamics ensemble, late localization is maximal at intermediate gaps
(`gap_ratio=30--60`), whereas population redistribution and complete
opposite-dominance relocations increase strongly with opening width.  The
matching hydrodynamic archives reproduce this qualitative structure while
substantially increasing redistribution.

This is evidence for **bottleneck-controlled collective trapping**, not just
generic school cohesion.  It is not evidence for stationary equilibrium
bistability: the narrow-gap occupation distributions retain initialization
memory, and no potential or switching rate is inferred.

## Data and compatibility

The analysis uses every compatible archived trajectory at `Ip=4`, `In=0.5`,
`N=100`, radius ratio 50, and duration 1000:

- no hydrodynamics: 40 seeds each at gap ratios 10 and 20, and 20 seeds each
  at 30, 40, 50, 60, 70, 80, and 85 (`220` trajectories);
- hydrodynamics: 20 seeds at each of the nine widths (`180` trajectories).

All trajectories use the production all-left initialization.  For a common
seed, the initial microscopic state is identical across no-hydrodynamic
widths and agrees between matching hydrodynamic and no-hydrodynamic files to
`1.78e-15`.  The hydrodynamic comparison holds `Ip`, `In`, `N`, geometry,
duration, and initialization fixed; its archived coupling is `If=1e-2`
instead of `If=0`.

Population counts are sampled every `0.1` time unit.  The primary late window
is `t=600--1000`; uncertainty resamples complete independent archive seeds.
The production sources are unchanged.  Their SHA256 hashes and the complete
archive audit are in `source_and_data_audit.json`.

## No-hydrodynamics width dependence

Late equal-seed occupation statistics are:

| Gap ratio | gap/R | `<q^2>` (95% seed CI) | `P(|q|>=0.8)` | `<|q(t+10)-q(t)|>` (95% CI) |
|---:|---:|---:|---:|---:|
| 10 | 0.2 | 0.505 [0.379, 0.635] | 0.421 | 0.0036 [0.0019, 0.0059] |
| 20 | 0.4 | 0.574 [0.450, 0.698] | 0.518 | 0.0313 [0.0226, 0.0405] |
| 30 | 0.6 | 0.862 [0.781, 0.932] | 0.838 | 0.0417 [0.0246, 0.0598] |
| 40 | 0.8 | 0.851 [0.805, 0.895] | 0.826 | 0.0946 [0.0653, 0.1266] |
| 50 | 1.0 | 0.883 [0.854, 0.911] | 0.869 | 0.1525 [0.1148, 0.1939] |
| 60 | 1.2 | 0.837 [0.811, 0.864] | 0.805 | 0.2423 [0.2040, 0.2821] |
| 70 | 1.4 | 0.807 [0.731, 0.857] | 0.772 | 0.3314 [0.2906, 0.3712] |
| 80 | 1.6 | 0.686 [0.587, 0.771] | 0.641 | 0.5556 [0.4986, 0.6095] |
| 85 | 1.7 | 0.674 [0.622, 0.707] | 0.620 | 0.8667 [0.7946, 0.9279] |

The localization curve is nonmonotonic.  In a paired comparison using the 20
seeds shared by all widths, gap 50 exceeds gap 10 in `<q^2>` by `0.352`
(`[0.161,0.532]`) and exceeds gap 85 by `0.209`
(`[0.163,0.262]`).  The analogous strong-localization differences are
`0.406` (`[0.190,0.622]`) and `0.250` (`[0.198,0.310]`).  Thus neither a
nearly closed bottleneck nor the broadest opening maximizes single-chamber
occupation.

Population redistribution is instead strongly increasing across the series.
The paired gap-85-minus-gap-10 difference is `0.862`
(`[0.790,0.924]`).  This separates an intermediate-width trapping regime
from a broad-opening mobile regime.

## The mobility trend is not an occupancy-range artifact

A raw mobility average can be small merely because a trajectory spends most
of its time near `|q|=1`, where motion in one direction is bounded.  We
therefore repeated the fixed-lag comparison conditional on the starting
population range, without introducing a new state coordinate:

- central, `|q|<0.4`: gap-85-minus-gap-10 mobility is `0.800`
  (`[0.755,0.840]`);
- intermediate, `0.4<=|q|<0.8`: `0.785` (`[0.626,0.895]`);
- strong, `|q|>=0.8`: `0.975` (`[0.896,1.059]`).

All contrasts are paired by archive seed.  The broad-gap enhancement remains
large in every fixed population band, so different occupancy of `q` cannot
by itself explain the geometric mobility trend.

## Temporal convergence and initialization memory

All runs start at `q=-1`.  Late signed means approach reflection symmetry at
many widths, but several ensembles retain appreciable seed trapping; for
example the late no-hydrodynamic mean is `-0.239` at gap 10 and `-0.267` at
gap 50.  Equal-seed distribution distances between `t=200--600` and
`600--1000` are `0.383`, `0.232`, and `0.212` at gaps 10, 20, and 30, then
fall to `0.06--0.14` at gaps 40--85.  Hydrodynamic distances show the same
pattern: `0.374` at gap 10 and `0.040--0.103` at gaps 40--85.

The wide-gap results are much closer to time-window convergence, but the
narrow-gap histograms are demonstrably finite-time occupation measures.
Consequently, the observed peaks in `P(q)` are not called stationary wells,
and `-log P(q)` is not interpreted as an energy.

## Complete opposite-dominance relocations and censoring

Using the established `q<=-0.6` to `q>=0.6` (or reflected) definition over
`t=200--1000` gives:

| Gap | No-hydro events (seeds) | Hydro events (seeds) |
|---:|---:|---:|
| 10 | 0 (0/40) | 22 (13/20) |
| 20 | 9 (6/40) | 75 (20/20) |
| 30 | 13 (10/20) | 198 (20/20) |
| 40 | 36 (16/20) | 307 (20/20) |
| 50 | 68 (19/20) | 400 (20/20) |
| 60 | 112 (20/20) | 497 (20/20) |
| 70 | 214 (20/20) | 618 (20/20) |
| 80 | 379 (20/20) | 1149 (20/20) |
| 85 | 631 (20/20) | 1353 (20/20) |

Every trajectory contributes a right-censored terminal interval.  At narrow
gaps, the event count is absent or sparse and most seeds show no relocation;
we therefore report counts and censoring rather than fitting a switching rate.
The increasing counts support the fixed-lag mobility result but do not turn
the nonstationary narrow-gap process into a stationary renewal process.

## Hydrodynamic comparison

Hydrodynamics preserves the nonmonotonic localization curve: late `<q^2>`
rises from `0.539` at gap 10 to `0.773` at gap 50 and falls to `0.610` at gap
85.  Paired seed contrasts resolve both sides of the peak: gap 50 minus gap 10
is `0.233` (`[0.123,0.345]`), and gap 50 minus gap 85 is `0.163`
(`[0.097,0.247]`).

Hydrodynamics strongly enhances redistribution at every width: late mobility
rises from `0.074` at gap 10 to `1.276` at gap 85.  The gap-85-minus-gap-10
paired difference is `1.202` (`[1.031,1.337]`) and remains resolved within
each fixed `q` band.  Thus hydrodynamic coupling shifts the transport scale
but does not create the geometric phenomenon or remove its bottleneck
dependence.

## Minimal physical hypothesis

The data support a three-regime aperture-encounter hypothesis:

1. A very narrow aperture is encountered too rarely within a collective
   persistence time, so relaxation from the initial population partition is
   extremely slow.
2. At intermediate width, coherent transfer is possible but return encounters
   remain infrequent, producing maximal single-chamber trapping.
3. At broad width, repeated encounters destroy long residence in either
   chamber and produce large population redistribution.

A simple dimensionless control suggested by the geometry is the number of
opening opportunities during one collective persistence time.  A school
travels `v tau_p/R` chamber radii during that time, while the open fraction of
the divider span is `w/(2R)`, giving the scaling variable

\[
\Lambda = {v\tau_p\over R}\,{w\over 2R},
\]

where `w` is physical gap width, `R` chamber radius, `v` collective swimming
speed, and `tau_p` the collective directional persistence time.  The
conjecture predicts frozen relaxation for `Lambda << 1`, strongest trapping
around `Lambda=O(1)`, and frequent redistribution for larger `Lambda`.
Varying `R`, `v`, or a microscopically derived `tau_p` should collapse the
crossover when plotted against `Lambda`, rather than against gap ratio alone.

This expression is a first-principles geometric hypothesis, not a derived
rate law: the effective collective persistence and the probability that an
opening encounter completes a coherent relocation remain to be derived.  In
particular, the earlier proposed `Ip_c proportional to In^2` boundary is not
used or treated as established.

## Scope

The result justifies a future predictive theory of bottleneck-controlled
trapping.  It does not justify a stationary potential, an equilibrium phase
transition, or a fitted switching rate.  No new microscopic simulations,
crossing-conditioned diagnostics, torque decompositions, or empirical return
maps were used.
