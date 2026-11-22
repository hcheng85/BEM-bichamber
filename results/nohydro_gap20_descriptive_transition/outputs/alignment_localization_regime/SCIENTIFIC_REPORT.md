# Alignment-induced macroscopic localization and transport

## Scientific decision

The bounded gate is positive for a **finite-time alignment-induced localization
regime**.  With social interactions removed, and also with alignment set to
zero while positional social interactions remain, population fluctuations
stay near the independent `1/N` scale and never reach `|q|>=0.8`.  Every
tested nonzero alignment strength produces order-one `q^2` and strongly
imbalanced occupation in at least part of the ensemble.  Population mobility
is largest near `Ip=1` and is smaller in the strongly trapped localized
ensembles at larger `Ip`.

This is not yet evidence for a stationary phase transition.  Three-seed scan
intervals are wide, finite-alignment occupation depends strongly on both time
window and initialization, and complete relocations are sparse.  The supported
phenomenon is therefore alignment-generated macroscopic localization with a
mobile onset regime, not a stationary effective potential or a measured
switching process.

This study changes level deliberately.  It does not condition on gap
encounters, decompose torques, fit a return map, or construct an effective
potential.  The observables are the chamber population imbalance

\[
q={N_R-N_L\over N},
\]

its time-occupied distribution, strong-localization occupation, finite-lag
population redistribution, and complete opposite-dominance relocations.

## Exact independent-particle benchmark

At stationary equilibrium in a reflection-symmetric bichamber, each of `N`
independent identical fish has right-chamber probability `1/2`.  Therefore

\[
 N_R\sim {\rm Binomial}(N,1/2),\qquad
 P_0(q)=2^{-N}{N\choose N_R},\quad q={2N_R-N\over N},
\]

and

\[
 \mathbb E_0 q=0,\qquad {\rm Var}_0(q)={1\over N}.
\]

For `N=100`, the exact variance is `0.01` and the exact probability of
`|q|>=0.8` is `3.0633e-17`.  The complete 101-point discrete law is provided
in `independent_binomial_benchmark.csv`; no Gaussian approximation is used.

The social-off simulation is a dynamical control, not the derivation of this
law.  It disables the entire visual social torque, including alignment and
the positional term, while retaining production wall avoidance, angular
noise, speed, clipping, and boundary collision handling.  `Ip=0` is kept as a
separate interacting control because its positional social term remains.

## Existing-data and source audit

The raw-data inventory contains 11 no-hydrodynamics ensembles.  All have
`Ip=4`: gap ratios 10--85 at noise `0.5`, plus gap ratios 10 and 20 at noise
`0.7`.  Only the 40 trajectories with gap ratio 20, `N=100`, noise `0.5`, and
no hydrodynamics match the requested scan.  Each has duration approximately
1000 and uses the production all-left localized initialization.  No compatible
stored trajectories exist at `Ip=0,1,2,8`, with balanced initialization, or
with all social interactions disabled.  The noise-0.7 data are inventoried
but not analyzed because a second control parameter lies beyond this gate.

The scan runner executes an extracted production step without modifying the
source.  It reproduces the archived localized initial state and first state
step to `1.78e-15`, and archived velocity to `4.44e-16`.  The balanced
initialization contains exactly 50 fish on each side.  The social-off torque
is exactly zero.  Fixed parameters are `N=100`, `In=0.5`, `Iw=0.94`, unit
speed, `dt=0.001`, radius ratio 50, and gap ratio 20.

## What the archived `Ip=4` data establish

Over the post-transient interval `t=200--1000`, the 40 archived localized
runs have

- `mean(q)=-0.0781`;
- `mean(q^2)=0.5853` with seed interval `[0.4711,0.6967]`;
- strong-localization occupation `P(|q|>=0.8)=0.5283`
  `[0.3999,0.6582]`;
- nine complete opposite-dominance relocations in six seeds over 32,000 total
  time units;
- 34 trajectories with no observed relocation; and
- median observed dominance-to-opposite-dominance transit duration 97.

This is an enormous departure from independent population fluctuations.  It
does not by itself prove stationary symmetry breaking.  The equal-seed
occupation distributions in `t=200--600` and `t=600--1000` have total-
variation distance 0.232, and complete relocations are sparse.  Although the
late ensemble mean and separate left/right strong-state occupations are
compatible with reflection symmetry, individual seeds remain trapped for long
periods.  The archived histogram is therefore a finite-time occupation law,
not an equilibrium energy landscape or a measured switching-rate process.

## Bounded scan and uncertainty

The new scan uses three independent seeds for each initialization and
condition.  It runs to `T=300`, saves only `q(t)` every 0.1, and compares the
common late window `t=200--300`.  The conditions are social-off and
`Ip=0,1,2,4,8`; archived localized `Ip=4` runs are reused, so only three new
balanced `Ip=4` trajectories are generated.  In total there are 33 new compact
runs.

Localized starts use the production all-left distribution.  Balanced starts
use the same production chamber distribution with exactly 50 fish on each
side and random headings.  Balanced is a neutral macroscopic initial
condition, not an equilibrium sample.

Uncertainty resamples complete trajectory seeds, never saved time points.
Signed and reflection-symmetrized distribution distances are both reported.
Complete relocation means first arrival in the opposite `|q|>=0.6` dominance
basin after the last sample in the source basin.  Every finite trajectory has
a right-censored terminal interval.  No transition rate is inferred.

## Alignment scan results

In the common `t=200--300` window, the independent social-off control gives
`mean(q^2)=0.00586` (`95%` seed bootstrap `[0.00129,0.01123]`) from localized
starts and `0.00572` (`[0.00238,0.00973]`) from balanced starts.  `Ip=0`
gives `0.01515` (`[0.00827,0.02731]`) and `0.00648`
(`[0.00425,0.00966]`), respectively.  None of these 12 runs ever occupies
`|q|>=0.8` in the comparison window.  Thus the positional social interaction
left at `Ip=0` does not produce the macroscopic localization seen with
alignment.

At `Ip=1`, `mean(q^2)` rises to `0.538` (`[0.293,0.729]`) for localized
starts and `0.317` (`[0.035,0.720]`) for balanced starts; strong-state
occupations are `0.391` and `0.248`.  At `Ip=2,4,8`, the point estimates of
`mean(q^2)` remain between `0.368` and `0.816`, compared with the exact
independent value `0.01`.  The intervals are broad because each new condition
has only three independent trajectories.  The 40-run archived `Ip=4`
localized ensemble is better resolved: over the same window,
`mean(q^2)=0.6107` (`[0.4785,0.7331]`) and strong-state occupation is
`0.5468` (`[0.4034,0.6881]`).

These results locate the onset within the deliberately coarse interval
`0<Ip<=1`; they do not resolve a critical value or the order of a transition.

## Initialization and stationarity

The independent control loses its initial localization: its localized-minus-
balanced `q^2` difference is `0.00014` (`[-0.00562,0.00590]`).  Finite
alignment does not generally forget initialization by `T=300`.  For example,
at `Ip=8` the difference is `0.428` (`[0.0737,0.6897]`), and the
reflection-symmetrized occupation distributions have total-variation distance
`0.640`.  At `Ip=1--4` the three-seed difference intervals include zero, but
symmetrized distribution distances remain `0.355--0.615`.

Adjacent-window symmetrized total-variation distances over `t=100--200`
versus `200--300` are `0.36--0.57` for most finite-alignment pilot groups.
Even the 40 archived `Ip=4` trajectories give distance `0.156` between
`t=200--600` and `600--1000`.  The scan occupation curves are consequently
reported as finite-time, equal-seed measures.  Balanced initialization is a
neutral macroscopic start, not a draw from a stationary distribution.

## Localization versus transport

The ten-time-unit mobility `mean(|q(t+10)-q(t)|)` is maximal at `Ip=1` in
both initialization groups: `0.131` (`[0.105,0.171]`) from localized starts
and `0.134` (`[0.092,0.161]`) from balanced starts.  The corresponding
social-off values are `0.0456` and `0.0610`.  For localized starts mobility
then falls to `0.0825`, `0.0330`, and `0.0209` at `Ip=2,4,8`, respectively,
while strong localization grows.  This supports a nonmonotonic competition
between alignment-enabled collective transport and alignment-induced
trapping.  The balanced `Ip=8` estimate rebounds to `0.107`; together with
its initialization sensitivity, this prevents claiming a universal
steady-state transport curve.

Only one complete opposite-dominance relocation occurs in all common
`t=200--300` scan windows (in the 40-run archived `Ip=4` group).  The longer
archived post-transient interval contains nine events in six seeds, with 34
of 40 seeds showing none.  Relocations seen before `t=200` in the short pilot
are initialization relaxation and are not used as stationary events.  These
counts are sufficient to document mobility and censoring, but not to estimate
an alignment-dependent switching rate.

## Physical interpretation and one next prediction

Alignment changes the population fluctuations from the `O(N^{-1/2})`
independent scale to chamber-scale collective occupation.  Near onset the
school can move fish coherently between chambers, producing the mobility peak;
stronger alignment makes the population state more persistent and suppresses
redistribution from localized starts.  This is the most compelling
macroscopic phenomenon in the present data.  Positional interactions alone
do not reproduce it.

One parameter-free regime-boundary scaling follows from the microscopic
angular law.  The orientation noise is `dtheta=In dW`, so its rotational
diffusivity is `D_r=In^2/2`, while the linearized alignment relaxation is
proportional to `Ip` times a dimensionless directed-network relaxation factor.
At fixed density and geometry, balancing these terms predicts

\[
 {Ip_c(In)\over Ip_c(0.5)} \simeq \left({In\over0.5}\right)^2.
\]

This is a testable boundary prediction, not a fitted critical law: a bounded
repeat near the onset at one second noise amplitude should shift the
localization boundary according to this ratio.  Stationary sampling from both
initializations must be demonstrated before fitting any `q` landscape.

## Scope and reproduction

The empirical `P(q)` curves are time-occupation distributions.  No scalar
potential, barrier, diffusion process, or switching rate is fitted.  Short
fish crossings and approach episodes are not used as substitutes for
collective relocation.

```bash
uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/run_alignment_localization_scan.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/run_alignment_localization_scan.py --stage all --workers 9

uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/analyze_alignment_localization_regime.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/analyze_alignment_localization_regime.py
```

Raw microscopic trajectories and resumable per-run caches are excluded from
version control.  The committed simulation product contains only compact
population-count trajectories and provenance metadata.
