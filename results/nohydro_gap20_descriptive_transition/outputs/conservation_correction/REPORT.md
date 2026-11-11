# Conservation-law correction to gap-20 plug flow

## Scope and definitions

This is one bounded test of why the parameter-free relation

\[
F_{\rm neck}\simeq {N J_a\over\ell}
\]

fails in the existing 40-seed, no-hydrodynamics gap-20 ensemble.  It uses the
original trajectories without modifying the simulator, adding transport
coefficients, or fitting a potential.

At each measurement window, the source chamber and longitudinal sign are
assigned using only that window's initial state:

- if the left chamber has more fish, $s=+1$;
- if the right chamber has more fish, $s=-1$;
- an exact population tie would use the sign of the instantaneous mean heading.

No tie occurred in the 1,061 analyzed windows.  The eventual crossing direction
is stored only as an audit field and is never used in a current or flux.  The
online sign agrees with that future direction in only 71--79% of the event
windows, independently confirming that future labels were not silently used.

Let $x=s x_{\rm lab}$ increase from the source chamber toward the gap, so the
approach rectangle is

\[
-\ell\leq x<0,\qquad |y|\leq w/2,\qquad
\ell=w=2.8635642126552705.
\]

The longitudinal density and current are the particle measures

\[
n(x,t)=\sum_i I_y(y_i)\,\delta(x-x_i),\qquad
m(x,t)=\sum_i I_y(y_i)v_{x,i}\,\delta(x-x_i),
\]

and $J_a=N^{-1}\int_{-\ell}^{0}m\,dx$.

## Exact finite-window identity

Starting from

\[
\partial_t n+\partial_x m+L=0,
\]

multiplication by $x+\ell$ and integration gives

\[
F_{\rm neck}(t)={N J_a(t)\over\ell}
-{1\over\ell}\int_{-\ell}^{0}(x+\ell)
\left[\partial_t n(x,t)+L(x,t)\right]dx.
\]

The numerical test does not estimate derivatives from spatial bins.  Every
native particle segment is clipped geometrically against the approach
rectangle.  For a window of duration $T$, define

\[
W(t)=\sum_{i\in A(t)}(x_i+\ell)
\]

and let $C_L^W$ be the signed lateral boundary-crossing count weighted by the
crossing position $x+\ell$.  Conservative particle accounting gives

\[
\bar F_{\rm neck}=
{1\over\ell T}\int_t^{t+T}\int_A m\,dx\,dt
-{W(t+T)-W(t)\over\ell T}
-{C_L^W\over\ell T}.
\]

The stored left-endpoint velocity definition of $J_a$ differs slightly from
the segment-clipped path integral when a fish crosses an approach boundary.
Keeping this separately gives

\[
\bar F_{\rm neck}=
{N\bar J_a\over\ell}
-A_W-L_W+\epsilon_{\rm path},
\]

where $A_W=\Delta W/(\ell T)$, $L_W=C_L^W/(\ell T)$, and
$\epsilon_{\rm path}$ is the known path-discretization correction.

Across 1,061 windows of length 0.2, the maximum absolute closure residual is
$6.04\times10^{-14}$ and the RMS residual is $1.09\times10^{-14}$.  The path
correction is typically below $10^{-3}$ fish/time.  Thus the conservation
identity closes to floating-point accuracy, and neither finite spatial bins nor
the stored-velocity convention explains the plug-flow failure.

## Phase-resolved mechanism

Phase labels select intervals retrospectively, but their current direction is
always assigned online as described above.  The pre-crossing interval ends one
native step before the first centerline crossing, so it contains no crossing
segment.  Values below are independent-trajectory seed means with seed-bootstrap
intervals.

| phase | neck flux | plug term | weighted accumulation $A_W$ | weighted lateral flux $L_W$ | total correction |
|---|---:|---:|---:|---:|---:|
| strictly pre-crossing | 0.000 | 1.684 (1.196, 2.226) | 8.975 (6.860, 11.155) | -7.290 (-9.466, -5.234) | 1.684 (1.196, 2.226) |
| sustained crossing | 4.118 (1.654, 6.654) | 2.977 (2.109, 3.859) | -4.275 (-7.384, -1.237) | 3.134 (0.452, 6.009) | -1.141 (-3.040, 0.704) |
| termination/reversal | 0.297 (-0.378, 1.000) | 0.841 (0.389, 1.310) | -7.124 (-9.312, -5.256) | 7.667 (5.687, 9.875) | 0.544 (-0.071, 1.098) |
| matched non-crossing | 0.000 | 0.071 (0.005, 0.145) | -0.091 (-0.738, 0.538) | 0.162 (-0.494, 0.838) | 0.071 (0.005, 0.145) |

The physical interpretation is:

1. **Before crossing, the apparent plug current fills the approach layer.**
   Downstream-weighted population storage is strongly positive in 35 of 37
   event seeds.  Simultaneous lateral inflow supplies much of that buildup, so
   $A_W$ and $L_W$ partly cancel.  Their remaining difference exactly removes
   the positive plug term, giving zero neck flux.
2. **During sustained crossing, depletion augments neck transport.**  Weighted
   storage decreases while lateral loss opposes transport.  On average the
   depletion is larger, so neck flux exceeds the plug estimate.
3. **At termination or reversal, lateral escape wins.**  The layer is still
   depleting, but outward lateral exchange is slightly larger and redirects
   the remaining motion away from the neck.
4. **Controls have neither large storage nor exchange.**  Both terms are small
   and their seed intervals contain zero.

The active-phase contributions are large and reproducible, but strongly
anticorrelated across seeds: the accumulation--lateral seed correlations are
-0.973 before crossing, -0.789 during sustained crossing, and -0.961 during
termination/reversal.  Consequently it is their phase-dependent cancellation,
not either term alone, that determines the net plug correction.  Weighted
population storage has the largest mean absolute contribution before and during
crossing; lateral exchange is slightly larger at termination.

## Held-out conditional closure at fixed $(q,J_a)$

The closure diagnostic uses seeds 0--29 to tabulate the mean total correction
in fixed, prespecified $(|q|,J_a)$ cells.  The use of $|q|$ only applies the
left--right symmetry and adds no state variable.  A cell requires at least
three training seeds.  Seeds 30--39 are held out.

- 308 of 340 held-out windows (90.6%) lie in 21 supported cells.
- Training and held-out cell-mean corrections have Spearman correlation 0.462.
- The plug prediction has seed-mean MAE 1.363; applying the training conditional
  correction gives 1.309.  The improvement is 0.053 with 95% interval
  (-0.155, 0.243), so it is not established.
- Mean bias improves from 0.676 to -0.037, but median within-seed rank
  correlation decreases from 0.429 to 0.387.

The correction therefore calibrates the ensemble mean but does not give a
controlled improvement in individual-window flux prediction.

The reason is explicit.  Across supported cells, train-to-test mean correlation
is only 0.168 for weighted accumulation and 0.196 for lateral exchange.
Held-out residual RMS values are respectively 8.92 and 8.28 fish/time, far
larger than the 2.52 RMS residual of their anticorrelated sum.  Neither large
term is determined by $(q,J_a)$; occasional cancellation only makes the total
look better behaved.

## Answers to the stopping questions

**A. Does the identity close?** Yes.  Conservative segment and boundary
accounting closes at approximately $10^{-14}$ without spatial binning.

**B. Why does plug flow fail?** The approach region is a renewing storage layer,
not a uniform pipe.  Longitudinally weighted population accumulation/depletion
and lateral entry/escape are each much larger than the net correction and
change roles across initiation, propagation, and termination.

**C. Does $(q,J_a)$ support a controlled approximate flux closure?** It supports
mean-bias calibration, but no statistically established held-out MAE or rank
improvement.  This is not a validated predictive closure and is not evidence
of Markovianity.

**D. Which neglected physics prevents closure?** The missing quantity is the
weighted storage-and-side-exchange term

\[
{1\over\ell}\int_{-\ell}^{0}(x+\ell)
[\partial_t n+L]\,dx.
\]

In practical terms, the strongest individual omission is longitudinally
weighted approach-population storage, while lateral exchange is inseparable
because it nearly cancels storage in a phase-dependent way.  Both require the
internal spatial distribution and boundary renewal of the approach layer;
neither can be recovered reproducibly from only $q$ and $J_a$.

## Reproducibility and limitations

- `analyze_conservation_correction.py` implements segment clipping, exact
  boundary accounting, phase selection, seed uncertainty, and held-out tests.
- `conservative_window_budgets.csv` contains all exact window balances.
- `phase_contribution_summary.csv` contains phase-wise seed estimates.
- `conditional_correction_cells.csv` and
  `heldout_conditional_closure.csv` contain the closure audit.
- `summary.json` and `figures/conservation_correction.png` are compact entry
  points.

All 40 trajectories were eligible inputs; 38 contributed at least one selected
event or matched-control interval because two seeds had neither a selected
event nor a selected control.  Phase selection is retrospective and therefore
descriptive, not an online phase classifier.  Controls reuse some raw
trajectories at distinct times, so control uncertainty is clustered by actual
trajectory seed.  No raw trajectory, simulator file, or cache is included.
