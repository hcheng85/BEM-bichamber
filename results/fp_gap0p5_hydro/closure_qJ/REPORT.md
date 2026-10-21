# Two-variable closure for hydrodynamic bichamber switching

## Conclusion: B. Partially successful closure

Adding the instantaneous neck current changes the result substantially. The
joint state `(q_l,J)` provides a useful **coarse-time Markov model** around lags
0.5--1.0 and predicts held-out chamber-switching times far better than occupancy
alone. It is not a scale-independent closed stochastic dynamics: the joint
Chapman--Kolmogorov test fails at lags 0.1--0.2, and a ballistic predictor
outperforms the binned joint conditional predictor at every tested horizon.

The physical result is therefore specific: directed transport through the neck
is the principal missing commitment variable in the one-dimensional model, but
the evolution of that current still contains unresolved acceleration,
orientation, spatial-organization, and hydrodynamic memory.

The analyzed archive is `R50/gap50`, with **gap/R = 1.0** and
**gap/(2R) = 0.5**. Twenty trajectories are used: 15 training and five held out,
with time 0--200 removed.

## Quantitative summary

| Test | Occupancy only | Occupancy + current | Interpretation |
|---|---:|---:|---|
| Corridor CK error, lag 0.1 | 0.0239 | 0.0387 | current state is worse at micro-lag |
| Corridor CK error, lag 0.2 | 0.0710 | 0.1047 | unresolved current acceleration |
| Corridor CK error, lag 0.5 | 0.1787 | 0.0595 | 66.7% reduction |
| Corridor CK error, lag 1.0 | 0.2646 | 0.1181 | 55.4% reduction |
| Held-out splitting Brier score | 0.2592 | 0.1785 | 31.2% improvement |
| Held-out splitting log loss | 0.7119 | 0.5168 | 27.4% improvement |
| LR first-passage prediction | 124.3 | 38.72 | observed restricted mean 36.27 |
| RL first-passage prediction | 145.9 | 38.26 | observed restricted mean 38.22 |
| LR survival maximum error | 0.411 | 0.172 | large improvement, residual remains |
| RL survival maximum error | 0.424 | 0.123 | large improvement, residual remains |

All 36 joint kernel states are supported. The switching-corridor comparison
uses all 12 joint origin states and about 144,000 origins at each tested lag.
At lag 0.5, the full joint CK error is 0.0745, below the independent train--test
kernel difference 0.0939. At lag 1.0 these values are 0.1337 and 0.1458. Thus
the 0.5--1.0 range is consistent with an approximate finite-resolution closure,
although the residual is not negligible.

## What the current explains

Inside `|q_l| <= 0.2`, the lowest current quintile reaches the right threshold
with probability 0.020 in training and 0.040 held out. The highest current
quintile has probabilities 0.947 and 0.983. This separation occurs across the
entire central occupancy range: the sign and magnitude of current determine
commitment much more strongly than occupancy itself.

The splitting surface generalizes to the five held-out realizations. Using
current reduces Brier score from 0.2592 (95% seed-bootstrap interval
0.2369--0.2716) to 0.1785 (0.1302--0.2338). It reduces log loss from 0.7119 to
0.5168. There are 353 completed training and 120 completed held-out corridor
excursions; one final excursion in each split is censored rather than assigned
an outcome.

## Why this is not merely a better short-time predictor

Current is strongly persistent. After removing its occupancy-conditioned mean,
its e-folding time is 3.0 and its positive integral correlation time is 2.54.
The median time from central-corridor entry to a dominance threshold is only
1.46 in training and 1.50 held out. Current therefore persists through the
physical commitment stage. It is not just an instantaneous indicator measured
after crossing.

However, the elementary ballistic relation explains most local accuracy. In
the central corridor at horizon 1.0, held-out RMSE is 0.2614 for occupancy only,
0.0768 for the binned `(q_l,J)` conditional predictor, and 0.0525 for
`q_l+tau J`. At horizons 0.1--0.5 the ballistic advantage is much larger. This
means the local predictive benefit of `J` is primarily persistent directed
motion, not evidence that all forces controlling the future of `J` have been
resolved.

## Switching validation with proper censoring

The lag-0.5 finite-state joint kernel gives LR MFPT 38.72
(training-seed bootstrap 35.17--42.45), versus held-out Kaplan--Meier restricted
mean 36.27 (validation-seed bootstrap 30.09--42.98). For RL it gives 38.26
(35.19--42.47), versus 38.22 (31.66--44.38). The occupancy-only kernel predicts
124.3 and 145.9.

The joint survival curves are not exact: maximum errors are 0.172 and 0.123.
But they reproduce both the transition scale and most of the distribution,
whereas occupancy alone leaves a very long artificial tail. This is genuine
held-out improvement rather than fitting the stationary histogram.

The five final incomplete validation residences are included as right-censored
observations: three LR and two RL. They are not discarded or counted as
completed passages.

## Correction of the two earlier validation issues

- Stationary comparisons now use training kernels against **test-seed** counts
  only. Held-out JS divergence is 0.000894 for the occupancy kernel and
  0.000599 for the joint kernel, evaluated on 400,005 test samples.
- First-passage survival and restricted means use Kaplan--Meier censoring, and
  uncertainty resamples whole realizations.

## Physical interpretation

1. **At fixed occupancy, current determines commitment strongly.** Opposite
   current quantiles change right-commitment probability from nearly zero to
   nearly one.
2. **Current memory is shorter than full residence time but comparable to the
   crossing-decision time.** Its e-folding time is 3.0, versus residence means
   near 36--38 and a median corridor commitment time near 1.5.
3. **Current explains the lag dependence of the 1-D reduction.** Projecting out
   `J` mixes configurations moving rapidly in opposite directions. At short
   lag they have distinct deterministic increments; at coarse lag their
   persistence is averaged into lag-dependent drift and diffusion.
4. **The evidence supports an effective persistent-current description.** It
   resembles an underdamped coordinate--velocity state more than overdamped
   barrier diffusion. This is an empirical analogy, not a derivation of mass or
   mechanical energy.
5. **The joint state is still incomplete at micro-lags.** The future change of
   current depends on fish orientations, neck crowding and arrangement,
   wall-induced turning, and hydrodynamic configuration. These variables are
   the likely source of the remaining short-lag CK failure and survival error.

## Bottom line

Directed collective transport through the neck is the main missing ingredient
in the one-dimensional potential picture. Adding it converts a poor
first-passage model into a quantitatively useful coarse Markov model. It does
not justify a continuous two-dimensional Fokker--Planck PDE or establish exact
Markov closure, so this investigation stops at conclusion B rather than adding
more phenomenological state variables.
