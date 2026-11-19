# Population-lattice kinetics and the landscape decision

## Decision

The polarization and approach-storage variables **cannot be eliminated to
obtain a useful one-dimensional Markov process for population imbalance**.
The supported coarse description is a reflection-symmetric hybrid
nonequilibrium process with

```text
n_R in {0,1,...,100}       discrete population coordinate,
p_x                         global signed axial polarization,
b_L, b_R                    chamber-resolved approach storage.
```

Within one directed transfer channel, the minimal tested jump intensity is
`lambda_+(n_R,p_x,b_L)` for left-to-right transfer and its reflected partner
`lambda_-(n_R,p_x,b_R)`.  Source-side storage carries the largest additional
rate information, while polarization contributes a smaller but independently
resolved increment.  Target-side storage does not improve the held-out
instantaneous crossing intensity after those variables are known.  Both
`b_L` and `b_R` remain in the global state because either chamber can be the
source and the majority label is singular at `q=0`.

A scalar birth--death potential is therefore not reported.  Computing one
from marginal `q` rates would force detailed balance onto dynamics that retain
resolved hidden-state dependence and reproducible projected probability
currents.

## Population record and exact state convention

The coordinate is `q=(N_R-N_L)/100=2n_R/100-1`, so one fish crossing changes
`q` by exactly `0.02`.  All 40 archived no-hydrodynamics, gap-20 trajectories
were analyzed after time 200.  No trajectory or microscopic simulation was
created.

The common population record has spacing 0.1 and contains 3,269
left-to-right and 3,116 right-to-left fish transfers (6,385 total).  Native
0.001 occupancy records exist for 20 seeds.  On those seeds the common record
retains 3,145 of 3,211 native crossings, or 97.95%.  Thus the rate calculation
uses a uniform 40-seed resolution while its small rapid-recrossing loss is
explicitly measured.

Reflection pooling uses

```text
(q,p_x,b_L,b_R, left-to-right)
    -> (-q,-p_x,b_R,b_L, right-to-left).
```

The tabulated rates use all partial fish-transfer events, not the nine complete
opposite-basin relocations.  At `q=0`, both marginal directional intensities
are 0.206 per time unit.  The exact-state rate curve is visibly rough away
from well-sampled states; no smoothing function or arbitrary potential was
fitted.

## Why a one-dimensional q process fails

### The nominal fast variables do not relax before population events

Microscopic configurations were sampled every 0.5 time unit.  After subtracting
the exact-`q` conditional mean and retaining only intervals with no population
crossing, the seed-mean conditional autocorrelations were:

| lag | `p_x` | `b_L` | `b_R` |
|---:|---:|---:|---:|
| 0.5 | 0.943 | 0.891 | 0.919 |
| 1.0 | 0.821 | 0.711 | 0.761 |
| 2.0 | 0.500 | 0.307 | 0.347 |
| 5.0 | -0.314 | -0.131 | -0.091 |

The negative long-lag correlations show oscillatory collective dynamics, not
instantaneous equilibration.  The variables remain strongly persistent over
the 0.5--2 time-unit interval in which partial transfers accumulate.

### A q-only kernel violates Chapman--Kolmogorov beyond sampling error

Transition kernels were learned from seeds 0--29 on the exact 101-state
population lattice.  For every tested lag, the weighted total-variation
difference between the directly measured `P(2 tau)` and `P(tau)^2` exceeds a
multinomial sampling-null envelope:

| `tau` | CK residual | sampling-null 95% interval | held-out one-lag TV |
|---:|---:|---:|---:|
| 0.1 | 0.0056 | 0.0025--0.0030 | 0.0146 |
| 0.5 | 0.0235 | 0.0058--0.0065 | 0.0417 |
| 1.0 | 0.0333 | 0.0080--0.0087 | 0.0588 |
| 2.0 | 0.0403 | 0.0103--0.0111 | 0.0817 |
| 5.0 | 0.0424 | 0.0138--0.0149 | 0.1477 |

At short lag the absolute residual is small but statistically resolved; at
long lag, where polarization and storage have crossed through zero
correlation, seed-transferable prediction becomes worse rather than
Markovian.  There is no demonstrated elimination window.

Recent population history also improves held-out prediction at the same
current `q`.  Its reduction of future-increment MSE is positive across the ten
test seeds at lags 0.1, 0.5, and 1.0; at lag 1 the improvement is
`2.31e-5`, with 95% seed interval `[1.18e-5,3.50e-5]`.  It is unresolved at
lag 2.  The history table overfits at lag 5, which does not rescue the q-only
model because the independently measured kernel and CK errors are then
largest.

## Which additional variables are necessary?

Crossing counts in each 0.5-time window were treated as directional count
data.  Fixed physical bins were populated from seeds 0--29; Poisson deviance
was evaluated on seeds 30--39.  These are out-of-seed predictions, not AUC
comparisons.

| candidate intensity state | held-out deviance | improvement over q-only, 95% seed interval |
|---|---:|---:|
| `q` | 0.5855 | -- |
| `(q,p_x)` | 0.5318 | 0.0537 `[0.0398,0.0712]` |
| `(q,b_source)` | 0.3191 | 0.2656 `[0.1972,0.3353]` |

On the common support, adding `p_x` to `(q,b_source)` gives a further deviance
reduction of 0.0204, with 95% seed interval `[0.0144,0.0261]`.  Adding target
storage after `(q,p_x,b_source)` changes deviance by -0.0035, interval
`[-0.0115,0.0043]`, and is unresolved.

The physical interpretation is direct.  `b_source` measures how much
population is already stored near the neck on the supplying side, while
`p_x` distinguishes whether collective axial motion is carrying that stored
population toward or away from the opening.  Population imbalance alone does
not specify either condition.

## Nonequilibrium landscape

The stationary occupancy from the preceding state-space audit remains a valid
statistical landscape, but it is not a dynamical potential.  At lag 0.5 the
antisymmetric transition flux in the `(q,p_x)` projection has an irreversible
edge fraction of 0.0205.  Its independently calculated train/test edge-current
cosine is 0.670, so the direction of the projected current is reproducible
across original seeds.  This is evidence for circulating nonequilibrium flux;
it is not a claim that the two-dimensional projection is closed.

Accordingly, a useful landscape representation must contain both stationary
mass and probability current in `(q,p_x,b_L,b_R)`.  A scalar `-log pi(q)` may
still be plotted as marginal occupancy, but it cannot be interpreted as the
generator's potential or used to infer a barrier and Kramers rate.

## Minimal stochastic theory to develop

The data support a hybrid generator of the form

```text
L f = L_cont(p_x,b_L,b_R | n_R) f
    + lambda_+(n_R,p_x,b_L) [f(n_R+1)-f(n_R)]
    + lambda_-(n_R,p_x,b_R) [f(n_R-1)-f(n_R)].
```

Reflection fixes `lambda_-(n,p_x,b_R)` from the reflected `lambda_+` law.
The continuous part must be constrained by the already derived polarization
evolution equation and exact approach-storage/neck-flux conservation laws.
It must retain circulating current; it should not be replaced by a fitted
gradient drift.

What has genuinely been predicted beyond the rate inputs is the crossing-count
distribution on ten independent seeds: source storage and polarization both
improve it with seed-bootstrap intervals excluding zero.  What has not yet
been predicted is the stationary joint density, the full population
relaxation spectrum, or the complete relocation rate from a closed
multidimensional generator.  Estimating the continuous conditional evolution
and then rerunning CK tests in this fixed physical state space is the next
first-principles step.  Complete relocations remain a final limited-sample
check, not the rate-estimation data set.

## Reproducibility and scope

`analyze_population_jump_kinetics.py` performs the reflection-symmetric rate
estimate, resolution audit, fixed-`q` relaxation test, CK calibration, history
test, held-out conditioned-intensity comparison, and projected-current audit.
The feature cache is excluded from version control.  No production source,
trajectory, crossing definition, empirical order parameter, fitted potential,
or microscopic simulation was introduced.
