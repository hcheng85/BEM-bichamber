# Conservation-constrained hybrid model for chamber transfer

## Decision

The state `X=(n_R,p_x,b_L,b_R)` supports an exact and physically interpretable
jump/conservation structure, but it does **not** support a closed Markov model
at the available resolution.  The proposed hybrid generator fails its
between-crossing memory and multi-lag Chapman--Kolmogorov tests.  It predicts
the coarse population component of partial-transfer pathways reasonably well,
but not their joint polarization/storage evolution or independent-seed
stationary occupancy.

Consequently there is not yet a dynamically validated collective landscape.
No scalar or multidimensional effective potential, stationary current field,
or whole-school switching rate is constructed from this failed closure.

The specific missing physics is the **transport phase of each approach
layer**: the same values of `(p_x,b_L,b_R)` can be filling or draining because
longitudinal supply and weighted lateral exchange retain their own temporal
persistence.  This diagnosis follows from the conservation law and held-out
memory; it is not a proposal to add another empirical order parameter.

## Exact coordinates and crossing maps

With approach length `ell` and `N=100`, the chamber storages are

```text
b_L = (1/N) sum I[-ell <= x_i < 0, |y_i| <= ell/2] (1+x_i/ell),
b_R = (1/N) sum I[0 <= x_i <= ell, |y_i| <= ell/2] (1-x_i/ell).
```

The weight is continuous and equal to one at the center opening.  Therefore a
single left-to-right crossing has the exact instantaneous map

```text
J_+(n_R,p_x,b_L,b_R)
    = (n_R+1, p_x, b_L-1/N, b_R+1/N),
```

and a right-to-left crossing has

```text
J_-(n_R,p_x,b_L,b_R)
    = (n_R-1, p_x, b_L+1/N, b_R-1/N).
```

Global `p_x=N^{-1} sum cos(theta_i)` has no membership jump because it includes
all fish on both sides.  The crossing transfers storage between chamber
accounts without changing `b_L+b_R` instantaneously.

For a finite native segment `x_0<0<x_1`, the crossing fish contributes

```text
Delta b_L = -1/N - x_0/(N ell),
Delta b_R = +1/N - x_1/(N ell).
```

The `+-1/N` terms are the centerline jump.  The remaining terms are continuous
travel from `x_0` to the opening and from the opening to `x_1`.  The analogous
right-to-left formulas reverse the jumps and directions.

This was checked for all 3,211 native crossings in the 20 available processed
occupancy files.  Nineteen seeds contain crossings.  There are no population
count mismatches, no crossings outside the physical opening, and the maximum
finite-segment jump-map residual is `1.60e-18`.

## Conservation between centerline jumps

Let `A_L` and `A_R` be the two approach rectangles and let `dL_L,dL_R` denote
weighted outward lateral exchange.  The exact piecewise-linear balances are

```text
db_L = [sum_(i in A_L) v_xi /(N ell)] dt - dL_L
       - dC_LR/N + dC_RL/N,

db_R = [-sum_(i in A_R) v_xi /(N ell)] dt - dL_R
       + dC_LR/N - dC_RL/N.
```

The upstream boundary has zero weight, so it contributes no storage jump.
Lateral entries and exits are discrete at microscopic resolution; at the
collective level they are part of the stochastic between-centerline evolution,
not part of `J_+` or `J_-`.

Two hundred half-time-unit windows, five from every original seed, were
evaluated for both chambers using clipped native trajectory segments.  All 400
budgets satisfy

```text
Delta b_c = longitudinal path contribution
            - weighted lateral outflux
            - net centerline outflux/N
```

with maximum residual `5.90e-17` and RMS residual `6.64e-18`.

## Candidate hybrid generator

Writing `z=(p_x,b_L,b_R)`, the conservation-consistent candidate is

```text
L f(X) = mu(X) . grad_z f
       + 1/2 C(X) : grad_z grad_z f
       + lambda_+(X) [f(J_+X)-f(X)]
       + lambda_-(X) [f(J_-X)-f(X)].
```

The storage components of `mu` are estimated only after subtracting the exact
centerline jumps above; they therefore contain longitudinal accumulation and
lateral exchange, as required by the conservation identity.  `C` is the
conditional covariance rate of those de-jumped increments.  Reflection acts
as

```text
R(n_R,p_x,b_L,b_R)=(N-n_R,-p_x,b_R,b_L),
```

and constrains `lambda_+(X)=lambda_-(R X)` and the corresponding drift and
covariance components.

To avoid a sparse four-dimensional table, coefficients and closure were
tested on only 60 physical cells:

- the five previously established population regimes;
- three axial regimes, `p_x<-0.25`, neutral, and `p_x>0.25`;
- low/high storage in each chamber at the established threshold `0.075`.

Only 28 cells have at least 100 reflection-pooled training increments, yet
they cover 99.04% of held-out configurations.  The population coordinate
remains the exact 101-state lattice in the jump maps; the five population bins
are only the supported projection used to test the continuous evolution.

## Are the 0.5-window crossing rates instantaneous intensities?

Not exactly.  Rates were re-estimated from common origins at horizons 0.1,
0.2, and 0.5 in 20 supported reflection-pooled cells of oriented population,
polarization, and source storage.

| horizon | conditional rate change from 0.1 | event-probability MAE | mean within-cell Fano factor |
|---:|---:|---:|---:|
| 0.1 | 0 | 0.0105 | 1.32 |
| 0.2 | 6.88% | 0.0269 | 1.56 |
| 0.5 | 11.56% | 0.0640 | 2.28 |

The seed-mean *total* rate changes are unresolved: at 0.5 the held-out
difference from the 0.1 estimate is `0.0021` fish/time, with 95% seed interval
`[-0.0105,0.0133]`.  Thus mean counts scale approximately with duration, but
the increasingly overdispersed counts and changing conditional rates reject a
frozen-rate Poisson interpretation of the 0.5-window estimate.  The 0.5 value
is an integrated conditional hazard.  The 0.1 estimate is the closest
available approximation to an instantaneous intensity and must be updated
after state changes.

## Between-crossing evolution retains memory

Training seeds 0--29 supply conditional means and covariance rates of the
de-jumped `(Delta p_x,Delta b_L,Delta b_R)`.  On seeds 30--39, residuals remain
strongly correlated with the preceding increment at the same current state:

| residual component | mean seed Spearman correlation | 95% seed interval |
|---|---:|---:|
| `p_x` | 0.475 | `[0.441,0.508]` |
| `b_L` | 0.551 | `[0.530,0.571]` |
| `b_R` | 0.466 | `[0.420,0.507]` |

All three intervals exclude zero by a wide margin.  The current state records
how much population is stored, but not whether longitudinal current and
lateral exchange are presently loading or unloading it.  Likewise `p_x`
retains turning/oscillation phase.  A white-noise diffusion for `z` therefore
does not close the generator.

## Joint operator and partial-transfer predictions

The reflection-symmetric empirical 0.5 operator on the 60 cells is the least
parametric Markov approximation to this state.  If even this saturated
operator fails, a more restrictive drift/diffusion parameterization on the
same state cannot repair closure.

| lag | Chapman--Kolmogorov TV | sampling-null 95% interval | held-out TV |
|---:|---:|---:|---:|
| 0.5 | 0 by construction | `0.0056--0.0081` | 0.0186 |
| 1 | 0.0248 | `0.0080--0.0112` | 0.0395 |
| 2 | 0.1101 | `0.0103--0.0139` | 0.0959 |
| 5 | 0.2270 | `0.0125--0.0158` | 0.1520 |

The CK residual exceeds sampling error from lag 1 onward and grows rapidly.
This is consistent with the direct residual-memory test.

Partial population transfers provide 1,042 independent-test pathway origins.
The model predicts their coarse `q`-regime distribution with TV error
`0.022--0.046` over the next 0.5--5 time units.  It does not predict the joint
state pathway: joint TV error is 0.102 at 0.5, 0.166 at 1, 0.234 at 2, and
0.212 at 5.  The useful prediction is therefore limited to population
redistribution; collective organization around that redistribution remains
descriptive.

## Stationary occupancy and landscape gate

The training operator has a numerically converged invariant distribution, but
it does not transfer to independent seeds.  Its five-bin `q` marginal is

```text
[0.321, 0.089, 0.180, 0.089, 0.321],
```

whereas the held-out occupancy is

```text
[0.178, 0.163, 0.319, 0.163, 0.178].
```

The q-marginal TV difference is 0.286 and the full joint difference is 0.287.
This large seed dependence reflects the slow, poorly mixed dominance basins;
it cannot be repaired by interpreting the training invariant distribution as
an energy landscape.

Because both the Markov and independent-stationarity gates fail, no
nonequilibrium landscape or probability-current interpretation is claimed.
The empirical occupancy can still be plotted descriptively, but `-log p` is
not a validated quasipotential and the antisymmetric flux of this unclosed
projection is not a closed collective current.

## What is established and what follows

Established exactly:

- the `+-1` population and `+-1/N` chamber-storage crossing maps;
- zero crossing jump in global `p_x`;
- exact chamber storage conservation with longitudinal, lateral, and
  centerline terms.

Predicted on independent trajectories:

- approximate short-horizon crossing count rates;
- the q-regime portion of frequent partial-transfer pathways.

Not dynamically predicted:

- polarization/storage pathways beyond one step;
- independent-seed stationary occupancy;
- a multidimensional probability landscape or current;
- the complete opposite-chamber relocation rate.

The nine complete relocations were not used to establish or tune the model.
Since the much larger partial-transfer cohort already rejects closure, those
nine events cannot validate a whole-school switching theory.

The next theoretical issue is now precise: derive a non-Markovian or
colored-transport representation of the existing approach-current and lateral
exchange terms, or demonstrate a controlled time-scale limit in which their
phase can be eliminated.  This report does not expand the state space or add a
phenomenological coordinate; it identifies why the requested four-variable
state is not yet a dynamically validated collective landscape.
