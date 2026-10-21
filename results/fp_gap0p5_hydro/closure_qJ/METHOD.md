# Method: empirical closure in chamber occupancy and neck current

## Scope and data

This extension uses the same 20 full-BEM hydrodynamic trajectories as the
one-dimensional audit. The archive label is `R50/gap50`; the stored geometry is
`gap/R=1.0`, equivalently `gap/(2R)=0.5`. It is not described here as literal
`gap/R=0.5`.

Seeds 0--14 are used for estimation and seeds 15--19 only for validation. The
first 200 time units are removed. Compact features are sampled every 0.01 time
unit. No trajectories are concatenated when increments, autocorrelations,
passages, or realization-level uncertainty are calculated.

The primary continuous occupancy is

\[
q_\ell=N^{-1}\sum_i\tanh(x_i/\ell),\qquad \ell=0.5.
\]

The current is calculated from the archived instantaneous velocity field,

\[
J=N^{-1}\sum_i (v_{x,i}/\ell)\operatorname{sech}^2(x_i/\ell).
\]

Thus no future occupancy is used to construct `J`.

## Transition-kernel discretization

The occupancy edges are

`[-1, -.9, -.75, -.5, -.4, -.2, 0, .2, .4, .5, .75, .9, 1]`.

They retain separate well, threshold, and switching-corridor regions. For the
joint kernel, each occupancy bin is divided into three current bins using
training-set conditional current tertiles. This gives 12 occupancy states and
36 joint states. All 12/12 occupancy states and 36/36 joint states meet the
pooled support threshold. In the switching corridor, all 4/4 occupancy origins
and 12/12 joint origins are retained; no sparse state is removed to improve the
score.

For each lag `tau = 0.1, 0.2, 0.5, 1.0`, transition counts are accumulated
within each realization. The Chapman--Kolmogorov comparison is between the
direct `P(2 tau)` and `P(tau)^2`. Two errors are retained:

1. Full-state total variation in the native state space.
2. Joint-kernel target states collapsed over final current, giving a projected
   occupancy error directly comparable with the occupancy-only model.

Origin rows require at least 500 training samples. The sampling/distribution
reference is the total-variation difference between independently estimated
training and held-out `P(2 tau)` kernels. Confidence intervals resample whole
training realizations 500 times.

The operational coarse-closure check requires, at both lags 0.5 and 1.0:

- smaller projected-occupancy CK error than the occupancy-only model;
- projected error no more than 15% above the train--test kernel difference;
- full joint-state error no larger than its full-state train--test difference.

This is a practical finite-resolution criterion, not proof of exact
Markovianity.

## Local prediction and current memory

Held-out future occupancy is predicted by:

1. the training conditional mean given occupancy;
2. the ballistic formula `q(t)+tau J(t)`;
3. the training conditional mean given occupancy and current.

Prediction uses the same 12 occupancy bins and 15 conditional current bins per
occupancy bin. All 180 joint prediction bins are supported, with minimum
training support above 1,000 samples. RMSE and MAE are first computed per held-
out realization and then bootstrapped over the five realization-level values.

Current autocorrelation is calculated after subtracting the training
occupancy-conditioned mean current. The reported e-folding and positive
integral times therefore measure residual directional persistence rather than
the trivial change of mean current across occupancy states.

## Splitting probability

A corridor excursion begins when a trajectory enters `|q_l| <= 0.2` from
outside and ends at the first hit of either dominance threshold `q_l = +/-0.4`.
States inside completed excursions are sampled every 0.1 time unit and inherit
the eventual threshold label. Final incomplete excursions are counted as
censored and are not assigned an outcome.

The committor grid has eight occupancy bins and five training-current quantile
bins. Beta(1,1) smoothing is used for empirical cell probabilities. Held-out
Brier score and log loss are calculated per seed and bootstrapped over seeds.

## First passage and censoring

Because the joint state passes the declared coarse-lag check, finite-state
kernels at lag 0.5 are tested for switching. The opposite threshold is made
absorbing, and the transient kernel predicts survival and MFPT from each actual
held-out initial state. The reported curve averages over those held-out starts.

Observed survival uses the Kaplan--Meier estimator. The final unfinished
residence interval in every held-out realization is retained as right-censored;
there are three censored LR episodes and two censored RL episodes. Observed
restricted means and model MFPT uncertainty are bootstrapped by resampling
whole validation and training realizations, respectively.

## Stationary-distribution correction

Stationary probabilities are obtained from training transition kernels, but
their JS divergences are evaluated against occupancy counts from test seeds
15--19 only (400,005 samples). No training histogram is labeled held out.

## Interpretation limits

Finite binning can hide memory and the train--test kernel difference includes
both sampling error and realization heterogeneity. A successful lag-0.5 kernel
is therefore evidence for a useful coarse Markov model, not a continuous
two-dimensional Fokker--Planck limit or a first-principles derivation. No 2-D
PDE, additional state variable, neural predictor, or effective potential is
introduced in this focused test.
