# R50/gap50 hydrodynamic one-dimensional FP validation

## Geometry normalization warning

The requested case was called `gap/R=0.5`, but the archived `R50/gap50` files
store physical `gap=7.15891` and `R=7.15891`.  Thus this ensemble has
**gap/R = 1.000** and **gap/(2R) = 0.500**.  No `R50/gap25`
hydrodynamic ensemble was found.  The results below answer the question for the
available gap-50 ensemble, which matches 0.5 only when normalized by chamber
diameter.

## Dataset and split

- Full BEM hydrodynamics: `If=0.01`, direct and boundary-mediated calls source-verified.
- Geometry: `R=7.15891`, `gap=7.15891`, gap/R = 1.000, gap/(2R) = 0.500.
- N=100, Ip=4.0, In=0.5, Iw=0.94, dt=0.001, T=1000.0.
- Twenty independent seeds; seeds 0--14 train, 15--19 held out; first 200 time units removed.
- Complete dominance switches at |q|=0.4: 316 train and 108 held out.
- The continuous coordinate uses ell=0.5, the rounded divider-tip radius; ell=0.25 and 1.0 are explicit sensitivity controls.

## Transition character

- Post-transient archived cluster events: 689.
- Absolute cluster mark: mean 62.478 fish, median 78.0, maximum 100.
- A nonzero native one-step occupancy change has mean magnitude 0.02014 in q and maximum 0.06000.
- Lag-one correlation of within-realization event waiting intervals: 0.111.

## One-dimensional discrete-q FP diagnostic

- Supported fitting-lag interval: None.  None means no interval met simultaneous increment-shape and diffusion-stability requirements.
- Reference maximum lag: 0.5, the first tested near-Gaussian increment scale; it is a coarse diagnostic, not a tau-to-zero estimate.
- Dynamic barrier at the reference lag: 4.099; barrier max/min ratio across tested lags: 1.32.
- Histogram barrier: 4.998.  Its existence confirms bistability, but not dynamical closure.
- Stationary JS divergence: 0.03011; stationary-current relative RMS residual: 0.8666.
- CK weighted-TV range: 0.002483--0.1381.
- Held-out LR MFPT: observed 34.207 [28.915, 39.302], predicted 65.822.
- Held-out RL MFPT: observed 37.548 [31.616, 44.178], predicted 77.410.

The lag-1.0 fit gives LR/RL MFPTs 33.017/35.503, close to the
held-out observations.  This agreement is post-hoc: at that lag the median
diffusion has changed by a factor 1.807 across the near-Gaussian lag range and the CK
error is large.  It therefore validates a coarse-time empirical propagator, not
a unique short-time diffusion coefficient.  Enforcing exact left-right symmetry
at the reference lag gives 70.992/70.992; it does not rescue the reference-lag error.

## Requested summary table

| Diagnostic | Result | Assessment |
|---|---:|---|
| Independent realizations | 15 train + 5 held out | adequate realization split |
| Complete dominance switches | 316 train + 108 test | much better sampled than gap 5 |
| Held-out mean transition time | LR 34.21; RL 37.55 | seed-bootstrap intervals reported above |
| Supported diffusion fitting interval | None | none |
| Near-Gaussian-lag diffusion variation | 1.807x | no coefficient plateau |
| Reference dynamic barrier | 4.099 at max lag 0.5 | finite but lag-dependent |
| CK weighted TV | 0.0025--0.1381 | deteriorates strongly with lag |
| Hidden-current RMSE reduction | 67.5% at horizon 1 | q is not closed |
| Reference MFPT relative error | LR 92.4%; RL 106.2% | fails held-out timing |
| Reference survival max error | LR 0.262; RL 0.266 | full FPT distribution also fails |

## Hidden-current closure test

Held-out prediction improvements after adding neck current:

- horizon 0.10: RMSE 0.017617 -> 0.010573 (39.98% improvement), conditional residual-current correlation 0.800.
- horizon 0.50: RMSE 0.07213 -> 0.022674 (68.57% improvement), conditional residual-current correlation 0.949.
- horizon 1.00: RMSE 0.13793 -> 0.044794 (67.52% improvement), conditional residual-current correlation 0.946.

## Interpretation rule

The potential is accepted as a predictive one-dimensional reduction only if its coefficients and barrier are stable over lag, CK discrepancies are small, adding neck current gives little held-out improvement, stationary current is small, and held-out MFPTs agree.  Failure of any of these checks is reported rather than repaired by interpolating unsupported states.

See the CSV tables and numbered figures in this directory for the complete diagnostics.

## Scientific conclusions

1. Relative to gap 5, gap 50 greatly improves transition counts (424 complete
   switches here), so insufficient data are no longer the main limitation.
2. The native process changes mostly in one-fish increments, but its increments
   become approximately Gaussian only after coarse aggregation around 0.5--1
   time unit.  Across that same range the inferred diffusion changes by
   1.807x, so there is no Huang-style short-lag plateau.
3. Chamber imbalance is not sufficient: neck current reduces held-out error by
   40--69% and remains strongly correlated with the q-only residual.
4. A finite two-well curve is reproducible as a phenomenological landscape,
   but the predeclared lag-0.5 model overpredicts held-out MFPT by roughly a
   factor of two and misses the held-out survival curve.  The good lag-1 fit is
   lag selection, not parameter-free prediction.
5. The dominant failure is hidden directional memory combined with temporal
   coarse-graining, not loss of metastability or too few switches.  Larger gap
   improves statistics but does not produce a closed one-dimensional continuous
   diffusion.
