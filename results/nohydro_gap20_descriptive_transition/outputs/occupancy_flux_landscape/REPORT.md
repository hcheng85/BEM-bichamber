# Population imbalance and past neck flux in no-hydrodynamics gap-20 trajectories

## Question and definitions

This is an empirical test using the existing 40 independent, N=100 trajectories. It does not fit a Fokker--Planck equation or interpret an occupancy histogram as energy. After the first 200 time units,

`q=(N_R-N_L)/N`, and `J_past=[N_R(t)-N_R(t-2)]/2` fish per unit time.

At prospective origins spaced by 5 time units, the outcome is a net movement of at least 5 fish from the chamber that is the majority *at measurement time* toward the minority chamber during `(t,t+5]`. Thus neither the predictor direction nor the label uses the eventual outcome. Predictor and outcome windows do not overlap.

## Numerical audit

The supplied synthetic conservation, direction, and risk-indexing tests pass. Direct reads of the original MAT trajectories gave 120 occupancy checks and 40 native-resolution crossing-count windows. Maximum errors were 1.11e-16 in q and 2.66e-15 fish/time in flux. The net-flux identity therefore agrees with both endpoint population changes and explicitly counted left-to-right minus right-to-left centerline crossings.

The horizontal structure in the joint plot is not a continuous-flow instability: q lies on a 0.02 population lattice and J_past on a 0.5 fish/time lattice. The mild zig-zag of the binned `-log p(q)` curve is correspondingly bin/lattice sensitive. It is not used as a physical barrier.

## Results

### Population imbalance is the dominant occupancy structure

The stationary sample contains 320,040 time points. Equal weighting of the 40 seeds gives `P(|q|>=0.8)=0.528` (seed-bootstrap 95% CI 0.396--0.659) and `P(|q|<0.5)=0.360` (0.255--0.471). The extreme populations near q=+/-1 are therefore real, but the broad central occupancy shows that partial redistribution is also common.

### Flux forms intermittent, persistent directional branches

Most samples lie on the zero-flux branch, with rarer positive and negative wings. Among nonoverlapping pairs for which both 2-unit windows contain nonzero net flux, the direction persists with probability 0.663 (seed-bootstrap 95% CI 0.586--0.736; 33 informative seeds). This supports short directional bursts, not a permanently flowing state.

Conditioning on the direction relative to current population balance makes the continuation effect clear:

| current abs(q) | past toward dominant | past zero | past toward balance |
|---|---:|---:|---:|
| 0--0.5 | 0/162 | 100/1927 | 26/161 |
| 0.5--0.8 | 0/54 | 46/597 | 8/67 |
| 0.8--1 | 0/30 | 28/3330 | 2/20 |

For `|q|<0.5`, the future-transfer probability is 0.161 after flux toward balance versus 0.052 after zero flux. For `0.5<=|q|<0.8`, the corresponding values are 0.119 and 0.077. The extreme-imbalance toward-balance cell has only 20 origins and 2 positives, so its apparent 0.100 probability is not precise; the seed-bootstrap interval is saved in `conditional_future_transfer.csv`.

### Held-out prediction improves, but only modestly

Seeds 0--29 train the fixed nine-cell probability table and seeds 30--39 are untouched test seeds. The held-out set has 1587 prospective windows and 85 positives. Mean seed-level Brier score is 0.05055 for q alone and 0.04960 for `(q,J_past)`. The improvement is 0.000944 (seed-bootstrap 95% CI 0.000271--0.001582), only 1.9% of the q-only score. Past flux therefore contains reproducible continuation information, but it is a small refinement rather than a replacement for population imbalance.

## Event counts and interpretation

The prior source-indexed analysis records 658 accepted short crossing episodes across 37 seeds. This analysis has 6348 nonoverlapping prospective windows, of which 210 in 31 seeds meet the 5-fish redistribution criterion. These are windows, not unique physical events. In contrast, only 9 complete opposite-basin relocations occur, across 6 seeds. That is far too few to validate a probability law for complete school relocation.

The simplest description supported here is a **discrete, directional transfer process**: q describes how population is stored between chambers, while nonzero past flux identifies whether a short crossing burst is already propagating and in which direction. J_past cannot explain first initiation during a crossing-free history, because it is then exactly zero. The data therefore support `(q,J_past)` for short-horizon continuation, not a static two-well potential and not a complete initiation theory.

## Files

- `figures/01_occupancy_landscape.png`: empirical p(q), seed-bootstrap bands, and statistical `-log p` only.
- `figures/02_joint_q_past_flux.png`: empirical joint distribution.
- `figures/03_conditional_transfer.png`: prospective conditional probabilities and counts.
- `heldout_validation.json`: fixed train/test split and seed bootstrap.
- `occupancy_flux_raw_validation.csv`: raw MAT checks.
- `conditional_future_transfer.csv`: cell counts and seed-cluster intervals.
- `summary.json` and `analysis_config.json`: event counts, definitions, and configuration.

No effective potential, diffusion closure, or Markov sufficiency is claimed.
