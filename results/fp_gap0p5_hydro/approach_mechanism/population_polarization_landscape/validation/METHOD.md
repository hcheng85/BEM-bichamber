# Method and reproducibility

## Sources and definitions

Raw trajectories:

`/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data/If1e-02_Ip40_In05_R50_num100_wid50_n{0..19}_BEM.mat`

Each file is about 4,000,009,536 bytes and contains:

- `state`, shape `(3,100,1000002)`: (x,y,\theta);
- `rdot`, shape `(2,100,1000002)`: total translational velocity;
- production parameters with native `dt=0.001`.

The synchronized extraction reuses the already validated stride-10 caches:

- `gap50_hydro_approach_features/gap50_hydro_approach_seedXX_stride10.npz`;
- `gap50_hydro_fp_validation/features/gap50_hydro_seedXX_stride10.npz`.

Those caches are outputs of `approach_mechanism/scripts/extract_approach_features.py` and the prior (q,J) extraction. Both derive from the same raw index, and this analysis verifies identical time axes and occupancy to (1.1\times10^{-16}).

Chambers are exhaustively assigned by `x<0` and `x>=0`, matching the existing event analysis. The opening has no finite (x)-length aperture class, so `N_gap=0`. Local polarization uses unit **headings**, not `rdot`. Empty-chamber polarization is undefined (`NaN`). Opening-directed mean heading is (h_L=\langle\cos\theta\rangle_L) and (h_R=-\langle\cos\theta\rangle_R). The diagnostic net rightward heading is

[
H=(N_L h_L-N_R h_R)/N.
]

The retrospective increment is strictly past-looking:

[
\Delta m_{past}(t)=m(t)-m(t-0.1).
]

The stored `J_neck` is retained solely for provenance/comparison and is never included in the proposed (z) state.

## Derived table

`tables/synchronized_features.npz` is a compressed columnar table. `seed` has shape `(20,)`, shared `time` has shape `(80001,)`, and every data column has shape `(20,80001)`. Columns are:

`N_L,N_R,N_gap,m,P_L,P_R,valid_L,valid_R,dm_past,heading_to_gap_L,heading_to_gap_R,net_heading_right,J_neck,q_l050,global_P,vertical_center`.

The table retains times 200.0--1000.0 at spacing 0.01 and occupies 34,042,404 bytes. The manifest records every raw path, byte size, sample count, empty-stratum count, and identity check. Vertical center is retained, but no phase is invented because phase was not a saved source field.

## Occupancy visualization

Training seeds are 0--14 and held-out seeds are 15--19. Mixed-state (P_L,P_R) histograms use three prespecified imbalance slices: `[-.75,-.2)`, `[-.2,.2)`, and `[.2,.75)`. The displayed 12-by-12 occupancy landscape masks cells with fewer than 100 training or 30 held-out observations. Robustness uses 8, 12 and 16 bins and reports train--test Jensen--Shannon divergence. Boundary states (m=\pm1) are plotted separately using only occupied-group polarization.

## Finite-state models

The imbalance-only model has 10 states: exact left/right full-occupancy boundaries plus eight mixed intervals with edges

`[-1,-.75,-.4,-.2,0,.2,.4,.75,1]`.

The proposed model has 38 states:

- each full-occupancy boundary is divided into training occupied-polarization tertiles (3+3 states);
- each mixed imbalance interval is divided at its training-conditional (P_L) and (P_R) medians (8 x 4 states).

All bin edges are learned from seeds 0--14 only. All 38 states exceed the support thresholds. No symmetry pooling is imposed.

For lags 0.1, 0.5 and 1.0, direct (P(2\tau)) is compared with (P(\tau)^2). Reported total variation first collapses the target over polarization to compare future imbalance; full 38-state error is also retained. Origin rows require 500 training samples. The train--test (P(2\tau)) difference, using at least 200 held-out samples per row, is the sampling/heterogeneity reference. The practical pass criterion is the same as the earlier (q,J) analysis: both projected and full errors must be no more than 15% above their corresponding train--test references.

## Splitting probability and memory

A corridor excursion starts when (|m|\leq0.2) is entered from outside and ends on first hitting (m\leq-0.4) or (m\geq0.4). Frames are sampled every 0.1 within the excursion and inherit its future boundary. The two final incomplete excursions are censored and receive no label.

The operator committor solves the discrete backward equation using each training transition kernel. Held-out Brier score and log loss are computed separately in each seed.

The hidden-memory diagnostic estimates beta-smoothed empirical probabilities on training excursions for:

1. (z) only;
2. (z) plus the sign of strictly past (Delta m);
3. (z) plus the sign of (H).

All three are evaluated on the same held-out observations for which every training conditional cell has at least 50 examples. These augmented variables are diagnostics of failure, not proposed additions to the landscape. Confidence intervals resample the five held-out seed scores.

## Nucleation

Only (m=\pm1) frames are used. Every 0.1 time unit, the outcome records whether any fish leaves the occupied chamber within the following 0.5. Training estimates compare occupied-polarization boundary state alone with the same state plus the sign of occupied-group heading toward the opening. Overlapping frames are not treated as independent; reported uncertainty is across held-out seeds.

## First passage and censoring

A dominance residence begins at (m\leq-0.4) or (m\geq0.4) and ends on first reaching the opposite threshold. The final incomplete residence in every trajectory is right-censored. Model survival propagates the empirical training kernel after making the opposite dominance set absorbing. Observed survival is Kaplan--Meier. Restricted means integrate from 0 to 180. Results are shown for lags 0.5 and 1.0; seed intervals resample the five validation trajectories.

The (q_\ell,J) benchmark is copied verbatim from `closure_qJ/tables/splitting_scores.csv` and `first_passage_validation.csv`. Its smooth occupancy definition produces slightly different residence counts, so it is a scale benchmark rather than a paired episode comparison.

## Commands

From the repository root:

```bash
UV_CACHE_DIR=/private/tmp/uv-cache uv run --offline \
  --with numpy --with scipy --with pandas --with matplotlib \
  python results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/validation/scripts/extract_population_polarization.py \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --approach-cache /Users/hao/Library/CloudStorage/Dropbox/code/bichamber/archive/archive_one_R_60/gap50_hydro_approach_features \
  --qj-cache /Users/hao/Library/CloudStorage/Dropbox/code/bichamber/archive/archive_one_R_60/gap50_hydro_fp_validation/features \
  --output results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/validation/tables/synchronized_features.npz

XDG_CACHE_HOME=/private/tmp/fontcache MPLCONFIGDIR=/private/tmp/mplconfig \
UV_CACHE_DIR=/private/tmp/uv-cache uv run --offline \
  --with numpy --with scipy --with pandas --with matplotlib \
  python results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/validation/scripts/analyze_population_polarization.py \
  --features results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/validation/tables/synchronized_features.npz \
  --output-dir results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/validation \
  --qj-results results/fp_gap0p5_hydro/closure_qJ
```

## Limits

- This is one ((I_w,I_n)) parameter point.
- Finite binning cannot prove exact Markovianity; failure is nevertheless robust because direction/history gives large held-out gains and first-passage times are wrong by a factor of roughly two.
- Boundary nucleation labels are rare and sampled from overlapping risk windows; only seed-level uncertainty is interpreted.
- The stationary landscape is descriptive and is never used as a potential.
