# Collective pathway and (q,a,b) closure diagnostic

## Scope
Read-only analysis of existing N=100, gap-20, no-hydrodynamics trajectories. No new microscopic simulation, fitted potential, fitted drift, or causal intervention.
337 fixed-common-cohort majority-to-minority first crossings; 185 transferred at least four ADDITIONAL fish during the next five units.
Independent original seeds represented: 37. Event-label groups are retrospective, not an online classifier.

## Physical coordinates
- `q=(N_R-N_L)/N` is signed population imbalance. Source side is fixed from q immediately before the first crossing (or from current q on prospective time grids); ties are excluded.
- `a=P_source . (-COM_source/|COM_source|)` is the gap-directed source polarization. When no source fish remain or its COM is undefined, a is missing.
- `b=(1/(N*ell)) sum_{source approach} (ell-u_i)`, with `u_i=source_sign*x_i` and `0<=u_i<=ell`, `|y_i|<=gap/2`. This equals the weighted-storage observable W/(N*ell) from the prior conservation identity.
- Event figures freeze the ORIGINAL source chamber so q, a, and b remain interpretable through changes of majority.

## Event-aligned dynamics
The figure compares continued and failed crossing episodes at fixed offsets on one common cohort: no prior crossing for at least five units and post-transient throughout. Uncertainty is clustered by original seed.
Chronology must be read from the curves: the code deliberately does not assert that a changes before b.

## Exact-q matched controls
Matched measurements: 658, at leads 2 and 1. Controls were selected using stable future no-flux conditions: retrospective case-control comparisons, NOT an unbiased crossing hazard.
- Lead 1.0: b(event)-b(control)=0.1488 (seed CI [0.12796, 0.17414]); 37 event seeds.
- Lead 2.0: b(event)-b(control)=0.05611 (seed CI [0.038931, 0.074699]); 37 event seeds.

## Prospective held-out diagnostic
Fixed time grid, no use of future event direction; seeds 0-29 are training and 30-39 are held out. Within each predeclared q/a cell, training b is split at its median; cells need adequate sample and independent-seed support. Unsupported cells fall back to q+a.
Held-out samples 1577, supported q/a/b cells 1577; supported coverage 1.
P(>=4 future net fish) Brier: q+a 0.052903, q+a+b 0.0511.
Seed-mean improvement (positive favors b): 0.0017988, CI [0.00060401, 0.0028547].
Future-net-fish MSE: q+a 21.059, q+a+b 20.506.
Seed-mean net-MSE improvement: 0.54697, CI [-0.086979, 1.4052].
These are nonparametric diagnostic conditional means, NOT a physically derived closure.

## Residual memory
Supported held-out states: 1577. Association between lagged a/b and future-net residual is a necessary closure check, not a sufficiency proof.
- a_prev: Spearman 0.18168, seed CI [0.11968814319792274, 0.2463011601862222].
- b_prev: Spearman 0.21485, seed CI [0.17966332489054995, 0.2482332764623463].

## Interpretation rule
Only describe an observed reorientation-accumulation-release sequence if the event-aligned ordering, control contrast and held-out storage test support it jointly. Even then, chronology alone is not causal evidence. If storage adds no held-out information, reject (q,a,b) as an adequate predictive state and revisit spatial/circulation coordinates. Full opposite-chamber relocation is not measured here.

## Reproducibility
Outputs include `event_paths.csv`, `event_outcomes.csv`, `matched_b_controls.csv`, `prospective_qab.csv`, `heldout_predictions.csv`, `qa_b_support.csv`, `summary.json`, `figures/` (including initial-imbalance-stratified pathways). Original MAT data and source files remain unchanged.
