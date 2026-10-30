# Codex follow-up: Population imbalance and *local* polarization landscape

Work in the BEM-bichamber repository. Start by reading:

- `results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/THEORY_m_PL_PR.md` **(primary hypothesis and definitions)**
- `results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/REPORT.md` **(existing evidence, not a finished potential)**
- `results/fp_gap0p5_hydro/approach_mechanism/REPORT.md`
- `results/fp_gap0p5_hydro/approach_mechanism/collective_coordinates/REPORT.md`
- `results/fp_gap0p5_hydro/closure_qJ/REPORT.md`
- The authoritative BEM simulator's saved-field definitions.

## Scientific question

Test the user's revised **three-coordinate proposal**, without imposing
current as a required state coordinate:

```
m = (N_R-N_L)/N,
P_L = |sum_{i in L} e(theta_i)/N_L|,
P_R = |sum_{i in R} e(theta_i)/N_R|.
```

The sign of a change in `m` identifies net transfer direction. Local
polarization magnitudes characterize *separate chamber schooling* rather
than cancellation in a global polarization. Test whether `(m,P_L,P_R)`
also suffices for *future* transition probability and switching times.
Earlier analysis of 512 nucleation episodes found that global `P` alone
poorly distinguishes failed from committed transitions (`AUC=0.553`),
whereas a signed heading toward the opening is informative (`AUC=0.793`
at first crossing). That comparison does not yet test local polarization.
Do not claim an effective energy barrier before validating the dynamics.

## Scope: one finite, falsifiable experiment

1. **Locate existing raw output** for the 20 `R50/gap50` full-hydrodynamics
   seeds in the configured archive or CARC mount. Do not download, clone, or
   reproduce the multi-GB microscopic data. If inaccessible, state the
   missing paths and stop; existing aggregate CSVs are not sufficient to
   reconstruct the joint `(m,P_L,P_R)` distribution.
2. Read the simulator and original analysis scripts to ascertain: chamber
   definition when fish occupy the aperture, the exact output cadence,
   polarization based on *heading* vs total velocity, the sign convention
   for opening-directed heading, and the definition of the prior neck current
   `J`/occupancy `q_l`. When a chamber is empty, its `P_c` is undefined:
   do not substitute zero. Track empty-chamber masks and any gap occupants,
   including a documented exhaustive population-assignment convention when
   used. Preserve time ordering and seed identity, and exclude the same
   first 200 time units as prior work.
3. Produce one compact time-resolved table, at the **existing** 0.01-time
   output cadence if present, containing `seed,t,N_L,N_R,N_gap,m,P_L,P_R,
   valid_L,valid_R,dm_past,heading_to_gap_L,heading_to_gap_R` (and, if
   already in the source, `J_neck,vertical_phase`). The *past* difference
   `dm_past=m(t)-m(t-h)` must be computed using frames strictly before or
   at `t`; do not use `m(t+h)` in a predictor for future commitment. Include
   optional global `P` only as a diagnostic. Compress the derived table,
   report file size and retained time sampling, and leave large raw archives
   in place.
4. Visualize `rho_ss(m,P_L,P_R)` separately for training seeds 0--14 and
   held-out seeds 15--19, with sample support and robustness to binning;
   plot `U_occ = -log(rho_ss)` only where density is supported. Label it a
   **stationary occupancy landscape**, not free energy or quasipotential.
   Because `P_L` or `P_R` is undefined when its chamber is empty, handle
   the `m=+/-1` boundary strata separately using the occupied group's `P`.
   Do not interpolate missing polarizations or drop the dominant metastable
   residence states. Consider 2D slices/conditional occupancy visualizations
   rather than an unreliable unsupported full 3D surface.
5. Test **predictive closure, not just correlation**:
   - Within supported, matched `(m,P_L,P_R)` bins, compare *future* splitting
     probability (left versus right chamber-dominance first) conditional on
     the **past** sign of `Delta m`. Also compare a split by signed
     chamber-mean heading toward the gap. These splits are diagnostic only:
     do not introduce `J` as a required 4th variable before the test.
   - Analyze nucleation (`|m|=1`, opposite group absent) separately from
     recruitment (at least one fish on each side). Avoid treating the
     unavailable empty-group `P` as a measured physical value.
   - Fit the finite-state transition operator in the proposed variables at
     several justified coarse lags (e.g., 0.1, 0.5, 1.0), where populated;
     test Chapman–Kolmogorov consistency, seed-held-out splitting/committor,
     first-passage survival and time scales. Compare state-only and
     past-increment/phase-conditioned predictors to detect hidden memory.
   - Use whole seeds for holdout and uncertainty. Do not randomize frames
     across train and test. Treat right-censored transition episodes correctly.
   - Compare side by side with the already validated `(q_l,J)` closure. No
     new fitted free-energy surface or neural surrogate.
6. End with a decisive statement: **(A)** autonomous `(m,P_L,P_R)` closure is
   predictive at a specified coarse lag; **(B)** useful visual/descriptive
   state but not dynamically sufficient; or **(C)** no meaningful added value
   over occupancy/current. For (A), show held-out passage statistics and
   calibration. For (B)/(C), state the missing physical variable and stop;
   do not keep adding arbitrary coarse coordinates.

## Analytical follow-through, only if justified

Use the exact microscopic angular noise amplitude and wall/social/hydrodynamic
equations to derive leading asymptotic approximations for (i) first-fish
angular entry and (ii) post-entry recruitment. In a well-defined stationary
event-counting regime, school-wide commitment frequency can be decomposed
into first-entry frequency times the conditional fraction that commits, but
do not assume a Poisson process or identify inverse stationary flux with MFPT
without validating renewal assumptions. Aim ultimately at an unfitted
prediction for the wall-strength × angular-noise optimal ridge; report
explicitly when the current data cannot support that prediction.

## Deliverables

- A new, separately named results `REPORT.md` with an unambiguous pass/fail
  and physical takeaway; preserve this hypothesis document and earlier report.
- Reproducible extraction/analysis scripts, compact derived tables, at most
  three essential figures and figure captions.
- The exact commands, seed split, numerical diagnostics, and raw-source paths.
- Preserve existing analyses. Work on a fresh branch, commit; push only if
  GitHub credentials/permissions already work. Never force-push.
