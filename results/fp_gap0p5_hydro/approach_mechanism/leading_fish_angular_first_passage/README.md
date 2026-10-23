# Reproduction

This directory contains derived tables and figures only; no raw trajectories
or LU matrices are committed.

Run the trajectory extraction from the repository root:

```bash
MPLCONFIGDIR=/tmp/mplconfig python3.10 \
  results/fp_gap0p5_hydro/approach_mechanism/leading_fish_angular_first_passage/scripts/analyze_leading_fish.py \
  --out results/fp_gap0p5_hydro/approach_mechanism/leading_fish_angular_first_passage \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --risk-set results/fp_gap0p5_hydro/approach_mechanism/reorientation_translation/tables/distance_matched_risk_set.csv \
  --geometry ../geo_data/BEM_geometry_R50_wid50.mat \
  --lu ../LU_py/lu_R_50_wid_50.npy \
  --piv ../LU_py/piv_R_50_wid_50.npy

MPLCONFIGDIR=/tmp/mplconfig python3.10 \
  results/fp_gap0p5_hydro/approach_mechanism/leading_fish_angular_first_passage/scripts/summarize_and_plot.py \
  results/fp_gap0p5_hydro/approach_mechanism/leading_fish_angular_first_passage
```

`tracked_fish_history.csv.gz` is a compact derived time series at 0.1 time-unit
spacing. `fixed_leading_fish.csv` contains one row per fixed fish and segment.
The primary numerical conclusions are in `mechanism_summary.csv`,
`gate_entry_budget.csv`, and `first_crosser_matched_rank_contrasts.csv`.

