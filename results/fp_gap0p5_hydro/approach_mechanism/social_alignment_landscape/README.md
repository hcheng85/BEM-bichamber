# Social-alignment angular landscape

This is an additive analysis of existing trajectories. No raw BEM data or
fitted potential is included.

From the repository root:

```bash
MPLCONFIGDIR=/tmp/mplconfig python3.10 \
  results/fp_gap0p5_hydro/approach_mechanism/social_alignment_landscape/scripts/analyze_alignment_landscape.py \
  --out results/fp_gap0p5_hydro/approach_mechanism/social_alignment_landscape \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --risk-set results/fp_gap0p5_hydro/approach_mechanism/reorientation_translation/tables/distance_matched_risk_set.csv \
  --fish-table results/fp_gap0p5_hydro/approach_mechanism/leading_fish_angular_first_passage/tables/fixed_leading_fish.csv

MPLCONFIGDIR=/tmp/mplconfig python3.10 \
  results/fp_gap0p5_hydro/approach_mechanism/social_alignment_landscape/scripts/summarize_and_plot.py \
  results/fp_gap0p5_hydro/approach_mechanism/social_alignment_landscape \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --risk-set results/fp_gap0p5_hydro/approach_mechanism/reorientation_translation/tables/distance_matched_risk_set.csv
```

Primary results are in `REPORT.md`, `summary.csv`, the three figures, and the
compact CSV tables. The compressed fish-level tables allow independent
recalculation without committing raw trajectories.

