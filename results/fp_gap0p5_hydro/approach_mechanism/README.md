# Bichamber approach mechanism

This directory tests the pre-transition physical mechanism in the existing 20
hydrodynamic `R50/gap50` trajectories. The scientific conclusion and limitations are
in [REPORT.md](REPORT.md); precise definitions are in [METHOD.md](METHOD.md).

From the repository root:

```bash
python results/fp_gap0p5_hydro/approach_mechanism/scripts/extract_approach_features.py \
  --raw-dir /path/to/BEM_data \
  --output-dir /path/to/local_feature_cache --stride 10

python results/fp_gap0p5_hydro/approach_mechanism/scripts/analyze_approach_mechanism.py \
  --feature-dir /path/to/local_feature_cache \
  --output-dir results/fp_gap0p5_hydro/approach_mechanism

python results/fp_gap0p5_hydro/approach_mechanism/scripts/plot_representative_episodes.py \
  --analysis-dir results/fp_gap0p5_hydro/approach_mechanism \
  --raw-dir /path/to/BEM_data --geometry-file /path/to/BEM_geometry_R50_wid50.mat

python results/fp_gap0p5_hydro/approach_mechanism/scripts/decompose_recruitment_interactions.py \
  --analysis-dir results/fp_gap0p5_hydro/approach_mechanism \
  --raw-dir /path/to/BEM_data --geometry-file /path/to/BEM_geometry_R50_wid50.mat \
  --lu-dir /path/to/LU_py
```

Dependencies: NumPy, SciPy and Matplotlib. Raw trajectories and local NPZ caches are
intentionally excluded from Git.
