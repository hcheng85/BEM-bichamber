# Occupancy--current closure extension

This directory tests whether `(q_l,J)` closes hydrodynamic bichamber switching
more successfully than occupancy alone. It reuses the compact feature caches
from the preceding analysis and does not require new BEM simulations.

Run from the repository root:

```bash
python results/fp_gap0p5_hydro/closure_qJ/scripts/analyze_qj_closure.py \
  --feature-dir /path/to/feature_cache \
  --output-dir results/fp_gap0p5_hydro/closure_qJ
```

Dependencies are the NumPy and Matplotlib requirements already documented in
the parent analysis. The generated `tables/` files contain all plotted values;
figures are provided in PDF and PNG. Raw MAT trajectories and NPZ feature
caches are intentionally not committed.
