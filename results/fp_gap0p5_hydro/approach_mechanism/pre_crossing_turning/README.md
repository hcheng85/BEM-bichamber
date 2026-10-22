# Pre-crossing turning mechanism

This directory contains a strictly pre-event replay of the deterministic angular
dynamics of the hydrodynamic gap-50 bichamber trajectories.  The main result is
in [REPORT.md](REPORT.md).

The raw 3.7-GB trajectory files and 401-MB LU factorization are intentionally
not committed.  With the archive mounted at the paths below, reproduce with:

```bash
MPLCONFIGDIR=/tmp/mplconfig python3.10 scripts/analyze_pre_crossing_turning.py \
  --analysis-dir . \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --cache-dir /path/to/archive_one_R_60/gap50_hydro_approach_features \
  --geometry-file /path/to/archive_one_R_60/geo_data/BEM_geometry_R50_wid50.mat \
  --lu-file /path/to/archive_one_R_60/LU_py/lu_R_50_wid_50.npy \
  --piv-file /path/to/archive_one_R_60/LU_py/piv_R_50_wid_50.npy
python3.10 scripts/summarize_results.py .
```

The second command applies the periodic-angle audit and generates the paired
seed contrasts and final figures.  No simulation and no fitted potential are
used.
