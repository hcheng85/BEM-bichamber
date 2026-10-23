# Collective-coordinate precursor analysis

The scientific conclusions are in [REPORT.md](REPORT.md), with definitions and
sampling details in [METHOD.md](METHOD.md).

Reproduce from the existing archive (no simulation is launched):

```bash
MPLCONFIGDIR=/tmp/mplconfig python3.10 scripts/analyze_collective_coordinates.py \
  --analysis-dir . \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --cache-dir /path/to/archive_one_R_60/gap50_hydro_collective_coordinates \
  --geometry-file /path/to/archive_one_R_60/geo_data/BEM_geometry_R50_wid50.mat \
  --lu-file /path/to/archive_one_R_60/LU_py/lu_R_50_wid_50.npy \
  --piv-file /path/to/archive_one_R_60/LU_py/piv_R_50_wid_50.npy

MPLCONFIGDIR=/tmp/mplconfig python3.10 scripts/summarize_collective_coordinates.py .
```

The raw 3.7-GB trajectories, BEM factorization, and disposable coordinate
caches are intentionally not committed. `--skip-replay` reuses existing compact
replay tables when only the global coordinate calculations need refreshing.
