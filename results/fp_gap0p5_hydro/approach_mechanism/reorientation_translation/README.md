# Reorientation-to-translation validation

See [REPORT.md](REPORT.md) for the physical conclusion and [METHOD.md](METHOD.md)
for definitions. The analysis uses only the existing 20 full-BEM trajectories.

Reproduce with:

```bash
MPLCONFIGDIR=/tmp/mplconfig python3.10 scripts/analyze_reorientation_translation.py \
  --analysis-dir . \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --cache-dir /path/to/archive_one_R_60/gap50_hydro_collective_coordinates

MPLCONFIGDIR=/tmp/mplconfig python3.10 scripts/summarize_and_plot.py .
```

Use `--skip-kinematics` to reuse the compact native-integration table. Raw
trajectories and disposable coordinate caches are intentionally not committed.
