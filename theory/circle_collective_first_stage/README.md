# Circular collective-motion first-stage analysis

This directory contains an additive, read-only analysis of the ten circular
full-BEM trajectories. Raw MAT files are not copied here.

Run from the repository checkout (replace the data path if needed):

```bash
python theory/circle_collective_first_stage/process_collective.py \
  --data-dir /Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data
python theory/circle_collective_first_stage/analyze_collective.py \
  --data-dir /Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data
python theory/circle_collective_first_stage/replay_torques.py \
  --data-dir /Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data \
  --source ../BEM_circular_parallel_instance.py
python theory/circle_collective_first_stage/analyze_drift_closure.py
```

The source-exact follow-up validation is run in the prescribed order with:

```bash
python theory/circle_collective_first_stage/check_stationarity_identity.py
python theory/circle_collective_first_stage/validate_drift_consistency.py \
  --data-dir /Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data
python theory/circle_collective_first_stage/analyze_polarization_response.py \
  --data-dir /Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data
```

Its compact outputs and report are under `outputs/drift_consistency/`.

The replay source is outside this publication checkout but is identified by
SHA-256 in `simulation_parameters.json` and `REPORT.md`.
