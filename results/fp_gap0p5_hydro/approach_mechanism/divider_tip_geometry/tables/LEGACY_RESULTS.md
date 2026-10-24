# Legacy observed-velocity results

The following files preserve the original calculation in which the tangent
velocity was taken from the observed trajectory while angular drift was
evaluated at a counterfactual heading:

- `tangent_snapshot_replay.csv`
- `class_frequencies.csv`
- `prospective_association.csv`
- `event_control_matched_summary.csv`
- `geometry_matched_class_pairs.csv`
- `event_control_geometry_matches.csv`
- `representative_traces_legacy_observed_gate_velocity.csv`
- `representative_trace_summary_legacy_observed_gate_velocity.csv`
- `representative_trace_selection_legacy_observed_gate_velocity.csv`
- `representative_crossings_legacy_observed_gate_velocity.csv`

These are comparison artifacts, not the primary results. Primary class and
association tables carry the `corrected_` prefix, and new trace files use each
trial heading's own replayed translational velocity.

The root-level `summary_legacy_observed_gate_velocity.csv` preserves the old
headline summary; `summary.csv` is now the corrected primary summary.
