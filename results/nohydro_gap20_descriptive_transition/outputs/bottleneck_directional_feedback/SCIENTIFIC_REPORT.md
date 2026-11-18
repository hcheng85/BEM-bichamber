# Directional bottleneck-feedback audit: stopped before intervention

## Outcome

The required geometry gate failed, so no directional branching simulations were run. In particular, the previous `full` and `crossgap_off` endpoints were not replayed and no source-receiver or target-receiver effect was estimated.

The supplied analyzer passed its synthetic self-test. The archived production SHA256 and restart-adapter SHA256 exactly matched the committed pilot:

- production: `5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`
- adapter: `60184cdd370bc06c5040a379ca526fe9128ba5d5053d01d16577edc1a43ebcd2`

Thus the stopping condition is geometric, not a provenance mismatch.

## Aperture-visibility mismatch

The audit examined the same 28 frozen microstates selected for the previous pilot. Native `GetVoronoi` plus `CheckVoronoi` admitted 276 directed links connecting fish with opposite signs of `x` (138 undirected links).

- 246 directed links intersected `x=0` inside the intended central aperture, `|y| <= gap/2 = 1.4317821063`.
- 30 directed links (15 unique links) intersected `x=0` outside that aperture.
- All 15 failures occurred in two pre-entry configurations: `pre_s1_k985000` contributed six unique links and `pre_s12_k240000` contributed nine.
- Their divider intersections were near the upper boundary, from `y=6.637` to `7.010`, rather than near the central opening.

The archived explicit vertical blocking segments terminate at `y=6.5833949804`, while the outer radius is `7.1589105316`. Consequently, these upper cross-divider links do not intersect a segment tested by `CheckVoronoi`; they also do not intersect either central tip circle. They are therefore admitted by the native visibility rule despite not passing through the intended central aperture.

This is not a small tolerance ambiguity: the worst link lies `5.5785` beyond `gap/2`. The source/target receiver torque partition itself closed exactly at fixed geometry, with maximum residual zero, but that algebra does not rescue the physical interpretation of the affected links.

## Scientific conclusion

The previous global `crossgap_off` treatment cannot yet be interpreted purely as removal of social feedback through the bottleneck. Its link set includes native interactions across an upper geometric bypass or excluded top-notch region in at least two audited frozen states. The present audit does not establish whether those particular links contributed materially during the 14 post-entry branches, because the registered protocol requires stopping on any geometry mismatch before simulation.

Accordingly, this run provides no evidence favoring source-side recruitment, target-side feedback, competing directional effects, or either mechanism being unresolved statistically. Those alternatives were not tested.

Before resuming the causal experiment, the physical domain and visibility convention should be reconciled explicitly: determine whether the upper region is intended fluid, an unintended second passage, or outside the admissible fish domain, and then define the cross-aperture link set from that source-grounded geometry. The production simulator and archived trajectories were not modified here.

Compact evidence is stored in `audit.json` and `aperture_visibility_audit.csv`. No raw trajectories, new simulations, or large intermediate data were produced.
