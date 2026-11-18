# Refined directional-feedback geometry review

## Decision

The revised audit resolves the earlier ambiguity as an inconsistency between the physical production boundary and the native social-visibility filter. The directional branching experiment was therefore stopped before any new restart simulations.

The nominal condition `|y| <= gap/2` is only a diagnostic classification and was not used as an extra wall predicate. Instead, the review reconstructed the exact boundary coordinates and curved components used by the archived simulator, overlaid every exceptional sightline, evaluated the native predicates, and audited all 658 indexed first-crossing frames.

## Provenance and fixed-cohort audit

- Production source SHA256: `5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`.
- Restart adapter SHA256 matched the committed pilot.
- The selected 28-state manifest matched exactly.
- Native `CheckVoronoi` and the optimized equivalent agreed on all 28 snapshots.
- Receiver-specific social-torque partition closed exactly.
- The attached script passed its self-test for neutral geometry classification, directional crossing counts, exact population balance, social masks, and 14-seed summaries.

## Actual boundary versus native visibility

Across the 28 frozen states, 276 directed native-admitted links connected fish on opposite signs of `x`. Of these, 246 intersected `x=0` within the nominal central interval and 30 directed links (15 unique) were off-central. All exceptional links occurred in `pre_s1_k985000` and `pre_s12_k240000`.

The exact production boundary contains rounded upper transition arcs centered at `circle_center2=(-1, 6.5834)` and `circle_center1=(1, 6.5834)`, each of radius `0.5`. All 15 unique exceptional sightlines intersect both of these physical boundary arcs: they leave the left chamber through its rounded cap, traverse exterior space, and enter the right chamber through its rounded cap.

The discrepancy is in native visibility. `CheckVoronoi` rejects links intersecting four straight wall segments and the two central tip circles `circle_center5` and `circle_center6`. It does not test the upper transition arcs (nor all other curved boundary components). Consequently, its segment and central-tip predicates return false for these 15 links even though the exact boundary overlay shows wall intersections. The fish endpoints themselves lie within their respective chambers; the problem is the connecting line of sight.

Therefore the earlier off-central diagnostic did identify physically obstructed links, but not merely because they were outside `gap/2`. The source-grounded reason is that they cross curved production-boundary arcs omitted from the native visibility predicate.

## Indexed first-crossing audit

All 658 indexed episode frames were read from the archived native trajectories and checked by interpolating each sign-changing particle segment at `x=0`.

- All 658 frames contained forward motion in the indexed direction.
- There were 659 individual particle crossings; one frame contained two simultaneous forward crossings.
- No indexed particle crossing was off-central.
- The largest indexed `|y|` was `1.40674`, below `gap/2=1.43178`.
- Population imbalance reconstructed from the native pre-frame agreed with indexed `q_pre` to `1.11e-16`.

Thus the event index and observed physical crossings support a single central opening. They do not explain or legitimize the through-wall visual links.

## Scientific consequence

No source-received or target-received directional ablation was run, and previous `full`/`crossgap_off` endpoint parity was not replayed because the physical-geometry gate precedes simulation. This run therefore provides no new estimate of recruitment, return suppression, retention, or distinct-fish recruitment.

The earlier global opposite-sign ablation cannot yet be assigned solely to social feedback through the bottleneck: native opposite-chamber coupling can include through-wall sightlines. Although all 14 post-entry starting snapshots had only central cross-links, unverified links could arise during five-unit branches. Correcting or redefining visibility would change the microscopic model, which this task explicitly forbids.

The appropriate next step is a separately authorized simulator correction and renewed parity/validation study, not reinterpretation of the present pilot or a parameter sweep. No production source, raw trajectory, or event definition was modified here.
