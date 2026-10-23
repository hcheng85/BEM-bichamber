# Method

## Scope

The analysis reuses the 456 distance- and phase-matched event/control pairs and
the ten fixed leading identities selected at A-2. It uses all 20 archived BEM
seeds and launches no simulation.

## Landscape construction

At each A-2 snapshot, the original Delaunay and obstruction rules define
visible neighbors. Their positions and headings are frozen only for the
conditional calculation. The four complex/real moments `A1`, `A2`, `B`, and
`C` are evaluated directly. Harmonic and direct microscopic torques are checked
on a 721-point grid over `[-pi,pi]`.

The right chamber is reflected by `(h_x,h_y)->(-h_x,h_y)`. Opening bounds are
the exact tangents to the radius-0.5 rounded tips. Stable orientations are
positive-to-negative torque zero crossings. The potential is integrated without
periodic detrending, preserving its cycle tilt.

For a fish outside the gate, the barrier is the minimum over both gate
boundaries and both angular directions of the maximum positive potential rise
along that path from the actual heading. It is undefined if no stable
orientation exists and zero if a downhill route exists.

## Prediction

The outcome is first angular-corridor entry within the following two time
units. Initially admissible fish are excluded. Models are trained on seeds
0--14 and tested on seeds 15--19. The baseline contains initial absolute angular
error and divider distance. Nested models add exact drift, stable-orientation
geometry, barrier height/existence, and cycle tilt. AUC, Brier score, and NLL
are reported per held-out seed; uncertainty is a paired bootstrap over seeds.

## Frozen-environment test

At 0, 0.1, 0.2, 0.5, and 1.0 time units after A-2, the actual neighbor graph
and harmonic coefficients are recomputed. Tests report neighbor Jaccard
retention, correlation and relative RMS change of the entire torque landscape,
and frozen/current torque agreement at the observed focal heading. Fish are
averaged within trajectories and uncertainty is clustered over 20 seeds.

