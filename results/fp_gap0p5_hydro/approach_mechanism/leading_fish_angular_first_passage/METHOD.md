# Method

## Selection and geometry

The input is the distance/phase-matched prospective risk set from
`reorientation_translation/tables/distance_matched_risk_set.csv`. The ten fish
with smallest centerline distance at A-2 are selected once. The same IDs are
followed at native `dt=0.001` until A.

For chamber sign `s=+1` (left) or `s=-1` (right), positions use
`x_local=s*x`, while headings use the reflection
`(h_x,h_y)->(s*h_x,h_y)`. The admissible interval is bounded by the exact
tangents from the fish position to the radius-0.5 rounded tip circles. Entry is
the first native sample whose reflected heading lies inside this interval.

## Update reconstruction

The archived simulator uses `RandomState(seed)`, then 100 heading uniforms, 100
radial uniforms and 100 polar-angle uniforms. Every update then consumes 100
standard normals. This recovers `In*sqrt(dt)*z` exactly unless the file was
resumed without its RNG state. The deterministic impulse is the saved heading
increment minus this stochastic increment, so simulator clipping and overrides
are included.

Wall, alignment, attraction and direct-hydrodynamic terms are replayed from the
source equations every 0.005 time units. The boundary-mediated effective term
is the residual required to close the exact deterministic budget; full BEM
solves at two snapshots per seed verify this attribution on unclipped,
non-overridden fish. Direct and BEM terms are reported separately, but their
sum is emphasized because near-singular opposite terms can cancel.

## Inference

Fish are never treated as independent replicates. Fish are averaged within
pair, pairs within seed, and the 20 seed values are bootstrapped (30,000
resamples). Event-control contrasts preserve the original matching. Eventual
first crossers are also compared with the fish of the same A-2 leading rank in
their paired control.

Seeds 10, 11, 13, 15, 17 and 19 have a verified RNG discontinuity at native
step 600,000. Segments extending beyond it are retained for geometry but marked
`rng_valid=0` and excluded from deterministic/noise contrasts.

