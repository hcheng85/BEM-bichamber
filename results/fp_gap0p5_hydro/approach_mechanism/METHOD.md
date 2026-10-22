# Method and reproducibility notes

## Dataset

- Existing full hydrodynamic BEM trajectories: `R50/gap50`, seeds 0--19.
- Stored parameters: `N=100`, `Ip=4`, `In=0.5`, `Iw=0.94`, `If=0.01`,
  native `dt=0.001`, duration 1000.
- Geometry: `R=g=7.158910531638177`, hence `g/R=1` and `g/(2R)=0.5`.
- Burn-in: 200 time units.
- Aggregate feature stride: 10 native steps, giving `dt=0.01`.
- Discovery seeds: 0--14; held-out confirmation: 15--19.

The feature extractor loads one raw MAT file at a time and saves only chamber-level
counts, centers, circular mean headings, polarization, standard deviations and
coordinate extrema. These caches are local working files and are not committed.

## Coordinates

For origin side sign `s=+1` for left and `s=-1` for right, the inward normal is
`(s,0)`. Center distance is `d=-s*Xc`; the inward heading component is
`s*hx/polarization`. The front is the maximum x fish for a left origin and minimum x
fish for a right origin. The geometric condition projects the chamber-conditioned
center to the divider midline.

## Episodes

An episode is armed only by exact occupancy `N_origin=100`. Event A is the next
sample with `N_opposite>0`. Success is reaching `N_opposite>=70` before returning to
zero; that threshold time is event B. A failed episode ends at the return to zero.
Aligned failed histories are censored there.

The 10-unit slope and pre-event geometric summaries require ten complete time units
with an empty opposite chamber. The matched plot is stricter: it requires 20 units
and matches successful and failed episodes without replacement using center distance,
vertical center and inward heading at -10.

## Oscillation phase

Within every contiguous exact-full-occupancy interval of at least 12 units, `Yc` is
Gaussian-smoothed with width 0.5. Top and bottom turning points are detected with
prominence 0.75 and minimum separation 3. Same-phase pairs with periods 3--40 define
`Delta d`. Each pair is prospectively labeled by whether a successful event A occurs
within the following measured period. No decreasing envelope is selected after the
outcome is known.

## Interaction replay

For one typical successful event from each held-out seed, the archived state is read
at the first samples reaching 1, 5, 20 and 50 opposite fish. The original formulas
are separated by source-fish side. The BEM forcing is linear, so same-side and
cross-side boundary responses are solved separately with the archived LU factors and
sum to the original response. Social Delaunay edges are filtered against the solid
divider segments and rounded tips as in the simulator.

For a torque contribution `omega`, its physically comparable inward-alignment rate is

`d(n_in dot e)/dt = -n_in,x sin(theta) omega`.

The interaction intervals bootstrap the five held-out seed values. Main event
intervals bootstrap independent seed means. The full episode and plotted-value tables
are retained for audit.

## Reproduction

Run the four scripts in order from the repository root. Raw and geometry paths are
environment-specific; see `README.md`. No script writes to the raw data directory.
