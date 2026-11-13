# Does school gap-directed polarization rise immediately before crossing?

Read-only analysis of existing no-hydro gap-20 data. The measured orientation is
`a = mean(e(theta)) dot unit(gap-COM_majority)`, evaluated within the
**fixed original majority chamber** over the two units before an event, and
in a previously matched same-occupancy control over the same duration.

## Data-level diagnostic

Eligible majority-to-minority paired events: 37.
Event-minus-control observed Δa, equal independent-seed weighting:
-0.312449 (95% seed-bootstrap CI [-0.4438209909962501, -0.1898250616967678]).
Source-chamber majority consistent at all sampled times in
37/37 pairs.

The screening result is opposite to the proposed rise: event schools begin the
window much more gap-directed than controls (event-minus-control start `a` =
0.586313, 95% CI [0.3773108459188384, 0.7982535991385529]), but that
advantage falls to 0.273864 (95% CI [0.020789218363015862, 0.5212014946406899])
immediately before crossing. Thus the known predictive orientation is already
present at least two time units before initiation; it is not produced as a
late two-unit increase.

## Exact endpoint heading/geometry split

The two-endpoint identity closes to a maximum absolute residual of
6.661e-16. In event windows, net heading
change contributes 0.092075
(95% CI [0.06208125762893842, 0.12247469125284019]), while rotation
of the school-COM-to-gap direction contributes
-0.428760
(95% CI [-0.5002565257377002, -0.35936509016295903]). Membership
contributes exactly zero because no fish crosses during these pre-first-crossing
windows. Thus gap-direction geometry is the statistically resolved dominant
contribution to the absolute event-window decline: the school headings turn
slightly more gapward, but motion of the school COM rotates the gap direction
away more strongly.

For the paired event-minus-control discrepancy, heading contributes
-0.147272
(95% CI [-0.23262609565580633, -0.0693654712319354]) and gap
geometry contributes -0.165177
(95% CI [-0.25193482189941807, -0.08088097727999863]). Their
difference is not resolved: geometry minus heading =
-0.017904
(95% CI [-0.12380739791549812, 0.09262804297314618]). Therefore
neither can be called the unique dominant explanation of the paired decline.

## Conditional source-exact decomposition

At each step the deterministic angle update is replayed from the archived
production source, including the unchanged visible Delaunay social law and
forward-ray wall law. The identity is

`Δa = wall + social + clip + Gaussian_noise_relaxation + realized_noise + membership + gap_direction_rotation`.

The membership term includes both fish joining/leaving the original source
chamber and its count normalization. The geometry term is the rotation of the
majority-COM-to-gap direction; it is **not another torque**. A collision
velocity override does not enter directly because a is defined from
headings and positions, although it can indirectly affect the COM movement.

The requested decision rule permits the eight-pair replay only if the 37-pair
screening supports a positive event-specific rise. It does not, so social,
wall, clip, noise-relaxation, and realized-noise terms were not computed on
the real trajectories. Exact native-step status: {
  "available": false,
  "reason": "not run because the preregistered positive screening condition failed"
}

## Interpretation

The preregistered positive screening condition failed, so the eight-pair source-exact replay was not run. The endpoint identity separates net heading change from school geometry, but the heading term still combines social, wall, and angular-noise effects and cannot attribute a torque source.

The social term remains the original combined visual rule. Noise relaxation
is the exact Gaussian conditional mean at each native step; realized noise is
zero-conditional-mean, but not necessarily zero in a selected event sample.
The event/control comparison has only matched observational meaning.
This establishes an identity and describes mechanisms, not a closed Markov
model for (q,a) or a whole-school transition-rate theory.

## Reproducibility

- `observed_intervals.csv`: event/control endpoints and exact endpoint terms.
- `observed_paired_differences.csv`: independent-seed paired contributions.
- `observed_time_series.csv` and `observed_time_series_summary.csv`: trajectory
  samples and seed-bootstrap bands.
- `observed_polarization_change.png`: pair-level endpoint changes.
- `observed_event_control_time_series.png`: mean time course and paired contrast.
- `observed_endpoint_decomposition.png`: heading, membership, and geometry split.
- `summary.json`: screening decision, uncertainty, and identity audit.
