# Exact social-alignment landscape for angular gap access

## Decision

**Outcome B:** the conditional angular landscape is mathematically exact at an
instant, and it contains modest new prospective information beyond initial
angular error and divider distance. It is not a closed static potential because
the visible-neighbor geometry and headings evolve on the same timescale as
angular entry.

An angular first-passage theory remains worthwhile only as a **nonautonomous
torque model with dynamically updated neighbor moments and a moving geometric
gate**. A frozen Kramers-style potential should not be developed.

## Exact derivation

For focal fish heading `theta`, let the frozen visible neighbors have headings
`theta_j` and bearings `phi_j` from the focal fish. The simulator's alignment
law is

```
Omega_align(theta) = Ip * sum_j sin(theta_j-theta)[1+cos(phi_j-theta)]
                          / sum_j [1+cos(phi_j-theta)].
```

Define

```
A1 = sum_j exp(i theta_j)
A2 = sum_j exp(i(theta_j+phi_j))
B  = sum_j exp(i phi_j)
C  = 1/2 sum_j sin(theta_j-phi_j).
```

Using `sin(a)cos(b)=[sin(a+b)+sin(a-b)]/2` gives exactly

```
Omega_align(theta) = Ip * {
    Im[A1 exp(-i theta)]
    + 1/2 Im[A2 exp(-2i theta)]
    + C
} / {n + Re[B exp(-i theta)]}.
```

The neighbor set is the simulator's Delaunay graph after exact divider-segment
and rounded-tip visibility obstruction. No neighbors, weights, harmonics, or
physical coefficients were fitted.

For all 9,120 prospective fixed-leading-fish configurations, the maximum error
between this expression and direct microscopic evaluation over 721 focal
headings was `4.71e-14`.

## Conditional landscape

The instantaneous landscape is the direct quadrature

```
V_i(theta) = - integral Omega_align_i(s) ds.
```

Because the constant term and heading-dependent denominator can give a nonzero
cycle integral, this is generally a tilted potential on the universal angular
cover rather than a periodic energy. The median absolute potential change over
one turn is 1.20. The analysis therefore retains the tilt and does not force
periodicity.

Stable orientations are zero crossings with negative torque slope. A
first-passage barrier is reported only when a stable orientation exists. Its
height is the least maximum uphill change in `V` encountered from the actual
initial heading to either rounded-tip gate boundary, considering both angular
directions. A downhill path has zero barrier; the calculation does not insert a
barrier merely because the fish is outside the gate.

Main landscape statistics:

- 99.62% have at least one stable social-alignment orientation;
- 17.60% have a stable orientation inside the opening corridor;
- among fish initially outside the corridor, 60.67% have a genuine positive
  first-passage barrier;
- their median barrier height is 0.0749;
- 35/9,120 configurations have no stable orientation and are treated as
  circulating rather than assigned a barrier.

## Prospective information

The outcome is future angular entry before A. Fish already inside the corridor
at A-2 are excluded from prediction. Models are trained on seeds 0--14 and
evaluated without refitting on seeds 15--19. Logistic coefficients are used
only as a validation tool; they are not interpreted as physical parameters.

| Initial predictors | Held-out AUC | Brier | NLL |
|---|---:|---:|---:|
| angular error + distance | 0.8262 | 0.1459 | 0.4558 |
| baseline + exact alignment drift | 0.8484 | 0.1427 | 0.4467 |
| baseline + drift + stable-orientation geometry | 0.8533 | 0.1381 | 0.4284 |
| full landscape, including genuine barrier and tilt | **0.8568** | **0.1302** | **0.4099** |
| landscape without initial angular error | 0.8323 | 0.1385 | 0.4313 |

The full landscape improves AUC in every held-out seed. The mean seed AUC gain
over the baseline is **0.0334 [0.0163, 0.0506]**. Brier score improves by
`-0.0158 [-0.0240, -0.0107]` and NLL by
`-0.0464 [-0.0671, -0.0304]`.

This is incremental information, not merely a renamed angular error. Initial
alignment closing rate is almost uncorrelated with initial absolute error
(`r=0.009`), yet adding it raises pooled AUC by 0.0222. Stable-orientation and
barrier variables add a further smaller gain. Conversely, removing initial
angular error makes the improvement inconsistent across held-out seeds. The
landscape complements rather than replaces heading geometry.

Prospectively imminent-event fish have, relative to matched controls:

- stable orientation 0.085 rad closer to the gate
  `[-0.117, -0.054]`;
- probability of a stable direction in the gate higher by 0.0399
  `[0.0272, 0.0525]`;
- genuine-barrier probability lower by 0.0598
  `[-0.0803, -0.0402]`;
- conditional barrier height lower by 0.178
  `[-0.248, -0.112]`.

The instantaneous alignment closing rate itself is not significantly different
between event and control configurations. Thus the useful information is the
shape and destination of the torque field, not simply a larger instantaneous
turning rate.

## Frozen-environment test

The relevant times are not separated:

- median observed gate-entry time: 0.429 (mean of seed medians 0.439,
  seed-bootstrap interval `[0.395, 0.481]`);
- median initial alignment-only reorientation scale
  `|eta|/Omega_closing`: 0.539.

At 0.2 time units, initial-neighbor retention is 0.683, frozen/current
landscape correlation is 0.903, relative RMS landscape change is 0.343, and
frozen/current torque signs agree for 0.763 of fish. At 0.5, near the typical
entry time, these become 0.566, 0.844, 0.459, and 0.680 respectively.

The broad landscape shape persists moderately well, explaining its held-out
predictive value, but the local torque that controls first passage is already
wrong in sign for roughly one third of fish by the typical entry time. This is
not adequate closure for a static first-passage or escape-rate calculation.

## Physical interpretation and recommendation

The exact instantaneous landscape reveals a meaningful mechanism: the social
configuration often creates a stable orientation nearer to, or inside, the
collision-free opening corridor, and lowers the angular uphill excursion needed
to access it. This is stronger than a correlation with current heading because
it predicts future entry beyond heading and distance on independent seeds.

It is not an autonomous energy landscape. Neighbor exchange, neighbor heading
evolution, focal translation, and the moving rounded-tip tangent directions
must remain dynamic. The simplest justified model is

```
d theta_i = Omega_align(theta_i; A1(t),A2(t),B(t),C(t)) dt
            + other microscopic torques dt + In dW_i,
```

with the coefficients and admissible angular interval updated from a rolling,
proximity-defined leading population on a timescale no coarser than about 0.1.
The present result does not justify a static potential, Kramers rate, or fitted
effective landscape.

