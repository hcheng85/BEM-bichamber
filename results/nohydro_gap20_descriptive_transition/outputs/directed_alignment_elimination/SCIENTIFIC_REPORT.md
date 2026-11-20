# Directed-alignment elimination in the polarized source school

## Decision

The production dynamics admit an accurate small-heading-deviation equation,
but the relative headings **cannot be adiabatically eliminated during a gap
approach**. The obstruction is not loss of polarization: median source-school
polarization is 0.971. Instead, the slow directed alignment mode relaxes on a
median time of 1.85 (held-out seeds), while the tested pre-entry response lasts
only 0.9 and the visible graph changes appreciably within 0.8. Consequently,
the left-null-vector response of a frozen directed graph is not predictive.

This route therefore does not yet give a closed low-dimensional turning law.
It does explain why nonreciprocity matters: it changes both the slowly evolving
relative-heading pattern and the distribution of influence over the extended
school. Those quantities retain memory throughout an approach. Replacing them
by their frozen-graph quasistatic limit removes the mechanism that the causal
intervention exposed.

## Production rule and exact collective projection

For a visible neighbour `j` of fish `i`, the production visual rule is

```text
w_ij = 1 + cos(beta_ij),                 Z_i = sum_j w_ij,
alpha_ij = w_ij / Z_i,
omega_i^vis = sum_j alpha_ij
              [Ip sin(theta_j-theta_i) + r_ij sin(beta_ij)].
```

Here `beta_ij` is the signed bearing of `j` from fish `i`; visibility is the
production Delaunay/wall-filtered directed graph. Thus `alpha_ij` is directed
both because `w_ij` is receiver-heading dependent and because each receiver
has its own full-neighbour denominator. The analysis replays this rule, the
production wall torque, its wall clip, and the final total-torque clip. No
production source was changed. The replayed production file has SHA-256
`5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`.

For the original source-school membership on a crossing-free interval, define

```text
P exp(i Theta) = (1/Ns) sum_i exp(i theta_i),
delta_i = wrap(theta_i-Theta),
c_i = cos(delta_i)/(Ns P).
```

The mean-heading identity is exact:

```text
dot Theta = sum_i c_i dot theta_i.
```

It remains the correct source-level quantity; an arithmetic heading mean is
not silently substituted for it.

## Small-deviation directed operator

Let `B_ij=alpha_ij` for internal source-source edges and zero otherwise, and

```text
L = diag(B 1) - B.
```

All source-external alignment, positional visual turning, wall turning, and
the exact clip correction are collected in `f`. At fixed positions and visible
graph,

```text
d theta = [-Ip L delta + f] dt + In dW + O(|delta|^2) dt.       (1)
```

The coefficient dependence on heading does not add a first-order alignment
term: its derivative multiplies `sin(theta_j-theta_i)=0` at the aligned state.
It does enter at higher order and through subsequent graph evolution.

Directedness matters. The right zero mode is

```text
L 1 = 0,
```

whereas the normalized left zero mode satisfies

```text
pi^T L = 0,      pi^T 1 = 1.
```

Generally `pi != 1/Ns`; its median effective support is only 6.26 fish in the
held-out data. Ordinary mean heading is therefore not conserved by internal
alignment. For a frozen, connected operator the directed collective phase
`Phi=pi^T theta` obeys

```text
d Phi = pi^T f dt + In pi^T dW.                                (2)
```

With `Q=I-1 pi^T`, the relative equation is

```text
d delta = [-Ip L delta + Q f] dt + In Q dW.                    (3)
```

If the nonzero modes were fast, the controlled leading approximation would be

```text
delta_* = Ip^(-1) L^# Q f,       dot Theta_ad = pi^T f,         (4)
```

where `L^#` is the inverse on the `pi`-orthogonal subspace. This is the tested
parameter-free approximation. Spatial extent is retained in `L` and `f`:
every coefficient uses the fish positions, directed visible edges, receiver
normalization, and spatially varying wall and positional inputs. Equation (4)
would cease to be low dimensional if those microscopic fields had to be
supplied at every time; the present test asks first whether they can at least
be eliminated dynamically.

For an evolving graph an additional term `dot(pi)^T theta` appears in
`d(pi^T theta)`. This term is negligible only if the graph changes slowly
compared with directed relaxation. Membership jumps at a crossing create a
separate discontinuity; all present tests stop before the first crossing.

## Data and validation design

The fixed 337-event prospective cohort is used without changing its event
definitions. Six crossing-free snapshots from lead 1.0 to lead 0.1 give 2,022
operator evaluations. Seeds 0--29 are the development set and seeds 30--39
are held out. The integrated test compares the observed source-heading change
over that interval with a quadrature of (2) evaluated along the archived path.
That quadrature is a stringent diagnostic, but is retrospective because it
uses later archived graphs and forcing; it is not presented as a prospective
forecast.

The directed-Laplacian self-test verifies distinct left and right null vectors.
Five of 2,022 production snapshots have more than one zero mode and are
reported rather than discarded. No new microscopic trajectories or fitted
coefficients are used.

## Small-angle accuracy

The strongly polarized expansion is locally valid. On held-out snapshots:

- median `P=0.971` and median relative-heading standard deviation is 0.242 rad;
- the exact-state linearized collective rate has correlation 0.997
  (seed interval 0.991--1.000) with the exact replay;
- its RMSE is 0.0126 (0.0048--0.0211), compared with an actual-rate standard
  deviation of 0.1466;
- 99% of absolute linearization errors are below 0.0842. The maximum over all
  2,022 snapshots is 0.310, so the approximation is not uniformly accurate in
  rare snapshots.

The successful exact-state test retains the measured `delta`. It validates
equation (1), not elimination of the relative state.

## Failure of time-scale separation

The nonzero eigenvalue with smallest real part gives the frozen-operator
asymptotic relaxation time `tau_rel=1/(Ip Re lambda_2)`. On held-out seeds:

```text
median tau_rel                         1.845
10--90% range                          1.154--3.683
fraction tau_rel < 0.5                 0
fraction tau_rel < 1.0                 0.030
```

Nonnormal transients can make a directed operator relax even less simply than
this eigenvalue estimate; they cannot create the missing clean scale
separation established by these numbers.

The frozen quasistatic relative pattern also fails directly. Its median
correlation with the observed relative headings is -0.152 on held-out data,
and its median RMSE is 3.61 times the observed relative-heading standard
deviation.

Meanwhile, from lead 1.0 to lead 0.2:

```text
median L1 change in pi                 0.508
median directed-edge Jaccard overlap   0.599
```

Thus the object presumed frozen in (2)--(4) changes on the same time scale as
the attempted elimination.

## Instantaneous and integrated response

At lead 1.0 on 129 held-out events:

| approximation | correlation (seed interval) | RMSE (seed interval) |
|---|---:|---:|
| frozen directed `pi^T f` | 0.062 (-0.043, 0.225) | 0.501 (0.285, 0.798) |
| uniform mean forcing | 0.703 (0.649, 0.744) | 0.108 (0.102, 0.112) |
| linearized rule retaining measured `delta` | 0.997 (0.991, 1.000) | 0.0126 (0.0048, 0.0211) |

The directed result is not caused solely by its largest near-reducible-graph
outlier: after removing the largest 5% of absolute directed predictions, its
held-out correlation is still only 0.149 and its RMSE is 0.273.

Over the archived lead-1.0 to lead-0.1 path:

| path quadrature | correlation (seed interval) | RMSE (seed interval) |
|---|---:|---:|
| directed `pi(t)^T f(t)` | 0.201 (0.054, 0.325) | 0.219 (0.188, 0.249) |
| uniform forcing | 0.555 (0.418, 0.644) | 0.112 (0.100, 0.121) |

The uniform forcing is a better empirical approximation, but it is not a
closed prediction. It is evaluated from the full microscopic forcing, its
integrated version uses the future archived path, and its RMSE remains large
relative to the observed integrated-response standard deviation (0.134).
Most importantly, equation (1) shows why it works better here: the extant
relative-heading field supplies an internal-alignment correction that the
quasistatic directed formula gets wrong. Dropping that field is not a
controlled derivation of a reciprocal collective law.

## Relation to the causal symmetrization experiment

If internal alignment were allowed to equilibrate under a frozen graph,
symmetrizing `alpha_ij` would replace the directed influence weights by the
uniform weights. The formal gap-alignment response at the initial state is

```text
Delta dot a_ad = -sin(psi) [mean(f) - pi^T f].                  (5)
```

For the exact 14 configurations used in the matched-noise intervention, (5)
predicts a mean first-unit symmetric-minus-baseline effect of -0.132. The
observed effect is -0.0285. Across states the correlation is 0.365 with a
seed-bootstrap interval of -0.088 to 0.736, sign agreement is 9/14, and RMSE
is 0.265. Equation (5) has the correct mean direction but greatly overstates
the response and does not resolve statewise prediction.

This does not contradict the causal result. The intervention proved that
nonreciprocal internal alignment contributes to turning and sustained
transfer. The present analysis shows that its effect is carried by a slowly
evolving relative-heading/visible-graph process, not by instantaneous
relaxation to the left-null-vector phase. During the branches, wall and
positional forcing, clipping status, graph edges, and relative headings all
respond before equilibration; the earlier experiment directly measured that
compensation.

## Physical conclusion and implication for reduced theory

The minimal supported mechanism remains

```text
directed visual interaction + persistent relative-heading organization
    -> collective turn
    -> coherent approach supply
    -> storage and neck flux.
```

The first arrow is not Markovian in the present global approach geometry. Its
memory time is comparable to the approach duration, and its operator is tied
to the spatially extended, rapidly rewiring visual graph. Therefore the
microscopic relative headings cannot be analytically eliminated while
retaining the measured nonreciprocal response.

This is a precise closure failure, not evidence for adding empirical moments.
A predictive low-dimensional transition-rate theory would require a separate
controlled limit in which the directed visual graph self-averages or remains
quasistationary for many alignment-relaxation times. The archived narrow-gap
approaches do not satisfy that limit. Without such a derivation, the honest
reduced description must retain a finite-memory turning response (for example
inside an attempt-level renewal theory) rather than claim an instantaneous
potential or fitted turning function.

## Reproduction

From `results/nohydro_gap20_descriptive_transition/`:

```bash
uv run --with numpy --with scipy --with matplotlib \
  python analyze_directed_alignment_elimination.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python analyze_directed_alignment_elimination.py \
  --raw-dir /Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data \
  --cohort outputs/collective_qab_pathway/event_outcomes.csv \
  --source ../../../BEM_bichamber_parallel_instance_no_hydro.py \
  --branch-effects outputs/nonreciprocal_alignment_intervention/per_seed_effects.csv \
  --output outputs/directed_alignment_elimination
```
