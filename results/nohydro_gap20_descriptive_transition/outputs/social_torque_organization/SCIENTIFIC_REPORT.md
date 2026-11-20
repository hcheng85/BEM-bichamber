# Organizational origin of collective social turning

## Physical conclusion

The collective turn is generated primarily by **nonreciprocal visual
alignment inside the source school**.  Small pairwise heading differences do
not cancel because the forward-looking visual rule assigns different weights
to the two directions of an otherwise symmetric Delaunay edge.  This converts
a weak, spatially organized heading mismatch into a net source-school turn.

The mechanism is not boundary-only, leader-only, or localized solely at the
opening.  Convex-hull fish carry disproportionate instantaneous torque, but
the continued-versus-failed separation is mainly carried by the school
interior and is shared by leading and trailing receivers.  Positional social
turning is subdominant.  Cross-chamber alignment and positional terms are each
substantial but oppose one another, leaving the internal nonreciprocal
alignment term as the robust commitment-specific contribution.

This gives a minimal physical chain:

```text
spatially organized small heading differences
  + directed forward visual weighting
    -> nonzero collective social turn
    -> mean heading remains directed toward the opening
    -> longitudinal approach current increases
    -> approach storage builds
    -> exact storage/flux law permits sustained neck transfer.
```

The first four arrows are source-grounded identities plus independently
reproduced temporal associations.  The final population transport is governed
by the previously established exact conservation law.  The analysis does not
establish that the pair-level torque can yet be predicted from a few global
coordinates: the necessary directed edge covariance is absent from `(d,psi)`
and from a perfectly coherent spatial approximation.

## Production rule and exact collective projection

For focal fish `i`, the unchanged production social torque is

```text
Omega_i^soc = (1/Z_i) sum_(j in V_i)
    [ Ip sin(theta_j-theta_i) + |rho_ij| sin(beta_ij) ] w_ij,

w_ij = 1 + cos(beta_ij),
Z_i  = sum_(j in V_i) w_ij,
rho_ij = r_j-r_i.
```

`V_i` is the wall-filtered Delaunay neighborhood and `beta_ij` is the signed
angle from heading `e_i` to `rho_ij`.  Therefore
`|rho_ij| sin(beta_ij)=e_i x rho_ij`.  The graph topology is reciprocal before
weighting, but `w_ij/Z_i` is directed.  A fish in front can have weight near 2
while the reverse, rearward view has weight near 0.

For source mean heading `h=P m`, the exact deterministic social contribution
to its angle is

```text
omega_soc = sum_(i in S) c_i Omega_i^soc,
c_i = cos(theta_i-theta_bar)/(N_s P).
```

This is the social term in the exact collective-heading equation derived in
the preceding analysis.  The numerical reconstruction matches native
`GetVision` fish torques to `1.98e-12` and the collective pair sum to
`2.36e-14` over 1,348 pre-crossing snapshots.

## Exact internal-pair identity

For an internal visible pair, let

```text
k_ij = c_i w_ij/Z_i,       k_ji = c_j w_ji/Z_j,
A_ij = Ip sin(theta_j-theta_i).
```

The pair alignment contribution is exactly

```text
T_align = A_ij (k_ij-k_ji).
```

It cancels only when the two directed coefficients match.  Writing
`a_i=c_i/Z_i` gives a second exact split:

```text
T_align = A_ij mean(a_i,a_j) (w_ij-w_ji)
        + A_ij mean(w_ij,w_ji) (a_i-a_j).
```

The first term is front--back visual asymmetry.  The second is receiver
projection/normalization asymmetry.

The positional pair torque also has the exact decomposition

```text
T_pos = mean(k_ij,k_ji) (e_i-e_j) x rho_ij
      + (k_ij-k_ji)/2 (e_i+e_j) x rho_ij.
```

These are respectively a heading-mismatch term and a directed-visibility
term.  Edges to fish outside the source school are unpaired external inputs
and are retained separately.

Pairing removes only `37.4%` `[37.2%,37.6%]` of the internal directed
`L1` torque magnitude.  Approximate polarization therefore does not imply
collective cancellation.

## Which term distinguishes commitment?

The fixed cohort contains 337 first majority-to-minority crossings from 37
original seeds.  Continued events transfer at least four additional fish over
five time units.  Every comparison below is made before the first crossing,
inside the same fixed population-imbalance, COM-distance, and relative-bearing
strata.  Uncertainty is resampled over original seeds, with seeds 30--39 held
out.

The torque terms are projected onto the reflection-invariant rate of improving
gap alignment, `-sin(psi) omega`.  On held-out seeds:

| lead | total social | internal visual-alignment asymmetry | normalization asymmetry | positional heading mismatch | positional visibility asymmetry |
|---:|---:|---:|---:|---:|---:|
| 2.0 | +0.043 `[0.000,0.084]` | +0.040 `[0.000,0.076]` | -0.000 `[-0.012,0.010]` | +0.002 `[0.000,0.004]` | -0.001 `[-0.004,0.005]` |
| 1.0 | +0.067 `[0.038,0.103]` | +0.056 `[0.026,0.090]` | +0.011 `[-0.006,0.022]` | +0.003 `[0.002,0.005]` | -0.003 `[-0.009,0.003]` |
| 0.5 | +0.026 `[-0.010,0.068]` | +0.047 `[0.011,0.088]` | -0.015 `[-0.030,-0.003]` | +0.002 `[0.000,0.004]` | +0.001 `[-0.003,0.007]` |
| 0.2 | +0.042 `[-0.013,0.084]` | +0.059 `[0.007,0.100]` | -0.008 `[-0.022,0.005]` | +0.004 `[0.001,0.005]` | -0.001 `[-0.006,0.010]` |

At lead 1, the internal visual-asymmetry term accounts for about 83% of the
held-out total social contrast.  Its training-seed contrast has the same sign
and is resolved.  The production alignment strength `Ip=4` amplifies small
`sin(theta_j-theta_i)` values, while directed visibility prevents their
antisymmetric cancellation.

External source-receiver terms do not supply the net effect at lead 1:

```text
external alignment contrast   -0.026 [-0.047,-0.009]
external positional contrast  +0.027 [ 0.007, 0.052].
```

They oppose one another almost exactly.  The commitment-specific social turn
is therefore primarily an internal source-school organization effect, not
simply cross-divider recruitment.

## Is it a leading-edge or boundary torque?

Three parameter-free receiver partitions were audited:

- leading/trailing halves are separated by the plane through the source COM
  normal to the instantaneous gap direction;
- boundary fish are vertices of the source-school convex hull;
- approach receivers use the already established physical approach rectangle.

Across held-out snapshots:

| receiver set | population share | absolute social-torque share | enrichment |
|---|---:|---:|---:|
| leading half | 0.469 | 0.531 | 1.13 |
| convex hull | 0.173 | 0.320 | 1.85 |
| approach region | 0.482 | 0.527 | 1.09 |

Boundary receivers are clearly more torque-active per fish.  They do not,
however, dominate the outcome separation.  At lead 1 on held-out seeds the
continued-minus-failed alignment contrasts are

```text
leading       +0.031 [0.009,0.056]
trailing      +0.036 [0.022,0.051]
hull          +0.001 [-0.010,0.016]
interior      +0.066 [0.041,0.092]
approach      +0.044 [0.012,0.078]
non-approach  +0.023 [0.013,0.035].
```

Thus the mechanism is distributed through the directed visual network.  The
opening-facing part is important, but neither a leader torque nor a hull-only
torque is an adequate closure.

## Controlled coherent-spatial approximation

A natural controlled approximation replaces every source heading by the mean
heading while preserving the exact positions, visible graph, and external
neighbors.  This removes internal alignment mismatch but retains positional
and directed-visibility geometry.

It correlates with the instantaneous exact social angular rate (`0.820` on
held-out seeds), but its RMSE is `0.0925`, compared with exact-rate standard
deviation `0.1415`.  More importantly, it fails the mechanism test.  At lead 1
its held-out continued-minus-failed alignment contrast is only `+0.009`
`[-0.005,0.025]`, compared with the exact `+0.067` `[0.038,0.103]`.

The approximation reproduces common spatial variation but deletes the small
heading differences that generate commitment-specific turning.  It is not a
quantitative coarse closure.

## Connection to approach storage and flux

For the existing source approach rectangle, define the exact stored-velocity
current `J_a` and weighted storage

```text
b = (1/(N ell)) sum_inside (ell-u_i).
```

The previously verified conservation law is

```text
F_neck = N J_a/ell - (1/ell) dW/dt - L_W + epsilon_path.
```

Social turning is not an additional term in this identity.  It changes the
headings that generate `J_a`; storage and lateral exchange then determine how
much of that supply becomes neck flux.

The chronology is consistent with that mechanism on held-out seeds.  After
the same `(q,d,psi)` controls:

| lead | continued-failed `J_a` | continued-failed `b` |
|---:|---:|---:|
| 2.0 | +0.005 `[-0.020,0.024]` | -0.021 `[-0.031,-0.012]` |
| 1.0 | +0.078 `[0.050,0.105]` | -0.001 `[-0.012,0.011]` |
| 0.5 | +0.093 `[0.068,0.121]` | +0.017 `[0.004,0.043]` |
| 0.2 | +0.107 `[0.072,0.127]` | +0.023 `[0.003,0.055]` |

The current separates before storage: directed supply is already larger at
lead 1, while the storage advantage appears by lead 0.5.  Within matched
geometry strata, lead-1 social alignment predicts the subsequent change in
`J_a` from lead 1 to 0.5.  The residual Spearman correlation is `0.516`
`[0.314,0.656]` on held-out seeds and `0.298` `[0.195,0.408]` on training
seeds.

This temporal association supports the physical chain from collective turn to
supply.  It is not a new conservation identity or a causal intervention.  The
exact identity begins with the supplied current and determines its conversion
into storage, lateral loss, and neck flux.

## Implication for a low-dimensional transition-rate theory

The analysis reduces the unresolved social torque to one dominant exact
structure:

```text
sum over internal visible edges of
  sin(pair heading difference)
  x front-back visual-weight difference
  x collective projection/normalization.
```

This is progress over an unspecified microscopic torque: pair symmetry,
positional turning, external input, leader dominance, and boundary dominance
have all been tested and found insufficient explanations.

However, the dominant edge covariance cannot be calculated from population
imbalance, COM distance, relative bearing, global polarization, or source COM
alone.  Exact polarization close to one is not enough because the required
small heading differences are spatially correlated with directed neighbor
geometry.  Introducing a sequence of empirical moments to encode that
covariance would violate the goal of a minimal theory.

Therefore no closed analytical transition rate is claimed.  A first-principles
closure would need to derive this edge covariance from a controlled
finite-extent polarization field or boundary-layer theory, and then couple its
predicted turn to the already exact approach-storage balance.  Until that is
done, the present result is a physically explicit mechanism and a precise
statement of the remaining closure problem, not a coarse-grained stochastic
law.

## Status of claims

- **Exact:** production pair decomposition and collective reconstruction.
- **Approximate cancellation tested:** only 37% of directed internal magnitude
  cancels on pairing.
- **Observed mechanism:** internal front--back visual weighting converts small
  heading mismatches into sustained collective turning.
- **Independent prediction:** this term and its link to later approach current
  reproduce on held-out seeds.
- **Rejected closures:** perfect-coherence spatial approximation, leader-only,
  hull-only, positional-only, and external-input-only explanations.
- **Not introduced:** fitted turning functions, new dynamical coordinates,
  transition kernels, potentials, or microscopic simulations.
