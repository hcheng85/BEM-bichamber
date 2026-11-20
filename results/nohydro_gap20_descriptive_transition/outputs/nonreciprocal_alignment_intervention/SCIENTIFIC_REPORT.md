# Causal suppression of nonreciprocal visual alignment

## Answers

**Is nonreciprocal visual alignment necessary for the observed collective
turning?**  It is necessary for the full turning amplitude, but it is not the
only source of turning.  Surgical symmetrization reduces first-unit turning
toward the opening by 44%, while residual wall, positional, external-social,
noise, and initial-condition effects retain a positive turn.

**Does suppression change sustained transfer probability?**  It substantially
reduces the *amount* of transfer but does not resolve a change in the binary
four-fish probability.  Mean five-unit transfer falls by `4.44` fish
(`-7.13,-1.76`), whereas the `>=4`-fish probability changes by only `-0.009`
(`-0.161,0.161`).  Thirteen paired clones lose threshold success and twelve
gain it.  Positive net transfer must therefore not be equated with a resolved
change in threshold commitment probability.

**Does this support a minimal physical theory?**  It establishes a causal
mechanistic link

```text
directed internal visual alignment
  -> gap-directed collective turning
  -> approach current and storage
  -> amount of population transfer.
```

It does not yet close a transition-rate theory because the directed pair
coefficient still depends on the evolving visible graph and small spatially
organized heading differences.  A coarse theory may include this torque only
after deriving that edge covariance; fitting it from global coordinates would
not be first-principles closure.

## Exact instantaneous counterfactual

For an internal source pair, the exact collective alignment contribution is

```text
T_ij = A_ij (k_ij-k_ji),
A_ij = Ip sin(theta_j-theta_i),
k_ij = c_i alpha_ij,
alpha_ij = w_ij/Z_i,
c_i = cos(theta_i-theta_bar)/(N_s P).
```

Replacing both effective coefficients by their symmetric part,

```text
k_ij -> (k_ij+k_ji)/2,
k_ji -> (k_ij+k_ji)/2,
```

cancels `T_ij` exactly while preserving the coefficient mean.  The analytical
and numerical residual is identically zero for all 337 established lead-1
snapshots.

This exact projected counterfactual is not a clean microscopic dynamics.
Implementing it fish by fish requires division by the instantaneous `c_i`, so
the pair rule depends on the global source mean heading and becomes singular
when a fish is perpendicular to that heading.  It would introduce artificial
global feedback that is absent from the production model.

## Feasible microscopic intervention

The implemented counterfactual instead symmetrizes the actual directed
fish-level alignment coefficients:

```text
alpha_ij, alpha_ji -> alpha_bar = (alpha_ij+alpha_ji)/2
```

for every currently internal pair in the original source chamber.  Fish `i`
then receives `alpha_bar A_ij` and fish `j` receives `-alpha_bar A_ij`.

At that instant the intervention preserves exactly:

- the internal positional social numerator;
- all external alignment and positional interactions;
- each original full-neighborhood visual denominator;
- the Delaunay and wall-visibility graph;
- the wall torque input and both production clipping rules;
- the collision rule and angular-noise realization.

The modified social input can change the numerical value produced by the final
total-torque clip if a state reaches that limit; the clipping operation itself
is unchanged.  Likewise, later collision outcomes can differ after trajectories
separate even though the collision rule is identical.

Because `c_i` and `c_j` need not be identical, the collective residual is
`alpha_bar A_ij(c_i-c_j)`.  Across all 337 lead-1 snapshots its mean absolute
value is `0.00566`, versus `0.08815` for the original internal collective
alignment: only 6.4% remains.  The intervention-delta reconstruction closes to
`2.1e-16`.

Once the arms evolve, their positions, visual graphs, collisions, and all
later torques can differ dynamically.  These are consequences of the surgical
change, not additional imposed modifications.

## Restart and experimental design

The production source remains read-only.  An AST-derived restart adapter
replaces only the explicit noise input and the `GetVision` call.  Baseline
one-step parity against archived states and velocities passes on all 14 initial
configurations, with maximum error `3.3e-16`.

The bounded experiment uses:

- 14 distinct original seeds;
- the first eligible fixed-cohort crossing configuration per seed, restarted
  one time unit before archived entry;
- selection independent of the archived continuation label;
- eight noise clones per microscopic state;
- identical initial configurations and noise arrays in baseline and symmetric
  arms;
- five time units per branch.

The cohort is crossing-conditioned.  It tests the fate of established lead-1
approaches, not the unconditional frequency with which arbitrary chamber
states initiate an approach.

Uncertainty below resamples the 14 original seeds.  Clones are paired
realizations, not independent bootstrap units.

## Causal effects on organization and supply

During the first time unit, symmetric minus baseline effects are:

| quantity | paired effect | 95% seed interval |
|---|---:|---:|
| turning contribution toward gap | -0.0285 | `[-0.0463,-0.0095]` |
| geometric alignment contribution | -0.0071 | `[-0.0115,-0.0023]` |
| mean gap-directed alignment | -0.0146 | `[-0.0247,-0.0043]` |
| mean approach current `J_a` | -0.0118 | `[-0.0187,-0.0049]` |
| mean weighted storage `b` | -0.00147 | `[-0.00230,-0.00068]` |
| mean polarization magnitude | +0.00247 | `[+0.00104,+0.00405]` |

The baseline turning contribution is `0.0652`; it falls to `0.0367`.  Thus
other interactions compensate for more than half of the turn, but not for the
removed nonreciprocal component.

The polarization increase is physically important.  Reciprocal alignment
makes headings slightly more parallel, yet the school turns less effectively
toward the opening.  Collective coherence and collective directional response
are therefore distinct.  This is the causal counterpart of the earlier
observational result that high polarization alone does not determine
commitment.

The reduction persists over five units:

```text
turning-alignment rate difference  -0.00651 [-0.01388,-0.00035]
mean approach-current difference   -0.01224 [-0.02172,-0.00294]
mean gap-alignment difference      -0.03812 [-0.06543,-0.01419]
mean polarization difference       +0.00431 [+0.00192,+0.00714].
```

The five-unit storage difference is unresolved (`-0.00109`,
`[-0.00313,0.00124]`).  Storage is a transient reservoir and is also depleted
by crossing; its endpoint-average need not remain separated even when supply
and integrated transfer differ.

## Entry, commitment threshold, and transfer amount

| outcome | baseline | symmetric | paired effect, 95% interval |
|---|---:|---:|---:|
| first entry probability | 0.839 | 0.875 | +0.036 `[-0.009,0.098]` |
| `>=4` net fish probability | 0.670 | 0.661 | -0.009 `[-0.161,0.161]` |
| mean net fish transferred | 10.60 | 6.16 | -4.44 `[-7.13,-1.76]` |

First entry is not reduced.  These initial states are already established
approaches, so geometry, momentum, wall turning, and noise can still carry the
first fish across.  The intervention primarily affects sustained directed
supply and the amount of subsequent propagation.

The four-fish indicator is insensitive here: 62 clone pairs succeed in both
arms, 25 fail in both, 13 lose success, and 12 gain it (paired discordant exact
`p=1`).  This does not contradict the resolved transfer reduction.  It shows
that a coarse threshold discards much of the causal transport response and is
strongly affected by trajectory-level compensation.

## Connection to exact transport conservation

The intervention does not alter the established identity

```text
F_neck = N J_a/ell - (1/ell)dW/dt - L_W + epsilon_path.
```

It acts upstream by reducing the collective turn that generates longitudinal
approach current.  The resolved decreases in `J_a`, early storage, and net
transfer are the expected ordered consequences.  Lateral exchange, storage
depletion, and later torque compensation determine why the binary threshold
does not change monotonically clone by clone.

## Implication for theory

The causal result supports a minimal *mechanism* with three layers:

1. a nonreciprocal internal visual-alignment torque sustains source-school
   turning toward the opening;
2. coherent translation converts that turn into directed approach current;
3. exact storage and flux conservation converts supply into neck transfer.

Only layers 2 and 3 currently have low-dimensional exact closure.  Layer 1
still requires

```text
sum_internal_edges
  sin(pair heading difference)
  x directed front-back visual-weight difference.
```

The intervention establishes that this term is causal, not merely correlated.
It does not make the term computable from `(q,d,psi)` or another existing global
state.  A predictive transition-rate theory therefore still requires a
controlled derivation of this directed edge covariance, rather than a fitted
turning function or additional empirical hierarchy.

No production source was modified and no parameter sweep was launched.  Only
the bounded five-unit restart branches described above were generated; no new
full-length microscopic trajectory ensemble was run or stored.
