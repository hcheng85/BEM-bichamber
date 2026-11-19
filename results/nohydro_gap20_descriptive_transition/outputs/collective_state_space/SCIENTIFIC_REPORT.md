# Minimal collective state for bichamber switching

## Decision

The archived dynamics support two long-lived, reflection-related
chamber-dominance states, but they do **not** support a scalar diffusive
reaction coordinate or a fitted double-well potential.

The smallest useful *diagram* is `(q,p_x)`, where `q=(N_R-N_L)/N` and
`p_x=N^{-1} sum cos(theta_i)` is global signed axial polarization.  It shows
the population loading and the direction of collective transport without
relabeling at a dominance reversal.  It is not dynamically closed.

The recommended next kinetic model is a reflection-symmetric hybrid process

```text
q in {-1,-0.98,...,1},       discrete population jumps,
p_x,                         persistent global axial motion,
b_L, b_R,                    left/right approach-layer storage.
```

Within a fixed dominance basin this reduces to `(q,p_x,b_s)`, with `b_s` the
current source reservoir.  Both `b_L` and `b_R` must be retained through
`q=0`, where majority-based source labels exchange.  The next model should be
a population-lattice jump process coupled to axial-motion and exact storage
budgets, not a Langevin diffusion in `q`.

## Relation to the requested coarse-graining framework

[Huang, Ling & Kanso (PNAS 2024)](https://doi.org/10.1073/pnas.2406293121)
first characterized confined-school states with global polar and milling order,
then projected onto a scalar coordinate only after the data occupied a nearly
one-dimensional curve in that order-parameter plane.  They separately showed
that `-log p_stationary` can imply a finite barrier even when dynamically
reconstructed drift and diffusion give an infinite one.

[Kolpas, Moehlis & Kevrekidis (PNAS 2007)](https://doi.org/10.1073/pnas.0608270104)
likewise required evidence for an attracting slow manifold and a meaningful
coarse observable before constructing an equation-free Fokker--Planck model.
Their lifting/healing procedure is not available here because the current task
forbids new microscopic simulations.  We therefore perform the necessary
state-space, memory, and jump audits but stop before drift/diffusion or
potential estimation.

## Physical coordinates

Post-transient microscopic states were sampled every one time unit from all 40
existing trajectories, giving 32,040 configurations.  The following quantities
were reconstructed directly from fish positions and headings:

- population imbalance `q=(N_R-N_L)/N`;
- global signed axial polarization `p_x` and transverse component `p_y`;
- global polar order `P=|(p_x,p_y)|`;
- signed milling and `Pi=|m|`, using normalized angular momentum about the
  instantaneous school COM;
- left- and right-chamber axial polarization, polar order, and milling;
- majority-source axial polarization;
- exact weighted approach storages `b_L` and `b_R`, whose source-side version
  is the previously validated `b` coordinate.

Reflection maps `(q,p_x,m,b_R-b_L)` to its negative and leaves
`(P,Pi,b_R+b_L)` unchanged.  The stationary diagram pools every configuration
with its reflected partner.  The unpooled finite ensemble has mean `q=-0.078`,
while mean `p_x=0.00012` and signed milling `-0.00138`; the left/right occupancy
difference is unresolved by seed bootstrap.  No physical left/right bias is
inferred.

## Collective-motion orders do not define a milling--schooling switch

The school is predominantly polarized:

```text
global P:  median 0.954, 10--90% range 0.590--0.981
global Pi: median 0.099, 10--90% range 0.017--0.324
corr(P,Pi) = -0.226
```

The joint `(P,Pi)` distribution is a broad wedge dominated by high polarization
and low milling, not two separated milling and schooling modes connected by a
thin curved manifold.  Chamber-resolved polar orders are also high, but
minority-chamber values are undefined when empty and biased upward when only a
few fish are present.  They are unsuitable primary landscape coordinates.

Consequently there is no data-based justification for constructing a nonlinear
scalar `K(P,Pi)` analogous to the PNAS study.  Milling is a useful phase
diagnostic but does not parameterize bichamber switching here.

## Metastability and actual relocations

Equal-seed stationary occupancy gives

```text
P(q >= 0.8) = 0.248, 95% seed CI 0.131--0.378
P(q <=-0.8) = 0.281, 95% seed CI 0.162--0.404
P(|q| < 0.5) = 0.360, 95% seed CI 0.245--0.470.
```

Thus two reflection-related dominance basins are substantially occupied, and
`q` is extremely persistent: seed-mean autocorrelation is 0.962 at lag 1,
0.894 at lag 5, and 0.724 at lag 20.  This supports metastable chamber
dominance.  The large central occupancy simultaneously shows that the system
is not accurately summarized as a two-state telegraph process.

Only nine complete passages between opposite `|q|>=0.6` basins occur, in six
of the 40 seeds.  Their `(q,p_x)` paths occupy a broad channel and do not
collapse onto a single curve.  Nine correlated transitions are insufficient
to estimate a switching rate, transition-path density, committor, or barrier.

## Slow dimension versus dynamical sufficiency

Lagged canonical correlations at lag 5 give:

| physical feature set | singular values |
|---|---|
| `q` | 0.996 |
| `(q,p_x)` | 0.996, 0.312 |
| `(q,p_x,P,Pi)` | 0.996, 0.388, 0.314, 0.165 |
| `(q,p_x,j_source,b,Pi)` | 0.996, 0.499, 0.372, 0.323, 0.310 |

There is one dominant slow coordinate, `q`.  Axial motion and storage are
faster: at lag 1 their autocorrelations are 0.825 (`p_x`) and 0.779 (`b`), and
they become oscillatory at longer lags.  The spectrum does not show a clean
two-dimensional slow manifold, and adding polar/milling order does not create
one.

Being slow is not the same as being Markov sufficient.  A held-out diagnostic
used seeds 0--29 to tabulate five-unit future population progress and seeds
30--39 to test whether residuals still correlate with the preceding two-unit
population history.  This is a closure diagnostic, not an AUC comparison:

| current state | held-out support | residual history Spearman, 95% seed CI | future-progress MSE |
|---|---:|---:|---:|
| `q` | 100% | 0.087 `[0.045,0.128]` | 0.00869 |
| `(q,p_x)` | 100% | 0.074 `[0.036,0.109]` | 0.00829 |
| `(q,b_s)` | 100% | 0.090 `[0.048,0.132]` | 0.00823 |
| `(q,p_x,Pi)` | 98.9% | 0.073 `[0.033,0.112]` | 0.00820 |
| `(q,p_x,b_s)` | 98.0% | 0.010 `[-0.022,0.039]` | 0.00608 |
| `(q,p_x,b_R-b_L,b_R+b_L)` | 93.3% | 0.004 `[-0.039,0.047]` | 0.00502 |

Neither axial motion, milling, nor storage alone closes the tested memory.
Axial motion and storage together do.  The source-only form is the smallest
tested description within one dominance basin.  The two chamber-resolved
storages are the smallest tested smooth representation through a full switch.
Passing this one history test is not proof of a Markov process; it identifies
the next state space on which stronger Chapman--Kolmogorov and renewal tests
should be performed.

## Stationary landscape and its limitation

The figure reports the reflection-pooled empirical occupancy in `(q,p_x)`.
Only cells containing at least 20 samples from at least five seeds are shown:
433 of 1,353 cells pass.  The displayed quantity is

```text
-log(empirical stationary probability mass) + constant.
```

It is a preliminary statistical landscape only.  It is not a coordinate-free
density, a free energy, or the potential of a reconstructed stochastic
dynamics.  In particular:

- `(q,p_x)` retains measurable history dependence;
- the dynamically better state includes storage and is at least three
  dimensional within a basin;
- sampling is sparse along complete transition paths;
- nine relocations cannot establish a finite barrier.

No barrier height or Kramers switching rate is reported.

## Population imbalance is a jump variable

At the native `dt=0.001`, 99.9875% of the sampled increments have `Delta q=0`;
every observed nonzero increment is exactly one fish, `|Delta q|=0.02`.  The
fourth-moment ratio to a Gaussian increment is 2,670.  Even at lag 0.1,
98.46% of increments remain zero and the ratio is 62.9.  The estimated second
Kramers--Moyal coefficient changes by more than an order of magnitude between
lags 0.001 and 2 rather than approaching a lag-independent diffusion value.

These facts reject a diffusion approximation for `q` at the available
resolution.  Multiple crossings accumulated over coarse windows do not change
the microscopic population-lattice character.

## Recommended next model

Develop a reflection-symmetric continuous-time jump/transport model:

```text
(q,p_x,b_L,b_R) -> (q +/- 2/N, p_x', b_L', b_R')
```

with chamber-crossing intensities conditioned on the current physical state.
Between jumps, use the already derived polarization and exact storage
evolution/conservation identities to constrain the continuous transport
variables.  Estimate jump intensities and deterministic/stochastic transport
only after checking lag convergence and Chapman--Kolmogorov consistency on
independent seeds.  A source-canonical `(q,p_x,b_s)` representation may be used
inside a dominance basin, but a full switching theory must keep both reservoirs
at ties and reversals.

Do not proceed next to a fitted potential.  The next scientifically controlled
step is a generator/renewal audit for this hybrid state, followed—only if it
closes—by equation-free short-burst calculations or a master equation for
switching kinetics.

## Reproducibility

`analyze_collective_state_space.py` reads only the archived trajectories and
existing occupancy caches.  It produces the single state diagram, supported
stationary cells, actual relocation paths, time correlations, held-out memory
diagnostics, and jump audit.  No microscopic simulation, nonlinear latent
coordinate, fitted potential, or parameter sweep was introduced.
