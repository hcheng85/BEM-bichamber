# Prospective entry and propagation intervention

## Decision

The bounded prospective experiment supports a causal division between first
entry and sustained transfer, with an important qualification.  Internal
visual nonreciprocity is **not resolved as necessary for binary first entry**:
symmetrization changes entry probability by exactly `0.000` under both
production and zero post-restart noise.  Suppressing angular noise also does
not reduce entry.  First entry can therefore be produced by the encounter
geometry and the remaining deterministic dynamics; it is not a simple
noise-activation event and it does not require the antisymmetric internal
alignment term in every attempt.

Propagation is different.  Under production noise, internal-alignment
symmetrization reduces complete-excursion net transfer by `3.125` fish
`[-8.215,-0.205]` and the probability of at least five net fish by `0.116`
`[0.018,0.250]`.  Removing noise from the original dynamics reduces the
five-fish probability by `0.098` `[0.009,0.223]`, while its `-1.813`-fish mean
transfer effect remains statistically unresolved.  The interaction is
physical: the symmetrization effect on transfer is `1.696` fish stronger with
production noise than without it `[0.018,4.206]`.  Individual angular
exploration and directed visual alignment therefore cooperate in propagation;
neither should be represented as diffusion of the school mean heading.

This establishes a minimal causal mechanism, not a closed analytical rate:

```text
recurrent finite-school encounter
  -> first fish enters through moving geometric aperture
  -> individual angular exploration + directed internal alignment
  -> sustained approach supply
  -> storage/neck-flux conservation converts supply to population transfer.
```

## Corrected angular convention

The aperture angles in `analyze_individual_heading_geometry.py` are local to
the inward-facing frame.  If `s=-1` for a left source and `s=+1` for a right
source, a local required increment must be converted before comparison with a
global heading change or global torque:

\[
  \Delta\Theta_{\rm required}^{\rm global}
    = \operatorname{wrap}[-s\,\Delta\phi_{\rm required}^{\rm local}].
\]

The corrected 69-event held-out missed-entry subset gives:

| corrected observational result | estimate | 95% source-seed interval |
|---|---:|---:|
| turns in required global direction | 97.1% | `[94.4%,100%]` |
| reaches the initial geometric angle | 47.8% | `[36.4%,59.1%]` |
| actual turn toward requirement | 0.483 rad | `[0.413,0.552]` |
| lead-1 wall rate toward requirement | -0.0160 | `[-0.0297,-0.0062]` |
| lead-1 social rate toward requirement | +0.1360 | `[+0.1225,+0.1509]` |

Left/right reflection reverses the required global increment, the global
heading change, and every global torque, leaving their projection toward the
opening invariant.  The self-test checks this identity numerically.  The
minimum-turn calculation now compares wrapped distances to both aperture
endpoints; a separate test crosses the `+/-pi` seam.  Thus the revised result
is not produced by an unwrapped endpoint comparison.

These corrected associations show that social turning has the required sign
near archived entry.  They remain event-conditioned observations; causality is
provided only by the restart experiment below.

## Prospective causal design

The experiment selects the first majority-source approach-storage onset after
the transient in each of 14 distinct original seeds.  An onset is the first
0.5-time sample with positive established weighted approach storage after an
empty sample.  Direction and majority status use only the present and previous
samples.  No archived crossing or later outcome enters selection.

Each microscopic state has eight clones and four arms:

| interaction rule | post-restart angular noise |
|---|---|
| production directed visual alignment | production noise |
| symmetrized internal alignment | the identical production-noise stream |
| production directed visual alignment | zero |
| symmetrized internal alignment | zero |

The symmetrization is the previously audited intervention: only the directed
internal alignment coefficients of the original source school are replaced by
their pairwise symmetric part.  Wall, positional-social, external-social,
collision, visibility, normalization, and torque-clipping rules are unchanged.
As arms evolve, their geometries and later interactions can differ; those are
consequences of the intervention.

Every branch ends at the first 0.5-grid return of *total* weighted approach
storage to exactly zero.  All 448 branches complete before the 20-time safety
limit, so no fixed future window combines distinct excursions.  One-step
baseline replay reproduces archived states and velocities with maximum error
`3.33e-16`.  Uncertainty resamples the 14 source seeds; the eight clones remain
paired and are not treated as independent bootstrap units.

For every branch, direction-aligned population change equals forward minus
backward crossing counts exactly.  The zero-noise clones are identical within
each state and interaction arm, providing an additional determinism audit.

## Four-arm results

Arm estimates and intervals are source-seed clustered:

| outcome | original + noise | symmetric + noise | original + zero | symmetric + zero |
|---|---:|---:|---:|---:|
| first-entry probability | 0.241 `[0.071,0.464]` | 0.241 `[0.089,0.429]` | 0.286 `[0.071,0.500]` | 0.286 `[0.071,0.571]` |
| excursion-end net transfer | 3.527 `[0.232,9.098]` | 0.402 `[-0.071,1.045]` | 1.714 `[0,4.714]` | 0.286 `[0,0.643]` |
| probability of at least five net fish | 0.170 `[0.036,0.348]` | 0.054 `[0,0.143]` | 0.071 `[0,0.214]` | 0.000 |
| mean longitudinal approach current | 0.0223 `[0.0077,0.0404]` | 0.0164 `[0.0078,0.0248]` | 0.0186 `[0.0032,0.0358]` | 0.0200 `[0.0101,0.0304]` |
| integrated approach current | 0.148 `[0.044,0.298]` | 0.096 `[0.044,0.148]` | 0.124 `[0.024,0.243]` | 0.120 `[0.059,0.183]` |

The paired causal contrasts are:

| intervention (`modified - reference`) | first entry | net transfer | `>=5` fish | integrated current |
|---|---:|---:|---:|---:|
| symmetrize, production noise | 0.000 `[-0.098,0.107]` | -3.125 `[-8.215,-0.205]` | -0.116 `[-0.250,-0.018]` | -0.052 `[-0.176,0.030]` |
| symmetrize, zero noise | 0.000 `[-0.214,0.214]` | -1.429 `[-4.143,0.071]` | -0.071 `[-0.214,0]` | -0.004 `[-0.089,0.074]` |
| suppress noise, original | +0.045 `[-0.018,0.143]` | -1.813 `[-4.607,0.107]` | -0.098 `[-0.223,-0.009]` | -0.024 `[-0.075,0.022]` |
| suppress noise, symmetric | +0.045 `[-0.071,0.170]` | -0.116 `[-0.473,0.179]` | -0.054 `[-0.143,0]` | +0.024 `[0.004,0.044]` |

The zero entry effect is an average causal effect, not clonewise identity.  With
production noise, symmetrization removes entry in five pairs and creates entry
in five others; without noise the discordant counts are eight and eight.
Thus other deterministic interactions and the evolving geometry compensate at
the probability level, even though individual microscopic outcomes change.

The symmetry-by-noise interaction is `-1.696` fish
`[-4.206,-0.018]` for transfer and `-0.0480` `[-0.0963,-0.0065]`
for integrated current.  These interaction results explain why symmetrization
has its clearest effect under the production fluctuations: angular noise
continually creates individual heading differences on which the directed
visual rule can act to organize net supply.

The propagation response is heterogeneous, as expected for rare collective
transfer.  For symmetrization with production noise, six source seeds have a
negative mean transfer effect, eight have zero effect, and none has a positive
effect; the two most responsive starts contribute 91% of the mean magnitude.
Dropping the single largest response changes the effect from -3.125 to -0.875
fish, with a seed-bootstrap interval `[-2.106,-0.106]`.  The causal direction
is therefore not a single-seed sign reversal, but its magnitude is strongly
encounter-state dependent.  Eight original/noise branches transfer at least
20 net fish (all from one highly responsive start), versus none after
symmetrization.  The largest transfer is 54 of the original 79 source fish:
it reverses chamber dominance but is not a complete whole-school relocation.

## Entry time and collective turning

Among clone pairs that enter in both interaction arms, symmetrization delays
entry by `0.274` `[0.100,0.451]` with production noise and by `0.228`
`[0.042,0.341]` without it.  This is a post-treatment concordant-entry subset
(22 pairs from five seeds and 24 pairs from three seeds); it cannot be read as
an unconditional entry-time effect.  It suggests that nonreciprocity can
accelerate an entry that still occurs under symmetric dynamics, while the
unconditional probability remains unchanged at this sample size.

Whole-excursion source-school angular change toward the *initial* geometric
requirement does not decrease under symmetrization: the production-noise
effect is +0.026 rad `[-0.0003,0.0613]`.  Turn-to-entry changes among concordant
entries by -0.083 rad `[-0.192,0.001]`.  Neither is resolved.  This does not
negate the source-exact near-entry torque result.  The initial aperture bearing
rotates as the extended school translates and deforms, and a whole-excursion
global heading displacement is not the time integral of torque toward that
moving bearing.  The causal transport result is carried by crossing and exact
current/storage outcomes, not by forcing this global angle to act as a new
state coordinate.

## Noise scale and physical interpretation

The production angular amplitude for an individual fish is `0.5`
rad/sqrt(time).  For the strongly polarized selected source schools (median
polarization `0.973`), the exact instantaneous mean-heading scale is only
`0.0553` rad/sqrt(time), a median factor of `9.05` smaller.  Setting the
post-restart increments to zero removes individual exploration and all of its
subsequent interaction-mediated consequences; it is not equivalent to merely
setting a scalar collective diffusivity to zero.

Zero noise does not suppress first entry, so a theory in which independent
angular diffusion alone activates the first crossing is rejected.  Yet zero
noise reduces the five-fish outcome under original interactions, and removing
nonreciprocity has its largest transfer effect when noise remains.  Noise is
therefore better interpreted as a source of microscopic angular diversity
that the directed visual network converts into organized propagation.

## Minimal analytical structure and remaining limitation

The evidence justifies a two-stage encounter description without adding a new
phenomenological coordinate:

1. **Entry first passage.**  The finite school supplies many moving
   fish-specific aperture intervals.  Entry is their first hit under the
   production deterministic wall/social dynamics plus individual angular
   diffusion.  The zero-noise result proves that the deterministic part alone
   has nonzero capture probability for this prospective ensemble.  The present
   experiment does not close its turning field from coarse geometry.
2. **Conditional propagation.**  After entry, directed internal visual
   alignment and individual angular diversity sustain longitudinal supply.
   The already exact approach-storage/lateral-exchange/neck-flux identity then
   converts that supply into population transfer.  Symmetrization causally
   reduces transfer without changing the binary entry probability.

Combined with the independently observed recurrent encounter clock, a renewal
rate could have the form `encounter frequency x entry probability x
propagation reward distribution`.  Only the transport identity and geometric
aperture are currently first-principles.  The encounter frequency remains
measured, and the evolving directed visual graph prevents a closed analytical
entry/propagation kernel.  Fitting that missing response would violate the
goal; deriving it from the production interaction law remains the necessary
theoretical step.

The causal division is thus supported but not absolute: nonreciprocity is not
shown to control whether the first fish ever enters, while it contributes
causally to the full population-transfer response under production noise.  The finite
entry intervals remain wide enough that a smaller entry effect cannot be
excluded.

## Reproduction

```bash
uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/analyze_individual_heading_geometry.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/run_prospective_entry_intervention.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/run_prospective_entry_intervention.py --stage all
```

The production simulator is unchanged.  No full trajectory, fitted turning
law, added collective coordinate, potential, return map, or parameter sweep
was introduced.
