# Prespecified fixed-q committor confirmation protocol

This protocol was written before inspecting any outcome from the new shots.

## Fixed design

- Microscopic configurations: the previously selected `q=-0.2` states from
  archive seeds 8, 12, 14, and 15. No new state selection is performed.
- Confirmatory sample: clone IDs 8--23, exactly 16 new angular-noise streams
  per state (64 shots total).
- Noise namespace: `202611270000 + 100*archive_seed + clone`, verified disjoint
  from the original pilot seeds.
- Dynamics: the source-extracted production no-hydrodynamics gap-20 step,
  including reconstructed previous-step wall/collision memory.
- Basins: A is `q<=-0.6`; B is `q>=+0.6`. Each trajectory terminates at its
  first native-step hit.
- Administrative horizons: 600 first, then continuation of unresolved shots
  to 1200. Censored shots are never classified as A.
- Checkpoint interval: 25 time units, storing microscopic state, complete
  collision memory, and exact random-generator state atomically.

## Confirmatory estimand and tests

The scientific estimand is the eventual B-before-A committor for each fixed
microscopic state. If all 64 new shots resolve, the primary test is the exact
conditional Fisher--Freeman--Halton 4x2 test of equal B probabilities. The
two-sided probability-ordering p-value is enumerated exactly conditional on
the total number of B outcomes and 16 shots in each state.

If any shots remain censored, resolved-only B fractions are not treated as
committor estimates. The primary reported test is instead a prespecified
100,000-permutation Pearson test of equality of the finite-horizon
A/B/censored distributions. Exact conditional 4x2 tests of B cumulative
incidence and unresolved survival are reported as secondary finite-horizon
diagnostics. These do not by themselves test equality of eventual committors.

For each state, empirical eventual-committor bounds are `B/16` to
`(B+censored)/16`. Conservative simultaneous confidence bounds use one-sided
Clopper--Pearson limits with Bonferroni control over the eight bound endpoints.
Ordinary 95% binomial intervals are reported only when a state is fully
resolved.

## Historical results

The original clone IDs 0--7 are excluded from every confirmatory test. After
the new analysis, they are combined with the new shots only in a separately
labelled descriptive 24-shot summary. The historical results motivated this
experiment and therefore are not independent confirmatory evidence.
