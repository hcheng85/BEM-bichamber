# Independent confirmation of fixed-`q` committor heterogeneity

## Scientific decision

The independent experiment rejects the hypothesis that the four microscopic
configurations have the same eventual forward committor, despite all having
exactly `q=-0.2`.

All 64 new shots resolved by time 1200. The prespecified exact conditional
4x2 test of equal B-before-A probabilities gives **p=0.0171**. The new
committor estimates span 0.125--0.5625, a difference of 0.4375. Population
imbalance therefore does not uniquely determine the transition state for
these configurations: microscopic school organization or motion retained in
the full state affects the probability of eventual relocation.

This result establishes committor insufficiency of `q` for the tested states.
It does not identify which additional physical variable is responsible, and
it does not by itself establish Markovian or kinetic closure in any enlarged
coordinate set.

## Fixed design and restart validation

The original archive configurations were retained without reselection:

| archive seed | state ID | initial q |
|---:|---|---:|
| 8 | `q-0.2_s8_k393600` | -0.2 |
| 12 | `q-0.2_s12_k296500` | -0.2 |
| 14 | `q-0.2_s14_k263900` | -0.2 |
| 15 | `q-0.2_s15_k583100` | -0.2 |

Each state received exactly 16 new noise realizations, clone IDs 8--23. The
64 seeds use the fixed namespace
`202611270000 + 100*archive_seed + clone`; uniqueness and nonoverlap with the
original clone-0--7 seeds were verified before shooting.

The source-extracted production dynamics were unchanged. Previous-step
wall/collision memory was reconstructed and included in every restart. Native
one-step replay passed with maximum error `2.22e-16`, and the reflection audit
passed with maximum error `1.41e-14`. The basin definitions remained
`A: q<=-0.6` and `B: q>=+0.6`, checked after every native `dt=0.001` step.

Checkpoints atomically store the full microscopic state, complete collision
memory, and RNG state every 25 time units. The 1200-stage run continued only
unresolved 600-stage paths; it did not replay their earlier evolution.

## Censoring-aware 600-unit result

At time 600, 10 of 64 new shots remained unresolved:

| seed | A | B | censored | empirical eventual-committor bounds |
|---:|---:|---:|---:|---:|
| 8 | 7 | 5 | 4 | 0.3125--0.5625 |
| 12 | 7 | 7 | 2 | 0.4375--0.5625 |
| 14 | 11 | 4 | 1 | 0.2500--0.3125 |
| 15 | 11 | 2 | 3 | 0.1250--0.3125 |

Resolved-only proportions were not interpreted as committors. The
prespecified 100,000-permutation test of equality of the finite-horizon
A/B/censored distributions gave `p=0.3677`. Exact secondary tests gave
`p=0.2875` for B cumulative incidence and `p=0.6639` for unresolved survival.
Thus the 600-unit experiment alone did not resolve state heterogeneity.

## Complete independent confirmation at 1200

Every censored path subsequently hit A or B; the latest first hit occurred at
time 1083.699.

| seed | A | B | new committor | exact 95% binomial interval |
|---:|---:|---:|---:|---:|
| 8 | 7 | 9 | 0.5625 | 0.2988--0.8025 |
| 12 | 7 | 9 | 0.5625 | 0.2988--0.8025 |
| 14 | 12 | 4 | 0.2500 | 0.0727--0.5238 |
| 15 | 14 | 2 | 0.1250 | 0.0155--0.3835 |

The global exact test conditions on 24 total B hits and enumerates all 2,445
4x2 tables with 16 shots per state. Its two-sided probability-ordering
p-value is 0.01705. This is the confirmatory result; it uses only the new
clone-8--23 shots.

The evidence supports a broad separation: the seed-8 and seed-12
configurations have higher B commitment than seed 14 and especially seed 15.
The original point ordering between seeds 8 and 12 did not reproduce exactly,
so the defensible conclusion is global fixed-`q` heterogeneity, not precise
ranking of every pair.

## Historical plus new shots: descriptive only

Combining the eight historical shots with the 16 new shots gives:

| seed | B / 24 | combined estimate | exact 95% binomial interval |
|---:|---:|---:|---:|
| 8 | 12 / 24 | 0.5000 | 0.2912--0.7088 |
| 12 | 16 / 24 | 0.6667 | 0.4468--0.8437 |
| 14 | 5 / 24 | 0.2083 | 0.0713--0.4215 |
| 15 | 3 / 24 | 0.1250 | 0.0266--0.3236 |

These 24-shot estimates improve descriptive precision, but they are not
independent confirmation because the original outcomes motivated the new
experiment. They are excluded from the confirmatory p-value.

## Interpretation

The result answers the narrow question decisively: `q=-0.2` identifies a
population count but not a unique progress probability. Two configurations at
that same count have estimated forward committors near 0.56, while another is
near 0.13. A one-dimensional committor `h(q)` cannot represent all four full
states within shooting uncertainty.

No new reaction coordinate is introduced here. Determining whether a specific
collective-motion variable explains this heterogeneity is a separate
hypothesis and is not inferred from these shootings. Likewise, committor
heterogeneity does not imply a fitted energy landscape or specify a transition
rate.

## Exact local commands

Run from `results/nohydro_gap20_descriptive_transition` in the repository:

```bash
uv run --with numpy --with scipy --with matplotlib \
  python run_qminus02_committor_confirmation.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python run_qminus02_committor_confirmation.py --stage audit

uv run --with numpy --with scipy --with matplotlib \
  python run_qminus02_committor_confirmation.py \
  --stage shoot --target-time 600 --checkpoint-interval 25 --workers 20

uv run --with numpy --with scipy --with matplotlib \
  python run_qminus02_committor_confirmation.py \
  --stage shoot --target-time 1200 --checkpoint-interval 25 --workers 10

uv run --with numpy --with scipy --with matplotlib \
  python run_qminus02_committor_confirmation.py \
  --stage analyze --target-time 1200
```

Rerunning either shooting command is safe: resolved checkpoints return
immediately, unresolved checkpoints continue from their stored RNG and
microscopic states, and provenance mismatches cause a hard failure. The local
`cache/` directory contains checkpoints and is excluded from GitHub.
