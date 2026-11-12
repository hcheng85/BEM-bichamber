# Event-aligned polarization and continuation after first crossing

## Question and exact accounting

This read-only analysis asks whether majority-school gap-directed polarization measured exactly 1 time unit before a previously indexed first crossing predicts subsequent propagation. It uses the existing 40 no-hydrodynamics gap-20 trajectories and changes neither simulator nor event index.

For event frame `k`, `q` and the source majority are defined from `state[k-1]`. Polarization is computed from `state[k-1000]`, using the direction from that majority school's COM to the gap center. The maximum measured lead-time error is 0.0e+00.

The follow-up endpoint is `state[k+5000]`. Additional transfer is

`s_majority * (N_R[k+5000] - N_R[k])`,

so the population change from `k-1` to `k` is explicitly excluded. The count identity closes with maximum error 0.0e+00. Of the 388 primary events, 1 had more than one net fish in the initial native timestep; that entire initial-step change is excluded. Additional net transfer ranges from -1 to 66 fish, including reversals.

## Sample

- Indexed episodes: 658.
- Population ties excluded: 9.
- First crossings opposite to the current majority-to-minority direction: 261 (counted, but excluded from the primary directional test).
- Lead-frame versus immediate-pre-crossing population mismatches: 0.
- Primary majority-to-minority crossings: 388 across 37 independent seeds.
- Continuations of at least 4 additional net fish within 5 units: 187/388.

Equal population at the lead and event endpoints does not mathematically exclude a compensating native cross-and-return between them. The existing event index and the earlier bounded native audit support the quiet-period interpretation, but this analysis does not claim an exhaustive native scan of every pre-event interval.

## Conditional continuation

At fixed population imbalance, continuation rises strongly with pre-crossing orientation:

| current abs(q) | `a<0.2` | `0.2<=a<0.6` | `a>=0.6` |
|---|---:|---:|---:|
| 0--0.5 | 6/82 = 0.073 | 56/108 = 0.519 | 47/54 = 0.870 |
| 0.5--0.8 | 1/27 = 0.037 | 32/51 = 0.627 | 16/16 = 1.000 |
| 0.8--1 | 3/11 = 0.273 | 7/16 = 0.438 | 19/23 = 0.826 |

Seed-bootstrap intervals are plotted and tabulated. The `0.5<=|q|<0.8, a>=0.6` cell is only 16 events despite its 16/16 continuation rate, so that exact value should not be generalized as certainty.

## Independent-seed prediction

Seeds 0--29 train the fixed binned probabilities; seeds 30--39 are held out. The held-out sample has 144 events and 75 continuations.

- AUC: 0.564 for q alone, 0.814 for `(q,a)`.
- Mean seed-level AUC gain: 0.192 (95% seed-bootstrap CI 0.077--0.284).
- Brier score: 0.25025 for q alone, 0.16850 for `(q,a)`.
- Mean seed-level Brier improvement: 0.07359 (95% CI 0.04364--0.10246).

The improvement is large and reproducible across independent seeds, although seed 33 is a small negative exception (five held-out events).

## Scientific conclusion

Gap-directed orientation predicts **propagation after initiation**, not initiation alone. A school already strongly directed toward the opening one unit before the first fish crosses is far more likely to carry additional fish through afterward. Because `a` is measured before the first crossing and the continuation count starts after it, this is not an after-the-fact flux signature.
The earlier initiation test increased held-out AUC from 0.619 to 0.785; the present continuation test increases it from 0.564 to 0.814. These are different conditioned outcomes and their numerical gains should not be treated as directly comparable effect sizes. Together they show that orientation is a pre-crossing **commitment coordinate** carrying information about both whether crossing starts and whether it grows.

This remains predictive rather than causal evidence. The test does not identify which wall, social, or stochastic torque generated the orientation, and it does not establish a whole-school relocation law, Markov closure, or effective potential.

## Reproducibility

- `event_aligned_samples.csv`: raw event-level accounting and timing audits.
- `conditional_continuation.csv`: counts, probabilities, and seed-bootstrap intervals.
- `heldout_validation.json`: independent-seed q versus `(q,a)` validation.
- `conditional_continuation.png`: conditional rates.
- `additional_transfer_distribution.png`: complete additional-transfer range; no events are clipped from the histogram.
- `summary.json`: sample and identity checks.
