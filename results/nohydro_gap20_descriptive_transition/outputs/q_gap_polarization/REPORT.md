# Population imbalance and majority-school gap-directed polarization

## Scope and definitions

This read-only analysis uses all 40 existing no-hydrodynamics gap-20 trajectories (`N=100`, `In=0.5`). No simulator code, raw data, dynamical equation, fitted potential, or parameter was changed.

At each prospective risk origin,

`q=(N_R-N_L)/N`.

The **majority school** is the fish in the currently more populated chamber. Its center of mass is their arithmetic mean position `C_majority`. The instantaneous direction to the center of the opening is `g=-C_majority/|C_majority|`, and

`a = mean(e(theta_i) dot g)` over majority-school fish.

Thus positive `a` means that the majority school points toward the gap. This definition uses only the current state; it is not signed by the later outcome. Origins are 5 time units apart, have constant coarse occupancy for the previous 2 units, and have no indexed initiation in that retrospective interval. The label is whether the first indexed crossing in `(t,t+5]` is from the current majority toward the minority chamber. Outcome windows do not overlap.

## Definition and label validation

- 228 raw snapshots reproduce cached q to 1.11e-16; majority counts have 0 mismatches.
- Every majority COM is on the expected chamber side (0 mismatches), and the direction vector is normalized to 2.22e-16.
- Vector projection and the independent scalar formula `mean(cos(theta-alpha_gap))` agree to 1.22e-15.
- All 658 existing native-refined event records show a raw population change at the stored frame, with 658/658 direction matches and zero timing error.
- In a bounded native-resolution audit of 40 risk windows, there were 0 hidden crossings during the nominally quiet retrospective intervals and 40/40 prospective label matches.

The last check is sampled rather than exhaustive, so extremely rare rapid cross-and-return events could still contaminate some nominal negatives. It finds no evidence that this drives the result.

## The three empirical distributions

### 1. `p(q)`

The population distribution is strongly concentrated near `q=+/-1`, while also retaining a broad central component. The bands resample independent seeds. Population imbalance is therefore a principal physical state variable, but it does not uniquely determine whether crossing starts.

### 2. `p(q,a)`

This is the empirical distribution over the 5792 prospective, crossing-free-history risk origins, not a stationary potential. Both toward-gap and away-from-gap headings occur at the same q. Therefore a contains configuration information not determined by chamber population alone.

### 3. Conditional first-crossing probability

The result is ordered in the physically expected direction. For `|q|<0.5`, the probability rises from 0.011 (4/361) for `a` in `[-1,-0.6)` to 0.301 (90/299) for `a` in `[0.6,1]`. For `0.5<=|q|<0.8`, it rises from 0.009 to 0.353. Even at `|q|>=0.8`, where initiation is rare, the corresponding values are 0.000 and 0.053. Seed-cluster intervals for every cell are in `conditional_first_crossing.csv`.

There are 5792 nonoverlapping risk origins, 343 positives, and 37 positive seeds. These are prospective first-entry windows, not complete school relocations.

## Held-out seeds

Seeds 0--29 define fixed binned probabilities; seeds 30--39 are held out. The held-out set contains 1370 origins and 129 positives. Pooled AUC increases from 0.619 for q alone to 0.785 for `(q,a)`. The mean seed-level AUC improvement is 0.191, with seed-bootstrap 95% CI 0.159--0.227.

Pooled Brier score decreases from 0.08399 to 0.07616. The mean seed-level improvement is 0.00826 (95% CI 0.00413--0.01270), about 9.6% of the mean q-only seed score. The gain is therefore reproducible and practically visible, not merely a tiny numerical change.

## Physical conclusion

Positive majority-school gap-directed polarization genuinely predicts first crossing **before any measured population redistribution**. The association survives conditioning on coarse population imbalance and independent-seed testing. The supported interpretation is: q determines how the population is loaded, while a measures whether the loaded school is oriented so that self-propulsion carries it toward the opening.

This is predictive evidence, not proof that polarization is the microscopic cause of the turn. In the no-hydrodynamics simulator, heading toward the gap is also kinematically tied to source-side self-propelled motion. The analysis does not identify which wall, social, or stochastic angular term created that heading, and first-fish entry is not the same as sustained collective relocation.

## Reproducibility

- `analyze_q_gap_polarization.py`: complete read-only analysis.
- `definition_validation_samples.csv` and `definition_and_label_validation.json`: raw definition and label audits.
- `prospective_q_a_windows.csv`: compact risk set.
- `conditional_first_crossing.csv`: counts, probabilities, and seed-cluster intervals.
- `heldout_validation.json`: independent-seed comparison.
- `figures/`: the three requested plots.
