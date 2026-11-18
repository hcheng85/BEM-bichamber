# Exact transport budget of continued and failed crossing attempts

## Physical conclusion

Failed crossing attempts suffer **both loss of directed supply and lateral escape**. Neither mechanism alone is statistically dominant during the first unit after entry: their contributions to the continued--failed neck-flux separation are comparable, and the bootstrap interval for their difference includes zero. The time ordering is nevertheless informative. Successful continuation is prepared by a source-side layer that is still being filled by longitudinal current and lateral inflow. Failed attempts have already switched to storage depletion and outward lateral loss before the first fish crosses. Immediately after entry, both groups release stored population, but failed events have weak replacement current and lose much more of that population through the sides instead of the neck.

This supports a **supply--retention--release** mechanism. It refines the previous persistent-orientation--accumulation picture without turning it into a closed stochastic model.

## Fixed cohort and comparison

The analysis reuses all 337 majority-to-minority episode crossings from the `(q,a,b)` common cohort: 185 events transferred at least four additional net fish during the following five time units and 152 did not. The events come from 37 independent original seeds. “Continued” therefore denotes short five-unit propagation, not complete relocation of the school.

For every event, the source chamber and forward direction are fixed from the population immediately before first entry. They are not reassigned if the chamber populations later become equal or reverse. Exact budgets are measured in nonoverlapping 0.2-unit windows from `t=-3` to `t=+2`. The crossing segment itself belongs to the window `(-0.2,0]`.

Continued and failed curves are directly standardized to the pooled cohort frequencies in the three previously declared initial-imbalance strata, `|q_pre|` in `[0,0.5)`, `[0.5,0.8)`, and `[0.8,1]`. Within each outcome and stratum, events are first averaged within original seed and then across seeds. Confidence intervals resample original seeds as clusters. The supply, lateral, and neck-flux differences reported below have the same sign and resolved seed-bootstrap intervals in each of the three imbalance strata.

## Exact identity and sign convention

The analysis reuses the validated native-segment conservation code. With local longitudinal coordinate increasing from the original source toward the gap,

\[
F_{\rm neck}=\frac{N\bar J_a}{\ell}-A_W-L_W+\epsilon_{\rm path}.
\]

Here $J_a$ is the archived, stored-velocity approach current; $A_W=\Delta W/(\ell T)$ is the change of downstream-weighted population storage; $L_W$ is the downstream-weighted lateral crossing flux, positive outward; and $\epsilon_{\rm path}$ is the audited correction between stored left-endpoint velocity and the segment-clipped path integral. The figure plots $-A_W$, so positive values mean reservoir release, and $-L_W$, so positive values mean lateral inflow while negative values mean lateral escape.

Across 8,425 event windows, the maximum absolute identity residual is `7.82e-14` fish/time and its RMS is `9.48e-15`. The mean absolute path correction is `7.88e-4` fish/time and its maximum is `9.70e-3`, negligible relative to the physical terms. These are exact accounting results to numerical precision. The explanation of why supply or lateral exchange differs between outcome groups remains statistical and is not implied by the identity.

## Time-dependent pathway

### Early approach

From `t=-3` to `-1`, directed supply is not resolved between outcomes: continued minus failed is `-0.257` fish/time, 95% CI `[-0.548, 0.033]`. Both groups are predominantly accumulating population. Continued events have greater lateral inflow (`+1.133`, CI `[0.300, 2.022]`), but there is no neck transport. Thus the outcome is not separated by an already active neck current several units in advance.

### Last unit ending at first entry

The budgets diverge strongly during `(-1,0]`:

| contribution, fish/time | continued | failed | continued - failed, 95% CI |
|---|---:|---:|---:|
| directed supply `N J_a/ell` | 5.004 | 2.715 | 2.289 `[1.733, 2.796]` |
| storage release `-A_W` | -6.998 | 1.042 | -8.040 `[-9.463, -6.673]` |
| lateral contribution `-L_W` | 2.997 | -2.757 | 5.754 `[4.248, 7.172]` |
| net neck flux | 1.003 | 1.000 | 0.003 `[0.000, 0.011]` |

The equal neck flux is expected because every selected event contributes its defining first crossing. Internally, however, the pathways are opposite. Continued events use strong forward supply plus lateral inflow to keep filling the weighted reservoir. Failed events already release storage while losing fish laterally. The first crossing therefore does not identify the state of the surrounding supply layer.

### First unit after entry

During `(0,1]`, the neck-flux difference becomes `6.550` fish/time, CI `[5.653, 7.381]`. Its exact decomposition is:

- lower directed supply in failures: `3.407`, CI `[2.962, 3.842]`;
- greater lateral escape in failures: `4.065`, CI `[2.620, 5.410]`;
- storage-release difference: `-0.921`, CI `[-2.519, 0.674]`;
- path-correction difference: `-0.0006` fish/time.

The positive supply and lateral contributions sum to `7.472`, CI `[6.095, 8.887]`; the unresolved storage difference partly offsets them. Lateral escape has the larger point contribution, but its advantage over directed supply is only `0.658`, CI `[-0.840, 2.175]`. The data therefore resolve **both mechanisms**, but do not resolve one as dominant over the other.

Notably, failed events release at least as much weighted storage immediately after entry as continued events. Their failure is not an absence of stored population release. The released population is insufficiently replenished and is preferentially lost through the lateral boundaries.

### One to two units after entry

Continuation then becomes reservoir-supported propagation. Continued events retain more directed supply (`+2.175`, CI `[1.739, 2.596]`) and more storage release (`+3.667`, CI `[2.435, 4.880]`), while the lateral difference is unresolved (`-0.549`, CI `[-1.616, 0.527]`). Continued neck flux remains `5.096` fish/time, whereas the failed mean is `-0.195`. Thus lateral escape distinguishes the initial branch, while sustained later transport additionally requires continued depletion of the larger loaded reservoir.

## Implication for reduced theory

Population imbalance, polarization, and scalar storage do not specify whether the approach layer is being replenished or leaking sideways. A physically controlled reduced description must retain the source-exact longitudinal supply and lateral boundary renewal associated with the conservation law. The key hidden state is not another fitted order parameter: it is the resolved transport state of the approach layer, including whether lateral exchange is inward or outward.

The budget does not identify the microscopic torque that creates persistent supply or turns lateral inflow into escape. Establishing that cause would require a separate source-level or controlled-intervention analysis. The present result establishes the mesoscopic mechanism: significant continuation occurs when storage release is accompanied by maintained forward supply and short-time lateral retention.

## Reproducibility

- `analyze_event_transport_budget.py` imports the already validated `segment_budget` implementation and changes no simulator code.
- `physical_budget.png` is the single compact physical-budget figure.
- `time_resolved_summary.csv`, `period_contrasts.csv`, and `q_stratum_sensitivity.csv` contain seed-bootstrap estimates.
- `event_window_budgets.csv` records the exact native-segment terms for audit; `summary.json` records identity closure.

No new trajectory or microscopic simulation was generated.
