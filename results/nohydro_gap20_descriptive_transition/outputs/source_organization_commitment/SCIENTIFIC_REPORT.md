# Source-chamber organization, approach, and commitment

## Physical conclusion

The data support a **two-stage organization → approach → commitment pathway**,
but not a transition from a polarized school into a persistent milling state.

1. Relative to exact-population no-crossing controls, future crossing episodes
   are already more polarized and more horizontally aligned toward the opening
   one to two time units before entry.  This is the organization of an
   **approach**.
2. Among schools that all subsequently make a first majority-to-minority
   crossing, sustained transfer is distinguished only during the final
   approximately 0.5 time unit.  Successful episodes develop stronger
   gap-directed alignment and a larger normalized milling component while
   retaining slightly *less* total polarization.  This is a late curved
   reorientation of an already highly polarized school, not loss of collective
   order.
3. Signed milling is less persistent in successful episodes than in failures.
   A stable, persistent circulation phase is therefore ruled out as the
   commitment mechanism.  The increased `|M|` is better interpreted as
   transient curved organization while the school turns and funnels toward the
   neck.

The minimal mechanism supported by these trajectories is thus:

```text
polarized source school
    -> horizontal approach toward the gap
    -> late curved/gap-aligned organization
    -> either reservoir-supported continuation or retreat.
```

This supplies an upstream collective-organization mechanism for the earlier
transport result: commitment requires a source school that remains directed
toward the opening while the approach reservoir is being supplied.  It does
not replace the exact approach-storage and lateral-flux conservation laws.

## Fixed cohort and definitions

The analysis reuses the established common cohort without relabeling:

- 337 first majority-to-minority crossings from 37 original seeds;
- at least five crossing-free time units before every event;
- 185 continued events, defined by at least four **additional** net fish
  transferred during the five units after first crossing;
- 152 failed events with fewer than four additional fish.

All organization measurements are made at negative times, from 5 to 0.1 time
units before the first fish crosses.  Future population transfer is used only
as the commitment label.  Complete opposite-chamber relocation is not the
outcome.

For the original source school, with `s=+1` for a right source and `s=-1` for
a left source, the three prespecified diagnostics are

```text
P = |mean_i h_i|,

M = mean_i [((r_i-r_COM) x h_i)_z / |r_i-r_COM|],

A = -s mean_i cos(theta_i).
```

`P` is source polarization magnitude, `|M|` is absolute normalized angular
momentum about the source-school centroid, and `A` is fixed horizontal
alignment toward the opening.  The attached starting script used a
radius-weighted angular momentum and a new three-fish outcome threshold; both
were corrected to match the established normalization and event definition.

Reflection leaves `P`, `|M|`, and `A` unchanged and reverses signed `M`; the
synthetic reflection test passes exactly.

## Commitment: continued versus failed crossings

Comparisons are stratified by the established initial `|q|` bins so different
population imbalance distributions cannot create the effect.  Uncertainty is
bootstrapped over original seeds.

At lead 2, continued episodes are not more organized toward commitment:
gap alignment is actually lower by `-0.0584`, 95% seed interval
`[-0.0981,-0.0196]`, and the milling difference is unresolved.  The pathways
then separate before any crossing:

| lead | success-failure `P` | success-failure `|M|` | success-failure `A` |
|---:|---:|---:|---:|
| 1.0 | -0.0138 `[-0.0166,-0.0108]` | +0.0401 `[0.0187,0.0587]` | +0.0504 `[0.0120,0.0836]` |
| 0.5 | -0.0138 `[-0.0166,-0.0112]` | +0.0560 `[0.0360,0.0744]` | +0.1031 `[0.0700,0.1346]` |
| 0.1 | -0.0187 `[-0.0217,-0.0159]` | +0.0761 `[0.0566,0.0939]` | +0.1250 `[0.0948,0.1540]` |

The total polarization remains very high: at lead 0.1 it is 0.961 for
continued events and 0.976 for failures.  The negative contrast is therefore
not disintegration.  It indicates that a small portion of almost perfectly
parallel motion has been converted into turning/curvature and axial motion
toward the opening.

The same chronology appears in changes measured from lead 2.  By lead 0.1,
continued-minus-failed change is

```text
Delta P       = -0.0144  [-0.0184,-0.0103]
Delta |M|     = +0.0631  [ 0.0450, 0.0803]
Delta A       = +0.1834  [ 0.1482, 0.2156].
```

Because both endpoints precede the first crossing, these changes cannot be a
consequence of chamber membership transfer.

## Independent-seed commitment prediction

Seeds 0--29 and 30--39 were kept separate.  No combined classifier was fit;
each number is a prespecified one-coordinate, q-stratified rank diagnostic.

At lead 0.5, the relationship is reproduced in both cohorts.  On the ten
held-out seeds (129 events, 75 continued), the q-stratified AUCs are

| diagnostic | held-out AUC | 95% seed interval |
|---|---:|---:|
| `P` | 0.164 | `[0.121,0.224]` |
| `|M|` | 0.748 | `[0.678,0.802]` |
| `A` | 0.781 | `[0.673,0.862]` |

The corresponding training AUCs are 0.187, 0.676, and 0.622.  Thus lower
polarization, higher milling magnitude, and higher gap alignment all replicate
in independent trajectories.  The direction is already fixed by the physical
definitions; no threshold or nonlinear combination was selected on test
seeds.

At lead 1, the earlier fixed-axis study found no useful trained fixed-axis
predictor.  The present result is consistent with it: the training
q-stratified AUC for `A` is 0.502 at lead 1 and becomes 0.622 only at lead 0.5.
The new evidence is a **late** commitment gate that the earlier lead-1
analysis did not sample.

## Approach initiation is not commitment

The existing exact-q matched controls distinguish a future crossing approach
from a configuration with no population flux for two units before and five
units afterward.  There are 329 pairs at each of leads 2 and 1; the held-out
subset contains 128 pairs from ten seeds.

Held-out event-minus-control differences are

| lead | `P` | `|M|` | `A` |
|---:|---:|---:|---:|
| 2.0 | +0.1036 `[0.0559,0.1508]` | -0.0090 `[-0.0317,0.0172]` | +0.2245 `[0.1334,0.3019]` |
| 1.0 | +0.0991 `[0.0692,0.1329]` | +0.0085 `[-0.0083,0.0248]` | +0.2842 `[0.1879,0.3650]` |

Polarization and gap alignment therefore characterize approach initiation.
Milling does not.  Within the crossing cohort, however, polarization has the
opposite relation to subsequent continuation, while milling and additional
alignment emerge only late.  Predicting an approach and predicting commitment
are physically and statistically distinct problems.

## Persistent circulation test

Signed milling is strongly persistent from lead 2 to lead 0.5 in both groups,
but more so in failed attempts:

```text
failed:    Spearman 0.804 [0.737,0.848], sign retention 0.846
continued: Spearman 0.650 [0.536,0.747], sign retention 0.813
difference in Spearman: -0.153 [-0.287,-0.024].
```

Thus persistent clockwise or counterclockwise circulation does not commit the
school.  Failures preserve their circulation phase more faithfully.  The late
increase in `|M|` during successful episodes accompanies reorientation and
partial loss of milling-phase persistence, consistent with a curved turn or
wall-guided funneling rather than a stable milling attractor.

## Implication for a minimal transition-rate theory

The results support a renewal decomposition rather than a new potential:

```text
rate of substantial transfer
    = rate of source-to-gap approach attempts
      x probability that an attempt commits.
```

Polarization and fixed gap alignment diagnose the attempt.  During the final
0.5 unit, additional gap alignment together with transient curved organization
diagnoses whether the attempt is likely to continue.  The exact storage
balance then governs how that organization becomes longitudinal supply,
lateral retention, and neck flux.

This does **not** establish a closed rate theory yet.  A renewal model would
still need a prospective definition of attempt start, a reset/retreat
criterion, inter-attempt waiting-time statistics, and validation that the
post-retreat source state loses sufficient memory.  The present analysis does
not promote `|M|` to an additional slow coordinate, fit a potential, or add a
combined empirical predictor.

## Scope

The organizational changes are held-out predictive associations and precede
crossing; they are not controlled causal interventions.  They may be generated
by wall-guided turning, collective geometry, or their coupling.  The analysis
does rule out two simpler mechanisms:

- commitment is not caused by a loss of collective order;
- commitment is not a persistent-milling-state transition.

No microscopic model was changed, no new trajectory was generated, and no
complete-relocation rate is inferred.
