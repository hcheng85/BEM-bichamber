# Fixed-axis versus COM-to-gap polarization

Both coordinates use precisely the same raw snapshots, eligibility rules, and future outcomes from the existing 40 no-hydrodynamics gap-20 trajectories. No event was relabeled and no simulation was run.

## Validation and cohorts

- First-entry risk origins: 5792; held-out: 1370 with 129 positives.
- Event-aligned continuation events: 388; held-out: 144 with 75 positives.
- Unconditional post-transient snapshots: 32040 (801 per seed at spacing 1).
- Raw-state q/a_COM cross-checks: 6180; maximum a_COM error 5.551e-16.
- Archived COM AUC and Brier values reproduce exactly for both outcomes.

## Same 37 matched pre-crossing trajectories

The moving COM coordinate decreases before events: event Δa_COM = -0.336685, 95% CI [-0.4266716773274447, -0.254812055394117]. The fixed horizontal coordinate instead increases: event Δa_fixed = 0.165288, 95% CI [0.1197679765519996, 0.21595198801630627].

Relative to unchanged controls, COM Δ = -0.312449, CI [-0.44836873272822364, -0.1908482734215461]; fixed-axis Δ = 0.090205, CI [-0.0016778786497714127, 0.19164534332106892]. The latter interval narrowly includes zero.

Therefore the apparent late decline of COM-to-gap polarization is largely reference-direction rotation, not fish turning away from the horizontal opening.

## Held-out prediction: seeds 0–29 train, 30–39 test

| Outcome | q AUC | q+a_COM AUC | q+a_fixed AUC | fixed-minus-COM seed-mean AUC difference (95% CI) | fixed-minus-COM Brier gain (95% CI) |
|---|---:|---:|---:|---:|---:|
| First entry | 0.619484 | 0.785410 | 0.699380 | -0.108936 [-0.16330796057493552, -0.06611352215561188] | -0.005797 [-0.00880944249960204, -0.003145746830780521] |
| Continuation | 0.564155 | 0.813720 | 0.502126 | -0.264373 [-0.3766073232323232, -0.1482823322510822] | -0.079779 [-0.10457051604018142, -0.05230598137881583] |

Positive Brier gain favors fixed-axis. Both fixed-minus-COM differences are negative and exclude zero: COM-to-gap is materially better, especially for continuation. Fixed-axis modestly improves first-entry prediction over q alone, but is worse than q alone for continuation.

## Unconditional empirical distributions

The equal-seed empirical p(q,a) landscapes use 32040 finite-time snapshots. q=0 occurs in 263 snapshots and has no majority coordinate. Left/right chambers are empty in 7483/8071 snapshots; these polarizations remain undefined rather than zero. Both chambers are occupied in 16486 snapshots (0.515).

The shifted -ln p panels are statistical occupation landscapes only. Fixed-axis retains horizontal heading; COM-to-gap also resolves off-axis location relative to the opening.

## Separate left/right and minority polarization

For initiation, the minority diagnostic retains 2593/5792 observations with at least five minority fish; held-out joint-cell coverage is 0.954. Minority Brier gain = 0.000646, CI [-0.0022877372808423504, 0.003516172290773991]; AUC gain = 0.016647, CI [-0.018063861327635718, 0.055112424659797814].

For continuation it retains 355/388 observations; coverage is 0.858. Minority Brier gain = -0.001407, CI [-0.015486737114493089, 0.013208771839434344]; AUC gain = 0.055787, CI [-0.0006145014964459366, 0.1206869262424818].

Neither outcome shows a statistically resolved benefit from minority polarization. Do not add it to the primary state; it is undefined for localized schools and produces sparse cells.

## Recommendation

**Retain majority-school COM-to-gap polarization a_COM for the eventual empirical collective landscape.** It is physically local to the opening and decisively outperforms fixed-horizontal polarization on independent seeds. Keep a_fixed as a mechanism diagnostic: headings turn horizontally toward the opening even while the moving COM projection falls.

## Scope and interpretation

These are finite-time empirical occupation distributions and held-out associations, not fitted potentials, rate barriers, Markov closure, or a whole-school relocation model. The rare complete relocations are not the outcomes tested here.
