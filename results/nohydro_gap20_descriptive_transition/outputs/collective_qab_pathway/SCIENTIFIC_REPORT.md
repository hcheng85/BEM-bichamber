# Collective q-a-b pathway: partial mechanism, inadequate closure

## Concrete conclusion

The data support a reproducible **persistent orientation → approach-layer loading → population release** sequence, but not a newly generated reorientation immediately before crossing. Weighted approach storage `b` contains information beyond coarse `(q,a)` for predicting whether at least four fish subsequently transfer. It does not give a resolved improvement for the magnitude of net transfer, and the tested `(q,a,b)` representation retains substantial dependence on earlier states. Therefore `(q,a,b)` is a useful physical description of the transition pathway, but it is not an adequate closed dynamical state for a transition-rate theory.

No simulator was run or modified. All results use the existing 40 no-hydrodynamics gap-20 trajectories and the committed event/control definitions.

## Cohorts and definitions

- Event chronology: 337 majority-to-minority first crossings from 37 original seeds, each with at least five crossing-free time units before entry; 185 events transferred at least four **additional** net fish during the next five units and 152 did not.
- Matched comparison: 658 exact-`q` event/control measurements at leads one and two, using the previously selected stable no-crossing controls.
- Prospective diagnostic: 6,308 fixed-grid observations. Seeds 0–29 trained the empirical conditional means and seeds 30–39 supplied 1,577 independent held-out observations, including 94 four-fish outcomes.
- `b=W/(N ell)` is the weighted source-side approach storage from the earlier exact conservation analysis. The present `j_heading` is only a heading-based kinematic indicator and is not substituted for the exact stored-velocity current.

The prior conservation law is an exact balance involving approach current, storage change, lateral exchange, neck flux, and the audited discretization correction. The pathway and prediction results below are statistical observations; they are not additional exact identities.

## Event-aligned pathway

Both outcome groups are already strongly gap-oriented five units before crossing: seed-equal mean `a=0.646` for continued events and `0.620` for failed events, while `b` is essentially zero. The distinguishing behavior is persistence, not a late upward reorientation. Because `a` projects mean heading onto the instantaneous COM-to-gap direction, its evolution can combine fish-heading changes with rotation of the school-geometry reference; this analysis does not attribute it to a particular turning torque.

- At `t=-3`, continued events retain `a=0.657` with `b=0.0047`; failed events have already fallen to `a=0.561` while loading earlier (`b=0.0285`).
- At `t=-1`, continued events still have `a=0.559`, compared with `0.291` for failed events. Storage is then similar (`b=0.151` versus `0.146`).
- Failed-event storage peaks earlier and lower: `b=0.164` at `t=-0.5`. Continued-event storage peaks at `b=0.228` at first entry.
- Population imbalance relative to the original source is unchanged before entry by construction. After entry it falls from `0.575` at `t=-5` to `0.397` at `t=+1` and `0.243` at `t=+5` for continued events. Failed events remain near `0.467` at `t=+5`.

The same qualitative separation appears within all three initial-imbalance strata. Thus the aggregate sequence is not produced only by different initial loading. The physical picture is that sustained source-school orientation supplies and maintains a larger near-gap reservoir; release of that reservoir follows first entry. Chronology alone does not prove that orientation causes loading or that storage causes transfer.

## Exact-q matched controls

Event configurations carry more weighted approach storage than their exact-`q` matched stable controls:

- Lead 2: `b_event-b_control = 0.0561`, seed-bootstrap 95% CI `[0.0389, 0.0747]`.
- Lead 1: `b_event-b_control = 0.1488`, CI `[0.1280, 0.1741]`.

Gap-directed polarization is also larger at both leads (`delta a=0.523` at lead 2 and `0.434` at lead 1). These retrospective case-control contrasts show that entry is approached through a jointly oriented and spatially loaded state. Because controls were selected for stable future no-flux behavior, these differences are not unbiased entry-hazard estimates.

## Held-out information contributed by b

All 1,577 held-out samples fell in adequately supported training cells. Adding a training-only median split in `b` to the predeclared `(q,a)` cells gave:

- Four-fish Brier score: `0.05290` for `(q,a)` and `0.05110` for `(q,a,b)`.
- Seed-mean Brier improvement: `0.00180`, CI `[0.00060, 0.00285]`; eight of ten test seeds improved. This is a small but resolved improvement, about 3.4% of the baseline Brier score.
- Future-net-fish MSE: `21.059` for `(q,a)` and `20.506` for `(q,a,b)`.
- Seed-mean MSE improvement: `0.547`, CI `[-0.087, 1.405]`; five seeds improved and five worsened. The approximately 2.6% aggregate change is statistically unresolved.

Therefore `b` adds reproducible information about entering the larger-transfer class, but does not robustly predict how many fish will transfer. It is more than a redundant accounting variable, yet it is not a sufficient slow coordinate.

## Closure and missing physics

The coarse `(q,a,b)` diagnostic fails the history test. On held-out seeds, future-net residuals remain associated with the state five units earlier:

- lagged `b`: Spearman `0.215`, seed-bootstrap CI `[0.180, 0.248]`;
- lagged `a`: Spearman `0.182`, CI `[0.120, 0.246]`.

This rejects approximate closure of the tested empirical representation. It does not prove that every continuous `(q,a,b)` description is non-Markovian, because the model uses coarse `q,a` bins and a binary storage split. Indeed, residual transport is strongly associated with continuous present `a` (Spearman `-0.369`, CI `[-0.415, -0.318]`) and weakly with continuous present `b` (`0.054`, CI `[0.018, 0.099]`), demonstrating material within-cell resolution loss.

The focused geometrical audit does not point to the global signed circulation component `tangent_p`: its residual association is unresolved (`-0.010`, CI `[-0.044, 0.025]`). In contrast, the localized approach heading current `j_heading` remains associated with residual transport (`0.151`, CI `[0.072, 0.223]`). This favors missing localized flux information over global circulation phase as the immediate explanation, while remaining descriptive because `j_heading` omits collision overrides.

## Physical mechanism and next theoretical step

The supported mechanism is conditional and mesoscopic: a school with persistent gap-directed organization loads a weighted approach layer and can then release multiple fish after first entry. Failed events lose that gap-directed organization earlier and their smaller reservoir drains without comparable population transfer. This does not distinguish fish-heading rotation from rotation or deformation of the COM-to-gap geometry. It describes short continuation episodes, not complete opposite-chamber relocation and not a whole-school committor.

The first-principles next step is not an effective potential in `(q,a,b)`. It is to augment the exact storage balance with the **source-exact stored-velocity approach current and weighted lateral boundary flux**, including collision/path corrections already defined by the conservation analysis, and test whether those localized fluxes remove the held-out history dependence. That directly targets the missing transport physics revealed here. A circulation coordinate should be reconsidered only if it improves this source-grounded flux closure; the present signed global circulation diagnostic does not support it as the leading omitted variable.

## Reproducibility and limits

Uncertainty is clustered by original simulation seed. Event-aligned success/failure comparisons are retrospective, whereas the held-out diagnostic is prospective on independent seeds. `focused_physical_diagnostics.json` records the seed-bootstrap chronology and the targeted residual-coordinate audit. No raw MAT file, simulator source, fitted potential, phenomenological drift, or large cache is included.
