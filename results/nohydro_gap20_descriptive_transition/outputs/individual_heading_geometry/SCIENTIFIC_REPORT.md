# Individual heading geometry and the necessity of encounter-time reorientation

## Decision

Initial angular dispersion explains a meaningful minority of first entries, but not most of them. At the established prospective approach onset, independent fish rays predict entry substantially better than a rigid mean-heading school. Nevertheless, they miss about half of all held-out entries, and their finite-two-unit version misses about 69%. Existing positions and headings therefore do not determine entry: dynamical reorientation during the encounter remains essential.

For the successful fixed-cohort entries missed by both frozen models, the median initial coherent turn required to expose any source fish to the aperture is 0.466 radians. The parameter-free independent-noise RMS accumulated before entry is only 0.107 radians, a median requirement-to-noise ratio of 4.37. After correcting the transformation from the local inward-facing angle to the global heading frame, the actual school mean heading turns in the required direction in 97.1% of held-out cases and reaches the full initial geometric angle in 47.8%. Fish positions, the aperture bearing, and individual headings evolve simultaneously, so full initial-angle attainment is not required for every entry.

The corrected source-exact torque projection strongly favors social rather than wall-guided turning during the final approach. At lead 1, the held-out social angular rate toward the required global turn is +0.1360 `[0.1225,0.1509]`, while the wall rate is -0.0160 `[-0.0297,-0.0062]`. The known internal nonreciprocal alignment term remains positive and resolved. This gives a defensible microscopic route to the required turn, but not a closed coarse theory: its magnitude still depends on the evolving directed visual graph and relative headings.

## Two parameter-free frozen predictions

The prospective surface is unchanged. An attempt starts at the first 0.5-time sample at which the established weighted approach storage becomes positive after an empty sample. The initiating side is selected from current storage without using a future crossing or outcome. The primary cohort contains 1,811 complete attempts initiated from the current majority chamber; seeds 0--29 are training and 30--39 are held out.

For fish `i` at actual position `(x_i,y_i)`, define local inward direction `D=-s`, where `s` is the source-chamber sign. The opening subtends the exact angular interval

\[
 \phi_i^- = \operatorname{atan2}(-h-y_i,\,s x_i),\qquad
 \phi_i^+ = \operatorname{atan2}( h-y_i,\,s x_i),
\]

with `h=gap/2` and `gap=2.8635642126552705`. Angles are measured from the inward horizontal direction. A ray targets the opening when its local heading lies in this interval and its segment is unobstructed by the production straight walls and rounded tip circles.

Two predictions use the production speed `v0=1`:

1. **Coherent ray:** every source fish is translated with the source mean-heading vector `P m`.
2. **Individual rays:** fish `i` is translated with its own instantaneous velocity direction `(cos theta_i,sin theta_i)`.

Both retain the full finite school geometry. Counts are computed both without a time cutoff, for eventual frozen targeting during the complete excursion, and with a fixed two-unit cutoff. Neither uses the eventual reset time to define the ray.

The minimum coherent turn is the shortest signed angular distance from the current mean heading to the union of all fish-specific aperture intervals:

\[
 \delta\Theta_{\min}=
 \mathop{\arg\min}_{i,\,\phi\in[\phi_i^-,\phi_i^+]}
 |\operatorname{wrap}(\phi-\Theta)|.
\]

This is a geometric requirement, not a fitted turning threshold. The aperture calculation is performed in the local inward-facing frame. For source sign `s`, its signed global-heading increment is exactly `-s deltaTheta_local`. The implementation now applies this transformation before comparing with global heading changes or torque rates. Reflection tests require the global increment to reverse sign when left and right chambers are exchanged. A separate seam test chooses the nearest wrapped aperture endpoint at `+/-pi`; both tests pass.

## Independent-seed entry validation

The held-out cohort contains 425 complete attempts, of which 159 eventually enter in their initiating direction.

| prospective predictor | entry AUC, 95% seed interval | sensitivity | specificity |
|---|---:|---:|---:|
| coherent frozen target | 0.527 `[0.507,0.543]` | 6.9% | 98.5% |
| individual frozen targets | 0.666 `[0.620,0.710]` | 50.9% | 80.8% |
| coherent rays within 2 units | 0.508 `[0.498,0.519]` | 3.1% | 98.5% |
| individual rays within 2 units | 0.636 `[0.599,0.672]` | 30.8% | 96.2% |

Thus individual angular dispersion is real physical information that the coherent ray discards. It explains roughly half of eventual entries, not most of them. The same ordering appears in training seeds; it is not produced only by the held-out ensemble.

Among held-out successful attempts, the median minimum coherent turn is 0.410 radians `[0.327,0.470]`, while the median fish-specific aperture width is 0.579 radians `[0.551,0.621]`. A wide angular opening as seen by nearby leading fish therefore coexists with a mean school heading that is usually outside every frozen aperture interval.

## Complete-excursion transfer, not a fixed future window

The primary transport reward is exactly the direction-aligned chamber-population change between the empty pre-attempt state and the physical storage reset. It obeys the existing crossing conservation identity and cannot include the next approach excursion.

On held-out seeds, the coherent target count has unresolved association with excursion reward: Spearman `rho=0.091` `[-0.049,0.192]`. Individual target count carries modest information, `rho=0.307` `[0.183,0.416]`, and its AUC for at least five net fish during the excursion is 0.603 `[0.519,0.689]`. The coherent AUC is 0.522 `[0.490,0.557]`. Initial individual geometry therefore influences how much population is available to enter, but it is far from a deterministic transfer law.

The requested fixed-five-unit audit confirms why excursion reward must be primary:

- 3.06% of held-out attempts reset and begin another approach before `onset+5`, so their fixed-five-unit population change can combine attempts.
- Fixed-five-unit and excursion-end fish counts differ in 17.65% of held-out attempts.
- The five-fish label changes in only 0.47%, but this numerical robustness does not make the underlying reward an attempt-local quantity.

All reported transport associations above therefore use excursion-end reward. Short transfer, sustained excursion transfer, and whole-school relocation remain distinct. The previously observed nine opposite-dominance relocations are not treated as a measured rate or as additional evidence for this entry mechanism.

## Required turn versus collective angular diffusion

For source mean-heading angle `Theta` and polarization `P`, the exact independent-noise quadratic variation from the production law `d theta_i = ... + I_n dW_i`, `I_n=0.5`, is

\[
 \sigma_\Theta^2 = {I_n^2\over N_s^2P^2}
 \sum_{i\in S}\cos^2(\theta_i-\Theta)
 \quad\text{per unit time}.
\]

No effective diffusivity is fitted. For held-out successful attempts, the median instantaneous scale is 0.0628 rad/sqrt(time) `[0.0607,0.0635]`.

Of the 315 fixed-cohort events uniquely associated with a prospective storage excursion, 169 are missed by both unlimited frozen-ray models. The held-out missed subset contains 69 events from all ten seeds. Its medians are

| quantity | estimate | 95% seed interval |
|---|---:|---:|
| initial minimum turn | 0.466 rad | `[0.434,0.493]` |
| elapsed onset-to-entry time | 3.064 | descriptive median |
| independent collective-noise RMS | 0.107 rad | `[0.100,0.115]` |
| required turn / noise RMS | 4.37 | `[4.04,4.74]` |
| actual mean-heading turn toward initial aperture | 0.483 rad | `[0.413,0.552]` |

Independent angular noise is therefore too small to be the typical source of the required coherent deflection. This comparison is a scale test, not a claim that entry is a Gaussian first-passage event: social coupling, wall torques, and the moving aperture bearing invalidate such a scalar Brownian model.

The actual mean heading turns in the initially required global direction in 97.1% `[94.4%,100%]` of held-out missed entries. It reaches the full initial minimum angle in 47.8% `[36.4%,59.1%]`. The latter is not a contradiction. As the spatially extended school translates and deforms, a different fish can become the leader, individual headings change, and the fish-specific aperture intervals rotate. Collective mean reorientation is nearly universal in direction but is not the entire kinematic path to entry.

## Wall-guided versus social turning

For the missed entries, the previously audited exact collective-heading projection was evaluated at fixed leads before crossing. Positive rates below point in the sign of the onset geometric requirement.

| held-out rate | lead 1 | lead 0.5 | lead 0.2 |
|---|---:|---:|---:|
| wall torque | -0.0160 `[-0.0297,-0.0062]` | -0.0321 `[-0.0541,-0.0164]` | -0.0322 `[-0.0583,-0.0137]` |
| total social torque | +0.1360 `[0.1225,0.1509]` | +0.1042 `[0.0912,0.1231]` | +0.0978 `[0.0883,0.1112]` |
| deterministic total | +0.1200 `[0.1072,0.1312]` | +0.0721 `[0.0594,0.0831]` | +0.0656 `[0.0401,0.0856]` |

The wall contribution is resolved opposite to the required global turn at every lead. It may shape the recurrent chamber trajectory and bring the school into the encounter region, but it does not supply the final deflection toward the aperture in this subset. The social term is resolved toward the requirement at all three held-out leads. This corrected sign separation is stronger than the originally committed result.

The internal visual-nonreciprocity contribution to improving gap alignment is independently positive in both splits. In the held-out missed subset it is +0.0735 `[0.0628,0.0862]` at lead 1, +0.0563 `[0.0442,0.0741]` at lead 0.5, and +0.0529 `[0.0392,0.0686]` at lead 0.2. Combined with the prior matched-noise intervention, this supports the physical route

```text
individual heading dispersion gives some immediate aperture rays
  + directed internal visual alignment turns the remaining coherent school
  -> maintained approach supply
  -> exact storage/lateral/neck-flux conversion.
```

The torque association here is retrospective; causality comes from symmetrization experiments. The new outcome-blind prospective restart test leaves first-entry probability unchanged (paired effect 0.000 with production noise, 95% seed interval `[-0.098,0.107]`) while reducing complete-excursion transfer by 3.125 fish `[0.205,8.215]`. Nonreciprocity is therefore a defensible source of reorientation and sustained propagation, but is not proven necessary for every first fish.

## Consequence for a minimal analytical encounter theory

A position-only or coherent frozen-ray theory is rejected. An independent-ray theory is a useful parameter-free geometric gate and explains part of initiation, but it does not close entry or transfer. Adding independent collective diffusion also fails by scale: the typical missed entry requires more than four noise RMS.

The minimal viable analytical structure would combine the exact fish-specific angular aperture with a derived collective angular response from the production social rule. The known nonreciprocal mechanism has the correct sign, magnitude scale, causal transfer effect, and persistence near entry. What remains unavailable is its value from coarse school geometry alone: it depends on directed visible-edge weights and organized relative-heading differences. Supplying those microscopic quantities dynamically would reproduce the many-fish model rather than produce low-dimensional closure.

The present stopping conclusion is therefore precise: geometry plus initial angular dispersion explains a substantial minority of archived entries; encounter-time dynamics is essential for the remainder; and the prospective intervention does not isolate internal nonreciprocity as necessary for binary first entry. It does establish that nonreciprocity supports subsequent population transfer. A low-dimensional analytical encounter theory is not closed until the collective turning response can be derived from coarse physical information without adding a phenomenological hierarchy.

## Reproduction

```bash
uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/analyze_individual_heading_geometry.py
```

No microscopic simulation was run, no production source or event definition was changed, and no fitted map, potential, or additional collective coordinate was introduced.
