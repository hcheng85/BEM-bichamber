# Ballistic encounter geometry versus collective reorientation

## Decision

A frozen coherent-school encounter is **not** a viable predictor of either initial entry or sustained population transfer from the beginning of a physical approach attempt. It is highly specific when it predicts a ray through the opening, but it misses nearly all entries because the source-school mean heading is usually still outward or oblique when the leading part first occupies the approach layer. Actual fish positions and finite school extent are already retained, and solid-wall obstruction eliminates only one otherwise eligible five-unit ray in 1,811 attempts. The failed assumption is therefore frozen collective orientation, not point-school geometry or an unmodeled static wall shadow.

Collective reorientation is needed even to bring the coherent school onto an aperture-intersecting path. Once entry has begun, evolving nonreciprocal visual alignment is causally important for the **amount** transferred: in the matched-noise intervention the initial geometric prediction is exactly identical in both arms, but suppressing internal visual nonreciprocity reduces mean five-unit transfer by 4.44 fish. The intervention does not resolve a change in first-entry probability or the binary four-fish threshold.

The recurrent approximately 15.5-unit sweep remains a separate, empirically measured encounter clock. This analysis neither predicts that period nor revises the failed Poincare map.

## Parameter-free coherent ballistic prediction

At a prospective attempt state, let `S` be the fish in the source chamber and let the production swimming speed be `v0=1`. The frozen coherent velocity is

\[
 V_0={v_0\over N_s}\sum_{i\in S}(\cos\theta_i,\sin\theta_i)
     =v_0 Pm .
\]

This is the instantaneous mean self-propulsion velocity, including the observed polarization magnitude `P`; unit COM speed is not imposed. Each actual source-fish position is translated rigidly,

\[
 r_i(t)=r_i(0)+V_0t.
\]

For a right or left source with sign `s`, the inward direction is `D=-s`. If `D V_{0x}>0`, fish `i` reaches the divider after

\[
 t_i={s x_i\over D V_{0x}},\qquad
 y_i^{\rm hit}=y_i+V_{0y}t_i.
\]

It is predicted to enter only if `0<t_i<=T`, `|y_i_hit|<=gap/2`, and the segment is unobstructed by the production wall segments and rounded tip circles. This uses the actual gap `2.8635642126552705`. The opening, horizons, and wall test contain no fitted distance, angle, or probability threshold.

The prospective surface is the established threshold-free storage-excursion onset: the first 0.5-time sample at which either approach layer becomes occupied after an empty sample. The initiating side is chosen from current storage, without future crossing information. The primary analysis retains the 1,811 attempts initiated from the current majority chamber. Seeds 0--29 are training and seeds 30--39 are held out. Outcomes are fixed independently of the predictor:

- entry within two time units;
- signed population change at `+2`;
- signed population change and at least five net fish at `+5`.

## Prospective held-out result

There are 425 majority-source attempts in the ten held-out seeds. Fifty enter within two units and 69 transfer at least five net fish by five units.

| target | coherent ballistic AUC (95% seed interval) | stored-current `Ja` AUC | difference |
|---|---:|---:|---:|
| entry by `+2` | 0.545 `[0.511,0.587]` | 0.770 `[0.706,0.827]` | -0.225 `[-0.287,-0.148]` |
| positive net transfer at `+2` | 0.555 `[0.526,0.579]` | 0.715 `[0.646,0.779]` | -0.160 `[-0.226,-0.084]` |
| at least five net fish at `+5` | 0.522 `[0.490,0.556]` | 0.644 `[0.571,0.710]` | -0.122 `[-0.186,-0.062]` |

The parameter-free binary ray prediction has 10.0% sensitivity and 98.9% specificity for two-unit entry. It predicts entry in 2.12% of attempts, compared with 11.76% observed. For five-unit sustained transfer its sensitivity is only 1.45%; it predicts the threshold in 0.47% of attempts, compared with 16.24% observed.

Ballistic count retains weak rank information about transferred amount: held-out Spearman correlations are 0.238 `[0.139,0.325]` at two units and 0.181 `[0.073,0.262]` at five units. It underpredicts by 0.285 fish `[-0.493,-0.111]` and 2.59 fish `[-3.59,-1.68]`, respectively. The training seeds show the same qualitative failure, with entry and sustained-transfer AUCs of 0.521 and 0.514. This is therefore not a held-out distribution accident.

These results do not contradict the earlier individual-heading aperture current. An instantaneous individual ray can be a useful final geometric gate close to entry. The present, stricter question starts at physical approach onset and asks whether a coherent school with its **current mean heading** already determines subsequent transfer. It does not.

## Exact transport budget and observational turning association

The fixed 337-event first-crossing cohort provides a second test at exactly one unit before entry. A frozen coherent ray through the defining entry time does not distinguish later continuation in held-out seeds: its AUC is 0.483 `[0.444,0.526]`. Allowing the same frozen ray two units, through one unit after entry, raises AUC only to 0.587 `[0.534,0.634]`, compared with 0.695 `[0.582,0.811]` for the actual lead-1 approach current.

For each event, the exact established budget is integrated from `t=-1` to `+1`:

\[
 \Delta n_{\rm neck}
 =\int {N J_a\over\ell}\,dt
  +\int(-A_W)\,dt
  +\int(-L_W)\,dt
  +\int\epsilon_{\rm path}\,dt .
\]

The maximum integrated closure residual is `1.04e-14` fish. Defining the ballistic residual as actual net neck transfer minus the number of frozen coherent rays predicted to arrive by `+1`, held-out seed-bootstrap associations are

| exact/physical contribution | Spearman correlation with residual | 95% interval |
|---|---:|---:|
| longitudinal supply | 0.731 | `[0.617,0.802]` |
| lateral contribution | 0.390 | `[0.177,0.578]` |
| storage contribution | -0.327 | `[-0.513,-0.102]` |
| lead-1 nonreciprocal alignment turn | 0.429 | `[0.256,0.549]` |

The conservation identity is exact; these correlations are not. They show that transfer missed by the frozen ray is realized primarily by evolving directed supply, with lateral exchange and storage determining its conversion to neck flux. The association with the previously derived internal visual-asymmetry term is reproduced in both training and held-out seeds, but observational association alone would not establish causality.

## Matched-noise causal comparison

The existing intervention begins from 14 archived lead-1 microstates with eight paired noise streams per state. The baseline and symmetrized arms have identical positions and headings at restart, so their ballistic predictions are exactly identical for all 112 clone pairs; the maximum geometry difference is zero.

Symmetric minus baseline effects are

| outcome | paired effect | 95% seed interval |
|---|---:|---:|
| first-unit turning toward the gap | -0.0285 | `[-0.0465,-0.0105]` |
| five-unit net transfer | -4.44 fish | `[-7.12,-1.85]` |
| first-entry probability | +0.036 | `[-0.009,+0.098]` |
| four-fish success probability | -0.009 | `[-0.161,+0.161]` |

Thus identical initial geometry can evolve into substantially different transfer amounts when the directed internal alignment rule is changed. Nonreciprocal collective turning is causally important for sustained supply and transfer amount. The already established, crossing-conditioned intervention does not demonstrate that it is necessary for the first fish to enter, nor does it resolve the coarse binary commitment threshold.

## Can a derived turning correction rescue the encounter theory?

The minimal parameter-free kinematic correction is clear. If a prescribed collective angular rate `Omega0` is frozen together with speed `V=v0 P`, then every actual fish position follows a common circular arc,

\[
 x_i(t)=x_i+{V\over\Omega_0}
 [\sin(\Theta_0+\Omega_0t)-\sin\Theta_0],
\]

\[
 y_i(t)=y_i-{V\over\Omega_0}
 [\cos(\Theta_0+\Omega_0t)-\cos\Theta_0],
\]

with the straight ray recovered as `Omega0 -> 0`. The production dynamics give an exact instantaneous collective rate by projecting the wall, social, clipping, and Ito terms onto the mean-heading angle. No fitted turning coefficient is required.

This is not yet a predictive low-dimensional correction. Previous source-exact tests found that a frozen microscopic rate is useful over roughly 0.1 time units, whereas entry and transfer here require 1--5 units. The directed visual graph and relative headings evolve on the approach time scale and cannot be adiabatically eliminated. Supplying their instantaneous microscopic torque at every step would reconstruct the original many-fish dynamics rather than close an analytical encounter model.

The supported minimal mechanism is therefore

```text
measured recurrent chamber sweep
  -> repeated finite-school encounters with the approach layer
  -> evolving wall and nonreciprocal social reorientation
  -> maintained longitudinal supply
  -> exact storage/lateral/neck-flux conversion.
```

Geometry alone does not predict entry from attempt onset, and a frozen turning correction is not controlled over the necessary flight time. A minimal analytical transfer-per-approach theory becomes viable only if the collective angular response can be derived over the full encounter without introducing an empirical state hierarchy. The present result stops at that precise missing closure.

## Reproduction and scope

```bash
uv run --with numpy --with scipy --with matplotlib \
  python results/nohydro_gap20_descriptive_transition/analyze_ballistic_encounter_transfer.py
```

The analysis reads the existing 40 trajectories and compact prior outputs. It runs no new microscopic simulations, changes no production code, fits no potential or predictor, and does not modify the Poincare section.
