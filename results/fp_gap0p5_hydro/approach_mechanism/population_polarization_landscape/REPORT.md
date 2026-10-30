# Population imbalance and polarization as a collective switching landscape

## Decision

The proposed pair `(m, P)` is a sensible **descriptive order-parameter plane**
for collective chamber switching, but is **not a predictive two-coordinate
energy landscape** for the present hydrodynamic, active bichamber model.

- `m = (N_R - N_L)/N` represents population occupancy. The two metastable
  chamber-dominated states lie near `m = +/-1`; switching proceeds through
  intermediate occupancy. This does **not** imply a universal inverse law
  between imbalance and transition rate. The reduction `|m|` occurs *during*
  a switch; it does not specify the hazard of initiating one while `|m| = 1`.
- `P = |N^-1 sum_i (cos(theta_i), sin(theta_i))|` is **global** polarization.
  It need not distinguish schools pointing toward versus away from the gap.
  Nor does a high global value guarantee two separately coherent schools.

For an exhaustive two-chamber assignment, conservation of fish gives the exact
event-counting identity

```
m(t) - m(0) = (2/N) [C_{L->R}(t) - C_{R->L}(t)].
```

Thus imbalance is an occupancy coordinate whose **increment** is governed
by *signed net crossing current*. A crossing rate is related to `dm/dt`
(interpreted as a jump current), not universally to `1/|m|`. If fish can
occupy a separate aperture/corridor class, the identity needs corresponding
flux terms or a consistent chamber-assignment convention.
- The existing episodes establish that signed collective heading at the first
  crossing is substantially more discriminating for successful recruitment
  than `P`. This is relevant to commitment *after the first crossing*, not a
  standalone prospective predictor of the first crossing.

The useful outcome is a sharper *two-stage* modeling target: an incoming,
directed leading-edge configuration determines nucleation; after a fish
crosses, directional persistence and recruitment control whether occupancy
changes substantially. `(m, P)` can visualize the recruitment trajectory,
but the direction of motion is needed for its dynamics.

## Exact relation: global and local polarization can disagree

Write `f_R = (1+m)/2`, `f_L = (1-m)/2` and let the chamber-conditional mean
heading vectors be `p_R = P_R e(psi_R)` and `p_L = P_L e(psi_L)`. Then

```
P^2 = f_R^2 P_R^2 + f_L^2 P_L^2
      + 2 f_R f_L P_R P_L cos(psi_R - psi_L).
```

In particular, with perfectly coherent but *oppositely* moving chamber
schools (`P_R = P_L = 1`, heading difference `pi`), **`P = |m|`**. At equal
populations, both schools can be individually perfectly polarized while
global `P = 0`. To evaluate stable local school organization, report `P_L`
and `P_R` or their population-weighted average, not just global `P`.

Conversely, at identical `m = +/-1` and `P = 1`, the occupied school's
heading can point toward or away from the gap, producing different escape
hazards. The mapping from microscopic configurations to `(m, P)` therefore
erases a dynamically causal sign, even before invoking hydrodynamic memory.

## What is actually supported by the committed archive

The `approach_mechanism/tables/nucleation_episodes.csv` archive contains
**512 first-fish-crossing episodes**, 353 of which subsequently reach the
70-fish commitment criterion and 159 of which return to zero opposite-chamber
fish; episodes come from **20 independent simulation seeds**. The comparison
below uses rank AUC and 3,000 bootstrap resamples of whole seeds, not frames:

| Observable at first crossing | Mean committed | Mean failed | AUC | Seed-block 95% CI |
|---|---:|---:|---:|---:|
| Global polarization `P` | 0.820 | 0.730 | **0.553** | [0.484, 0.633] |
| Signed heading toward opening `cos(psi)` | 0.524 | 0.104 | **0.793** | [0.747, 0.836] |
| Wall-normal fish extent `std(x)` | 0.936 | 0.846 | **0.638** | [0.583, 0.689] |

An AUC of 0.5 means chance-level discrimination of success from failure when
comparing an arbitrary committed and failed episode. `P` exhibits a modest
mean difference but ranks success little better than chance; its seed-level
interval includes 0.5.

For the 372 initially highly polarized crossing episodes with `P >= 0.8`,
**96 fail** (25.8%). Crucially, conditioning on this high-`P` subset, the
directional heading still has **AUC 0.788**. Its empirical commitment fraction
is 35/78 = 44.9% for `-0.2 <= cos(psi) < 0.2`, versus 167/186 = 89.8% for
`0.5 <= cos(psi) <= 1`. These are episode fractions at the selected model
parameters, not predictive phase-diagram formulae; unequal bin counts and
seed clustering should be respected in any formal regression.

![Diagnostic, not an inferred potential](figures/imbalance_polarization_diagnostics.png)

### What the population trajectory shows — and does not show

Archived `event_aligned_histories.csv` gives *group-wise mean* polarization
and opposite-chamber counts around 33 strictly matched, censored event pairs.
The diagnostic figure plots the **pair of marginal group means**
`(mean[1-2*N_opp/N], mean[P])` at different relative times. Successful
episodes on average progress toward population balance and continue across;
failed episodes largely return to source occupancy.

**This is not a joint trajectory sample and cannot be converted into**
`rho_ss(m,P)`, a quasipotential, a saddle, or a Markov transition kernel.
It merely documents the projected mean course of a selected episode cohort.

These results also agree with the independent
`approach_mechanism/collective_coordinates/REPORT.md`: angular alignment is a
prospective precursor of nucleation, whereas spatial extension increases
nearer that event. Crucially, that `(psi, Lambda)` pair failed short-lag
closure due to phase and extension-rate memory. Another independent analysis,
`closure_qJ/REPORT.md`, found that population occupancy plus signed neck
current substantially improves held-out splitting probabilities and passage
times compared with occupancy alone, but is not closed at arbitrarily short
lags. That is further evidence against replacing direction by `|P|`.

## Statistical landscape versus physical energy

If full, synchronized long-time samples are available, a histogram can define
an **occupancy landscape**

```
U_occ(m,P) = -log[rho_ss(m,P)] + constant.
```

Its minima indicate frequently visited states, not necessarily dynamical
attractors or the least-action escape path. The microscopic model is actively
driven, hydrodynamically coupled, noisy and generally has nonzero steady
probability currents. Accordingly, `-log(rho_ss)` alone does **not** give a
Kramers barrier or switching rate. A physically predictive quasipotential
would require a defensible closed stochastic generator (or a minimum-action
construction from the actual dynamics) with diffusion/memory handled correctly.

This caveat is the key lesson from Kolpas, Moehlis and Kevrekidis (PNAS 104,
5931, 2007, DOI 10.1073/pnas.0608270104): their effective potential was
obtained **after** validating a dynamically meaningful coarse observable and
estimating its drift/diffusion from microscopic dynamics. See also Brackston,
Wynn and Stumpf (Phys. Rev. E 98, 022136, 2018, DOI
10.1103/PhysRevE.98.022136) for the distinction between a quasipotential
and a simple stationary-probability analogy.

## Next falsifiable test

1. Extract a small, *joint and time-synchronized* derived table from the
   20 available full-BEM trajectories: `seed,time,N_L,N_R,m,P,P_L,P_R,
   cos(psi),J_neck,phase` at the existing output cadence; note how each
   observable is computed and how gap fish are classified. Save only the
   derived table and metadata, not the multi-GB originals.
2. Compute `rho_ss(m,P)`, **clearly labeled an occupancy landscape**. Check
   geometry/symmetry of left and right basins, effect of bin size, and
   seed-to-seed uncertainty. If the already computed `q_l` is exactly related
   to `m`, reuse it to avoid incompatible definitions.
3. At fixed `(m,P)` bins, split by `cos(psi)` and by signed `J_neck`. Measure
   *future*, not simultaneous, first-passage/commitment probabilities to the
   two chamber-dominance sets. If these differ, a 2D Markov closure fails.
   Validate with whole-seed held-out data and ensure matched event selection
   does not determine the outcome by construction.
4. Test Markov memory with several lags and phase/history conditioning; if
   occupancy plus `P` fails, stop treating its landscape as a rate model.
   Compare with the existing `(q_l,J)` baseline rather than inventing
   another fitted double-well or declaring a saddle from visualization.
5. For a *predictive analytical paper*, seek first-principles asymptotic
   descriptions of the pre-cross leading-edge angular escape and the
   post-cross directed recruitment, then connect their event frequencies.
   Do not fit an arbitrary `a(m,P)(1-m^2)^2` potential to reproduce the ridge.

An executable next-stage task for Codex is in `NEXT_CODEX_TEST.md`.

## Reproduce this diagnostic

From the repository root:

```
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/bichamber_p_mpl \
  python3 results/fp_gap0p5_hydro/approach_mechanism/population_polarization_landscape/analyze_existing_evidence.py
```

Outputs: three aggregate CSV tables plus a high-`P` directional split, and
the PDF/PNG diagnostic in `figures/`. No microscopic trajectories were
rerun, no theory coefficients were fitted, and no 2D landscape was inferred.
