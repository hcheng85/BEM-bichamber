# Proposed collective landscape: population imbalance and chamber-local polarization

**Status:** A physics-motivated *hypothesis to test*, not a derived escape
potential or a confirmed Markov closure. This document supersedes the
coordinate choice in `NEXT_CODEX_TEST.md` as originally committed, while the
earlier `REPORT.md` remains the evidence baseline.

## 1. Physical definitions

For `N` swimmers, let `N_L` and `N_R` denote the number assigned to left and
right chambers. Define a consistent convention for aperture occupants; if
aperture fish form a distinct class, track `N_gap` and fluxes separately.

The population imbalance is

```
m(t) = [N_R(t)-N_L(t)]/N,     -1 <= m <= 1.
```

Thus `m=+1` means the right chamber is fully occupied, `m=-1` the left, and
`m=0` equal populations. The *sign of m* identifies the currently dominant
chamber. The *sign of the increment*

```
Delta m(t;h) = m(t+h)-m(t)
```

identifies net rightward/leftward population transfer over the interval. With
an exhaustive two-chamber assignment, the following is an exact counting
identity:

```
m(t+h)-m(t) = (2/N) [C_{L->R}(t,t+h)-C_{R->L}(t,t+h)].
```

No separately recorded neck-current variable is *necessary* to determine the
net transfer direction retrospectively from a resolved `m(t)` series.
Nevertheless `m(t)` alone contains neither its own increment nor the
incoming fish headings, so a projected model may retain memory.

For each nonempty chamber, let `e_i = (cos(theta_i),sin(theta_i))` and define

```
p_L = (1/N_L) sum_{i in L} e_i;   P_L = |p_L|,
p_R = (1/N_R) sum_{i in R} e_i;   P_R = |p_R|.
```

`P_c=1` signals perfect heading alignment of the fish in chamber `c`;
`P_c~0` signals broadly dispersed headings. **Local coherence does not
require the two chamber groups to point in the same direction.** Importantly,
`P_c` is *undefined* when `N_c=0`, rather than physically zero. Record a
missing-value mask and both populations. Near the one-chamber basins, plot
the coherence of the occupied group; do not fabricate the empty group's
polarization, and do not silently remove all `|m|=1` configurations.

## 2. Hypothesized order-parameter landscape

We first examine the state

```
z(t) = (m(t), P_L(t), P_R(t)),
```

where both polarizations are defined, and its natural lower-dimensional
one-chamber boundary strata (`m=+/-1`, occupied-group polarization only).
Chamber exchange under a perfectly symmetric model maps
`(m,P_L,P_R)` to `(-m,P_R,P_L)`; measure any simulation asymmetry rather than
assuming it away. The chamber-dominated metastable states should occur near
the two occupancy boundaries **if** the observations show persistent
residences there. Intermediates describe progression and possible temporary
loss of coherence during recruitment.

The empirically measurable *occupancy landscape* is

```
U_occ(z) = -log rho_ss(z) + constant,
```

on statistically supported bins of each occupancy stratum. For the active,
hydrodynamic, non-equilibrium school, this is **not** automatically a
mechanical energy, Freidlin-Wentzell quasipotential, or a rate barrier. Its
stationary marginal leaves out probability currents and motion history.

No inverse scaling law `switching_rate ~ 1/|m|` follows: `|m|` mostly tracks
how much of an already-started transfer has occurred, whereas the waiting
hazard at `|m|=1` describes *initiation* of a transfer. Treat the latter and
post-entry recruitment separately if the reduced state cannot capture both.

## 3. Critical information-loss test

Two configurations can share exactly the same `(m,P_L,P_R)` but have the
coherent population swim toward or away from the gap. The local scalar
polarizations do not encode the signs of the chamber heading vectors.
Likewise `m=0.4` can be moving toward `m=1` or toward `m=-1`.

Therefore a **closed Markov model in these three magnitudes is a falsifiable
approximation, not a consequence of their definitions**. Test it directly:

1. At fixed or closely matched `(m,P_L,P_R)`, compare held-out probabilities
   of reaching the right or left chamber-dominated set first.
2. Stratify the *same* bins by the **past**, observable sign of `Delta m`,
   and separately by chamber-conditional signed projections
   `h_L=p_L dot n_(L->gap)` and `h_R=p_R dot n_(R->gap)` (when defined).
   Only use information available at the prediction time. A substantial
   change in the future commitment probability means the proposed three
   coordinates are not an autonomous kinetic state.
3. Verify time-lag stability, approximate Chapman-Kolmogorov consistency and
   held-out first-passage times/survival curves with proper censoring.
4. Use the earlier occupancy-plus-neck-current `(q_l,J)` study as an
   independent benchmark. This is a comparison/diagnostic, **not** a demand
   to include `J` in the proposed landscape from the outset.

The existing evidence makes this a serious test: across 512 first-crossing
episodes, global `P` discriminates commitment only weakly (`AUC=0.553`),
whereas opening-directed heading gives `AUC=0.793`. Chamber-local `P_L,P_R`
are better measures of *separate school coherence* but their magnitudes
still omit the sign of each group's direction. This prior result neither
confirms nor rejects the three-coordinate proposal, which has not yet been
measured on joint synchronized time series.

## 4. A possible coarse dynamics, conditioned on validation

For `N=100`, the population jumps occur in increments `2/N=0.02` under
the exhaustive assignment. A finite-state transition kernel or a
jump-process description of `m` is initially more faithful than assuming
smooth overdamped diffusion. Polarization can vary continuously and may
possess inertial/oscillatory memory from the full swimming dynamics.

If a coarse-lag Markov test *passes*, estimate

```
K_tau(z'|z) = Pr[z(t+tau) in dz' | z(t)=z],
```

with seed-held-out validation and uncertainty. Its committor solves the
appropriate discrete backward equation on the chosen chamber-dominance
sets. First-passage distributions can then be predicted from that operator.
Estimate drift/diffusion and a Fokker-Planck limit only after jump size and
time-lag convergence have been demonstrated; do not infer a Kramers rate
from a fitted double-well plot.

If Markov closure *fails*, preserve `(m,P_L,P_R)` as scientifically useful
order parameters, but do not present `U_occ` as a predictive theory. The
first hypothesis to check is that **directional heading/history** is missing;
then seek a physical reduction of initiation near the gap and directional
post-entry recruitment from the microscopic wall/social/hydrodynamic laws.

The longer-term paper goal is a parameter-free, or explicitly
controlled-asymptotic, explanation of the *wall-strength × angular-noise*
switching-rate diagram and optimal ridge. A single-parameter-point density
plot or empirical transition surface does not meet that goal.

## 5. Required evidence and limitation

The committed `nucleation_episodes.csv` contains first-crossing event-level
global polarization, but **not the synchronized local `(P_L,P_R)` time
series**. The separate event-aligned mean population and global polarization
curves are not joint measurements. Reconstructing a legitimate three-state
occupancy landscape requires the underlying 20 full-BEM trajectories (or
precomputed, *synchronous* derived variables with the same time stamps).

Keep the original 20-seed data split (training seeds 0--14, held-out seeds
15--19), discard the same initial 200 time units, and use the same transition
event definitions as `approach_mechanism/REPORT.md`: first fish entering the
opposite chamber is A; 70 opposite-chamber fish before returning to zero is
successful commitment B. Respect episode censoring and whole-seed
dependence. The implementation requirements are in `NEXT_CODEX_TEST.md`.

### Methodological references

- Kolpas, Moehlis & Kevrekidis, *PNAS* **104**, 5931 (2007), DOI:
  `10.1073/pnas.0608270104`: dynamical coarse-observable testing and
  microscopic estimates of effective drift/diffusion before constructing a
  collective potential.
- Brackston, Wynn & Stumpf, *Physical Review E* **98**, 022136 (2018), DOI:
  `10.1103/PhysRevE.98.022136`: distinction between stationary occupancy
  analogies and the dynamics-based quasipotential.
