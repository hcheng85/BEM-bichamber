# Exact planar-wall reduction and a falsification of a simple tip-optimum mechanism

Status: 24 October 2026. **Verified results for two restricted problems, not yet a full-school explanation or a prediction of its rate optimum.** No data-fitted constants.

## 1. What the archived model actually says at the rounded tips

The full BEM simulator uses a heading-dependent ray intersection, then applies
`omega_wall = Iw * sign(wrapped_theta_minus_normal) / ray_distance`, with an
additional near-wall factor and clipping where applicable.
The exact BEM focal replay in
`../divider_tip_geometry/scripts/exact_tip_dynamics.py` provides the one-sided
angular drift at the geometrical upper/lower tip corridor.

`test_wall_drift_continuation.py` reconstructs the archived one-sided drifts at
`Iw=.94` to a maximum absolute error of **3.55e-15** and varies only the wall
strength in each *fixed* configuration. Every one of the 5,012 eligible
`epsilon=.001` tip evaluations is in an unclipped branch through `Iw=2`.

On these frozen configurations, the derivative of the **outward-edge drift**
toward the corridor with respect to `Iw` is **negative in all 5,012 cases**:
min `-4.80469`, median `-0.24458`, max `-0.09047` (in the normalized
outside-positive convention). The percentage with locally inward edge drift
increases monotonically `52.04%, 54.59%, 57.16%` at `Iw=0, .94, 2`.
The stricter transmission classification increases `51.98%, 52.67%, 53.33%`
while pinning also grows `0.06%, 1.92%, 3.83%`.

Thus the **immediate, frozen-state wall-induced angle change at the rounded tip
cannot produce an interior optimum over 0<=Iw<=2**. It can *increase* pinning,
but the reported transmission fraction does not fall. This does **not** rule
out a dynamic tip mechanism once the whole school's state distribution is
allowed to change, or effects beyond the checked strength range. The 5,012
evaluations share seeds, fish and episodes; they are not independent replicates.

## 2. Exact first-principles reduction on the flat divider

Use dimensionless self-propulsion speed `U=1`. In the *left chamber*, let the
flat solid divider be `x=-r`, the wall distance be `d=-r-x>0`, and take a
heading `theta` directed toward the wall, `-pi/2<theta<pi/2`. Define
`beta=abs(theta)` on an interval where the sign of `theta` is fixed.
Provided the **selected ray** intersects this particular straight divider,
is not amplified and is not clipped, its intersection distance is

    ell_ray = d / cos(beta).

The exact original wall rule and swimming kinematics then give, in the
**wall-only model**,

    dd = -cos(beta) dt,
    d(beta) = [Iw cos(beta)/d] dt + In dW_beta.

There is no translation Brownian noise in this model. Crucially, the change
of variable

    H = beta + Iw log(d)

cancels the deterministic nonlinear wall torque **exactly**:

    dH = In dW_beta.                              [1]

The sign change of the driving Wiener increment from `theta` to `beta`
does not affect its Brownian law while the angle sign is fixed.
No Ito correction is needed because `d` has finite variation.
This is a *stochastic first integral of the locally ray-cast wall dynamics*.

For `In=0`, the invariant predicts

    beta(d) = beta0 + Iw log(d0/d).

For a nonzero wall strength, the ideal infinite planar branch has an
asymptotic tangential stand-off distance

    d_star = d0 exp[-(pi/2-beta0)/Iw].             [2]

For increasing `Iw`, `d_star` increases strictly: the swimmer becomes
wall-parallel farther from the divider. This is **not** a finite-time
turnaround, and for a finite wall the actual ray eventually switches to a
rounded tip or the outer boundary.

If there are other deterministic angular torques `Omega_other` and an
additional wall-normal translational velocity component `u_x` (social,
direct and BEM hydrodynamics), then before angular clipping,

    dH = [sgn(theta)*Omega_other - Iw*u_x/d]dt
         + In*sgn(theta)dW.                     [3]

Consequently equation [1] is exact only for the isolated planar wall and
provides a **specific, source-derived residual** to test in the collective
simulation. Do not assume it is small because `Ip=4` and hydrodynamics are
active in the archived system.

## 3. A rigorous noise-dependent upper bound, not an escape probability

Consider the planar-branch subproblem started at `(d0,beta0)`. Let `dc<d0`
be a prescribed distance required to enter an upstream layer, with the
flat-wall ray and angle branch intact until first reaching `dc`. If the fish
ever reaches `dc` in time `T` with `0<=beta<pi/2`, then necessarily

    H <= Hc = pi/2 + Iw log(dc).

The required *downward* change of the Brownian label is

    Delta = H0-Hc = Iw log(d0/dc) - (pi/2-beta0).

When `Delta>0`, the reflection principle applied to [1] proves

    P(reach dc by T while retaining the planar branch)
      <= erfc[Delta / (sqrt(2)*In*sqrt(T))].    [4]

When `Delta<=0`, [4] gives only the trivial upper bound `1`.
For deterministic `In=0`, a **necessary** condition for finite-time approach
is

    Iw < Iw_crit = (pi/2-beta0) / log(d0/dc).  [5]

At equality the local planar solution approaches `dc` asymptotically;
for larger strengths it cannot reach `dc` without leaving the branch.
Equation [4] is the exact *first-passage probability of the reduced label*
and an **upper bound** on physical proximity, NOT an equality for arrival at
the tip. Merely crossing `Hc` does not place a fish at `dc`: position and
remaining travel time still matter. If collective terms, branch changes or
clipping arise first, the physical use of this bound ceases.

This yields a candidate scaling for the flank of a blocked-access region:
`In sqrt(T)` versus `Iw log(d0/dc) - (pi/2-beta0)`. It has no adjustable
coefficients but it is **not** the observed school's two-parameter
transition-rate phase diagram.

## 4. Verification

Commands, from repository root:

    python3 results/fp_gap0p5_hydro/approach_mechanism/deterministic_gate_threshold/test_wall_drift_continuation.py
    MPLCONFIGDIR=/tmp/bichamber_mpl python3 results/fp_gap0p5_hydro/approach_mechanism/deterministic_gate_threshold/analyze_planar_access.py

The second command independently re-evaluates the true source wall/raycast
for >=30 geometrically eligible points to machine accuracy, confirms the
exact cancellation, verifies all 5,012 tip slope signs and compares a
fixed-time Gaussian-label random walk to the reflection-principle result.
The archived wall-only Euler trajectory
`../divider_tip_geometry/tables/planar_wall_identity.csv` independently
verifies the same planar invariant, max error about `1.6e-6` with Euler
time discretization.

Figure: `figures/competing_wall_mechanisms.pdf` and `.png`. Data:
`tables/counterfactual_one_sided_Iw.csv`.

## 5. Physics and explicit falsification criteria

- **Established:** Tip-wall angular bias is locally facilitating and
  increases monotonically with `Iw` on the fixed archived near-gap states.
- **Established:** On the flat straight-divider ray branch, a conserved
  deterministic heading-distance label becomes pure Brownian motion under
  rotational noise. Stronger wall interaction increases the local stand-off
  distance and can block approach to a chosen near-wall layer without noise.
- **Not established:** Farther stand-off means less gap crossing. Stronger
  wall steering also accelerates motion *along* the divider toward its tip.
  These are opposing geometric consequences whose net effect depends on
  finite geometry, arrival phase, neighbors and hydrodynamics. A physical
  passage rate cannot be inferred by multiplying [4] by a frozen tip fraction.
- **Not established:** The source of the first school-directed turn; whether
  the residual in [3] is a deterministic alignment/hydrodynamic impulse,
  noise crossing a threshold, or a finite-geometry ray switch.
- **Not established:** The observed wall-strength optimum, the
  wall/noise optimal ridge, or a full-school transition MFPT.

The decisive next test is *not* another projection of swimming direction:
calculate the time-resolved budget [3] at a prospectively selected section
of the oscillating school's cycle, before knowing whether it crosses, and
use independent short counterfactual rollouts to check whether changes in
`H` actually trigger an approach and subsequent school-wide transition.
The full BEM states and LU matrices required for that test are not present
in this repository; see `NEXT_CODEX_TEST.md`.

## Literature connection (method, not a claim of equivalence)

The derivation follows exact coordinate reduction and first-passage methods,
not an inferred effective potential. The **collective** decision about the
gap should then be checked against the bifurcation/collective-variable
strategy developed for noisy schooling transitions (Chan and Kanso,
*Bioinspiration & Biomimetics* **21**, 016024, 2026,
DOI `10.1088/1748-3190/ae39bd`), and against confinement-induced
collective escape (Aranson and Pikovsky, *Phys. Rev. Lett.* **128**, 108001,
2022). These papers study different models; neither provides a validated
closure for the present bichamber's hydrodynamic 100-fish dynamics.
