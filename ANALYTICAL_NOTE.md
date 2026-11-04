# Exact circular reduction of the original forward-ray wall law

## Scientific purpose and scope

Replace the bichamber boundary by a disk of radius `R`, while preserving the
*physical forward-ray wall avoidance* in the existing BEM fish simulator. Start
with the **single-fish, wall-only subsystem** before reintroducing the actual
Delaunay-neighbor social alignment/attraction and BEM hydrodynamics. This is a
mechanistic benchmark, **not** a claim that independent fish describe the
interacting school or its two-chamber transitions.

**Prior-work baseline:** Huang, Ling and Kanso, *Collective phase transitions
in confined fish schools*, PNAS **121** (2024), DOI
`10.1073/pnas.2406293121`, already studied the **full interacting fish model
in a circle** (including geometric boundary and hydrodynamics), found
schooling/milling and other collective states, and used an empirically inferred
one-dimensional stochastic collective-coordinate potential to characterize
intermittency. We must not present circular confinement or the existence of
collective milling as new. The complementary analytical objective here is a
**microscopically derived wall-induced orbit, its linear stability and noise
covariance**, followed by justified collective moment/phase reductions. The
exact single-fish orbit does not yet explain the full school's phase diagram.

The original simulator has unit self-propulsion speed, angular noise amplitude
`In`, and wall strength `Iw`. Its unmodified, away-from-contact rule is

\[
\dot{\mathbf x}=\mathbf e(\theta),\qquad
d\theta=\frac{I_w\,\mathrm{sgn}(\theta-\beta)}{d}\,dt+I_n\,dB_t,
\]

where `d` is distance to the *first forward ray intersection*, and `beta` is
the outward normal there. Interpret the signed angular difference wrapped into
`(-pi,pi]`. In the disk it always lies in `[-pi/2,pi/2]` for interior points.

**Caveat:** the actual code multiplies wall torque by
`1 + 1/d**(high_order-1)` inside its `|d| <= |safedis|` safety layer, and clips
angular velocities at `pi/dt`. These special near-contact/numerical safeguards
are NOT part of the exact smooth-branch reduction below. The source implementation
is `results/fp_gap0p5_hydro/approach_mechanism/pre_crossing_turning/scripts/replay_dynamics.py`
(`ray_wall_intersection`, `wall_avoidance`). The circular expressions are exact
for the ordinary branch outside that layer.

## Step 1: remove the spatial azimuth exactly

Let `r=|x|`, `varphi=arg(x)`, and `alpha=theta-varphi` be the heading relative
to the radial direction. The ray-wall distance, found from
`|x+d e(theta)|^2=R^2`, is

\[
d(r,\alpha)=-r\cos\alpha+\sqrt{R^2-r^2\sin^2\alpha}.
\]

Angular momentum geometry gives
`R sin(theta-beta)=r sin(alpha)` and hence
`sgn(theta-beta)=sgn(sin(alpha))` for `r>0`. The angular position of the fish
obeys `dot(varphi)=sin(alpha)/r`. Thus the **exact** two-dimensional Markov
subsystem, with neither a fitted torque nor an effective potential, is

\[
dr=\cos\alpha\,dt,\qquad
d\alpha=\left[\frac{I_w\,\mathrm{sgn}(\sin\alpha)}{d(r,\alpha)}
-\frac{\sin\alpha}{r}\right]dt+I_n\,dB_t.
\tag{1}
\]

This coordinate chart is singular at `r=0`; the full Cartesian ray rule is
regular there except that its wall torque vanishes for an exactly radial
center departure (`theta-beta=0`). Such exact radial initial conditions are
non-generic and require care with the simulation's boundary safeguard.

## Step 2: exact wall-induced circulating solutions

With `In=0`, two circular solutions have constant radius and tangent heading:

\[
\alpha_*=\pm\frac\pi2,\qquad
r_* = \frac{R}{\sqrt{1+I_w^2}},\qquad
d_* = \frac{R I_w}{\sqrt{1+I_w^2}},\qquad
\dot\varphi_* = \dot\theta_* = \pm\frac{\sqrt{1+I_w^2}}R.
\tag{2}
\]

These are parameter-free predictions of the original wall rule. Stronger
repulsion moves the milling orbit inward. The smooth rule is applicable at
the orbit only when `d_* > |safedis|`; the simulator's clipping must also be
inactive (`sqrt(1+Iw**2)/R < pi/dt`). For small `Iw`, the formal orbit approaches
the wall and safety regularization can change it.

## Step 3: local stability is exactly solvable

Let `a=Iw>0`, `s=sqrt(1+a**2)`, and perturb a counterclockwise orbit by
`u=r-r_*` and `v=alpha-pi/2`. To first order,

\[
\begin{pmatrix}\dot u\\\dot v\end{pmatrix}
=\begin{pmatrix}0&-1\\C&-\gamma\end{pmatrix}
\begin{pmatrix}u\\v\end{pmatrix},\qquad
C=\frac{(1+a^2)^2}{R^2a^2},\quad
\gamma=\frac{\sqrt{1+a^2}}{Ra}.
\]

The clockwise version changes the two off-diagonal signs; the eigenvalues
are the same:

\[
\lambda_\pm=-\frac\gamma2\pm\frac{i}{2}\sqrt{4C-\gamma^2}
=-\frac{\sqrt{1+a^2}}{2Ra}
\pm i\frac{\sqrt{(1+a^2)(4a^2+3)}}{2Ra}.
\tag{3}
\]

Both circulating branches are therefore **locally attracting spirals for
every `a>0`**, within the smooth forward-ray model. This is a wall-generated
mill, not an assumed double-well potential. No global-basin or uniqueness
claim is made; the radial configurations and center need special treatment.

## Step 4: small-angular-noise predictions without fitting

Near either circulating branch the linearized stochastic problem is a 2D
Ornstein-Uhlenbeck process driven only in the `v` component. Its stationary
covariance follows from the Lyapunov equation, giving

\[
\operatorname{Var}(\alpha-\alpha_*)
=\frac{I_n^2}{2\gamma}
=\frac{I_n^2 Ra}{2\sqrt{1+a^2}},\qquad
\operatorname{Var}(r-r_*)
=\frac{I_n^2}{2\gamma C}
=\frac{I_n^2R^3a^3}{2(1+a^2)^{5/2}},\qquad
\operatorname{Cov}(r,\alpha)=0.
\tag{4}
\]

These variances are **local weak-noise predictions**, not global moments if
noise drives switching between the two chiral branches or into the near-wall
regularization. The angular azimuth remains a free global phase.

## Step 5: the correct 2D stochastic first-passage problem

Write `f(r,alpha)=cos(alpha)`,
`g(r,alpha)=Iw*sgn(sin(alpha))/d(r,alpha)-sin(alpha)/r` and
`Dtheta=In**2/2`. For a density with respect to `dr dalpha`,

\[
\partial_t p=-\partial_r(fp)-\partial_\alpha(gp)
+D_\theta\partial_\alpha^2p.
\tag{5}
\]

For a physically specified target set `A` (for example a change in signed
circulation, defined away from the radial/center singularities), the backward
mean first-passage problem is

\[
f\,\partial_r T+g\,\partial_\alpha T
+D_\theta\partial_\alpha^2T=-1,\qquad T|_A=0.
\tag{6}
\]

Boundary conditions must be supplied **from the actual wall rule**, including
its near-contact regularization, rather than silently imposing specular
reflection. The noise is angular, so a scalar radial Kramers potential is not
generally justified.

In a weak-noise WKB approximation `p~exp[-V/Dtheta]`, the leading-order
quasipotential satisfies the **two-dimensional** Hamilton-Jacobi equation

\[
f\,\partial_r V+g\,\partial_\alpha V
+(\partial_\alpha V)^2=0.
\tag{7}
\]

Even this is a *proposed* rare-transition calculation: the non-smooth
chirality boundary and simulator's contact rule require special handling.
We should not infer a barrier from stationary population imbalance.

## How to extend this to the actual interacting school

There is an **exact, unclosed moment hierarchy** for the actual school in the
circle; it provides the correct bridge from one fish to `N` fish and avoids
introducing a neural-network surrogate. Let

\[
\mathbf p=\frac1N\sum_i\mathbf e(\theta_i),\qquad
M=\frac1N\sum_i\sin\alpha_i,\qquad
\alpha_i=\theta_i-\arg\mathbf x_i.
\]

Here `M` is *signed milling about the fixed circle center* (not necessarily
the school's instantaneous center of mass). With the full simulator's
`d theta_i=Omega_i dt+In dB_i`, where `Omega_i` includes wall, social and
hydrodynamic angular terms, Itô's formula gives

\[
d\mathbf p=\left[\frac1N\sum_i \mathbf e_\perp(\theta_i)\Omega_i
-\frac{I_n^2}{2}\mathbf p\right]dt
+\frac{I_n}{N}\sum_i\mathbf e_\perp(\theta_i)dB_i,
\tag{8}
\]

and, with the **actual** translational velocity `v_i` (including the BEM
velocity if retained),

\[
dM=\left[\frac1N\sum_i\cos\alpha_i
\left(\Omega_i-\frac{\mathbf x_i\times\mathbf v_i}{r_i^2}\right)
-\frac{I_n^2}{2}M\right]dt
+\frac{I_n}{N}\sum_i\cos\alpha_i\,dB_i.
\tag{9}
\]

These are exact away from the circle center and discrete safety overrides,
assuming independent Brownian heading noises. They show explicitly why
polarization and milling are *not* determined by population count or wall
strength alone: the unresolved torque--heading and velocity--geometry
correlations matter. The terms can be calculated directly from simulated
microstates; a subsequent closure must justify their approximation.

1. Keep the same wall function for **every fish** in a disk. Reinstate the
   simulator's exact social alignment and attraction on its moving visible
   Delaunay graph. Turn off direct and boundary hydrodynamics initially, then
   add them separately; do not substitute an unnormalized Kuramoto network.
2. Test whether individual wall-induced rotating branches persist as a
   *collective* milling state. The mean polarization `P` and signed milling
   `M=(1/N)sum_i [(x_i/|x_i|) cross e(theta_i)]` are distinct:
   `P` detects aligned movement, while `M` detects coherent circling.
3. Examine the competition of **wall-induced curvature**, coherence/attraction
   from true social torques, and angular noise. Derive linear modes of any
   symmetric collective solution before guessing an effective potential.
4. Only after the circular test identifies genuine stable macrostates and
   switching mechanisms return to bichamber transitions. A single disk has
   **no left/right population-imbalance transition**; its appropriate basic
   targets are milling stability, chirality reversal, and schooling/milling
   transitions.

## Reproducibility / falsification

`validate_orbits.py` independently integrates the **Cartesian forward-ray
wall rule**, compares orbit radius and angular frequency with (2), and
confirms linear covariance (4) by an algebraic Lyapunov solver (not by noisy
Monte Carlo). It writes `orbit_validation.png`. This validates the smooth
single-fish branch; it does not validate a collective transition mechanism.
