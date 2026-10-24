# Phase-conditioned crossing: exact renewal structure and a solvable angular-gate limit

**Status:** analytical iteration, **not** a calibrated explanation of the BEM switching-rate phase diagram. All statements labeled *exact* below hold within their specified state/process or reduced-model definitions. Inferences from the archived simulations are explicitly distinguished from further modeling assumptions.

## 1. What the archived analysis resolves, and what it does not

The archived `R50/gap50` full-BEM trajectories have `N=100`, swimming speed 1, `I_p=4`, angular-noise amplitude `I_n=0.5`, wall strength `I_w=0.94`, `I_f=0.01`, and rounded-tip radius `r=0.5`. **The actual geometry is `gap/R=1`, not `gap/R=0.5`** (`gap/(2R)=0.5`). Only one `(I_w,I_n)` is available in that 20-seed BEM ensemble; a two-parameter phase diagram cannot be established from it.

Earlier audits found (a) a roughly 24--28 time-unit school-height oscillation, but no monotonically decreasing divider distance between cycles; (b) a moving tangent corridor that is strongly associated with imminent crossing; (c) extensive leading-fish identity turnover; (d) no stable wall-induced angular potential barrier; and (e) a precise inward-current identity that establishes propagation by self-swimming, but does not identify an initiating torque. At the first-fish event A, 353 of 512 nucleations eventually reached the operational 70-fish commitment B before the opposite-side population returned to zero. Neither the observed flux nor its torque budget by itself is the desired causal trigger.

There are two scientifically distinct tasks: **(i)** obtain an exact phase-conditioned description of escape and recruitment; **(ii)** find a tractable model that explains whether moderate angular noise facilitates geometric access but excessive noise destroys directional travel. They must not be conflated.

## 2. Exact full-state backward and renewal formulation

Write the full microscopic state as `z=(r_1,...,r_N,theta_1,...,theta_N)`. In each smooth branch of the production wall, visibility, and BEM rules,

\[
d\mathbf r_i=\mathbf v_i(z;I_w)\,dt,\qquad
d\theta_i=\Omega_i(z;I_w)\,dt+I_n\,dW_i,
\]

with \(\mathbf v_i=\mathbf e(\theta_i)+\mathbf U_i^{\rm direct}+\mathbf U_i^{\rm BEM}\) whenever no simulator override is active. The full-state generator is

\[
\mathcal L=\sum_i\mathbf v_i\cdot\nabla_{\mathbf r_i}
+\sum_i\Omega_i\partial_{\theta_i}
+\frac{I_n^2}{2}\sum_i\partial_{\theta_i}^2.
\]

Collision/override surfaces and changes of Delaunay or wall-ray branch require the actual *hybrid* update or its continuous-time interpretation; the smooth-generator expression must not be imposed across these surfaces without matching conditions.

Let \(\Sigma\) be a **prospective** school-height Poincare section, for example an upward crossing of a fixed height with all fish initially in one chamber. Define the *next complete return* to \(\Sigma\) only after visiting a specified opposite-height section; this hysteresis avoids counting infinitesimal/noisy recrossings as distinct cycles. Let \(\tau_B\) be the first 70-fish commitment and \(\tau_\Sigma\) the next section return. For any point \(z\in\Sigma\), define

\[
m(z)=\mathbb E_z[\tau_B\wedge\tau_\Sigma],\qquad
K(z,dz')=\mathbb P_z\{\tau_\Sigma<\tau_B,\ Z_{\tau_\Sigma}\in dz'\}.
\]

The **exact Markov-renewal MFPT equation** is

\[
T(z)=m(z)+\int_\Sigma K(z,dz')\,T(z').
\]

No independent-cycle assumption has been made. This is a high-dimensional equation, **not a solved scalar prediction**. A scalar estimate \(T\approx\langle\tau_{\rm cycle}\rangle/p_B\) needs a demonstrably mixing phase-section distribution and rare, approximately independent successful cycles.

To distinguish initiation from collective recruitment, let A be the first appearance of an opposite-side fish, let \(\tau_0\) be the *subsequent* return to zero opposite fish, and define the exact post-A commitment committor

\[
c(z_A)=\mathbb P_{z_A}(\tau_B<\tau_0).
\]

For the *first* A during a cycle, the probability that this particular nucleation commits before retreat is

\[
p_{A\to B}(z)=\mathbb E_z\!\left[
\mathbf1_{\{\tau_A<\tau_\Sigma\}}c(Z_{\tau_A})\right].
\]

This is an exact strong-Markov identity for the designated first excursion. It need **not** equal the probability of any B before the next return to \(\Sigma\) if unsuccessful A events may be followed by further attempts. Correlated candidate fish must not be multiplied as if independent.

## 3. Exact geometry for a deliberately reduced angular problem

Reflect the left/right origin chambers so positive \(x\) points toward the divider. For a focal fish at normal distance \(d>0\) and height \(y\), set the two rounded-tip centers to

\[
c_-=-g/2-r,\qquad c_+=g/2+r.
\]

The collision-free inward angular interval \([\alpha_-,\alpha_+]\) is bounded by the tip tangents:

\[
\alpha_- =\operatorname{atan2}(c_--y,d)
+\arcsin\frac{r}{\sqrt{d^2+(c_--y)^2}},
\]

\[
\alpha_+ =\operatorname{atan2}(c_+-y,d)
-\arcsin\frac{r}{\sqrt{d^2+(c_+-y)^2}}.
\]

The corridor is **closed** when \(\alpha_-\geq\alpha_+\), and must not then be used. In a temporarily frozen corridor define half-width \(a=(\alpha_+-\alpha_-)/2\), midpoint \(\alpha_c=(\alpha_++\alpha_-)/2\), and centered heading \(\vartheta=\theta-\alpha_c\). Under pure unit-speed self-propulsion, if \(\max(|\alpha_-|,|\alpha_+|)<\pi/2\), continuous residence in this cone until

\[
t_f=\frac{d}{\min_{\theta\in[\alpha_-,\alpha_+]}\cos\theta}
\]

is a **sufficient, conservative** geometric condition for crossing the divider without hitting either tip (the trajectory remains in the cone wedge and has positive inward speed). It is generally **not necessary**: other curved trajectories cross while leaving this frozen cone, and translation from hydrodynamics invalidates the unit-speed guarantee. This bound is only for the isolated angular-geometric model.

## 4. Exact solution of the two-stage, finite-window angular model

Assume angular preparation for a specified positive time \(t_p\) with prescribed constant drift \(b_p\), angular Brownian diffusivity \(D=I_n^2/2\), and initial centered orientation \(\vartheta_0\). At the start of a crossing opportunity, the wrapped-Gaussian density on the circle is **exact** for this preparation model:

\[
p_{\rm prep}(\vartheta;t_p)
=\frac{1}{\sqrt{4\pi D t_p}}
\sum_{k\in\mathbb Z}
\exp\!\left[-\frac{(\vartheta-\vartheta_0-b_pt_p+2\pi k)^2}
{4Dt_p}\right].
\]

For analytical transparency only, take zero angular drift during the subsequent flight and kill the angular process on leaving \((-a,a)\). The **exact** finite-flight angular survival probability solves \(\partial_tS=D\partial_\vartheta^2S\), with \(S(\pm a,t)=0\), \(S(\vartheta,0)=1\), and is

\[
S(\vartheta;t_f)=\frac4\pi\sum_{k=0}^\infty
\frac{(-1)^k}{2k+1}
\cos\!\left(\frac{(2k+1)\pi\vartheta}{2a}\right)
\exp\!\left[-\frac{D(2k+1)^2\pi^2t_f}{4a^2}\right].
\]

Hence the following quantity is exactly solvable by rapidly convergent series and one quadrature:

\[
\boxed{P_{\rm gate}(D)=\int_{-a}^{a}
p_{\rm prep}(\vartheta;t_p)S(\vartheta;t_f)\,d\vartheta.}
\]

Zero drift during the flight is **not essential for solvability**. If a constant flight torque produces angular drift \(b_f\), set \(\beta=b_f/(2D)\) and \(\nu_n=n\pi/(2a)\). The killed survival kernel instead is

\[
S_{b_f}(\vartheta;t_f)=e^{-\beta\vartheta-b_f^2t_f/(4D)}
\sum_{n=1}^\infty A_n\sin\big(\nu_n(\vartheta+a)\big)
e^{-D\nu_n^2t_f},
\]

\[
A_n=\frac{\nu_n\big[e^{-\beta a}-(-1)^ne^{\beta a}\big]}
{a(\beta^2+\nu_n^2)}.
\]

This is an **exact solution only for prescribed constant** \(b_f\) and a stationary cone. Strongly varying social/hydrodynamic torques and the actual moving tip rays invalidate that constant-drift assumption. At very low \(D\), direct evaluation of this eigenseries with nonzero \(b_f\) can be ill-conditioned because of large cancelling exponentials; the code is intended for moderate diffusivities in that case.

For a single self-propelled swimmer in the frozen geometric wedge, this is the probability of a **sufficient angular-success event**, not the full spatial first-passage probability. For the actual school it is a conceptual building block, not a calibrated event rate: candidate fish interact, tips and school geometry move, and the total torque has phase-dependent drift and branch switches.

### Why the optimum exists, and when it disappears

Let \(\mu=\operatorname{wrap}(\vartheta_0+b_pt_p)\). If \(|\mu|>a\) and \(t_p,t_f>0\), then

\[
\lim_{D\to0}P_{\rm gate}(D)=0,
\qquad \lim_{D\to\infty}P_{\rm gate}(D)=0.
\]

The first limit follows because zero-noise preparation misses the cone; the second follows because absorbing angular survival decays at least with the leading rate \(D\pi^2/(4a^2)\). Continuity and positivity for finite \(D\) **prove the existence of an interior optimum** without invoking a potential barrier.

If zero-noise preparation instead places the heading **strictly inside the cone** (\(|\mu|<a\)) and there is no flight drift, then \(P_{\rm gate}(D)\to1\) as \(D\to0\). In that regime an interior maximum is not required; the toy model favors zero noise. Thus a wall-dependent *preparation drift* can move the preferred noise strength or eliminate the optimum altogether.

When zero-noise preparation misses the cone by \(\Delta=|\mu|-a>0\), balancing the small-noise Gaussian-entry exponent \(\Delta^2/(4Dt_p)\) against the large-noise principal-survival exponent \(D\pi^2t_f/(4a^2)\) suggests

\[
D_{\rm balance}\sim\frac{a\Delta}{\pi\sqrt{t_pt_f}},\qquad
I_{n,\rm balance}\sim\sqrt{2D_{\rm balance}}.
\]

**This is an uncontrolled balance estimate**, not an asymptotically exact expression for the maximizer of the full series. It gives a geometric interpretation of the trend: a narrower gate, greater initial mismatch, slower flight or shorter preparation changes the preferred noise amplitude. The numerical series, not this estimate, must be used for quantitative toy-model predictions.

### Numerical evaluation (mathematical example, not data fitting)

The accompanying `finite_window_gate.py` evaluates the exact wrapped heat kernel and absorbing eigenfunction expansion. Using the archive's geometric lengths `g=7.1589105`, `r=0.5`, and an **illustrative** isolated-fish location `(d,y)=(2.34,0)` gives `a=0.94349` rad and the conservative pure-swimming flight time `t_f=3.99`. Assign an illustrative preparation duration `t_p=2` and initial angular mismatch of 0.5 rad *outside* the cone. These timing and mismatch choices are **not extracted as microscopic torque data or a measured leading-fish distribution**.

| assigned constant preparation drift \(b_p\) | exact-series maximizing \(I_n\) | maximum \(P_{\rm gate}\) | interpretation |
|---:|---:|---:|---|
| 0.0 rad/time | 0.4262 | 0.04394 | finite-noise optimum |
| -0.2 rad/time | 0.2958 | 0.13981 | deterministic steering nearly reaches the cone |
| -0.4 rad/time | tends to 0 | tends to 1 | zero-noise heading is safely inside the cone |

For the first row, independently increasing the angular quadrature resolution from 128 to 1024 nodes and eigenmode count from 30 to 150 changes `P_max=0.0439372709319` by less than displayed precision. The balance estimate is `I_n=0.3261`, whereas the exact series gives `0.4262`, illustrating why the balance estimate must not be reported as a precise ridge law. The numerical closeness of any illustrative result to the archive's `I_n=0.5` is **not validation**; the actual phase-conditioned drift, candidate location, opportunity time, and recruitment probability were not measured here.

## 5. Connection to a full first-principles theory

The microscopic wall rule and the moving tip tangents can determine the *instantaneous* angular gate \(\alpha_\pm(z)\) and angular drift \(\Omega_i(z;I_w)\) with **no fitted coefficient**. The mathematically appropriate spatial-angular reduction would include at least \((d,y,\theta,\phi)\), with school oscillation phase \(\phi\) and the actual deterministic velocity and torque evaluated in a controlled background-state closure:

\[
\dot d=-v_x(d,y,\theta,\phi;I_w),\quad
\dot y=v_y(d,y,\theta,\phi;I_w),\quad
d\theta=b(d,y,\theta,\phi;I_w)dt+I_n dW,
\quad \dot\phi=\omega(\phi).
\]

The backward first-passage equation in the smooth part of the reduced domain is

\[
-v_x\partial_d q+v_y\partial_y q+b\partial_\theta q
+\omega\partial_\phi q+\frac{I_n^2}{2}\partial_\theta^2q=0,
\]

with **hitting** boundary `q=1` on a valid opposite-chamber crossing and `q=0` on an explicitly defined unsuccessful phase reset. This equation is exact for a specified four-coordinate process, but reducing the full interactive school to that process is a **closure approximation**, and an imposed phase trajectory may fail in real crossings. A static potential or scalar net gap flux is neither necessary nor justified.

In a small-noise regime *conditional on a deterministic orbit that does not already cross*, a more faithful non-gradient landscape is the controlled **path action**

\[
\mathcal S_* =\inf_{T,u_i}\left\{\frac12\int_0^T
\sum_i u_i(t)^2\,dt:\quad
\dot{\theta}_i=\Omega_i(Z;I_w)+u_i,\quad
\dot{\mathbf r}_i=\mathbf v_i(Z;I_w),\quad
\text{first crossing before phase reset}\right\}.
\]

On regular smooth branches, Freidlin-Wentzell reasoning gives a conditional rare-event exponent of order \(e^{-\mathcal S_*/I_n^2}\). This is a **path-dependent action, not a potential of one coordinate**. The action is zero if the unforced deterministic trajectory already crosses; ray-switch and neighbor-graph discontinuities require particular care. High-noise persistence loss and post-A recruitment must still be modeled separately to obtain a full optimum ridge.

## 6. Single decisive falsification before further theory growth

The archived full-BEM ensemble has only one wall/noise pair, and previous close-to-A retrospective replays do not reveal initiation. Rather than another fitted coordinate search or large BEM torque replay, conduct **one bounded prospective check** on previously saved trajectories:

1. Anchor candidates at the same upward or downward full-occupancy height section (before knowing whether they will cross), never at `A-2` selected from future events.
2. For every *actual first crosser*, track the time-dependent rounded-tip cone, the angular residence and the remaining physical travel time from its first cone entry until A. Compare with same-phase at-risk noncrossers using independent seeds for uncertainty.
3. Check whether an appreciable fraction of actual first crossings meet a *persistent angular corridor* condition with sufficient time to traverse. If few do, the solvable frozen-cone model is at most a severe lower bound and must **not** be promoted to the mechanism. Move directly to the full kinetic `(d,y,theta,phase)` first-passage calculation.
4. Independently, in a small predeclared sample of phase-section states, integrate the **unmodified deterministic production dynamics** with `I_n=0` for one next cycle. Count whether genuine first crossings occur. If they commonly do, there is no positive small-noise activation action from those states; don't invoke Kramers-type scaling.

Only if these checks support a limited set of slowly varying preparation/flight variables should \(b_p(I_w)\), \(a(I_w)\), \(t_p(I_w)\), and \(t_f(I_w)\) be computed from microscopic laws and inserted into the series **without adjusting them to reproduce the phase diagram**. The resulting noise optimum and wall-strength ridge can then be prospectively confronted with simulation. The recruitment committor must be measured or derived separately for a 70-fish switch; one fish passing through the opening is not the whole-school transition.
