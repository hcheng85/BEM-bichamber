# Independent aperture geometry and a parameter-free first-exit problem

## Scientific conclusion

The aperture filter contributes information beyond the removal of negative velocity contributions. Positive inward supply alone is not sufficient to explain the previously observed advantage of the capture current.

On held-out seeds, `Jcap` outperforms `Jin` for future first entry and four-fish transfer. The gain remains resolved after restricting comparisons to the existing `(q,a,b)` cells and then further to outcome-blind, training-only `Jin` quartiles within those cells. In the most controlled comparison:

- first entry: `AUC(Jcap)-AUC(Jin)=0.01863`, seed-bootstrap 95% CI `[0.00958, 0.02891]`;
- four-fish transfer: gain `0.01532`, CI `[0.00953, 0.02492]`.

By contrast, `Jin-Ja` is unresolved in every corresponding comparison. The supported physical quantity is therefore **positive inward supply filtered by actual aperture accessibility**, not merely a rectified longitudinal current.

This supports a parameter-free first-exit description of a single approach attempt. It does not establish a closed collective theory because the social turning rate evolves with the full many-fish configuration.

## The three currents

All quantities use the same source-side approach region and the current majority chamber; no future direction enters their definitions. With local forward velocity `u_i=s cos(theta_i)`:

\[
J_a={1\over N}\sum_{i\in A}s\,v_{x,i},
\]

where `v_x` is the archived stored velocity, including signed outward contributions and collision overrides;

\[
J_{\rm in}={1\over N}\sum_{i\in A}\max(u_i,0);
\]

and

\[
J_{\rm cap}={1\over N}\sum_{i\in A}\max(u_i,0)
I\!\left(\left|y_i-z_i{\sin\theta_i\over u_i}\right|\le h\right).
\]

Thus `Jin` and `Jcap` use the same heading velocity and differ only by the exact straight-ray aperture indicator. Their comparison isolates geometry from sign rectification.

## Held-out empirical comparison

Seeds 0--29 define only the pre-existing `b` split and the outcome-blind `Jin` quartiles. Seeds 30--39 remain held out. Uncertainty resamples the ten original held-out seeds as clusters. No coefficient or outcome probability is fitted.

### Future first entry

| comparison | `Ja` AUC | `Jin` AUC | `Jcap` AUC | `Jcap-Jin`, 95% CI |
|---|---:|---:|---:|---:|
| all held-out origins | 0.657 | 0.677 | 0.707 | 0.0293 `[0.0081, 0.0532]` |
| same `(q,a,b)` cell | 0.554 | 0.560 | 0.599 | 0.0387 `[0.0196, 0.0557]` |
| same `(q,a,b,Jin)` stratum | 0.506 | 0.510 | 0.529 | 0.0186 `[0.0096, 0.0289]` |

The final comparison uses 10,628 positive-negative pairs across 19 supported held-out cells. Pair-weighted positive-minus-negative `Jin` is only `0.00160`, confirming that the geometric contrast is evaluated at closely comparable inward supply.

### Four-fish transfer

| comparison | `Ja` AUC | `Jin` AUC | `Jcap` AUC | `Jcap-Jin`, 95% CI |
|---|---:|---:|---:|---:|
| all held-out origins | 0.732 | 0.720 | 0.762 | 0.0428 `[0.0254, 0.0620]` |
| same `(q,a,b)` cell | 0.645 | 0.643 | 0.667 | 0.0236 `[0.0071, 0.0462]` |
| same `(q,a,b,Jin)` stratum | 0.532 | 0.534 | 0.549 | 0.0153 `[0.0095, 0.0249]` |

The final comparison uses 6,790 pairs across 18 cells; its pair-weighted residual `Jin` difference is `0.00630`. The aperture advantage is small but independently reproducible.

`Jin-Ja` is unresolved for initiation (`0.0207`, CI `[-0.0085, 0.0424]`) and transfer (`-0.0127`, CI `[-0.0348, 0.0057]`) before conditioning, and remains unresolved after both conditioning stages. Therefore simply deleting negative current is not the supported explanation.

These AUCs quantify rank information, not a fitted transition law. Training-only supply quartiles are used solely to form comparable held-out pairs, not as a new collective coordinate or predictor.

## Backward first-exit formulation

### State and microscopic inputs

Use local source coordinates with the physical opening centered at the divider. For one fish during one approach attempt,

\[
D=(-\ell,0)\times(-h,h),\qquad h=\ell/2,
\]

and let `phi` be heading relative to the positive local longitudinal direction. Freeze or prescribe the instantaneous microscopic turning rate `Omega` and retain the production angular noise `I_n=0.5`:

\[
dZ_t=\cos\Phi_t\,dt,
\qquad dY_t=\sin\Phi_t\,dt,
\qquad d\Phi_t=\Omega(Z_t,Y_t,\Phi_t)\,dt+I_n\,dW_t.
\]

For a fish-specific frozen-torque approximation, `Omega` is the instantaneous production value

\[
\Omega_i^0=\operatorname{clip}
(\omega_{{\rm wall},i}+\omega_{{\rm social},i}),
\]

held constant over the short first-exit calculation. No coefficient is introduced.

Define the first exit

\[
\tau=\inf\{t>0:(Z_t,Y_t)\notin D\}.
\]

Capture occurs when the exit is through the physical aperture `z=0`, rather than through the lateral approach boundaries or the upstream boundary `z=-ell`. The lateral and upstream exits represent failure of this approach attempt, not removal of the fish from the microscopic system; later re-entry is a separate renewal event.

### Backward generator and boundary data

The backward generator is

\[
\mathcal L=cos\phi\,\partial_z+sin\phi\,\partial_y
+\Omega(z,y,\phi)\,\partial_\phi
+{I_n^2\over2}\,\partial_{\phi\phi}.
\]

The aperture committor

\[
q(z,y,\phi)=P_{z,y,\phi}\{\text{exit through the opening}\}
\]

solves the hypoelliptic backward problem

\[
\mathcal L q=0\quad\text{in }D\times S^1,
\]

with periodicity in `phi` and kinetic boundary conditions only on outgoing phase-space boundaries:

\[
q=1\quad\text{at }z=0,\ \cos\phi>0,
\]

\[
q=0\quad\text{at }z=-\ell,\ \cos\phi<0,
\]

\[
q=0\quad\text{at }y=h,\ \sin\phi>0,
\qquad
q=0\quad\text{at }y=-h,\ \sin\phi<0.
\]

Incoming orientations require no independently imposed boundary values; they are determined by the interior kinetic solution.

The Laplace transform of the capture-time distribution,

\[
G_\lambda=E[e^{-\lambda\tau}I_{\rm capture}],
\]

satisfies

\[
(\mathcal L-\lambda)G_\lambda=0
\]

with the same aperture/failure values. Hence `q=G_0`. The capture-weighted first moment

\[
m=E[\tau I_{\rm capture}]
\]

satisfies

\[
\mathcal L m=-q,
\]

with zero outgoing boundary data, and the conditional mean capture time is `m/q` where `q>0`.

A persistence-aware, parameter-free supply functional is then

\[
J_{\rm FP}={1\over N}\sum_{i\in A}\max(\cos\phi_i,0)\,
q(z_i,y_i,\phi_i;\Omega_i^0,I_n).
\]

This is the theoretical extension logically supported by the empirical geometry test.

## Rigorous deterministic limit

Set `I_n=0` and `Omega=0`. The backward characteristics are straight:

\[
Z_t=z+t\cos\phi,\qquad Y_t=y+t\sin\phi.
\]

For `u=cos(phi)>0`, the divider is reached at

\[
t_g=-z/u,
\]

and the transverse intercept is

\[
y_g=y-z\tan\phi.
\]

Since the approach rectangle is convex, the characteristic reaches the opening before either lateral boundary exactly when

\[
|y-z\tan\phi|\le h.
\]

For `u<=0`, the straight characteristic cannot exit through the forward opening before an upstream or lateral exit. Therefore

\[
q_0(z,y,\phi)=I(u>0)I(|y-z\tan\phi|\le h),
\]

and substituting `q_0` into `J_FP` gives exactly `Jcap`. The deterministic straight-line capture condition is thus the zero-noise, zero-turning characteristic solution of the backward first-exit problem—not an independently fitted rule.

For nonzero prescribed constant `Omega`, the zero-noise characteristics are circular arcs,

\[
\Phi_t=\phi+\Omega t,
\]

\[
Z_t=z+{\sin(\phi+\Omega t)-\sin\phi\over\Omega},
\qquad
Y_t=y+{\cos\phi-\cos(\phi+\Omega t)\over\Omega},
\]

and first contact with the same physical boundaries supplies the deterministic finite-turn correction.

## Why this is not an exact collective closure

The wall torque is determined by the fish position and fixed boundary, but the social torque is not a prescribed one-particle field. It depends on every other fish position and heading, the instantaneous Delaunay neighborhood, wall-blocked visibility, and visual normalization. Those quantities evolve while the candidate fish approaches the opening and can change discontinuously when neighbor topology changes.

Consequently `(z,y,phi)` with frozen `Omega_i^0` is not Markovian under the full school dynamics. The exact backward equation would live in the full `3N`-dimensional microscopic state and its capture set would couple individual exits to the evolving school. The frozen-torque first-exit problem is a short-time conditional approximation, appropriate for the empirically supported aperture gate; it cannot by itself predict collective continuation, repeated lateral renewal, or whole-school relocation.

The next theoretical step is therefore to solve or asymptotically approximate this parameter-free single-attempt problem and test it against archived short-time capture outcomes. Any longer-time collective theory would additionally require a source-grounded description of how the many-fish torque environment evolves, not a fitted potential.

## Reproducibility

- `analyze_aperture_geometry_independence.py` recomputes `Jin`, verifies exact `Jcap` parity, constructs outcome-blind supply strata, and performs seed-cluster held-out comparisons.
- `geometry_independence.png` is the compact empirical comparison.
- `heldout_auc.csv`, `paired_auc_contrasts.csv`, and `training_supply_strata.csv` contain the validation results.
- `summary.json` records all estimates and pair counts.

No microscopic simulation, simulator modification, fitted predictor, potential, or new empirical coordinate was introduced.
