# Exact collective evolution of source-school gap-directed polarization

## Physical conclusion

Collective gap-directed motion is sustained primarily by the production
model's **visual social turning**, while motion of the school center of mass
rotates the gap reference direction and is the dominant term destroying the
measured polarization.  Continued-transfer events retain a much larger
positive social contribution and suffer less destructive reference rotation
before first crossing than failed attempts.  Wall turning does not maintain
the late precursor, angular diffusion relaxes rather than creates its positive
polarization, and production torque clipping is negligible in this regime.

The net gap-directed polarization still declines in both outcome groups.  The
mechanism is therefore persistence by partial compensation, not growth of a
newly forming order parameter: social turning offsets part of the stronger
COM-geometry loss for long enough that the source keeps supplying the opening.
This is an exact source decomposition combined with a retrospective outcome
association.  It does not by itself prove that changing social torque would
cause a selected event to continue.

## Production dynamics and definitions

Fix the original source chamber for each event.  At time `t`, let `S(t)` be the
fish currently in that physical chamber, `N_s=|S|`, and

```text
R_s = (1/N_s) sum_(i in S) r_i,
g_s = -R_s/|R_s|,
P_s = (1/N_s) sum_(i in S) e_i,
a   = P_s . g_s,
e_i = (cos(theta_i), sin(theta_i)).
```

The archived production update has unit free-swimming speed, collision
overrides that can set the stored translational velocity to zero, and

```text
d theta_i = Omega_i dt + In dW_i,        In = 0.5,
D_theta = In^2/2 = 0.125.
```

The actual deterministic torque is not an unsaturated sum.  If
`omega_w^0` is the raw forward-ray wall torque and `omega_soc` is the visual
social torque, the source applies

```text
c_w   = clip(omega_w^0, +/-pi/dt) - omega_w^0,
omega_w = omega_w^0 + c_w,
c_tot = clip(omega_w + omega_soc, +/-pi/dt)
        - (omega_w + omega_soc),
Omega = omega_w^0 + omega_soc + c_w + c_tot.
```

Both clipping stages are retained below.  The social rule itself is the
production visual-weighted combination of heading alignment and relative
position/bearing turning; “social” should not be read as pure polar alignment.

## Exact crossing-free stochastic equation

Between chamber crossings, membership is fixed.  Ito's formula gives

```text
d e_i = e_i_perp Omega_i dt - D_theta e_i dt
        + In e_i_perp dW_i.
```

The source COM moves with the actual stored source velocity
`V_s=(1/N_s) sum v_i`, including collision overrides.  Therefore

```text
d g_s/dt = -(I - g_s g_s^T) V_s / |R_s|.
```

Combining these identities gives the exact crossing-free semimartingale

```text
da = [ (1/N_s) sum_(i in S) (g_s . e_i_perp) Omega_i
       - D_theta a
       - P_s . (I-g_s g_s^T) V_s/|R_s| ] dt
     + (In/N_s) sum_(i in S) (g_s . e_i_perp) dW_i.
```

Substitution of the actual production torque separates the angular drift into
raw wall, social, internal wall clipping, and final total clipping terms.  The
third deterministic term is rotation of the source-COM-to-gap direction; it
is geometry, not a torque.  The continuous quadratic variation is

```text
d[a]_t = In^2/N_s^2 sum_(i in S) (g_s . e_i_perp)^2 dt.
```

Thus even the noise amplitude is not determined by `a` alone.

At the native finite timestep, angular diffusion was split exactly into its
conditional Ito relaxation and observed realization:

```text
noise relaxation = [exp(-In^2 dt/2)-1] times deterministic projection,
realized noise    = observed projection minus its conditional mean.
```

The torque increments use the exact secant slope of cosine, so the discrete
budget is algebraic rather than an Euler approximation to `da`.

## Membership changes at chamber crossings

At a crossing, the original chamber remains the source but `S`, `N_s`, `P_s`,
`R_s`, and `g_s` jump.  With positions and headings evaluated at the crossing,

```text
Delta a_jump = (P_s^+ - P_s^-) . g_s^-
               + P_s^+ . (g_s^+ - g_s^-).
```

The first term is the heading-moment and normalization change caused by the
entering/leaving fish.  The second is the instantaneous COM-reference change.
In the finite-step audit, fixed-membership COM motion is applied before these
two membership terms; all three pieces are separately exact.  The complete
evolution is the crossing-free stochastic equation plus the sum of these
jumps.

## Source-grounded evaluation

The fixed cohort contains 337 majority-to-minority first-crossing episodes
(185 continued by at least four additional net fish over five time units and
152 failed that threshold) from 37 independent seeds.  The outcome definition
describes short propagation, not complete school relocation.

Production torques were replayed at 10,447 prospective native steps from
`t=-2` to `+1` at spacing 0.1.  All 1,929 chamber-crossing steps in the same
windows were additionally replayed for the jump balance.  Results were
standardized to the pooled initial-imbalance strata and uncertainty resampled
the 37 original seeds.  No fish or event was treated as an independent
bootstrap cluster.

The maximum absolute finite-step identity residual was `1.38e-15`.  Across all
replayed steps, only one fish-step activated the internal wall cap and one
activated the final total-torque cap.  The period-mean clipping contribution
was consequently zero at displayed precision.

### Polarization retained before crossing

At `t=-2`, the q-standardized source polarization was 0.612 in continued
events and 0.460 in failed events.  Immediately after the defining crossing it
was 0.357 versus -0.079.  Both groups lose polarization, but failed events lose
it much faster.  In the final pre-crossing unit, the exact smooth source terms
in polarization per unit time were:

| term | continued | failed | continued - failed, 95% seed CI |
|---|---:|---:|---:|
| raw social turning | +0.1524 | +0.0584 | +0.0939 `[0.0721, 0.1147]` |
| wall turning | -0.0123 | -0.0239 | +0.0117 `[-0.0001, 0.0233]` |
| torque clipping | 0 | 0 | 0 |
| Ito noise relaxation | -0.0598 | -0.0138 | -0.0460 `[-0.0542, -0.0372]` |
| fixed-membership COM geometry | -0.2687 | -0.3371 | +0.0684 `[0.0462, 0.0919]` |
| smooth total, excluding realized noise | -0.1884 | -0.3164 | +0.1280 `[0.0946, 0.1635]` |

Social turning is the only large positive angular drift in both groups.  Its
continued-minus-failed separation is larger than the geometry separation in
the final pre-crossing unit.  COM-reference rotation is nevertheless the
largest absolute term and makes the net drift negative.  Noise relaxation is
more negative in continued events because their positive `a` is larger; this
is the exact `-D_theta a` effect, not evidence that noise selectively targets
successful schools.  The realized-noise contrast is unresolved
(`-0.0378`, CI `[-0.1482, 0.0803]`).

Wall turning contributes a small positive separation earlier (`+0.0202`, CI
`[0.0076, 0.0319]` over `t=-2` to `-1`) but becomes negative in both groups in
the last unit, with an unresolved contrast.  It is not the mechanism sustaining
the late gap-directed precursor.

### Immediately after first crossing

Social turning remains strongly separated: +0.0970 in continued events versus
+0.0042 in failures, difference +0.0928 `[0.0717, 0.1149]`.  However, continued
events then experience more destructive COM-reference rotation and stronger
diffusive relaxation, so the net smooth-drift difference is unresolved.  Their
advantage persists because they enter the post-crossing interval with much
higher `a`, rather than because their total post-crossing drift becomes
positive.

Crossing membership also destroys polarization during continuation.  In the
first post-crossing unit, the heading/normalization jump is -0.0451 per unit
time in continued events versus -0.0085 in failures.  The associated COM jump
partly compensates it (+0.0269 versus +0.0042), leaving net membership effects
of -0.0182 and -0.0043.  This is why membership changes cannot be represented
by simply continuing the crossing-free equation with a varying `N_s`.

## Connection to exact approach and population conservation

Let `W` be the established downstream-weighted source approach storage,
`J_a` the archived stored-velocity approach current, `L_W` weighted lateral
exchange, and `F_neck` positive net transfer out of the source.  The previously
validated identity is

```text
F_neck = N J_a/ell - (1/ell) dW/dt - L_W + epsilon_path,
dN_s/dt = -F_neck                         (between discrete notation changes).
```

Equivalently, for `b=W/(N ell)`, population release is constrained by directed
supply, lateral exchange, and neck flux.  The polarization equation supplies
the missing angular mechanism behind `J_a`: positive social torque maintains
source headings with a gap-directed component, while COM-reference rotation
erodes that organization.

The two laws do not collapse into one scalar equation.  `J_a` samples only
fish in the approach rectangle and projects on the fixed chamber axis, whereas
`a` averages the whole current source population and projects on a moving
COM-to-gap direction.  Moreover, `F_neck` determines only the number crossing.
The polarization jump also needs each crossing fish's heading and its effect
on the source COM.

Numerically, the linked transport budget gives continued-minus-failed net neck
flux `+6.550` fish/time in the first post-crossing unit.  The identical source
count-loss contrast from the membership audit is `+6.550` (seed CI
`[5.680, 7.345]`).  This directly connects the two exact accounts.  The prior
transport identity attributes that flux separation to maintained longitudinal
supply plus reduced lateral loss; the present identity shows that persistent
positive social turning is the principal microscopic angular term associated
with that maintained organization.

## Why an analytical collective closure is still impossible

The derivation identifies the missing information precisely:

- The social drift contains the correlation
  `mean[(g_s.e_i_perp) omega_soc,i]`.  It depends on the full visible Delaunay
  neighbor graph, relative headings, distances, and bearings, not on the means
  `q`, `a`, or `b`.
- The wall drift requires each fish's forward-ray intersection distance and
  side, and clipping depends on the distribution of the summed torques.
- COM-reference rotation requires `|R_s|` and the transverse component of the
  actual source velocity, including collision overrides.
- Rotational-noise variance requires a second angular moment, even though its
  conditional mean contribution is exactly `-D_theta a`.
- Membership jumps require the headings and positions of the particular fish
  crossing in each direction.  Scalar neck flux supplies only their signed
  count.

Therefore no exact drift or diffusion for `a`, or for `(q,a,b)`, follows from
the current variables.  Closing the theory requires deriving conditional
microscopic moment and boundary-flux laws from the production interactions;
fitting an empirical drift or potential would hide rather than solve this
problem.

## Interpretation boundary and reproducibility

The equation, saturation split, finite-step identities, and connection between
source count loss and neck flux are exact.  Successful-versus-failed term
differences are retrospective associations standardized by initial imbalance.
They establish which terms coexist with sustained transfer, not an
interventional causal effect.  In particular, the aggregate social term does
not distinguish within-source alignment from cross-divider feedback.

`analyze_collective_polarization_evolution.py` performs the source replay,
identity audit, seed bootstrap, and figure generation.  The unchanged
production source SHA-256 is
`5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`.
No new microscopic simulation, fitted drift, effective potential, or new
collective coordinate was introduced.
