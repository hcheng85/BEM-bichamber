# Branch-free heading mode and exact power budget

> **Follow-up qualification:** The [connected-circulation audit](../connected_circulation/REPORT.md) resamples the same `t=100–160` mill blocks much more densely and tests native-step conditional balance. It shows substantial grid sensitivity in the mean expected power. The historical tables below are preserved, but “damped mode” must not be read as a measured restoring eigenvalue: negative mean deterministic power can balance stochastic injection, and does not establish linear stability. Use the follow-up for the corrected stationarity interpretation.

## Scientific verdict

The branch-free result is decisive for the proposed modal mechanism:

1. **In stable mills, `H2` is a damped but noise-excited mode.** Its deterministic drift power is negative in both independent milling seeds. Finite-`N` angular noise injects power and maintains fluctuations.
2. **The mode is not milling-specific.** Mean `|H2|` is larger in the two schooling controls (0.117 and 0.110) than in the two sustained mills (0.068 and 0.078). The `m=1` control is also larger than `m=2` in the mills.
3. **No growing `H2` precursor survives.** In seed 9, mean `|H2|` falls from 0.110 during `t=250--300` to 0.087 during `t=306--312.4`, while deterministic `H2` power is strongly negative. It falls again during the polarization collapse.
4. **The earlier wrapped `V2` fluctuations were not entirely artificial, but their torque interpretation was coordinate-dependent.** In mills, `|H2|` and wrapped `|V2|` are strongly correlated (about 0.85), yet the source signs differ: alignment damps the wrapped coordinate but supplies branch-free `H2` power. `H2` is the defensible observable.

Therefore the traveling two-lobed heading instability should not be pursued as the initiating schooling--milling mechanism for these data. One observed transition cannot establish causality, but it directly contradicts the required leading positive `H2` production.

## Exact observable

For COM-relative polar angle `phi_i` and `alpha_i=theta_i-phi_i`, define

```text
A_m = mean(exp(-i m phi_i))
B_0 = mean(exp(i alpha_i))
B_m = mean(exp(i alpha_i-i m phi_i))
H_m = B_m-B_0 A_m
```

Equivalently, with `u_i=exp(i alpha_i)` and `f_i=exp(-i m phi_i)`,

```text
H_m = mean[u_i (f_i-A_m)].
```

This is smooth in heading and removes the angular-density contribution. For a perfectly tangent mill, `u_i=+i` or `-i` for every fish, so `B_m=B_0 A_m` for any nonuniform angular density. A deliberately nonuniform synthetic configuration gave exactly zero `H1` and `H2` to floating-point precision.

No fish was exactly at the COM in the replayed snapshots; the minimum radius was `0.000543 R`. Unlike the prior wrapped analysis, no branch exclusion or regularization is required.

## Production update order

The relevant production lines are in `BEM_circular_parallel_instance.py`:

- `GetNoise` creates only an orientation increment, `In*sqrt(dt)*Normal(0,1)`; its two position rows are zero.
- Velocity and angular torque use `state[n]`, including the old heading.
- Total angular rate is clipped to `[-pi/dt,+pi/dt]` and collision/extreme-velocity overrides are applied.
- `state[n+1]` is formed with deterministic position increment `rDot*dt`, deterministic angular increment, and the new angular noise.
- The wall-safety retry may set selected position velocities to zero, but reuses the same deterministic angular increment and angular-noise draw.
- Stored `rdot[:,:,n+1]` advances `state[:,:,n]`.

The noise random numbers are literally drawn before the final update, so it would be inaccurate to say that the draw occurs after positions are stored. However, the draw never enters velocity, collision detection, or the positional retry, and its positional components are identically zero. Consequently `x[n+1]`, its COM, `phi_i'`, and `A_m'` are conditionally deterministic given the current state. This is exactly the condition needed below.

## Exact one-step Gaussian identities

Let

```text
theta_i* = theta_i[n] + Omega_i,eff dt
c_i = exp(i(theta_i*-phi_i')) (exp(-i m phi_i')-A_m').
```

Independence of the new Gaussian angular increments gives

```text
E[H_m[n+1] | F_n] = exp(-In^2 dt/2) mean(c_i)

Var(H_m[n+1] | F_n)
 = (1-exp(-In^2 dt)) (1-|A_m'|^2)/N,
```

where complex variance means `E|H-EH|^2`.

These identities were tested at six fixed snapshots spanning both mill chiralities, seed-9 schooling, pre-collapse, collapse, and early milling. For each snapshot and `m=1,2`, 20,000 independent angular-noise redraws used the production `In=0.5` and `dt=0.001`:

- largest complex-mean discrepancy: 1.42 Monte Carlo standard errors;
- largest relative variance discrepancy: 0.694%;
- no fitted coefficient or altered noise convention was used.

The numerical validation is in `outputs/conditional_redraw_validation.csv`.

## Exact continuous-time budget

For angular source `k`,

```text
D_k = i mean[u_i (f_i-A_m) Omega_i^(k)].
```

The moving COM-polar basis contributes

```text
D_geometry = -i(m+1) mean(u_i f_i omega_i^C)
             +i A_m mean(u_i omega_i^C)
             +i m B_0 mean(omega_i^C f_i).
```

Angular diffusion contributes `D_diffusion=-(In^2/2)H_m`. Therefore

```text
d|H_m|^2/dt drift = sum_k 2 Re(conj(H_m) D_k)
                    +2 Re(conj(H_m) D_geometry)
                    -In^2 |H_m|^2
                    +In^2 (1-|A_m|^2)/N.
```

The last term is finite-`N` stochastic injection, not deterministic steering. Alignment, attraction, wall, direct hydrodynamics, BEM hydrodynamics, and clipping are retained separately.

### Finite-step closure

At every one of 2,658 source-replayed snapshots, the exact conditional next-step power was compared with the first-order continuous-time budget. Actual power change is decomposed algebraically as

```text
actual change = dt * continuous expected power
              + finite-step discretization residual
              + realized zero-mean noise residual.
```

- algebraic closure error: `2.71e-20` maximum;
- median absolute discretization residual: `5.70e-7`;
- 95th percentile: `4.12e-6`;
- mean absolute residual was 5.49% of mean absolute exact conditional change;
- realized-noise residual mean: `(6.37 +/- 4.11)e-6` (mean ± ordinary SE), consistent with zero.

The maximum discretization residual, `0.00111`, occurs in a highly nonlinear angular update. The exact conditional formula, not the first-order expansion, is the reference at such steps.

Velocity replay agreed with stored `rdot` to `3.55e-15`; the alignment-plus-attraction split agreed with the production visual torque to `2.23e-12`.

## Stable-state source budget

The uncertainty units are independent seeds, not the three time blocks within a seed. With only two milling seeds, the table reports the range of the two seed means rather than a misleading narrow confidence interval.

| `H2` power source | stable-mill mean | seed-mean range | interpretation |
|---|---:|---:|---|
| alignment | +0.003692 | [+0.003476, +0.003907] | supplies power in both seeds |
| attraction | -0.000313 | [-0.000524, -0.000101] | weak damping |
| wall | -0.000206 | [-0.000221, -0.000192] | weak damping |
| direct hydrodynamics | -0.002512 | [-0.005646, +0.000621] | large and seed-dependent |
| BEM hydrodynamics | -0.00000048 | [-0.00000075, -0.00000020] | negligible damping |
| COM geometry | -0.003602 | [-0.003651, -0.003553] | consistent damping |
| angular diffusion | -0.001786 | [-0.002080, -0.001491] | consistent damping |
| **deterministic drift total** | **-0.004745** | **[-0.007633, -0.001857]** | damped in both seeds |
| finite-`N` injection | +0.002428 | [+0.002412, +0.002445] | noise excitation |
| expected total | -0.002317 | [-0.005221, +0.000588] | seed-dependent after noise |

Thus visual alignment does inject branch-free `H2` power, consistent with a real shape--heading coupling. It is overcome by COM geometric transport, angular diffusion, and—in seed 7—direct hydrodynamics. The resulting deterministic mode is damped, not unstable.

Schooling controls behave differently: alignment power is even larger (`+0.01040` mean), and `|H2|` itself is larger. Positive alignment production is therefore not a milling-instability signature.

## Angular clipping

Clipping occurred for only four fish in two of the 2,658 snapshots, both in seed 8:

- `t=112.4`: two fish; maximum clip-residual rate `6994.2`; exact discrete clip contribution to next-step `|H2|^2` was `+2.82e-4`.
- `t=522.8`: two fish; maximum clip-residual rate `1243.8`; exact discrete contribution was `+6.07e-4`.

These events explain why clipping looked substantial in block-averaged wrapped-coordinate budgets. At the individual steps, the continuous `Omega_clip` projection is not quantitatively reliable because the turn is order one over a single timestep. The report therefore uses the exact counterfactual difference between clipped and unclipped conditional next-step power. Clipping is locally important but rare, occurs neither during seed 9's transition nor persistently in the mills, and cannot explain the transition.

## Seed 9: chronology and matched control

The transition is at approximately `t=312.41`.

| interval | mean `P` | mean `M_C` | mean `|H2|` | deterministic `H2` power | expected power incl. injection |
|---|---:|---:|---:|---:|---:|
| schooling 250--300 | 0.900 | +0.216 | 0.110 | +0.00430 | +0.00645 |
| late school 300--306 | 0.914 | +0.204 | 0.121 | +0.01062 | +0.01288 |
| pre-collapse 306--312.4 | 0.908 | +0.086 | 0.087 | **-0.01496** | **-0.01272** |
| collapse 312.4--315 | 0.497 | -0.694 | 0.061 | -0.00834 | -0.00618 |
| early mill 315--330 | 0.273 | -0.841 | 0.097 | -0.00096 | +0.00132 |
| mill 330--360 | 0.200 | -0.859 | 0.083 | -0.01186 | -0.00943 |

During the critical pre-collapse interval, alignment power is `-0.00168`, direct-hydrodynamic power is `-0.01024`, and geometry is also negative. There is no positive deterministic precursor hidden by the wrapped coordinate.

A prospective control was selected from nonoverlapping earlier 6.4-unit seed-9 windows using only initial `P` and `M_C`, with the same 0.2 replay interval. The closest control was `t=282.0--288.4` (standardized initial-state distance 0.59):

- control `|H2|=0.0965`, deterministic power `+0.0190`;
- pre-collapse `|H2|=0.0870`, deterministic power `-0.0150`.

This is descriptive within one seed, not an independent replicate. Nevertheless, its direction is opposite to a growing modal trigger. Neither `Q2` nor wrapped `V2` shows a clean leading growth sequence, and `H1` remains larger than `H2` through much of the trajectory.

## Answers to the requested decisions

**A.** `H2` is a damped but noise-excited mode in the two sustained mills. It is not a growing instability and is not unique to milling.

**B.** Alignment supplies `H2` power, while COM geometric transport and angular diffusion consistently remove it. Attraction and wall turning weakly remove power. Direct hydrodynamics is large but changes sign between seeds; BEM turning is negligible. Rare clipped turns have measurable exact one-step effects but no sustained or transitional role.

**C.** `H2` does not causally precede the observed transition. Its amplitude and deterministic production decrease before polarization collapses. With one event, causality cannot be established in general; the available event specifically rejects the required positive precursor.

**D.** Stop pursuing the modal-instability hypothesis for this dataset. The smallest next first-principles test is a bounded joint drift audit of existing snapshots in `(P^2,M_C)` conditioned on two COM-geometry quantities: COM displacement from the tank center and the angle between polarization and the COM radial direction. The exact source decompositions for `P^2` and `M_C` already exist. Testing whether their restoring signs change before seed 9's transition, against the same prospective control, would distinguish loss of translational-school geometry from angular-circulation capture without fitting a surrogate or running new simulations.

## Files

- `analyze_branch_free_mode.py`: complete streaming, replay, redraw, budget, and plotting pipeline.
- `outputs/branch_free_snapshot_budget.csv`: source-exact snapshot budget.
- `outputs/interval_source_budget.csv`: disjoint block summaries.
- `outputs/state_source_budget.csv`: seed-level state comparison.
- `outputs/seed9_transition_windows.csv`: transition chronology and matched control.
- `outputs/conditional_redraw_validation.csv`: independent Gaussian-redraw audit.
- `outputs/validation_audit.json`: numerical closure and provenance.
- `outputs/seed9_branch_free_timeline.png`, `stable_mill_H2_source_budget.png`, and `conditional_identity_validation.png`: targeted figures.

The compact NPZ files contain only streamed reduced observables; raw trajectories remain outside the repository and unchanged.
