# Targeted validation of an oscillatory shape--heading mechanism

## Verdict

**The proposed weakly unstable, traveling two-lobed feedback is not supported as a persistent mechanism in the full model at this parameter point.** A COM-relative `m=2` deformation is present, but it is neither specific to milling nor characterized by a reproducible period near 11.4. In the two sustained mills, the source-exact visual-alignment projection removes `m=2` heading-mode power in every one of six time blocks. Direct hydrodynamic turning is much more variable and can transiently inject power, including shortly before the only observed transition, but it does not produce a common growing oscillation.

The bounded result is therefore negative for the ideal-ring instability, not for all shape--heading coupling. The full dynamics contain instantaneous coupling and an event-specific hydrodynamic impulse, but the data do not establish it as the cause of schooling--milling bistability or a stability boundary.

## Scope and provenance

This read-only analysis used the original `state` and `rdot` arrays for seeds 0, 7, 8, and 9. No simulator, safety condition, MAT file, interaction coefficient, or trajectory was changed. Stage A streamed frames every 0.05 time units. Stage B replayed 2,658 source-exact snapshots every 0.2 time units in three 60-unit blocks from each sustained mill (seeds 7 and 8), the seed-9 interval 250--360, and one schooling-control block from seed 0.

The replayed production source is `BEM_circular_parallel_instance.py`, SHA-256 `2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`. The source-exact implementation was reused from `theory/circle_collective_first_stage/com_relative_turning/analyze_com_relative_turning.py`; the new driver is `analyze_modes_and_feedback.py`.

## Exact definitions and identities

With school center `C`, relative position `y_i=x_i-C`, radius `rho_i`, angle `phi_i`, and mean radius `rho_bar`, the measured modes are

```text
Q_m = mean[(rho_i-rho_bar) exp(-i m phi_i)],
A_m = mean[exp(-i m phi_i)].
```

For circulation sign `s=sign(mean(sin(theta_i-phi_i)))`, the wrapped tangent-heading deviation and its mode are

```text
v_i = wrap(theta_i-phi_i-s*pi/2),
V_m = mean[v_i exp(-i m phi_i)].
```

The exact kinematic shape derivative used here is

```text
dot Q_m = mean[((u_r,i-mean(u_r))-i m (rho_i-rho_bar) omega_i^C)
               exp(-i m phi_i)].
```

The deterministic, moving-basis `m=2` heading derivative is

```text
dot V_2 = mean[((Omega_i-omega_i^C)-2 i v_i omega_i^C) exp(-2 i phi_i)].
```

`Omega_i` is split exactly into alignment, attraction, wall, direct-hydrodynamic, BEM-mediated, and clipping-residual terms. The geometric contribution contains both `-omega_i^C` and the Fourier-basis term. Modal power from source `k` is the physically interpretable instantaneous quantity

```text
G_k = Re(conj(V_2) D_k),
```

where `D_k` is that source's contribution to `dot V_2`. Positive `G_k` amplifies the heading-deformation amplitude; negative `G_k` damps it. This is an exact projection, not a fitted energy or potential.

To remove rigid orbital rotation, `psi(t)=integral mean_i(omega_i^C) dt` and `Z_m^co=Z_m exp(+i m psi)` for `Z=Q,A,V`. Thus a laboratory-frame phase rotation is not counted as an internal oscillation.

## Numerical validation

- Stored velocity replay: maximum error `3.55e-15`.
- Alignment-plus-attraction reconstruction: maximum error `2.23e-12`.
- Branch-safe native-step `V_2` identity: maximum residual `1.07e-4`.
- Wrapped-coordinate branch crossings: 60 fish-snapshots, explicitly flagged rather than used for that derivative audit.
- COM cutoff: `rho_i/R <= 0.005`. Stage A excluded 2,063 of 7,839,600 fish-frames (0.0263%); Stage B excluded 73 of 265,800 (0.0275%). Modes are averaged over the remaining fish in those rare frames.
- Halving the Stage-A interval from 0.05 to 0.025 changed mean `|Q_m|` by at most 0.0306%; common-frame values agreed exactly.

The small nonzero native-step residual is consistent with finite-step polar-coordinate geometry. It is not a mismatch of the replayed simulator velocity or social torque.

## Measurements in sustained mills

The simplified ring predicts `r*=2.52318`, `Re(lambda_2)=0.00701`, `Im(lambda_2)=0.55169`, and period 11.4. The actual sustained mills are substantially different base states: mean COM-relative radii are 1.454 and 1.464 (about 42% below `r*`), and their mean orbital periods are 7.35 and 7.59. Orbital period is not itself the internal-mode period, but the discrepancy shows why the ideal-ring eigenvalue cannot be transferred directly.

Across the six prespecified mill blocks, the observed `m=2` spectral peak periods were 20.02, 7.51, 10.01, 12.01, 6.67, and 10.01. Power at the reference frequency ranged from 1.5% to 100% of each block's peak. Complex `Q_2`, `V_2`, and torque phases likewise showed no repeatable phase lag across blocks. There is therefore no single resolved full-model oscillation at the predicted frequency.

The deformation is modest and not milling-specific. Mean `|Q_2|/rho_bar` is 0.068 and 0.062 over seeds 7 and 8, while the schooling control is 0.150. Seed 9 changes from 0.154 before the transition to 0.078 in the mill. An `m=2` signal alone cannot be interpreted as the milling instability.

Source-resolved heading-mode power, treating the six time blocks as uncertainty units:

| contribution | mean `G_k` | 95% t interval | signs by block |
|---|---:|---:|---:|
| alignment | -0.001210 | [-0.001756, -0.000663] | 0 positive / 6 negative |
| attraction | -0.000195 | [-0.000313, -0.000077] | 0 / 6 |
| wall | -0.000394 | [-0.000563, -0.000224] | 0 / 6 |
| direct hydrodynamics | -0.002277 | [-0.007243, 0.002689] | 2 / 4 |
| BEM hydrodynamics | -0.00000057 | [-0.00000098, -0.00000015] | 0 / 6 |
| COM geometry | +0.000715 | [-0.000071, 0.001501] | 5 / 1 |
| total | -0.000794 | [-0.006372, 0.004784] | 3 / 3 |

The corresponding mean exact `Q_2` amplitude derivative is `-0.000429`, with 95% interval `[-0.002036, 0.001178]`. Thus neither the deformation nor total heading-mode power shows reproducible positive growth. Most importantly for the proposed nonreciprocal visual mechanism, alignment is consistently damping in these full-model mills.

## Seed 9 transition chronology

The dashed line in the timeline marks the archived transition time 312.41. Before the polarization collapse, `M_C` begins changing sign and there is a short interval (306--311.5) in which direct hydrodynamic power is positive (`+0.0233`), total `V_2` power is positive (`+0.0102`), and the exact `Q_2` amplitude derivative is positive (`+0.0228`). This is direct evidence for a transient hydrodynamically supported shape--heading episode before the collapse.

It is not evidence for the proposed growing eigenmode:

1. Mean `|Q_2|` is already lower in 306--311.5 (0.172) than during 250--300 (0.204) and 300--306 (0.220).
2. During the actual collapse interval 311.5--315, alignment, direct hydrodynamics, and total heading power are all negative.
3. `Q_2` does not grow monotonically or reach an amplitude absent from the long schooling interval.
4. Comparable hydrodynamic power fluctuations occur in sustained states, and there is only one transition event.

The defensible temporal statement is therefore: a direct-hydrodynamic `m=2` impulse accompanies the loss of circulation and precedes the sharp polarization collapse, but a coherent two-lobed disturbance does not demonstrably initiate or destabilize the state.

## Physical interpretation and relation to prior work

[Huang, Ling & Kanso (2024)](https://doi.org/10.1073/pnas.2406293121) established confined-school bistability and inferred stochastic reduced dynamics. [Hang et al. (2026)](https://doi.org/10.1038/s41467-026-70569-y) derived directional information transport from nonreciprocal visual interactions. The present test asks a narrower bridge question: whether that visual transport feeds a traveling `m=2` instability in the original full model.

The answer at this parameter point is **no**. Visual alignment does participate in the exact angular dynamics, but its `m=2` projection is restoring on average. COM geometry tends to return some power, while direct hydrodynamics supplies large, state-dependent fluctuations and one pre-collapse impulse. This supports a noisy, dynamically reorganizing full school rather than the particular weakly unstable two-neighbor-ring mode.

This does not refute the ideal-ring calculation: its base radius, fixed neighbor graph, absence of hydrodynamics, and uniform spacing do not match the observed mill. Nor does one negative parameter-point test establish the full stability boundary.

## Smallest meaningful next calculation

If the stability mechanism must be pursued, the smallest informative calculation is **not** a parameter sweep or an effective-potential fit. Numerically differentiate the exact deterministic simulator right-hand side around a few archived stable-mill snapshots using prescribed `m=1,2,3,4` radial and heading perturbations, while retaining each snapshot's actual Delaunay graph and BEM response. The resulting small projected tangent operator would directly show which full-model mode is least damped and whether hydrodynamics changes the sign of the `m=2` eigenvalue. A handful of snapshots from both mill chiralities is enough for this bounded check; no new long trajectories are required.

## Reproduction and outputs

Run from the repository root:

```bash
MPLCONFIGDIR=/tmp/mplconfig uv run --with numpy --with scipy --with matplotlib \
  python theory/circle_collective_first_stage/oscillatory_shape_heading/analyze_modes_and_feedback.py \
  --data-dir /path/to/BEM_data
```

Decisive outputs are:

- `outputs/shape_modes_stable_and_transition.png`
- `outputs/corotating_Q2_V2.png`
- `outputs/stable_mill_modal_power.png`
- `outputs/seed9_mode_feedback_timeline.png`
- `outputs/stable_mill_feedback_summary.csv`
- `outputs/seed9_feedback_windows.csv`
- `outputs/mode_spectral_summary.csv`
- `outputs/cross_phase_summary.csv`
- `outputs/sampling_convergence.csv`
- `outputs/validation_audit.json`

The compressed NPZ files contain only the streamed reduced modes, not raw trajectories.
