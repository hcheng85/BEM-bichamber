# Bounded same-noise counterfactual protocol

The baseline validation gate must pass before any causal perturbation is run.

## Immutable source and archived data

- Seed 9: `If1e-02_Ip40_In05_R50_num100_wid100_n9_BEM.mat` in `/Volumes/data/BEM_bichamber/archive_one_R_60/BEM_data/`.
- Source: archive-parent `BEM_circular_parallel_instance.py`, SHA256 `2a415d26ad83cf42a2f3cb0f4ba156b8dd6c11b32972ce8305499eebced22da4`.
- Original `dt=0.001`, N=100, R=7.158910531638177, Ip=4, In=0.5, Iw=0.94, If=0.01; no altered interactions, safety conditions, or hydrodynamic approximations.
- Validation runtime: Python 3.14.5, NumPy 2.5.4, SciPy 1.18.1, macOS 15.6.1 arm64. BLAS thread counts are set to one.

## Gate 1: recover angular increments, 309–320

For steps n=309000,…,319999, apply the actual production deterministic update to the archived state. Recover

`xi_i[n] = theta_i[n+1] - (theta_i[n] + dt * Omega_i,effective[n])`.

The parentheses follow the simulator's floating-point addition order. These are angular increments `In*dW`, not dW itself. Do not wrap them, smooth them, rescale them, or replace them with new Gaussian samples. Recovery from rounded headings does not guarantee recovery of the original RNG bits.

The adapter extracts the actual statements from `step_Computation` using the Python AST. The only substituted physical input is `wNoise=GetNoise(para)`, replaced by the supplied increment. Original source files remain unchanged. Allocation is limited to a two-frame buffer; checkpoint I/O and the outer simulation loop are omitted. Every interaction call, summation order, clipping operation, collision override, and positional retry is retained. The extracted code is saved for inspection.

The positional-safety memory `dis_last` is the previous **uncorrected proposed** signed distance, not necessarily the distance of the corrected archived position. Its value at the start is reconstructed from the immediately preceding archived state. Subsequently it is carried forward exactly as in production.

Check one-step positions independently against the archive and velocities against `rdot[:,:,n+1]`. Heading round-trip equality alone is not independent validation because noise was inferred from the headings.

## Gate 2: autonomous baseline

Start once from the archived t=309 state, then evolve without re-anchoring using the recovered per-fish increments. Predeclared tolerances at every native timestep are:

- maximum per-fish Euclidean position error ≤ 10^-6 R;
- maximum circular heading error ≤ 10^-6 radians.

Stop at the first violation and report it. Do not loosen tolerances, reset states, alter the recovered noise, or proceed to perturbed trajectories to obtain an apparently successful counterfactual.

## Conditional next stage, only if both gates pass

At t=309 and 310, solve the local minimum-Euclidean-norm heading change subject to exact complex-p preservation and exact K change ±0.02. Positions remain fixed. Verify nonlinear constraints and optimizer stationarity; a linear projected direction alone is not sufficient to preserve p exactly.

For each target amplitude, construct an equal-heading-RMS control with exact p preservation and zero first-order K change. Controls must be chosen without observing their outcomes. Integrate these short trajectories only to t=320 with the same corresponding increment sequence and the exact source update. Compare to the passing baseline, reporting K, P, signed COM circulation and COM geometry. No long continuation or parameter sweep is authorized by this protocol.

If the baseline fails, these downstream steps are **not performed**. This failure supplies no evidence for or against K being causal.
