# Prospective source-level test of the planar invariant and a genuine turn trigger

Use this prompt in the full local BEM environment, where the original 100-fish
archive, geometry, and LU matrices exist. **Do not require copying the 3.7 GB
trajectory archive to this repository.** This is a targeted counterfactual,
not another broad phase-diagram fit.

> Investigate the *causal trigger for approach to the bichamber opening* using
> the newly derived exact planar-wall coordinate, not correlations of the
> instantaneous swimming velocity with the gap direction.
>
> 1. Read the full `RESULT_AND_LIMITS.md` and scripts in
>    `results/fp_gap0p5_hydro/approach_mechanism/deterministic_gate_threshold`.
>    Verify every wall equation against `BEM_bichamber_parallel_instance.py`
>    and the actual angular-noise discretization, including ray switching,
>    near-wall amplification, clipping and hydrodynamic overrides. Do not
>    silently replace the production wall law.
>
> 2. From the 20 archived `Iw=.94, In=.5, Ip=4, If=.01` trajectories, select a
>    **prospective** upstream section without using any future crossing outcome.
>    For example, take the first point per school oscillation at which an
>    approaching focal fish is in the solid flat-divider ray branch,
>    `d=|x|-r` is in the predeclared range `[1.5,2.5]`, its heading points
>    into the divider and vertically toward the nearest tip, and
>    `|y|>gap/2+r`. Report how many sections these fixed rules identify.
>    If this creates too few samples, diagnose and propose one transparent,
>    pre-outcome change, then rerun. Do not select windows based on their
>    subsequently known successes.
>
> 3. For each selected history, compute the full-state decomposition while
>    the **same** straight-wall ray, heading sign and unclipped rule persist:
>    `H=abs(theta_local)+Iw*log(d)` and
>    `dH=[sgn(theta_local)*Omega_other_local-Iw*u_normal/d]dt
>        + In*sgn(theta_local)dW_local`, where `u_normal` is translational
>    velocity **toward the divider** (not radial velocity toward the gap).
>    Check left/right sign conventions explicitly; e.g. for left chamber
>    `u_normal=u_x`, and for right `u_normal=-u_x`. Distinguish the sign of
>    `Omega_other_local` in each chamber. Compare measured `dH` to the
>    integrated drift and reconstructed angular noise at per-step resolution.
>    Report exact-branch duration and the fraction of valid time; do not
>    extrapolate the identity across branch changes.
>
> 4. For a prespecified approach distance `dc=r` and a fixed short horizon
>    `T=4` simulation-time units, compute the first reach to `dc`, first tip
>    ray switch, first gap crossing and the sign/time decomposition of
>    cumulative `Delta H`. Test the exact wall-only first-passage bound in
>    `RESULT_AND_LIMITS.md` **only on the wall-only isolated control**.
>    For the full school, report how large the explicit other-torque and
>    hydrodynamic residual terms are relative to `In*sqrt(T)`. This decides
>    whether a Brownian streamline mechanism is relevant or quantitatively
>    invalidated by collective drift.
>
> 5. Before any expensive grid, run a small **paired full-school** experiment:
>    roughly 8-12 eligible prospective 100-fish states, `Iw` in
>    `{0.4,0.94,1.5}`, `In` in `{0,0.5}`, horizon `T=4`, retaining the
>    actual social, visibility, translational/rotational hydrodynamic and
>    collision/ray override updates. Use common initial states and matched
>    angular-noise streams where noise is on. Change parameters in the
>    source, not in frozen torque postprocessing. Record the **first
>    turning impulse**, its terms, approach to `dc` and first crossing.
>    Only if a systematic mechanism survives this check, extend a few paired
>    histories to at least one complete oscillation and to the 70-fish
>    school-wide commitment criterion.
>
> 6. Do not fit an effective potential, threshold, rate or noise coefficient.
>    If the first-principles planar reduction fails for the actual school,
>    clearly falsify it, quantify why, and identify the minimal alternative
>    coordinate directly from the source equations. In particular, do not
>    call a swimming-current decomposition a mechanism for turning.
>
> 7. Commit **small** scripts, CSVs and figures, a data provenance note, and
>    one concise report containing derivations and controls. The report must
>    distinguish an exact isolated-branch theorem, a numerical diagnostic,
>    a causal paired-rollout result, and an unproven proposal. Do not push
>    giant raw trajectories or numerical LU objects. Include command lines,
>    parameters, SHA, units and failure cases.
>
> **Stop condition:** If full-state trajectories/LU files are unavailable,
> state the missing path(s) and stop. If the residual dominates `dH`, do not
> iterate more scalar Fokker-Planck approximations or declare a paper theory;
> report that the wall-only closure is inapplicable and the triggering
> interactions must be resolved first.
