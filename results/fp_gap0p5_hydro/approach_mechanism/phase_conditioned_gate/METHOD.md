# Method and reproducibility

## Scope and provenance

This is a read-only analysis of the existing 20 full-BEM files

    /Volumes/storage2/BEM_bichamber/archive/archive_one_R_60/BEM_data/If1e-02_Ip40_In05_R50_num100_wid50_n{seed}_BEM.mat

for seeds 0--19. The accepted first-crossing episode definitions are from:

    results/fp_gap0p5_hydro/approach_mechanism/tables/nucleation_episodes.csv

The analysis code is:

    results/fp_gap0p5_hydro/approach_mechanism/phase_conditioned_gate/scripts/analyze_phase_conditioned_gate.py

It reuses the geometry conventions validated in:

- 'results/fp_gap0p5_hydro/approach_mechanism/divider_tip_geometry/scripts/exact_tip_dynamics.py'
- 'results/fp_gap0p5_hydro/approach_mechanism/divider_tip_geometry/REPORT.md'
- 'results/fp_gap0p5_hydro/approach_mechanism/inward_current_decomposition/scripts/analyze_inward_current.py'

No BEM solve, simulator run, fitted coefficient, potential, or outcome-trained classifier is used. The corrected counterfactual tangent-velocity implementation remains the authoritative source for one-sided drift classifications; the present bounded test does not recompute or use a one-sided drift.

## Imported analytical files

The working branch starts from 'analysis/divider-tip-one-sided-drift' at '0daeb3f'. Exactly these files were imported from theory commit '161f777':

- 'PHASE_CONDITIONED_GATE_THEORY.md'
- 'finite_window_gate.py'
- 'test_finite_window_gate.py'

All five numerical tests were run directly as assert functions because 'pytest' was unavailable in the offline environment:

    UV_CACHE_DIR=/private/tmp/uv-cache uv run --offline \
      --with numpy --with scipy --with matplotlib \
      python -c "import test_finite_window_gate as t; \
    names=[n for n in dir(t) if n.startswith('test_')]; \
    [(getattr(t,n)(),print('PASS',n)) for n in names]; print('TOTAL',len(names))"

The five passed checks are:

1. tangent reflection symmetry and closed-cone rejection;
2. probability bounds and limiting behavior;
3. constant-flight-drift reflection symmetry;
4. agreement of the general and zero-drift eigenseries;
5. quadrature convergence.

These tests verify the numerical implementation of the stated toy model, not its closure for BEM trajectories.

## Geometry and conventions

The archive has 'gap=R', so 'gap/R=1' and 'gap/(2R)=0.5'. Saved nondimensional lengths are \(g=R=7.158910531638177\) and tip radius \(r=0.5\). For origin sign \(s=+1\) in the left chamber and \(s=-1\) in the right chamber,

\[
d=-s x,\qquad
\theta_{\rm local}=\operatorname{atan2}(\sin\theta,s\cos\theta).
\]

The rounded-tip centers are \(c_-=-g/2-r\) and \(c_+=g/2+r\). The tangent interval is

\[
\alpha_-=\operatorname{atan2}(c_--y,d)+\arcsin(r/\sqrt{d^2+(c_--y)^2}),
\]
\[
\alpha_+=\operatorname{atan2}(c_+-y,d)-\arcsin(r/\sqrt{d^2+(c_+-y)^2}).
\]

A corridor is valid only for \(d>0\), focal-tip distance greater than \(r\), and \(\alpha_-<\alpha_+\). Its angular margin is

\[
m={\alpha_+-\alpha_-\over2}
-\left|\operatorname{wrap}\!\left(\theta_{\rm local}
-{\alpha_++\alpha_-\over2}\right)\right|.
\]

Positive margin means inside. The three prespecified shells are \(d=r=0.5\), \(d=2r=1\) (primary), and \(d=4r=2\).

## Actual-crossing analysis

For each of the 512 accepted nucleation episodes:

1. Resolve the first native saved index with \(s x_i\geq0\) within the archived event-time tolerance. If fish tie at that step, use the most inward fish and record the tie count.
2. Define the full-occupancy episode start from 'full_dwell_before'.
3. For the eventual first crosser, find the first inward shell crossing \(d_{n-1}>d_{\rm ref}\), \(d_n\leq d_{\rm ref}\).
4. If it is already inside at episode start, mark 'reference_left_censored'; do not substitute a later outward/re-entry passage.
5. From the shell index through the last pre-crossing state, evaluate both:
   - the corridor frozen at the reference state;
   - the instantaneous corridor at every native step.

For each definition, record fraction of time inside, total/initial/final residence, inside-to-outside exit count, and persistence at all steps. The native interval is \(\Delta t=0.001\).

Self-propelled normal displacement is

\[
\Delta x_{\rm self}=\sum_n s\cos\theta_n\,\Delta t.
\]

The saved-update residual

\[
\Delta x_{\rm nonself}=s[x(T)-x(0)]-\Delta x_{\rm self}
\]

is exact for the realized trajectory. It is direct+BEM translation on ordinary production updates and also absorbs collision/override corrections. It is consequently reported as an **effective non-self translational residual**, despite the compatibility column name 'translational_hydrodynamics'. A complete source split would require long microscopic BEM replay and is outside the bounded task.

The isolated safe-flight time is \(d/\min_{\theta\in[\alpha_-,\alpha_+]}\cos\theta\). It is a sufficient bound only for a unit-speed isolated swimmer continuously remaining in the frozen cone.

## Prospective Poincaré risk set

Native indices after time 200 are sampled every 0.01 time unit only for prospective section/opportunity construction.

A section occurs when:

- all fish occupy one chamber;
- the school-center height crosses \(Y_c=0\) upward;
- chamber identity is unchanged across that sampled step.

After an accepted section, another same-origin upward crossing is accepted only if the school first reaches \(Y_c\leq-r=-0.5\). This is the hysteresis reset. Opposite-origin full-occupancy sections are also admissible; in practice a cycle terminates at its first fish crossing before whole-school relocation.

For each accepted section, the prospective outcome is first fish crossing before the next accepted section return. No future event time or fish identity is used to select the starting state.

Within a cycle and for each fish, the first sampled valid state satisfying \(0<d\leq d_{\rm ref}\) is its shell opportunity. Fish already within the shell at the section are retained and flagged 'censored_inside_shell'. This lets leading membership renew naturally. The primary cycle-level observation is the earliest fish opportunity; exact ties select the closest fish. Initial angular margin, corridor width, distance, height, candidate crossing, cycle crossing, and whether that candidate became the first crosser are stored. This is a shell opportunity, not a claim that the heading is already admissible.

The fish-level table is intentionally git-ignored because fish within a school are correlated and it is a 5.2 MB regeneration artifact. Inference uses one earliest observation per cycle and independent-seed aggregation.

## Uncertainty and descriptive discrimination

Every estimate first averages within seed. A 95% interval resamples the seed means with replacement 20,000 times. For quantities requiring both outcome classes, only seeds containing both classes contribute and their number is reported.

The AUC is computed independently in each seed from the prespecified angular margin, then averaged across seeds. It is a descriptive prospective discrimination statistic, not a fitted classifier. The inside-versus-outside contrast is likewise computed inside each seed before bootstrap. Frames, fish, and cycles are never treated as independent replicates for confidence intervals.

## Reproduction

    UV_CACHE_DIR=/private/tmp/uv-cache \
    MPLCONFIGDIR=/private/tmp/mplconfig \
    uv run --offline --with numpy --with scipy --with matplotlib \
    python results/fp_gap0p5_hydro/approach_mechanism/phase_conditioned_gate/scripts/analyze_phase_conditioned_gate.py \
      --raw-dir /Volumes/storage2/BEM_bichamber/archive/archive_one_R_60/BEM_data \
      --episodes results/fp_gap0p5_hydro/approach_mechanism/tables/nucleation_episodes.csv \
      --out results/fp_gap0p5_hydro/approach_mechanism/phase_conditioned_gate \
      --force

Per-seed checkpoints are ignored by git. Remove '--force' to rebuild aggregate tables and figures from existing checkpoints.

Committed compact outputs are:

- 'tables/seed_level_summary.csv': the independent-seed quantities used for uncertainty;
- 'tables/crossing_persistence_summary.csv': pooled seed-mean estimates and intervals;
- 'tables/prospective_summary.csv': prospective seed-mean estimates and intervals;
- 'tables/actual_crossing_persistence.csv': one row per event and shell, including censoring;
- 'tables/prospective_cycles.csv': section and stopping indices;
- 'tables/prospective_cycle_opportunities.csv': one earliest opportunity per cycle and shell.

## Limitations

- Only one \((I_w,I_n)\) point exists in this ensemble.
- Reference-shell left censoring rises from 15.2% at \(r\) to 43.8% at \(4r\).
- A finite observed state sampling interval cannot prove continuous-time cone residence between saved states.
- Exit counts are sensitive to native-step boundary chatter; persistence and residence fraction are the primary measures.
- The prospective crossing fraction is high, especially conditional on reaching \(d=r\), which reduces precision of the near-gap inside/outside contrast.
- First-fish nucleation is not collective commitment. Recruitment is intentionally not modeled here.
