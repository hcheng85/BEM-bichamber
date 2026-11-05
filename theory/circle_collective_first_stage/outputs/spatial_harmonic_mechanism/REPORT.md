# Spatial-harmonic test of the circular schooling--milling transition

## Verdict

The proposed decomposition is **geometrically exact**, but spatial first-harmonic
asymmetry alone does not explain the loss of polarization.

Across stationary regimes, polarization `P` and fixed-center `|A1|` are strongly
associated because the schooling state is an off-center, wall-following school,
whereas the milling state is a more centered ring.  This is not a closed
mechanism.  Within either regime their framewise correlations are only 0.200
(schooling) and 0.149 (milling), and the correction `H1` is essential.  In
milling, `H1` is comparable in magnitude to `A1` and usually points against the
tangent contribution.  The low polarization is produced by vector cancellation,
not by `A1` alone.

The single observed transition in seed 9 strengthens this negative conclusion:
`|H1|` changes and coherent signed circulation weakens before `P` collapses,
while fixed-center `|A1|` remains high through the collapse and changes only
later.  Thus fixed-center spatial `m=1` asymmetry distinguishes the stationary
states but does not lead this transition.  These ten trajectories at one
parameter point do not identify a stability boundary.

## Immutable inputs and streaming procedure

The analysis used the ten original files

```text
If1e-02_Ip40_In05_R50_num100_wid100_n{0..9}_BEM.mat
```

with `state[0]=x`, `state[1]=y`, `state[2]=theta`, `N=100`,
`R=7.158910531638177`, and `dt=0.001`.  The large uncompressed MATLAB-v5 arrays
were opened read-only with
`theory/circle_collective_first_stage/matv5_memmap.py`.  Batches of at most
1,000 sampled frames were copied; complete trajectories were never loaded.
The simulator, wall law, safety conditions, MAT files, and previous analyses
were not modified.

After the 20-time-unit formation transient, every 100th native frame was used:
98,000 frames in total.  Persistent schooling/milling labels use the previously
documented one-second order-parameter smoothing and ten-second, 90% persistence
rule.  Seed 9 was separately read every 0.01 time unit over `t=280..340`.
Uncertainty below is a bootstrap over independent seed/state intervals, not
over correlated frames.  There are eight schooling and only three milling
units (seeds 7, 8, and the post-transition part of seed 9), so milling
uncertainty remains limited by trajectory count.

## Exact identity and the role of `H1`

For each frame, about the **fixed circle center**, the code evaluates

```text
phi_i   = atan2(y_i, x_i)
alpha_i = wrap(theta_i - phi_i)
p       = mean(exp(i theta_i))
M       = mean(sin(alpha_i))
A1      = mean(exp(i phi_i))
A2      = mean(exp(2 i phi_i)).
```

With `s=sign(M)`, direct substitution gives

```text
p = i s A1 + H1,
H1 = mean[exp(i phi_i) (exp(i alpha_i) - i s)].
```

The maximum complex residual over all 98,000 sampled frames is
`9.74e-15`; the RMS residual is `1.28e-15`.  This is floating-point precision.
For the 0.14% of frames with `|M|<0.2`, the same algebra remains true after an
arbitrary sign choice, but the tangent-chirality interpretation is explicitly
classified as ambiguous.  Changing the eligibility threshold to 0.1 or 0.3
does not change the conclusions.

The stationary seed-level comparison is:

| quantity | schooling mean (95% seed bootstrap) | milling mean (95% seed bootstrap) |
|---|---:|---:|
| `P` | 0.905 [0.903, 0.907] | 0.196 [0.185, 0.206] |
| `|M|` | 0.904 [0.901, 0.906] | 0.763 [0.750, 0.777] |
| fixed-center `|A1|` | 0.929 [0.926, 0.931] | 0.403 [0.381, 0.422] |
| fixed-center `|A2|` | 0.751 [0.744, 0.757] | 0.174 [0.168, 0.177] |
| `|H1|` | 0.166 [0.160, 0.174] | 0.375 [0.356, 0.396] |
| `|H1|/P` | 0.184 [0.177, 0.193] | 3.08 [2.85, 3.51] |

The global eligible-frame correlation of `P` with `|A1|` is 0.905, but this
mostly reflects pooling two different regimes.  The tangent-only RMS error
`RMS(|H1|)/RMS(P)` is 0.220 in schooling and 1.90 in milling.  In milling the
mean cosine between `i s A1` and `H1` is -0.790, demonstrating systematic
cancellation.  Therefore silently dropping `H1` would give the wrong physical
interpretation precisely in the low-polarization state.

## Fixed circle center versus school center of mass

The fixed center is required for the exact orbital/tangent construction, but it
must not be confused with the instantaneous school center.  The stationary
school center lies at mean radius `0.476 R`; the milling center lies at
`0.133 R`.  When azimuth is instead measured about the instantaneous center of
mass, mean `|A1|` is small in both states: 0.0634 in schooling and 0.0493 in
milling.  No fish was excluded by the `1e-8 R` center-of-mass radial tolerance.

Thus the large fixed-center `A1` and `A2` in schooling primarily quantify the
off-center spatial arrangement of a compact school.  The second harmonic
actually separates the stationary regimes more strongly than the first
(`mill-school = -0.578`, bootstrap interval `[-0.585,-0.569]`), but neither is
by itself an angular-alignment mechanism.

## Seed-9 chronology

The one transition is referenced to the previously detected persistent-state
boundary at `t=312.41`.  With a prespecified rule--the first one-second window
containing at least 70 of 100 values beyond three standard deviations from the
seed-9 `t=280..300` baseline--the onsets are:

| diagnostic | first sustained change |
|---|---:|
| `|H1|` increases | 306.08 |
| `|M|` decreases | 306.57 |
| `P` decreases | 311.81 |
| fixed-center `|A1|` decreases | 321.67 |
| center-of-mass radius decreases | 322.32 |
| fixed-center `|A2|` sustained increase | not detected |

The exact onset times depend on this descriptive threshold, and one transition
cannot establish a universal precursor.  The robust ordering visible without
threshold fitting is that the tangent-deviation term grows and circulation
coherence is lost while `A1` is still large; polarization then collapses; only
later does the school become spatially centered.  Spatial `m=1` asymmetry is
therefore not the initiating transition coordinate in seed 9.

## Orbital-speed identity

The stored convention was used exactly: `rdot[:,:,step+1]` advances
`state[:,:,step]`.  A separate finite-difference audit gives position-update
RMS errors below `1.3e-13`.  For fish away from the fixed-center singularity,

```text
omega_i = (x_i v_yi - y_i v_xi) / r_i^2,
dA1/dt  = i mean(omega_i exp(i phi_i)).
```

For each differencing interval, a fixed set of fish that remained beyond
`0.01 R` was used on both sides.  At native lag `dt=0.001`, two schooling
fish-steps out of 219,000 and 83 milling fish-steps out of 81,000 were excluded.
The complex relative RMSE is `2.69e-4` in schooling and `3.59e-3` in milling.
Across lags 0.001--0.01 it remains below `4.9e-4` and `6.2e-3`, respectively.
The larger milling error is the expected conditioning problem of azimuth near
the fixed center, not a velocity-index mismatch.

Writing `omega_bar=mean(omega_i)`, the right-hand side separates exactly into

```text
i omega_bar A1
+ i mean[(omega_i-omega_bar) exp(i phi_i)].
```

The mean magnitudes of the rigid and differential terms are 0.304 and 0.0128
in schooling; differential motion contributes only 3.5% of their summed
magnitudes.  In milling they are 0.241 and 0.210, with a 45.5% differential
fraction.  Their vector sum has mean magnitude only 0.135.  A milling ring is
therefore neither perfectly centered nor rigidly rotating: differential orbital
speeds are dynamically important and partly cancel rigid rotation in `dA1/dt`.

## Scientific interpretation and limitations

Directly established:

1. The `p = i s A1 + H1` identity is exact.
2. The two stationary regimes have very different fixed-center spatial
   harmonics and center-of-mass locations.
3. `H1` cannot be neglected, especially in milling.
4. Seed 9 loses tangent-heading coherence before fixed-center spatial `A1`
   changes appreciably.
5. The orbital `A1` derivative identity and the `step+1` velocity convention
   pass numerical validation; milling contains large differential orbital
   motion.

Not established:

- `A1` does not close or predict the polarization dynamics by itself.
- The single natural transition cannot determine causality, transition rates,
  or a stability boundary.
- No effective potential, surrogate, or fitted transition model is justified
  by this test.

The defensible next first-principles state must retain heading relative to the
radial/tangent direction (the information in `H1` or its lower moments) together
with spatial centering/orbital geometry.  Controlled data at additional
parameter points and multiple transitions would be required to predict a
stability boundary.

## Reproducibility map

- `analyze_spatial_harmonic_mechanism.py`: complete read-only analysis.
- `analysis_config.json`: sampling, thresholds, velocity indexing, and
  bootstrap configuration.
- `identity_audit.json`, `state_relationship_summary.csv`: exact identity and
  tangent-only diagnostics.
- `stationary_seed_summary.csv`, `stationary_comparison.csv`: seed-level
  statistics and uncertainty.
- `seed9_transition_timeline.csv`, `seed9_transition_chronology.json`: compact
  high-resolution transition diagnostics.
- `orbital_identity_samples.csv`, `orbital_identity_summary.csv`: derivative,
  rigid/differential, exclusion, and discretization audit.
- `harmonic_sampled_frames.npz`: compact sampled reductions only; no raw
  trajectories or per-fish histories.

