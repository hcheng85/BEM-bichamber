# Full-state committor pilot at no-hydrodynamics gap 20

## Decision

This bounded experiment does **not resolve whether population imbalance `q` is
a sufficient reaction coordinate**.  The obstruction is specific: basin
first-passage times are too long for the documented 100-time-unit safety
horizon.  Of 128 independent shots, 26 first hit basin A, 12 first hit basin
B, and 90 (70.3%) remained strictly between the basins.  Censored shots were
not called returns to A.

There is suggestive fixed-`q` heterogeneity at `q=0`, but it is not established:
only two of four states supplied at least two resolved shots, the
shooting-variance-corrected state variance was 0.375, and its original-seed
bootstrap interval was 0--0.375.  The predeclared selective-expansion criterion
was not met at any `q`.  More 100-unit shots would mostly refine a
finite-horizon unresolved fraction, rather than identify the eventual
committor.  The additional orientation test was therefore not opened.

## Definition and prospective state inventory

The forward committor is

\[
h(X)=\Pr_X(\tau_B<\tau_A),\qquad
A=\{q\leq-0.6\},\quad B=\{q\geq+0.6\}.
\]

All selected configurations were strictly between A and B.  Four states were
selected at each exact population-lattice value `q=-0.4,-0.2,0,+0.2`, using 16
distinct original seeds.  For each level the rule selected the earliest exact
`q` frame after time 200 from the first unused eligible seeds, restricted to
times 200--900.  It used no archived future outcome.  Archive coverage was
ample before selection: 17--19 eligible seeds and 1,058--1,999 frames per
level.

Each state received eight independent angular-noise realizations.  A basin was
tested after every native `dt=0.001` step and a shot stopped at its first hit.
The 100-unit horizon is a computational safety boundary, not a physical return
definition.  The compact shot table records A, B, or right-censored status.

## Restart and symmetry audits

The shooting code reuses the validated production-step extractor.  It
reconstructs previous-step wall/collision memory from archived step `k-1` and
retains that memory at restart.  Native `k` to `k+1` replay passed on all 16
states, including collision cases, with maximum state/velocity error
`3.33e-16`.  Population remained exactly 100 in every shot.

The microscopic one-step left/right reflection test passed with maximum error
`7.11e-15`.  The much weaker statistical check,
`mean h(-0.2)+mean h(+0.2)-1`, was 0.229 with a seed-bootstrap interval
-0.125--0.667, so the shooting data neither detect nor contradict reflection
symmetry.

## First-hitting results

| initial q | states | shots | A | B | censored | mean resolved fraction | equal-state empirical committor bounds |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -0.4 | 4 | 32 | 8 | 0 | 24 | 0.250 | 0.000--0.750 |
| -0.2 | 4 | 32 | 9 | 2 | 21 | 0.344 | 0.063--0.719 |
| 0.0 | 4 | 32 | 8 | 3 | 21 | 0.344 | 0.094--0.750 |
| +0.2 | 4 | 32 | 1 | 7 | 24 | 0.250 | 0.219--0.969 |

The last column assigns every censored realization first to A and then to B;
it is an empirical extreme-outcome bound for these shots, not a confidence
interval for the underlying probability.  Binomial intervals in
`state_committors.csv` apply only to resolved shots.  Because resolution can
depend on destination, resolved-only proportions need not estimate the full
committor without bias.

The data show the expected directional tendency with `q`: no B hits were seen
from `q=-0.4`, whereas seven B and one A hit were seen from `q=+0.2`.  They do
not locate an intermediate-committor surface precisely.  At `q=0`, one state
resolved all eight shots with one B hit, while another produced two B hits and
six censored shots.  That contrast motivates the heterogeneity diagnostic but
does not survive the required seed-level precision gate.

## Sufficiency, closure, and next inference

This pilot validates the full-state restart and first-hitting machinery, but it
does not establish a full-state committor surface.  Consequently:

- `q` sufficiency is **unresolved**, not accepted.
- Collective orientation is not justified as an additional reaction
  coordinate by this experiment; it was intentionally not tested after the
  gate failed.
- Even a future finding `h(X)\approx h(q)` would establish committor sufficiency
  only.  It would not prove that `q(t)` is Markovian or kinetically closed.

The limiting quantity is the first-hitting-time distribution relative to the
safety horizon, not the number of archived configurations or native restart
accuracy.  Resolving the scientific question would require a separately
justified longer-horizon/censoring design that preserves terminal states or a
valid survival/competing-risk treatment.  This run stops at the requested
decision gate rather than expanding shots that fail to resolve either basin.

## Reproduction and compact outputs

From `results/nohydro_gap20_descriptive_transition`:

```bash
uv run --with numpy --with scipy --with matplotlib \
  python run_full_state_committor.py --self-test
uv run --with numpy --with scipy --with matplotlib \
  python run_full_state_committor.py --stage audit
```

The complete pilot is reproducible with `--stage pilot`; completed shots are
cached locally and skipped on rerun.  Raw branch trajectories and the local
shot cache are excluded from version control.  Committed outputs contain the
audit, prospective selection and coverage, per-shot terminal summaries,
state-level censoring bounds, fixed-`q` variance diagnostic, and the compact
figure.
