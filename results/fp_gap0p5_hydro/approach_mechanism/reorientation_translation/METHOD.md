# Method

## Geometric decomposition

`sigma_perp` is reconstructed exactly from the archived coordinates as
`Lambda*d`. Log differences are formed within each event/control pair, averaged
within seed, and bootstrapped over the 20 seed means.

## Logistic comparison

The previous samples and fixed split are unchanged: seeds 0--14 train and
15--19 test. No model uses outcome-derived features. `psi` is encoded by its
sine and cosine; all other features are linear and standardized using training
data. Pooled metrics, metrics for each held-out seed, and paired differences of
seed AUCs are reported.

## Native-step kinematics

For each saved segment, update origins are `A-2 <= t < A-1`. The archived
`rdot[:,:,k+1]` advances `state[:,:,k]` to `state[:,:,k+1]`. Direct dipolar
translation is evaluated at every native origin using the original complex
velocity equation. Saved updates equal to zero are wall-freeze overrides;
updates equal to the heading vector are self-only hydrodynamic overrides.
Contributions suppressed by those rules are set to zero. Effective BEM velocity
is the saved velocity minus the applied self and direct contributions. This is
the boundary-mediated velocity actually retained by the simulator, rather than
the raw pre-override field.

The identity is checked per native update and after integration. Direct
velocity was independently compared with the previous equation replay at three
snapshots; maximum component error was `3.33e-16`.

## Distance-matched risk set

Risk time is `t0=A-2`; outcome is first crossing within the next two time
units. Controls need remain in the origin chamber only through this horizon.
Matching is within seed and chamber using `(sin phase, cos phase, d)`, with no
heading variable. Controls are unique and separated by at least 0.5 time units.
The heading-only model retains the same 15/5 seed split.
