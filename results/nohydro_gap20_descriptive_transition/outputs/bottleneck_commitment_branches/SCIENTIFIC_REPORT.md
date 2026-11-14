# Controlled bottleneck commitment pilot

## Scope

This bounded experiment restarted 14 prospective pre-entry configurations and 14 post-first-crossing configurations from distinct original no-hydrodynamics gap-20 seeds. Each microscopic state used eight independent angular-noise arrays over five time units. The same explicit array was applied to every arm of a state/clone pair. No full-length simulation was run, no archived event definition was changed, and selection did not use future crossing or continuation labels.

The endpoints are short-time entry and propagation outcomes. They are not a whole-school switching committor.

## Restart and structural audits

The unchanged archived source has SHA256 `5e57b9dca2cb962e2c35f540f3982d5f4018eb719067223e6b319866b4482284`. Its public production function does not accept a restart state or exogenous noise. The additive adapter extracts the native timestep statements and supplies only the restart state, previous provisional wall-distance memory, explicit angular-noise increment, and registered social intervention.

- Native parameters were reproduced: `dt=0.001`, `N=100`, gap ratio 20 (`gap=2.8635642126552705`), radius ratio 50 (`R=7.158910531638177`), `Ip=4`, `In=0.5`, `Iw=0.94`, and `If=0`.
- Twenty archived one-step transitions were replayed. Maximum state error was `4.44e-16`; maximum velocity error was `3.33e-16`. Four audited frames exercised collision overrides, affecting 26 fish in total.
- The inferred angular increment had mean framewise standard deviation `0.01572`, consistent with `In*sqrt(dt)=0.01581`.
- The Delaunay-edge visibility optimization was exactly identical to the archived all-pairs `CheckVoronoi` result on all 20 frames. Wall avoidance, social turning, torque clipping, collision handling, and previous provisional-distance memory therefore remained the archived rules.
- Cross-divider social turning was structurally present: 81 admissible cross-gap pairs across the 20 frozen frames (37 pre, 44 post). Ten frames had a nonzero cross-divider torque contribution. Removing only opposite-chamber numerator terms while retaining the original full-neighborhood visual normalization changed social torque by mean absolute `0.0197` and maximum absolute `2.20`.

The parity gate passed and the cross-gap treatment was therefore run.

## Stage I: imposed orientation before entry

The baseline produced 15 entries in 112 clone trials, a seed-equal probability of `0.1339` (95% seed-bootstrap CI `0.0357–0.2502`). The heading interventions changed only source-majority headings at fixed positions.

| Arm | Entries / 112 | Entry probability difference vs baseline | Mean net transfer difference |
|---|---:|---:|---:|
| toward | 22 | `+0.0625` (`0.0179, 0.1250`) | `+0.366` (`0.098, 0.696`) |
| away | 12 | `-0.0268` (`-0.0536, 0.0000`) | `-0.491` (`-1.009, -0.054`) |

The imposed cue changed polarization monotonically: toward minus baseline was `+0.0805` (`0.0636, 0.0973`), and away minus baseline was `-0.0858` (`-0.1030, -0.0683`). Under common noise, seven trials entered only after the toward intervention; three baseline entries were lost under the away intervention. No paired trial changed in the opposite direction.

Thus gap-directed orientation is causally relevant under this controlled perturbation. This does **not** establish that native social torque is the natural source of the spontaneous orientation precursor; wall turning, school geometry, and other native dynamics were not source-ablated in Stage I.

## Stage II: propagation after first crossing

Under full dynamics, 68 of 112 clone trials transferred at least four additional net fish. The seed-equal propagation probability was `0.6071` (95% CI `0.3482–0.8304`), and mean additional transfer was `9.82` fish (`4.87–15.14`).

The surgical cross-gap-off arm produced 67 successes. Its probability difference was `-0.0089` (`-0.0268, 0.0000`): only one paired clone changed threshold outcome, so an effect on the four-fish success probability is not statistically resolved by this pilot. However, mean additional transfer decreased by `1.46` fish (`-2.66, -0.44`); 12 of 14 seed means were negative and two positive. Cross-divider social input therefore has a resolved positive effect on the *amount* of five-unit propagation in these selected post-entry states, while evidence for changing the binary four-fish commitment threshold is insufficient.

The supported physical picture is two-part: orientation toward the opening affects geometric initiation, and cross-divider social feedback modestly amplifies subsequent transport. It is not the sole propagation mechanism, because the full success rate remains high and the ablation changes only one binary threshold outcome. Within-chamber recruitment, approach-zone storage, wall geometry, and self-propulsion remain active.

## Interpretation limits

- Stage I demonstrates sensitivity to an imposed heading cue, not the natural dynamical origin of that cue.
- Stage II conditions on first crossing and tests at least four additional net fish over five time units. It does not measure full-school relocation, basin arrival, return probability, or a complete switching rate.
- Confidence intervals use the 14 original trajectory seeds as independent units; the eight noise clones are paired conditional replicates.
- The pilot was not enlarged after observing outcomes. Wider independent-seed studies would be needed to resolve small changes in the binary propagation probability.

Compact machine-readable results are in `summary.json`, `per_microstate.csv`, `branches.csv`, `simulator_audit.json`, and `selection_manifest.csv`.
