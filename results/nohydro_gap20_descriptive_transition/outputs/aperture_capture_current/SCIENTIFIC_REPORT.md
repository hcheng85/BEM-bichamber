# Straight-line aperture capture as a near-term transition gate

## Conclusion

The parameter-free aperture-capture current is supported as a **near-term geometric gate**, not as a closed long-horizon transition coordinate. It improves on the ordinary longitudinal approach current for held-out first entry and four-fish transfer, including comparisons made within the pre-existing `(q,a,b)` cells. It also separates continued and failed crossing attempts during the final unit before entry.

It does not replace collective polarization. Gap-directed polarization remains the best five-unit first-entry discriminator and separates continued from failed attempts earlier than aperture capture. The physical sequence supported by the combined results is therefore:

1. collective orientation supplies fish toward the approach layer;
2. approach geometry selects the subset whose instantaneous rays reach the opening;
3. maintained capture-directed supply supports continued transfer.

Individual ballistic rays are accurate over short travel times and progressively fail as angular turning accumulates. The straight-line current should consequently be interpreted as instantaneous aperture eligibility, not a prediction that the same fish must follow an unturned path for five units.

## Exact geometric construction

Fix the source from the current majority chamber, without using the future crossing direction. Let `source_sign=+1` for a right-chamber source and `-1` for a left-chamber source, and define the forward local coordinate

\[
z=sx,\qquad s=-\mathrm{source\_sign}.
\]

The source-side approach rectangle is

\[
-\ell\le z<0,\qquad |y|\le h,\qquad h=\ell/2,
\]

with the audited physical gap width `ell=2.8635642126552705`. For heading `theta`, the inward velocity is

\[
u=s\cos\theta.
\]

When `u>0`, the unit-speed straight ray reaches the divider after

\[
\tau=-z/u
\]

at transverse coordinate

\[
y_{\rm hit}=y-z\,{\sin\theta\over u}.
\]

The capture condition is exactly

\[
I_{\rm cap}=I(-\ell\le z<0)I(|y|\le h)I(u>0)I(|y_{\rm hit}|\le h).
\]

Because the starting point and divider intercept both lie in the convex interval `[-h,h]`, this condition also guarantees that the straight segment reaches the opening before crossing either lateral edge of the approach rectangle. It is an exact ray-geometry statement; it is not an exact statement about the subsequently turning fish.

The requested current is

\[
J_{\rm cap}={1\over N}\sum_i I_{{\rm cap},i}\,u_i.
\]

It contains no threshold, fitted coefficient, future direction, or outcome-dependent sign. The comparison current `J_a` uses the archived stored velocity of every fish in the same approach rectangle, including outward contributions and collision overrides.

## Prospective evidence on independent seeds

The first-entry cohort contains 5,792 crossing-free-history risk origins, including 343 future entries. The transfer cohort contains 6,308 fixed-grid origins, including 225 outcomes with at least four net fish moving toward the current minority chamber. Seeds 30--39 are held out; seeds 0--29 are used only to reproduce the previous training-only `b` split inside the declared `(q,a)` cells.

Held-out rank discrimination is:

| coordinate/current | first-entry AUC | four-fish-transfer AUC |
|---|---:|---:|
| `-|q|` | 0.604 `[0.514, 0.687]` | 0.554 `[0.460, 0.649]` |
| `a` | 0.787 `[0.762, 0.807]` | 0.772 `[0.731, 0.813]` |
| `b` | 0.674 `[0.633, 0.707]` | 0.702 `[0.656, 0.738]` |
| ordinary `J_a` | 0.657 `[0.634, 0.681]` | 0.732 `[0.704, 0.768]` |
| aperture `J_cap` | 0.707 `[0.681, 0.733]` | 0.762 `[0.730, 0.796]` |

Relative to `J_a`, the `J_cap` AUC gain is `0.0499`, seed-bootstrap CI `[0.0250, 0.0728]`, for first entry and `0.0301`, CI `[0.0117, 0.0474]`, for transfer. Within the existing `(q,a,b)` cells, the conditional-concordance gains remain resolved: `0.0453` `[0.0258, 0.0624]` and `0.0223` `[0.0064, 0.0395]`, respectively. Thus the aperture filter adds geometric information beyond the ordinary current and the coarse existing coordinates.

For five-unit first entry, `a` exceeds `J_cap` by `0.0801`, CI `[0.0475, 0.1072]`. For transfer their difference is unresolved: `J_cap-a=-0.0091`, CI `[-0.0614, 0.0468]`. This distinction is physically sensible: `a` measures collective orientation before fish necessarily enter the approach rectangle, whereas `J_cap` measures a localized, later gate.

The gate is sparse. Only 217 of 5,792 first-entry risk origins have nonzero `J_cap`, but 139 of those 217 lead to an entry within five units. In the transfer cohort, 112 of the 322 nonzero origins produce a four-fish outcome. These are prospective associations, not probabilities of individual ballistic arrival.

## Matched entry and continued-versus-failed attempts

At exactly matched population imbalance, future-entry configurations have more capture-directed current than stable no-entry controls:

- two units before entry: `Delta J_cap=0.01033`, CI `[0.00807, 0.01264]`, and `Delta n_cap=1.57` fish;
- one unit before entry: `Delta J_cap=0.03062`, CI `[0.02620, 0.03548]`, and `Delta n_cap=4.95` fish.

These exact-`q` pairs are retrospective case-control contrasts; they verify the geometric precursor but are not unbiased prospective hazard estimates.

For the fixed 337-event cohort, continued and failed curves are standardized to the same initial-`|q|` distribution. The continued-minus-failed `J_cap` difference is slightly negative at lead 3 (`-0.00240`, CI `[-0.00473, -0.00046]`), unresolved at lead 2 (`0.00163`, CI `[-0.00241, 0.00539]`), then becomes positive at lead 1 (`0.03122`, CI `[0.02445, 0.03796]`), lead 0.5 (`0.05190`, CI `[0.04369, 0.06015]`), and lead 0.2 (`0.06272`, CI `[0.05463, 0.07135]`). The number of ray-captured fish differs by `4.25`, `7.10`, and `8.51` fish at those final three leads.

Polarization already differs at lead 3, while positive capture-current separation appears only within the last unit. Aperture capture is therefore the geometric conversion of earlier collective organization into immediate opening access, not the source of that organization.

## Is the straight-line approximation physically reasonable?

The audit follows every ray-qualified fish from the prospective first-entry origins to its own ballistic arrival time. No outcome label enters the ray definition.

| ballistic time `tau` | fish | seeds | median absolute heading change | median endpoint error |
|---|---:|---:|---:|---:|
| 0--0.25 | 28 | 11 | 0.162 rad | 0.014 |
| 0.25--0.5 | 27 | 14 | 0.117 rad | 0.036 |
| 0.5--1 | 153 | 24 | 0.190 rad | 0.093 |
| 1--2 | 383 | 31 | 0.253 rad | 0.194 |
| 2--4.1 | 194 | 32 | 0.485 rad | 0.724 |

As a complementary retrospective audit, the fish that actually makes the next first crossing is already ray-qualified at the prospective origin in 100% of events within 0.5 time units, 85% within 0.5--1, 69% within 1--2, and only 7.9% within 2--5. The sample is modest at the shortest times, but the degradation is unambiguous.

Therefore straight-line capture is physically reasonable for short approach flights, useful but approximate through about two units, and not a defensible individual-path approximation over the full five-unit outcome horizon. Its longer-horizon predictive value reflects the current organization of supply, not literal persistence of a particular ray.

## Minimal finite-angular-persistence correction

The archived microscopic law is

\[
d\theta=\Omega_i\,dt+I_n\,dW_t,\qquad
\Omega_i=\operatorname{clip}(\omega_{\rm wall,i}+\omega_{\rm social,i}),
\qquad I_n=0.5,
\]

with the production clipping bounds. For a short flight, hold the instantaneous `Omega_i` locally constant and linearize about the initial heading. If `u=cos(phi)>0` in local source coordinates and `tau=-z/u`, the transverse intercept perturbation is

\[
\delta y_{\rm hit}={1\over u}\int_0^\tau\delta\theta(t)\,dt.
\]

Consequently,

\[
E[\delta y_{\rm hit}]={\Omega_i\tau^2\over2u},\qquad
\operatorname{Var}(\delta y_{\rm hit})={I_n^2\tau^3\over3u^2}.
\]

Writing `mu_i` for the ballistic intercept plus the mean shift and `sigma_i` for the square root of this variance gives the parameter-free aperture probability

\[
p_i^{\rm cap}=\Phi\!\left({h-\mu_i\over\sigma_i}\right)
-\Phi\!\left({-h-\mu_i\over\sigma_i}\right),
\]

and the minimal persistence-corrected current

\[
J_{\rm cap}^{\rm pers}={1\over N}\sum_{u_i>0}u_i p_i^{\rm cap}.
\]

The supplied analysis implements this formula but does not fit or tune it. It uses only the native wall/social torque, noise strength, and geometry. It is a controlled local correction when angular changes are small; it is not valid once the state-dependent social and wall torques evolve substantially. Extending beyond that regime requires a first-passage calculation for the microscopic `(y,theta)` dynamics with state-dependent torque, not another empirical collective coordinate.

## Scope and reproducibility

- `analyze_aperture_capture_current.py` contains the ray condition, current, cohort tests, persistence audit, and analytical correction.
- `aperture_capture_summary.png` is the compact geometry-and-evidence figure.
- `predictive_auc.csv`, `matched_and_event_summary.csv`, `persistence_summary.csv`, and `first_fish_recall.csv` are the compact statistical summaries.
- `summary.json` records the held-out comparisons and correction formula.

The analysis used the existing 40 trajectories. It did not modify the simulator, run new trajectories, fit a potential, create a new phenomenological coordinate, or tune a coefficient.
