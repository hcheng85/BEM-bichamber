# Encounter geometry versus collective reorientation

## Physical conclusion

Commitment is not explained by a school merely encountering the opening with
an unchanged heading, and it is not explained by heading rotation at a fixed
position.  The supported mechanism is a **coupled curved encounter**:

```text
coherent source school
    -> COM moves toward the opening and toward its centreline
    -> mean heading rotates with the changing COM-to-opening bearing
    -> first entry either develops into sustained transfer or retreats.
```

The position-aware alignment becomes predictive earlier than the old fixed
horizontal alignment.  An exact finite-change identity then assigns roughly
half of its late separation to rotation of the mean-heading direction and half
to rotation of the COM-to-opening bearing.  The coherence-magnitude term is
small and negative.  Thus the earlier milling signal is most consistently
interpreted as transient turning during an oblique encounter, not a persistent
circulation state and not simple depolarization.

Geometry is nevertheless not a complete replacement coordinate at the start
of a prospectively defined attempt: on held-out seeds its improvement over
fixed-axis alignment is unresolved there.  The evidence supports geometry plus
reorientation as a mechanism during the attempt, rather than a static
geometric gate that fixes the outcome at attempt onset.

## Definitions and audits

For source-chamber fish, let

```text
R = source-school centre of mass,
h = mean heading vector,
g = -R / |R|,
A_fixed = -s h_x,
A_COM = h . g,
```

where `s=+1` for a right source and `s=-1` for a left source.  `A_COM` uses the
actual direction from the school COM to the central opening; it is a proposed
replacement for `A_fixed`, not an additional fitted coordinate.  The analysis
is reflection covariant and the synthetic reflection test passes.

The established cohort is unchanged: 337 first majority-to-minority crossing
episodes from 37 independent seeds, including 185 episodes with at least four
additional net fish transferred in the following five time units and 152
failures.  All measurements discussed below precede the first crossing.
Results are controlled in prespecified `|q|` and COM-distance bins, with
uncertainty resampled over original seeds.  Seeds 0--29 and 30--39 remain the
training and held-out cohorts.

No microscopic simulations were run and no outcome was used to define an
attempt.  The code reads the existing 40 trajectories.

## Earlier held-out prediction

Within the fixed crossing cohort, the position-aware direction distinguishes
later continuation well before the fixed horizontal direction does:

| lead before crossing | held-out fixed AUC | held-out COM AUC | COM - fixed, 95% seed interval |
|---:|---:|---:|---:|
| 3.0 | 0.454 | 0.604 | +0.150 `[0.036, 0.325]` |
| 2.0 | 0.509 | 0.713 | +0.203 `[0.073, 0.401]` |
| 1.0 | 0.691 | 0.834 | +0.143 `[0.054, 0.288]` |
| 0.5 | 0.810 | 0.824 | +0.014 `[-0.076, 0.189]` |
| 0.2 | 0.834 | 0.824 | -0.010 `[-0.103, 0.163]` |

These are one-variable, `q`- and proximity-stratified rank diagnostics, not a
fitted classifier.  They show that actual encounter bearing contains earlier
information.  Near entry, when the school is already funneled toward the
neck, the fixed and position-aware definitions become empirically equivalent.
Because this cohort is selected to contain a future crossing, the table is a
pre-crossing commitment prediction, not a prospective estimate of attempt
incidence.

## Exact translation--turning decomposition

For two pre-crossing times, the finite change obeys the exact identity

```text
Delta(h . g)
  = Delta h . mean(g) + mean(h) . Delta g.
```

Writing `h=P m`, where `P=|h|` and `m` is the unit mean-heading direction,
gives a second exact split:

```text
Delta h . mean(g)
  = Delta P mean(m) . mean(g)
    + mean(P) Delta m . mean(g).
```

Both numerical identities close to machine precision (`3.4e-16` and
`1.7e-16` maximum residual).  From lead 2 to lead 0.5, the standardized
continued-minus-failed changes are

| contribution | difference | 95% seed interval |
|---|---:|---:|
| mean-heading vector | +0.0984 | `[0.0708, 0.1232]` |
| coherence magnitude `P` | -0.0062 | `[-0.0077, -0.0047]` |
| mean-heading direction `m` | +0.1046 | `[0.0777, 0.1297]` |
| COM-to-opening bearing `g` | +0.0975 | `[0.0830, 0.1133]` |
| total `A_COM` | +0.1959 | `[0.1671, 0.2234]` |

Over the same interval, successful episodes reduce COM distance to the
opening by an additional `0.376` length units (`[-0.458,-0.297]`).  The axial
COM difference is unresolved at this point, so most of the early approach is
oblique/centreline convergence rather than purely horizontal translation.
By lead 0.1 the directional and bearing contributions remain comparable
(`+0.123` and `+0.108`), and the additional distance closure reaches `0.524`.

The identities establish *how* the alignment changes, not *why* the school
turns.  The outcome association is empirical; it does not by itself identify
wall torque or social torque as the cause.

## Prospective approach attempts

An attempt begins when the source COM enters the previously established
conservation-law approach rectangle after at least five crossing-free time
units.  Attempts have non-overlapping two-unit outcome windows.  Only current
geometry and past population history enter this definition; future crossings
only label the already-detected attempt.

Across all 40 seeds there are 1,080 attempts: 220 have a target-directed first
entry within two time units, and 138 reach the established four-additional-fish
commitment threshold.  Thirty-three attempts instead first contain a crossing
toward the majority chamber and are counted among the 860 attempts without a
target entry.

On the ten held-out seeds (293 attempts), after matching `q` and proximity:

| outcome | fixed AUC | COM AUC | inward-COM-speed AUC | COM - fixed |
|---|---:|---:|---:|---:|
| target entry | 0.656 | 0.659 | 0.646 | +0.003 `[-0.083,0.093]` |
| commitment | 0.696 | 0.764 | 0.745 | +0.069 `[-0.063,0.194]` |

Thus encounter geometry, fixed alignment, and translation all carry
prospective information, but the independent advantage of `A_COM` at attempt
onset is not resolved.  The strong earlier advantage in the event-aligned
cohort develops during the approach.  This is the key reason not to interpret
`A_COM` as a static commitment coordinate.

## Renewal interpretation

The empirical equal-seed estimates are an attempt rate of `0.0438` per time
unit `[0.0402,0.0474]`, a target-entry probability of `0.193`
`[0.146,0.243]`, and a four-fish commitment probability per detected attempt
of `0.119` `[0.086,0.152]`.  The resulting committed-attempt rate is `0.00625`
per time unit `[0.00437,0.00847]`.  Equal-seed summaries use different
weighting, so their displayed means are not required to multiply exactly.

Inter-attempt intervals have mean `28.84`, coefficient of variation `0.624`,
and no detected lag-one correlation (`-0.004`, `[-0.060,0.055]`).  Commitment
labels retain weak positive lag-one dependence (`0.109`, `[0.003,0.210]`).
Consequently a Poisson process with an iid constant success probability is not
yet justified.

The simplest analytical next model is a renewal/first-passage description of
a coherent aggregate using only its COM position and **relative bearing**
`psi = angle(h)-angle(g)`.  A return flux through the approach boundary would
set the attempt frequency; coupled translation and rotational persistence in
`(R,psi)` would set target-entry and commitment probabilities; retreat would
reset an attempt.  This is a concrete first-principles target, not an
established theory.  It must reproduce the non-exponential waiting law and
explain or remove the residual outcome-to-outcome memory before it can predict
a collective transition rate.

## What is established and what is not

- **Exact:** the two finite-change identities and their numerical closure.
- **Observed:** successful transfer combines COM convergence and genuine
  mean-heading rotation; loss of polarization magnitude does not generate the
  signal.
- **Predicted on independent seeds:** position-aware alignment separates
  commitment earlier within the common crossing cohort.
- **Unresolved prospectively:** replacing fixed alignment by position-aware
  alignment at attempt onset does not give a resolved held-out improvement.
- **Mechanistic interpretation:** a curved geometric encounter requiring
  active reorientation is supported; pure encounter geometry and persistent
  circulation are not.
- **Not established:** causal origin of the turning, an iid renewal process,
  or an analytical collective transition-rate theory.

No potential, transition matrix, additional empirical state hierarchy, or
new microscopic simulation was introduced.
