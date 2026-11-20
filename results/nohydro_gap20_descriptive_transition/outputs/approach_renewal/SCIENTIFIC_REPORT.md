# Approach excursions, recovery, and renewal

## Decision

The established approach-storage geometry provides a clean, physically
interpretable definition of an attempt, but the extracted attempts are **not
an ordinary renewal process**. Complete emptying of both approach layers does
not erase the memory of the preceding school circulation and transfer.

A minimal four-outcome semi-Markov organization—no entry, entry followed by
return, partial transfer, or sustained transfer—predicts held-out recovery
timing and the direction of the next attempt better than an iid renewal model.
It does not resolve prediction of the next sustained transfer, and substantial
waiting-time memory remains after conditioning. It is therefore a useful
episode-level description, not a dynamically closed transition-rate theory.

The strongest physical result is a separation between attempt generation and
commitment. Attempt frequency transfers reasonably from seeds 0--29 to seeds
30--39, but entry, partial transfer, sustained transfer, and transferred fish
do not. Geometry and coherent chamber circulation appear to organize repeated
encounters; the unresolved nonreciprocal turning environment controls their
marks and prevents a universal renewal reward law.

## A threshold-free physical attempt surface

For the already established left and right approach rectangles, let

```text
b_L = (1/N) sum_{left approach}  (ell+x_i)/ell,
b_R = (1/N) sum_{right approach} (ell-x_i)/ell.
```

These are the same weighted storage variables used in the exact conservation
analysis, with `ell=2.8635642126552705`. Define total approach storage
`B=b_L+b_R`.

At the available 0.5-time feature resolution:

1. the last sample with `B=0` is the recovered state;
2. an attempt begins when `B` next becomes positive;
3. it remains active through the complete positive excursion;
4. reset occurs when `B` returns to zero.

The initiating direction is the side with greater storage at the first
occupied sample. This definition uses neither a future crossing direction nor
the outcome. It detects the leading portion of an extended school and avoids
the previously rejected source-COM entry surface. Apart from a `1e-14`
floating-point zero tolerance, no storage, distance, current, or performance
threshold is selected.

The primary extraction contains 2,338 complete excursions from all 40 seeds,
with 2,298 successive-attempt transitions. Excursions touching the beginning
or end of an archived trajectory are excluded because their recovery is
censored. Both approach sides are occupied at onset in only 116 attempts; the
larger storage gives an unambiguous direction in every case. The initiating
side is the majority chamber in 1,811 attempts; 24 onsets have equal chamber
populations. Minority-side attempts are retained rather than silently
relabelled using their future outcome.

At 1.0-time resolution, 2,194 attempts are recovered instead of 2,338. The
main conclusions persist: lag-one waiting correlation is 0.350 rather than
0.407, and the excess direction persistence is 0.026 rather than 0.023. The
0.5 grid therefore affects boundary timing and merges some short empty gaps,
but it does not create the renewal failure.

## Stages and exact population reward

Within each complete excursion, centerline counts from the existing population
cache define:

- **first entry:** at least one crossing in the initiating direction;
- **return:** first entry occurs, but the net direction-aligned population
  change at reset is nonpositive;
- **partial transfer:** net direction-aligned change is 1--4 fish;
- **sustained transfer:** net change is at least 5 fish.

The five-fish boundary is the established “first plus at least four
additional” scale. Here it is evaluated at the physical reset, not at an
arbitrary new future horizon. Thus it is an attempt-complete transfer mark and
is not identical to the previous fixed five-time continuation label.

The four mutually exclusive outcome states are consequently

```text
no_entry, returned, partial, sustained.
```

Across all attempts their counts are `1640, 114, 306, 278`. The identity

```text
R_k = direction_k (n_R(reset)-n_R(pre))
    = forward crossings - backward crossings
```

holds exactly for every attempt; the maximum residual is zero fish.

This population identity connects directly to the exact weighted-storage
law. Since the approach storage is zero at both endpoints, its integrated
storage-change term vanishes over a complete excursion. Net neck transfer is
therefore exactly balanced by integrated longitudinal supply, weighted lateral
exchange, and the already audited finite-step path correction. Conservation
constrains each reward `R_k`; it does not determine the stochastic distribution
of supply and lateral escape that generates that reward.

## Renewal and semi-Markov formulations

Let attempt `k` have outcome `M_k`, signed initiating direction `D_k`, reward
`R_k`, and recovery wait `W_k` from reset to the next detected onset.

An ordinary marked renewal description requires successive tuples to be iid:

```text
(W_k, D_{k+1}, M_{k+1}, R_{k+1}) independent of earlier attempts.
```

It need not be Poisson; an arbitrary nonexponential waiting distribution is
allowed. Consequently writing

```text
transfer rate = attempt frequency x mean reward per attempt
```

is only a renewal-reward identity if this independence is established.

The smallest tested alternative is the four-state Markov-renewal kernel

```text
K_rs(t,r,d) = P(W_k <= t, M_{k+1}=s,
                R_{k+1}=r, D_{k+1}/D_k=d | M_k=r_state).
```

No microscopic coordinate or fitted drift is added. Training seeds estimate
only empirical outcome-conditioned probabilities and means. The complete
four-by-four count/probability table is provided for transparency; held-out
predictions use the previous outcome only.

## Ordinary renewal is rejected

Waiting times are strongly nonexponential. On held-out seeds their mean is
7.82, median 9.5, and coefficient of variation 0.468. The exponential
distribution with the same mean has Kolmogorov distance 0.357. The empirical
survival curve has a strong concentration near 9--11 time units, consistent
with a circulation/return time rather than a constant attempt hazard.

Nonexponential waiting alone would not reject renewal. Serial dependence does:

| held-out diagnostic | estimate | 95% seed interval |
|---|---:|---:|
| lag-one waiting-time Spearman correlation | 0.386 | 0.299, 0.461 |
| lag-one transfer Spearman correlation | 0.059 | 0.008, 0.109 |
| same-direction excess over within-seed iid marks | 0.059 | 0.015, 0.101 |

The same waiting-memory result is present in training seeds (`0.392`, interval
`0.300--0.454`) and remains at the coarser extraction resolution. An approach
layer becoming empty is therefore a geometric reset, but not a statistical
reset of the circulating, spatially extended school.

## Minimal semi-Markov test

Previous outcome has a reproducible effect on recovery. Mean recovery waits in
training/held-out data are:

| previous outcome | training | held out |
|---|---:|---:|
| no entry | 9.74 | 9.06 |
| returned | 4.74 | 5.39 |
| partial | 5.14 | 4.83 |
| sustained | 6.95 | 7.73 |

An encounter with no crossing is followed by nearly a full circulation before
the next approach. Entry and partial transfer often leave fish near the neck
or activate the opposite-side layer, shortening recovery and reducing
directional persistence.

On seeds 30--39, previous-outcome conditioning changes prediction as follows:

| target | ordinary error | semi-Markov error | improvement, 95% seed interval |
|---|---:|---:|---:|
| recovery wait MSE | 13.872 | 10.655 | +3.217 `[2.051,4.487]` |
| next transfer MSE | 64.037 | 63.538 | +0.499 `[0.089,0.919]` |
| next first-entry Brier | 0.2707 | 0.2627 | +0.0080 `[-0.0018,0.0204]` |
| next partial-transfer Brier | 0.1482 | 0.1436 | +0.0046 `[-0.0031,0.0124]` |
| next sustained-transfer Brier | 0.1522 | 0.1517 | +0.00046 `[-0.00128,0.00218]` |
| same-direction Brier | 0.2410 | 0.2232 | +0.0178 `[0.0085,0.0265]` |

Thus the minimal history state genuinely improves encounter timing and
direction. Its improvement for transfer amount is statistically resolved but
less than 1% of the ordinary transfer MSE. Entry, partial-transfer, and
sustained-transfer probability improvements are unresolved.

Most importantly, the held-out lag-one correlation of the outcome-conditioned
waiting residual is still `0.140 [0.069,0.209]`. The four-state kernel has not
produced a semi-Markov closure. Adding more empirical states to absorb that
memory would violate the purpose of this analysis, so the state hierarchy is
not expanded.

## Independent-seed population-transfer prediction

An ordinary renewal model trained on seeds 0--29 gives a reasonable attempt
frequency but fails to transfer the mark distribution to seeds 30--39:

| equal-seed quantity | training prediction | held-out observation | held-out minus train, 95% interval |
|---|---:|---:|---:|
| attempts per 100 time | 7.18 | 7.68 | +0.49 `[-0.04,0.97]` |
| first-entry probability | 0.228 | 0.415 | +0.187 `[0.071,0.299]` |
| partial-transfer probability | 0.101 | 0.172 | +0.071 `[0.004,0.135]` |
| sustained-transfer probability | 0.088 | 0.176 | +0.088 `[0.039,0.136]` |
| mean direction-aligned fish transferred | 1.30 | 2.65 | +1.35 `[0.44,2.27]` |
| mean gross crossings | 2.00 | 4.06 | +2.06 `[0.89,3.20]` |

This is not repaired by the previous-outcome model: its held-out sustained
probability gain is unresolved. The data therefore support a relatively
repeatable encounter clock but not a seed-independent per-attempt commitment
law. Unequal event counts are not treated as independent evidence; all quoted
intervals resample original seeds.

## Partial transfer is not complete relocation

The existing opposite-basin definition—successive visits to opposite
`|q|>=0.6` basins—finds only nine relocations in six seeds. Their median
duration is 97 time units, and they overlap a median of 10 complete approach
attempts (range 2--45). Every relocation center and target-basin arrival occurs
during an approach excursion, but usually not the same excursion.

The attempt model therefore addresses first entry and partial/sustained
short-transfer rewards. A complete relocation is an accumulated, multi-attempt
population process. Nine events cannot validate a first-passage law for that
accumulation, and a five-fish attempt must not be called a whole-school switch.

## Physical interpretation and route to theory

The supported episode-level mechanism is

```text
coherent chamber circulation
  -> leading-edge loading of the approach layer
  -> entry or retreat
  -> conservation-controlled partial transfer
  -> geometric emptying, but incomplete dynamical reset.
```

The following pieces are exact or potentially derivable from first principles:

- the approach geometry and weighted storage;
- the endpoint population reward and integrated storage/flux identity;
- coherent aggregate translation during an encounter;
- in principle, the circulation return time from chamber geometry, swimming
  speed, wall turning, angular diffusion, and the production interaction law.

The following remain empirical:

- the nonexponential return-time distribution;
- its serial phase memory after an empty-layer reset;
- the four-state outcome kernel and transfer-reward distribution;
- the probability that an approach sustains five-fish transfer;
- accumulation of many rewards into an opposite-basin relocation.

The previous directed-alignment analysis explains why the commitment mark is
not yet derivable: relative headings and the visual graph retain memory over an
attempt. The present result adds that the chamber-scale circulation phase also
survives a geometric retreat. A defensible analytical theory should therefore
derive a phase/age-dependent return kernel for the coherent extended school
and couple its attempt rewards to exact storage conservation. It should not
assume Poisson attempts or fit another instantaneous potential.

At present, an attempt-level semi-Markov notation is useful bookkeeping, but a
predictive collective transition-rate theory has not been established.

## Reproduction

From `results/nohydro_gap20_descriptive_transition/`:

```bash
uv run --with numpy --with scipy --with matplotlib \
  python analyze_approach_renewal.py --self-test

uv run --with numpy --with scipy --with matplotlib \
  python analyze_approach_renewal.py \
  --features outputs/initiation_propagation/cache/population_kinetics_features_dt0p5.npz \
  --q-cache-dir outputs/initiation_propagation/cache \
  --output outputs/approach_renewal
```

The feature cache and population indices are compact read-only representations
of the existing 40 archived trajectories. No microscopic simulation, fitted
potential, additional collective coordinate, or parameter sweep is used.
