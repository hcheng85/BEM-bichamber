# Time-sampled transition law for (population imbalance, gap-directed polarization)

Scope: post-transient raw no-hydro gap-20 trajectories, 40 independent seeds; origins separated by their horizon τ. Train seeds 0–29; held-out seeds 30–39.

## What is measured

The starting state is the EXACT chamber imbalance q (integer |N_R-N_L|, not a q bin), and the majority-school COM-to-gap heading projection a_COM. A neutral q=0 has no majority and is excluded only from the two-coordinate predictor.
For an initial majority sign s=sign(q), the next population is represented by q_next_aligned=s*q(t+τ); because s is observed this is equivalent to q(t+τ). The reported future-q kernel has 20 outcome bins spanning [-1,1]. The directional net change is Δn=-s[N_R(t+τ)-N_R(t)]; positive values move fish toward the original minority, and 7 outcome classes distinguish backward transfers, zero transfers, and different forward transfers. A positive short transfer is Δn≥5.

Training uses empirical conditional frequencies with fixed nested prior-strength smoothing κ=20 (q broad bin → exact q → exact q+a → optional q+a+history). This smoothing is a statistical regularizer, not a dynamical force or fitted physical parameter. History codes use only whether q changed and whether a increased over the preceding τ; they are diagnostic, not additional proposed landscape coordinates.

## Held-out scores

| τ | n test | transfer≥5 | q′ log loss q | q′ log loss q+a | Δn log loss q | Δn log loss q+a | five-fish Brier q | five-fish Brier q+a |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7917 | 103 | 0.2685 | 0.2557 | 0.4212 | 0.4259 | 0.0128 | 0.0129 |
| 5 | 1577 | 83 | 0.7715 | 0.7450 | 0.7970 | 0.7954 | 0.0496 | 0.0483 |
| 10 | 783 | 80 | 1.2000 | 1.1942 | 1.1571 | 1.1724 | 0.0905 | 0.0894 |

A **positive** seed-level loss difference means the more complex predictor is better.

## Independent-seed gains and confidence intervals

| τ | q+a vs q: q′ log loss | q+a vs q: Δn log loss | q+a vs q: five-fish Brier | history vs q+a: q′ log loss | history vs q+a: Δn log loss |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.012744 [0.005126, 0.019879] | -0.004731 [-0.010431, 0.000539] | -0.000064 [-0.000150, 0.000014] | 0.013990 [0.010258, 0.018525] | 0.020406 [0.012908, 0.027734] |
| 5 | 0.026464 [0.015203, 0.038506] | 0.001482 [-0.011574, 0.013743] | 0.001273 [0.000532, 0.001973] | 0.002767 [-0.003357, 0.008808] | -0.000925 [-0.005721, 0.003953] |
| 10 | 0.005789 [-0.013115, 0.024816] | -0.015389 [-0.030103, -0.000254] | 0.001072 [0.000138, 0.002021] | -0.010061 [-0.020580, 0.000852] | -0.016006 [-0.024805, -0.007242] |

## Answers to the scientific questions

### 1. Does polarization improve redistribution prediction beyond q?

For the complete binned future-population distribution, yes at τ=1 and 5: both q′ log loss and Brier score improve reproducibly. At τ=10 the q′ Brier gain remains positive, but the log-loss gain is unresolved. For the harder seven-class signed-transfer distribution, Δn log-loss gains are unresolved at τ=1 and 5 and significantly negative at τ=10. Thus polarization contains useful population-redistribution information, but not a uniformly better full transfer law at every lag.

For the threshold Δn≥5, polarization improves held-out AUC at all three horizons and Brier score at τ=5 and 10; τ=1 Brier is unresolved. Threshold prediction should not be confused with modeling the entire signed transfer distribution.

### 2. Is the improvement reproducible across independent trajectories?

The q′ improvements at τ=1 and 5 and the five-fish improvements at τ=5 and 10 have seed-bootstrap intervals excluding zero. The signed-transfer log-loss improvement does not reproduce at τ=1 or 5 and reverses at τ=10. The evidence is therefore outcome- and horizon-dependent, rather than a universal gain from adding a.

### 3. Does additional temporal history expose a limitation?

Yes. At τ=1, the predeclared past-change flags improve q′ log loss, q′ Brier, Δn log loss, Δn Brier, five-fish Brier, and AUC with independent-seed intervals excluding zero. Current (q,a) is therefore not predictively sufficient even for this limited history test. At τ=5 the log-loss gains are unresolved despite Brier improvements; at τ=10 history worsens the full Δn scores. This mixed longer-lag behavior does not repair the short-lag closure failure.

### 4. Is this sufficient for a transition-rate theory?

No. The coordinates are useful descriptive predictors, especially for binned future q and short-transfer risk, but the reproducible τ=1 history dependence rejects the tested memoryless description. The majority coordinate is also undefined at q=0 and changes reference through a reversal. Quantitative kernels remain sensitive to regularization/binning, and the held-out five-fish counts are only 103 at τ=1, 83 at τ=5, 80 at τ=10. These are short redistributions, not complete relocations.

## Conservation, parity, and sample audit

- τ=1: 7917 held-out origins, 691 nonzero transfers, 103 with Δn≥5, 263 q=0 origins excluded across all seeds, 63 future sign reversals, maximum conservation residual 1.110e-16.
- τ=5: 1577 held-out origins, 295 nonzero transfers, 83 with Δn≥5, 52 q=0 origins excluded across all seeds, 80 future sign reversals, maximum conservation residual 1.110e-16.
- τ=10: 783 held-out origins, 247 nonzero transfers, 80 with Δn≥5, 26 q=0 origins excluded across all seeds, 77 future sign reversals, maximum conservation residual 1.110e-16.

All 5,792 archived COM-polarization origins are reproduced by the raw-state parity audit; exact values and errors are in summary.json.

## Interpretation and limits

This is a **necessary predictive-sufficiency test**, not proof that (q,a) is a Markov process. The q′ conditional alone does not characterize the full joint (q′,a′) dynamics. History is tested only through predeclared signs and a q-change indicator; an unresolved history effect cannot establish Markovianity. Conversely, a reproducible held-out history improvement is evidence against the tested state being sufficient at that lag.
The data contain very few complete opposite-basin relocations; a five-fish directional change is NOT a full school relocation. Strong finite-time chamber trapping and lack of confirmed mixing prohibit interpreting -log p(q,a) as a rate barrier.
As the majority definition flips when q changes sign (and is undefined at q=0), the present two-coordinate description is not globally smooth through a complete relocation. Results may depend on smoothing strength and bin size; report sensitivity before using these empirical kernels in a rate calculation.

## Reproducibility

- summary.json and heldout_tau*.json: outcomes and seed-cluster intervals
- heldout_seed_tau*.csv: per-seed proper scores
- conditional_transfer_tau*.csv: pooled counts, seed-bootstrap intervals, mean Δn
- qnext_kernel_tau*.csv: broad-cell empirical binned future-population kernel
- heldout_transition_scores.png and conditional_transfer_tau5.png: diagnostics
- run configuration: analysis_config.json
