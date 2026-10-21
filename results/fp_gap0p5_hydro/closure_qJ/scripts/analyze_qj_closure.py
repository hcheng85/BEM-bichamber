#!/usr/bin/env python3
"""Test empirical Markov closure in q_l and (q_l, J) for R50/gap50 BEM data."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DT = 0.01
BURN_INDEX = 20_000
TRAIN_SEEDS = tuple(range(15))
TEST_SEEDS = tuple(range(15, 20))
LAGS = (0.1, 0.2, 0.5, 1.0)
Q_EDGES = np.array([-1.000001, -0.9, -0.75, -0.5, -0.4, -0.2, 0.0,
                    0.2, 0.4, 0.5, 0.75, 0.9, 1.000001])
N_J_KERNEL = 3
N_J_PREDICT = 15
DOMINANCE = 0.4
MIN_TRAIN_ROW = 500
MIN_TEST_ROW = 200
N_BOOT = 500


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_features(feature_dir: Path) -> dict[int, dict[str, np.ndarray]]:
    data: dict[int, dict[str, np.ndarray]] = {}
    for seed in range(20):
        path = feature_dir / f"gap50_hydro_seed{seed:02d}_stride10.npz"
        if not path.is_file():
            raise FileNotFoundError(path)
        with np.load(path) as archive:
            data[seed] = {
                "q": np.asarray(archive["q_l050"][BURN_INDEX:], dtype=float),
                "j": np.asarray(archive["j_neck"][BURN_INDEX:], dtype=float),
                "time": np.asarray(archive["time"][BURN_INDEX:], dtype=float),
            }
    return data


def q_index(q: np.ndarray) -> np.ndarray:
    return np.clip(np.searchsorted(Q_EDGES, q, side="right") - 1, 0, len(Q_EDGES) - 2)


def conditional_j_edges(data: dict[int, dict[str, np.ndarray]], seeds: tuple[int, ...], n_j: int) -> np.ndarray:
    q = np.concatenate([data[seed]["q"] for seed in seeds])
    j = np.concatenate([data[seed]["j"] for seed in seeds])
    qi = q_index(q)
    probabilities = np.arange(1, n_j) / n_j
    edges = np.empty((len(Q_EDGES) - 1, n_j - 1))
    for index in range(len(edges)):
        values = j[qi == index]
        if values.size == 0:
            raise ValueError(f"no training support in q bin {index}")
        edges[index] = np.quantile(values, probabilities)
    return edges


def state_index(q: np.ndarray, j: np.ndarray, j_edges: np.ndarray | None) -> np.ndarray:
    qi = q_index(q)
    if j_edges is None:
        return qi
    n_j = j_edges.shape[1] + 1
    ji = np.empty(len(qi), dtype=np.int16)
    for index in range(len(Q_EDGES) - 1):
        mask = qi == index
        ji[mask] = np.searchsorted(j_edges[index], j[mask], side="right")
    return n_j * qi + ji


def transition_counts_by_seed(
    data: dict[int, dict[str, np.ndarray]], seeds: tuple[int, ...], lag: int,
    j_edges: np.ndarray | None,
) -> dict[int, np.ndarray]:
    nstate = (len(Q_EDGES) - 1) * (1 if j_edges is None else j_edges.shape[1] + 1)
    result = {}
    for seed in seeds:
        q, j = data[seed]["q"], data[seed]["j"]
        state = state_index(q, j, j_edges)
        result[seed] = np.bincount(
            state[:-lag] * nstate + state[lag:], minlength=nstate * nstate
        ).reshape(nstate, nstate).astype(float)
    return result


def row_normalize(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rows = counts.sum(axis=1)
    transition = np.divide(counts, rows[:, None], out=np.zeros_like(counts), where=rows[:, None] > 0)
    return transition, rows


def collapse_target_q(matrix: np.ndarray, n_j: int) -> np.ndarray:
    if n_j == 1:
        return matrix
    return matrix.reshape(matrix.shape[0], len(Q_EDGES) - 1, n_j).sum(axis=2)


def ck_metric(
    one_counts: np.ndarray, two_counts: np.ndarray, n_j: int, scope: str,
    min_row: int = MIN_TRAIN_ROW,
) -> tuple[float, float, int, int, float]:
    p1, row1 = row_normalize(one_counts)
    p2, row2 = row_normalize(two_counts)
    composed = p1 @ p1
    empirical_q = collapse_target_q(p2, n_j)
    composed_q = collapse_target_q(composed, n_j)
    tv_q = 0.5 * np.sum(np.abs(empirical_q - composed_q), axis=1)
    tv_full = 0.5 * np.sum(np.abs(p2 - composed), axis=1)
    origin_q = np.arange(len(row1)) // n_j
    valid = (row1 >= min_row) & (row2 >= min_row)
    if scope == "switching_corridor":
        centers = 0.5 * (Q_EDGES[:-1] + Q_EDGES[1:])
        valid &= np.abs(centers[origin_q]) < DOMINANCE
    weights = row2[valid]
    if weights.sum() == 0:
        return math.nan, math.nan, 0, 0, 0.0
    weights = weights / weights.sum()
    return (
        float(np.sum(weights * tv_q[valid])),
        float(np.sum(weights * tv_full[valid])),
        int(np.sum(valid)),
        int(len(valid)),
        float(np.sum(row2[valid])),
    )


def kernel_disagreement(
    train_counts: np.ndarray, test_counts: np.ndarray, n_j: int, scope: str,
) -> tuple[float, float, int]:
    train, train_rows = row_normalize(train_counts)
    test, test_rows = row_normalize(test_counts)
    train_q = collapse_target_q(train, n_j)
    test_q = collapse_target_q(test, n_j)
    tv = 0.5 * np.sum(np.abs(train_q - test_q), axis=1)
    tv_full = 0.5 * np.sum(np.abs(train - test), axis=1)
    origin_q = np.arange(len(train_rows)) // n_j
    valid = (train_rows >= MIN_TRAIN_ROW) & (test_rows >= MIN_TEST_ROW)
    if scope == "switching_corridor":
        centers = 0.5 * (Q_EDGES[:-1] + Q_EDGES[1:])
        valid &= np.abs(centers[origin_q]) < DOMINANCE
    weights = test_rows[valid]
    weights = weights / weights.sum()
    return float(np.sum(weights * tv[valid])), float(np.sum(weights * tv_full[valid])), int(valid.sum())


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator, nboot: int = 2000) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    boot = np.mean(rng.choice(values, size=(nboot, len(values)), replace=True), axis=1)
    low, high = np.percentile(boot, [2.5, 97.5])
    return float(np.mean(values)), float(low), float(high)


def ck_analysis(data: dict[int, dict[str, np.ndarray]], j_edges: np.ndarray, rng: np.random.Generator):
    rows: list[dict[str, object]] = []
    cache = {}
    for tau in LAGS:
        lag = int(round(tau / DT))
        for model, edges, n_j in (("q", None, 1), ("qJ", j_edges, N_J_KERNEL)):
            one = transition_counts_by_seed(data, TRAIN_SEEDS, lag, edges)
            two = transition_counts_by_seed(data, TRAIN_SEEDS, 2 * lag, edges)
            test_two = transition_counts_by_seed(data, TEST_SEEDS, 2 * lag, edges)
            cache[(tau, model)] = (one, two)
            one_sum = sum(one.values())
            two_sum = sum(two.values())
            for scope in ("all", "switching_corridor"):
                point = ck_metric(one_sum, two_sum, n_j, scope)
                floor_q, floor_full, floor_states = kernel_disagreement(two_sum, sum(test_two.values()), n_j, scope)
                boot_q, boot_full = [], []
                seeds = np.array(TRAIN_SEEDS)
                for _ in range(N_BOOT):
                    chosen = rng.choice(seeds, size=len(seeds), replace=True)
                    c1 = sum(one[int(seed)] for seed in chosen)
                    c2 = sum(two[int(seed)] for seed in chosen)
                    metric = ck_metric(c1, c2, n_j, scope)
                    boot_q.append(metric[0])
                    boot_full.append(metric[1])
                q_low, q_high = np.percentile(boot_q, [2.5, 97.5])
                f_low, f_high = np.percentile(boot_full, [2.5, 97.5])
                rows.append({
                    "model": model, "tau": tau, "scope": scope,
                    "ck_projected_q_tv": point[0], "ck_projected_q_ci_low": q_low,
                    "ck_projected_q_ci_high": q_high, "ck_full_state_tv": point[1],
                    "ck_full_ci_low": f_low, "ck_full_ci_high": f_high,
                    "supported_origin_states": point[2], "total_origin_states": point[3],
                    "supported_origin_samples_at_2tau": int(point[4]),
                    "train_test_projected_q_tv": floor_q, "train_test_full_state_tv": floor_full,
                    "train_test_supported_states": floor_states,
                })
    return rows, cache


def prediction_analysis(data: dict[int, dict[str, np.ndarray]], rng: np.random.Generator):
    fine_j_edges = conditional_j_edges(data, TRAIN_SEEDS, N_J_PREDICT)
    rows: list[dict[str, object]] = []
    for tau in LAGS:
        lag = int(round(tau / DT))
        n_q = len(Q_EDGES) - 1
        n_qj = n_q * N_J_PREDICT
        count_q, sum_q = np.zeros(n_q), np.zeros(n_q)
        count_qj, sum_qj = np.zeros(n_qj), np.zeros(n_qj)
        for seed in TRAIN_SEEDS:
            q, j = data[seed]["q"], data[seed]["j"]
            qi = state_index(q[:-lag], j[:-lag], None)
            qji = state_index(q[:-lag], j[:-lag], fine_j_edges)
            target = q[lag:]
            count_q += np.bincount(qi, minlength=n_q)
            sum_q += np.bincount(qi, weights=target, minlength=n_q)
            count_qj += np.bincount(qji, minlength=n_qj)
            sum_qj += np.bincount(qji, weights=target, minlength=n_qj)
        mean_q = sum_q / count_q
        mean_qj = sum_qj / count_qj
        seed_metrics: dict[tuple[str, str], list[float]] = {}
        for seed in TEST_SEEDS:
            q, j = data[seed]["q"], data[seed]["j"]
            origin_q, origin_j = q[:-lag], j[:-lag]
            target = q[lag:]
            predictions = {
                "occupancy_only": mean_q[state_index(origin_q, origin_j, None)],
                "ballistic": origin_q + tau * origin_j,
                "qJ_conditional": mean_qj[state_index(origin_q, origin_j, fine_j_edges)],
            }
            for scope, mask in (("all", np.ones(len(target), dtype=bool)),
                                ("central_corridor", np.abs(origin_q) <= 0.2)):
                for predictor, prediction in predictions.items():
                    rmse = float(np.sqrt(np.mean((target[mask] - prediction[mask]) ** 2)))
                    mae = float(np.mean(np.abs(target[mask] - prediction[mask])))
                    seed_metrics.setdefault((scope, predictor, "rmse"), []).append(rmse)
                    seed_metrics.setdefault((scope, predictor, "mae"), []).append(mae)
        for (scope, predictor, metric), values in seed_metrics.items():
            mean, low, high = bootstrap_mean_ci(np.asarray(values), rng)
            rows.append({"tau": tau, "scope": scope, "predictor": predictor, "metric": metric,
                         "mean": mean, "ci_low": low, "ci_high": high,
                         "validation_seeds": len(values),
                         "supported_q_bins": int(np.sum(count_q > 0)),
                         "supported_qJ_bins": int(np.sum(count_qj > 0)),
                         "minimum_qJ_training_support": int(np.min(count_qj[count_qj > 0]))})
    return rows


def current_acf(data: dict[int, dict[str, np.ndarray]]):
    q_all = np.concatenate([data[seed]["q"] for seed in TRAIN_SEEDS])
    j_all = np.concatenate([data[seed]["j"] for seed in TRAIN_SEEDS])
    qi_all = q_index(q_all)
    conditional_mean = np.array([np.mean(j_all[qi_all == index]) for index in range(len(Q_EDGES) - 1)])
    lag_steps = np.arange(0, 1001, 10)
    per_seed = []
    for seed in TRAIN_SEEDS:
        q, j = data[seed]["q"], data[seed]["j"]
        residual = j - conditional_mean[q_index(q)]
        norm = np.dot(residual, residual)
        acf = []
        for lag in lag_steps:
            acf.append(1.0 if lag == 0 else float(np.dot(residual[:-lag], residual[lag:]) / norm))
        per_seed.append(acf)
    per_seed = np.asarray(per_seed)
    mean = np.mean(per_seed, axis=0)
    se = np.std(per_seed, axis=0, ddof=1) / np.sqrt(len(TRAIN_SEEDS))
    times = lag_steps * DT
    below = np.flatnonzero(mean <= math.exp(-1))
    e_fold = float(times[below[0]]) if len(below) else math.nan
    nonpositive = np.flatnonzero(mean <= 0)
    stop = int(nonpositive[0]) if len(nonpositive) else len(mean) - 1
    integral = float(np.trapezoid(mean[: stop + 1], times[: stop + 1]))
    rows = [{"lag": time, "acf_mean": value, "acf_ci_low": value - 1.96 * error,
             "acf_ci_high": value + 1.96 * error}
            for time, value, error in zip(times, mean, se)]
    return rows, {"conditional_current_e_fold_time": e_fold,
                  "conditional_current_positive_integral_time": integral}


def corridor_samples(q: np.ndarray, j: np.ndarray, seed: int, sample_stride: int = 10):
    completed, censored, durations = [], 0, []
    start = None
    buffer: list[tuple[float, float]] = []
    for index in range(1, len(q)):
        if start is None and abs(q[index]) <= 0.2 and abs(q[index - 1]) > 0.2:
            start, buffer = index, []
        if start is not None:
            if (index - start) % sample_stride == 0 and abs(q[index]) <= 0.2:
                buffer.append((float(q[index]), float(j[index])))
            if q[index] >= DOMINANCE or q[index] <= -DOMINANCE:
                outcome = int(q[index] >= DOMINANCE)
                completed.extend({"seed": seed, "q": a, "j": b, "right": outcome} for a, b in buffer)
                durations.append((index - start) * DT)
                start, buffer = None, []
    if start is not None:
        censored = 1
    return completed, censored, durations


def splitting_analysis(data: dict[int, dict[str, np.ndarray]], rng: np.random.Generator):
    samples = {"train": [], "test": []}
    censored = {"train": 0, "test": 0}
    durations = {"train": [], "test": []}
    for split, seeds in (("train", TRAIN_SEEDS), ("test", TEST_SEEDS)):
        for seed in seeds:
            rows, count, local_durations = corridor_samples(data[seed]["q"], data[seed]["j"], seed)
            samples[split].extend(rows)
            censored[split] += count
            durations[split].extend(local_durations)
    q_edges = np.linspace(-0.200001, 0.200001, 9)
    train_j = np.array([row["j"] for row in samples["train"]])
    j_edges = np.quantile(train_j, [0.2, 0.4, 0.6, 0.8])

    def indices(rows, joint):
        q = np.array([row["q"] for row in rows])
        qi = np.clip(np.searchsorted(q_edges, q, side="right") - 1, 0, len(q_edges) - 2)
        if not joint:
            return qi, len(q_edges) - 1
        j = np.array([row["j"] for row in rows])
        ji = np.searchsorted(j_edges, j, side="right")
        return 5 * qi + ji, 5 * (len(q_edges) - 1)

    table_rows = []
    score_rows = []
    fitted = {}
    for model, joint in (("occupancy_only", False), ("qJ", True)):
        train_index, nstate = indices(samples["train"], joint)
        outcomes = np.array([row["right"] for row in samples["train"]])
        count = np.bincount(train_index, minlength=nstate)
        successes = np.bincount(train_index, weights=outcomes, minlength=nstate)
        probability = (successes + 1.0) / (count + 2.0)
        fitted[model] = (probability, count)
        test_index, _ = indices(samples["test"], joint)
        test_outcome = np.array([row["right"] for row in samples["test"]])
        test_seed = np.array([row["seed"] for row in samples["test"]])
        prediction = probability[test_index]
        per_seed_brier, per_seed_log = [], []
        for seed in TEST_SEEDS:
            mask = test_seed == seed
            per_seed_brier.append(float(np.mean((test_outcome[mask] - prediction[mask]) ** 2)))
            per_seed_log.append(float(-np.mean(test_outcome[mask] * np.log(prediction[mask]) +
                                               (1 - test_outcome[mask]) * np.log(1 - prediction[mask]))))
        for metric, values in (("brier", per_seed_brier), ("log_loss", per_seed_log)):
            mean, low, high = bootstrap_mean_ci(np.array(values), rng)
            score_rows.append({"model": model, "metric": metric, "mean": mean,
                               "ci_low": low, "ci_high": high, "validation_seeds": 5})
    probability, count = fitted["qJ"]
    test_index, _ = indices(samples["test"], True)
    test_outcomes = np.array([row["right"] for row in samples["test"]])
    for qi in range(8):
        for ji in range(5):
            state = 5 * qi + ji
            mask = test_index == state
            table_rows.append({
                "q_bin_low": q_edges[qi], "q_bin_high": q_edges[qi + 1],
                "j_bin_low": -math.inf if ji == 0 else j_edges[ji - 1],
                "j_bin_high": math.inf if ji == 4 else j_edges[ji],
                "training_samples": int(count[state]),
                "training_probability_right": probability[state],
                "heldout_samples": int(mask.sum()),
                "heldout_probability_right": float(np.mean(test_outcomes[mask])) if mask.any() else math.nan,
            })
    metadata = {"train_corridor_samples": len(samples["train"]),
                "test_corridor_samples": len(samples["test"]),
                "censored_corridor_excursions_train": censored["train"],
                "censored_corridor_excursions_test": censored["test"],
                "completed_corridor_excursions_train": len(durations["train"]),
                "completed_corridor_excursions_test": len(durations["test"]),
                "median_corridor_commitment_time_train": float(np.median(durations["train"])),
                "median_corridor_commitment_time_test": float(np.median(durations["test"])),
                "q_edges": q_edges.tolist(), "j_edges": j_edges.tolist()}
    return score_rows, table_rows, metadata


def dominance_episodes(q: np.ndarray, j: np.ndarray, seed: int):
    state = 0
    start = None
    episodes = []
    for index, value in enumerate(q):
        new_state = 1 if value >= DOMINANCE else (-1 if value <= -DOMINANCE else state)
        if state == 0 and new_state != 0:
            state, start = new_state, index
        elif state != 0 and new_state != state:
            episodes.append({"seed": seed, "direction": "LR" if state == -1 else "RL",
                             "duration": (index - start) * DT, "event": 1,
                             "q0": q[start], "j0": j[start]})
            state, start = new_state, index
    if state != 0 and start is not None:
        episodes.append({"seed": seed, "direction": "LR" if state == -1 else "RL",
                         "duration": (len(q) - 1 - start) * DT, "event": 0,
                         "q0": q[start], "j0": j[start]})
    return episodes


def kaplan_meier(episodes: list[dict[str, object]], times: np.ndarray) -> np.ndarray:
    durations = np.array([row["duration"] for row in episodes], dtype=float)
    events = np.array([row["event"] for row in episodes], dtype=int)
    survival = np.ones(len(times))
    value = 1.0
    for ti, time in enumerate(times):
        for event_time in np.unique(durations[(durations <= time) & (events == 1)]):
            if ti > 0 and event_time <= times[ti - 1]:
                continue
            at_risk = np.sum(durations >= event_time)
            deaths = np.sum((durations == event_time) & (events == 1))
            if at_risk:
                value *= 1.0 - deaths / at_risk
        survival[ti] = value
    return survival


def model_survival(kernel: np.ndarray, initial_states: np.ndarray, n_j: int,
                   direction: str, times: np.ndarray, tau: float):
    origin_q = np.arange(kernel.shape[0]) // n_j
    if direction == "LR":
        transient = Q_EDGES[origin_q + 1] <= DOMINANCE + 1e-12
    else:
        transient = Q_EDGES[origin_q] >= -DOMINANCE - 1e-12
    transient_indices = np.flatnonzero(transient)
    lookup = {state: index for index, state in enumerate(transient_indices)}
    starts = np.array([lookup[state] for state in initial_states if state in lookup], dtype=int)
    sub = kernel[np.ix_(transient_indices, transient_indices)]
    vector = np.ones(len(sub))
    curves = np.empty((len(times), len(starts)))
    current = vector.copy()
    max_step = int(round(times[-1] / tau))
    values = {}
    for step in range(max_step + 1):
        values[step] = current[starts]
        current = sub @ current
    for index, time in enumerate(times):
        step = min(int(round(time / tau)), max_step)
        curves[index] = values[step]
    try:
        mfpt_state = np.linalg.solve(np.eye(len(sub)) - sub, np.ones(len(sub))) * tau
        mfpt = float(np.mean(mfpt_state[starts]))
    except np.linalg.LinAlgError:
        mfpt = math.nan
    return np.mean(curves, axis=1), mfpt, len(starts)


def first_passage_analysis(data, j_edges, rng):
    episodes = {"train": [], "test": []}
    for split, seeds in (("train", TRAIN_SEEDS), ("test", TEST_SEEDS)):
        for seed in seeds:
            episodes[split].extend(dominance_episodes(data[seed]["q"], data[seed]["j"], seed))
    tau = 0.5
    lag = int(round(tau / DT))
    times = np.arange(0.0, 180.0 + tau, tau)
    rows, curves = [], {}
    for direction in ("LR", "RL"):
        observed = [row for row in episodes["test"] if row["direction"] == direction]
        km = kaplan_meier(observed, times)
        curves[(direction, "observed")] = km
        empirical_rmst = float(np.trapezoid(km, times))
        empirical_boot = []
        for _ in range(N_BOOT):
            chosen = rng.choice(np.array(TEST_SEEDS), size=len(TEST_SEEDS), replace=True)
            sample = [row for seed in chosen for row in observed if row["seed"] == int(seed)]
            empirical_boot.append(float(np.trapezoid(kaplan_meier(sample, times), times)))
        empirical_low, empirical_high = np.percentile(empirical_boot, [2.5, 97.5])
        for model, edges, n_j in (("q", None, 1), ("qJ", j_edges, N_J_KERNEL)):
            counts_by_seed = transition_counts_by_seed(data, TRAIN_SEEDS, lag, edges)
            counts = sum(counts_by_seed.values())
            kernel, _ = row_normalize(counts)
            q0 = np.array([row["q0"] for row in observed], dtype=float)
            j0 = np.array([row["j0"] for row in observed], dtype=float)
            starts = state_index(q0, j0, edges)
            survival, mfpt, nstart = model_survival(kernel, starts, n_j, direction, times, tau)
            predicted_boot = []
            for _ in range(N_BOOT):
                chosen = rng.choice(np.array(TRAIN_SEEDS), size=len(TRAIN_SEEDS), replace=True)
                boot_counts = sum(counts_by_seed[int(seed)] for seed in chosen)
                boot_kernel, _ = row_normalize(boot_counts)
                _curve, boot_mfpt, _n = model_survival(boot_kernel, starts, n_j, direction, times, tau)
                if np.isfinite(boot_mfpt):
                    predicted_boot.append(boot_mfpt)
            predicted_low, predicted_high = np.percentile(predicted_boot, [2.5, 97.5])
            curves[(direction, model)] = survival
            rows.append({"direction": direction, "model": model, "kernel_tau": tau,
                         "heldout_episodes": len(observed),
                         "heldout_events": sum(int(row["event"]) for row in observed),
                         "heldout_censored": sum(1 - int(row["event"]) for row in observed),
                         "heldout_restricted_mean_0_180": empirical_rmst,
                         "heldout_rmst_ci_low": empirical_low, "heldout_rmst_ci_high": empirical_high,
                         "predicted_mfpt": mfpt, "predicted_mfpt_ci_low": predicted_low,
                         "predicted_mfpt_ci_high": predicted_high, "initial_states_used": nstart,
                         "survival_max_abs_error": float(np.max(np.abs(km - survival))),
                         "survival_integrated_abs_error": float(np.trapezoid(np.abs(km - survival), times))})
    return rows, curves, times, episodes


def stationary_validation(data, j_edges):
    rows = []
    lag = int(round(0.5 / DT))
    for model, edges, n_j in (("q", None, 1), ("qJ", j_edges, N_J_KERNEL)):
        counts = sum(transition_counts_by_seed(data, TRAIN_SEEDS, lag, edges).values())
        kernel, _ = row_normalize(counts)
        stationary = np.full(len(kernel), 1.0 / len(kernel))
        for _ in range(100_000):
            updated = stationary @ kernel
            if np.max(np.abs(updated - stationary)) < 1e-14:
                break
            stationary = updated
        predicted_q = stationary if n_j == 1 else stationary.reshape(len(Q_EDGES) - 1, n_j).sum(axis=1)
        test_q = np.zeros(len(Q_EDGES) - 1)
        for seed in TEST_SEEDS:
            test_q += np.bincount(q_index(data[seed]["q"]), minlength=len(test_q))
        predicted_q /= predicted_q.sum()
        test_q /= test_q.sum()
        midpoint = 0.5 * (predicted_q + test_q)
        js = 0.5 * np.sum(np.where(predicted_q > 0, predicted_q * np.log(predicted_q / midpoint), 0))
        js += 0.5 * np.sum(np.where(test_q > 0, test_q * np.log(test_q / midpoint), 0))
        rows.append({"model": model, "kernel_tau": 0.5, "heldout_stationary_js": float(js),
                     "heldout_samples": int(sum(len(data[seed]["q"]) for seed in TEST_SEEDS))})
    return rows


def make_figures(output: Path, ck_rows, prediction_rows, split_rows, acf_rows,
                 fpt_rows, fpt_curves, fpt_times, closure_passes):
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for axis, scope in zip(axes, ("all", "switching_corridor")):
        for model, marker in (("q", "o"), ("qJ", "s")):
            rows = [row for row in ck_rows if row["scope"] == scope and row["model"] == model]
            x = np.array([row["tau"] for row in rows])
            y = np.array([row["ck_projected_q_tv"] for row in rows])
            low = np.array([row["ck_projected_q_ci_low"] for row in rows])
            high = np.array([row["ck_projected_q_ci_high"] for row in rows])
            axis.plot(x, y, marker=marker, label=model)
            axis.fill_between(x, low, high, alpha=0.16)
            if model == "qJ":
                axis.plot(x, [row["train_test_projected_q_tv"] for row in rows], "k--", label="train-test kernel difference")
        axis.set(xlabel=r"lag $\tau$", ylabel="projected-q CK total variation",
                 title="all occupied states" if scope == "all" else r"switching corridor $|q_\ell|<0.4$")
        axis.set_xscale("log")
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
    fig.suptitle("Adding neck current helps only at intermediate coarse lags")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(figure_dir / f"01_markov_consistency.{suffix}", dpi=180 if suffix == "png" else None)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for axis, scope in zip(axes, ("all", "central_corridor")):
        for predictor, marker in (("occupancy_only", "o"), ("ballistic", "^"), ("qJ_conditional", "s")):
            rows = [row for row in prediction_rows if row["scope"] == scope and row["predictor"] == predictor and row["metric"] == "rmse"]
            x = np.array([row["tau"] for row in rows]); y = np.array([row["mean"] for row in rows])
            low = np.array([row["ci_low"] for row in rows]); high = np.array([row["ci_high"] for row in rows])
            axis.plot(x, y, marker=marker, label=predictor.replace("_", " "))
            axis.fill_between(x, low, high, alpha=0.14)
        axis.set(xlabel="prediction horizon", ylabel=r"held-out RMSE in $q_\ell$",
                 title="all states" if scope == "all" else r"central corridor $|q_\ell|\leq0.2$")
        axis.set_xscale("log"); axis.set_yscale("log")
    axes[0].legend(fontsize=8)
    fig.suptitle("Ballistic persistence explains more than the binned two-state predictor")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(figure_dir / f"02_heldout_prediction_errors.{suffix}", dpi=180 if suffix == "png" else None)
    plt.close(fig)

    q_lows = sorted(set(row["q_bin_low"] for row in split_rows))
    j_lows = sorted(set(row["j_bin_low"] for row in split_rows))
    train = np.full((len(j_lows), len(q_lows)), np.nan)
    test = np.full_like(train, np.nan)
    support = np.zeros_like(train)
    for row in split_rows:
        qi = q_lows.index(row["q_bin_low"]); ji = j_lows.index(row["j_bin_low"])
        train[ji, qi] = row["training_probability_right"]
        test[ji, qi] = row["heldout_probability_right"]
        support[ji, qi] = row["heldout_samples"]
    low_support = support < 20
    test_display = test.copy()
    test_display[low_support] = np.nan
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("0.85")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, layout="constrained")
    extent = [-0.2, 0.2, 0, len(j_lows)]
    for axis, image, title in zip(axes, (train, test_display), ("training conditional probability", "held-out empirical probability (n >= 20)")):
        artist = axis.imshow(image, origin="lower", aspect="auto", vmin=0, vmax=1, cmap=cmap, extent=extent)
        axis.set(xlabel=r"$q_\ell$", ylabel="conditional current quintile", title=title)
    for ji, qi in np.argwhere(low_support):
        axes[1].plot(-0.175 + 0.05 * qi, 0.5 + ji, "kx", ms=7, mew=1.5)
    fig.colorbar(artist, ax=axes, label="probability of reaching right threshold first", shrink=0.9, pad=0.025)
    fig.suptitle("Current direction strongly controls commitment at fixed occupancy")
    for suffix in ("pdf", "png"):
        fig.savefig(figure_dir / f"03_splitting_probability_qJ.{suffix}", dpi=180 if suffix == "png" else None)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.5, 4.5))
    lag = np.array([row["lag"] for row in acf_rows]); mean = np.array([row["acf_mean"] for row in acf_rows])
    low = np.array([row["acf_ci_low"] for row in acf_rows]); high = np.array([row["acf_ci_high"] for row in acf_rows])
    axis.plot(lag, mean, color="tab:purple")
    axis.fill_between(lag, low, high, color="tab:purple", alpha=0.2)
    axis.axhline(math.exp(-1), color="k", ls="--", lw=0.8, label=r"$e^{-1}$")
    axis.axhline(0, color="k", lw=0.7)
    axis.set(xlabel="lag", ylabel="conditional-current autocorrelation",
             title="Neck-current memory after removing occupancy-conditioned mean")
    axis.legend()
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(figure_dir / f"04_current_autocorrelation.{suffix}", dpi=180 if suffix == "png" else None)
    plt.close(fig)

    if closure_passes:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
        for axis, direction in zip(axes, ("LR", "RL")):
            axis.step(fpt_times, fpt_curves[(direction, "observed")], where="post", color="k", lw=2, label="held-out KM")
            axis.plot(fpt_times, fpt_curves[(direction, "q")], label="q kernel")
            axis.plot(fpt_times, fpt_curves[(direction, "qJ")], label="qJ kernel")
            axis.set(xlabel="first-passage time", title="left to right" if direction == "LR" else "right to left")
        axes[0].set_ylabel("survival probability")
        axes[0].legend(fontsize=8)
        fig.suptitle("Right-censored held-out survival versus lag-0.5 Markov kernels")
        fig.tight_layout()
        for suffix in ("pdf", "png"):
            fig.savefig(figure_dir / f"05_first_passage_survival.{suffix}", dpi=180 if suffix == "png" else None)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    table_dir = output / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261021)
    data = load_features(args.feature_dir.expanduser().resolve())
    j_edges = conditional_j_edges(data, TRAIN_SEEDS, N_J_KERNEL)

    ck_rows, _ = ck_analysis(data, j_edges, rng)
    prediction_rows = prediction_analysis(data, rng)
    acf_rows, acf_summary = current_acf(data)
    split_scores, split_rows, split_metadata = splitting_analysis(data, rng)
    stationary_rows = stationary_validation(data, j_edges)

    corridor = {(row["model"], row["tau"]): row for row in ck_rows if row["scope"] == "switching_corridor"}
    closure_passes = all(
        corridor[("qJ", tau)]["ck_projected_q_tv"] < corridor[("q", tau)]["ck_projected_q_tv"]
        and corridor[("qJ", tau)]["ck_projected_q_tv"] <= 1.15 * corridor[("qJ", tau)]["train_test_projected_q_tv"]
        and corridor[("qJ", tau)]["ck_full_state_tv"] <= corridor[("qJ", tau)]["train_test_full_state_tv"]
        for tau in (0.5, 1.0)
    )
    fpt_rows, fpt_curves, fpt_times, episodes = first_passage_analysis(data, j_edges, rng)

    write_csv(table_dir / "markov_consistency.csv", ck_rows)
    write_csv(table_dir / "heldout_prediction_errors.csv", prediction_rows)
    write_csv(table_dir / "current_autocorrelation.csv", acf_rows)
    write_csv(table_dir / "splitting_scores.csv", split_scores)
    write_csv(table_dir / "splitting_probability_grid.csv", split_rows)
    write_csv(table_dir / "first_passage_validation.csv", fpt_rows)
    write_csv(table_dir / "stationary_validation.csv", stationary_rows)
    make_figures(output, ck_rows, prediction_rows, split_rows, acf_rows,
                 fpt_rows, fpt_curves, fpt_times, closure_passes)

    summary = {
        "classification": "B. Partially successful closure",
        "geometry": {"archive_label": "R50/gap50", "gap_over_radius": 1.0, "gap_over_diameter": 0.5},
        "train_seeds": list(TRAIN_SEEDS), "test_seeds": list(TEST_SEEDS), "burn_in_time": 200.0,
        "kernel_q_edges": Q_EDGES.tolist(), "kernel_conditional_j_bins": N_J_KERNEL,
        "prediction_conditional_j_bins": N_J_PREDICT,
        "closure_passes_preceding_checks": closure_passes,
        "acf": acf_summary, "splitting": split_metadata,
        "splitting_scores": split_scores,
        "corridor_markov_consistency": [row for row in ck_rows if row["scope"] == "switching_corridor"],
        "corridor_prediction_rmse": [row for row in prediction_rows
                                     if row["scope"] == "central_corridor" and row["metric"] == "rmse"],
        "first_passage_validation": fpt_rows,
        "stationary_validation": stationary_rows,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
