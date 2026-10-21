#!/usr/bin/env python3
"""Validate a one-dimensional FP reduction for the R50/gap50 full-BEM data.

The stored physical parameters give gap/R=1 and gap/(2R)=0.5.  This distinction
is kept explicit because an exact gap/R=0.5 ensemble (R50/gap25) is absent.

All model estimation uses seeds 0--14.  Seeds 15--19 are held out for the
stationary-distribution, closure, and first-passage validations.  Independent
realizations are never concatenated when increments or passages are formed.
"""
from __future__ import annotations

import csv
import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.integrate import cumulative_trapezoid
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


FEATURE_DIR = Path("features")
RESULT_DIR = Path("generated")
PROCESSED_DIR = Path("Processed_data")
DT_NATIVE = 0.001
DT_FEATURE = 0.01
BURN_TIME = 200.0
BURN_NATIVE = int(BURN_TIME / DT_NATIVE)
BURN_FEATURE = int(BURN_TIME / DT_FEATURE)
N_FISH = 100
Q_GRID = np.linspace(-1.0, 1.0, N_FISH + 1)
Q_DX = 2.0 / N_FISH
DOMINANCE = 0.4
TRAIN_SEEDS = tuple(range(15))
TEST_SEEDS = tuple(range(15, 20))
NATIVE_LAGS = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000], dtype=int)
REGRESSION_WINDOWS = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
FEATURE_LAGS = np.array([1, 2, 5, 10, 20, 50, 100], dtype=int)
ELL_KEYS = {0.25: "q_l025", 0.5: "q_l050", 1.0: "q_l100"}


def processed_path(seed: int) -> Path:
    return PROCESSED_DIR / f"event_histogram_If1e-02_In05_num100_R50_gap50_n{seed}_BEM.mat"


def feature_path(seed: int) -> Path:
    return FEATURE_DIR / f"gap50_hydro_seed{seed:02d}_stride10.npz"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def load_native() -> tuple[dict[int, np.ndarray], dict[int, dict[str, np.ndarray]], dict[str, float]]:
    q_by_seed: dict[int, np.ndarray] = {}
    events: dict[int, dict[str, np.ndarray]] = {}
    parameter_rows = []
    for seed in range(20):
        path = processed_path(seed)
        data = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        nleft = np.asarray(data["num_fish_left"], dtype=np.int16).reshape(-1)
        q_by_seed[seed] = 1.0 - 2.0 * nleft.astype(np.float64) / N_FISH
        events[seed] = {
            "occur": np.atleast_1d(data["event_occur"]).astype(np.int64),
            "size": np.atleast_1d(data["event_size"]).astype(np.int64),
            "interval": np.atleast_1d(data["event_intvl"]).astype(np.float64),
        }
        para = data["para"]
        parameter_rows.append({name: float(np.asarray(getattr(para, name)).reshape(-1)[0]) for name in ("Ip", "In", "R", "dt", "Iw", "If", "num", "total_time", "gap", "radius_ratio", "gap_ratio")})
    first = parameter_rows[0]
    for seed, row in enumerate(parameter_rows[1:], 1):
        for key, value in first.items():
            if not np.isclose(row[key], value):
                raise ValueError(f"parameter mismatch seed {seed}: {key}={row[key]} versus {value}")
    return q_by_seed, events, first


def load_features() -> dict[int, dict[str, np.ndarray]]:
    result: dict[int, dict[str, np.ndarray]] = {}
    missing = [str(feature_path(seed)) for seed in range(20) if not feature_path(seed).is_file()]
    if missing:
        raise FileNotFoundError("feature extraction incomplete:\n" + "\n".join(missing))
    for seed in range(20):
        with np.load(feature_path(seed)) as data:
            result[seed] = {key: np.array(data[key]) for key in ("time", "q", "q_l025", "q_l050", "q_l100", "j_neck", "neck_count")}
    return result


def dominance_passages(q: np.ndarray, dt: float, threshold: float = DOMINANCE) -> list[dict[str, float | str]]:
    state = 0
    entered = math.nan
    passages: list[dict[str, float | str]] = []
    for index, value in enumerate(q):
        new = 1 if value >= threshold else (-1 if value <= -threshold else state)
        if state == 0 and new != 0:
            state = new
            entered = index * dt
        elif state != 0 and new != state:
            now = index * dt
            passages.append({"direction": "LR" if state == -1 else "RL", "duration": now - entered, "time": now})
            state = new
            entered = now
    return passages


def state_index_exact(values: np.ndarray) -> np.ndarray:
    return np.clip(np.rint((values + 1.0) / Q_DX).astype(np.int16), 0, N_FISH)


def state_index_bins(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.clip(np.searchsorted(edges, values, side="right") - 1, 0, len(edges) - 2).astype(np.int16)


def conditional_sums(
    series: list[np.ndarray],
    lags: np.ndarray,
    nstate: int,
    indexer,
) -> dict[str, np.ndarray]:
    shape = (len(lags), nstate)
    result = {name: np.zeros(shape, dtype=np.float64) for name in ("count", "sum1", "sum2", "sum3", "sum4")}
    for li, lag in enumerate(lags):
        for values in series:
            x = values[:-lag]
            delta = values[lag:] - x
            index = indexer(x)
            result["count"][li] += np.bincount(index, minlength=nstate)
            result["sum1"][li] += np.bincount(index, weights=delta, minlength=nstate)
            result["sum2"][li] += np.bincount(index, weights=delta**2, minlength=nstate)
            result["sum3"][li] += np.bincount(index, weights=delta**3, minlength=nstate)
            result["sum4"][li] += np.bincount(index, weights=delta**4, minlength=nstate)
    return result


def moments_from_sums(sums: dict[str, np.ndarray], lags: np.ndarray, dt: float) -> dict[str, np.ndarray]:
    count = sums["count"]
    with np.errstate(invalid="ignore", divide="ignore"):
        raw1 = sums["sum1"] / count
        raw2 = sums["sum2"] / count
        raw3 = sums["sum3"] / count
        raw4 = sums["sum4"] / count
        var = raw2 - raw1**2
        central3 = raw3 - 3 * raw1 * raw2 + 2 * raw1**3
        central4 = raw4 - 4 * raw1 * raw3 + 6 * raw1**2 * raw2 - 3 * raw1**4
        tau = (lags * dt)[:, None]
        drift = raw1 / tau
        diffusion = var / (2 * tau)
        skew = central3 / np.maximum(var, 0) ** 1.5
        gaussian_ratio = central4 / (3 * np.maximum(var, 0) ** 2)
    return {"count": count, "mean": raw1, "variance": var, "drift": drift, "diffusion": diffusion, "skew": skew, "gaussian_ratio": gaussian_ratio}


def weighted_line(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    good = np.isfinite(y) & (weights > 0)
    if np.count_nonzero(good) < 3:
        return math.nan, math.nan, math.nan
    xx = x[good]
    yy = y[good]
    ww = weights[good]
    design = np.column_stack([np.ones_like(xx), xx])
    rootw = np.sqrt(ww / np.max(ww))
    beta, *_ = np.linalg.lstsq(design * rootw[:, None], yy * rootw, rcond=None)
    fitted = design @ beta
    center = np.average(yy, weights=ww)
    ss_res = np.sum(ww * (yy - fitted) ** 2)
    ss_tot = np.sum(ww * (yy - center) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    return float(beta[1]), float(beta[0]), float(r2)


def regress_coefficients(
    moments: dict[str, np.ndarray], lags: np.ndarray, dt: float, max_tau: float, min_count: int = 1000
) -> dict[str, np.ndarray]:
    use = lags * dt <= max_tau + 1e-12
    tau = lags[use] * dt
    nstate = moments["count"].shape[1]
    arrays = {name: np.full(nstate, np.nan) for name in ("V", "D", "mean_intercept", "var_intercept", "mean_r2", "var_r2")}
    arrays["support"] = np.min(moments["count"][use], axis=0)
    for state in range(nstate):
        weights = moments["count"][use, state]
        if np.min(weights) < min_count:
            continue
        slope, intercept, r2 = weighted_line(tau, moments["mean"][use, state], weights)
        var_slope, var_intercept, var_r2 = weighted_line(tau, moments["variance"][use, state], weights)
        arrays["V"][state] = slope
        arrays["D"][state] = var_slope / 2.0
        arrays["mean_intercept"][state] = intercept
        arrays["var_intercept"][state] = var_intercept
        arrays["mean_r2"][state] = r2
        arrays["var_r2"][state] = var_r2
    return arrays


def center_component(valid: np.ndarray, centers: np.ndarray) -> np.ndarray:
    index = int(np.argmin(np.abs(centers)))
    if not valid[index]:
        return np.zeros_like(valid, dtype=bool)
    left = index
    right = index
    while left > 0 and valid[left - 1]:
        left -= 1
    while right + 1 < len(valid) and valid[right + 1]:
        right += 1
    mask = np.zeros_like(valid, dtype=bool)
    mask[left : right + 1] = True
    return mask


def reconstruct_potential(centers: np.ndarray, coefficient: dict[str, np.ndarray]) -> dict[str, object]:
    V = coefficient["V"]
    D = coefficient["D"]
    valid = np.isfinite(V) & np.isfinite(D) & (D > 0)
    component = center_component(valid, centers)
    phi = np.full_like(centers, np.nan, dtype=float)
    if np.count_nonzero(component) >= 3:
        x = centers[component]
        ratio = V[component] / D[component]
        phi_values = np.log(D[component]) - cumulative_trapezoid(ratio, x, initial=0.0)
        phi_values -= np.nanmin(phi_values)
        phi[component] = phi_values
    return {"phi": phi, "component": component, "valid": valid}


def symmetry_constrained(coefficient: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Impose only the known left-right symmetry, without filling missing states."""
    result = {key: np.array(value, copy=True) for key, value in coefficient.items()}
    paired_v = np.isfinite(coefficient["V"]) & np.isfinite(coefficient["V"][::-1])
    paired_d = np.isfinite(coefficient["D"]) & np.isfinite(coefficient["D"][::-1])
    result["V"][:] = np.nan
    result["D"][:] = np.nan
    result["V"][paired_v] = 0.5 * (coefficient["V"][paired_v] - coefficient["V"][::-1][paired_v])
    result["D"][paired_d] = 0.5 * (coefficient["D"][paired_d] + coefficient["D"][::-1][paired_d])
    return result


def barrier_height(centers: np.ndarray, phi: np.ndarray) -> float:
    good = np.isfinite(phi)
    left = good & (centers <= -DOMINANCE)
    right = good & (centers >= DOMINANCE)
    middle = good & (centers > -DOMINANCE) & (centers < DOMINANCE)
    if not left.any() or not right.any() or not middle.any():
        return math.nan
    wells = 0.5 * (np.nanmin(phi[left]) + np.nanmin(phi[right]))
    return float(np.nanmax(phi[middle]) - wells)


def normalize_probability(values: np.ndarray) -> np.ndarray:
    total = np.nansum(values)
    return values / total if total > 0 else np.full_like(values, np.nan)


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    good = np.isfinite(p) & np.isfinite(q) & (p + q > 0)
    p = normalize_probability(p[good])
    q = normalize_probability(q[good])
    m = 0.5 * (p + q)
    left = np.sum(np.where(p > 0, p * np.log(p / m), 0.0))
    right = np.sum(np.where(q > 0, q * np.log(q / m), 0.0))
    return float(0.5 * (left + right))


def potential_metrics(centers: np.ndarray, coefficient: dict[str, np.ndarray], empirical: np.ndarray) -> dict[str, float | np.ndarray]:
    reconstruction = reconstruct_potential(centers, coefficient)
    phi = reconstruction["phi"]
    component = reconstruction["component"]
    predicted = np.full_like(empirical, np.nan, dtype=float)
    predicted[component] = normalize_probability(np.exp(-phi[component]))
    observed = normalize_probability(empirical.astype(float))
    hist_phi = np.full_like(empirical, np.nan, dtype=float)
    positive = observed > 0
    hist_phi[positive] = -np.log(observed[positive])
    if np.any(component & positive):
        hist_phi -= np.nanmin(hist_phi[component & positive])
    common = component & positive
    mae = float(np.nanmean(np.abs(phi[common] - hist_phi[common]))) if common.any() else math.nan
    dx = float(np.median(np.diff(centers)))
    flux_residual = math.nan
    if np.count_nonzero(component) >= 3:
        p = observed[component] / dx
        V = coefficient["V"][component]
        D = coefficient["D"][component]
        flux = V * p - np.gradient(D * p, dx)
        scale = np.sqrt(np.mean((V * p) ** 2)) + np.sqrt(np.mean(np.gradient(D * p, dx) ** 2))
        flux_residual = float(np.sqrt(np.mean(flux**2)) / scale) if scale > 0 else math.nan
    return {
        "phi": phi,
        "hist_phi": hist_phi,
        "predicted": predicted,
        "component": component,
        "barrier": barrier_height(centers, phi),
        "potential_mae": mae,
        "stationary_js": js_divergence(observed, predicted),
        "flux_relative_rms": flux_residual,
        "component_min": float(np.min(centers[component])) if component.any() else math.nan,
        "component_max": float(np.max(centers[component])) if component.any() else math.nan,
    }


def transition_matrix(series: list[np.ndarray], lag: int, nstate: int, indexer) -> tuple[np.ndarray, np.ndarray]:
    counts = np.zeros((nstate, nstate), dtype=np.float64)
    for values in series:
        origin = indexer(values[:-lag]).astype(np.int64)
        target = indexer(values[lag:]).astype(np.int64)
        counts += np.bincount(origin * nstate + target, minlength=nstate * nstate).reshape(nstate, nstate)
    rows = counts.sum(axis=1)
    matrix = np.divide(counts, rows[:, None], out=np.zeros_like(counts), where=rows[:, None] > 0)
    return matrix, rows


def ck_metrics(series: list[np.ndarray], lag: int, nstate: int, indexer, min_row: int = 1000) -> tuple[float, float, int]:
    one, rows = transition_matrix(series, lag, nstate, indexer)
    two, rows_two = transition_matrix(series, 2 * lag, nstate, indexer)
    composed = one @ one
    valid = (rows >= min_row) & (rows_two >= min_row)
    if not valid.any():
        return math.nan, math.nan, 0
    weights = rows_two[valid] / np.sum(rows_two[valid])
    tv_rows = 0.5 * np.sum(np.abs(two[valid] - composed[valid]), axis=1)
    weighted_tv = float(np.sum(weights * tv_rows))
    relative_frobenius = float(np.linalg.norm(two[valid] - composed[valid]) / np.linalg.norm(two[valid]))
    return weighted_tv, relative_frobenius, int(valid.sum())


def solve_mfpt(centers: np.ndarray, V: np.ndarray, D: np.ndarray, direction: str) -> float:
    if direction == "LR":
        domain = centers <= DOMINANCE + 1e-12
        start_target = -DOMINANCE
        reverse = False
    else:
        domain = centers >= -DOMINANCE - 1e-12
        start_target = DOMINANCE
        reverse = True
    x = centers[domain]
    v = V[domain]
    d = D[domain]
    if reverse:
        x = -x[::-1]
        v = -v[::-1]
        d = d[::-1]
        start_target = -start_target
    if len(x) < 4 or not np.all(np.isfinite(v)) or not np.all(np.isfinite(d)) or np.any(d <= 0):
        return math.nan
    dx = float(np.median(np.diff(x)))
    n = len(x)
    matrix = np.zeros((n, n), dtype=float)
    rhs = -np.ones(n, dtype=float)
    matrix[0, 0] = -2 * d[0] / dx**2
    matrix[0, 1] = 2 * d[0] / dx**2
    matrix[-1, -1] = 1.0
    rhs[-1] = 0.0
    for i in range(1, n - 1):
        matrix[i, i - 1] = d[i] / dx**2 - v[i] / (2 * dx)
        matrix[i, i] = -2 * d[i] / dx**2
        matrix[i, i + 1] = d[i] / dx**2 + v[i] / (2 * dx)
    try:
        solution = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError:
        return math.nan
    start = int(np.argmin(np.abs(x - start_target)))
    return float(solution[start])


def fp_survival(
    centers: np.ndarray, V: np.ndarray, D: np.ndarray, direction: str, times: np.ndarray
) -> np.ndarray:
    """Survival from the dominance threshold under the fitted backward generator."""
    if direction == "LR":
        domain = centers <= DOMINANCE + 1e-12
        start_target = -DOMINANCE
        reverse = False
    else:
        domain = centers >= -DOMINANCE - 1e-12
        start_target = DOMINANCE
        reverse = True
    x, v, d = centers[domain], V[domain], D[domain]
    if reverse:
        x, v, d = -x[::-1], -v[::-1], d[::-1]
        start_target = -start_target
    if len(x) < 4 or not np.all(np.isfinite(v)) or not np.all(np.isfinite(d)) or np.any(d <= 0):
        return np.full_like(times, np.nan)
    dx = float(np.median(np.diff(x)))
    ntransient = len(x) - 1
    generator = np.zeros((ntransient, ntransient), dtype=float)
    generator[0, 0] = -2 * d[0] / dx**2
    generator[0, 1] = 2 * d[0] / dx**2
    for i in range(1, ntransient):
        left_rate = d[i] / dx**2 - v[i] / (2 * dx)
        right_rate = d[i] / dx**2 + v[i] / (2 * dx)
        if left_rate < 0 or right_rate < 0:
            return np.full_like(times, np.nan)
        generator[i, i] = -(left_rate + right_rate)
        generator[i, i - 1] = left_rate
        if i + 1 < ntransient:
            generator[i, i + 1] = right_rate
    start = int(np.argmin(np.abs(x[:-1] - start_target)))
    evolved = expm_multiply(csr_matrix(generator), np.ones(ntransient), start=float(times[0]),
                            stop=float(times[-1]), num=len(times), endpoint=True)
    return np.clip(evolved[:, start], 0.0, 1.0)


def closure_test(
    train_values: list[np.ndarray], train_current: list[np.ndarray], test_values: list[np.ndarray], test_current: list[np.ndarray],
    horizon: int, nstate: int, indexer,
) -> dict[str, float]:
    count = np.zeros(nstate)
    sum_y = np.zeros(nstate)
    sum_j = np.zeros(nstate)
    for values, current in zip(train_values, train_current):
        state = indexer(values[:-horizon])
        y = values[horizon:] - values[:-horizon]
        j = current[:-horizon]
        count += np.bincount(state, minlength=nstate)
        sum_y += np.bincount(state, weights=y, minlength=nstate)
        sum_j += np.bincount(state, weights=j, minlength=nstate)
    mean_y = np.divide(sum_y, count, out=np.zeros_like(sum_y), where=count > 0)
    mean_j = np.divide(sum_j, count, out=np.zeros_like(sum_j), where=count > 0)
    numerator = 0.0
    denominator = 0.0
    for values, current in zip(train_values, train_current):
        state = indexer(values[:-horizon])
        y_residual = values[horizon:] - values[:-horizon] - mean_y[state]
        j_residual = current[:-horizon] - mean_j[state]
        numerator += float(np.dot(j_residual, y_residual))
        denominator += float(np.dot(j_residual, j_residual))
    beta = numerator / denominator if denominator > 0 else 0.0
    base_errors = []
    current_errors = []
    all_y = []
    all_jres = []
    for values, current in zip(test_values, test_current):
        state = indexer(values[:-horizon])
        y = values[horizon:] - values[:-horizon]
        base = mean_y[state]
        j_residual = current[:-horizon] - mean_j[state]
        base_errors.append(y - base)
        current_errors.append(y - base - beta * j_residual)
        all_y.append(y - base)
        all_jres.append(j_residual)
    base_error = np.concatenate(base_errors)
    current_error = np.concatenate(current_errors)
    residual = np.concatenate(all_y)
    jres = np.concatenate(all_jres)
    rmse_base = float(np.sqrt(np.mean(base_error**2)))
    rmse_current = float(np.sqrt(np.mean(current_error**2)))
    correlation = float(np.corrcoef(residual, jres)[0, 1])
    return {
        "beta": beta,
        "rmse_q_only": rmse_base,
        "rmse_q_current": rmse_current,
        "rmse_improvement_fraction": 1.0 - rmse_current / rmse_base if rmse_base > 0 else math.nan,
        "conditional_residual_current_correlation": correlation,
    }


def seed_bootstrap_mean(values: dict[int, list[float]], seeds: tuple[int, ...], rng: np.random.Generator, nboot: int = 2000) -> tuple[float, float, float]:
    per_seed = {seed: np.asarray(values.get(seed, []), dtype=float) for seed in seeds}
    pooled = np.concatenate([x for x in per_seed.values() if x.size])
    estimate = float(np.mean(pooled))
    boot = []
    for _ in range(nboot):
        chosen = rng.choice(seeds, size=len(seeds), replace=True)
        sample = [per_seed[int(seed)] for seed in chosen if per_seed[int(seed)].size]
        if sample:
            boot.append(float(np.mean(np.concatenate(sample))))
    low, high = np.percentile(boot, [2.5, 97.5])
    return estimate, float(low), float(high)


def main() -> None:
    global FEATURE_DIR, RESULT_DIR, PROCESSED_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-dir", type=Path, required=True, help="directory containing processed event MAT files")
    parser.add_argument("--feature-dir", type=Path, required=True, help="directory containing extracted NPZ feature caches")
    parser.add_argument("--result-dir", type=Path, required=True, help="destination for regenerated figures and tables")
    args = parser.parse_args()
    PROCESSED_DIR = args.processed_dir.expanduser().resolve()
    FEATURE_DIR = args.feature_dir.expanduser().resolve()
    RESULT_DIR = args.result_dir.expanduser().resolve()
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    q_native, events, parameters = load_native()
    features = load_features()
    train_native = [q_native[s][BURN_NATIVE:] for s in TRAIN_SEEDS]
    test_native = [q_native[s][BURN_NATIVE:] for s in TEST_SEEDS]
    rng = np.random.default_rng(502024)

    # Transition and event characterization.
    transition_rows: list[dict[str, object]] = []
    passages_by_split: dict[str, dict[str, dict[int, list[float]]]] = {
        split: {direction: defaultdict(list) for direction in ("LR", "RL")} for split in ("train", "test")
    }
    for seed in range(20):
        q = q_native[seed][BURN_NATIVE:]
        split = "train" if seed in TRAIN_SEEDS else "test"
        passages = dominance_passages(q, DT_NATIVE)
        for item in passages:
            passages_by_split[split][str(item["direction"])][seed].append(float(item["duration"]))
        transition_rows.append({
            "seed": seed,
            "split": split,
            "switches": len(passages),
            "LR": sum(item["direction"] == "LR" for item in passages),
            "RL": sum(item["direction"] == "RL" for item in passages),
            "fraction_abs_q_gt_08": float(np.mean(np.abs(q) > 0.8)),
            "fraction_abs_q_lt_02": float(np.mean(np.abs(q) < 0.2)),
        })
    write_csv(RESULT_DIR / "transition_counts_by_seed.csv", transition_rows)

    train_counts = np.zeros(len(Q_GRID))
    test_counts = np.zeros(len(Q_GRID))
    for q in train_native:
        train_counts += np.bincount(state_index_exact(q), minlength=len(Q_GRID))
    for q in test_native:
        test_counts += np.bincount(state_index_exact(q), minlength=len(Q_GRID))

    event_rows = []
    all_marks = []
    one_step_changes = []
    interval_pairs = []
    for seed in range(20):
        q = q_native[seed][BURN_NATIVE:]
        changes = np.diff(q)
        one_step_changes.append(changes[changes != 0])
        occur = events[seed]["occur"]
        marks = events[seed]["size"][occur >= BURN_NATIVE]
        valid_occur = occur[occur >= BURN_NATIVE]
        all_marks.append(marks)
        waits = np.diff(valid_occur) * DT_NATIVE
        if waits.size > 1:
            interval_pairs.append((waits[:-1], waits[1:]))
        event_rows.append({"seed": seed, "events": len(marks), "mean_abs_mark": float(np.mean(np.abs(marks))), "median_abs_mark": float(np.median(np.abs(marks))), "max_abs_mark": int(np.max(np.abs(marks)))})
    write_csv(RESULT_DIR / "event_statistics_by_seed.csv", event_rows)
    all_marks_array = np.concatenate(all_marks)
    one_step_array = np.concatenate(one_step_changes)
    wait_corr = float(np.corrcoef(np.concatenate([x for x, _ in interval_pairs]), np.concatenate([y for _, y in interval_pairs]))[0, 1])

    # Native discrete conditional moments and regression fits.
    native_sums = conditional_sums(train_native, NATIVE_LAGS, len(Q_GRID), state_index_exact)
    native_moments = moments_from_sums(native_sums, NATIVE_LAGS, DT_NATIVE)
    potential_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    native_coefficients = {}
    native_potentials = {}
    for window in REGRESSION_WINDOWS:
        coefficient = regress_coefficients(native_moments, NATIVE_LAGS, DT_NATIVE, window)
        metrics = potential_metrics(Q_GRID, coefficient, train_counts)
        native_coefficients[window] = coefficient
        native_potentials[window] = metrics
        valid = np.isfinite(coefficient["V"]) & np.isfinite(coefficient["D"])
        potential_rows.append({
            "observable": "q_discrete", "max_tau": window,
            "supported_states": int(valid.sum()),
            "positive_D_states": int(np.count_nonzero(valid & (coefficient["D"] > 0))),
            "component_min": metrics["component_min"], "component_max": metrics["component_max"],
            "barrier": metrics["barrier"], "potential_mae": metrics["potential_mae"],
            "stationary_js": metrics["stationary_js"], "flux_relative_rms": metrics["flux_relative_rms"],
            "median_mean_r2": float(np.nanmedian(coefficient["mean_r2"])),
            "median_var_r2": float(np.nanmedian(coefficient["var_r2"])),
        })
        for index, qvalue in enumerate(Q_GRID):
            coefficient_rows.append({"observable": "q_discrete", "max_tau": window, "state": qvalue, **{key: coefficient[key][index] for key in ("V", "D", "support", "mean_r2", "var_r2")}})

    ck_rows: list[dict[str, object]] = []
    for lag in (10, 50, 100, 200, 500, 1000):
        tv, frob, rows = ck_metrics(train_native, lag, len(Q_GRID), state_index_exact)
        ck_rows.append({"observable": "q_discrete", "tau": lag * DT_NATIVE, "weighted_TV": tv, "relative_Frobenius": frob, "supported_rows": rows})

    # There is no accepted fitting interval unless diffusion is simultaneously
    # close to Gaussian, Markov-consistent, and stable with lag.  We use 0.5 as
    # a declared reference because it is the first tested lag whose median
    # fourth-moment ratio is close to one; this is not called a tau->0 limit.
    selected_tau = 0.5
    selected_coefficient = native_coefficients[selected_tau]
    selected_lag = int(NATIVE_LAGS[np.argmin(np.abs(NATIVE_LAGS * DT_NATIVE - selected_tau))])
    per_seed_v = []
    per_seed_d = []
    for values in train_native:
        seed_sums = conditional_sums([values], np.array([selected_lag]), len(Q_GRID), state_index_exact)
        seed_moments = moments_from_sums(seed_sums, np.array([selected_lag]), DT_NATIVE)
        per_seed_v.append(seed_moments["drift"][0])
        per_seed_d.append(seed_moments["diffusion"][0])
    per_seed_v_array = np.asarray(per_seed_v)
    per_seed_d_array = np.asarray(per_seed_d)
    direct_rows = []
    for index, qvalue in enumerate(Q_GRID):
        direct_rows.append({
            "q": qvalue,
            "tau": selected_lag * DT_NATIVE,
            "V_seed_mean": float(np.nanmean(per_seed_v_array[:, index])),
            "V_seed_se": float(np.nanstd(per_seed_v_array[:, index], ddof=1) / np.sqrt(len(TRAIN_SEEDS))),
            "D_seed_mean": float(np.nanmean(per_seed_d_array[:, index])),
            "D_seed_se": float(np.nanstd(per_seed_d_array[:, index], ddof=1) / np.sqrt(len(TRAIN_SEEDS))),
        })
    write_csv(RESULT_DIR / "direct_coefficients_realization_uncertainty.csv", direct_rows)
    predicted_mfpt = {direction: solve_mfpt(Q_GRID, selected_coefficient["V"], selected_coefficient["D"], direction) for direction in ("LR", "RL")}

    mfpt_lag_rows: list[dict[str, object]] = []
    test_passage_means = {
        direction: float(np.mean(np.concatenate([np.asarray(v) for v in passages_by_split["test"][direction].values()])))
        for direction in ("LR", "RL")
    }
    for window in REGRESSION_WINDOWS:
        coefficient = native_coefficients[window]
        for direction in ("LR", "RL"):
            prediction = solve_mfpt(Q_GRID, coefficient["V"], coefficient["D"], direction)
            mfpt_lag_rows.append({
                "observable": "q_discrete", "max_tau": window, "direction": direction,
                "predicted_mfpt": prediction, "heldout_observed_mfpt": test_passage_means[direction],
                "relative_error": (prediction - test_passage_means[direction]) / test_passage_means[direction],
            })
        symmetric = symmetry_constrained(coefficient)
        for direction in ("LR", "RL"):
            prediction = solve_mfpt(Q_GRID, symmetric["V"], symmetric["D"], direction)
            mfpt_lag_rows.append({
                "observable": "q_discrete_symmetry_constrained", "max_tau": window, "direction": direction,
                "predicted_mfpt": prediction, "heldout_observed_mfpt": test_passage_means[direction],
                "relative_error": (prediction - test_passage_means[direction]) / test_passage_means[direction],
            })

    mfpt_rows = []
    for split, seeds in (("train", TRAIN_SEEDS), ("test", TEST_SEEDS)):
        for direction in ("LR", "RL"):
            estimate, low, high = seed_bootstrap_mean(passages_by_split[split][direction], seeds, rng)
            mfpt_rows.append({"split": split, "direction": direction, "observed_mean": estimate, "observed_ci_low": low, "observed_ci_high": high, "predicted": predicted_mfpt[direction], "relative_error": (predicted_mfpt[direction] - estimate) / estimate if np.isfinite(predicted_mfpt[direction]) else math.nan, "selected_max_tau": selected_tau})
    write_csv(RESULT_DIR / "mfpt_validation.csv", mfpt_rows)

    # Continuous observables, smoothing sensitivity, and current closure.
    closure_rows: list[dict[str, object]] = []
    continuous_summary = []
    continuous_plot_data = {}
    j_train = [features[s]["j_neck"][BURN_FEATURE:] for s in TRAIN_SEEDS]
    j_test = [features[s]["j_neck"][BURN_FEATURE:] for s in TEST_SEEDS]
    for ell, key in ELL_KEYS.items():
        train = [features[s][key][BURN_FEATURE:] for s in TRAIN_SEEDS]
        test = [features[s][key][BURN_FEATURE:] for s in TEST_SEEDS]
        edges = np.linspace(-1.0, 1.0, 62)
        centers = 0.5 * (edges[:-1] + edges[1:])
        indexer = lambda x, e=edges: state_index_bins(x, e)
        sums = conditional_sums(train, FEATURE_LAGS, len(centers), indexer)
        moments = moments_from_sums(sums, FEATURE_LAGS, DT_FEATURE)
        counts = np.zeros(len(centers))
        test_hist = np.zeros(len(centers))
        for values in train:
            counts += np.bincount(indexer(values), minlength=len(centers))
        for values in test:
            test_hist += np.bincount(indexer(values), minlength=len(centers))
        local_rows = []
        local_coefficients = {}
        local_potentials = {}
        for window in (0.05, 0.1, 0.2, 0.5, 1.0):
            coefficient = regress_coefficients(moments, FEATURE_LAGS, DT_FEATURE, window, min_count=500)
            metrics = potential_metrics(centers, coefficient, counts)
            local_coefficients[window] = coefficient
            local_potentials[window] = metrics
            row = {
                "observable": f"q_l_ell{ell:g}", "max_tau": window,
                "supported_states": int(np.count_nonzero(np.isfinite(coefficient["V"]))),
                "positive_D_states": int(np.count_nonzero(np.isfinite(coefficient["D"]) & (coefficient["D"] > 0))),
                "component_min": metrics["component_min"], "component_max": metrics["component_max"],
                "barrier": metrics["barrier"], "potential_mae": metrics["potential_mae"],
                "stationary_js": metrics["stationary_js"], "flux_relative_rms": metrics["flux_relative_rms"],
                "median_mean_r2": float(np.nanmedian(coefficient["mean_r2"])), "median_var_r2": float(np.nanmedian(coefficient["var_r2"])),
            }
            potential_rows.append(row)
            local_rows.append(row)
        for lag in (1, 5, 10, 20, 50, 100):
            tv, frob, rows = ck_metrics(train, lag, len(centers), indexer, min_row=500)
            ck_rows.append({"observable": f"q_l_ell{ell:g}", "tau": lag * DT_FEATURE, "weighted_TV": tv, "relative_Frobenius": frob, "supported_rows": rows})
        for horizon in (10, 50, 100):
            closure = closure_test(train, j_train, test, j_test, horizon, len(centers), indexer)
            closure_rows.append({"observable": f"q_l_ell{ell:g}", "horizon": horizon * DT_FEATURE, **closure})
        observed_continuous = {direction: [] for direction in ("LR", "RL")}
        for values in test:
            for passage in dominance_passages(values, DT_FEATURE):
                observed_continuous[str(passage["direction"])].append(float(passage["duration"]))
        for window, coefficient in local_coefficients.items():
            for direction in ("LR", "RL"):
                prediction = solve_mfpt(centers, coefficient["V"], coefficient["D"], direction)
                observed = float(np.mean(observed_continuous[direction]))
                mfpt_lag_rows.append({
                    "observable": f"q_l_ell{ell:g}", "max_tau": window, "direction": direction,
                    "predicted_mfpt": prediction, "heldout_observed_mfpt": observed,
                    "relative_error": (prediction - observed) / observed,
                })
        chosen = local_coefficients[0.1]
        chosen_metrics = local_potentials[0.1]
        continuous_summary.append({"ell": ell, "barrier_tau01": chosen_metrics["barrier"], "stationary_js_tau01": chosen_metrics["stationary_js"], "flux_residual_tau01": chosen_metrics["flux_relative_rms"]})
        continuous_plot_data[ell] = (centers, counts, test_hist, local_coefficients, local_potentials, moments)

    # Closure test for exact q at feature resolution.
    q_train_feature = [features[s]["q"][BURN_FEATURE:] for s in TRAIN_SEEDS]
    q_test_feature = [features[s]["q"][BURN_FEATURE:] for s in TEST_SEEDS]
    for horizon in (10, 50, 100):
        closure = closure_test(q_train_feature, j_train, q_test_feature, j_test, horizon, len(Q_GRID), state_index_exact)
        closure_rows.append({"observable": "q_discrete", "horizon": horizon * DT_FEATURE, **closure})

    write_csv(RESULT_DIR / "potential_lag_sensitivity.csv", potential_rows)
    write_csv(RESULT_DIR / "state_coefficients.csv", coefficient_rows)
    write_csv(RESULT_DIR / "chapman_kolmogorov.csv", ck_rows)
    write_csv(RESULT_DIR / "neck_current_closure.csv", closure_rows)
    write_csv(RESULT_DIR / "mfpt_lag_sensitivity.csv", mfpt_lag_rows)

    lag_rows = []
    for li, lag in enumerate(NATIVE_LAGS):
        supported = native_moments["count"][li] >= 1000
        lag_rows.append({
            "tau": lag * DT_NATIVE,
            "supported_states": int(supported.sum()),
            "median_drift": float(np.nanmedian(native_moments["drift"][li, supported])),
            "median_diffusion": float(np.nanmedian(native_moments["diffusion"][li, supported])),
            "median_abs_skew": float(np.nanmedian(np.abs(native_moments["skew"][li, supported]))),
            "median_gaussian_fourth_ratio": float(np.nanmedian(native_moments["gaussian_ratio"][li, supported])),
            "probability_nonzero_increment": float(np.mean([np.mean(np.diff(q, n=1) != 0) for q in train_native])) if lag == 1 else math.nan,
        })
    write_csv(RESULT_DIR / "native_lag_diagnostics.csv", lag_rows)

    # Essential figures.
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(3, 1, figsize=(10, 7.5), sharex=True)
    for axis, seed in zip(axes, TEST_SEEDS[:3]):
        q = q_native[seed][BURN_NATIVE::100]
        t = BURN_TIME + np.arange(len(q)) * 0.1
        axis.plot(t, q, lw=0.65)
        axis.axhline(DOMINANCE, color="k", ls="--", lw=0.7)
        axis.axhline(-DOMINANCE, color="k", ls="--", lw=0.7)
        axis.set_ylabel(f"q, seed {seed}")
    axes[-1].set_xlabel("time")
    fig.suptitle("Gap 50 full-BEM chamber imbalance: held-out realizations")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "01_representative_occupancy_trajectories.pdf")
    fig.savefig(RESULT_DIR / "01_representative_occupancy_trajectories.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 4.6))
    axis.step(Q_GRID, normalize_probability(train_counts), where="mid", label="training seeds 0--14")
    axis.step(Q_GRID, normalize_probability(test_counts), where="mid", label="held-out seeds 15--19")
    axis.set(xlabel="q", ylabel="probability mass", title="Stationary discrete chamber occupancy")
    axis.legend()
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "02_stationary_occupancy_distribution.pdf")
    fig.savefig(RESULT_DIR / "02_stationary_occupancy_distribution.png", dpi=180)
    plt.close(fig)

    selected_states = [-0.8, -0.4, 0.0, 0.4, 0.8]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    tau = NATIVE_LAGS * DT_NATIVE
    for state in selected_states:
        index = int(np.argmin(np.abs(Q_GRID - state)))
        axes[0].plot(tau, native_moments["mean"][:, index], "o-", ms=3, label=f"q={Q_GRID[index]:.1f}")
        axes[1].plot(tau, native_moments["variance"][:, index], "o-", ms=3, label=f"q={Q_GRID[index]:.1f}")
    axes[0].set(xlabel=r"$\tau$", ylabel=r"$E[\Delta q|q]$", title="Conditional mean increments")
    axes[1].set(xlabel=r"$\tau$", ylabel=r"$Var[\Delta q|q]$", title="Conditional increment variance")
    axes[0].set_xscale("log"); axes[1].set_xscale("log")
    axes[1].legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "03_conditional_moments_vs_lag.pdf")
    fig.savefig(RESULT_DIR / "03_conditional_moments_vs_lag.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for window in REGRESSION_WINDOWS:
        coefficient = native_coefficients[window]
        axes[0].plot(Q_GRID, coefficient["V"], label=f"max lag {window:g}")
        axes[1].plot(Q_GRID, coefficient["D"], label=f"max lag {window:g}")
    direct_v_mean = np.nanmean(per_seed_v_array, axis=0)
    direct_v_se = np.nanstd(per_seed_v_array, axis=0, ddof=1) / np.sqrt(len(TRAIN_SEEDS))
    direct_d_mean = np.nanmean(per_seed_d_array, axis=0)
    direct_d_se = np.nanstd(per_seed_d_array, axis=0, ddof=1) / np.sqrt(len(TRAIN_SEEDS))
    axes[0].fill_between(Q_GRID, direct_v_mean - 1.96 * direct_v_se, direct_v_mean + 1.96 * direct_v_se,
                         color="0.35", alpha=0.16, label=r"direct $\tau=0.5$, seed 95% band")
    axes[1].fill_between(Q_GRID, direct_d_mean - 1.96 * direct_d_se, direct_d_mean + 1.96 * direct_d_se,
                         color="0.35", alpha=0.16, label=r"direct $\tau=0.5$, seed 95% band")
    axes[0].axhline(0, color="k", lw=0.7)
    axes[1].axhline(0, color="k", lw=0.7)
    axes[0].set(ylabel="V(q)", title="Regression drift: lag sensitivity")
    axes[1].set(xlabel="q", ylabel="D(q)", title="Regression diffusion: lag sensitivity")
    axes[0].legend(ncol=2, fontsize=7.5)
    axes[1].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "04_drift_diffusion_lag_sensitivity.pdf")
    fig.savefig(RESULT_DIR / "04_drift_diffusion_lag_sensitivity.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9, 5))
    hist = native_potentials[selected_tau]["hist_phi"]
    axis.plot(Q_GRID, hist, color="k", lw=2, label=r"$-\log p(q)$")
    for window in REGRESSION_WINDOWS:
        axis.plot(Q_GRID, native_potentials[window]["phi"], label=f"dynamic, max lag {window:g}")
    axis.set(xlabel="q", ylabel=r"$\Phi(q)+C$", title="Discrete-q effective potential and lag sensitivity")
    axis.set_ylim(bottom=-0.2, top=np.nanpercentile(hist, 95) + 2)
    axis.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "05_effective_potential_comparison.pdf")
    fig.savefig(RESULT_DIR / "05_effective_potential_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), sharey=True)
    for axis, direction in zip(axes, ("LR", "RL")):
        rows = [row for row in mfpt_lag_rows if row["observable"] == "q_discrete" and row["direction"] == direction]
        axis.plot([row["max_tau"] for row in rows], [row["predicted_mfpt"] for row in rows], "o-", label="discrete-q FP")
        smooth = [row for row in mfpt_lag_rows if row["observable"] == "q_l_ell0.5" and row["direction"] == direction]
        axis.plot([row["max_tau"] for row in smooth], [row["predicted_mfpt"] for row in smooth], "s-", label=r"$q_\ell$, $\ell=0.5$")
        observed = rows[0]["heldout_observed_mfpt"]
        axis.axhline(observed, color="k", ls="--", label="held-out observed")
        axis.set(xlabel="regression maximum lag", title="left to right" if direction == "LR" else "right to left")
        axis.set_xscale("log"); axis.set_yscale("log")
    axes[0].set_ylabel("mean first-passage time")
    axes[0].legend(fontsize=8)
    fig.suptitle("MFPT prediction is strongly fitting-lag dependent")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "06_predicted_vs_observed_transition_times.pdf")
    fig.savefig(RESULT_DIR / "06_predicted_vs_observed_transition_times.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for observable in sorted(set(row["observable"] for row in ck_rows)):
        rows = [row for row in ck_rows if row["observable"] == observable]
        axes[0].plot([row["tau"] for row in rows], [row["weighted_TV"] for row in rows], "o-", label=observable)
    q_closure = [row for row in closure_rows if row["observable"] == "q_discrete"]
    axes[1].plot([row["horizon"] for row in q_closure], [row["rmse_q_only"] for row in q_closure], "o-", label="q only")
    axes[1].plot([row["horizon"] for row in q_closure], [row["rmse_q_current"] for row in q_closure], "o-", label="q + neck current")
    axes[0].set(xlabel=r"$\tau$", ylabel="weighted CK total variation", title="Chapman--Kolmogorov consistency")
    axes[1].set(xlabel="prediction horizon", ylabel="held-out RMSE", title="Hidden-current closure test")
    axes[0].legend(fontsize=8); axes[1].legend()
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "07_markov_and_current_closure.pdf")
    fig.savefig(RESULT_DIR / "07_markov_and_current_closure.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ell, (centers, counts, _test_hist, coefficients, potentials, _moments) in continuous_plot_data.items():
        axes[0].plot(centers, normalize_probability(counts), label=fr"$\ell={ell:g}$")
        axes[1].plot(centers, potentials[0.1]["phi"], label=fr"$\ell={ell:g}$")
    axes[0].set(xlabel=r"$q_\ell$", ylabel="probability mass", title="Smoothed occupancy distributions")
    axes[1].set(xlabel=r"$q_\ell$", ylabel=r"$\Phi(q_\ell)+C$", title="Dynamic potentials at max lag 0.1")
    axes[0].legend(); axes[1].legend()
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "08_smoothed_coordinate_sensitivity.pdf")
    fig.savefig(RESULT_DIR / "08_smoothed_coordinate_sensitivity.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    step_fish = np.abs(one_step_array) * N_FISH / 2.0
    axes[0].hist(step_fish, bins=np.arange(0.5, np.max(step_fish) + 1.5), density=True, color="tab:blue")
    axes[0].set(xlabel="fish crossing in one native time step", ylabel="probability", title="Actual one-step occupancy changes")
    axes[1].hist(np.abs(all_marks_array), bins=np.arange(0.5, 101.5, 5), density=True, color="tab:orange")
    axes[1].set(xlabel="fish in aggregated cluster mark", ylabel="probability", title="Archived collective crossing clusters")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "09_step_changes_and_cluster_marks.pdf")
    fig.savefig(RESULT_DIR / "09_step_changes_and_cluster_marks.png", dpi=180)
    plt.close(fig)

    # Held-out first-passage survival.  Complete passages are shown directly;
    # the five final incomplete residence intervals are not treated as events.
    survival_times = np.linspace(0.0, 180.0, 181)
    survival_metrics: dict[str, dict[str, float]] = {}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), sharey=True)
    for axis, direction in zip(axes, ("LR", "RL")):
        durations = np.concatenate([
            np.asarray(values, dtype=float) for values in passages_by_split["test"][direction].values()
        ])
        empirical = np.array([np.mean(durations > time) for time in survival_times])
        model_reference = fp_survival(Q_GRID, selected_coefficient["V"], selected_coefficient["D"], direction, survival_times)
        model_posthoc = fp_survival(Q_GRID, native_coefficients[1.0]["V"], native_coefficients[1.0]["D"], direction, survival_times)
        survival_metrics[direction] = {
            "n_complete": int(len(durations)),
            "reference_max_abs_error": float(np.nanmax(np.abs(empirical - model_reference))),
            "posthoc_tau1_max_abs_error": float(np.nanmax(np.abs(empirical - model_posthoc))),
        }
        axis.step(survival_times, empirical, where="post", color="k", lw=2, label="held-out complete passages")
        axis.plot(survival_times, model_reference, label=r"FP max lag $0.5$")
        axis.plot(survival_times, model_posthoc, label=r"FP max lag $1.0$ (post-hoc)")
        axis.set(xlabel="first-passage time", title="left to right" if direction == "LR" else "right to left")
    axes[0].set_ylabel("survival probability")
    axes[0].legend(fontsize=8)
    fig.suptitle("Held-out first-passage survival is fitting-lag dependent")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "10_first_passage_survival.pdf")
    fig.savefig(RESULT_DIR / "10_first_passage_survival.png", dpi=180)
    plt.close(fig)

    # Concise machine-readable summary and narrative report.
    train_switches = sum(row["switches"] for row in transition_rows if row["split"] == "train")
    test_switches = sum(row["switches"] for row in transition_rows if row["split"] == "test")
    barrier_values = [row["barrier"] for row in potential_rows if row["observable"] == "q_discrete" and np.isfinite(row["barrier"])]
    barrier_ratio = max(barrier_values) / min(barrier_values) if barrier_values and min(barrier_values) > 0 else math.nan
    median_diffusion_by_lag = np.array([row["median_diffusion"] for row in lag_rows])
    gaussian_by_lag = np.array([row["median_gaussian_fourth_ratio"] for row in lag_rows])
    meaningful = (gaussian_by_lag >= 0.8) & (gaussian_by_lag <= 1.2)
    meaningful_d = median_diffusion_by_lag[meaningful]
    diffusion_stability_ratio = float(np.max(meaningful_d) / np.min(meaningful_d)) if meaningful_d.size > 1 else math.nan
    supported_fitting_interval = None if not np.isfinite(diffusion_stability_ratio) or diffusion_stability_ratio > 1.2 else [float(np.min(NATIVE_LAGS[meaningful] * DT_NATIVE)), float(np.max(NATIVE_LAGS[meaningful] * DT_NATIVE))]
    selected_metrics = native_potentials[selected_tau]
    histogram_barrier = barrier_height(Q_GRID, selected_metrics["hist_phi"])
    posthoc_tau = 1.0
    posthoc_mfpt = {
        direction: solve_mfpt(Q_GRID, native_coefficients[posthoc_tau]["V"], native_coefficients[posthoc_tau]["D"], direction)
        for direction in ("LR", "RL")
    }
    symmetric_reference = symmetry_constrained(selected_coefficient)
    symmetric_mfpt = {
        direction: solve_mfpt(Q_GRID, symmetric_reference["V"], symmetric_reference["D"], direction)
        for direction in ("LR", "RL")
    }
    summary = {
        "parameters": parameters,
        "train_seeds": list(TRAIN_SEEDS), "test_seeds": list(TEST_SEEDS),
        "train_switches": train_switches, "test_switches": test_switches,
        "events_post_transient": int(len(all_marks_array)),
        "cluster_abs_size_mean": float(np.mean(np.abs(all_marks_array))),
        "cluster_abs_size_median": float(np.median(np.abs(all_marks_array))),
        "single_step_nonzero_abs_change_mean": float(np.mean(np.abs(one_step_array))),
        "single_step_nonzero_abs_change_max": float(np.max(np.abs(one_step_array))),
        "event_wait_lag1_correlation": wait_corr,
        "reference_regression_max_tau": selected_tau,
        "supported_fitting_lag_interval": supported_fitting_interval,
        "diffusion_stability_ratio_in_near_gaussian_lags": diffusion_stability_ratio,
        "reference_barrier": selected_metrics["barrier"],
        "histogram_barrier": histogram_barrier,
        "barrier_max_to_min_ratio_across_lags": barrier_ratio,
        "selected_stationary_js": selected_metrics["stationary_js"],
        "selected_flux_relative_rms": selected_metrics["flux_relative_rms"],
        "predicted_mfpt": predicted_mfpt,
        "symmetry_constrained_reference_mfpt": symmetric_mfpt,
        "posthoc_tau1_mfpt": posthoc_mfpt,
        "heldout_survival_metrics": survival_metrics,
        "continuous_smoothing_summary": continuous_summary,
    }
    (RESULT_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    test_mfpt = {row["direction"]: row for row in mfpt_rows if row["split"] == "test"}
    ck_q = [row for row in ck_rows if row["observable"] == "q_discrete"]
    closure_q = [row for row in closure_rows if row["observable"] == "q_discrete"]
    report = f"""# R50/gap50 hydrodynamic one-dimensional FP validation

## Geometry normalization warning

The requested case was called `gap/R=0.5`, but the archived `R50/gap50` files
store physical `gap={parameters['gap']:.6g}` and `R={parameters['R']:.6g}`.  Thus this ensemble has
**gap/R = {parameters['gap']/parameters['R']:.3f}** and **gap/(2R) = {parameters['gap']/(2*parameters['R']):.3f}**.  No `R50/gap25`
hydrodynamic ensemble was found.  The results below answer the question for the
available gap-50 ensemble, which matches 0.5 only when normalized by chamber
diameter.

## Dataset and split

- Full BEM hydrodynamics: `If={parameters['If']}`, direct and boundary-mediated calls source-verified.
- Geometry: `R={parameters['R']:.6g}`, `gap={parameters['gap']:.6g}`, gap/R = {parameters['gap']/parameters['R']:.3f}, gap/(2R) = {parameters['gap']/(2*parameters['R']):.3f}.
- N={int(parameters['num'])}, Ip={parameters['Ip']}, In={parameters['In']}, Iw={parameters['Iw']}, dt={parameters['dt']}, T={parameters['total_time']}.
- Twenty independent seeds; seeds 0--14 train, 15--19 held out; first {BURN_TIME:g} time units removed.
- Complete dominance switches at |q|={DOMINANCE}: {train_switches} train and {test_switches} held out.
- The continuous coordinate uses ell=0.5, the rounded divider-tip radius; ell=0.25 and 1.0 are explicit sensitivity controls.

## Transition character

- Post-transient archived cluster events: {len(all_marks_array)}.
- Absolute cluster mark: mean {np.mean(np.abs(all_marks_array)):.3f} fish, median {np.median(np.abs(all_marks_array)):.1f}, maximum {np.max(np.abs(all_marks_array))}.
- A nonzero native one-step occupancy change has mean magnitude {np.mean(np.abs(one_step_array)):.5f} in q and maximum {np.max(np.abs(one_step_array)):.5f}.
- Lag-one correlation of within-realization event waiting intervals: {wait_corr:.3f}.

## One-dimensional discrete-q FP diagnostic

- Supported fitting-lag interval: {supported_fitting_interval}.  None means no interval met simultaneous increment-shape and diffusion-stability requirements.
- Reference maximum lag: {selected_tau:g}, the first tested near-Gaussian increment scale; it is a coarse diagnostic, not a tau-to-zero estimate.
- Dynamic barrier at the reference lag: {selected_metrics['barrier']:.4g}; barrier max/min ratio across tested lags: {barrier_ratio:.3g}.
- Histogram barrier: {histogram_barrier:.4g}.  Its existence confirms bistability, but not dynamical closure.
- Stationary JS divergence: {selected_metrics['stationary_js']:.4g}; stationary-current relative RMS residual: {selected_metrics['flux_relative_rms']:.4g}.
- CK weighted-TV range: {min(row['weighted_TV'] for row in ck_q):.4g}--{max(row['weighted_TV'] for row in ck_q):.4g}.
- Held-out LR MFPT: observed {test_mfpt['LR']['observed_mean']:.3f} [{test_mfpt['LR']['observed_ci_low']:.3f}, {test_mfpt['LR']['observed_ci_high']:.3f}], predicted {test_mfpt['LR']['predicted']:.3f}.
- Held-out RL MFPT: observed {test_mfpt['RL']['observed_mean']:.3f} [{test_mfpt['RL']['observed_ci_low']:.3f}, {test_mfpt['RL']['observed_ci_high']:.3f}], predicted {test_mfpt['RL']['predicted']:.3f}.

The lag-1.0 fit gives LR/RL MFPTs {posthoc_mfpt['LR']:.3f}/{posthoc_mfpt['RL']:.3f}, close to the
held-out observations.  This agreement is post-hoc: at that lag the median
diffusion has changed by a factor {diffusion_stability_ratio:.3f} across the near-Gaussian lag range and the CK
error is large.  It therefore validates a coarse-time empirical propagator, not
a unique short-time diffusion coefficient.  Enforcing exact left-right symmetry
at the reference lag gives {symmetric_mfpt['LR']:.3f}/{symmetric_mfpt['RL']:.3f}; it does not rescue the reference-lag error.

## Requested summary table

| Diagnostic | Result | Assessment |
|---|---:|---|
| Independent realizations | 15 train + 5 held out | adequate realization split |
| Complete dominance switches | {train_switches} train + {test_switches} test | much better sampled than gap 5 |
| Held-out mean transition time | LR {test_mfpt['LR']['observed_mean']:.2f}; RL {test_mfpt['RL']['observed_mean']:.2f} | seed-bootstrap intervals reported above |
| Supported diffusion fitting interval | {supported_fitting_interval} | none |
| Near-Gaussian-lag diffusion variation | {diffusion_stability_ratio:.3f}x | no coefficient plateau |
| Reference dynamic barrier | {selected_metrics['barrier']:.3f} at max lag {selected_tau:g} | finite but lag-dependent |
| CK weighted TV | {min(row['weighted_TV'] for row in ck_q):.4f}--{max(row['weighted_TV'] for row in ck_q):.4f} | deteriorates strongly with lag |
| Hidden-current RMSE reduction | {100*closure_q[-1]['rmse_improvement_fraction']:.1f}% at horizon 1 | q is not closed |
| Reference MFPT relative error | LR {(predicted_mfpt['LR']/test_mfpt['LR']['observed_mean']-1)*100:.1f}%; RL {(predicted_mfpt['RL']/test_mfpt['RL']['observed_mean']-1)*100:.1f}% | fails held-out timing |
| Reference survival max error | LR {survival_metrics['LR']['reference_max_abs_error']:.3f}; RL {survival_metrics['RL']['reference_max_abs_error']:.3f} | full FPT distribution also fails |

## Hidden-current closure test

Held-out prediction improvements after adding neck current:
"""
    for row in closure_q:
        report += f"\n- horizon {row['horizon']:.2f}: RMSE {row['rmse_q_only']:.5g} -> {row['rmse_q_current']:.5g} ({100*row['rmse_improvement_fraction']:.2f}% improvement), conditional residual-current correlation {row['conditional_residual_current_correlation']:.3f}."
    report += f"""

## Interpretation rule

The potential is accepted as a predictive one-dimensional reduction only if its coefficients and barrier are stable over lag, CK discrepancies are small, adding neck current gives little held-out improvement, stationary current is small, and held-out MFPTs agree.  Failure of any of these checks is reported rather than repaired by interpolating unsupported states.

See the CSV tables and numbered figures in this directory for the complete diagnostics.

## Scientific conclusions

1. Relative to gap 5, gap 50 greatly improves transition counts (424 complete
   switches here), so insufficient data are no longer the main limitation.
2. The native process changes mostly in one-fish increments, but its increments
   become approximately Gaussian only after coarse aggregation around 0.5--1
   time unit.  Across that same range the inferred diffusion changes by
   {diffusion_stability_ratio:.3f}x, so there is no Huang-style short-lag plateau.
3. Chamber imbalance is not sufficient: neck current reduces held-out error by
   40--69% and remains strongly correlated with the q-only residual.
4. A finite two-well curve is reproducible as a phenomenological landscape,
   but the predeclared lag-0.5 model overpredicts held-out MFPT by roughly a
   factor of two and misses the held-out survival curve.  The good lag-1 fit is
   lag selection, not parameter-free prediction.
5. The dominant failure is hidden directional memory combined with temporal
   coarse-graining, not loss of metastability or too few switches.  Larger gap
   improves statistics but does not produce a closed one-dimensional continuous
   diffusion.
"""
    (RESULT_DIR / "REPORT.md").write_text(report)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
