#!/usr/bin/env python3
"""Macroscopic localization and transport across social-alignment strength."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

from analyze_approach_current_balance import RawReader


HERE = Path(__file__).resolve().parent
RAW_DEFAULT = Path("/Volumes/storage2/BEM_bichamber/archive_one_R_60/BEM_data")
CACHE_DEFAULT = HERE / "outputs/initiation_propagation/cache"
OUT_DEFAULT = HERE / "outputs/alignment_localization_regime"
N = 100
DT = 0.1
LATE_START = 200.0
STRONG = 0.8
BASIN = 0.6
Q_VALUES = np.arange(-N, N + 1, 2) / N
FILE_RE = re.compile(
    r"If(?P<If>[^_]+)_Ip(?P<Ip>\d+)_In(?P<In>\d+)_R(?P<R>\d+)_num(?P<N>\d+)_wid(?P<gap>\d+)_n(?P<seed>\d+)_BEM(?P<nohydro>_no_hydro)?\.mat$"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def independent_probability() -> np.ndarray:
    nright = np.arange(N + 1)
    probability = binom.pmf(nright, N, 0.5)
    return probability


def inventory(raw_dir: Path) -> list[dict]:
    grouped = defaultdict(list)
    for path in raw_dir.glob("*.mat"):
        match = FILE_RE.match(path.name)
        if not match:
            continue
        d = match.groupdict(); key = (d["If"], int(d["Ip"]), int(d["In"]), int(d["R"]), int(d["N"]), int(d["gap"]), bool(d["nohydro"]))
        grouped[key].append((int(d["seed"]), path))
    rows = []
    for key, files in sorted(grouped.items()):
        fluid, ip10, in10, radius, number, gap, nohydro = key
        path = min(files)[1]
        reader = RawReader(path); nframes = int(reader.info["state"].shape[-1]); reader.close()
        rows.append(dict(If=fluid, Ip=ip10 / 10, In=in10 / 10, radius_ratio=radius, N=number,
                         gap_ratio=gap, dynamics="no_hydrodynamics" if nohydro else "hydrodynamics",
                         n_files=len(files), seeds=json.dumps(sorted(seed for seed, _ in files)),
                         n_frames=nframes, native_dt=0.001, duration=(nframes - 1) * 0.001,
                         initialization="all fish localized in left chamber" if nohydro else "not source-audited here"))
    return rows


def archived_ip4(cache_dir: Path) -> list[dict]:
    runs = []
    for seed in range(40):
        with np.load(cache_dir / f"seed{seed:02d}_stride100_qscan_v4.npz") as data:
            time = np.asarray(data["step"], dtype=float) * 0.001
            q = np.asarray(data["q"], dtype=float)
        runs.append(dict(run_id=f"archived_Ip4_localized_seed{seed}", condition="alignment", ip=4.0,
                         initialization="localized", seed=seed, time=time, q=q, source="archived"))
    return runs


def new_scan(output: Path) -> list[dict]:
    manifest = read_csv(output / "scan_manifest.csv")
    with np.load(output / "scan_q_trajectories.npz") as data:
        time = np.asarray(data["time"], dtype=float)
        qrows = np.asarray(data["q_count_difference"], dtype=float) / N
    if len(manifest) != len(qrows):
        raise RuntimeError("Scan manifest/trajectory mismatch")
    runs = []
    for row, q in zip(manifest, qrows):
        runs.append(dict(run_id=row["run_id"], condition=row["condition"], ip=float(row["ip"]),
                         initialization=row["initialization"], seed=int(row["seed"]), time=time, q=q, source="new"))
    return runs


def window(run: dict, lo: float, hi: float | None = None) -> np.ndarray:
    keep = run["time"] >= lo - 1e-12
    if hi is not None:
        keep &= run["time"] <= hi + 1e-12
    return run["q"][keep]


def relocation_events(q: np.ndarray, dt: float = DT, bound: float = BASIN) -> list[dict]:
    current = 0; last_source = None; events = []
    for index, value in enumerate(q):
        basin = 1 if value >= bound else (-1 if value <= -bound else 0)
        if basin == current and basin:
            last_source = index
        elif basin and current and basin != current:
            events.append(dict(source_basin=current, target_basin=basin,
                               departure_index=last_source, arrival_index=index,
                               duration=(index - last_source) * dt))
            current = basin; last_source = index
        elif basin and not current:
            current = basin; last_source = index
    return events


def seed_boot(values: np.ndarray, seed: int, nboot: int = 5000) -> tuple[float, list[float]]:
    rng = np.random.default_rng(seed)
    draw = rng.choice(values, (nboot, len(values)), replace=True).mean(axis=1)
    return float(np.mean(values)), np.percentile(draw, [2.5, 97.5]).tolist()


def distribution(runs: list[dict], lo: float, hi: float | None = None) -> np.ndarray:
    per_seed = []
    for run in runs:
        q = window(run, lo, hi)
        counts = np.bincount(np.rint((q + 1) * N / 2).astype(int), minlength=N + 1).astype(float)
        per_seed.append(counts / counts.sum())
    return np.mean(per_seed, axis=0)


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    return float(0.5 * np.abs(a - b).sum())


def reflection_symmetrize(probability: np.ndarray) -> np.ndarray:
    return 0.5 * (probability + probability[::-1])


def summarize_group(runs: list[dict], lo: float, hi: float | None, tag: str, boot_seed: int) -> dict:
    probability = independent_probability()
    values = [] ; all_events = []
    for run in runs:
        q = window(run, lo, hi)
        lag = max(1, round(10 / DT))
        events = relocation_events(q)
        all_events.extend(dict(run_id=run["run_id"], **event) for event in events)
        values.append([
            np.mean(q), np.mean(np.abs(q)), np.mean(q * q), np.mean(np.abs(q) >= STRONG), np.mean(q >= STRONG), np.mean(q <= -STRONG),
            np.mean(np.abs(q[lag:] - q[:-lag])) if len(q) > lag else np.nan,
            len(events),
        ])
    values = np.asarray(values)
    names = ("mean_q", "mean_abs_q", "mean_q2", "strong_fraction", "right_strong_fraction", "left_strong_fraction",
             "lag10_abs_delta_q", "relocations_per_run")
    out = dict(window=tag, n_runs=len(runs), n_samples_per_run_min=min(len(window(r, lo, hi)) for r in runs),
               exposure_time=float(sum(max(0, len(window(r, lo, hi)) - 1) * DT for r in runs)),
               relocation_count=len(all_events), seeds_with_relocation=len({r["run_id"] for r in all_events}),
               relocation_duration_median=float(np.median([r["duration"] for r in all_events])) if all_events else None,
               runs_without_relocation=sum(not relocation_events(window(r, lo, hi)) for r in runs),
               right_censored_terminal_intervals=len(runs))
    for index, name in enumerate(names):
        estimate, ci = seed_boot(values[:, index], boot_seed + index)
        out[name] = estimate; out[name + "_ci95"] = json.dumps(ci)
    empirical = distribution(runs, lo, hi)
    out["tv_from_binomial"] = total_variation(empirical, probability)
    out["reflection_symmetrized_tv_from_binomial"] = total_variation(reflection_symmetrize(empirical), probability)
    out["variance_q"] = out["mean_q2"] - out["mean_q"] ** 2
    rng = np.random.default_rng(boot_seed + 99)
    picked = values[rng.integers(len(values), size=(5000, len(values)))].mean(axis=1)
    variance_draw = picked[:, 2] - picked[:, 0] ** 2
    out["variance_q_ci95"] = json.dumps(np.percentile(variance_draw, [2.5, 97.5]).tolist())
    out["variance_over_independent"] = out["variance_q"] / (1 / N)
    out["q2_over_independent"] = out["mean_q2"] / (1 / N)
    return out


def grouped_runs(runs: list[dict]) -> dict[tuple[str, float, str], list[dict]]:
    grouped = defaultdict(list)
    for run in runs:
        grouped[(run["condition"], run["ip"], run["initialization"])].append(run)
    return grouped


def group_summaries(runs: list[dict]) -> list[dict]:
    rows = []
    for count, (key, part) in enumerate(sorted(grouped_runs(runs).items()), 1):
        condition, ip, initialization = key
        max_time = min(float(r["time"][-1]) for r in part)
        lo = LATE_START; hi = min(300.0, max_time)
        rows.append(dict(condition=condition, Ip=ip, initialization=initialization,
                         **summarize_group(part, lo, hi, f"{lo:g}-{hi:g}", 20261130 + 20 * count)))
    return rows


def stationarity_audit(runs: list[dict]) -> list[dict]:
    rows = []
    groups = grouped_runs(runs)
    for count, (key, part) in enumerate(sorted(groups.items()), 1):
        condition, ip, initialization = key
        end = min(float(r["time"][-1]) for r in part)
        windows = [(0., 100.), (100., 200.), (200., min(300., end))]
        if end >= 1000 - 1e-9:
            windows += [(200., 1000.), (200., 600.), (600., 1000.)]
        for index, (lo, hi) in enumerate(windows):
            if hi <= lo:
                continue
            rows.append(dict(condition=condition, Ip=ip, initialization=initialization,
                             **summarize_group(part, lo, hi, f"{lo:g}-{hi:g}", 20261210 + 40 * count + index)))
        if end >= 300:
            early = distribution(part, 100, 200); late = distribution(part, 200, 300)
            for row in rows[::-1]:
                if row["condition"] == condition and row["Ip"] == ip and row["initialization"] == initialization and row["window"] == "200-300":
                    row["tv_previous_window"] = total_variation(early, late)
                    row["reflection_symmetrized_tv_previous_window"] = total_variation(reflection_symmetrize(early), reflection_symmetrize(late)); break
        if end >= 1000 - 1e-9:
            a = distribution(part, 200, 600); b = distribution(part, 600, 1000)
            for row in rows[::-1]:
                if row["condition"] == condition and row["Ip"] == ip and row["initialization"] == initialization and row["window"] == "600-1000":
                    row["tv_previous_window"] = total_variation(a, b)
                    row["reflection_symmetrized_tv_previous_window"] = total_variation(reflection_symmetrize(a), reflection_symmetrize(b)); break
    return rows


def initialization_sensitivity(runs: list[dict]) -> list[dict]:
    groups = grouped_runs(runs); rows = []
    for condition, ip in sorted({(key[0], key[1]) for key in groups}):
        left = groups.get((condition, ip, "localized")); balanced = groups.get((condition, ip, "balanced"))
        if not left or not balanced:
            continue
        lp = distribution(left, 200, 300); bp = distribution(balanced, 200, 300)
        lq2 = np.asarray([np.mean(window(r, 200, 300) ** 2) for r in left]); bq2 = np.asarray([np.mean(window(r, 200, 300) ** 2) for r in balanced])
        ls = np.asarray([np.mean(np.abs(window(r, 200, 300)) >= STRONG) for r in left]); bs = np.asarray([np.mean(np.abs(window(r, 200, 300)) >= STRONG) for r in balanced])
        rng = np.random.default_rng(20261220 + int(ip * 10) + (100 if condition == "independent" else 0))
        q2_draw = rng.choice(lq2, (5000, len(lq2)), replace=True).mean(axis=1) - rng.choice(bq2, (5000, len(bq2)), replace=True).mean(axis=1)
        strong_draw = rng.choice(ls, (5000, len(ls)), replace=True).mean(axis=1) - rng.choice(bs, (5000, len(bs)), replace=True).mean(axis=1)
        tv_draw = []; sym_tv_draw = []
        for _ in range(2000):
            lr = [left[i] for i in rng.integers(len(left), size=len(left))]
            br = [balanced[i] for i in rng.integers(len(balanced), size=len(balanced))]
            ldraw=distribution(lr,200,300);bdraw=distribution(br,200,300)
            tv_draw.append(total_variation(ldraw,bdraw))
            sym_tv_draw.append(total_variation(reflection_symmetrize(ldraw),reflection_symmetrize(bdraw)))
        rows.append(dict(condition=condition, Ip=ip, n_localized=len(left), n_balanced=len(balanced),
                         distribution_tv=total_variation(lp, bp), reflection_symmetrized_distribution_tv=total_variation(reflection_symmetrize(lp),reflection_symmetrize(bp)),
                         mean_q2_localized=float(lq2.mean()), mean_q2_balanced=float(bq2.mean()),
                         strong_localized=float(ls.mean()), strong_balanced=float(bs.mean()),
                         q2_difference=float(lq2.mean() - bq2.mean()), q2_difference_ci95=json.dumps(np.percentile(q2_draw,[2.5,97.5]).tolist()),
                         strong_difference=float(ls.mean() - bs.mean()), strong_difference_ci95=json.dumps(np.percentile(strong_draw,[2.5,97.5]).tolist()),
                         distribution_tv_bootstrap_ci95=json.dumps(np.percentile(tv_draw,[2.5,97.5]).tolist()),
                         reflection_symmetrized_distribution_tv_bootstrap_ci95=json.dumps(np.percentile(sym_tv_draw,[2.5,97.5]).tolist())))
    return rows


def relocation_rows(runs: list[dict]) -> list[dict]:
    rows = []
    for run in runs:
        for event_id, event in enumerate(relocation_events(run["q"]), 1):
            rows.append(dict(run_id=run["run_id"], condition=run["condition"], Ip=run["ip"], initialization=run["initialization"],
                             seed=run["seed"], event_id=event_id, departure_time=run["time"][event["departure_index"]],
                             arrival_time=run["time"][event["arrival_index"]], **event))
    return rows


def figure(runs: list[dict], summary: list[dict], output: Path) -> None:
    prob = independent_probability(); groups = grouped_runs(runs)
    conditions = [("independent", 0.0)] + [("alignment", x) for x in (0., 1., 2., 4., 8.)]
    fig, axes = plt.subplots(3, 2, figsize=(11, 11))
    colors = {"localized": "#b2182b", "balanced": "#2166ac"}
    for row, (condition, ip) in enumerate(conditions):
        ax = axes[0, 0]
        for initialization, marker in (("localized", "o"), ("balanced", "s")):
            part = groups.get((condition, ip, initialization))
            if not part: continue
            x = -1 if condition == "independent" else ip
            val = np.asarray([np.mean(window(r, 200, 300) ** 2) for r in part]); mean, ci = seed_boot(val, 20261230 + row)
            ax.errorbar(x + (-.06 if initialization == "localized" else .06), mean, yerr=[[mean-ci[0]],[ci[1]-mean]],
                        fmt=marker, color=colors[initialization], capsize=3, label=initialization if row == 0 else None)
    axes[0, 0].axhline(1/N, color="k", ls=":", label="independent exact"); axes[0, 0].set(xticks=[-1,0,1,2,4,8], xticklabels=["social\noff","0","1","2","4","8"], ylabel=r"$\langle q^2\rangle$", title="Macroscopic localization"); axes[0,0].legend(frameon=False)
    for metric, ax, title in (("strong_fraction", axes[0,1], r"Fraction $|q|\geq0.8$"), ("lag10_abs_delta_q", axes[1,0], r"Population mobility $\langle|\Delta q|\rangle_{10}$")):
        for initialization in ("localized", "balanced"):
            selected = [r for r in summary if r["initialization"] == initialization]
            x = [-1 if r["condition"] == "independent" else r["Ip"] for r in selected]
            y = [r[metric] for r in selected]
            ci = np.asarray([json.loads(r[metric + "_ci95"]) for r in selected]); x=np.asarray(x); y=np.asarray(y); order=np.argsort(x)
            error=np.vstack((y-ci[:,0],ci[:,1]-y))
            ax.errorbar(x[order],y[order],yerr=error[:,order],fmt="o-",capsize=3,color=colors[initialization],label=initialization)
        ax.set(xticks=[-1,0,1,2,4,8], xticklabels=["off","0","1","2","4","8"], title=title); ax.legend(frameon=False)
    ax = axes[1,1]
    archived = groups[("alignment",4.0,"localized")]
    for lo, hi, color in ((0,100,"#f4a582"),(200,600,"#92c5de"),(600,1000,"#2166ac")):
        ax.plot(Q_VALUES, distribution(archived,lo,hi), color=color, label=f"Ip=4, t={lo}-{hi}")
    ax.plot(Q_VALUES, prob, "k:", lw=2, label="independent binomial"); ax.set(xlabel="q", ylabel="P(q)", title="Archived Ip=4 initialization/stationarity"); ax.legend(frameon=False, fontsize=8)
    for col, initialization in enumerate(("localized","balanced")):
        ax = axes[2,col]
        for condition, ip in conditions:
            part = groups.get((condition,ip,initialization))
            if not part: continue
            label="social off" if condition=="independent" else f"Ip={ip:g}"
            ax.plot(Q_VALUES, distribution(part,200,300), label=label)
        ax.plot(Q_VALUES, prob,"k:",lw=1.5,label="binomial"); ax.set(xlabel="q",ylabel="P(q)",title=f"Late occupancy: {initialization}"); ax.legend(frameon=False,fontsize=7,ncol=2)
    fig.tight_layout(); fig.savefig(output, dpi=220); plt.close(fig)


def self_test() -> None:
    p = independent_probability(); assert abs(p.sum() - 1) < 1e-12
    assert abs(np.sum(Q_VALUES * p)) < 1e-12
    assert abs(np.sum(Q_VALUES**2 * p) - 1/N) < 1e-12
    q = np.r_[np.full(5,-.8), [-.4,0,.4], np.full(5,.8), [.4,0,-.4], np.full(5,-.8)]
    events = relocation_events(q); assert len(events) == 2 and all(e["duration"] > 0 for e in events)
    print("PASS: exact binomial moments and opposite-basin relocation state machine")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=RAW_DEFAULT)
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DEFAULT)
    parser.add_argument("--output", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    args.output.mkdir(parents=True, exist_ok=True)
    inventory_rows = inventory(args.raw_dir); write_csv(args.output / "existing_data_inventory.csv", inventory_rows)
    runs = archived_ip4(args.cache_dir) + new_scan(args.output)
    summaries = group_summaries(runs); stationarity = stationarity_audit(runs); sensitivity = initialization_sensitivity(runs); relocations = relocation_rows(runs)
    write_csv(args.output / "macroscopic_summary.csv", summaries)
    write_csv(args.output / "stationarity_windows.csv", stationarity)
    write_csv(args.output / "initialization_sensitivity.csv", sensitivity)
    write_csv(args.output / "complete_relocations.csv", relocations)
    p = independent_probability()
    write_csv(args.output / "independent_binomial_benchmark.csv", [dict(n_right=n, q=Q_VALUES[n], probability=p[n]) for n in range(N+1)])
    figure(runs, summaries, args.output / "alignment_localization_regime.png")
    relevant = [r for r in inventory_rows if r["dynamics"] == "no_hydrodynamics" and r["gap_ratio"] == 20 and r["N"] == 100]
    summary = dict(
        exact_independent_benchmark=dict(N=N, mean_q=float(np.sum(Q_VALUES*p)), variance_q=float(np.sum(Q_VALUES**2*p)),
                                         strong_fraction=float(np.sum(p[np.abs(Q_VALUES)>=STRONG])), law="n_R~Binomial(N,1/2), q=(2n_R-N)/N"),
        existing_relevant_inventory=relevant, analysis_window="equal-seed occupancy over t=200--300 for scan comparability",
        archived_ip4_long_window="separate 200--600 versus 600--1000 stationarity audit",
        macroscopic_summary=summaries, initialization_sensitivity=sensitivity,
        relocation_definition="first arrival in opposite |q|>=0.6 dominance basin after last sample in source basin",
        relocation_count=len(relocations), relocation_seeds=len({r["run_id"] for r in relocations}),
        caveat="occupancy is not called stationary where adjacent-window TV or initialization sensitivity remains substantial; no effective potential or switching rate is fitted",
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
