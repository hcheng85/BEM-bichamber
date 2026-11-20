#!/usr/bin/env python3
"""Test a storage-excursion renewal description of collective transfer.

An attempt is a complete excursion of the established total weighted approach
storage away from zero.  The definition has no optimized threshold and detects
the leading portion of an extended school.  Existing 0.5-time trajectory
features and 0.1-time population indices are read only; no model is simulated.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kstest, spearmanr


DT = 0.5
TRANSIENT = 200.0
N = 100
TRANSFER_THRESHOLD = 5
OUTCOMES = ("no_entry", "returned", "partial", "sustained")
TARGETS = ("waiting_time", "next_transfer", "next_first_entry",
           "next_partial", "next_sustained", "same_direction")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def load_features(path: Path) -> dict[int, dict[str, np.ndarray]]:
    with np.load(path) as z:
        arrays = {key: np.asarray(z[key]) for key in z.files}
    result = {}
    for seed in range(40):
        take = arrays["seed"] == seed
        order = np.argsort(arrays["index"][take])
        result[seed] = {key: value[take][order] for key, value in arrays.items()
                        if key not in ("seed", "index")}
    return result


def classify(forward: int, transfer: int) -> str:
    if forward == 0:
        return "no_entry"
    if transfer <= 0:
        return "returned"
    if transfer < TRANSFER_THRESHOLD:
        return "partial"
    return "sustained"


def extract_attempts(features: dict[int, dict[str, np.ndarray]], stride: int = 1) -> list[dict]:
    rows = []
    for seed, a in features.items():
        sampled = np.arange(0, len(a["q"]), stride, dtype=int)
        active = (a["b_left"][sampled] + a["b_right"][sampled]) > 1e-14
        starts = np.flatnonzero(active & np.r_[True, ~active[:-1]])
        stops = np.flatnonzero(active & np.r_[~active[1:], True])
        complete = (starts > 0) & (stops < len(sampled) - 1)
        seed_rows = []
        for attempt_id, (js, je) in enumerate(zip(starts[complete], stops[complete])):
            start, stop = int(sampled[js]), int(sampled[je])
            pre, reset = int(sampled[js - 1]), int(sampled[je + 1])
            # The initiating side is the side with greater storage at the first
            # occupied sample.  No outcome or future direction is used.
            direction = 1 if a["b_left"][start] > a["b_right"][start] else -1
            nright_pre = int(round((float(a["q"][pre]) + 1) * N / 2))
            nright_reset = int(round((float(a["q"][reset]) + 1) * N / 2))
            transfer = int(direction * (nright_reset - nright_pre))
            forward_series = a["plus"] if direction == 1 else a["minus"]
            backward_series = a["minus"] if direction == 1 else a["plus"]
            forward = int(np.sum(forward_series[pre:reset]))
            backward = int(np.sum(backward_series[pre:reset]))
            crossing_bins = np.flatnonzero(forward_series[pre:reset] > 0)
            first_entry_upper = (TRANSIENT + (pre + crossing_bins[0] + 1) * DT
                                 if len(crossing_bins) else np.nan)
            outcome = classify(forward, transfer)
            qpre = float(a["q"][pre])
            majority_direction = 1 if qpre < 0 else (-1 if qpre > 0 else 0)
            seed_rows.append(dict(
                seed=seed, attempt_id=attempt_id, resolution=DT * stride,
                pre_index=pre, start_index=start, end_index=stop, reset_index=reset,
                onset_lower=TRANSIENT + pre * DT,
                onset_upper=TRANSIENT + start * DT,
                reset_time=TRANSIENT + reset * DT,
                duration_bracket=(reset - pre) * DT,
                direction=direction, q_pre=qpre, q_reset=float(a["q"][reset]),
                b_left_onset=float(a["b_left"][start]),
                b_right_onset=float(a["b_right"][start]),
                n_right_pre=nright_pre, n_right_reset=nright_reset,
                population_transfer=transfer, forward_crossings=forward,
                backward_crossings=backward, gross_crossings=forward + backward,
                first_entry=int(forward > 0), first_entry_time_upper=first_entry_upper,
                returned=int(forward > 0 and transfer <= 0),
                partial=int(0 < transfer < TRANSFER_THRESHOLD),
                sustained=int(transfer >= TRANSFER_THRESHOLD), outcome=outcome,
                majority_source=int(direction == majority_direction),
                max_b_left=float(np.max(a["b_left"][start:stop + 1])),
                max_b_right=float(np.max(a["b_right"][start:stop + 1])),
                max_b_total=float(np.max(a["b_left"][start:stop + 1] +
                                         a["b_right"][start:stop + 1])),
                crossing_conservation_residual=transfer - (forward - backward)))
        for current, following in zip(seed_rows[:-1], seed_rows[1:]):
            current["recovery_wait"] = following["onset_upper"] - current["reset_time"]
            current["cycle_time"] = following["onset_upper"] - current["onset_upper"]
        if seed_rows:
            seed_rows[-1]["recovery_wait"] = np.nan
            seed_rows[-1]["cycle_time"] = np.nan
        rows.extend(seed_rows)
    return rows


def transitions(attempts: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for row in attempts:
        by[row["seed"]].append(row)
    result = []
    for seed, seq in by.items():
        seq.sort(key=lambda r: r["attempt_id"])
        for current, following in zip(seq[:-1], seq[1:]):
            result.append(dict(
                seed=seed, attempt_id=current["attempt_id"],
                previous_outcome=current["outcome"], next_outcome=following["outcome"],
                waiting_time=current["recovery_wait"], cycle_time=current["cycle_time"],
                previous_transfer=current["population_transfer"],
                next_transfer=following["population_transfer"],
                next_first_entry=following["first_entry"],
                next_partial=following["partial"],
                next_sustained=following["sustained"],
                same_direction=int(current["direction"] == following["direction"])))
    return result


def cluster_boot(rows: list[dict], metric, seed: int, nboot: int = 4000) -> dict:
    by = defaultdict(list)
    for row in rows:
        by[int(row["seed"])].append(row)
    seeds = sorted(by); point = float(metric(rows)); rng = np.random.default_rng(seed); draw = []
    for _ in range(nboot):
        sample = []
        for replicate, chosen in enumerate(rng.choice(seeds, len(seeds), replace=True)):
            for row in by[int(chosen)]:
                copied = dict(row); copied["seed"] = replicate; sample.append(copied)
        value = metric(sample)
        if np.isfinite(value): draw.append(float(value))
    return dict(estimate=point, ci95=np.percentile(draw, [2.5, 97.5]).tolist(),
                n=len(rows), n_seeds=len(seeds))


def serial_pairs(rows: list[dict], field: str) -> list[dict]:
    by = defaultdict(list)
    for row in rows: by[int(row["seed"])].append(row)
    out = []
    for seed, seq in by.items():
        seq.sort(key=lambda r: r["attempt_id"])
        values = [r[field] for r in seq if np.isfinite(float(r[field]))]
        out.extend(dict(seed=seed, x=float(x), y=float(y)) for x, y in zip(values[:-1], values[1:]))
    return out


def corr_metric(rows: list[dict]) -> float:
    x = np.asarray([r["x"] for r in rows]); y = np.asarray([r["y"] for r in rows])
    return float(spearmanr(x, y).statistic) if np.std(x) and np.std(y) else np.nan


def direction_excess(rows: list[dict]) -> float:
    by = defaultdict(list)
    for row in rows: by[int(row["seed"])].append(row)
    values = []
    for seq in by.values():
        if len(seq) < 2: continue
        direction = np.asarray([r["direction"] for r in seq]); p = np.mean(direction == 1)
        observed = np.mean(direction[:-1] == direction[1:]); values.append(observed - (p*p + (1-p)**2))
    return float(np.mean(values))


def fit_models(train: list[dict]) -> tuple[dict, dict]:
    ordinary = {field: float(np.mean([r[field] for r in train])) for field in TARGETS}
    semi = {}
    for outcome in OUTCOMES:
        part = [r for r in train if r["previous_outcome"] == outcome]
        semi[outcome] = {field: float(np.mean([r[field] for r in part])) for field in TARGETS}
        semi[outcome]["n_train"] = len(part)
    return ordinary, semi


def kernel_table(rows: list[dict]) -> list[dict]:
    result = []
    for cohort, pred in (("train", lambda s:s<30), ("heldout", lambda s:s>=30)):
        selected=[r for r in rows if pred(r["seed"])]
        for previous in OUTCOMES:
            origin=[r for r in selected if r["previous_outcome"]==previous]
            for following in OUTCOMES:
                cell=[r for r in origin if r["next_outcome"]==following]
                result.append(dict(cohort=cohort,previous_outcome=previous,next_outcome=following,
                                   n=len(cell),probability=(len(cell)/len(origin) if origin else np.nan),
                                   mean_wait=(float(np.mean([r["waiting_time"] for r in cell])) if cell else np.nan),
                                   median_wait=(float(np.median([r["waiting_time"] for r in cell])) if cell else np.nan),
                                   mean_next_transfer=(float(np.mean([r["next_transfer"] for r in cell])) if cell else np.nan)))
    return result


def heldout_scores(test: list[dict], ordinary: dict, semi: dict) -> tuple[list[dict], list[dict]]:
    per_seed = []
    for seed in sorted({int(r["seed"]) for r in test}):
        part = [r for r in test if int(r["seed"]) == seed]
        for target in TARGETS:
            actual = np.asarray([r[target] for r in part], float)
            p0 = np.full(len(part), ordinary[target])
            p1 = np.asarray([semi[r["previous_outcome"]][target] for r in part])
            per_seed.append(dict(seed=seed, target=target, n=len(part),
                                 ordinary_mse=float(np.mean((actual-p0)**2)),
                                 semi_markov_mse=float(np.mean((actual-p1)**2)),
                                 improvement=float(np.mean((actual-p0)**2-(actual-p1)**2))))
    summary = []
    for target in TARGETS:
        part = [r for r in per_seed if r["target"] == target]
        values = np.asarray([r["improvement"] for r in part]); rng = np.random.default_rng(20261201 + TARGETS.index(target))
        pick = rng.integers(len(values), size=(10000, len(values)))
        summary.append(dict(target=target, n_test=sum(r["n"] for r in part), n_test_seeds=len(part),
                            ordinary_mse=float(np.mean([r["ordinary_mse"] for r in part])),
                            semi_markov_mse=float(np.mean([r["semi_markov_mse"] for r in part])),
                            mse_improvement=float(values.mean()),
                            improvement_ci95=json.dumps(np.percentile(values[pick].mean(axis=1), [2.5,97.5]).tolist())))
    return per_seed, summary


def residual_memory(test: list[dict], semi: dict) -> dict:
    by = defaultdict(list)
    for row in test: by[int(row["seed"])].append(row)
    pairs = []
    for seed, seq in by.items():
        seq.sort(key=lambda r: r["attempt_id"])
        residual = [r["waiting_time"] - semi[r["previous_outcome"]]["waiting_time"] for r in seq]
        pairs.extend(dict(seed=seed, x=x, y=y) for x, y in zip(residual[:-1], residual[1:]))
    return cluster_boot(pairs, corr_metric, 20261220)


def per_seed_aggregate(attempts: list[dict], features: dict[int,dict[str,np.ndarray]]) -> list[dict]:
    rows = []
    for seed in sorted({int(r["seed"]) for r in attempts}):
        part = [r for r in attempts if int(r["seed"]) == seed]
        exposure = (len(features[seed]["q"])-1)*DT
        rows.append(dict(seed=seed, n_attempts=len(part), exposure=exposure,
                         attempt_rate=100*len(part)/exposure,
                         first_entry_probability=float(np.mean([r["first_entry"] for r in part])),
                         partial_probability=float(np.mean([r["partial"] for r in part])),
                         sustained_probability=float(np.mean([r["sustained"] for r in part])),
                         mean_population_transfer=float(np.mean([r["population_transfer"] for r in part])),
                         mean_gross_crossings=float(np.mean([r["gross_crossings"] for r in part]))))
    return rows


def aggregate_validation(seed_rows: list[dict]) -> list[dict]:
    train = [r for r in seed_rows if r["seed"] < 30]; test = [r for r in seed_rows if r["seed"] >= 30]
    result = []
    fields = ("attempt_rate", "first_entry_probability", "partial_probability",
              "sustained_probability", "mean_population_transfer", "mean_gross_crossings")
    for j, field in enumerate(fields):
        x = np.asarray([r[field] for r in train]); y = np.asarray([r[field] for r in test]); difference = y.mean()-x.mean()
        rng = np.random.default_rng(20261240+j); draw = []
        for _ in range(10000):
            draw.append(rng.choice(y,len(y),replace=True).mean()-rng.choice(x,len(x),replace=True).mean())
        result.append(dict(metric=field, train_prediction=float(x.mean()), heldout_observed=float(y.mean()),
                           heldout_minus_train=float(difference), difference_ci95=json.dumps(np.percentile(draw,[2.5,97.5]).tolist()),
                           n_train_seeds=len(x), n_test_seeds=len(y)))
    return result


def relocation_events(cache_dir: Path, attempts: list[dict], bound: float = .6) -> list[dict]:
    by = defaultdict(list)
    for row in attempts: by[int(row["seed"])].append(row)
    result = []
    for seed in range(40):
        with np.load(cache_dir/f"seed{seed:02d}_stride100_qscan_v4.npz") as z:
            step=np.asarray(z["step"]); q=np.asarray(z["q"])
        take=step>=round(TRANSIENT/.001); step,q=step[take],q[take]
        previous=0; last_index=None; event_id=0
        for i,value in enumerate(q):
            state=1 if value>=bound else (-1 if value<=-bound else 0)
            if state:
                if previous and state!=previous:
                    lo=int(last_index); segment=q[lo:i+1]; center=lo+int(np.argmin(np.abs(segment)))
                    t0=float(step[lo]*.001); tc=float(step[center]*.001); t1=float(step[i]*.001)
                    overlap=[r for r in by[seed] if not (r["reset_time"]<t0 or r["onset_lower"]>t1)]
                    center_attempt=[r for r in overlap if r["onset_lower"]<=tc<=r["reset_time"]]
                    result.append(dict(seed=seed, relocation_id=event_id, source_basin=previous, target_basin=state,
                                       departure_time=t0, center_time=tc, arrival_time=t1, duration=t1-t0,
                                       n_overlapping_attempts=len(overlap), center_attempt_outcome=(center_attempt[0]["outcome"] if center_attempt else ""),
                                       center_attempt_transfer=(center_attempt[0]["population_transfer"] if center_attempt else "")))
                    event_id+=1
                previous=state; last_index=i
    return result


def diagnostics(attempts: list[dict], trans: list[dict], semi: dict) -> dict:
    result = {}
    for name, pred in (("train",lambda s:s<30),("heldout",lambda s:s>=30),("all",lambda s:True)):
        a=[r for r in attempts if pred(r["seed"])]; t=[r for r in trans if pred(r["seed"])]
        waits=np.asarray([r["waiting_time"] for r in t]); mu=float(waits.mean())
        wait_pairs=serial_pairs(t,"waiting_time"); transfer_pairs=serial_pairs(a,"population_transfer")
        result[name]=dict(n_attempts=len(a),n_transitions=len(t),outcomes=Counter(r["outcome"] for r in a),
                          waiting_mean=mu,waiting_median=float(np.median(waits)),waiting_cv=float(waits.std()/mu),
                          exponential_ks=float(kstest(waits,"expon",args=(0,mu)).statistic),
                          lag1_wait=cluster_boot(wait_pairs,corr_metric,20261260),
                          lag1_transfer=cluster_boot(transfer_pairs,corr_metric,20261261),
                          direction_persistence_excess=cluster_boot(a,direction_excess,20261262))
    result["heldout_semi_markov_wait_residual_lag1"]=residual_memory([r for r in trans if r["seed"]>=30],semi)
    return result


def sensitivity(primary: list[dict], coarse: list[dict]) -> dict:
    def compact(rows):
        trans=transitions(rows);wait_pairs=serial_pairs(trans,"waiting_time")
        return dict(n_attempts=len(rows),median_duration=float(np.median([r["duration_bracket"] for r in rows])),
                    sustained_fraction=float(np.mean([r["sustained"] for r in rows])),
                    lag1_wait_spearman=corr_metric(wait_pairs),direction_persistence_excess=direction_excess(rows))
    return {"dt_0p5":compact(primary),"dt_1p0":compact(coarse)}


def figure(attempts, trans, score_summary, relocations, out: Path) -> None:
    fig,axes=plt.subplots(2,2,figsize=(10,7.5))
    ax=axes[0,0]
    for name,pred,color in (("train",lambda s:s<30,"#327C91"),("held-out",lambda s:s>=30,"#B66A50")):
        w=np.sort([r["waiting_time"] for r in trans if pred(r["seed"])]);surv=1-np.arange(1,len(w)+1)/len(w);ax.step(w,surv,where="post",label=name,color=color)
    mean=np.mean([r["waiting_time"] for r in trans if r["seed"]<30]);x=np.linspace(0,20,200);ax.plot(x,np.exp(-x/mean),"k:",label="exponential, train mean");ax.set(xlabel="recovery waiting time",ylabel="survival",title="Waiting times are not exponential");ax.legend(frameon=False)
    ax=axes[0,1];x=np.arange(len(OUTCOMES));width=.36
    for offset,(name,pred,color) in zip((-.18,.18),(("train",lambda s:s<30,"#327C91"),("held-out",lambda s:s>=30,"#B66A50"))):
        values=[np.mean([r["waiting_time"] for r in trans if pred(r["seed"]) and r["previous_outcome"]==o]) for o in OUTCOMES];ax.bar(x+offset,values,width,label=name,color=color)
    ax.set_xticks(x,OUTCOMES,rotation=20);ax.set(ylabel="mean recovery wait",title="Outcome-conditioned recovery");ax.legend(frameon=False)
    ax=axes[1,0];names=[r["target"].replace("next_","").replace("same_direction","direction") for r in score_summary];values=[r["mse_improvement"] for r in score_summary];lo=[];hi=[]
    for row,value in zip(score_summary,values):
        ci=json.loads(row["improvement_ci95"]);lo.append(value-ci[0]);hi.append(ci[1]-value)
    ax.bar(np.arange(len(names)),values,yerr=[lo,hi],capsize=3);ax.axhline(0,color="k",ls=":");ax.set_xticks(np.arange(len(names)),names,rotation=30);ax.set(ylabel="held-out MSE/Brier improvement",title="Previous-outcome semi-Markov vs renewal")
    ax=axes[1,1];ax.scatter([r["duration"] for r in relocations],[r["n_overlapping_attempts"] for r in relocations],s=35);ax.set(xlabel="opposite-basin relocation duration",ylabel="overlapping approach attempts",title="Complete relocations are multi-attempt")
    fig.tight_layout();fig.savefig(out/"approach_renewal_validation.png",dpi=220);plt.close(fig)


def self_test() -> None:
    # Two complete total-storage excursions; the first transfers five fish.
    features={0:dict(q=np.asarray([-.8,-.8,-.7,-.7,-.8,-.8,-.8]),b_left=np.asarray([0,.1,.2,0,0,.1,0]),b_right=np.zeros(7),plus=np.asarray([0,3,2,0,0,0,0]),minus=np.asarray([0,0,0,5,0,0,0]))}
    rows=extract_attempts(features);assert len(rows)==2 and rows[0]["sustained"] and rows[0]["crossing_conservation_residual"]==0 and rows[1]["outcome"]=="no_entry"
    print("PASS: threshold-free storage excursions, transfer conservation, and outcome stages")


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("--features",type=Path);parser.add_argument("--q-cache-dir",type=Path);parser.add_argument("--output",type=Path);parser.add_argument("--self-test",action="store_true");args=parser.parse_args()
    if args.self_test:self_test();return
    if not all((args.features,args.q_cache_dir,args.output)):parser.error("--features, --q-cache-dir and --output are required")
    args.output.mkdir(parents=True,exist_ok=True);features=load_features(args.features);attempts=extract_attempts(features);trans=transitions(attempts);train=[r for r in trans if r["seed"]<30];test=[r for r in trans if r["seed"]>=30];ordinary,semi=fit_models(train);kernel=kernel_table(trans);per_seed_scores,score_summary=heldout_scores(test,ordinary,semi);seed_rows=per_seed_aggregate(attempts,features);aggregate=aggregate_validation(seed_rows);relocations=relocation_events(args.q_cache_dir,attempts);diag=diagnostics(attempts,trans,semi);coarse=extract_attempts(features,2);sens=sensitivity(attempts,coarse)
    write_csv(args.output/"attempt_episodes.csv",attempts);write_csv(args.output/"attempt_transitions.csv",trans);write_csv(args.output/"semi_markov_kernel.csv",kernel);write_csv(args.output/"heldout_score_by_seed.csv",per_seed_scores);write_csv(args.output/"heldout_prediction_summary.csv",score_summary);write_csv(args.output/"attempt_statistics_by_seed.csv",seed_rows);write_csv(args.output/"heldout_aggregate_validation.csv",aggregate);write_csv(args.output/"opposite_basin_relocations.csv",relocations);figure(attempts,trans,score_summary,relocations,args.output)
    summary=dict(configuration=dict(attempt_definition="complete excursion of b_left+b_right from zero back to zero",feature_dt=DT,transfer_threshold_fish=TRANSFER_THRESHOLD,train_seeds="0-29",heldout_seeds="30-39",new_simulations=0,fitted_potentials=0),audits=dict(max_abs_crossing_conservation_residual=max(abs(r["crossing_conservation_residual"]) for r in attempts),n_complete_attempts=len(attempts),n_transitions=len(trans),n_relocations=len(relocations),n_both_sides_occupied_at_onset=sum(r["b_left_onset"]>0 and r["b_right_onset"]>0 for r in attempts),n_majority_source_attempts=sum(r["majority_source"] for r in attempts),n_equal_population_onsets=sum(r["q_pre"]==0 for r in attempts)),ordinary_renewal_diagnostics=diag,semi_markov_model=dict(states=list(OUTCOMES),ordinary=ordinary,conditional=semi,heldout_scores=score_summary),heldout_aggregate_validation=aggregate,resolution_sensitivity=sens,relocations=dict(n=len(relocations),n_seeds=len({r["seed"] for r in relocations}),median_duration=float(np.median([r["duration"] for r in relocations])),median_overlapping_attempts=float(np.median([r["n_overlapping_attempts"] for r in relocations])),minimum_overlapping_attempts=min(r["n_overlapping_attempts"] for r in relocations),maximum_overlapping_attempts=max(r["n_overlapping_attempts"] for r in relocations)))
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else float(x))+"\n")


if __name__=="__main__":main()
