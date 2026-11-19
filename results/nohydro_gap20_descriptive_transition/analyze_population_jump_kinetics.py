#!/usr/bin/env python3
"""Test elimination of polarization/storage from the gap-20 population process.

The population coordinate is kept on its exact 101-state lattice.  The common
0.1-time occupancy record is used for all 40 trajectories; native occupancy
records, where available, are used only to quantify temporal-resolution loss.
No microscopic dynamics are run or modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from analyze_approach_current_balance import ELL, N, raw_path
from analyze_initiation_propagation import RawReader


DT_Q = 0.1
DT_FEATURE = 0.5
TRANSIENT = 200.0
N_FISH = 100
LAGS = (0.1, 0.5, 1.0, 2.0, 5.0)
Q_EDGES = np.asarray([-1.000001, -.6, -.2, .2, .6, 1.000001])
P_EDGES = np.asarray([-1.000001, -.25, 0., .25, 1.000001])
B_EDGES = np.asarray([0., .025, .075, .15, 1.000001])


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def bin_index(x: float, edges: np.ndarray) -> int:
    return int(np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2))


def load_q(cache_dir: Path, seed: int) -> tuple[np.ndarray, np.ndarray]:
    with np.load(cache_dir / f"seed{seed:02d}_stride100_qscan_v4.npz") as data:
        steps = np.asarray(data["step"], int); q = np.asarray(data["q"], float)
    keep = steps >= round(TRANSIENT / .001)
    return steps[keep], np.rint((q[keep] + 1) * N_FISH / 2).astype(np.int16)


def reflection_rates(series: dict[int, np.ndarray], seeds: list[int]) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    exposure = np.zeros(N_FISH + 1); birth = np.zeros(N_FISH + 1)
    raw_plus = np.zeros(N_FISH + 1); raw_minus = np.zeros(N_FISH + 1)
    raw_exp = np.zeros(N_FISH + 1)
    for seed in seeds:
        n = series[seed]; start = n[:-1]; delta = np.diff(n)
        raw_exp += np.bincount(start, minlength=N_FISH + 1)
        raw_plus += np.bincount(start[delta > 0], weights=delta[delta > 0], minlength=N_FISH + 1)
        raw_minus += np.bincount(start[delta < 0], weights=-delta[delta < 0], minlength=N_FISH + 1)
    # A rightward transfer at n and its reflected leftward transfer at 100-n
    # estimate the same rate without imposing any stationarity beyond symmetry.
    exposure = DT_Q * (raw_exp + raw_exp[::-1])
    birth = raw_plus + raw_minus[::-1]
    lam_plus = np.divide(birth, exposure, out=np.full_like(birth, np.nan), where=exposure > 0)
    lam_minus = lam_plus[::-1]
    rows = [dict(n_right=n, q=2*n/N_FISH-1, exposure_time=exposure[n],
                 rightward_count=birth[n], lambda_right=lam_plus[n], lambda_left=lam_minus[n])
            for n in range(N_FISH + 1) if exposure[n] > 0]
    return lam_plus, lam_minus, rows


def native_resolution_audit(processed_dir: Path, series: dict[int, np.ndarray]) -> dict:
    rows = []
    for seed in range(40):
        path = processed_dir / f"event_histogram_If00_In05_num100_R50_gap20_n{seed}_BEM_no_hydro.mat"
        if not path.exists():
            continue
        nl = np.asarray(loadmat(path, variable_names=["num_fish_left"])["num_fish_left"]).reshape(-1)
        nl = nl[round(TRANSIENT/.001):]
        native = float(np.abs(np.diff(nl)).sum())
        coarse = float(np.abs(np.diff(series[seed])).sum())
        rows.append(dict(seed=seed, native_crossings=native, coarse_crossings=coarse,
                         retained_fraction=coarse/native if native else np.nan))
    return dict(n_native_seeds=len(rows), native_crossings=sum(r["native_crossings"] for r in rows),
                coarse_crossings=sum(r["coarse_crossings"] for r in rows),
                retained_fraction=sum(r["coarse_crossings"] for r in rows)/sum(r["native_crossings"] for r in rows),
                rows=rows)


def empirical_kernel(series: dict[int, np.ndarray], seeds: list[int], lag: int) -> tuple[np.ndarray, np.ndarray]:
    count = np.zeros((N_FISH + 1, N_FISH + 1))
    for seed in seeds:
        n = series[seed]; np.add.at(count, (n[:-lag], n[lag:]), 1)
        # Reflection augmentation enforces the exact left-right model symmetry.
        np.add.at(count, (N_FISH-n[:-lag], N_FISH-n[lag:]), 1)
    support = count.sum(axis=1)
    kernel = np.zeros_like(count)
    ok = support > 0; kernel[ok] = count[ok] / support[ok, None]
    kernel[~ok, np.flatnonzero(~ok)] = 1
    return kernel, support


def weighted_tv(a: np.ndarray, b: np.ndarray, weight: np.ndarray, threshold: int = 30) -> tuple[float, int]:
    ok = weight >= threshold
    if not np.any(ok): return np.nan, 0
    w = weight[ok] / weight[ok].sum()
    return float(np.sum(w * .5 * np.abs(a[ok]-b[ok]).sum(axis=1))), int(ok.sum())


def kernel_tests(series: dict[int, np.ndarray]) -> list[dict]:
    out = []; rng=np.random.default_rng(20261120)
    for tau in LAGS:
        lag = round(tau/DT_Q)
        p, support = empirical_kernel(series, list(range(30)), lag)
        p2, support2 = empirical_kernel(series, list(range(30)), 2*lag)
        ck, nrow = weighted_tv(p2, p@p, support2)
        ok=support2>=30;weight=support2[ok]/support2[ok].sum();null=[];target=(p@p)[ok]
        for _ in range(300):
            row_tv=[]
            for probability,count in zip(target,support2[ok].astype(int)):
                sampled=rng.multinomial(count,probability)/count
                row_tv.append(.5*np.abs(sampled-probability).sum())
            null.append(float(np.dot(weight,row_tv)))
        test, test_support = empirical_kernel(series, list(range(30,40)), lag)
        general, ng = weighted_tv(test, p, test_support)
        test2, test2_support = empirical_kernel(series, list(range(30,40)), 2*lag)
        predicted2, npred = weighted_tv(test2, p@p, test2_support)
        out.append(dict(tau=tau, ck_total_variation=ck, ck_supported_q_states=nrow,
                        ck_sampling_null_ci95=json.dumps(np.percentile(null,[2.5,97.5]).tolist()),
                        ck_to_null_median_ratio=ck/np.median(null),
                        heldout_kernel_total_variation=general, heldout_supported_q_states=ng,
                        heldout_two_step_total_variation=predicted2, heldout_two_step_supported_q_states=npred))
    return out


def table_means(rows: list[dict], key_fields: tuple[str, ...], min_count: int = 20) -> dict[tuple, float]:
    by = defaultdict(list)
    for row in rows: by[tuple(row[k] for k in key_fields)].append(row["future_delta"])
    return {k: float(np.mean(v)) for k,v in by.items() if len(v) >= min_count}


def history_test(series: dict[int, np.ndarray], tau: float) -> dict:
    lag = round(tau/DT_Q); train=[]; test=[]
    for seed,n in series.items():
        target = train if seed < 30 else test
        for i in range(lag, len(n)-lag):
            past = int(np.sign(n[i]-n[i-lag])); future = float((n[i+lag]-n[i])/N_FISH)
            target.append(dict(seed=seed, q=int(n[i]), history=past, future_delta=future))
    qmean = table_means(train,("q",)); hmean = table_means(train,("q","history"))
    stats=[]
    for seed in range(30,40):
        part=[r for r in test if r["seed"]==seed and (r["q"],) in qmean and (r["q"],r["history"]) in hmean]
        if not part: continue
        e0=np.mean([(r["future_delta"]-qmean[(r["q"],)])**2 for r in part])
        e1=np.mean([(r["future_delta"]-hmean[(r["q"],r["history"])])**2 for r in part])
        stats.append((seed,e0,e1,len(part)))
    arr=np.asarray([[x[1],x[2]] for x in stats]); rng=np.random.default_rng(20261118+lag)
    pick=rng.integers(len(arr),size=(3000,len(arr))); gain=arr[:,0]-arr[:,1]; boot=gain[pick].mean(axis=1)
    return dict(tau=tau,n_test=sum(x[3] for x in stats),n_test_seeds=len(stats),q_only_mse=float(arr[:,0].mean()),
                q_history_mse=float(arr[:,1].mean()),history_mse_improvement=float(gain.mean()),
                history_improvement_ci95=json.dumps(np.percentile(boot,[2.5,97.5]).tolist()))


def sample_features(raw_dir: Path, cache_dir: Path, series: dict[int,np.ndarray], cache: Path) -> list[dict]:
    if cache.exists():
        with np.load(cache) as z:
            fields=list(z.files); count=len(z[fields[0]])
            return [{k:z[k][i].item() for k in fields} for i in range(count)]
    out=[]; stride=round(DT_FEATURE/.001); qstride=round(DT_FEATURE/DT_Q)
    for seed in range(40):
        print(f"feature snapshots {seed+1}/40 seed={seed}",flush=True)
        reader=RawReader(raw_path(raw_dir,seed)); n=series[seed]; steps=np.arange(round(TRANSIENT/.001),reader.info["state"].shape[-1]-1,stride,dtype=int)
        length=min(len(steps),len(n[::qstride]))
        try:
            for j,step in enumerate(steps[:length]):
                st=reader.frame("state",int(step)); x,y,theta=st; right=x>=0; e=np.cos(theta)
                px=float(e.mean()); q=float(2*right.sum()/N_FISH-1)
                la=(~right)&((-x)<=ELL)&(np.abs(y)<=ELL/2); ra=right&(x<=ELL)&(np.abs(y)<=ELL/2)
                bl=float(np.sum(np.clip((ELL+x[la])/ELL,0,1))/N_FISH)
                br=float(np.sum(np.clip((ELL-x[ra])/ELL,0,1))/N_FISH)
                i=j*qstride; segment=n[i:min(i+qstride+1,len(n))]; d=np.diff(segment)
                plus=int(d[d>0].sum()); minus=int(-d[d<0].sum())
                out.append(dict(seed=seed,index=j,q=q,px=px,b_left=bl,b_right=br,plus=plus,minus=minus))
        finally: reader.close()
    arrays={k:np.asarray([r[k] for r in out]) for k in out[0]};cache.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(cache,**arrays)
    return out


def directed_rows(features: list[dict]) -> list[dict]:
    out=[]
    for r in features:
        for direction in (-1,1):
            source = r["b_left"] if direction==1 else r["b_right"]
            target = r["b_right"] if direction==1 else r["b_left"]
            count = r["plus"] if direction==1 else r["minus"]
            out.append(dict(seed=r["seed"],qbin=bin_index(direction*r["q"],Q_EDGES),
                            pbin=bin_index(direction*r["px"],P_EDGES),
                            sbin=bin_index(source,B_EDGES),tbin=bin_index(target,B_EDGES),count=count))
    return out


def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    mu=np.maximum(mu,1e-10);term=mu.copy();positive=y>0
    term[positive]=y[positive]*np.log(y[positive]/mu[positive])-(y[positive]-mu[positive])
    return float(2*np.mean(term))


def intensity_models(features: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    rows=directed_rows(features);train=[r for r in rows if r["seed"]<30];test=[r for r in rows if r["seed"]>=30]
    specs=(("q",("qbin",)),("q+px",("qbin","pbin")),("q+source_storage",("qbin","sbin")),
           ("q+px+source_storage",("qbin","pbin","sbin")),
           ("q+px+source+target_storage",("qbin","pbin","sbin","tbin")))
    model_results=[]; detail=[]; fitted={};rng=np.random.default_rng(20261119)
    q_by=defaultdict(list)
    for r in train:q_by[(r["qbin"],)].append(r["count"])
    qmean={k:np.mean(v) for k,v in q_by.items() if len(v)>=50}
    for name,fields in specs:
        by=defaultdict(list)
        for r in train:by[tuple(r[f] for f in fields)].append(r["count"])
        mean={k:float(np.mean(v)) for k,v in by.items() if len(v)>=50}
        fitted[name]=(fields,mean)
        seedstats=[]
        for seed in range(30,40):
            part=[r for r in test if tuple(r[f] for f in fields) in mean and (r["qbin"],) in qmean]
            part=[r for r in part if r["seed"]==seed]
            if not part:continue
            y=np.asarray([r["count"] for r in part]);pred=np.asarray([mean[tuple(r[f] for f in fields)] for r in part]);base=np.asarray([qmean[(r["qbin"],)] for r in part])
            seedstats.append((seed,poisson_deviance(y,base),poisson_deviance(y,pred),len(part)))
        a=np.asarray([[x[1],x[2]] for x in seedstats]);gain=a[:,0]-a[:,1];pick=rng.integers(len(a),size=(3000,len(a)))
        model_results.append(dict(model=name,n_train=len(train),n_test_supported=sum(x[3] for x in seedstats),support_fraction=sum(x[3] for x in seedstats)/len(test),n_cells=len(mean),
                                  q_only_deviance_same_support=float(a[:,0].mean()),model_deviance=float(a[:,1].mean()),deviance_improvement=float(gain.mean()),
                                  improvement_ci95=json.dumps(np.percentile(gain[pick].mean(axis=1),[2.5,97.5]).tolist())))
        for key,value in mean.items():detail.append(dict(model=name,cell=json.dumps(key),n_train=len(by[key]),crossing_intensity=value/DT_FEATURE))
    comparisons=[]
    for reference,candidate in (("q","q+px"),("q","q+source_storage"),("q+source_storage","q+px+source_storage"),("q+px+source_storage","q+px+source+target_storage")):
        rf,rm=fitted[reference];cf,cm=fitted[candidate];stats=[]
        for seed in range(30,40):
            part=[r for r in test if r["seed"]==seed and tuple(r[f] for f in rf) in rm and tuple(r[f] for f in cf) in cm]
            y=np.asarray([r["count"] for r in part]);pr=np.asarray([rm[tuple(r[f] for f in rf)] for r in part]);pc=np.asarray([cm[tuple(r[f] for f in cf)] for r in part])
            stats.append((poisson_deviance(y,pr),poisson_deviance(y,pc),len(part)))
        a=np.asarray([[x[0],x[1]] for x in stats]);gain=a[:,0]-a[:,1];pick=rng.integers(len(a),size=(3000,len(a)))
        comparisons.append(dict(reference=reference,candidate=candidate,n_test_supported=sum(x[2] for x in stats),reference_deviance=float(a[:,0].mean()),candidate_deviance=float(a[:,1].mean()),
                                deviance_improvement=float(gain.mean()),improvement_ci95=json.dumps(np.percentile(gain[pick].mean(axis=1),[2.5,97.5]).tolist())))
    return model_results,detail,comparisons


def fixed_q_relaxation(features: list[dict], series: dict[int,np.ndarray]) -> list[dict]:
    # Subtract exact-q conditional means before correlating, so slow variation
    # of q cannot masquerade as relaxation of a fast variable.
    means={}
    for field in ("px","b_left","b_right"):
        by=defaultdict(list)
        for r in features:by[round((r["q"]+1)*50)].append(r[field])
        means[field]={k:float(np.mean(v)) for k,v in by.items()}
    byseed=defaultdict(list)
    for r in features:byseed[r["seed"]].append(r)
    out=[]
    for tau in LAGS[1:]:
        lag=round(tau/DT_FEATURE); correlations=defaultdict(list);fractions=[]
        for seed,rows in byseed.items():
            qfine=series[seed]; qstride=round(DT_FEATURE/DT_Q);valid=[]
            for i in range(len(rows)-lag):
                lo=i*qstride;hi=(i+lag)*qstride
                valid.append(bool(np.all(qfine[lo:hi+1]==qfine[lo])))
            valid=np.asarray(valid);fractions.append(valid.mean())
            for field in correlations.keys() | {"px","b_left","b_right"}:
                x=np.asarray([r[field]-means[field][round((r["q"]+1)*50)] for r in rows])
                ok=valid & np.isfinite(x[:-lag]) & np.isfinite(x[lag:])
                if ok.sum()>10 and np.std(x[:-lag][ok])>0 and np.std(x[lag:][ok])>0:
                    correlations[field].append(float(np.corrcoef(x[:-lag][ok],x[lag:][ok])[0,1]))
        for field,vals in correlations.items():
            out.append(dict(tau=tau,coordinate=field,mean_seed_fixed_q_autocorrelation=float(np.nanmean(vals)),n_seeds=len(vals),fraction_intervals_no_crossing=float(np.mean(fractions))))
    return out


def current_field(features: list[dict]) -> tuple[list[dict],dict]:
    # Current on the interpretable (q,px) projection.  It is descriptive of
    # nonequilibrium circulation, not claimed to close the dynamics.
    qe=np.linspace(-1.000001,1.000001,11);pe=np.linspace(-1.000001,1.000001,9);shape=(10,8);centers=np.array([(i,j) for i in range(10) for j in range(8)])
    seed_counts=[]
    byseed=defaultdict(list)
    for r in features:byseed[r["seed"]].append(r)
    for seed,rows in byseed.items():
        c=np.zeros((80,80))
        for a,b in zip(rows[:-1],rows[1:]):
            i=bin_index(a["q"],qe)*8+bin_index(a["px"],pe);j=bin_index(b["q"],qe)*8+bin_index(b["px"],pe);c[i,j]+=1
        seed_counts.append(c)
    train=sum(seed_counts[:30]);test=sum(seed_counts[30:]);antisym=train-train.T;testa=test-test.T
    iu=np.triu_indices(80,1);x=antisym[iu];y=testa[iu];mask=(np.abs(x)+np.abs(y))>0
    cosine=float(np.dot(x[mask],y[mask])/(np.linalg.norm(x[mask])*np.linalg.norm(y[mask]))) if np.any(mask) else np.nan
    rows=[]
    qc=(qe[:-1]+qe[1:])/2;pc=(pe[:-1]+pe[1:])/2
    for idx,(i,j) in enumerate(centers):
        vx=vy=0.
        for dst,(k,l) in enumerate(centers):
            flux=antisym[idx,dst];vx+=flux*(qc[k]-qc[i]);vy+=flux*(pc[l]-pc[j])
        rows.append(dict(q_center=qc[i],px_center=pc[j],occupancy=float(train[idx].sum()),current_q=vx,current_px=vy))
    irr=float(np.abs(x).sum()/max(train.sum(),1))
    return rows,dict(train_irreversible_edge_fraction=irr,train_test_edge_current_cosine=cosine,lag=DT_FEATURE)


def make_figure(rate_rows, kernels, histories, relax, currents, output: Path) -> None:
    fig,ax=plt.subplots(2,2,figsize=(10,7.7))
    q=np.asarray([r["q"] for r in rate_rows]);rp=np.asarray([r["lambda_right"] for r in rate_rows]);rm=np.asarray([r["lambda_left"] for r in rate_rows]);
    ax[0,0].plot(q,rp,"o-",ms=3,label="left to right");ax[0,0].plot(q,rm,"o-",ms=3,label="right to left");ax[0,0].set(xlabel="q",ylabel="crossing intensity",title="reflection-symmetric population rates");ax[0,0].legend(frameon=False)
    kt=[r["tau"] for r in kernels];null_hi=[json.loads(r["ck_sampling_null_ci95"])[1] for r in kernels]
    ax[0,1].plot(kt,[r["ck_total_variation"] for r in kernels],"o-",label="CK residual");ax[0,1].plot(kt,null_hi,"k--",label="sampling-null 97.5%")
    ax[0,1].plot(kt,[r["heldout_kernel_total_variation"] for r in kernels],"s-",label="held-out kernel");ax[0,1].set(xscale="log",xlabel="coarse lag",ylabel="weighted total variation",title="q-kernel consistency");ax[0,1].legend(frameon=False,fontsize=8)
    for field,marker in (("px","o"),("b_left","s"),("b_right","^")):
        part=[r for r in relax if r["coordinate"]==field];ax[1,0].plot([r["tau"] for r in part],[r["mean_seed_fixed_q_autocorrelation"] for r in part],marker+"-",label=field)
    ax[1,0].axhline(0,color="k",lw=.6);ax[1,0].set(xlabel="lag while q remains fixed",ylabel="conditional autocorrelation",title="candidate fast-variable relaxation");ax[1,0].legend(frameon=False)
    c=[r for r in currents if r["occupancy"]>=30];scale=max([math.hypot(r["current_q"],r["current_px"]) for r in c]+[1]);ax[1,1].quiver([r["q_center"] for r in c],[r["px_center"] for r in c],[r["current_q"]/scale for r in c],[r["current_px"]/scale for r in c],angles="xy",scale_units="xy",scale=1,color="tab:red")
    ax[1,1].set(xlim=(-1,1),ylim=(-1,1),xlabel="q",ylabel="global axial polarization px",title="empirical probability-current projection")
    fig.suptitle("Population jump kinetics and elimination test",fontsize=11);fig.tight_layout();fig.savefig(output,dpi=220);plt.close(fig)


def self_test() -> None:
    q={0:np.asarray([20,20,21,20,19,19,20],dtype=np.int16)}
    lp,lm,rows=reflection_rates(q,[0]);assert np.isclose(lp[20],lm[80],equal_nan=True)
    assert sum(r["rightward_count"] for r in rows)==4
    p,s=empirical_kernel(q,[0],1);assert np.allclose(p.sum(axis=1),1) and s.sum()==12
    assert bin_index(0,Q_EDGES)==2 and bin_index(.1,B_EDGES)==2
    print("PASS: reflection pooling, population lattice, kernels, and physical bins")


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path);p.add_argument("--processed-dir",type=Path);p.add_argument("--cache-dir",type=Path);p.add_argument("--output",type=Path);p.add_argument("--self-test",action="store_true");args=p.parse_args()
    if args.self_test:self_test();return
    if not all((args.raw_dir,args.processed_dir,args.cache_dir,args.output)):p.error("data and output paths are required")
    args.output.mkdir(parents=True,exist_ok=True);series={s:load_q(args.cache_dir,s)[1] for s in range(40)}
    audit=native_resolution_audit(args.processed_dir,series);lp,lm,rate_rows=reflection_rates(series,list(range(40)));train_lp,train_lm,train_rates=reflection_rates(series,list(range(30)))
    kernels=kernel_tests(series);histories=[history_test(series,t) for t in LAGS]
    features=sample_features(args.raw_dir,args.cache_dir,series,args.cache_dir/"population_kinetics_features_dt0p5.npz")
    models,intensity_detail,nested=intensity_models(features);relax=fixed_q_relaxation(features,series);currents,current_summary=current_field(features)
    write_csv(args.output/"population_jump_rates.csv",rate_rows);write_csv(args.output/"kernel_tests.csv",kernels);write_csv(args.output/"history_tests.csv",histories);write_csv(args.output/"conditioned_intensity_models.csv",models);write_csv(args.output/"nested_intensity_comparisons.csv",nested);write_csv(args.output/"conditioned_intensities.csv",intensity_detail);write_csv(args.output/"fixed_q_relaxation.csv",relax);write_csv(args.output/"projected_probability_currents.csv",currents);write_csv(args.output/"native_resolution_audit.csv",audit.pop("rows"))
    make_figure(rate_rows,kernels,histories,relax,currents,args.output/"population_kinetic_decision.png")
    total_right=sum(np.diff(n)[np.diff(n)>0].sum() for n in series.values());total_left=sum(-np.diff(n)[np.diff(n)<0].sum() for n in series.values())
    summary=dict(configuration=dict(seeds=list(range(40)),train_seeds=list(range(30)),test_seeds=list(range(30,40)),q_dt=DT_Q,feature_dt=DT_FEATURE,transient=TRANSIENT,new_simulations=0,fitted_potential=0),
                 crossing_sample=dict(left_to_right=int(total_right),right_to_left=int(total_left),total=int(total_right+total_left),resolution_audit=audit),
                 kernel_tests=kernels,history_tests=histories,conditioned_intensity_models=models,nested_intensity_comparisons=nested,fixed_q_relaxation=relax,projected_current=current_summary,
                 potential_policy="derive a birth-death stationary potential only if CK and history tests support elimination")
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (args.output/"analysis_configuration.json").write_text(json.dumps(dict(raw_dir=str(args.raw_dir),processed_dir=str(args.processed_dir),cache_dir=str(args.cache_dir),common_population_resolution=DT_Q,reflection_symmetry="(q,px,b_left,b_right) -> (-q,-px,b_right,b_left)",new_simulations=0),indent=2)+"\n")


if __name__=="__main__":main()
