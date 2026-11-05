#!/usr/bin/env python3
"""Read-only validation of spatial harmonics in circular full-BEM trajectories.

The large MATLAB arrays are opened with read-only NumPy memory maps.  Only
bounded frame batches are copied.  No simulator function or trajectory file is
modified.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_collective import classify
from matv5_memmap import open_state_rdot
from process_collective import FILE_PATTERN, seed_from_name
from replay_torques import persistent_label


RADIUS = 7.158910531638177
DT = 0.001
N_FISH = 100
CHIRALITY_THRESHOLD = 0.2


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalise_state_batch(array: np.ndarray) -> np.ndarray:
    """Return a copied state batch as (component, fish, time)."""
    out = np.asarray(array, dtype=np.float64)
    if out.ndim == 2:
        return out[:, :, None]
    if out.shape[0] == 3:
        return out
    if out.shape[-1] == 3:
        return np.moveaxis(out, -1, 0)
    raise ValueError(f"Unexpected state batch shape {out.shape}")


def frame_metrics(state_batch: np.ndarray) -> dict[str, np.ndarray]:
    x, y, theta = state_batch
    phi = np.arctan2(y, x)
    alpha = np.arctan2(np.sin(theta-phi), np.cos(theta-phi))
    p = np.mean(np.exp(1j*theta), axis=0)
    M = np.mean(np.sin(alpha), axis=0)
    A1 = np.mean(np.exp(1j*phi), axis=0)
    A2 = np.mean(np.exp(2j*phi), axis=0)
    # The identity is exact for either s=+/-1.  sign(0) is assigned +1 only
    # for the algebraic audit; the tangent interpretation is flagged below.
    s = np.where(M < 0, -1.0, 1.0)
    H1 = np.mean(
        np.exp(1j*phi)*(np.exp(1j*alpha)-1j*s[None, :]), axis=0
    )
    residual = np.abs(p-(1j*s*A1+H1))

    comx, comy = np.mean(x, axis=0), np.mean(y, axis=0)
    xc, yc = x-comx, y-comy
    rc = np.hypot(xc, yc)
    valid_com = rc > 1e-8*RADIUS
    phi_com = np.arctan2(yc, xc)
    exp_com = np.where(valid_com, np.exp(1j*phi_com), 0.0)
    n_com = np.sum(valid_com, axis=0)
    A1_com = np.sum(exp_com, axis=0)/np.maximum(n_com, 1)
    M_com = np.sum(
        np.where(valid_com, np.sin(theta-phi_com), 0.0), axis=0
    )/np.maximum(n_com, 1)
    tangent = 1j*s*A1
    phase_cos = np.real(tangent*np.conj(H1))/np.maximum(np.abs(tangent)*np.abs(H1), 1e-15)
    return {
        "P": np.abs(p), "p": p, "M": M, "A1": A1, "A2": A2,
        "H1": H1, "identity_residual": residual,
        "com_radius": np.hypot(comx, comy)/RADIUS,
        "A1_com": A1_com, "M_com": M_com,
        "com_excluded": N_FISH-n_com,
        "vector_coherence": np.abs(p)/np.maximum(np.abs(A1)+np.abs(H1), 1e-15),
        "tangent_h1_phase_cos": phase_cos,
    }


def load_labels(outputs: Path, seed: int) -> tuple[np.lib.npyio.NpzFile, np.ndarray]:
    z = np.load(outputs/f"metrics_seed{seed}.npz")
    metric_dt = float(z["time"][1]-z["time"][0])
    raw, _, _ = classify(z["P"], z["Pi"], metric_dt)
    return z, persistent_label(raw, metric_dt)


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator, n_boot: int=20000) -> tuple[float, float, float]:
    values = np.asarray(values, float)
    draws = np.mean(values[rng.integers(0, len(values), (n_boot, len(values)))], axis=1)
    lo, hi = np.quantile(draws, [.025, .975])
    return float(np.mean(values)), float(lo), float(hi)


def bootstrap_difference(a: np.ndarray, b: np.ndarray, rng: np.random.Generator, n_boot: int=20000) -> tuple[float, float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    da = np.mean(a[rng.integers(0, len(a), (n_boot, len(a)))], axis=1)
    db = np.mean(b[rng.integers(0, len(b), (n_boot, len(b)))], axis=1)
    d = db-da
    lo, hi = np.quantile(d, [.025, .975])
    return float(np.mean(b)-np.mean(a)), float(lo), float(hi)


def orbital_sample(state, rdot, steps: np.ndarray, lags: tuple[int, ...],
                   radial_cutoff: float) -> list[dict]:
    rows: list[dict] = []
    nt = state.shape[2]
    for step in steps:
        if step+max(lags) >= nt-1:
            continue
        for lag in lags:
            positions=[np.asarray(state[:2, :, step+q], float) for q in range(lag+1)]
            valid=np.ones(positions[0].shape[1],dtype=bool)
            for pos in positions:
                valid &= np.hypot(pos[0],pos[1]) > radial_cutoff*RADIUS
            if not np.any(valid):
                continue
            rhs_values = []
            rigid_values = []
            diff_values = []
            excluded_total = 0
            for q in range(lag):
                st = np.asarray(state[:, :, step+q], float)
                vel = np.asarray(rdot[:, :, step+q+1], float)
                x, y = st[0], st[1]
                rr2 = x*x+y*y
                excluded_total += int(np.count_nonzero(~valid))
                phi = np.arctan2(y[valid], x[valid])
                ephi = np.exp(1j*phi)
                omega = (x[valid]*vel[1, valid]-y[valid]*vel[0, valid])/rr2[valid]
                A1 = np.mean(ephi)
                omega_bar = np.mean(omega)
                rigid = 1j*omega_bar*A1
                differential = 1j*np.mean((omega-omega_bar)*ephi)
                rhs_values.append(rigid+differential)
                rigid_values.append(rigid)
                diff_values.append(differential)
            first = np.asarray(state[:, :, step], float)
            last = np.asarray(state[:, :, step+lag], float)
            A_first = np.mean(np.exp(1j*np.arctan2(first[1,valid], first[0,valid])))
            A_last = np.mean(np.exp(1j*np.arctan2(last[1,valid], last[0,valid])))
            lhs = (A_last-A_first)/(lag*DT)
            rhs = np.mean(rhs_values)
            rigid = np.mean(rigid_values)
            differential = np.mean(diff_values)
            rows.append({
                "step": int(step), "time": step*DT, "lag_steps": lag,
                "lhs_real": lhs.real, "lhs_imag": lhs.imag,
                "rhs_real": rhs.real, "rhs_imag": rhs.imag,
                "rigid_real": rigid.real, "rigid_imag": rigid.imag,
                "differential_real": differential.real,
                "differential_imag": differential.imag,
                "abs_error": abs(lhs-rhs), "rhs_abs": abs(rhs),
                "excluded_fish_steps": excluded_total,
                "n_included_fish":int(np.count_nonzero(valid)),
                "radial_cutoff_over_R":radial_cutoff,
            })
    return rows


def summarise_orbital(rows: list[dict]) -> list[dict]:
    result=[]
    for state_name in ("schooling", "milling"):
        for lag in sorted({int(r["lag_steps"]) for r in rows}):
            use=[r for r in rows if r["macrostate"]==state_name and int(r["lag_steps"])==lag]
            if not use:
                continue
            lhs=np.array([complex(r["lhs_real"],r["lhs_imag"]) for r in use])
            rhs=np.array([complex(r["rhs_real"],r["rhs_imag"]) for r in use])
            rigid=np.array([complex(r["rigid_real"],r["rigid_imag"]) for r in use])
            diff=np.array([complex(r["differential_real"],r["differential_imag"]) for r in use])
            rmse=float(np.sqrt(np.mean(np.abs(lhs-rhs)**2)))
            result.append({
                "macrostate":state_name,"lag_steps":lag,"delta_t":lag*DT,
                "n_samples":len(use),"complex_rmse":rmse,
                "relative_rmse":rmse/max(float(np.sqrt(np.mean(np.abs(rhs)**2))),1e-15),
                "mean_abs_rhs":float(np.mean(np.abs(rhs))),
                "mean_abs_rigid":float(np.mean(np.abs(rigid))),
                "mean_abs_differential":float(np.mean(np.abs(diff))),
                "differential_fraction_of_sum":float(np.mean(np.abs(diff)/(np.abs(rigid)+np.abs(diff)+1e-15))),
                "excluded_fish_steps":int(sum(int(r["excluded_fish_steps"]) for r in use)),
            })
    return result


def sustained_change(time: np.ndarray, value: np.ndarray, baseline: np.ndarray, direction: str,
                     zscore: float=3.0, window: int=100, required: int=70) -> float | None:
    mean, sd=float(np.mean(baseline)),float(np.std(baseline,ddof=1))
    flag=value < mean-zscore*sd if direction=="down" else value > mean+zscore*sd
    kernel=np.ones(window,dtype=int)
    counts=np.convolve(flag.astype(int),kernel,mode="valid")
    found=np.flatnonzero(counts>=required)
    return None if not len(found) else float(time[int(found[0])])


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--data-dir",type=Path,required=True)
    p.add_argument("--outputs",type=Path,default=Path(__file__).resolve().parent/"outputs")
    p.add_argument("--outdir",type=Path,default=Path(__file__).resolve().parent/"outputs"/"spatial_harmonic_mechanism")
    p.add_argument("--sample-stride",type=int,default=100,help="raw steps; default 0.1 time units")
    p.add_argument("--orbital-per-seed",type=int,default=300)
    args=p.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    files=sorted(args.data_dir.glob(FILE_PATTERN),key=seed_from_name)
    if len(files)!=10:
        raise FileNotFoundError(f"Expected 10 files, found {len(files)} in {args.data_dir}")

    frame_arrays=[]; orbital_rows=[]; rng=np.random.default_rng(20261104)
    for file in files:
        seed=seed_from_name(file); z,label=load_labels(args.outputs,seed)
        state,rdot,_=open_state_rdot(file)
        steps=np.arange(20000,state.shape[2]-2,args.sample_stride,dtype=np.int64)
        collected={k:[] for k in ("P","M","A1","A2","H1","identity_residual","com_radius","A1_com","M_com","com_excluded","vector_coherence","tangent_h1_phase_cos")}
        for start in range(0,len(steps),1000):
            st=steps[start:start+1000]
            batch=normalise_state_batch(state[:,:,st])
            metrics=frame_metrics(batch)
            for key in collected:
                collected[key].append(metrics[key])
        data={k:np.concatenate(v) for k,v in collected.items()}
        metric_i=steps//10
        lab=label[metric_i]
        frame_arrays.append({"seed":np.full(len(steps),seed,np.int8),"step":steps,"time":steps*DT,
                             "label":lab,**data})

        eligible=np.flatnonzero((steps>=20000)&(steps<=999000)&(lab!=0))
        if len(eligible):
            pick=eligible[np.linspace(0,len(eligible)-1,min(args.orbital_per_seed,len(eligible))).round().astype(int)]
            rows=orbital_sample(state,rdot,steps[pick],(1,2,5,10),.01)
            for row in rows:
                mi=min(int(round(row["step"]/10)),len(label)-1)
                row["seed"]=seed
                row["macrostate"]={-1:"milling",0:"intermediate",1:"schooling"}[int(label[mi])]
            orbital_rows.extend(rows)
        print(f"seed {seed}: {len(steps)} harmonic frames",flush=True)

    merged={key:np.concatenate([d[key] for d in frame_arrays]) for key in frame_arrays[0]}
    # Compact sampled diagnostics.  Complex arrays are retained directly.
    np.savez_compressed(args.outdir/"harmonic_sampled_frames.npz",**{
        k:(v.astype(np.complex64) if np.iscomplexobj(v) else v.astype(np.float32)
           if np.issubdtype(v.dtype,np.floating) else v)
        for k,v in merged.items()
    })

    stationary_rows=[]
    for seed in range(10):
        for value,name in ((1,"schooling"),(-1,"milling")):
            use=(merged["seed"]==seed)&(merged["label"]==value)
            if not np.any(use):
                continue
            p=merged["P"][use]; a1=np.abs(merged["A1"][use]); a2=np.abs(merged["A2"][use]); h=np.abs(merged["H1"][use])
            elig=np.abs(merged["M"][use])>=CHIRALITY_THRESHOLD
            stationary_rows.append({
                "seed":seed,"macrostate":name,"n_frames":int(np.sum(use)),
                "eligible_chirality_fraction":float(np.mean(elig)),
                "mean_P":float(np.mean(p)),"mean_abs_M":float(np.mean(np.abs(merged["M"][use]))),
                "mean_abs_A1":float(np.mean(a1)),"mean_abs_A2":float(np.mean(a2)),
                "mean_abs_H1":float(np.mean(h)),"mean_H1_over_P":float(np.mean(h/np.maximum(p,1e-12))),
                "mean_tangent_fraction":float(np.mean(a1/np.maximum(a1+h,1e-12))),
                "mean_vector_coherence":float(np.mean(merged["vector_coherence"][use])),
                "mean_com_radius":float(np.mean(merged["com_radius"][use])),
                "mean_abs_A1_com":float(np.mean(np.abs(merged["A1_com"][use]))),
                "mean_abs_M_com":float(np.mean(np.abs(merged["M_com"][use]))),
                "corr_P_abs_A1":float(np.corrcoef(p,a1)[0,1]),
                "max_identity_residual":float(np.max(merged["identity_residual"][use])),
                "com_center_exclusions":int(np.sum(merged["com_excluded"][use])),
            })
    write_csv(args.outdir/"stationary_seed_summary.csv",stationary_rows)

    comparison_rows=[]
    quantities=("mean_P","mean_abs_M","mean_abs_A1","mean_abs_A2","mean_abs_H1","mean_H1_over_P","mean_com_radius","mean_abs_A1_com")
    for q in quantities:
        school=np.array([float(r[q]) for r in stationary_rows if r["macrostate"]=="schooling"])
        mill=np.array([float(r[q]) for r in stationary_rows if r["macrostate"]=="milling"])
        sm,slo,shi=bootstrap_mean(school,rng); mm,mlo,mhi=bootstrap_mean(mill,rng)
        diff,dlo,dhi=bootstrap_difference(school,mill,rng)
        comparison_rows.append({"quantity":q,"school_mean":sm,"school_ci_low":slo,"school_ci_high":shi,
            "mill_mean":mm,"mill_ci_low":mlo,"mill_ci_high":mhi,"mill_minus_school":diff,
            "difference_ci_low":dlo,"difference_ci_high":dhi,"n_school_units":len(school),"n_mill_units":len(mill)})
    write_csv(args.outdir/"stationary_comparison.csv",comparison_rows)

    identity={
        "n_sampled_frames":int(len(merged["P"])),
        "max_identity_residual":float(np.max(merged["identity_residual"])),
        "rms_identity_residual":float(np.sqrt(np.mean(merged["identity_residual"]**2))),
        "chirality_threshold":CHIRALITY_THRESHOLD,
        "ambiguous_fraction":float(np.mean(np.abs(merged["M"])<CHIRALITY_THRESHOLD)),
        "threshold_sensitivity":{},
    }
    for threshold in (.1,.2,.3):
        use=np.abs(merged["M"])>=threshold
        identity["threshold_sensitivity"][str(threshold)]={
            "eligible_fraction":float(np.mean(use)),
            "mean_H1_over_P":float(np.mean(np.abs(merged["H1"][use])/np.maximum(merged["P"][use],1e-12))),
            "corr_P_abs_A1":float(np.corrcoef(merged["P"][use],np.abs(merged["A1"][use]))[0,1]),
        }
    (args.outdir/"identity_audit.json").write_text(json.dumps(identity,indent=2)+"\n")

    relationship_rows=[]
    for value,name in ((1,"schooling"),(-1,"milling")):
        use=merged["label"]==value
        p=merged["P"][use]; a1=np.abs(merged["A1"][use]); h=np.abs(merged["H1"][use])
        relationship_rows.append({
            "macrostate":name,"n_frames":int(np.sum(use)),
            "corr_P_abs_A1":float(np.corrcoef(p,a1)[0,1]),
            "rms_tangent_only_error":float(np.sqrt(np.mean(h*h))),
            "rms_error_over_rms_P":float(np.sqrt(np.mean(h*h))/np.sqrt(np.mean(p*p))),
            "mean_H1_over_P":float(np.mean(h/np.maximum(p,1e-12))),
            "mean_tangent_H1_phase_cos":float(np.mean(merged["tangent_h1_phase_cos"][use])),
            "mean_vector_coherence":float(np.mean(merged["vector_coherence"][use])),
            "ambiguous_chirality_fraction":float(np.mean(np.abs(merged["M"][use])<CHIRALITY_THRESHOLD)),
        })
    write_csv(args.outdir/"state_relationship_summary.csv",relationship_rows)

    # Seed 9 high-resolution transition window.
    file=files[9]; state,_,_=open_state_rdot(file); z,label=load_labels(args.outputs,9)
    steps=np.arange(280000,340001,10,dtype=np.int64)
    hm=frame_metrics(normalise_state_batch(state[:,:,steps]))
    transition_rows=[]
    for j,step in enumerate(steps):
        transition_rows.append({"step":int(step),"time":step*DT,"persistent_label":int(label[step//10]),
            "P":hm["P"][j],"M":hm["M"][j],"abs_A1":abs(hm["A1"][j]),"abs_A2":abs(hm["A2"][j]),
            "abs_H1":abs(hm["H1"][j]),"com_radius":hm["com_radius"][j],
            "abs_A1_com":abs(hm["A1_com"][j]),"M_com":hm["M_com"][j]})
    write_csv(args.outdir/"seed9_transition_timeline.csv",transition_rows)
    tt=steps*DT; baseline=(tt>=280)&(tt<300)
    chronology={"baseline":"seed 9, t in [280,300)","transition_reference_time":312.41,"rule":"first 1 s window with at least 70 of 100 samples beyond 3 baseline SD"}
    for name,direction in (("P","down"),("abs_M","down"),("abs_A1","down"),("abs_A2","up"),("abs_H1","up"),("com_radius","down")):
        values=np.abs(hm["M"]) if name=="abs_M" else np.abs(hm[name[4:]]) if name.startswith("abs_") else hm[name]
        chronology[name]=sustained_change(tt,values,values[baseline],direction)
    (args.outdir/"seed9_transition_chronology.json").write_text(json.dumps(chronology,indent=2)+"\n")

    write_csv(args.outdir/"orbital_identity_samples.csv",orbital_rows)
    orbital_summary=summarise_orbital(orbital_rows)
    write_csv(args.outdir/"orbital_identity_summary.csv",orbital_summary)

    # Decisive figures.
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
    axs[0].plot(tt,hm["P"],label="P",lw=1); axs[0].plot(tt,np.abs(hm["M"]),label="|M|",lw=1)
    axs[0].plot(tt,np.abs(hm["A1"]),label="|A1|",lw=1); axs[0].plot(tt,np.abs(hm["H1"]),label="|H1|",lw=1)
    axs[1].plot(tt,np.abs(hm["A2"]),label="|A2|",lw=1,color="#f58518")
    axs[1].plot(tt,hm["com_radius"],label="COM radius / R",lw=1,color="#54a24b")
    axs[1].plot(tt,np.abs(hm["A1_com"]),label="|A1| about COM",lw=1,color="#b279a2")
    for ax in axs: ax.axvline(312.41,color="black",ls="--",lw=1); ax.grid(alpha=.2); ax.legend(ncol=4,fontsize=8)
    axs[1].set_xlabel("time"); fig.suptitle("Seed 9: spatial harmonics through the schooling-to-milling event")
    fig.tight_layout(); fig.savefig(args.outdir/"seed9_spatial_harmonic_timeline.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.2,5.2))
    keep=(merged["label"]!=0); select=np.flatnonzero(keep)[::10]
    colors=np.where(merged["label"][select]==1,"#4c78a8","#e45756")
    ax.scatter(np.abs(merged["A1"][select]),merged["P"][select],c=colors,s=5,alpha=.2,rasterized=True)
    lim=max(float(np.max(np.abs(merged["A1"][select]))),float(np.max(merged["P"][select])))
    ax.plot([0,lim],[0,lim],"k--",lw=1,label="tangent-only P=|A1|")
    ax.set(xlabel="|A1| about fixed circle center",ylabel="polarization P",title="H1 prevents a tangent-only identification")
    ax.grid(alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(args.outdir/"P_vs_A1_stationary.png",dpi=220); plt.close(fig)

    fig,axs=plt.subplots(1,3,figsize=(12,4))
    for ax,key,title in zip(axs,("mean_abs_A1","mean_abs_A2","mean_abs_H1"),("|A1|","|A2|","|H1|")):
        for x,name,color in ((0,"schooling","#4c78a8"),(1,"milling","#e45756")):
            vals=[float(r[key]) for r in stationary_rows if r["macrostate"]==name]
            ax.scatter(np.full(len(vals),x)+np.linspace(-.05,.05,len(vals)),vals,color=color,s=35)
            ax.plot([x-.13,x+.13],[np.mean(vals)]*2,color="black",lw=2)
        ax.set_xticks([0,1],["school","mill"]); ax.set_title(title); ax.grid(alpha=.2)
    fig.suptitle("Independent seed/interval means"); fig.tight_layout(); fig.savefig(args.outdir/"stationary_harmonic_comparison.png",dpi=220); plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.5,4.5))
    for state_name,color in (("schooling","#4c78a8"),("milling","#e45756")):
        use=[r for r in orbital_summary if r["macrostate"]==state_name]
        ax.plot([r["delta_t"] for r in use],[r["relative_rmse"] for r in use],"o-",label=state_name,color=color)
    ax.set(xlabel="finite-difference interval",ylabel="complex relative RMSE",title="Orbital A1 identity discretization check")
    ax.grid(alpha=.2); ax.legend(); fig.tight_layout(); fig.savefig(args.outdir/"orbital_identity_discretization.png",dpi=220); plt.close(fig)

    config={"data_pattern":FILE_PATTERN,"radius":RADIUS,"dt":DT,"n_fish":N_FISH,
        "uniform_sample_stride_steps":args.sample_stride,"uniform_sample_dt":args.sample_stride*DT,
        "formation_transient":20.0,"seed9_window":[280.0,340.0],"seed9_dt":.01,
        "chirality_threshold":CHIRALITY_THRESHOLD,"chirality_sensitivity":[.1,.2,.3],
        "com_radial_exclusion":"r_com <= 1e-8 R","orbital_radial_exclusion":"r <= 0.01 R at any point in differencing interval",
        "orbital_velocity_index":"rdot[:,:,step+1] advances state[:,:,step]",
        "orbital_lag_steps":[1,2,5,10],"orbital_samples_per_seed":args.orbital_per_seed,
        "bootstrap_units":"independent seed/state intervals","bootstrap_replicates":20000,"random_seed":20261104}
    (args.outdir/"analysis_config.json").write_text(json.dumps(config,indent=2)+"\n")
    print(f"saved {args.outdir}",flush=True)


if __name__=="__main__":
    main()
