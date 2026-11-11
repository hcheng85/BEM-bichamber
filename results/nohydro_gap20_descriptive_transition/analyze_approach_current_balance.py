#!/usr/bin/env python3
"""Source-exact balance, plug-flow audit, and closure test for approach current.

The simulator and raw MAT files are immutable inputs.  Gap-20 geometry is
constructed in memory from the archived production source by changing only its
top-level geometry parameter from gap_ratio=5 to gap_ratio=20.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import os
import sys
import types
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "theory" / "circle_collective_first_stage"))
sys.path.insert(0, str(HERE.parent))
from matv5_memmap import inspect_numeric_arrays  # noqa: E402


N = 100
DT = 0.001
IN = 0.5
ELL = 2.8635642126552705
WIDTH = ELL


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def numeric_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Recover numeric CSV fields while retaining categorical strings."""
    out = []
    for row in rows:
        converted: dict[str, Any] = {}
        for key, value in row.items():
            try:
                number = float(value)
                converted[key] = int(number) if number.is_integer() else number
            except (TypeError, ValueError):
                converted[key] = value
        out.append(converted)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


class RawReader:
    def __init__(self, path: Path):
        self.path = path; self.info = inspect_numeric_arrays(path); self.fd = os.open(path, os.O_RDONLY)

    def block(self, name: str, first: int, last: int) -> np.ndarray:
        info = self.info[name]; shape = info.shape[:-1]
        frame_bytes = int(np.prod(shape)) * info.dtype.itemsize
        count = last - first + 1
        raw = os.pread(self.fd, count * frame_bytes, info.data_offset + first * frame_bytes)
        return np.frombuffer(raw, dtype=info.dtype).reshape((*shape, count), order="F")

    def close(self) -> None:
        os.close(self.fd)


def load_gap20_production(source: Path):
    """Execute the unchanged source logic with only the archived gap parameter set to 20."""
    text = source.read_text()
    old = "parameters_list = [50, 5, 20]"
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one exact geometry assignment {old!r}")
    text = text.replace(old, "parameters_list = [50, 20, 20]", 1)
    if "tqdm" not in sys.modules:
        shim = types.ModuleType("tqdm"); shim.tqdm = lambda iterable, *a, **k: iterable; sys.modules["tqdm"] = shim
    module = types.ModuleType("gap20_nohydro_production_replay")
    module.__file__ = str(source)
    cwd = os.getcwd(); os.chdir(source.parent)
    try:
        exec(compile(text, str(source) + "::gap20_geometry", "exec"), module.__dict__)
    finally:
        os.chdir(cwd)
    if not np.isclose(module.gap, ELL, rtol=0, atol=1e-12):
        raise RuntimeError((module.gap, ELL))
    return module


def fast_check_voronoi(mod, vn: np.ndarray, location: np.ndarray) -> np.ndarray:
    """Loop only over Delaunay edges; predicates are exactly the production predicates."""
    out = vn.copy()
    for i, j in np.argwhere(vn):
        a = location[:, i]; b = location[:, j]
        if mod.check_intersections(mod.segments, [a, b]):
            out[i, j] = False; continue
        if mod.check_inersect_circle(a, b, mod.circle_center6, mod.r):
            out[i, j] = False; continue
        if mod.check_inersect_circle(a, b, mod.circle_center5, mod.r):
            out[i, j] = False
    return out


def replay_torques(mod, state: np.ndarray, para: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vn0 = mod.GetVoronoi(state[:2], para)
    vn = fast_check_voronoi(mod, vn0, state[:2])
    # GetVision prints once per step if any fish has zero visual weight; its
    # documented fallback is retained while only the diagnostic print is muted.
    with np.errstate(invalid="ignore", divide="ignore"), contextlib.redirect_stdout(io.StringIO()):
        social = np.asarray(mod.GetVision(vn, state, para)).reshape(-1)
    wall = np.asarray(mod.GetWallAvoid(state, para)).reshape(-1)
    raw = wall + social
    effective = np.clip(raw, -np.pi / para["dt"], np.pi / para["dt"])
    return wall, social, effective - raw, effective


def indicator(state: np.ndarray, direction: int) -> np.ndarray:
    x, y = state[0], state[1]; s = float(direction)
    return (s * x < 0) & ((-s * x) <= ELL) & (np.abs(y) <= WIDTH / 2)


def current(state: np.ndarray, velocity: np.ndarray, direction: int) -> tuple[float, int]:
    inside = indicator(state, direction)
    return float(np.sum(direction * velocity[0] * inside) / N), int(np.sum(inside))


def exact_balance(mod, reader: RawReader, start: int, stop: int, direction: int, para: dict[str, Any]) -> tuple[list[dict], dict]:
    """Native-step exact finite-step balance on [start, stop]."""
    state = reader.block("state", start, stop)
    velocity = reader.block("rdot", start + 1, stop + 1)
    components = {name: np.zeros(stop - start) for name in
                  ("boundary", "wall", "social", "clip", "noise_relaxation", "realized_noise", "collision_override")}
    measured = np.zeros(stop - start)
    noise_z = []
    exp_noise = math.exp(-0.5 * IN * IN * DT)
    maximum_visibility_error = 0
    clip_count = 0
    for k in range(stop - start):
        st = state[:, :, k]; nxt = state[:, :, k + 1]
        vel = velocity[:, :, k]; vel_next = velocity[:, :, k + 1]
        inside = indicator(st, direction); inside_next = indicator(nxt, direction)
        theta = st[2]; theta_next = nxt[2]
        wall, social, clip, omega = replay_torques(mod, st, para)
        clip_count += int(np.count_nonzero(clip))
        delta = omega * DT
        small = np.abs(delta) < 1e-12
        slope = np.empty_like(delta)
        slope[small] = -np.sin(theta[small])
        slope[~small] = (np.cos(theta[~small] + delta[~small]) - np.cos(theta[~small])) / delta[~small]
        factor = direction * inside / N
        components["wall"][k] = np.sum(factor * slope * wall * DT)
        components["social"][k] = np.sum(factor * slope * social * DT)
        components["clip"][k] = np.sum(factor * slope * clip * DT)
        cos_det = np.cos(theta + delta)
        components["noise_relaxation"][k] = np.sum(factor * (exp_noise - 1) * cos_det)
        components["realized_noise"][k] = np.sum(factor * (np.cos(theta_next) - exp_noise * cos_det))
        components["boundary"][k] = np.sum(direction * (inside_next.astype(int) - inside.astype(int)) * np.cos(theta_next) / N)
        ideal = np.sum(direction * inside * np.cos(theta)) / N
        ideal_next = np.sum(direction * inside_next * np.cos(theta_next)) / N
        actual = np.sum(direction * inside * vel[0]) / N
        actual_next = np.sum(direction * inside_next * vel_next[0]) / N
        components["collision_override"][k] = (actual_next - ideal_next) - (actual - ideal)
        measured[k] = actual_next - actual
        noise_z.extend(((theta_next - theta - delta) / (IN * math.sqrt(DT))).tolist())
    reconstructed = sum(components.values())
    closure = measured - reconstructed
    half = (stop - start) // 2
    rows = []
    for label, lo, hi in (("first_half", 0, half), ("second_half", half, stop - start), ("full", 0, stop - start)):
        row = {"interval": label, "measured_delta_Ja": float(np.sum(measured[lo:hi]))}
        for name, values in components.items(): row[name] = float(np.sum(values[lo:hi]))
        row["reconstructed_delta_Ja"] = float(np.sum(reconstructed[lo:hi]))
        row["closure_residual"] = float(np.sum(closure[lo:hi]))
        rows.append(row)
    audit = dict(max_abs_step_closure=float(np.max(np.abs(closure))), full_closure=float(np.sum(closure)),
                 noise_mean=float(np.mean(noise_z)), noise_variance=float(np.var(noise_z, ddof=1)),
                 angular_clip_count=clip_count, n_native_steps=stop-start,
                 unused_visibility_error=maximum_visibility_error)
    return rows, audit


def choose_pairs(events: list[dict[str, str]], matched: list[dict[str, str]]) -> list[dict[str, Any]]:
    event_map = {(int(r["seed"]), int(r["event_id"])): r for r in events}
    lead2 = [r for r in matched if float(r["lead"]) == 2.0]
    candidates = []
    for r in lead2:
        key = (int(r["event_seed"]), int(r["event_id"]))
        ev = event_map[key]
        candidates.append(dict(event_seed=key[0], event_id=key[1], event_time=float(ev["time"]),
                               first_native=int(ev["first_native"]), direction=int(r["direction"]), q=float(r["q"]),
                               control_seed=int(r["control_seed"]), control_time=float(r["control_time"])))
    chosen = []
    for seed in sorted({r["event_seed"] for r in candidates}):
        options = sorted((r for r in candidates if r["event_seed"] == seed), key=lambda r: (r["event_time"], r["event_id"]))
        chosen.append(options[0])
    return chosen


def raw_path(raw_dir: Path, seed: int) -> Path:
    return raw_dir / f"If00_Ip40_In05_R50_num100_wid20_n{seed}_BEM_no_hydro.mat"


def coarse_series(reader: RawReader, start: int, stop: int, direction: int, stride: int = 10) -> dict[str, np.ndarray]:
    states = reader.block("state", start, stop)
    velocities = reader.block("rdot", start + 1, stop + 1)
    indices = np.arange(0, stop-start+1, stride)
    q=[]; ja=[]; na=[]
    for k in indices:
        st=states[:,:,k]; vel=velocities[:,:,k]
        q.append(1-2*np.sum(st[0]<0)/N)
        value,count=current(st,vel,direction); ja.append(value); na.append(count)
    return dict(step=start+indices,time=(start+indices)*DT,q=np.asarray(q),Ja=np.asarray(ja),Na=np.asarray(na),
                state=states[:,:,indices],velocity=velocities[:,:,indices])


def cluster_bootstrap_difference(rows: list[dict], key: str, nboot: int = 5000) -> tuple[float,float,float]:
    seeds=sorted({int(r["seed"]) for r in rows}); means=[]
    for seed in seeds: means.append(np.mean([float(r[key]) for r in rows if int(r["seed"])==seed]))
    values=np.asarray(means); rng=np.random.default_rng(20261111)
    draws=rng.choice(values,(nboot,len(values)),replace=True).mean(axis=1)
    return float(values.mean()),*map(float,np.quantile(draws,[.025,.975]))


def seed_bootstrap_mean(rows: list[dict], key: str, nboot: int = 5000) -> tuple[float,float,float]:
    values=np.asarray([float(r[key]) for r in rows if np.isfinite(float(r[key]))])
    if not len(values): return (float("nan"),)*3
    rng=np.random.default_rng(20261111);draws=rng.choice(values,(nboot,len(values)),replace=True).mean(axis=1)
    return float(values.mean()),*map(float,np.quantile(draws,[.025,.975]))


def balance_analysis(mod, para, raw_dir: Path, pairs: list[dict], out: Path, n_exact: int) -> tuple[list[dict], dict]:
    positions=np.linspace(0,len(pairs)-1,min(n_exact,len(pairs))).round().astype(int)
    selected=[pairs[i] for i in positions]
    rows=[]; audits=[]; visibility_checks=[]
    for pair in selected:
        print(f"native balance event seed {pair['event_seed']} control {pair['control_seed']}",flush=True)
        for kind,seed,start in (("event",pair["event_seed"],pair["first_native"]-2000),
                                ("control",pair["control_seed"],round(pair["control_time"]/DT))):
            reader=RawReader(raw_path(raw_dir,seed))
            # Verify optimized visibility against production at three fixed frames.
            for offset in (0,1000,1999):
                st=reader.block("state",start+offset,start+offset)[:,:,0]
                v0=mod.GetVoronoi(st[:2],para); a=mod.CheckVoronoi(v0.copy(),st[:2],para); b=fast_check_voronoi(mod,v0,st[:2])
                visibility_checks.append(int(np.count_nonzero(a!=b)))
            part,audit=exact_balance(mod,reader,start,start+2000,pair["direction"],para); reader.close()
            for r in part: rows.append({"seed":pair["event_seed"],"trajectory_seed":seed,"event_id":pair["event_id"],
                                        "kind":kind,"direction":pair["direction"],"q":pair["q"],**r})
            audits.append({"seed":pair["event_seed"],"kind":kind,**audit})
    write_csv(out/"native_balance_intervals.csv",rows); write_csv(out/"native_balance_audit.csv",audits)
    full=[r for r in rows if r["interval"]=="full"]
    comparison=[]
    for pair_seed in sorted({int(r["seed"]) for r in full}):
        ev=next(r for r in full if int(r["seed"])==pair_seed and r["kind"]=="event")
        ct=next(r for r in full if int(r["seed"])==pair_seed and r["kind"]=="control")
        row={"seed":pair_seed}
        for key in ("measured_delta_Ja","boundary","wall","social","clip","noise_relaxation","realized_noise","collision_override"):
            row[f"event_{key}"]=float(ev[key]);row[f"control_{key}"]=float(ct[key]);row[f"delta_{key}"]=float(ev[key])-float(ct[key])
        comparison.append(row)
    write_csv(out/"native_balance_event_control.csv",comparison)
    summary={}
    for key in ("measured_delta_Ja","boundary","wall","social","clip","noise_relaxation","realized_noise","collision_override"):
        summary[key]=cluster_bootstrap_difference(comparison,f"delta_{key}")
    audit=dict(n_pairs=len(comparison),pair_event_seeds=[r["seed"] for r in comparison],
               max_visibility_graph_disagreement=max(visibility_checks,default=0),
               max_abs_step_closure=max(float(r["max_abs_step_closure"]) for r in audits),
               max_abs_full_closure=max(abs(float(r["full_closure"])) for r in audits),
               pooled_noise_mean=float(np.mean([float(r["noise_mean"]) for r in audits])),
               pooled_noise_variance=float(np.mean([float(r["noise_variance"]) for r in audits])),
               paired_difference_seed_bootstrap={k:{"estimate":v[0],"ci95_low":v[1],"ci95_high":v[2]} for k,v in summary.items()})
    (out/"native_balance_summary.json").write_text(json.dumps(audit,indent=2)+"\n")
    return comparison,audit


def closure_analysis(raw_dir: Path, pairs: list[dict], out: Path) -> tuple[list[dict],dict]:
    rows=[]
    for pair in pairs:
        for kind,seed,start in (("event",pair["event_seed"],pair["first_native"]-2000),
                                ("control",pair["control_seed"],round(pair["control_time"]/DT))):
            reader=RawReader(raw_path(raw_dir,seed)); z=coarse_series(reader,start,start+2000,pair["direction"],stride=100);reader.close()
            for k in range(2,len(z["Ja"])-2):
                rows.append(dict(seed=pair["event_seed"],kind=kind,q=float(z["q"][k]),Ja=float(z["Ja"][k]),
                                 history_dJa=float((z["Ja"][k]-z["Ja"][k-2])/.2),
                                 future_dJa=float((z["Ja"][k+2]-z["Ja"][k])/.2)))
    # Fixed quantile bins are descriptive, not fitted dynamics.
    qedges=np.unique(np.quantile([r["q"] for r in rows],[0,.25,.5,.75,1])); jedges=np.unique(np.quantile([r["Ja"] for r in rows],[0,.2,.4,.6,.8,1]))
    for r in rows:
        r["q_bin"]=int(np.clip(np.searchsorted(qedges,float(r["q"]),side="right")-1,0,len(qedges)-2))
        r["Ja_bin"]=int(np.clip(np.searchsorted(jedges,float(r["Ja"]),side="right")-1,0,len(jedges)-2))
        r["history_sign"]="rising" if float(r["history_dJa"])>=0 else "falling"
    cell_rows=[]
    for qb in range(len(qedges)-1):
        for jb in range(len(jedges)-1):
            cell=[r for r in rows if r["q_bin"]==qb and r["Ja_bin"]==jb]
            for kind in ("event","control"):
                part=[r for r in cell if r["kind"]==kind]
                if part: cell_rows.append(dict(q_bin=qb,Ja_bin=jb,group=kind,n=len(part),n_seeds=len({r['seed'] for r in part}),mean_future_dJa=float(np.mean([r['future_dJa'] for r in part]))))
            for sign in ("rising","falling"):
                part=[r for r in cell if r["history_sign"]==sign]
                if part: cell_rows.append(dict(q_bin=qb,Ja_bin=jb,group=sign,n=len(part),n_seeds=len({r['seed'] for r in part}),mean_future_dJa=float(np.mean([r['future_dJa'] for r in part]))))
    write_csv(out/"closure_samples.csv",rows);write_csv(out/"closure_binned_drift.csv",cell_rows)
    # History effect: seed-level, within-cell rising minus falling, cells requiring both histories.
    seed_effect=[]
    for seed in sorted({r["seed"] for r in rows}):
        effects=[];weights=[]
        for qb in range(len(qedges)-1):
            for jb in range(len(jedges)-1):
                up=[r["future_dJa"] for r in rows if r["seed"]==seed and r["q_bin"]==qb and r["Ja_bin"]==jb and r["history_sign"]=="rising"]
                dn=[r["future_dJa"] for r in rows if r["seed"]==seed and r["q_bin"]==qb and r["Ja_bin"]==jb and r["history_sign"]=="falling"]
                if up and dn: effects.append(np.mean(up)-np.mean(dn));weights.append(min(len(up),len(dn)))
        if effects: seed_effect.append(dict(seed=seed,history_effect=float(np.average(effects,weights=weights)),n_cells=len(effects)))
    write_csv(out/"closure_seed_history_effect.csv",seed_effect)
    estimate,lo,hi=cluster_bootstrap_difference(seed_effect,"history_effect")
    # Reproducibility across lower/upper seed halves on identical bins.
    cell_half={}
    for half,seeds in (("low",set(range(20))),("high",set(range(20,40)))):
        for qb in range(len(qedges)-1):
            for jb in range(len(jedges)-1):
                vals=[r["future_dJa"] for r in rows if r["seed"] in seeds and r["q_bin"]==qb and r["Ja_bin"]==jb]
                if len(vals)>=10: cell_half[(half,qb,jb)]=float(np.mean(vals))
    common=[(qb,jb) for qb in range(len(qedges)-1) for jb in range(len(jedges)-1) if ("low",qb,jb) in cell_half and ("high",qb,jb) in cell_half]
    low=np.array([cell_half[("low",*c)] for c in common]);high=np.array([cell_half[("high",*c)] for c in common])
    summary=dict(q_edges=qedges.tolist(),Ja_edges=jedges.tolist(),n_samples=len(rows),n_seed_pairs=len(pairs),
                 history_effect=estimate,history_effect_ci95=[lo,hi],history_effect_n_seeds=len(seed_effect),
                 split_half_common_cells=len(common),split_half_drift_correlation=float(np.corrcoef(low,high)[0,1]) if len(common)>1 else None,
                 split_half_drift_rmse=float(np.sqrt(np.mean((low-high)**2))) if len(common) else None,
                 pooled_future_drift_sd=float(np.std([r["future_dJa"] for r in rows],ddof=1)))
    return rows,augment_closure_summary(rows,summary,out)


def augment_closure_summary(rows: list[dict], summary: dict, out: Path) -> dict:
    """Test hidden outcome/history dependence without fitting a drift model."""
    seed_kind=[]
    for seed in sorted({int(r["seed"]) for r in rows}):
        effects=[];weights=[]
        cells=sorted({(int(r["q_bin"]),int(r["Ja_bin"])) for r in rows if int(r["seed"])==seed})
        for qb,jb in cells:
            ev=[float(r["future_dJa"]) for r in rows if int(r["seed"])==seed and int(r["q_bin"])==qb and int(r["Ja_bin"])==jb and r["kind"]=="event"]
            ct=[float(r["future_dJa"]) for r in rows if int(r["seed"])==seed and int(r["q_bin"])==qb and int(r["Ja_bin"])==jb and r["kind"]=="control"]
            if ev and ct: effects.append(np.mean(ev)-np.mean(ct));weights.append(min(len(ev),len(ct)))
        if effects: seed_kind.append(dict(seed=seed,event_control_effect=float(np.average(effects,weights=weights)),n_cells=len(effects)))
    write_csv(out/"closure_seed_event_control_effect.csv",seed_kind)
    kind_est,kind_lo,kind_hi=seed_bootstrap_mean(seed_kind,"event_control_effect")

    # Remove the descriptive cell means, then ask whether the previous and next
    # 0.2-time-unit changes remain correlated within independent seed pairs.
    cell_means={}
    for cell in {(int(r["q_bin"]),int(r["Ja_bin"])) for r in rows}:
        part=[r for r in rows if (int(r["q_bin"]),int(r["Ja_bin"]))==cell]
        cell_means[cell]=(np.mean([float(r["history_dJa"]) for r in part]),np.mean([float(r["future_dJa"]) for r in part]))
    seed_memory=[]
    for seed in sorted({int(r["seed"]) for r in rows}):
        part=[r for r in rows if int(r["seed"])==seed]
        h=np.asarray([float(r["history_dJa"])-cell_means[(int(r["q_bin"]),int(r["Ja_bin"]))][0] for r in part])
        f=np.asarray([float(r["future_dJa"])-cell_means[(int(r["q_bin"]),int(r["Ja_bin"]))][1] for r in part])
        if len(part)>=8 and np.std(h)>0 and np.std(f)>0:
            seed_memory.append(dict(seed=seed,n=len(part),residual_spearman=float(spearmanr(h,f).statistic)))
    write_csv(out/"closure_seed_residual_memory.csv",seed_memory)
    mem_est,mem_lo,mem_hi=seed_bootstrap_mean(seed_memory,"residual_spearman")
    kind_memory=[];kind_stats={}
    for kind in ("event","control"):
        subset=[r for r in rows if r["kind"]==kind]
        kind_cell={}
        for cell in {(int(r["q_bin"]),int(r["Ja_bin"])) for r in subset}:
            part=[r for r in subset if (int(r["q_bin"]),int(r["Ja_bin"]))==cell]
            kind_cell[cell]=(np.mean([float(r["history_dJa"]) for r in part]),np.mean([float(r["future_dJa"]) for r in part]))
        current=[]
        for seed in sorted({int(r["seed"]) for r in subset}):
            part=[r for r in subset if int(r["seed"])==seed]
            h=np.asarray([float(r["history_dJa"])-kind_cell[(int(r["q_bin"]),int(r["Ja_bin"]))][0] for r in part])
            f=np.asarray([float(r["future_dJa"])-kind_cell[(int(r["q_bin"]),int(r["Ja_bin"]))][1] for r in part])
            if len(part)>=8 and np.std(h)>0 and np.std(f)>0:
                row=dict(seed=seed,kind=kind,n=len(part),residual_spearman=float(spearmanr(h,f).statistic));kind_memory.append(row);current.append(row)
        est,lo,hi=seed_bootstrap_mean(current,"residual_spearman")
        kind_stats[kind]=dict(estimate=est,ci95=[lo,hi],n_seeds=len(current))
    write_csv(out/"closure_seed_kind_specific_memory.csv",kind_memory)
    summary.update(event_control_effect=kind_est,event_control_effect_ci95=[kind_lo,kind_hi],
                   event_control_effect_n_seeds=len(seed_kind),residual_history_spearman=mem_est,
                   residual_history_spearman_ci95=[mem_lo,mem_hi],residual_history_n_seeds=len(seed_memory),
                   kind_specific_residual_history=kind_stats)
    (out/"closure_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary


def boundary_flux_window(states: np.ndarray, direction: int) -> tuple[float,float,float,float]:
    """Net rates across approach boundaries; exactly closes coarse membership."""
    s=float(direction); dt=(states.shape[2]-1)*DT
    x,y=states[0],states[1]
    inside=(s*x<0)&((-s*x)<=ELL)&(np.abs(y)<=WIDTH/2)
    changed=inside[:,:-1]!=inside[:,1:]
    signs=inside[:,:-1].astype(float)-inside[:,1:].astype(float)  # positive outward
    neck_cross=(((s*x[:,:-1])<0)&((s*x[:,1:])>=0))|(((s*x[:,:-1])>=0)&((s*x[:,1:])<0))
    d0=-s*x[:,:-1];d1=-s*x[:,1:]
    upstream_cross=((d0<ELL)&(d1>=ELL))|((d0>=ELL)&(d1<ELL))
    neck=float(np.sum(signs*changed*neck_cross))
    upstream=float(np.sum(signs*changed*(~neck_cross)*upstream_cross))
    lateral=float(np.sum(signs*changed*(~neck_cross)*(~upstream_cross)))
    accumulation=float(np.sum(inside[:,-1])-np.sum(inside[:,0]))/dt
    # dN/dt = -(out_neck + out_upstream + out_lateral)
    return neck/dt,upstream/dt,lateral/dt,accumulation


def plug_flow_analysis(raw_dir: Path, events: list[dict[str,str]], out: Path) -> tuple[list[dict],dict]:
    selected=[r for r in events if int(r["seed"])>=30 and int(r["continuation"])==1 and float(r["duration"])>=1.0]
    rows=[]
    for event_index,ev in enumerate(selected):
        if event_index % 20 == 0: print(f"plug audit {event_index}/{len(selected)}",flush=True)
        seed=int(ev["seed"]);direction=int(ev["direction"]);first=int(ev["first_native"])
        reader=RawReader(raw_path(raw_dir,seed)); states=reader.block("state",first,first+5000);vel=reader.block("rdot",first+1,first+5001);reader.close()
        for offset in range(0,3001,100):
            st=states[:,:,offset];vv=vel[:,:,offset];inside=indicator(st,direction);vdir=direction*vv[0]
            ja=float(np.sum(vdir*inside)/N);na=int(inside.sum())
            positive=inside&(vdir>0)
            if ja<=0 or not np.any(positive): continue
            distance=-direction*st[0]
            travel=float(np.sum(distance[positive])/np.sum(vdir[positive]))
            lag_steps=int(round(travel/DT)); window=200
            begin=offset+lag_steps;end=begin+window
            if end>=states.shape[2]: continue
            segment=states[:,:,begin:end+1]
            neck,upstream,lateral,accum=boundary_flux_window(segment,direction)
            neck_now,_,_,_=boundary_flux_window(states[:,:,offset:offset+window+1],direction)
            predicted=N*ja/ELL
            rows.append(dict(seed=seed,event_id=int(ev["event_id"]),sample_time=offset*DT,Ja=ja,Na=na,
                             plug_flux=predicted,travel_time=travel,actual_neck_flux=neck,actual_neck_flux_unshifted=neck_now,
                             upstream_outflux=upstream,lateral_outflux=lateral,approach_accumulation=accum,
                             conservation_residual=accum+neck+upstream+lateral,
                             plug_error=predicted-neck))
    write_csv(out/"plug_flow_heldout_samples.csv",rows)
    seed_rows=[]
    for seed in sorted({r["seed"] for r in rows}):
        part=[r for r in rows if r["seed"]==seed];p=np.array([r["plug_flux"] for r in part]);a=np.array([r["actual_neck_flux"] for r in part]);a0=np.array([r["actual_neck_flux_unshifted"] for r in part])
        seed_rows.append(dict(seed=seed,n=len(part),bias=float(np.mean(p-a)),mae=float(np.mean(np.abs(p-a))),
                              correlation=float(spearmanr(p,a).statistic) if len(np.unique(a))>1 else np.nan,
                              unshifted_bias=float(np.mean(p-a0)),unshifted_mae=float(np.mean(np.abs(p-a0))),
                              unshifted_correlation=float(spearmanr(p,a0).statistic) if len(np.unique(a0))>1 else np.nan,
                              conservation_rms=float(np.sqrt(np.mean(np.array([r['conservation_residual'] for r in part])**2)))))
    write_csv(out/"plug_flow_heldout_seed_summary.csv",seed_rows)
    p=np.array([r["plug_flux"] for r in rows]);a=np.array([r["actual_neck_flux"] for r in rows]);a0=np.array([r["actual_neck_flux_unshifted"] for r in rows]);e=p-a
    accumulation=np.array([r["approach_accumulation"] for r in rows]);lateral=np.array([r["lateral_outflux"] for r in rows])
    upstream=np.array([r["upstream_outflux"] for r in rows])
    bias_ci=seed_bootstrap_mean(seed_rows,"bias");mae_ci=seed_bootstrap_mean(seed_rows,"mae")
    now_bias_ci=seed_bootstrap_mean(seed_rows,"unshifted_bias");now_mae_ci=seed_bootstrap_mean(seed_rows,"unshifted_mae")
    summary=dict(heldout_seeds=sorted({r["seed"] for r in rows}),n_sustained_events=len(selected),n_samples=len(rows),
                 plug_bias=float(np.mean(e)),plug_mae=float(np.mean(np.abs(e))),plug_rmse=float(np.sqrt(np.mean(e**2))),
                 plug_spearman=float(spearmanr(p,a).statistic),median_travel_time=float(np.median([r["travel_time"] for r in rows])),
                 unshifted_plug_spearman=float(spearmanr(p,a0).statistic),unshifted_plug_bias=float(np.mean(p-a0)),
                 unshifted_plug_mae=float(np.mean(np.abs(p-a0))),
                 residual_vs_accumulation_spearman=float(spearmanr(e,accumulation).statistic),
                 residual_vs_lateral_outflux_spearman=float(spearmanr(e,lateral).statistic),
                 residual_vs_upstream_outflux_spearman=float(spearmanr(e,upstream).statistic),
                 mean_plug_flux=float(np.mean(p)),mean_delayed_neck_flux=float(np.mean(a)),mean_unshifted_neck_flux=float(np.mean(a0)),
                 delayed_zero_flux_fraction=float(np.mean(a==0)),unshifted_zero_flux_fraction=float(np.mean(a0==0)),
                 maximum_conservation_residual=float(np.max(np.abs([r["conservation_residual"] for r in rows]))),
                 seed_mean_bias=bias_ci[0],seed_mean_bias_ci95=list(bias_ci[1:]),seed_mean_mae=mae_ci[0],seed_mean_mae_ci95=list(mae_ci[1:]),
                 seed_median_correlation=float(np.nanmedian([r["correlation"] for r in seed_rows])),
                 seed_mean_unshifted_bias=now_bias_ci[0],seed_mean_unshifted_bias_ci95=list(now_bias_ci[1:]),
                 seed_mean_unshifted_mae=now_mae_ci[0],seed_mean_unshifted_mae_ci95=list(now_mae_ci[1:]),
                 seed_median_unshifted_correlation=float(np.nanmedian([r["unshifted_correlation"] for r in seed_rows])))
    (out/"plug_flow_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return rows,summary


def figures(out: Path, balance: list[dict], closure: list[dict], plug: list[dict]) -> None:
    keys=("boundary","wall","social","noise_relaxation","realized_noise","collision_override")
    fig,axes=plt.subplots(2,2,figsize=(10.5,8.2));axes=axes.ravel()
    x=np.arange(len(keys));ev=[];ct=[]
    for key in keys:
        ev.append(np.mean([r[f"event_{key}"] for r in balance]));ct.append(np.mean([r[f"control_{key}"] for r in balance]))
    axes[0].bar(x-.18,ev,.36,label="future crossing");axes[0].bar(x+.18,ct,.36,label="same-q control");axes[0].set_xticks(x,keys,rotation=45,ha="right");axes[0].set_ylabel("integrated contribution to $\\Delta J_a$");axes[0].legend(frameon=False)
    for sign,color in (("rising","tab:red"),("falling","tab:blue")):
        part=[r for r in closure if r["history_sign"]==sign];axes[1].scatter([r["Ja"] for r in part],[r["future_dJa"] for r in part],s=8,alpha=.25,label=sign,color=color)
    axes[1].axhline(0,color=".5",ls=":");axes[1].set(xlabel="$J_a(t)$",ylabel="future $\\Delta J_a/0.2$");axes[1].legend(frameon=False,title="previous 0.2")
    for ax,key,title in ((axes[2],"actual_neck_flux_unshifted","same 0.2 window"),(axes[3],"actual_neck_flux","after mean travel time")):
        ax.scatter([r["plug_flux"] for r in plug],[r[key] for r in plug],s=12,alpha=.3)
        lim=max(np.quantile([r["plug_flux"] for r in plug],.99),np.quantile([r[key] for r in plug],.99));ax.plot([0,lim],[0,lim],"k:");ax.set(xlabel="$N J_a/\\ell$",ylabel="observed net neck flux",xlim=(0,lim),ylim=(-.2*lim,lim),title=title)
    fig.tight_layout();fig.savefig(out/"figures"/"approach_current_mechanism.png",dpi=220);plt.close(fig)


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument("--raw-dir",type=Path,required=True);p.add_argument("--source",type=Path,required=True)
    p.add_argument("--previous",type=Path,required=True);p.add_argument("--output",type=Path,required=True);p.add_argument("--n-exact",type=int,default=8)
    p.add_argument("--reuse-balance-closure",action="store_true",help="Resume the plug audit from completed compact outputs")
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True);(out/"figures").mkdir(exist_ok=True)
    events=read_csv(args.previous/"crossing_episodes.csv");matched=read_csv(args.previous/"matched_precursor_pairs.csv");pairs=choose_pairs(events,matched)
    if args.reuse_balance_closure:
        balance=numeric_rows(read_csv(out/"native_balance_event_control.csv"));balance_summary=json.loads((out/"native_balance_summary.json").read_text())
        closure=numeric_rows(read_csv(out/"closure_samples.csv"));closure_summary=json.loads((out/"closure_summary.json").read_text())
        closure_summary=augment_closure_summary(closure,closure_summary,out)
    else:
        mod=load_gap20_production(args.source.resolve())
        para=dict(Ip=4.0,In=IN,R=mod.R,light=1.0,dt=DT,Iw=.94,If=0.0,num=N,delta=.01,total_time=1000,gap=mod.gap,radius_ratio=50,gap_ratio=20)
        balance,balance_summary=balance_analysis(mod,para,args.raw_dir,pairs,out,args.n_exact)
        closure,closure_summary=closure_analysis(args.raw_dir,pairs,out)
    plug,plug_summary=plug_flow_analysis(args.raw_dir,events,out)
    figures(out,balance,closure,plug)
    config=dict(raw_dir=str(args.raw_dir),production_source=str(args.source.resolve()),source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                geometry_override="in-memory parameters_list [50,5,20] -> [50,20,20]; production file unchanged",N=N,dt=DT,In=IN,
                approach_length=ELL,approach_width=WIDTH,n_all_pairs=len(pairs),n_native_exact_pairs=args.n_exact,
                native_balance_interval="event A-2 to A; same-q control t to t+2",closure_sampling_dt=.1,
                heldout_plug_seeds="30--39",plug_evaluation_window=.2)
    (out/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    combined=dict(balance=balance_summary,closure=closure_summary,plug_flow=plug_summary)
    (out/"summary.json").write_text(json.dumps(combined,indent=2)+"\n")


if __name__=="__main__": main()
