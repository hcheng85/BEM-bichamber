#!/usr/bin/env python3
"""Exact finite-step polarization budget for the existing branching ensemble.

This replays exactly the already-prespecified transport branches and noise
sequences.  It adds diagnostics to the unchanged production step; it does not
create a new intervention ensemble.
"""
from __future__ import annotations

from pathlib import Path
import argparse, ast, copy, csv, hashlib, json, multiprocessing as mp, sys, time
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
TRANSPORT = HERE.parent / "stochastic_branching_transport"
sys.path[:0] = [str(BASE), str(BASE / "com_relative_turning"),
                str(HERE.parent / "same_noise_counterfactual"),
                str(HERE.parent / "stochastic_branching"), str(TRANSPORT)]
from analyze_com_relative_turning import COMReplay
from replay_torques import load_production, social_split
from run_branching import coordinate, metric, save_csv

DT = .001
IN = .5
N = 100
STEPS = 1000
TIMES = (309., 310.)
REPS = 4
BRANCHES = ("baseline", "baseline_duplicate") + tuple(
    f"{s}_d{j}" for s in ("plus", "minus") for j in range(3))
METRICS = ("K", "P", "MC", "com_x", "com_y", "com_radius_R",
           "mean_radius", "radial_std", "sigma_major", "sigma_minor", "aspect")
SOURCES = ("alignment", "attraction", "wall", "direct_hydro", "bem_hydro",
           "effective_update_residual", "finite_step_curvature",
           "angular_diffusion", "realized_stochastic")
TAUS = (.05, .1, .25, 1.)
_MOD = _PARA = _STEP = None


def extract_step(mod, source: Path, out: Path):
    """Extract the production update and return already-computed torque arrays."""
    tree = ast.parse(source.read_text())
    fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef)
              and x.name == "step_Computation")
    loop = next(x for x in fn.body if isinstance(x, ast.For))
    first = next(i for i, x in enumerate(loop.body)
                 if isinstance(x, ast.Assign) and isinstance(x.targets[0], ast.Name)
                 and x.targets[0].id == "CurState")
    body = copy.deepcopy(loop.body[first:])
    replacements = 0
    for node in body:
        if (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == "wNoise"):
            node.value = ast.Name(id="noise_increment", ctx=ast.Load())
            replacements += 1
    assert replacements == 1
    prefix = ast.parse('''def archived_step_polarization(current, para, noise_increment, dis_last):
    state = np.empty((3, para["num"], 2))
    state[:,:,0] = current
    rDotList = np.empty((2, para["num"], 2))
    count = 0
''').body[0]
    suffix = ast.parse('''
angular = np.vstack([np.asarray(wHydro).ravel(), np.asarray(w_BEM).ravel(),
                     np.asarray(wWall).ravel(), np.asarray(wVision).ravel(),
                     np.asarray(thetaDot).ravel()])
flags = dict(wall_clip=int(np.sum(abs(wWall) >= np.pi/para["dt"])),
             angular_clip=int(np.sum(abs(thetaDot) >= np.pi/para["dt"])),
             collision=int(if_jump.sum()),
             safety_candidates=int(np.sum(dis >= safedis)))
return state[:,:,count].copy(), rDot.copy(), dis_last.copy(), flags, angular, voro.copy()
''').body
    prefix.body += body + suffix
    module = ast.fix_missing_locations(ast.Module(body=[prefix], type_ignores=[]))
    text = ast.unparse(module) + "\n"
    path = out / "extracted_production_step_polarization.py"
    path.write_text(text)
    env = dict(vars(mod))
    exec(compile(module, str(source) + "::archived_step_polarization", "exec"), env)
    return env["archived_step_polarization"], hashlib.sha256(text.encode()).hexdigest()


def prepare(replay, step):
    """Use the exact perturbations already saved by the transport experiment."""
    sm, _ = replay.arrays(9)
    vectors = np.load(TRANSPORT / "outputs" / "perturbation_vectors.npz")
    states, memory = {}, {}
    zero = np.zeros((3, N))
    for t0 in TIMES:
        n = round(t0 / DT)
        st = np.asarray(sm[:, :, n], float)
        _, _, dis0, _, _, _ = step(np.asarray(sm[:, :, n-1], float),
                                    replay.para, zero, None)
        branches = {"baseline": st.copy(), "baseline_duplicate": st.copy()}
        for sign in ("plus", "minus"):
            for direction in range(3):
                name = f"{sign}_d{direction}"
                x = st.copy()
                x[2] += vectors[f"t{int(t0)}_{name}"]
                branches[name] = x
        states[t0], memory[t0] = branches, dis0
    return states, memory


def init_worker(source, extracted):
    global _MOD, _PARA, _STEP
    _MOD = load_production(Path(source))
    _PARA = dict(Ip=4., In=IN, R=_MOD.R, light=1., dt=DT, Iw=.94, If=.01,
                 num=N, delta=.01, total_time=1000, gap=_MOD.gap,
                 radius_ratio=50, gap_ratio=50)
    env = dict(vars(_MOD))
    exec(Path(extracted).read_text(), env)
    _STEP = env["archived_step_polarization"]


def task_worker(task):
    t0, rep, noise_seed, states, dis0, expected_file = task
    rng = np.random.default_rng(noise_seed)
    noise = IN * np.sqrt(DT) * rng.normal(size=(STEPS, N))
    current = {k: v.copy() for k, v in states.items()}
    memory = {k: dis0.copy() for k in BRANCHES}
    data = np.empty((len(BRANCHES), STEPS + 1, len(METRICS)))
    p = np.empty((len(BRANCHES), STEPS + 1), complex)
    dp = np.empty((len(BRANCHES), STEPS, len(SOURCES)), complex)
    closure = np.empty((len(BRANCHES), STEPS), complex)
    flags = {k: dict(wall_clip=0, angular_clip=0, collision=0,
                     safety_candidates=0) for k in BRANCHES}
    for j, name in enumerate(BRANCHES):
        data[j, 0] = metric(current[name])
        p[j, 0] = np.mean(np.exp(1j * current[name][2]))
    duplicate_max = 0.
    start = time.monotonic()
    decay = np.exp(-.5 * IN * IN * DT)
    for n in range(STEPS):
        inc = np.zeros((3, N))
        inc[2] = noise[n]
        for j, name in enumerate(BRANCHES):
            old = current[name]
            new, _, memory[name], fl, angular, voro = _STEP(
                old, _PARA, inc, memory[name])
            direct, bem, wall, vision, effective = angular
            align, attraction = social_split(_MOD, voro, old, _PARA["Ip"])
            if np.max(np.abs(align + attraction - vision)) > 1e-10:
                raise RuntimeError("social split failed")
            omega = (align, attraction, wall, direct, bem)
            z = np.exp(1j * old[2])
            first = np.array([DT * np.mean(1j * z * w) for w in omega])
            raw = sum(omega)
            update_residual = DT * np.mean(1j * z * (effective - raw))
            p_old = z.mean()
            p_star = np.mean(np.exp(1j * (old[2] + effective * DT)))
            curvature = p_star - p_old - DT * np.mean(1j * z * effective)
            diffusion = (decay - 1.) * p_star
            p_next = np.mean(np.exp(1j * new[2]))
            stochastic = p_next - decay * p_star
            pieces = np.r_[first, update_residual, curvature, diffusion, stochastic]
            dp[j, n] = pieces
            closure[j, n] = p_next - p_old - pieces.sum()
            current[name] = new
            data[j, n+1] = metric(new)
            p[j, n+1] = p_next
            for key in flags[name]:
                flags[name][key] += fl[key]
        duplicate_max = max(duplicate_max, float(np.max(np.abs(
            current["baseline"] - current["baseline_duplicate"]))))
    expected = np.load(expected_file, allow_pickle=False)["data"]
    metric_replay_max = float(np.max(np.abs(data - expected)))
    return dict(t0=t0, replicate=rep, noise_seed=noise_seed, data=data, p=p,
                dp=dp, closure=closure, duplicate_max=duplicate_max,
                metric_replay_max=metric_replay_max, flags=flags,
                seconds=time.monotonic()-start)


def _complex_parts(row, prefix, value):
    row[prefix + "_real"] = float(np.real(value))
    row[prefix + "_imag"] = float(np.imag(value))


def summarize(out: Path):
    runs = {}
    audit_rows, interval_rows, transport_rows = [], [], []
    for path in sorted(out.glob("trajectory_t*_r*.npz")):
        z = np.load(path, allow_pickle=False)
        key = (float(z["t0"]), int(z["replicate"]))
        runs[key] = (z["data"], z["p"], z["dp"])
        audit_rows.append(dict(time=key[0], replicate=key[1],
            max_step_closure=float(np.max(np.abs(z["closure"])))) )

    for t0 in TIMES:
        for rep in range(REPS):
            data, p, dp = runs[(t0, rep)]
            cumulative = np.zeros((len(BRANCHES), STEPS+1, len(SOURCES)), complex)
            cumulative[:, 1:] = np.cumsum(dp, axis=1)
            for j, branch in enumerate(BRANCHES[2:], 2):
                sign = branch.split("_")[0]
                direction = int(branch[-1])
                er = complex(data[0, 0, 3], data[0, 0, 4])
                er /= max(abs(er), 1e-30)
                heading = p[0, 0] / max(abs(p[0, 0]), 1e-30)
                for tau in TAUS:
                    n = int(round(tau/DT))
                    response = p[j, n] - p[0, n]
                    r = dict(time=t0, replicate=rep, branch=branch, sign=sign,
                             direction=direction, tau=tau,
                             inward_total=float(-np.real(np.conj(er)*response)),
                             longitudinal_total=float(np.real(np.conj(heading)*response)),
                             transverse_total=float(np.imag(np.conj(heading)*response)),
                             delta_P=float(abs(p[j,n])-abs(p[0,n])))
                    _complex_parts(r, "observed_delta_p", response)
                    for q, source in enumerate(SOURCES):
                        value = cumulative[j,n,q] - cumulative[0,n,q]
                        _complex_parts(r, source, value)
                        r[source+"_inward"] = float(-np.real(np.conj(er)*value))
                    interval_rows.append(r)

                    # Exact vector swimming response as a double integral of dp.
                    dc = cumulative[j, :n] - cumulative[0, :n]
                    swim_parts = DT * dc.sum(axis=0)
                    swim_observed = DT * np.sum(p[j,:n] - p[0,:n])
                    # Exact chord projection used by the validated COM budget:
                    # g_n=(C_n+C_{n+1})/(|C_n|+|C_{n+1}|).
                    cb0=data[j,:n][:,[3,4]];cb1=data[j,1:n+1][:,[3,4]]
                    gb=(cb0+cb1)/np.maximum((np.linalg.norm(cb0,axis=1)+np.linalg.norm(cb1,axis=1))[:,None],1e-30)
                    g_branch=gb[:,0]+1j*gb[:,1]
                    c00=data[0,:n][:,[3,4]];c01=data[0,1:n+1][:,[3,4]]
                    g0=(c00+c01)/np.maximum((np.linalg.norm(c00,axis=1)+np.linalg.norm(c01,axis=1))[:,None],1e-30)
                    g_base=g0[:,0]+1j*g0[:,1]
                    radial_parts = DT*np.real(np.conj(g_branch[:,None])*dc).sum(axis=0)
                    geometry = DT*np.real(np.conj(g_branch-g_base)*p[0,:n]).sum()
                    radial_observed = DT*np.sum(np.real(np.conj(g_branch)*p[j,:n])-
                                                np.real(np.conj(g_base)*p[0,:n]))
                    tr = dict(time=t0,replicate=rep,branch=branch,sign=sign,
                              direction=direction,tau=tau,
                              swimming_radial_observed=float(radial_observed),
                              geometry_projection_residual=float(geometry),
                              radial_reconstruction_error=float(radial_observed-geometry-radial_parts.sum()))
                    _complex_parts(tr,"swimming_vector_observed",swim_observed)
                    for q,source in enumerate(SOURCES):
                        _complex_parts(tr,source+"_swimming_vector",swim_parts[q])
                        tr[source+"_swimming_radial"] = float(radial_parts[q])
                    transport_rows.append(tr)
    save_csv(out/"finite_step_audit.csv", audit_rows)
    save_csv(out/"polarization_interval_budget.csv", interval_rows)
    save_csv(out/"source_resolved_swimming_transport.csv", transport_rows)

    summary=[]
    keys=["inward_total"]+[s+"_inward" for s in SOURCES]
    for t0 in TIMES:
        for sign in ("plus","minus"):
            for tau in TAUS:
                use=[r for r in interval_rows if r["time"]==t0 and r["sign"]==sign and r["tau"]==tau]
                # Average the three prespecified directions inside each noise block.
                row=dict(time=t0,sign=sign,tau=tau,n_noise=REPS,n_directions=3)
                for k in keys:
                    vals=np.array([np.mean([r[k] for r in use if r["replicate"]==rep]) for rep in range(REPS)])
                    row[k+"_mean"]=float(vals.mean());row[k+"_noise_sd"]=float(vals.std(ddof=1))
                    dmeans=np.array([np.mean([r[k] for r in use if r["direction"]==d]) for d in range(3)])
                    row[k+"_direction_range"]=float(dmeans.max()-dmeans.min())
                summary.append(row)
    save_csv(out/"polarization_source_summary.csv",summary)

    transport_summary=[]
    tkeys=[s+"_swimming_radial" for s in SOURCES]+["geometry_projection_residual","swimming_radial_observed"]
    for t0 in TIMES:
        for sign in ("plus","minus"):
            for tau in TAUS:
                use=[r for r in transport_rows if r["time"]==t0 and r["sign"]==sign and r["tau"]==tau]
                row=dict(time=t0,sign=sign,tau=tau,n_noise=REPS,n_directions=3)
                for k in tkeys:
                    vals=np.array([np.mean([r[k] for r in use if r["replicate"]==rep]) for rep in range(REPS)])
                    row[k+"_mean"]=float(vals.mean());row[k+"_noise_sd"]=float(vals.std(ddof=1))
                transport_summary.append(row)
    save_csv(out/"swimming_transport_source_summary.csv",transport_summary)

    # Decisive figures: early p response and the induced swimming displacement.
    colors=dict(alignment="C0",attraction="C2",wall="C3",direct_hydro="C1",
                bem_hydro="C4",effective_update_residual="C7",
                finite_step_curvature="C8",angular_diffusion="C6",
                realized_stochastic="C5")
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for col,t0 in enumerate(TIMES):
        for row,sign in enumerate(("plus","minus")):
            ax=axes[row,col]
            for source in SOURCES:
                yy=[next(x for x in summary if x["time"]==t0 and x["sign"]==sign and x["tau"]==tau)[source+"_inward_mean"] for tau in TAUS]
                ax.plot(TAUS,yy,"o-",label=source if (row,col)==(0,0) else None,color=colors[source],lw=1)
            yy=[next(x for x in summary if x["time"]==t0 and x["sign"]==sign and x["tau"]==tau)["inward_total_mean"] for tau in TAUS]
            ax.plot(TAUS,yy,"ko-",lw=2,label="observed total" if (row,col)==(0,0) else None)
            ax.axhline(0,color="k",lw=.5);ax.grid(alpha=.2);ax.set_title(f"t0={t0:g}, {sign}K")
            ax.set_ylabel("inward component of Delta p")
    for ax in axes[1]:ax.set_xlabel("elapsed time")
    axes[0,0].legend(fontsize=7,ncol=2);fig.tight_layout()
    fig.savefig(out/"polarization_source_response.png",dpi=180);plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
    for col,t0 in enumerate(TIMES):
        for row,sign in enumerate(("plus","minus")):
            ax=axes[row,col]
            for source in SOURCES:
                yy=[next(x for x in transport_summary if x["time"]==t0 and x["sign"]==sign and x["tau"]==tau)[source+"_swimming_radial_mean"] for tau in TAUS]
                ax.plot(TAUS,yy,"o-",label=source if (row,col)==(0,0) else None,color=colors[source],lw=1)
            geom=[next(x for x in transport_summary if x["time"]==t0 and x["sign"]==sign and x["tau"]==tau)["geometry_projection_residual_mean"] for tau in TAUS]
            obs=[next(x for x in transport_summary if x["time"]==t0 and x["sign"]==sign and x["tau"]==tau)["swimming_radial_observed_mean"] for tau in TAUS]
            ax.plot(TAUS,geom,"--",color="0.5",label="projection geometry" if (row,col)==(0,0) else None)
            ax.plot(TAUS,obs,"ko-",lw=2,label="observed swimming" if (row,col)==(0,0) else None)
            ax.axhline(0,color="k",lw=.5);ax.grid(alpha=.2);ax.set_title(f"t0={t0:g}, {sign}K")
            ax.set_ylabel("radial swimming displacement")
    for ax in axes[1]:ax.set_xlabel("elapsed time")
    axes[0,0].legend(fontsize=7,ncol=2);fig.tight_layout()
    fig.savefig(out/"source_resolved_swimming_response.png",dpi=180);plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-dir",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--stage",choices=("all","setup","simulate","summary"),default="all")
    args=ap.parse_args()
    out=HERE/"outputs";out.mkdir(parents=True,exist_ok=True)
    source=BASE.parents[1].parent/"BEM_circular_parallel_instance.py"
    source_sha=hashlib.sha256(source.read_bytes()).hexdigest()
    expected_sha=json.loads((BASE/"simulation_parameters.json").read_text())["source_sha256"]
    if source_sha != expected_sha: raise RuntimeError("production source hash changed")
    replay=COMReplay(args.data_dir,source,BASE/"outputs")
    step,step_sha=extract_step(replay.mod,source,out)
    states,memory=prepare(replay,step)
    config=dict(source_sha256=source_sha,extracted_step_sha256=step_sha,dt=DT,
                steps=STEPS,initial_times=TIMES,replicates=REPS,branches=BRANCHES,
                sources=SOURCES,noise_seeds="identical to stochastic_branching_transport",
                scope="source accounting of the existing branch ensemble; no new interventions")
    (out/"analysis_configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    if args.stage=="setup":return
    tasks=[]
    for t0 in TIMES:
        for rep in range(REPS):
            expected=TRANSPORT/"outputs"/f"trajectory_t{int(t0)}_r{rep}.npz"
            tasks.append((t0,rep,960000+int(t0)*100+rep,states[t0],memory[t0],expected))
    if args.stage in ("all","simulate"):
        todo=[x for x in tasks if not (out/f"trajectory_t{int(x[0])}_r{x[1]}.npz").exists()]
        results=[]
        if todo:
            start=time.monotonic();ctx=mp.get_context("spawn")
            with ctx.Pool(min(args.workers,len(todo)),initializer=init_worker,
                          initargs=(str(source),str(out/"extracted_production_step_polarization.py"))) as pool:
                for z in pool.imap_unordered(task_worker,todo):
                    path=out/f"trajectory_t{int(z['t0'])}_r{z['replicate']}.npz"
                    np.savez_compressed(path,t0=z["t0"],replicate=z["replicate"],
                        noise_seed=z["noise_seed"],data=z["data"],p=z["p"],dp=z["dp"],
                        closure=z["closure"],branches=np.array(BRANCHES),sources=np.array(SOURCES),
                        metrics=np.array(METRICS))
                    results.append({k:v for k,v in z.items() if k not in ("data","p","dp","closure","flags")}|{"flags":z["flags"]})
                    print(f"finished t={z['t0']}, rep={z['replicate']}, replay={z['metric_replay_max']:.3g}",flush=True)
            (out/"simulation_audit.json").write_text(json.dumps(dict(workers=args.workers,
                wall_seconds=time.monotonic()-start,runs=results),indent=2)+"\n")
        audit=json.loads((out/"simulation_audit.json").read_text())
        if max(x["duplicate_max"] for x in audit["runs"])>1e-13:raise RuntimeError("duplicate failure")
        if max(x["metric_replay_max"] for x in audit["runs"])>1e-12:raise RuntimeError("stored-trajectory replay failure")
    if args.stage=="simulate":return
    summarize(out)


if __name__=="__main__":main()
