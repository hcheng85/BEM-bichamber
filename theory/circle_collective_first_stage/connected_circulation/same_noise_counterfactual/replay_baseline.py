#!/usr/bin/env python3
"""First gate of a bounded counterfactual: source-exact same-noise baseline.

Extracts the actual production one-step statements via AST, changing only the
noise provider to the supplied archived increment. Does not edit production.
"""
from pathlib import Path
import ast,copy,csv,hashlib,json,sys,time,argparse
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
sys.path[:0]=[str(BASE),str(BASE/'com_relative_turning')]
from analyze_com_relative_turning import COMReplay

DT=.001;R=7.158910531638177;START=309000;END=320000
POSITION_TOL=1e-6*R;HEADING_TOL=1e-6

def save(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def source_step(mod,source,out):
    tree=ast.parse(source.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='step_Computation')
    loop=next(x for x in fn.body if isinstance(x,ast.For))
    first=next(i for i,x in enumerate(loop.body) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='CurState')
    body=copy.deepcopy(loop.body[first:])
    replacements=0
    for node in body:
        if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='wNoise':
            assert isinstance(node.value,ast.Call) and node.value.func.id=='GetNoise'
            node.value=ast.Name(id='noise_increment',ctx=ast.Load());replacements+=1
    assert replacements==1
    prefix=ast.parse('''def archived_step(current, para, noise_increment, dis_last):
    state = np.empty((3, para["num"], 2))
    state[:,:,0] = current
    rDotList = np.empty((2, para["num"], 2))
    count = 0
''').body[0]
    suffix=ast.parse('''return state[:,:,count].copy(), rDot.copy(), thetaDot.copy(), dis_last.copy(), dict(wall_clip=int(np.sum(abs(wWall)>=np.pi/para["dt"])), angular_clip=int(np.sum(abs(thetaDot)>=np.pi/para["dt"])), collision=int(if_jump.sum()), safety_candidates=int(np.sum(dis>=safedis)))''').body
    prefix.body+=body+suffix
    module=ast.fix_missing_locations(ast.Module(body=[prefix],type_ignores=[]))
    text=ast.unparse(module)+'\n';(out/'extracted_production_step.py').write_text(text)
    env=dict(vars(mod));exec(compile(module,str(source)+'::archived_step','exec'),env)
    return env['archived_step'],hashlib.sha256(text.encode()).hexdigest()

def coordinates(st):
    xy=st[:2]-st[:2].mean(axis=1)[:,None];phi=np.arctan2(xy[1],xy[0]);e=np.exp(1j*st[2]);f=np.exp(-1j*phi)
    p=e.mean();B=(e*f).mean();C=B-p*f.mean()
    return dict(K=float(C.imag),P=float(abs(p)),MC=float(B.imag),com_radius_R=float(np.linalg.norm(st[:2].mean(axis=1))/R),mean_school_radius=float(np.linalg.norm(xy,axis=0).mean()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,required=True);args=ap.parse_args()
    out=HERE/'outputs';out.mkdir(parents=True,exist_ok=True)
    source=BASE.parents[1].parent/'BEM_circular_parallel_instance.py'
    sha=hashlib.sha256(source.read_bytes()).hexdigest()
    assert sha==json.loads((BASE/'simulation_parameters.json').read_text())['source_sha256']
    replay=COMReplay(args.data_dir,source,BASE/'outputs');sm,rm=replay.arrays(9)
    step,step_sha=source_step(replay.mod,source,out);para=replay.para;zero=np.zeros((3,100));cost=time.monotonic()
    # dis_last stores the previous PRE-safety proposed signed distances.
    _,_,_,initial_dis,_=step(np.asarray(sm[:,:,START-1],float),para,zero,None)
    noise_path=out/'reconstructed_noise.npz';audit_path=out/'reconstruction_audit.json'
    if noise_path.exists() and audit_path.exists():
        z=np.load(noise_path);assert int(z['start'])==START and int(z['end'])==END
        old=json.loads(audit_path.read_text());assert old['source_sha256']==sha and old['step_sha256']==step_sha
        noise=z['increments'];print('Reusing validated reconstruction cache',flush=True)
    else:
        noise=np.empty((END-START,100));rows=[];dis_last=initial_dis.copy();maxpos=maxvel=maxheadinground=0.
        for index,n in enumerate(range(START,END)):
            current=np.asarray(sm[:,:,n],float);target=np.asarray(sm[:,:,n+1],float)
            predicted,vel,omega,dis_last,flags=step(current,para,zero,dis_last)
            # Match the production addition order: theta_next=(theta+dt*Omega)+noise.
            noise[index]=target[2]-predicted[2]
            maxheadinground=max(maxheadinground,float(np.max(abs(predicted[2]+noise[index]-target[2]))))
            pe=float(np.max(abs(predicted[:2]-target[:2])));ve=float(np.max(abs(vel-rm[:,:,n+1])))
            maxpos=max(maxpos,pe);maxvel=max(maxvel,ve)
            if index%10==0 or pe>1e-10 or ve>1e-10:
                rows.append(dict(time=n*DT,max_position_error=pe,max_velocity_error=ve,noise_rms=float(np.sqrt(np.mean(noise[index]**2))),**flags))
            if (index+1)%500==0:print(f'noise {index+1}/{END-START}, elapsed {time.monotonic()-cost:.1f}s, position error {maxpos:.3g}',flush=True)
        np.savez_compressed(noise_path,increments=noise,start=START,end=END,dt=DT,initial_dis=initial_dis)
        save(out/'anchored_reconstruction.csv',rows)
        audit=dict(source_sha256=sha,step_sha256=step_sha,steps=END-START,
            max_anchored_position_error=maxpos,max_anchored_velocity_error=maxvel,max_heading_roundtrip_error=maxheadinground,
            noise_mean=float(noise.mean()),noise_std=float(noise.std()),expected_noise_std=.5*np.sqrt(DT),
            reconstruction_seconds=time.monotonic()-cost,position_tolerance=POSITION_TOL,heading_tolerance=HEADING_TOL)
        audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    # Autonomous replay: no re-anchoring and no altered noise after this point.
    current=np.asarray(sm[:,:,START],float).copy();dis_last=initial_dis.copy();rows=[];failure=None
    start_replay=time.monotonic();maxpos=maxhead=0.
    for index,n in enumerate(range(START,END)):
        increment=zero.copy();increment[2]=noise[index]
        current,vel,omega,dis_last,flags=step(current,para,increment,dis_last)
        target=np.asarray(sm[:,:,n+1],float)
        pe=float(np.max(np.linalg.norm(current[:2]-target[:2],axis=0)))
        he=float(np.max(abs(np.angle(np.exp(1j*(current[2]-target[2]))))))
        maxpos=max(maxpos,pe);maxhead=max(maxhead,he)
        base=coordinates(current);truth=coordinates(target)
        rows.append(dict(time=(n+1)*DT,max_position_error=pe,max_heading_error=he,
            **{'replay_'+k:v for k,v in base.items()},**{'archive_'+k:v for k,v in truth.items()},**flags))
        if pe>POSITION_TOL or he>HEADING_TOL:
            failure=dict(step=n+1,time=(n+1)*DT,position_error=pe,heading_error=he);print('BASELINE FAILED',failure,flush=True);break
        if (index+1)%500==0:print(f'baseline {index+1}/{END-START}, errors {pe:.3g}, {he:.3g}',flush=True)
    save(out/'baseline_replay.csv',rows)
    result=dict(passed=failure is None,failure=failure,baseline_steps=len(rows),max_position_error=maxpos,max_heading_error=maxhead,
        baseline_seconds=time.monotonic()-start_replay,total_this_run_seconds=time.monotonic()-cost,
        position_tolerance=POSITION_TOL,heading_tolerance=HEADING_TOL,counterfactuals_run=False)
    (out/'baseline_verdict.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,ax=plt.subplots(2,1,figsize=(9,6),sharex=True)
    ax[0].semilogy([r['time'] for r in rows],[max(r['max_position_error']/R,1e-17) for r in rows],label='position / R')
    ax[0].semilogy([r['time'] for r in rows],[max(r['max_heading_error'],1e-17) for r in rows],label='heading [rad]')
    ax[0].axhline(1e-6,color='k',ls='--',label='predeclared threshold');ax[0].legend();ax[0].set_ylabel('max fish error')
    for k in ('K','P','MC'):
        ax[1].plot([r['time'] for r in rows],[r['replay_'+k]-r['archive_'+k] for r in rows],label=k)
    ax[1].legend();ax[1].set(xlabel='time',ylabel='replay − archive')
    fig.tight_layout();fig.savefig(out/'baseline_validation.png',dpi=180);plt.close(fig)
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
