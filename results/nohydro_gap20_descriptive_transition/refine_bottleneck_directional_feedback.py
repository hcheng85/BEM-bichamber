#!/usr/bin/env python3
"""Native-wall directional feedback experiment, with gross crossing accounting.

Place beside run_bottleneck_commitment_branches.py in the archived repo.
No production source or data are changed.  The x=0 aperture classification is
DIAGNOSTIC, never an extra wall predicate.  See companion protocol.
"""
from __future__ import annotations
import argparse
import ast
import contextlib
import hashlib
import io
import json
import multiprocessing as mp
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PILOT = HERE / 'outputs' / 'bottleneck_commitment_branches'
OUT = HERE / 'outputs' / 'bottleneck_directional_refined'
G = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def geometry(state, vn, gap):
    """Native-admitted opposite-sign links; no assertion of wall obstruction."""
    x, y = np.asarray(state[:2], float)
    v = np.asarray(vn, bool)
    cross = v & (x[:, None]*x[None, :] < 0)  # exact legacy ablation mask
    den = x[None, :] - x[:, None]
    frac = np.divide(-x[:, None], den, out=np.zeros_like(den), where=den != 0)
    yy = y[:, None] + frac*(y[None, :]-y[:, None])
    central = cross & (np.abs(yy) <= gap/2 + 1e-10)
    return cross, central, cross & ~central, yy


def source_mask(state, sign):
    if sign not in (-1, 1):
        raise ValueError('source sign must be +/-1')
    return state[0] >= 0 if sign == 1 else state[0] < 0


def crossing_counts(prev, nxt, sign, gap):
    """Exact per-fish discrete sign changes, locations by step interpolation."""
    old, new = source_mask(prev, sign), source_mask(nxt, sign)
    forward, backward = old & ~new, ~old & new
    change = forward | backward
    ids = np.flatnonzero(change)
    dx = nxt[0, ids]-prev[0, ids]
    if np.any(dx == 0):
        raise RuntimeError('Sign change without displacement')
    yy = prev[1, ids] - prev[0, ids]*(nxt[1, ids]-prev[1, ids])/dx
    iscentral = np.abs(yy) <= gap/2 + 1e-10
    return forward, backward, ids, yy, iscentral


def arms(scope):
    # The original all-off arm is always kept for end-point parity.
    if scope == 'native':
        return ('full', 'crossgap_off', 'source_receives_off', 'target_receives_off')
    return ('full', 'crossgap_off', 'central_off', 'central_source_off', 'central_target_off')


def masked_social(base, mod, vn, state, para, mask, full):
    """Subset numerator deletion; keep original full visual denominator."""
    if not np.any(mask):
        return full.copy()
    kept = vn.copy(); kept[mask] = False
    with np.errstate(invalid='ignore', divide='ignore'), contextlib.redirect_stdout(io.StringIO()):
        kept_torque = np.asarray(mod.GetVision(kept, state, para)).reshape(-1)
    old_den = base.vision_denominator(mod, vn, state, para)
    new_den = base.vision_denominator(mod, kept, state, para)
    reduced = np.divide(kept_torque*new_den, old_den, out=np.zeros_like(full), where=old_den != 0)
    # Only receivers with removed edges may change: retain native zero-degree fallback exactly.
    return np.where(mask.any(axis=0) | mask.any(axis=1), reduced, full)


def social(base, mod, vn, state, para, mode, sign, stats=None):
    cross, central, other, _ = geometry(state, vn, float(mod.gap))
    if stats is not None:
        stats['native_link_steps'] += int(cross.sum())
        stats['central_link_steps'] += int(central.sum())
        stats['other_link_steps'] += int(other.sum())
        stats['steps_with_other_links'] += int(other.any())
        stats['zero_plane_links'] += int(np.sum(vn & ((state[0,:,None]==0) | (state[0,None,:]==0))))
    full = base.social_with_filter(mod, vn, state, para, False)
    if mode == 'full':
        return full
    if mode in ('crossgap_off', 'source_receives_off', 'target_receives_off'):
        off = base.social_with_filter(mod, vn, state, para, True)
    elif mode in ('central_off', 'central_source_off', 'central_target_off'):
        # Use legacy implementation when the masks coincide exactly.
        off = (base.social_with_filter(mod, vn, state, para, True)
               if np.array_equal(cross, central) else masked_social(base, mod, vn, state, para, central, full))
    else:
        raise ValueError(mode)
    if mode in ('crossgap_off', 'central_off'):
        return off
    src = source_mask(state, sign)
    receiver = src if mode in ('source_receives_off', 'central_source_off') else ~src
    return np.where(receiver, off, full)


def structural_audit(base, mod, para, origins, source, out):
    """Check actual native predicates; give Codex exact segment endpoints to plot.

    This is not independent verification of completeness of physical wall model.
    A native-admitted off-central line is not labelled solid.
    """
    rows, links, points = [], [], []
    for t in origins:
        st = t['state']; vn0 = mod.GetVoronoi(st[:2], para)
        native = mod.CheckVoronoi(vn0.copy(), st[:2], para)
        fast = base.fast_check_voronoi(mod, vn0.copy(), st[:2])
        if not np.array_equal(native, fast):
            raise RuntimeError('Native/fast visibility mismatch '+t['sample_id'])
        cross, cen, other, yy = geometry(st, native, float(mod.gap))
        full = social(base,mod,native,st,para,'full',t['source_sign'])
        off = social(base,mod,native,st,para,'crossgap_off',t['source_sign'])
        src = social(base,mod,native,st,para,'source_receives_off',t['source_sign'])
        tar = social(base,mod,native,st,para,'target_receives_off',t['source_sign'])
        clo = float(np.max(np.abs((full-src)+(full-tar)-(full-off))))
        if not np.all(np.isfinite([full, off, src, tar])) or clo > 1e-11:
            raise RuntimeError('Torque partition failed '+t['sample_id'])
        rows.append(dict(sample_id=t['sample_id'],phase=t['phase'],seed=t['seed'],
                         native_cross_directed=int(cross.sum()),central_directed=int(cen.sum()),
                         other_directed=int(other.sum()),receiver_closure=clo))
        if not other.any():
            continue
        for k in range(st.shape[1]):
            points.append(dict(sample_id=t['sample_id'],fish=k,x=st[0,k],y=st[1,k],theta=st[2,k]))
        for i,j in np.argwhere(other):
            a,b=st[:2,i],st[:2,j]
            wall_seg=bool(mod.check_intersections(mod.segments,[a,b]))
            wall5=bool(mod.check_inersect_circle(a,b,mod.circle_center5,mod.r))
            wall6=bool(mod.check_inersect_circle(a,b,mod.circle_center6,mod.r))
            links.append(dict(sample_id=t['sample_id'],phase=t['phase'],seed=t['seed'],fish_i=int(i),fish_j=int(j),
                              xi=float(a[0]),yi=float(a[1]),xj=float(b[0]),yj=float(b[1]),y_at_x0=float(yy[i,j]),
                              native_visible=True,segment_predicate=wall_seg,tip5_predicate=wall5,tip6_predicate=wall6))
            if wall_seg or wall5 or wall6:
                raise RuntimeError('Admitted link conflicts with native predicate')
    from analyze_directional_bottleneck_feedback import csv_write
    csv_write(out/'native_visibility_audit.csv',rows)
    csv_write(out/'offcentral_link_endpoints.csv',links)
    csv_write(out/'offcentral_snapshots.csv',points)
    # Archive exact native geometry values and relevant source functions for review.
    geom={}
    for name,v in vars(mod).items():
        if name == 'segments' or name in ('R','r','gap') or name.startswith('circle_center'):
            try: geom[name]=np.asarray(v,dtype=float).tolist()
            except (ValueError,TypeError): pass
    dump_json(out/'native_geometry_values.json',geom)
    txt=source.read_text();tree=ast.parse(txt)
    wanted={'CheckVoronoi','check_intersections','check_inersect_circle','GetWallAvoid','step_Computation','GetVoronoi','GetVision'}
    pieces=[ast.get_source_segment(txt,n) for n in tree.body if isinstance(n,ast.FunctionDef) and (n.name in wanted or any(s in n.name.lower() for s in ('wall','distance','collis','inside','intersect')))]
    (out/'native_geometry_source_excerpt.txt').write_text('\n\n'.join(p for p in pieces if p)+'\n')
    summary=dict(n_snapshots=len(rows),n_cross_directed=sum(r['native_cross_directed'] for r in rows),
                 n_central_directed=sum(r['central_directed'] for r in rows),n_other_directed=sum(r['other_directed'] for r in rows),
                 offcentral_examples=sorted({r['sample_id'] for r in links}),native_predicates_consistent=True,
                 independent_physical_wall_check='REQUIRES source-grounded local review; off-central is not synonymous with obstructed')
    dump_json(out/'native_geometry_audit.json',summary)
    return summary


def audit_indexed_crossings(base, raw, gap, out):
    """Bounded read-only check of archived indexed first-entry frames (not all events)."""
    from analyze_directional_bottleneck_feedback import csv_read,csv_write
    path=HERE/'outputs'/'initiation_propagation'/'crossing_episodes.csv'
    if not path.is_file():raise FileNotFoundError(path)
    events=csv_read(path);groups=defaultdict(list)
    for r in events:groups[int(r['seed'])].append(r)
    rows=[];exceptions=[]
    for seed,rs in sorted(groups.items()):
        rr=base.RawReader(base.raw_path(raw,seed))
        try:
            for r in sorted(rs,key=lambda r:int(r['first_native'])):
                k=int(r['first_native']);sign=-int(r['direction'])
                block=rr.block('state',k-1,k);pre,nxt=block[:,:,0],block[:,:,1]
                f,b,ids,yy,cen=crossing_counts(pre,nxt,sign,gap)
                qpre=2*int((pre[0]>=0).sum())/pre.shape[1]-1
                rec=dict(seed=seed,event_id=r['event_id'],native_step=k,source_sign=sign,
                         forward=int(f.sum()),backward=int(b.sum()),other=int((~cen).sum()),
                         q_pre_error=float(abs(qpre-float(r['q_pre']))),
                         max_abs_y=float(np.max(np.abs(yy))) if len(yy) else None)
                rows.append(rec)
                for fish,y0 in zip(ids[~cen],yy[~cen]):
                    exceptions.append(dict(seed=seed,event_id=r['event_id'],native_step=k,fish=int(fish),y_at_x0=float(y0),
                                           x_before=float(pre[0,fish]),y_before=float(pre[1,fish]),
                                           x_after=float(nxt[0,fish]),y_after=float(nxt[1,fish])))
        finally:rr.close()
    csv_write(out/'indexed_first_crossing_audit.csv',rows)
    csv_write(out/'indexed_offcentral_motion.csv',exceptions)
    meta=dict(n_indexed_events=len(rows),n_forward_events=sum(r['forward']>0 for r in rows),
              n_offcentral_particle_crossings=sum(r['other'] for r in rows),
              n_total_particle_crossings=sum(r['forward']+r['backward'] for r in rows),
              max_q_pre_error=max(r['q_pre_error'] for r in rows),
              n_direction_mismatches=sum(r['forward']<=r['backward'] for r in rows),
              scope='indexed first-crossing frames only, not exhaustive trajectory audit')
    dump_json(out/'indexed_first_crossing_audit.json',meta)
    return meta


def init_worker(source,adapter_text):
    import run_bottleneck_commitment_branches as base
    from analyze_approach_current_balance import load_gap20_production,fast_check_voronoi
    mod=load_gap20_production(Path(source));mod.CheckVoronoi=lambda vn,loc,pa:fast_check_voronoi(mod,vn,loc)
    env=dict(vars(mod));para=base.parameters(mod)
    env['GetVision_intervention']=lambda vn,st,pa,mode:social(base,mod,vn,st,pa,mode,G['sign'],G['stats'])
    exec(adapter_text,env)
    G.update(base=base,step=env['archived_step'],para=para,gap=float(mod.gap))


def run_task(t):
    b, step, para, gap = [G[k] for k in ('base','step','para','gap')]
    sign=int(t['source_sign']);G['sign']=sign
    start=t['state'];nr0=int((start[0]>=0).sum());scope=t['scope'];aa=arms(scope)
    states={a:start.copy() for a in aa};memory={a:t['dis_last'].copy() for a in aa}
    counts={a:dict(forward=0,backward=0,forward_central=0,backward_central=0,forward_other=0,backward_other=0,
                   native_link_steps=0,central_link_steps=0,other_link_steps=0,steps_with_other_links=0,zero_plane_links=0) for a in aa}
    seen_forward={a:np.zeros(b.N,bool) for a in aa};seen_backward={a:np.zeros(b.N,bool) for a in aa}
    initial_source=source_mask(start,sign)
    rng=np.random.default_rng(int(t['noise_seed']))
    noise=0.5*np.sqrt(b.DT)*rng.normal(size=(b.HORIZON_STEPS,b.N));inc=np.zeros((3,b.N))
    for z in noise:
        inc[2]=z
        for a in aa:
            prev=states[a];G['stats']=counts[a]
            nxt,_,_,memory[a],_=step(prev,para,inc,memory[a],a)
            f,back,ids,yy,cen=crossing_counts(prev,nxt,sign,gap)
            c=counts[a];c['forward']+=int(f.sum());c['backward']+=int(back.sum())
            c['forward_central']+=int(np.sum(f[ids]&cen));c['backward_central']+=int(np.sum(back[ids]&cen))
            c['forward_other']+=int(np.sum(f[ids]&~cen));c['backward_other']+=int(np.sum(back[ids]&~cen))
            seen_forward[a]|=f;seen_backward[a]|=back;states[a]=nxt
    rows=[]
    for a in aa:
        nr=int((states[a][0]>=0).sum());net=-sign*(nr-nr0);c=counts[a]
        if c['forward']-c['backward']!=net or c['forward_central']+c['forward_other']!=c['forward'] or c['backward_central']+c['backward_other']!=c['backward']:
            raise RuntimeError('Exact gross crossing balance failed')
        rows.append(dict(phase='post',sample_id=t['sample_id'],seed=int(t['seed']),clone=int(t['clone']),arm=a,
                         source_sign=sign,nright0=nr0,nrightT=nr,additional_net=int(net),success4=int(net>=4),
                         unique_initial_source_forward=int(np.sum(seen_forward[a]&initial_source)),
                         unique_initial_target_backward=int(np.sum(seen_backward[a]&~initial_source)),**c))
    return rows


def summary(rows,scope):
    from analyze_directional_bottleneck_feedback import mean_ci
    aa=arms(scope);index={(r['sample_id'],int(r['clone']),r['arm']):r for r in rows}
    if len(index)!=len(rows):raise ValueError('Duplicate branch row')
    pairs=sorted({(r['sample_id'],int(r['clone'])) for r in rows})
    seeds=defaultdict(list)
    outcomes=('additional_net','forward','backward','unique_initial_source_forward','unique_initial_target_backward','success4')
    for sid,clone in pairs:
        full=index[sid,clone,'full'];row={}
        for a in aa:
            r=index[sid,clone,a]
            if int(r['seed'])!=int(full['seed']) or int(r['nright0'])!=int(full['nright0']) or int(r['source_sign'])!=int(full['source_sign']):
                raise ValueError('Pair mismatch')
            for k in outcomes:row[f'{a}_{k}']=float(r[k])
        seeds[int(full['seed'])].append(row)
    if len(pairs)!=112 or len(seeds)!=14 or any(len(p)!=8 for p in seeds.values()):raise ValueError('Incomplete fixed pilot cohort')
    per=[]
    for seed,rs in sorted(seeds.items()):
        d=dict(seed=seed,n_clones=len(rs))
        for a in aa:
            for k in outcomes:
                key=f'{a}_{k}';d[key]=float(np.mean([r[key] for r in rs]))
                if a!='full':d[f'effect_{key}']=float(np.mean([r[f'full_{k}']-r[key] for r in rs]))
        for a in aa[1:]:
            assert abs(d[f'effect_{a}_additional_net']-(d[f'effect_{a}_forward']-d[f'effect_{a}_backward']))<1e-10
        per.append(d)
    estimates={k:mean_ci([r[k] for r in per]) for k in per[0] if k not in ('seed','n_clones')}
    exposure={a:{k:sum(int(r[k]) for r in rows if r['arm']==a) for k in ('native_link_steps','central_link_steps','other_link_steps','steps_with_other_links','zero_plane_links','forward_other','backward_other')} for a in aa}
    return per,dict(scope=scope,n_seeds=len(seeds),n_pairs=len(pairs),effects_and_baselines=estimates,exposure=exposure,
                    effect_definition='full minus ablated; positive forward means full increases gross forward crossings, negative backward means full suppresses returns')


def self_test():
    # Three distinct changes: forward through central, backwards outside central,
    # and an exact-plane endpoint.  The last follows x>=0 membership convention.
    pre=np.array([[1,-1,0,1],[0,2,0,0],[0,0,0,0]],float)
    nxt=np.array([[-1,1,-1,1],[0,2,0,0],[0,0,0,0]],float)
    f,b,ids,yy,cen=crossing_counts(pre,nxt,1,1.)
    assert f.sum()==2 and b.sum()==1 and cen.tolist()==[True,False,True]
    assert f.sum()-b.sum()==-((nxt[0]>=0).sum()-(pre[0]>=0).sum())
    f2,b2,*_=crossing_counts(pre,nxt,-1,1.)
    assert np.array_equal(f,b2) and np.array_equal(b,f2)
    st=np.array([[1,-1,1,-1],[0,0,2,2],[0,0,0,0]],float)
    vn=np.zeros((4,4),bool);vn[0,1]=vn[1,0]=vn[2,3]=vn[3,2]=True
    cross,cen,other,_=geometry(st,vn,1.)
    assert cross.sum()==4 and cen.sum()==2 and other.sum()==2
    # Mock normalized social kernel. Central deletion must leave the other
    # native-admitted cross-side link unchanged and preserve the denominator.
    class Mod:
        gap=1.
        @staticmethod
        def GetVision(v, st, pa):
            d=v.sum(axis=1)
            return np.divide(v@np.arange(1,5.),d,out=np.zeros(4),where=d>0)
    class Base:
        @staticmethod
        def vision_denominator(m,v,s,p):return v.sum(axis=1)
        @staticmethod
        def social_with_filter(m,v,s,p,off):
            full=m.GetVision(v,s,p)
            if not off:return full
            mask=v&(s[0,:,None]*s[0,None,:]<0)
            return masked_social(Base,m,v,s,p,mask,full)
    full=social(Base,Mod,vn,st,None,'full',1)
    no=social(Base,Mod,vn,st,None,'central_off',1)
    assert np.array_equal(no[2:],full[2:]) and np.all(no[:2]==0)
    for sign in (-1,1):
        a=social(Base,Mod,vn,st,None,'central_source_off',sign)
        d=social(Base,Mod,vn,st,None,'central_target_off',sign)
        assert np.allclose((full-a)+(full-d),full-no)
    # Whole analysis smoke test on the expected 14 x 8 paired cohort.
    rows=[]
    for seed in range(14):
        for clone in range(8):
            for arm in arms('native'):
                fwd,back={'full':(12,3),'crossgap_off':(9,4),'source_receives_off':(9,3),'target_receives_off':(12,4)}[arm]
                rows.append(dict(sample_id=f's{seed}',seed=seed,clone=clone,arm=arm,source_sign=1,nright0=80,
                                 additional_net=fwd-back,forward=fwd,backward=back,unique_initial_source_forward=fwd,
                                 unique_initial_target_backward=back,success4=int(fwd-back>=4),
                                 native_link_steps=8,central_link_steps=7,other_link_steps=1,steps_with_other_links=1,
                                 zero_plane_links=0,forward_other=0,backward_other=0))
    per,res=summary(rows,'native')
    assert len(per)==14 and res['effects_and_baselines']['effect_crossgap_off_additional_net']['mean']==4
    assert res['effects_and_baselines']['effect_source_receives_off_forward']['mean']==3
    assert res['effects_and_baselines']['effect_target_receives_off_backward']['mean']==-1
    print('PASS: neutral geometry classification; both crossing directions; exact count balance; scoped social masks; matched 14-seed summary')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true');p.add_argument('--stage',choices=('audit','simulate'),default='audit')
    p.add_argument('--source',type=Path);p.add_argument('--raw-dir',type=Path);p.add_argument('--pilot-dir',type=Path,default=PILOT)
    p.add_argument('--output',type=Path,default=OUT);p.add_argument('--workers',type=int,default=6)
    p.add_argument('--scope',choices=('native','central'),default='native')
    p.add_argument('--geometry-review',type=Path,help='A source-grounded local review JSON with status=consistent and explanation; needed to simulate')
    args=p.parse_args()
    if args.self_test:return self_test()
    import run_bottleneck_commitment_branches as base
    from analyze_approach_current_balance import load_gap20_production
    from analyze_directional_bottleneck_feedback import csv_read,csv_write,previous_pilot_parity
    source=(args.source or base.SOURCE).resolve();raw=(args.raw_dir or base.RAW_DEFAULT).resolve()
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    old=json.loads((args.pilot_dir/'simulator_audit.json').read_text())
    mod=load_gap20_production(source);para=base.parameters(mod);step,adapter,adapter_sha=base.extract_step(mod,source)
    source_sha=sha(source)
    if source_sha!=old['source_sha256'] or adapter_sha!=old['adapter_sha256']:raise RuntimeError('Source/adapter hash gate failed')
    origins=base.setup_origins(raw,step,para,base.choose_pre(),base.choose_post())
    manifest=csv_read(args.pilot_dir/'selection_manifest.csv')
    expected={(r['phase'],r['sample_id'],int(r['seed']),int(r['native_step']),int(r['source_sign']),int(r['nright0'])) for r in manifest}
    current={(r['phase'],r['sample_id'],int(r['seed']),int(r['native_step']),int(r['source_sign']),int(r['nright0'])) for r in origins}
    if current!=expected or len(current)!=28:raise RuntimeError('Pilot manifest parity failed')
    audit=structural_audit(base,mod,para,origins,source,out)
    audit['indexed_motion']=audit_indexed_crossings(base,raw,float(mod.gap),out)
    print(json.dumps(audit,indent=2),flush=True)
    if args.stage=='audit':return
    if audit['indexed_motion']['n_direction_mismatches'] or audit['indexed_motion']['max_q_pre_error']>1e-10:
        raise RuntimeError('Archived initial crossing indexing mismatch; resolve before branching')
    if args.geometry_review is None:raise RuntimeError('Source-grounded geometry review required, not just x=0 classification')
    review=json.loads(args.geometry_review.read_text())
    if review.get('source_sha256')!=source_sha or review.get('status')!='consistent' or not str(review.get('explanation','')).strip():
        raise RuntimeError('Unresolved or inconsistent actual geometry: stop branching')
    out=out/args.scope;out.mkdir(parents=True,exist_ok=True)
    config=dict(source_sha256=source_sha,adapter_sha256=adapter_sha,scope=args.scope,manifest_sha256=sha(args.pilot_dir/'selection_manifest.csv'),review_sha256=sha(args.geometry_review),script_sha256=sha(Path(__file__)))
    cpath=out/'config.json'
    if cpath.exists() and json.loads(cpath.read_text())!=config:raise RuntimeError('Resume configuration mismatch; use a fresh output directory')
    dump_json(cpath,config)
    tasks=[dict(t,clone=c,scope=args.scope,noise_seed=202611130000+2*100000+int(t['seed'])*100+c) for t in origins if t['phase']=='post' for c in range(base.N_CLONES)]
    partial=out/'branches_partial.csv';rows=[];complete=set()
    if partial.exists():
        groups=defaultdict(list)
        for r in csv_read(partial):groups[(r['sample_id'],int(r['clone']))].append(r)
        valid={(t['sample_id'],int(t['clone'])) for t in tasks}
        for key,rr in groups.items():
            if key in valid and len(rr)==len(arms(args.scope)) and {r['arm'] for r in rr}==set(arms(args.scope)):
                rows.extend(rr);complete.add(key)
    pending=[t for t in tasks if (t['sample_id'],int(t['clone'])) not in complete]
    if pending:
        with mp.get_context('spawn').Pool(min(args.workers,len(pending)),initializer=init_worker,initargs=(str(source),adapter)) as pool:
            for n,rr in enumerate(pool.imap_unordered(run_task,pending),1):
                rows.extend(rr)
                if n%5==0 or n==len(pending):csv_write(partial,rows);print(f'{n}/{len(pending)} blocks',flush=True)
    parity=previous_pilot_parity(rows,csv_read(args.pilot_dir/'branches.csv'))
    dump_json(out/'parity.json',parity)
    if not parity['passed']:raise RuntimeError('Legacy full/all-off endpoint parity failed')
    per,res=summary(rows,args.scope);csv_write(out/'branches.csv',sorted(rows,key=lambda r:(r['sample_id'],int(r['clone']),r['arm'])))
    csv_write(out/'per_seed.csv',per);dump_json(out/'summary.json',res)
    if partial.exists():partial.unlink()
    print(json.dumps(res,indent=2))

if __name__=='__main__':main()
