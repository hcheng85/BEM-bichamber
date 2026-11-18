#!/usr/bin/env python3
"""Directional cross-aperture social feedback experiment (read-only archive).

Run from the repository results/nohydro_gap20_descriptive_transition folder.
This script imports the audited archived restart adapter, preserves its geometry,
noisy timestep and original two arms, and adds two receiver-specific social ablations.
No simulator file or archived MAT is modified.

Scientifically essential gates:
  * wall-visible Delaunay links crossing x=0 must intersect the physical aperture;
  * source and adapter SHA256 must match the already-audited pilot;
  * full and crossgap_off replays must reproduce the committed 112 paired outputs.

See accompanying README. Use --self-test before running on local data.
"""
from __future__ import annotations
import argparse
import contextlib
import csv
import hashlib
import io
import json
import multiprocessing as mp
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT_DEFAULT = HERE / 'outputs' / 'bottleneck_directional_feedback'
PILOT_DEFAULT = HERE / 'outputs' / 'bottleneck_commitment_branches'
ARMS = ('full','crossgap_off','source_receives_off','target_receives_off')
GLOBAL = {}


def csv_read(p):
    with open(p, newline='') as f: return list(csv.DictReader(f))


def csv_write(p, rows):
    if not rows: return
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    keys=list(dict.fromkeys(k for row in rows for k in row))
    tmp=p.with_name(p.name+'.tmp')
    with open(tmp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n')
        w.writeheader();w.writerows(rows)
    os.replace(tmp,p)


def crossing_link_geometry(state, vn, gap, tol=1e-10):
    """Classify directed Delaunay/visual links crossing divider as aperture/solid.

    This audits *visual links already admitted by the native wall filtering*.
    No inference about social torque signs is made here.
    """
    x=np.asarray(state[0],float); y=np.asarray(state[1],float)
    v=np.asarray(vn,bool)
    if v.shape!=(len(x),len(x)): raise ValueError('Expected square adjacency')
    cross=(x[:,None]*x[None,:]<0) & v
    denom=x[None,:]-x[:,None]
    frac=np.divide(-x[:,None],denom,out=np.zeros_like(denom),where=np.abs(denom)>1e-14)
    y_div=y[:,None]+frac*(y[None,:]-y[:,None])
    opening=cross & (np.abs(y_div)<=gap/2+tol)
    solid=cross & ~opening
    return dict(cross=cross,opening=opening,solid=solid,y_intersection=y_div)


def selective_social(base, mod, vn, state, para, mode, source_sign):
    """Remove cross-side social input for specified receiver chamber only.

    Reuses the existing audited full/crossgap_off code and its original visual
    denominator, so no new numerator or collision law is introduced.
    """
    if mode=='full': return base.social_with_filter(mod,vn,state,para,False)
    if mode=='crossgap_off': return base.social_with_filter(mod,vn,state,para,True)
    if mode not in ('source_receives_off','target_receives_off'):
        raise ValueError('Unknown social arm '+str(mode))
    full=base.social_with_filter(mod,vn,state,para,False)
    no_cross=base.social_with_filter(mod,vn,state,para,True)
    if int(source_sign) not in (-1,1): raise ValueError('Bad majority sign')
    x=state[0]
    source=(x>=0) if source_sign==1 else (x<0)
    recipient=source if mode=='source_receives_off' else ~source
    return np.where(recipient,no_cross,full)


def structural_audit(base,mod,para,origins):
    rows=[]; exceed=0; worst=0.; examples=[]
    for item in origins:
        state=item['state'];
        vn0=mod.GetVoronoi(state[:2],para)
        native=mod.CheckVoronoi(vn0.copy(),state[:2],para)
        fast=base.fast_check_voronoi(mod,vn0.copy(),state[:2])
        if not np.array_equal(native,fast):raise RuntimeError('Voronoi parity failed '+item['sample_id'])
        link=crossing_link_geometry(state,native,float(mod.gap))
        n_cross=int(link['cross'].sum());n_open=int(link['opening'].sum());n_solid=int(link['solid'].sum())
        n_unique=int(np.triu(link['cross']|link['cross'].T,1).sum())
        if n_solid:
            exceed+=n_solid
            yy=abs(link['y_intersection'][link['solid']])
            worst=max(worst,float(np.max(yy)-mod.gap/2))
            examples.append(item['sample_id'])
        with np.errstate(invalid='ignore',divide='ignore'),contextlib.redirect_stdout(io.StringIO()):
            full=base.social_with_filter(mod,native,state,para,False)
            no_cross=base.social_with_filter(mod,native,state,para,True)
            no_source=selective_social(base,mod,native,state,para,'source_receives_off',item['source_sign'])
            no_target=selective_social(base,mod,native,state,para,'target_receives_off',item['source_sign'])
        # At fixed geometry, social torque removal partitions exactly by receiver.
        closure=np.max(np.abs((full-no_source)+(full-no_target)-(full-no_cross)))
        if closure>1e-11:raise RuntimeError(f'Receiver decomposition failed {item["sample_id"]}: {closure:g}')
        rows.append(dict(phase=item['phase'],sample_id=item['sample_id'],seed=int(item['seed']),
                         n_cross_directed=n_cross,n_cross_unique=n_unique,n_open_directed=n_open,
                         n_solid_directed=n_solid,mean_abs_cross_torque=float(np.mean(np.abs(full-no_cross))),
                         mean_abs_source_receive_torque=float(np.mean(np.abs(full-no_source))),
                         mean_abs_target_receive_torque=float(np.mean(np.abs(full-no_target))),
                         max_abs_receiver_closure=float(closure)))
    return rows, dict(n_frames=len(rows),cross_directed=sum(r['n_cross_directed'] for r in rows),
                      cross_unique=sum(r['n_cross_unique'] for r in rows),
                      opening_directed=sum(r['n_open_directed'] for r in rows),
                      solid_directed=sum(r['n_solid_directed'] for r in rows),
                      worst_solid_excess=float(worst),solid_examples=examples,
                      all_cross_links_pass_through_gap=(exceed==0))


def init_worker(source,adapter_text):
    # One process owns its independent compiled production timestep.
    import run_bottleneck_commitment_branches as base
    from analyze_approach_current_balance import load_gap20_production,fast_check_voronoi
    mod=load_gap20_production(Path(source))
    mod.CheckVoronoi=lambda vn,loc,pa:fast_check_voronoi(mod,vn,loc)
    para=base.parameters(mod)
    env=dict(vars(mod))
    env['GetVision_intervention']=lambda vn,st,pa,mode:selective_social(base,mod,vn,st,pa,mode,GLOBAL['sign'])
    exec(adapter_text,env)
    GLOBAL.update(base=base,step=env['archived_step'],para=para,sign=1)


def run_task(t):
    b=GLOBAL['base']; step=GLOBAL['step'];para=GLOBAL['para']
    sign=int(t['source_sign']); GLOBAL['sign']=sign
    start=t['state']; nright0=int((start[0]>=0).sum())
    states={arm:start.copy() for arm in ARMS}
    memory={arm:t['dis_last'].copy() for arm in ARMS}
    rng=np.random.default_rng(int(t['noise_seed']))
    noises=0.5*np.sqrt(b.DT)*rng.normal(size=(b.HORIZON_STEPS,b.N))
    inc=np.zeros((3,b.N))
    for z in noises:
        inc[2]=z
        for arm in ARMS:
            out,_,_,memory[arm],_=step(states[arm],para,inc,memory[arm],arm)
            states[arm]=out
    return [dict(phase='post',sample_id=t['sample_id'],seed=int(t['seed']),clone=int(t['clone']),
                 arm=arm,source_sign=sign,nright0=nright0,nrightT=int((state[0]>=0).sum()),
                 additional_net=int(-sign*((state[0]>=0).sum()-nright0)),
                 success4=int(-sign*((state[0]>=0).sum()-nright0)>=4))
            for arm,state in states.items()]


def mean_ci(v,boot=8000,seed=20261113):
    a=np.asarray(v,float)
    if len(a)<2:return dict(mean=float(np.mean(a)),ci95=[None,None],n=len(a))
    rng=np.random.default_rng(seed)
    sample=rng.choice(a,size=(boot,len(a)),replace=True).mean(axis=1)
    return dict(mean=float(np.mean(a)),ci95=[float(x) for x in np.percentile(sample,[2.5,97.5])],n=len(a))


def summarize(rows):
    index={(r['sample_id'],int(r['clone']),r['arm']):r for r in rows}
    keys=sorted({(r['sample_id'],int(r['clone'])) for r in rows})
    by_seed=defaultdict(list)
    discordant={arm:dict(gained=0,lost=0) for arm in ARMS if arm!='full'}
    for sid,clone in keys:
        b=index[(sid,clone,'full')]
        for arm in ARMS[1:]:
            x=index[(sid,clone,arm)]
            if int(x['success4'])>int(b['success4']):discordant[arm]['gained']+=1
            if int(x['success4'])<int(b['success4']):discordant[arm]['lost']+=1
        rec=dict(seed=int(b['seed']),sample_id=sid,clone=clone)
        for arm in ARMS:
            x=index[(sid,clone,arm)]
            rec['net_'+arm]=int(x['additional_net'])
            rec['p_'+arm]=int(x['success4'])
        by_seed[int(b['seed'])].append(rec)
    seeds=sorted(by_seed)
    if len(seeds)!=14 or len(keys)!=112:raise ValueError(f'Expected 14 original seeds and 112 pairs, got {len(seeds)} and {len(keys)}')
    state_rows=[]
    for seed in seeds:
        rs=by_seed[seed]
        d=dict(seed=seed,sample_id=rs[0]['sample_id'],n_clones=len(rs))
        for arm in ARMS:
            d['net_'+arm]=np.mean([r['net_'+arm] for r in rs])
            d['p_'+arm]=np.mean([r['p_'+arm] for r in rs])
            if arm!='full':
                d['effect_net_'+arm]=d['net_full']-d['net_'+arm]
                d['effect_p_'+arm]=d['p_full']-d['p_'+arm]
        d['nonadditivity_net']=d['effect_net_source_receives_off']+d['effect_net_target_receives_off']-d['effect_net_crossgap_off']
        state_rows.append(d)
    terms=['effect_net_crossgap_off','effect_net_source_receives_off','effect_net_target_receives_off',
           'effect_p_crossgap_off','effect_p_source_receives_off','effect_p_target_receives_off','nonadditivity_net']
    estimates={key:mean_ci([r[key] for r in state_rows],seed=20261113+i) for i,key in enumerate(terms)}
    baselines={arm:dict(net=mean_ci([r['net_'+arm] for r in state_rows]),
                        success4=mean_ci([r['p_'+arm] for r in state_rows])) for arm in ARMS}
    return state_rows,dict(n_seeds=len(seeds),n_paired_trials=len(keys),baselines=baselines,
                           effects=estimates,discordant_threshold4=discordant,
                           positive_effect_counts={key:dict(positive=int(sum(r[key]>0 for r in state_rows)),negative=int(sum(r[key]<0 for r in state_rows))) for key in terms})


def previous_pilot_parity(newrows,prior_rows):
    orig={(r['sample_id'],int(r['clone']),r['arm']):r for r in prior_rows if r['phase']=='post'}
    count=0;mismatches=[]
    for r in newrows:
        if r['arm'] not in ('full','crossgap_off'):continue
        key=(r['sample_id'],int(r['clone']),r['arm'])
        old=orig.get(key)
        if old is None or int(old['nrightT'])!=int(r['nrightT']):mismatches.append(key)
        else:count+=1
    return dict(n_compared=count,n_mismatches=len(mismatches),sample_mismatches=mismatches[:8],passed=(count==224 and not mismatches))


def report(meta):
    a=meta['structural'],meta.get('results');s,res=a
    lines=['# Directional social-feedback test after first crossing','',
           'Fixed 14 post-first-crossing microstates, eight paired noise streams each. Production dynamics and positions are unchanged except for receiver-side cross-divider social torque ablations.','',
           '## Geometry and reproducibility gates',
           f'- Previously audited production SHA/adapter parity: {meta["hash_gate_passed"]}.',
           f'- Admissible visual links crossing the divider: {s["cross_directed"]} directed, of which {s["opening_directed"]} pass through the opening and {s["solid_directed"]} cross solid wall.',
           f'- All cross-divider visual links cross the aperture: {s["all_cross_links_pass_through_gap"]}.',
           f'- Original pilot arm replay: {meta.get("pilot_parity",{}).get("n_compared","not run")} paired arm endpoints matched, {meta.get("pilot_parity",{}).get("n_mismatches","not run")} mismatches.','']
    if res:
        lines+=['## Causal effect estimates','',
                'Each effect is **full minus ablated**, so positive net means the removed coupling supported additional transfer. 95% intervals are bootstrap over the 14 original simulation seeds.','',
                '| Removed social input received by | Mean reduction of additional net fish | 95% bootstrap interval | Threshold ≥4 outcome changes |',
                '|---|---:|---|---|']
        d=res['effects']; disc=res['discordant_threshold4']
        for arm,label in [('source_receives_off','Source-chamber fish'),('target_receives_off','Target-chamber fish'),('crossgap_off','Both sides')]:
            x=d['effect_net_'+arm];p=disc[arm]
            lines.append(f'| {label} | {x["mean"]:+.3f} | [{x["ci95"][0]:+.3f}, {x["ci95"][1]:+.3f}] | {p["gained"]} gains / {p["lost"]} losses |')
        ni=d['nonadditivity_net'];lines +=['',f'- Nonadditivity of final transfer effects: {ni["mean"]:+.3f} (CI {ni["ci95"]}). This is dynamical interaction, not an algebraic error.',
                '', '## Interpretation','',
                'The two receiver-targeted ablations determine whether cross-divider social coupling changes the motion of fish still in the source chamber, of fish already across the gap, or both. Neither a positive net count nor a threshold-four event by itself establishes whole-school relocation.',
                'Selectively removing terms from GetVision preserves its original full-neighborhood denominator. Changes in final fish counts remain nonlinear in the perturbation, so directional effects need not sum to the both-sides ablation.',
                'Source/target designations here refer to **current spatial chamber** relative to the original majority side at each time step; fish identities can cross.',
                '']
    else:lines+=['## Stopping conclusion','','No directional simulation was run because the geometry/source audit failed. Identify the structural mismatch before interpreting the previous cross-divider result.','']
    return '\n'.join(lines)


def self_test():
    # Link (0,1) intersects y=0 through gap; (2,3) crosses solid y=2.
    st=np.array([[-1,1,-1,1],[0,0,2,2],[0,0,0,0]],float)
    vn=np.zeros((4,4),bool);vn[0,1]=vn[1,0]=vn[2,3]=vn[3,2]=True
    g=crossing_link_geometry(st,vn,1.)
    assert g['opening'].sum()==2 and g['solid'].sum()==2
    class Stub:
        @staticmethod
        def social_with_filter(mod, vn, st, pa, off):
            return np.array([1,2,3,4] if off else [6,7,8,9],float)
    for sign in (-1,1):
        f=selective_social(Stub,None,vn,st,None,'full',sign)
        b=selective_social(Stub,None,vn,st,None,'crossgap_off',sign)
        s=selective_social(Stub,None,vn,st,None,'source_receives_off',sign)
        t=selective_social(Stub,None,vn,st,None,'target_receives_off',sign)
        assert np.allclose((f-s)+(f-t),f-b)
        assert np.count_nonzero(f-s)==2 and np.count_nonzero(f-t)==2
    # Paired seeded mock endpoint example: 14 seeds, 8 clones, four arms.
    rows=[]
    for seed in range(14):
        for clone in range(8):
            for arm in ARMS:
                v={'full':9,'crossgap_off':7,'source_receives_off':8,'target_receives_off':7}[arm]
                rows.append(dict(phase='post',sample_id=f's{seed}',seed=seed,clone=clone,
                                 arm=arm,additional_net=v,success4=int(v>=4)))
    per,summary=summarize(rows)
    assert len(per)==14 and summary['effects']['effect_net_crossgap_off']['mean']==2
    assert summary['effects']['effect_net_source_receives_off']['mean']==1
    print('PASS: crossing geometry, receiver-specific torque mask, paired-seed summary')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--source',type=Path,default=None)
    parser.add_argument('--raw-dir',type=Path,default=None)
    parser.add_argument('--pilot-dir',type=Path,default=PILOT_DEFAULT)
    parser.add_argument('--output',type=Path,default=OUT_DEFAULT)
    parser.add_argument('--workers',type=int,default=min(6,os.cpu_count() or 1))
    parser.add_argument('--stage',choices=['audit','all'],default='all')
    args=parser.parse_args()
    if args.self_test:return self_test()
    import run_bottleneck_commitment_branches as base
    from analyze_approach_current_balance import load_gap20_production
    source=(args.source or base.SOURCE).resolve();raw_dir=(args.raw_dir or base.RAW_DEFAULT).resolve()
    if not source.is_file():raise FileNotFoundError(source)
    original=json.loads((args.pilot_dir/'simulator_audit.json').read_text())
    mod=load_gap20_production(source);para=base.parameters(mod)
    step,adapter_text,adapter_sha=base.extract_step(mod,source)
    source_sha=hashlib.sha256(source.read_bytes()).hexdigest()
    hash_pass=(source_sha==original['source_sha256'] and adapter_sha==original['adapter_sha256'])
    if not hash_pass:raise RuntimeError('Production source or restart adapter changed: stop before running')
    before=base.choose_pre();after=base.choose_post()
    origins=base.setup_origins(raw_dir,step,para,before,after)
    audit_rows,audit_sum=structural_audit(base,mod,para,origins)
    args.output.mkdir(parents=True,exist_ok=True)
    csv_write(args.output/'aperture_visibility_audit.csv',audit_rows)
    meta=dict(structural=audit_sum,hash_gate_passed=True,pilot_parity={},results=None,
              source_sha256=source_sha,adapter_sha256=adapter_sha,
              analysis='same 14 post-entry states and same 8 noise streams as previous pilot')
    (args.output/'audit.json').write_text(json.dumps(meta,indent=2)+'\n')
    if not audit_sum['all_cross_links_pass_through_gap']:
        (args.output/'REPORT.md').write_text(report(meta)+'\n')
        raise RuntimeError('Some eligible crossdivider links cross solid wall; audit physical visibility before intervention')
    if args.stage=='audit':
        (args.output/'REPORT.md').write_text(report(meta)+'\n')
        print(json.dumps(audit_sum,indent=2));return
    post=[x for x in origins if x['phase']=='post']
    tasks=[]
    for item in post:
        for clone in range(base.N_CLONES):
            noise_seed=202611130000+2*100000+int(item['seed'])*100+clone
            tasks.append(dict(item,clone=clone,noise_seed=noise_seed))
    rows=[];partial=args.output/'branches_partial.csv'
    prior={}
    if partial.exists():
        partial_rows=csv_read(partial)
        prior=defaultdict(list)
        for r in partial_rows:prior[(r['sample_id'],int(r['clone']))].append(r)
        rows=[r for key,group in prior.items() if len(group)==4 and {g['arm'] for g in group}==set(ARMS) for r in group]
    done={(r['sample_id'],int(r['clone'])) for r in rows}
    pending=[t for t in tasks if (t['sample_id'],t['clone']) not in done]
    if pending:
        context=mp.get_context('spawn')
        with context.Pool(min(args.workers,len(pending)),initializer=init_worker,
                          initargs=(str(source),adapter_text)) as pool:
            for n,part in enumerate(pool.imap_unordered(run_task,pending),1):
                rows.extend(part)
                if n%5==0 or n==len(pending):
                    csv_write(partial,rows)
                    print(f'completed {n}/{len(pending)} paired original-state restarts',flush=True)
    rows.sort(key=lambda r:(r['sample_id'],int(r['clone']),r['arm']))
    old=csv_read(args.pilot_dir/'branches.csv')
    parity=previous_pilot_parity(rows,old)
    meta['pilot_parity']=parity
    if not parity['passed']:
        (args.output/'audit.json').write_text(json.dumps(meta,indent=2)+'\n')
        raise RuntimeError(f'Original pilot arm parity failed: {parity}')
    csv_write(args.output/'branches.csv',rows)
    if partial.exists():partial.unlink()
    per,result=summarize(rows)
    csv_write(args.output/'per_microstate.csv',per)
    meta['results']=result
    (args.output/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
    (args.output/'REPORT.md').write_text(report(meta)+'\n')
    print(json.dumps({k:v for k,v in result['effects'].items()},indent=2))

if __name__=='__main__':main()
