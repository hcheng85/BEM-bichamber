#!/usr/bin/env python3
"""Recompute precursor chronology on one common event cohort across all leads."""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np
from scipy.stats import rankdata

LEADS = (5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.2, 0.1)

def auc(y, x):
    y=np.asarray(y,bool); x=np.asarray(x,float); n1=y.sum(); n0=(~y).sum()
    return float((rankdata(x)[y].sum()-n1*(n1+1)/2)/(n1*n0))

def write(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    rows=list(csv.DictReader(a.input.open())); by={}
    for r in rows: by.setdefault((r['event_seed'],r['event_id']),set()).add(float(r['lead']))
    common={key for key,val in by.items() if set(LEADS)<=val}
    rows=[r for r in rows if (r['event_seed'],r['event_id']) in common]
    rng=np.random.default_rng(20261110); out=[]
    for lead in LEADS:
        sub=[r for r in rows if float(r['lead'])==lead]; seeds=sorted({int(r['event_seed']) for r in sub})
        for metric in ('a','j_source','j_approach'):
            byseed={s:np.mean([float(r[f'delta_{metric}']) for r in sub if int(r['event_seed'])==s]) for s in seeds}
            means=np.array(list(byseed.values())); draws=rng.choice(means,(4000,len(means)),replace=True).mean(1); lo,hi=np.quantile(draws,[.025,.975])
            ev=np.array([float(r[f'event_{metric}']) for r in sub]); ct=np.array([float(r[f'control_{metric}']) for r in sub]); estimate_auc=auc(np.r_[np.ones(len(ev)),np.zeros(len(ct))],np.r_[ev,ct])
            auc_draw=[]
            for _ in range(3000):
                chosen=rng.choice(seeds,len(seeds),replace=True); selected=[]
                for s in chosen: selected.extend(r for r in sub if int(r['event_seed'])==s)
                ee=np.array([float(r[f'event_{metric}']) for r in selected]); cc=np.array([float(r[f'control_{metric}']) for r in selected]); auc_draw.append(auc(np.r_[np.ones(len(ee)),np.zeros(len(cc))],np.r_[ee,cc]))
            alo,ahi=np.quantile(auc_draw,[.025,.975])
            out.append({'lead':lead,'metric':metric,'mean_paired_difference':means.mean(),'difference_ci95_low':lo,'difference_ci95_high':hi,'auc':estimate_auc,'auc_ci95_low':alo,'auc_ci95_high':ahi,'n_events':len(sub),'n_event_seeds':len(seeds)})
    write(a.output,out)

if __name__=='__main__': main()
