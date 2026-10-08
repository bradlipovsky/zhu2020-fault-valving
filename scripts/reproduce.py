#!/usr/bin/env python3
"""Run the independent calculation. No reference document or old result is read."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import time
from common import MAIN_CASES, save_json

VALIDATION_CASES=['baseline_coarse','baseline_fine','baseline_tight','initialization','normal_stress']

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--simulate-only',action='store_true')
    parser.add_argument('--postprocess-only',action='store_true')
    parser.add_argument('--cases',nargs='+')
    parser.add_argument('--jobs',type=int,default=12)
    parser.add_argument('--status-file',default='data/run_status.json',help='Separate progress file for concurrent case batches.')
    parser.add_argument('--fresh',action='store_true',help='Move existing generated cases aside before regenerating.')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent.parent;os.chdir(str(root))
    Path('data').mkdir(exist_ok=True)
    cases=args.cases or MAIN_CASES+VALIDATION_CASES
    if not args.postprocess_only:
        subprocess.run(['bash','scripts/bootstrap.sh'],check=True)
        with open('data/verification.txt','w') as f:subprocess.run(['build/valving','--test'],stdout=f,check=True)
        with open('data/coupled_verification.txt','w') as f:subprocess.run(['build/coupled_test'],stdout=f,check=True)
        with open('data/steady_verification.txt','w') as f:subprocess.run(['build/steady_test'],stdout=f,check=True)
        with open('data/graded_verification.txt','w') as f:subprocess.run(['build/graded_test'],stdout=f,check=True)
        with open('data/analysis_verification.txt','w') as f:subprocess.run(['python3','tests/analysis_boundaries.py'],stdout=f,check=True)
        subprocess.run(['build/valving','--laws','configs/baseline.cfg','data/laws'],check=True)
        source_hash=digest('src/model.cpp'); executable_hash=digest('build/valving')
        pending=[]
        for case in cases:
            folder=Path('data',case); config_hash=digest(Path('configs',case+'.cfg'))
            if folder.exists() and args.fresh:
                backup=Path('.tmp','previous-{}-{}'.format(case,int(time.time()*1e9)))
                folder.rename(backup)
            if (folder/'completed.json').exists():
                prior=json.loads((folder/'provenance.json').read_text())
                matching=prior['source_sha256']==source_hash and prior['config_sha256']==config_hash
                if matching and (folder/'fields.bin').exists():
                    print('Verified cached independent run:',case,flush=True);continue
                if matching:
                    folder.rename(Path('.tmp','metadata-{}-{}'.format(case,int(time.time()*1e9))))
                else:raise RuntimeError('Source/config changed for '+case+'; use --fresh')
            if (folder/'fields.bin').exists():
                raise RuntimeError('Incomplete existing data for '+case+'; inspect its process before restarting')
            pending.append(case)
        active={}; outcomes=[]; started=time.time()
        cpus=sorted(os.sched_getaffinity(0)); slots=max(1,min(args.jobs,len(cpus)//4))
        available=list(range(slots))
        while pending or active:
            while pending and available:
                case=pending.pop(0);slot=available.pop(0);folder=Path('data',case);folder.mkdir(exist_ok=True)
                log=(folder/'run.log').open('w')
                env=dict(os.environ,OMP_PROC_BIND='true',OMP_PLACES='cores',
                         TMPDIR=str(root/'.tmp'),MPLBACKEND='Agg')
                assigned=cpus[slot*4:(slot+1)*4]
                def affinity(assigned=assigned):os.sched_setaffinity(0,set(assigned))
                command=['build/valving','configs/'+case+'.cfg',str(folder)]
                proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,preexec_fn=affinity)
                meta=dict(case=case,command=command,pid=proc.pid,start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                    source_sha256=source_hash,executable_sha256=executable_hash,
                    config_sha256=digest('configs/'+case+'.cfg'),cpu_affinity=assigned,
                    scientific_input='published equations and explicitly configured assumptions; no imported simulation state')
                save_json(folder/'provenance.json',meta)
                active[case]=(proc,log,slot,meta);print('START',case,'pid',proc.pid,flush=True)
            running=[]
            for case,(proc,log,slot,meta) in list(active.items()):
                code=proc.poll()
                if code is not None:
                    log.close();meta.update(exit_code=code,end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
                    save_json(Path('data',case,'provenance.json'),meta)
                    outcomes.append(dict(case=case,exit_code=code));available.append(slot);del active[case]
                    print('EXIT',case,code,flush=True)
                else:running.append(dict(case=case,pid=proc.pid))
            save_json(args.status_file,dict(runner_pid=os.getpid(),running=running,pending=pending,
                outcomes=outcomes,elapsed_s=time.time()-started,source_sha256=source_hash))
            if active:time.sleep(10)
        if any(o['exit_code'] for o in outcomes):raise RuntimeError('A model run failed; inspect its run.log')
    if not args.simulate_only:
        subprocess.run(['python3','scripts/analyze.py']+cases,check=True)
        subprocess.run(['python3','scripts/prepare_panels.py'],check=True)
        subprocess.run(['python3','scripts/plot.py','all'],check=True)
        subprocess.run(['python3','scripts/report_tables.py'],check=True)
        subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','reproduction.tex'],cwd='report',check=True)
        subprocess.run(['python3','scripts/package_report.py','--verify'],check=True)
        subprocess.run(['python3','scripts/verify_artifacts.py'],check=True)

if __name__=='__main__':main()
