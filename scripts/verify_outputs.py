#!/usr/bin/env python3
"""Check complete raw runs, including output outside the selected figure windows."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import numpy as np
from common import YEAR, MAIN_CASES, configuration, fields, physical_depths, save_json
from reproduce import VALIDATION_CASES


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def verify_case(case):
    folder=Path('data',case)
    done=json.loads((folder/'completed.json').read_text())
    provenance=json.loads((folder/'provenance.json').read_text())
    cfg=configuration('configs/'+case+'.cfg')
    assert done['completed'] and provenance['exit_code']==0,case+' incomplete run'
    for key,path in [('source','src/model.cpp'),('executable','build/valving'),
                     ('config','configs/'+case+'.cfg')]:
        assert provenance[key+'_sha256']==sha(path),case+' '+key+' changed'
    assert (folder/'config.cfg').read_bytes()==Path('configs',case+'.cfg').read_bytes(),case+' archived configuration'
    assert np.isclose(done['time_s'],cfg['years']*YEAR,rtol=1e-14,atol=0),case+' end time'
    assert done['n']==cfg['n'] and done['dynamic_nodes']==len(physical_depths(cfg)),case+' mesh'

    with (folder/'history.csv').open() as stream:
        names=stream.readline().strip().split(',')
        values=np.loadtxt(stream,delimiter=',',ndmin=2)
    assert values.shape==(done['accepted_steps']+1,len(names)),case+' history length'
    assert np.all(np.isfinite(values)),case+' nonfinite history'
    h={name:values[:,i] for i,name in enumerate(names)}
    assert np.array_equal(h['step'],np.arange(len(values))),case+' missing history steps'
    assert h['time_s'][0]==0 and h['dt_s'][0]==0,case+' initial history'
    assert np.all(np.diff(h['time_s'])>0) and np.all(h['dt_s'][1:]>0),case+' history ordering'
    assert np.isclose(h['time_s'][-1],done['time_s'],rtol=1e-14,atol=0),case+' final history'
    assert np.allclose(np.diff(h['time_s']),h['dt_s'][1:],rtol=0,
                       atol=4*np.finfo(float).eps*h['time_s'][-1]),case+' step durations'
    assert np.all((h['error']>=0)&(h['error']<=1)),case+' accepted error'
    assert np.all(h['min_effective_pa']>0),case+' nonpositive effective stress'
    for name in names:
        if name.startswith('effective_'):assert np.all(h[name]>0),case+' '+name
        if name.startswith('permeability_'):
            assert np.all(h[name]>=cfg['kmin']*(1-1e-12)) and np.all(h[name]<=cfg['kmax']*(1+1e-12)),case+' '+name
    for key in ['Lb','hstar']:
        minimum=float(h['min_'+key+'_cells'].min())
        assert minimum>0 and np.isclose(minimum,done['minimum_'+key+'_cells'],rtol=1e-12,atol=0),case+' resolution minimum'

    z,records=fields(case)
    assert len(records)>1 and np.all(np.isfinite(z)) and np.all(np.diff(z)>0),case+' field coordinates'
    assert np.all(np.isin(z,physical_depths(cfg))),case+' field depths outside mesh'
    expected_size=16+8*len(z)+len(records)*records.dtype.itemsize
    assert (folder/'fields.bin').stat().st_size==expected_size,case+' truncated field record'
    steps=records['step']
    assert steps[0]==0 and steps[-1]==done['accepted_steps'] and np.all(np.diff(steps)>0),case+' saved steps'
    assert np.array_equal(records['time'],h['time_s'][steps]),case+' saved times differ from history'
    # Float32 snapshots need a rounding allowance; the accepted history is float64.
    rounding=8*np.finfo(np.float32).eps
    for first in range(0,len(records),256):
        chunk=records['fields'][first:first+256]
        assert np.all(np.isfinite(chunk)),case+' nonfinite saved fields'
        assert np.all(chunk[:,2]>0),case+' nonpositive saved effective stress'
        for field in [3,4]:
            assert np.all(chunk[:,field]>=cfg['kmin']*(1-rounding)) and np.all(chunk[:,field]<=cfg['kmax']*(1+rounding)),case+' saved permeability bounds'
        assert np.all(chunk[:,3]<=chunk[:,4]*(1+rounding)),case+' permeability exceeds reference permeability'
        if cfg['fixed_pressure']:
            assert np.all(chunk[:,2]==records['fields'][0,2]),case+' fixed pressure changed'

    result=dict(case=case,checks_passed=True,command=['python3','scripts/verify_outputs.py',case],
        scope='All accepted-step history rows and all saved field profiles; does not establish convergence or agreement with the paper.',
        source_sha256=provenance['source_sha256'],executable_sha256=provenance['executable_sha256'],
        config_sha256=provenance['config_sha256'],verifier_sha256=sha('scripts/verify_outputs.py'),
        years=done['years'],accepted_steps=done['accepted_steps'],saved_profiles=len(records),
        maximum_accepted_error=float(h['error'].max()),minimum_effective_pa=float(h['min_effective_pa'].min()),
        minimum_Lb_cells=float(h['min_Lb_cells'].min()),minimum_hstar_cells=float(h['min_hstar_cells'].min()),
        raw_files={name:dict(bytes=(folder/name).stat().st_size,sha256=sha(folder/name))
                   for name in ['history.csv','fields.bin']},
        checks=['Successful completed run, source/executable/configuration hashes, archived configuration, end time and mesh.',
                'Finite complete history, consecutive accepted steps, increasing time, consistent step durations and accepted errors at most one.',
                'Positive effective stress, permeability bounds and whole-history resolution minima.',
                'Complete binary records, mesh coordinates, finite fields and exact saved time/step agreement with history.',
                'Positive saved stress, constitutive permeability bounds and permeability no greater than reference permeability.',
                'Effective stress constant in every fixed-pressure saved profile when that option is enabled.'])
    save_json(folder/'output_verification.json',result)
    print('VERIFIED',case,done['accepted_steps'],'accepted steps;',len(records),'saved profiles',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases',nargs='*')
    args=parser.parse_args()
    os.chdir(str(Path(__file__).resolve().parent.parent))
    for case in args.cases or MAIN_CASES+VALIDATION_CASES:verify_case(case)
