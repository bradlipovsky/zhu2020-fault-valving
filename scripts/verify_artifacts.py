#!/usr/bin/env python3
"""Audit completeness, numerical panel provenance, and the compiled deliverable."""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
from common import MAIN_CASES, save_json, configuration
from inventory import PANELS, GROUPS
from reproduce import VALIDATION_CASES

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def main():
    source=sha('src/model.cpp');checked=[]
    assert len(PANELS)==49 and len(GROUPS)==12
    for case in MAIN_CASES+VALIDATION_CASES:
        folder=Path('data',case)
        done=json.loads((folder/'completed.json').read_text());provenance=json.loads((folder/'provenance.json').read_text())
        assert done['completed'] and done['accepted_steps']>0,case
        assert provenance['source_sha256']==source,case+' source changed'
        assert provenance['config_sha256']==sha('configs/'+case+'.cfg'),case+' configuration changed'
        assert provenance['exit_code']==0,case+' failed'
        assert abs(done['time_s']/31557600-done['years'])<1e-10
        configured=configuration('configs/'+case+'.cfg')
        assert abs(done['years']-configured['years'])<1e-9,case+' incorrect end time'
        assert done['n']==int(configured['n']),case+' incorrect mesh'
        assert 8<=done['dynamic_nodes']<=done['n'],case+' invalid physical mesh'
        assert done['minimum_Lb_cells']>0 and done['minimum_hstar_cells']>0,case+' invalid resolution diagnostic'
        checked.extend([folder/'completed.json',folder/'analysis.json',folder/'provenance.json'])
    for name,spec in PANELS.items():
        path=Path(spec['data_file']);meta=json.loads(path.with_suffix('.json').read_text())
        assert meta['model_source_sha256']==source,name
        for script,digest in meta['pipeline_source_sha256'].items():assert sha(script)==digest,name+' script changed: '+script
        assert meta['data_command']==spec['data_command'] and meta['plot_command']==spec['plot_command']
        assert meta['status'] in ['independently reproduced','partial','not reproduced']
        with np.load(str(path),allow_pickle=False) as d:
            for key in d.files:
                values=d[key]
                if values.dtype.kind in 'f':assert np.all(np.isfinite(values)),name+' '+key
            if 'values' in d:
                assert d['values'].shape==(len(d['time_s']),len(d['depth_m'])),name
                assert len(d['time_s'])>1 and np.all(np.diff(d['time_s'])>=0),name
                assert d['step'].shape==d['time_s'].shape and np.all(np.diff(d['step'])>=0),name+' inconsistent step coordinate'
                assert np.all(np.diff(d['depth_m'])>0),name
                field=int(d['field'])
                if field==2:assert d['values'].min()>0,name+' nonpositive effective stress'
                if field==3:
                    assert d['values'].min()>1e-19*(1-1e-6) and d['values'].max()<1e-15*(1+1e-6),name
        checked.extend([path,path.with_suffix('.json')])
        for extension in ['pdf','png']:
            f=Path('figures/panels',name+'.'+extension);assert f.stat().st_size>1000,str(f);checked.append(f)
    for group in GROUPS:
        for extension in ['pdf','png']:
            f=Path('figures',group+'.'+extension);assert f.stat().st_size>1000;checked.append(f)
    pdf=Path('report/reproduction.pdf');assert pdf.read_bytes()[:4]==b'%PDF'
    text=subprocess.check_output(['pdftotext',str(pdf),'-']).decode()
    for phrase in ['Independent calculations','Panel provenance','Quantitative comparisons','Supplementary']:
        assert phrase in text,'Missing report content: '+phrase
    for name in PANELS:assert name in text,'Panel missing from PDF provenance: '+name
    for forbidden in ['TODO','PLACEHOLDER','awaiting simulation']:
        assert forbidden not in text,'Unfinished report marker: '+forbidden
    log=Path('report/reproduction.log').read_text(errors='replace')
    for error in ['LaTeX Error','There were undefined references','Citation `']:
        assert error not in log,error
    checked.extend(Path('report').glob('*.tex'));checked.append(pdf)
    checked.extend([Path('data/verification.txt'),Path('data/coupled_verification.txt'),Path('data/steady_verification.txt'),Path('data/graded_verification.txt'),
        Path('data/quantitative_comparisons.json'),Path('data/validation_summary.json')])
    checked.extend(Path('data/laws').glob('*.csv'))
    checked.extend(Path('configs',case+'.cfg') for case in MAIN_CASES+VALIDATION_CASES)
    checked.extend(Path('scripts').glob('*.py'));checked.extend(Path('src').glob('*.cpp'))
    bundle=json.loads(Path('data/report_bundle_verification.json').read_text())
    assert bundle['isolated_compile'] and bundle['sha256']==sha(bundle['archive'])
    checked.extend([Path(bundle['archive']),Path('data/report_bundle_verification.json'),Path('data/report_bundle_build.log')])
    manifest={str(p):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(checked))}
    save_json('data/artifact_manifest.json',dict(model_source_sha256=source,panels=49,figure_groups=12,
        complete_cases=len(MAIN_CASES+VALIDATION_CASES),files=manifest))
    print('VERIFIED: 12 completed numerical cases, 49 panel datasets, 49 individual panel figures, 12 figure groups, compiled report.')

if __name__=='__main__':main()
