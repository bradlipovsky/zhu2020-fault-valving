#!/usr/bin/env python3
"""Check that the raw-output audit rejects damaged copies of a software fixture."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import struct
import subprocess
import tempfile

root=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('fixture',type=Path,help='Directory made by scripts/smoke_pipeline.py')
args=parser.parse_args()
fixture=args.fixture.resolve()
assert (fixture/'NOT_SCIENTIFIC_RESULTS.txt').exists(),'Use a software fixture, not production output.'
outcomes=[]
with tempfile.TemporaryDirectory(prefix='output-audit-',dir=str(root/'.tmp')) as temporary:
    work=Path(temporary)
    for name in ['scripts','src','build','configs','data/reference']:(work/name).mkdir(parents=True,exist_ok=True)
    for name in ['verify_outputs.py','common.py','reproduce.py']:
        shutil.copy(str(root/'scripts'/name),str(work/'scripts'/name))
    for name in ['src/model.cpp','build/valving','configs/reference.cfg']:
        shutil.copy(str(fixture/name),str(work/name))
    for name in ['config.cfg','completed.json','provenance.json','history.csv','fields.bin']:
        shutil.copy(str(fixture/'data/reference'/name),str(work/'data/reference'/name))
    folder=work/'data/reference'
    original={name:(folder/name).read_bytes() for name in ['history.csv','fields.bin']}
    nz=struct.unpack_from('<Q',original['fields.bin'],8)[0]
    second_record=16+8*nz+16+24*nz
    command=['python3','scripts/verify_outputs.py','reference']
    def run(label,expected):
        result=subprocess.run(command,cwd=str(work),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              env=dict(os.environ,TMPDIR=str(root/'.tmp')))
        output=result.stdout.decode()
        assert (result.returncode==0)==(expected is None),(label,output)
        if expected:assert expected in output,(label,output)
        outcomes.append(dict(check=label,expected_failure=expected,exit_code=result.returncode,output=output))
        for name,content in original.items():(folder/name).write_bytes(content)
    run('Intact independently calculated short fixture',None)
    (folder/'fields.bin').write_bytes(original['fields.bin'][:-1])
    run('Truncated final binary record','truncated field record')
    rows=original['history.csv'].splitlines(keepends=True)
    (folder/'history.csv').write_bytes(b''.join(rows[:2]+rows[3:]))
    run('Missing accepted history row','history length')
    for label,column,value,expected in [
        ('Nonfinite history value',5,b'nan','nonfinite history'),
        ('Nonpositive effective stress',5,b'-1','nonpositive effective stress'),
        ('Accepted error above one',7,b'2','accepted error')]:
        damaged=list(rows);values=damaged[2].rstrip(b'\n').split(b',');values[column]=value
        damaged[2]=b','.join(values)+b'\n';(folder/'history.csv').write_bytes(b''.join(damaged))
        run(label,expected)
    for label,offset,fmt,value,expected in [
        ('Inconsistent saved time',second_record,'<d',123.0,'saved times differ from history'),
        ('Permeability outside constitutive bounds',second_record+16+3*nz*4,'<f',2e-15,'saved permeability bounds'),
        ('Changing fixed-pressure profile',second_record+16+2*nz*4,'<f',1.0,'fixed pressure changed')]:
        damaged=bytearray(original['fields.bin']);struct.pack_into(fmt,damaged,offset,value)
        (folder/'fields.bin').write_bytes(damaged);run(label,expected)

record=dict(scope='Software audit validation only; no damaged or fixture output is used as a scientific result.',
    command=['python3','tests/output_audit.py',str(args.fixture)],fixture=str(fixture),
    verifier_sha256=hashlib.sha256((root/'scripts/verify_outputs.py').read_bytes()).hexdigest(),
    test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks_passed=True,outcomes=outcomes)
(root/'data/output_audit_verification.json').write_text(json.dumps(record,indent=2)+'\n')
print('VERIFIED intact fixture and',len(outcomes)-1,'deliberate output corruptions')
