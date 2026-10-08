#!/usr/bin/env python3
"""Exercise the full artifact pipeline in isolation; these are never scientific results."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import time
from common import save_json

def main():
    root=Path(__file__).resolve().parent.parent
    test=root/'.tmp'/('pipeline-smoke-'+str(int(time.time()*1e9)));test.mkdir(parents=True)
    for name in ['src','tests','scripts','configs','documentation']:
        shutil.copytree(str(root/name),str(test/name),ignore=shutil.ignore_patterns('__pycache__'))
    (test/'report').mkdir()
    shutil.copy(str(root/'report/reproduction.tex'),str(test/'report/reproduction.tex'))
    for name in ['CMakeLists.txt','requirements.txt']:shutil.copy(str(root/name),str(test/name))
    (test/'.deps').symlink_to(root/'.deps',target_is_directory=True)
    # These intentionally coarse, short runs test software integration only.
    overrides=dict(n=512,threads=1,years=.002,output_years=.0001,output_spacing=1000)
    for path in (test/'configs').glob('*.cfg'):
        lines=[]
        for line in path.read_text().splitlines():
            key=line.split('=')[0].strip()
            lines.append('{} = {}'.format(key,overrides[key]) if key in overrides else line)
        path.write_text('\n'.join(lines)+'\n')
    (test/'NOT_SCIENTIFIC_RESULTS.txt').write_text('Software smoke test only. Do not publish these generated figures as reproduction results.\n')
    command=['python3','scripts/reproduce.py','--jobs','3'];started=time.time()
    with (root/'data/pipeline_smoke.log').open('w') as log:
        result=subprocess.run(command,cwd=str(test),stdout=log,stderr=subprocess.STDOUT)
    record=dict(purpose='software integration only; no smoke figure used in scientific report',
        isolated_directory=str(test),configuration_overrides=overrides,command=command,
        source_sha256=hashlib.sha256((root/'src/model.cpp').read_bytes()).hexdigest(),
        exit_code=result.returncode,wall_s=time.time()-started)
    save_json(root/'data/pipeline_smoke_verification.json',record)
    print(json.dumps(record,indent=2))
    if result.returncode:raise SystemExit(result.returncode)

if __name__=='__main__':main()
