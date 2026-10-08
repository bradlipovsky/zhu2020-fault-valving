"""Compare one and eight OpenMP threads through one saved seismic second.

The deposited initial checkpoint is from the T=1e10 s calculation. The two
processes run concurrently on the same host; timings describe this benchmark.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
import struct
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
state=ROOT/'results/threading_serial/phase1_state.bin'
with state.open('rb') as stream:
    stream.seek(-8,2);start=struct.unpack('d',stream.read(8))[0]
n=(state.stat().st_size-8)//32-1

def run(threads):
    name='threading_serial' if threads==1 else 'threading_parallel'
    path=ROOT/'results'/name
    command=[str(ROOT/'build/valving'),'--n',str(n),'--T','1e10',
        '--start-year',repr(start/31536000),'--years',repr((start+1)/31536000),
        '--state',str(state),'--profile','inputs/T1e10_restart.csv','--output',str(path)]
    before=time.monotonic()
    with path.with_suffix('.log').open('w') as log:
        subprocess.run(command,cwd=str(ROOT),env=dict(os.environ,OMP_NUM_THREADS=str(threads)),
                       stdout=log,stderr=log,check=True)
    elapsed=time.monotonic()-before
    metadata=json.loads((path/'metadata.json').read_text())
    assert metadata['worker_threads']==threads,'Build with OpenMP to run this comparison'
    return elapsed

with ThreadPoolExecutor(max_workers=2) as pool:
    serial,parallel=list(pool.map(run,[1,8]))
result={'start_year':start/31536000,'end_year':(start+1)/31536000,
        'serial_wall_s':serial,'parallel_wall_s':parallel,'speedup':serial/parallel,
        'threads':8,'bitwise_identical':True,
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
        'source_sha256':hashlib.sha256((ROOT/'src/valving.cpp').read_bytes()).hexdigest(),
        'initial_state_sha256':hashlib.sha256(state.read_bytes()).hexdigest()}
for name in ['fields.bin','final_state.bin','history.csv']:
    a=(ROOT/'results/threading_serial'/name).read_bytes()
    b=(ROOT/'results/threading_parallel'/name).read_bytes()
    assert a==b,'Thread count changed '+name
    result[name+'_sha256']=hashlib.sha256(a).hexdigest()
(ROOT/'results/threading_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
