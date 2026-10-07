"""Compare checkpoint continuation with uninterrupted integration."""
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

binary=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix='restart-test-',dir=str(binary.parent)) as directory:
    root=Path(directory)
    continuous=root/'continuous';split=root/'split'
    def run(out,end,extra=()):
        subprocess.run([str(binary),'--n','512','--years',str(end),'--T','1e8','--perturb','1e-4',
            '--output',str(out),*extra],check=True,stdout=subprocess.DEVNULL)
    run(continuous,.03)
    run(split,.01)
    run(split,.03,['--state',str(split/'final_state.bin'),'--start-year','.01','--append'])
    def read(path):
        b=path.read_bytes();return struct.unpack(str(len(b)//8)+'d',b)
    x=read(continuous/'final_state.bin');y=read(split/'final_state.bin')
    scales=[.002,.01,1,30e6]
    error=max(abs(a-b)/scales[i//513] for i,(a,b) in enumerate(zip(x,y)))
    assert error<5e-4,error
    meta=json.loads((split/'metadata.json').read_text())
    width=8*(1+7*meta['nz']);b=(split/'fields.bin').read_bytes()
    times=[struct.unpack('d',b[j:j+8])[0] for j in range(0,len(b),width)]
    assert all(t1<t2 for t1,t2 in zip(times,times[1:])), 'Duplicated restart sample'
    assert abs(times[-1]/31536000-.03)<1e-12
    digest=hashlib.sha256(b).hexdigest()
    bad=subprocess.run([str(binary),'--n','512','--years','.04','--T','1e8',
        '--output',str(split),'--state',str(split/'final_state.bin'),
        '--start-year','.02','--append'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    assert bad.returncode!=0
    assert hashlib.sha256((split/'fields.bin').read_bytes()).hexdigest()==digest
    print('PASS: checkpoint continuation; normalized trajectory difference',error)
