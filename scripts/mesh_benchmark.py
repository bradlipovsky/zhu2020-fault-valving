#!/usr/bin/env python3
"""Repeat the controlled interseismic mesh comparisons; never supplies panel data."""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from configurations import BASE

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',default='.tmp/mesh-benchmark')
    parser.add_argument('--jobs',type=int,default=3)
    args=parser.parse_args();folder=Path(args.directory)
    folder.mkdir(parents=True,exist_ok=False)
    subprocess.run(['bash','scripts/bootstrap.sh'],check=True)
    subprocess.run(['cmake','--build','build','--target','mesh_compare'],check=True)
    cases={
        'uniform':dict(n=32768,years=5,fine_depth=0,surface_ratio=1,tolerance=.001),
        'graded':dict(n=32768,years=5,fine_depth=30000,surface_ratio=1,tolerance=.001),
        'fine_deep':dict(n=262144,years=1,fine_depth=30000,surface_ratio=1,tolerance=.001),
        'fine_surface':dict(n=262144,years=1,fine_depth=30000,surface_ratio=8,tolerance=.001),
        'fine_deep_tight':dict(n=262144,years=1,fine_depth=30000,surface_ratio=1,tolerance=1e-5),
        'fine_surface_tight':dict(n=262144,years=1,fine_depth=30000,surface_ratio=8,tolerance=1e-5)}
    def run(item):
        name,change=item
        values=dict(BASE,threads=1,maximum_spacing=2000,stretch_scale=10000,surface_scale=1000)
        values.update(change);config=folder/(name+'.cfg')
        config.write_text(''.join('{} = {:.16g}\n'.format(k,v) for k,v in values.items()))
        command=['build/valving',str(config),str(folder/name)]
        with (folder/(name+'.log')).open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        return name,dict(command=command,configuration=values,
            completion=json.loads((folder/name/'completed.json').read_text()))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        runs=dict(pool.map(run,cases.items()))
    comparisons={}
    for a,b in [('uniform','graded'),('fine_deep','fine_surface'),
                ('fine_deep_tight','fine_surface_tight'),('fine_surface','fine_surface_tight')]:
        command=['build/mesh_compare',str(folder/(a+'.cfg')),str(folder/a/'checkpoint.bin'),
                 str(folder/(b+'.cfg')),str(folder/b/'checkpoint.bin')]
        comparisons[a+'__'+b]=dict(command=command,results=json.loads(subprocess.check_output(command).decode()))
    result=dict(source_sha256=hashlib.sha256(Path('src/model.cpp').read_bytes()).hexdigest(),
        scope='Interseismic method check only; excluded from figure data.',runs=runs,comparisons=comparisons)
    (folder/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(folder/'results.json')

if __name__=='__main__':main()
