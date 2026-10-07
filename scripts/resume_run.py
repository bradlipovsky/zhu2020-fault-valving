"""Resume a stopped integration from its latest periodic checkpoint.

Run only after the simulation has stopped. Any output after the checkpoint is
preserved in build/interrupted before truncation to the accepted saved state.
The time embedded in the atomically replaced binary checkpoint is authoritative.
"""
from pathlib import Path
import argparse
import datetime
import json
import os
import shutil
import struct
import subprocess

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory',type=Path)
args=parser.parse_args()
path=args.directory.resolve();meta=json.loads((path/'metadata.json').read_text())
state=path/'checkpoint_state.bin'
assert state.stat().st_size==4*meta['n']*8+8, 'Checkpoint does not match metadata'
with state.open('rb') as stream:
    stream.seek(-8,2);time=struct.unpack('d',stream.read(8))[0]
assert time<meta['years']*31536000, 'Run already completed'
fields=path/'fields.bin';width=8*(1+7*meta['nz']);keep=None
with fields.open('rb') as stream:
    for j in range(fields.stat().st_size//width):
        stream.seek(j*width);t=struct.unpack('d',stream.read(8))[0]
        if abs(t-time)<1e-5:keep=(j+1)*width;break
assert keep is not None,'No field record at the checkpoint time'
stamp=datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S%f')
backup=ROOT/'build'/'interrupted'/(path.name+'-'+stamp)
backup.mkdir(parents=True)
for name in ['fields.bin','history.csv','metadata.json','checkpoint_state.bin','checkpoint.json']:
    if (path/name).exists():shutil.copy2(str(path/name),str(backup/name))
with fields.open('r+b') as stream:stream.truncate(keep)
lines=(path/'history.csv').read_text().splitlines()
lines=[lines[0]]+[line for line in lines[1:] if float(line.split(',')[0])<=time]
(path/'history.csv').write_text('\n'.join(lines)+'\n')
command=[str(ROOT/'build'/'valving'),'--output',str(path),'--n',str(meta['n']-1),
    '--T',str(meta['T_s']),'--years',str(meta['years']),'--height',str(meta['height_m']),
    '--rtol',str(meta['rtol']),'--max-dt',str(meta['max_dt_s']),
    '--stride',str(meta.get('stride',10)),
    '--perturb',str(meta['initial_state_perturbation']),
    '--state',str(state),'--start-year',repr(time/31536000),'--append']
if meta['initial_profile']:command+=['--profile',meta['initial_profile']]
if not meta['coupled']:command+=['--fixed']
print('Preserved interrupted output in',backup,flush=True)
print('Resuming at year',time/31536000,flush=True)
os.chdir(str(ROOT))
subprocess.run(command,check=True)
