"""Read only output written by the independent C++ executable."""
from pathlib import Path
import json
import struct
import numpy as np

YEAR = 365.25*86400
FIELD_NAMES = ['slip', 'velocity', 'effective', 'permeability', 'kstar', 'flux']
MAIN_CASES = ['reference','baseline','short','long','long_reference','verylong','verylong_reference']

def configuration(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        line = line.split('#')[0]
        if '=' in line:
            key, value = line.split('=',1)
            result[key.strip()] = float(value)
    return result

def fields(case, complete=True):
    folder = Path('data',case)
    if complete and not (folder/'completed.json').exists():
        raise RuntimeError('Incomplete simulation: '+case)
    with (folder/'fields.bin').open('rb') as stream:
        if stream.read(8) != b'ZHUIND01':
            raise ValueError('Not output of the independent model')
        nz = struct.unpack('<Q',stream.read(8))[0]
        z = np.fromfile(stream,dtype='<f8',count=nz)
    dtype = np.dtype([('time','<f8'),('step','<i8'),('fields','<f4',(6,nz))])
    count = ((folder/'fields.bin').stat().st_size-16-8*nz)//dtype.itemsize
    records = np.memmap(folder/'fields.bin',dtype=dtype,mode='r',offset=16+8*nz,shape=(count,))
    return z, records

def history(case):
    return np.genfromtxt(Path('data',case,'history.csv'),delimiter=',',names=True)

def save_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
