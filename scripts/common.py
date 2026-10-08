"""Read only output written by the independent C++ executable."""
from pathlib import Path
import json
import math
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

def physical_depths(cfg):
    """Reconstruct the configured physical mesh for output-coordinate metadata."""
    cfg={'n':32768,'Lz':500000.,'surface_ratio':1.,'surface_scale':1000.,
         'fine_depth':0.,'maximum_spacing':2000.,'stretch_scale':10000.,**cfg}
    n=int(cfg['n']);h=cfg['Lz']/n;nodes=[]
    j=int((cfg['surface_ratio']-1)/2)
    while True:
        nodes.append(j)
        if j==n-1:break
        z=(j+.5)*h
        ratio=max(1.,cfg['surface_ratio']*math.exp(-z/cfg['surface_scale']))
        if cfg['fine_depth']>0 and z>cfg['fine_depth']:
            ratio=math.exp(min(math.log(cfg['maximum_spacing']/h),(z-cfg['fine_depth'])/cfg['stretch_scale']))
        # The positive round-to-nearest convention matches std::round in C++.
        stride=max(1,int(math.floor(ratio+.5)));remaining=n-1-j
        j+=remaining if remaining<1.5*stride else stride
    return (np.asarray(nodes)+.5)*h

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
