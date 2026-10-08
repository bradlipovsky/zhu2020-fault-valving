"""Independently check the km^2 permeability conversion using deep Darcy flux."""
from pathlib import Path
import json
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
results={}
for case in ['T1e7','T1e8','T1e9','T1e10']:
    with h5py.File(ROOT/'reference'/'data'/(case+'.mat')) as f:
        z=f['depth'][:]*1000
        # Below the taper, pressure and flux are nearly depth-independent in
        # effective-stress coordinates. Avoid upper/lower boundary stencils.
        indices=np.flatnonzero((z>100000)&(z<400000))
        errors=[]
        for j in [0,f['pore_pressure'].shape[1]//2,-1]:
            p=f['pore_pressure'][:,j]*1e6
            k=f['permeability'][:,j]*1e6
            q=f['flux'][:,j]
            predicted=k/1e-4*(np.gradient(p,z)-9800)
            errors.extend(np.abs(predicted[indices]/q[indices]-1))
        error=float(max(errors))
        assert error<1e-6,(case,error)
        results[case]={'maximum_relative_deep_Darcy_flux_error':error}
(ROOT/'results'/'archive_units_check.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS: Darcy flux confirms the archived pressure, depth, and permeability units')
