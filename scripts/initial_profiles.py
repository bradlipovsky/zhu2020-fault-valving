"""Recover approximate restart states from the authors' partial archive.

The aging-law state and shear traction were not archived. Reconstruct theta
using piecewise constant slip rates from archived slip increments, starting
with theta=1 s. State memory after 4 years is quantified in the metadata.
These profiles support archive-assisted reruns, not from-scratch replication.
"""
from pathlib import Path
import json
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
YEAR=31536000

def make(case,offset=4):
    with h5py.File(ROOT/'reference'/'data'/(case+'.mat')) as f:
        z=f['depth'][:]*1000;t=f['time'][:].ravel();j=np.argmin(abs((t-t[0])/YEAR-offset))
        slip=f['slip'][:,:j+1];theta=np.ones(len(z));memory=np.ones(len(z))
        for i in range(1,j+1):
            ds=slip[:,i]-slip[:,i-1]
            if np.any(ds<0):raise ValueError('Reverse slip requires a signed aging-state reconstruction')
            u=ds/.002;decay=np.exp(-u)
            integral=np.ones_like(u);sel=u>1e-10;integral[sel]=-np.expm1(-u[sel])/u[sel]
            theta=theta*decay+(t[i]-t[i-1])*integral;memory*=decay
        a=np.interp(z,[0,14900,27600,60000],[.0105,.03,.07,.173])
        b=np.interp(z,[0,13600,14900,27600,60000],[.02,.0378,.0356,.0375,.0424])
        psi=.6+b*np.log(1e-6*theta/.002)
        p=f['pore_pressure'][:,j]*1e6;k=f['permeability'][:,j]*1e6
        ks=1e-19+(k-1e-19)*np.exp((22050*z-p)/30e6)
        v=f['slip_velocity'][:,j];eff=f['effective_normal_stress'][:,j]*1e6
        w=np.log(v/(2e-6))+psi/a
        friction=a*np.where(w>30,w+np.log(2),np.arcsinh(np.exp(np.minimum(w,30))))
        tau=eff*friction+4.68e6*v
        path=ROOT/'inputs'/(case+'_restart.csv')
        np.savetxt(path,np.c_[z,psi,ks,p,tau],delimiter=',',header='z_m,psi,kstar_m2,p_Pa,tau_Pa',comments='')
        meta=dict(source=case+'.mat',source_index=int(j),offset_years=float((t[j]-t[0])/YEAR),
            absolute_time_s=float(t[j]),theta_initial_s=1,
            max_relative_theta_change_if_initial_2e6_s=float(np.max(1999999*memory/theta)),
            state_method='Exact constant-rate aging update per archived slip increment',
            pressure_units='Pa',permeability_units='m^2',stress_units='Pa')
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        print(case,meta,flush=True)

if __name__=='__main__':
    for name in ['T1e7','T1e8','T1e9','T1e10']:make(name)
