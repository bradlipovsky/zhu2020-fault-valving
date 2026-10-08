#!/usr/bin/env python3
"""Synthetic catalog-boundary check; these arrays never supply scientific panels."""
from pathlib import Path
import json
import os
import struct
import sys
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'scripts'))
from analyze import analyze,fit_migration

def main():
    # Known translating front, followed by the same front interrupted by a quake.
    front=[[float(t),9000-4200*t,1e-8] for t in np.arange(60)*.025]
    fitted=fit_migration(front,2000,10000,.001)
    assert len(fitted)==1 and abs(fitted[0]['speed_m_per_year']-4200)<1e-8
    assert not fit_migration(front,13000,20000,.001)
    front[30][2]=.01;fitted=fit_migration(front,2000,10000,.001)
    assert len(fitted)==2 and all(abs(s['speed_m_per_year']-4200)<1e-8 for s in fitted)
    (ROOT/'.tmp').mkdir(exist_ok=True)
    original=Path.cwd()
    with tempfile.TemporaryDirectory(prefix='catalog-boundary-',dir=str(ROOT/'.tmp')) as temporary:
        try:
            os.chdir(temporary);folder=Path('data/fixture');folder.mkdir(parents=True)
            z=np.array([1000,5000,10000,15000,20000],dtype='<f8');t=np.arange(12,dtype=float)
            # Two closed events, followed by a third that remains active at the output boundary.
            active=np.isin(np.arange(12),[2,3,6,7,10,11]);speed=np.where(active,.01,1e-10)
            dtype=np.dtype([('time','<f8'),('step','<i8'),('fields','<f4',(6,len(z)))])
            records=np.zeros(len(t),dtype=dtype);records['time']=t;records['step']=np.arange(len(t))
            records['fields'][:,0,:]=(.5*np.cumsum(active))[:,None]
            records['fields'][:,1,:]=speed[:,None];records['fields'][:,2,:]=2e7
            records['fields'][:,3:5,:]=1e-16;records['fields'][:,5,:]=3e-9
            maximum_depth=np.full(len(t),15000.)
            def write_fixture():
                with (folder/'fields.bin').open('wb') as f:
                    f.write(struct.pack('<8sQ',b'ZHUIND01',len(z)));z.tofile(f);records.tofile(f)
                np.savetxt(str(folder/'history.csv'),np.column_stack([t,speed,maximum_depth,np.full(len(t),2e7)]),
                    delimiter=',',header='time_s,vmax_m_s,z_vmax_m,min_effective_pa',comments='')
            write_fixture()
            (folder/'config.cfg').write_text('Vp = 1e-9\nseismic_threshold = 0.001\nT = 1e8\n')
            (folder/'completed.json').write_text(json.dumps(dict(completed=True,n=16,dynamic_nodes=16,
                dz_m=1000,minimum_Lb_cells=5,minimum_hstar_cells=20)))
            result=analyze('fixture')
            assert [e['complete'] for e in result['events']]==[True,True,False]
            assert result['large_event_count']==2 and result['truncated_event_count']==1
            assert result['complete_cycles_s']==[[4.,8.]] and result['selected_cycle_s']==[4.,8.]
            # A file beginning during an event cannot manufacture a complete initial event either.
            records['fields'][:2,1,:]=.01;speed[:2]=.01
            write_fixture()
            result=analyze('fixture')
            assert [e['complete'] for e in result['events']]==[False,True,False]
            assert result['large_event_count']==1 and not result['selected_cycle_complete']
            assert result['truncated_event_count']==2
            # Now the second closed event slips only at depth, and a later pressure minimum
            # gives a known drainage diagnostic. The final open event remains excluded.
            speed[:2]=1e-10;records['fields'][:2,1,:]=1e-10
            records['fields'][:,0,:]=(.5*np.cumsum(np.isin(np.arange(12),[2,3])))[:,None]
            records['fields'][:,0,2:]+=(.5*np.cumsum(np.isin(np.arange(12),[6,7])))[:,None]
            records['fields'][:,0,:]+=(.5*np.cumsum(np.isin(np.arange(12),[10,11])))[:,None]
            records['fields'][9,2,:]=2.5e7
            write_fixture()
            result=analyze('fixture');phase=result['phase_diagnostics'];event=phase['first_partial_ruptures'][0]
            assert len(phase['first_partial_ruptures'])==1
            assert event['footprint_top_m']==10000 and event['footprint_bottom_m']==20000
            assert event['max_slip_m']==1 and phase['pressure_drop_from_window_start_mpa']==5
            assert phase['maximum_effective_stress_since_window_start_years']==9/31557600
            # A later, larger partial event must not replace either of the first
            # two chronological events in the phase-comparison table.
            speed[:]=1e-10;speed[[2,3,6,7,10]]=.01
            records['fields'][:,1,:]=speed[:,None];records['fields'][:,0,:]=0
            increments=np.zeros(len(t));increments[[2,3]]=.1;increments[[6,7]]=.2;increments[10]=2
            records['fields'][:,0,2:]=np.cumsum(increments)[:,None]
            write_fixture();result=analyze('fixture')
            assert len(result['events'])==3 and all(e['complete'] for e in result['events'])
            partial=result['phase_diagnostics']['first_partial_ruptures']
            assert len(partial)==2 and np.allclose([e['max_slip_m'] for e in partial],[.2,.4])
            # A seismic event outside the saved depth range still invalidates a slow-slip
            # episode. The accepted-step history carries that whole-fault maximum.
            t*=31557600;records['time']=t;speed[:]=1e-10;speed[2:6]=2e-9
            records['fields'][:,1,:]=speed[:,None]
            records['fields'][:,0,:]=(np.cumsum(speed)*31557600)[:,None]
            speed[4]=.01;maximum_depth[4]=40000;write_fixture()
            result=analyze('fixture')
            assert all(not d['episodes'] for d in result['slow_slip'].values())
            speed[4]=2e-9;maximum_depth[4]=15000;write_fixture()
            result=analyze('fixture')
            assert all(len(d['episodes'])==1 for d in result['slow_slip'].values())
        finally:os.chdir(str(original))
    print('CATALOG BOUNDARY CHECKS PASSED: truncated events cannot close complete cycles.')
    print('PHASE DIAGNOSTIC CHECKS PASSED: connected rupture depths and pressure-minimum timing.')
    print('CATALOG ORDER CHECKS PASSED: first partials retained despite a larger later event.')
    print('ASEISMIC CLASSIFICATION CHECKS PASSED: use the whole-fault accepted-step maximum.')
    print('MIGRATION FIT CHECKS PASSED: speed, depth selection, and seismic interruption.')

if __name__=='__main__':main()
