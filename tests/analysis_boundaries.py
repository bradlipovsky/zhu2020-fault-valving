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
from analyze import analyze,fit_migration,velocity_crossings
from prepare_panels import run_panel

def integrated(increments):
    # Increment i belongs to the interval from saved time i to time i+1.
    return np.r_[0.,np.cumsum(increments)[:-1]]

def main():
    assert np.allclose(velocity_crossings(np.arange(3,14,2)*1000,
        np.array([1e-10,-1e-8,1e-10,1e-8,1e-10,1e-8]),1e-9),[4000,8000,12000])
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
            records['fields'][:,0,:]=(.5*integrated(active))[:,None]
            records['fields'][:,1,:]=speed[:,None];records['fields'][:,2,:]=2e7
            records['fields'][:,3:5,:]=1e-16;records['fields'][:,5,:]=3e-9
            maximum_depth=np.full(len(t),15000.)
            def write_fixture():
                with (folder/'fields.bin').open('wb') as f:
                    f.write(struct.pack('<8sQ',b'ZHUIND01',len(z)));z.tofile(f);records.tofile(f)
                np.savetxt(str(folder/'history.csv'),np.column_stack([t,speed,maximum_depth,np.full(len(t),2e7)]),
                    delimiter=',',header='time_s,vmax_m_s,z_vmax_m,min_effective_pa',comments='')
            write_fixture()
            (folder/'config.cfg').write_text('Vp = 1e-9\nseismic_threshold = 0.001\nT = 1e8\ninflux = 3e-9\n')
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
            records['fields'][:,0,:]=(.5*integrated(np.isin(np.arange(12),[2,3])))[:,None]
            records['fields'][:,0,2:]+=(.5*integrated(np.isin(np.arange(12),[6,7])))[:,None]
            records['fields'][:,0,:]+=(.5*integrated(np.isin(np.arange(12),[10,11])))[:,None]
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
            records['fields'][:,0,2:]=integrated(increments)[:,None]
            write_fixture();result=analyze('fixture')
            assert len(result['events'])==3 and all(e['complete'] for e in result['events'])
            partial=result['phase_diagnostics']['first_partial_ruptures']
            assert len(partial)==2 and np.allclose([e['max_slip_m'] for e in partial],[.2,.4])
            # A brief threshold excursion with sub-centimetre slip remains in
            # the catalog, but does not count as another resolved rupture.
            increments[[2,3]]=.001;records['fields'][:,0,2:]=integrated(increments)[:,None]
            write_fixture();result=analyze('fixture')
            assert len(result['events'])==3 and result['selected_event_count']==3
            assert result['selected_small_ruptures']==2 and result['selected_unresolved_intervals']==1
            partial=result['phase_diagnostics']['first_partial_ruptures']
            assert np.allclose([e['max_slip_m'] for e in partial],[.4,2])
            # Earlier creep must not inflate the slip of a subsequent brief
            # threshold excursion. Only the crossing profiles bound its slip.
            speed[:]=1e-10;speed[[2,3]]=.002;records['fields'][:,1,:]=speed[:,None]
            increments[:]=0;increments[0]=.02;increments[2:4]=.0002
            records['fields'][:,0,:]=integrated(increments)[:,None]
            write_fixture();result=analyze('fixture');event=result['events'][0]
            assert len(result['events'])==1 and result['selected_unresolved_intervals']==1
            assert event['onset_bracket_s']==[1.,2.] and event['slip_start_s']==2 and event['slip_end_s']==4
            assert abs(event['max_slip_m']-.0004)<1e-8 and not event['rupture_intervals_m']
            # A seismic event outside the saved depth range still invalidates a slow-slip
            # episode. The accepted-step history carries that whole-fault maximum.
            t*=31557600;records['time']=t;speed[:]=1e-10;speed[2:6]=2e-9
            records['fields'][:,1,:]=speed[:,None]
            records['fields'][:,0,:]=(integrated(speed)*31557600)[:,None]
            speed[4]=.01;maximum_depth[4]=40000;write_fixture()
            result=analyze('fixture')
            assert all(not d['episodes'] for d in result['slow_slip'].values())
            speed[4]=2e-9;maximum_depth[4]=15000;write_fixture()
            result=analyze('fixture')
            assert all(len(d['episodes'])==1 for d in result['slow_slip'].values())
            # An interior pulse must not hide a simultaneous leading boundary.
            records['fields'][:,1,:]=[1e-10,1e-10,1e-8,1e-10,1e-8]
            speed[:]=1e-8;write_fixture();result=analyze('fixture')
            fronts=np.genfromtxt(str(folder/'front.csv'),delimiter=',',names=True)
            assert np.allclose(fronts['leading_front_depth_m'],7500,atol=.001)
            assert np.allclose(fronts['deepest_front_depth_m'],17500,atol=.001)
            # Three closed large events delimit two cycles. A slip panel must
            # retain both intervening intervals, including their closing ruptures.
            t=np.arange(14,dtype=float);active=np.isin(np.arange(14),[2,3,6,7,10,11])
            speed=np.where(active,.01,1e-10);maximum_depth=np.full(len(t),15000.)
            records=np.zeros(len(t),dtype=dtype);records['time']=t;records['step']=np.arange(len(t))
            records['fields'][:,0,:]=(.5*integrated(active))[:,None]
            records['fields'][:,1,:]=speed[:,None];records['fields'][:,2,:]=2e7
            records['fields'][:,3:5,:]=1e-16;records['fields'][:,5,:]=3e-9
            (folder/'provenance.json').write_text('{"purpose":"synthetic window verification"}')
            write_fixture();analyze('fixture')
            spec=dict(cases=['fixture'],kind='slip_profiles',field=0)
            panel,metadata=run_panel('F2d',spec)
            assert metadata['window_s']==[4.,12.] and metadata['complete_cycle_count']==2
            assert np.array_equal(panel['cycle_bounds_s'],[[4.,8.],[8.,12.]])
            assert np.all(panel['values'][0]==0) and np.all(panel['values'][-1]==2)
            active[10:]=False;speed[10:]=1e-10
            records['fields'][:,1,:]=speed[:,None]
            records['fields'][:,0,:]=(.5*integrated(active))[:,None]
            write_fixture();analyze('fixture');panel,metadata=run_panel('F2d',spec)
            assert metadata['window_s']==[4.,8.] and metadata['complete_cycle_count']==1
            assert metadata['requested_cycle_count']==2 and np.all(panel['values'][-1]==1)
        finally:os.chdir(str(original))
    print('CATALOG BOUNDARY CHECKS PASSED: truncated events cannot close complete cycles.')
    print('PHASE DIAGNOSTIC CHECKS PASSED: connected rupture depths and pressure-minimum timing.')
    print('CATALOG ORDER CHECKS PASSED: first partials retained despite a larger later event.')
    print('THRESHOLD CATALOG CHECKS PASSED: retain intervals without resolved slip footprints.')
    print('CROSSING SNAPSHOT CHECKS PASSED: preceding creep excluded from event slip.')
    print('ASEISMIC CLASSIFICATION CHECKS PASSED: use the whole-fault accepted-step maximum.')
    print('MIGRATION FIT CHECKS PASSED: speed, depth selection, and seismic interruption.')
    print('FRONT BRANCH CHECKS PASSED: preserve leading and deeper crossings separately.')
    print('SLIP PROFILE WINDOW CHECKS PASSED: two cycles retained; missing cycles counted.')

if __name__=='__main__':main()
