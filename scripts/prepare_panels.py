#!/usr/bin/env python3
"""Select and store each panel's numerical inputs from fresh C++ output."""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
from common import fields, history, configuration, save_json, YEAR
from inventory import PANELS
from analyze import analyze

def analysis(case):
    path=Path('data',case,'analysis.json')
    return json.loads(path.read_text()) if path.exists() else analyze(case)

def run_panel(name,spec):
    case=spec['cases'][0];z,r=fields(case);a=analysis(case);cfg=a['configuration']
    t=np.asarray(r['time']);window=list(a['selected_cycle_s']);kind=spec['kind'];field=spec['field']
    reason=a['selection_rule'];complete=a['selected_cycle_complete']; depth=(0,25000)
    if name.startswith('S2'):
        offset=2 if name.endswith('a') else 3
        cycles=a['complete_cycles_s']
        if len(cycles)>=offset:window=cycles[-offset];reason='Earlier complete cycle at fixed reverse index {}'.format(offset)
        else:reason='Insufficient earlier cycles; full simulation shown';window=[float(t[0]),float(t[-1])];complete=False
    if name.startswith('F4'):
        # Fixed two-year pre-rupture window; no search for the closest published image.
        window=[max(window[0],window[1]-2*YEAR),window[1]];depth=(2000,10000)
        reason='Final two years preceding the selected closing large earthquake'
    if kind=='depth_histories':
        cycles=a['complete_cycles_s']
        window=[cycles[-2][0],cycles[-1][1]] if len(cycles)>=2 else [float(t[0]),float(t[-1])]
        reason='Last two complete cycles; full simulation if fewer than two'
    ti=np.flatnonzero((t>=window[0])&(t<=window[1]))
    zi=np.flatnonzero((z>=depth[0])&(z<=depth[1]))
    result=dict(depth_m=z[zi],time_s=t[ti],step=np.asarray(r['step'][ti]),
        origin_s=np.array(window[0]),values=np.asarray(r['fields'][ti,field,:][:,zi]),
        steady_values=np.asarray(r['fields'][0,field,zi]),
        case=np.array(case),field=np.array(field),Vp=np.array(cfg['Vp']),influx=np.array(cfg['influx']))
    if kind=='depth_histories':
        # Use accepted-step histories to preserve brief seismic maxima.
        h=history(case);mask=(h['time_s']>=window[0])&(h['time_s']<=window[1]);hh=h[mask]
        lookup={1:'velocity',2:'effective',3:'permeability',5:'flux'}
        vals=np.column_stack([hh[lookup[field]+'_'+str(d)] for d in [5000,10000,15000,20000]])
        # Store all accepted samples; no interpolation or copied reference curves.
        result.update(depth_m=np.array([5000,10000,15000,20000]),time_s=hh['time_s'],step=hh['step'],values=vals,
            steady_values=np.array([h[lookup[field]+'_'+str(d)][0] for d in [5000,10000,15000,20000]]))
    if kind in ['profiles','slip_profiles']:
        interval=1.5 if cfg['T']<=1e8 else 4
        requested=np.arange(window[0],window[1],interval*YEAR)
        selected=[];labels=[]
        # Regular interseismic contours.
        for target in requested:
            j=int(np.clip(np.searchsorted(t,target),0,len(t)-1))
            selected.append(j);labels.append('interseismic')
        if field==2:
            for target in np.arange(window[0],min(window[0]+10*YEAR,window[1]),(.5 if cfg['T']<=1e8 else 1)*YEAR):
                selected.append(int(np.clip(np.searchsorted(t,target),0,len(t)-1)));labels.append('early')
        for event in a['events']:
            if event['complete'] and window[0]<event['peak_s']<=window[1]:
                for target in np.arange(event['start_s'],event['end_s'],1.):
                    selected.append(int(np.clip(np.searchsorted(t,target),0,len(t)-1)));labels.append('coseismic')
        result.update(profile_time_s=t[selected],profile_values=np.asarray(r['fields'][selected,field,:][:,zi]),
            profile_class=np.asarray(labels))
    if field==0:
        # Cumulative slip relative to cycle onset, not relative to a published restart.
        j=int(np.clip(np.searchsorted(t,window[0]),0,len(t)-1));zero=np.asarray(r['fields'][j,0,zi])
        result['values']=result['values']-zero
        if 'profile_values' in result:result['profile_values']=result['profile_values']-zero
    return result,dict(case=case,config_file='configs/'+case+'.cfg',parameters=cfg,
        window_s=window,selection=reason,complete_cycle=complete,
        source_provenance=json.loads(Path('data',case,'provenance.json').read_text()))

def prepare(name):
    spec=PANELS[name];kind=spec['kind'];metadata=dict(spec)
    if not spec['cases']:
        profile=np.genfromtxt('data/laws/steady_profile.csv',delimiter=',',names=True)
        result={key:profile[key] for key in profile.dtype.names}
        if kind=='stress_law':
            a=np.genfromtxt('data/laws/stress_law.csv',delimiter=',',names=True)
            result={key:a[key] for key in a.dtype.names}
        if kind=='evolution_law':
            for file in ['healing_law','slip_law']:
                a=np.genfromtxt('data/laws/'+file+'.csv',delimiter=',',names=True)
                result.update({key:a[key] for key in a.dtype.names})
        metadata.update(parameters=configuration('configs/baseline.cfg'),config_file='configs/baseline.cfg',
            selection='Direct constitutive laws or independently calculated steady state')
    elif kind=='fronts':
        result={};selections={}
        for case in spec['cases']:
            a=analysis(case);z,r=fields(case);window=a['selected_cycle_s']
            ti=np.flatnonzero((r['time']>=window[0])&(r['time']<=window[1]));zi=np.flatnonzero(z<=25000)
            result[case+'_depth_m']=z[zi];result[case+'_time_years']=(r['time'][ti]-window[0])/YEAR
            result[case+'_velocity']=np.asarray(r['fields'][ti,1,:][:,zi]);result[case+'_Vp']=np.array(a['configuration']['Vp'])
            selections[case]=dict(window_s=window,configuration=a['configuration'],
                complete_cycle=a['selected_cycle_complete'],selection=a['selection_rule'])
        metadata['selections']=selections
    else:
        result,details=run_panel(name,spec);metadata.update(details)
    status='partial'; discrepancy='Exact friction and initialization are not fully specified; numerical trajectory differs.'
    if name in ['F1b','F1c']:
        status='independently reproduced';discrepancy='Analytic laws agree with numerical limiting-case verification.'
    elif name=='F1a':
        discrepancy='New schematic, not a simulation reproduction; friction profiles are declared approximations.'
    elif name=='F1d':discrepancy='Normal stress gradient is assumed; steady Darcy balance is independently verified.'
    elif name.startswith('F4'):
        a=analysis('baseline');small=[e for e in a['events'] if e['complete'] and not e['large']
            and metadata['window_s'][0]<e['peak_s']<=metadata['window_s'][1]
            and any(top<10000 and bottom>2000 for top,bottom in e['rupture_intervals_m'])]
        metadata['middepth_small_events']=small
        discrepancy='{} small ruptures intersect 2--10 km in the fixed two-year window; exact swarm sequence differs.'.format(len(small))
        if len(small)<3:
            status='not reproduced'
            discrepancy+=' Fewer than three such events; generated attempt shown.'
    elif name.startswith('S2') and not metadata.get('complete_cycle'):
        status='not reproduced';discrepancy='Too few complete cycles to produce this alternate-cycle comparison.'
    elif name.startswith('F5'):
        pulses=analysis('short')['slow_slip'];depths=['15000','18000','20000']
        counts=[str(len(pulses[d]['episodes'])) for d in depths]
        intervals=[pulses[d]['median_interval_years'] for d in depths]
        discrepancy='Complete aseismic episodes above 1.1 Vp at 15/18/20 km: '+ '/'.join(counts)+'. '
        discrepancy+='Median intervals: '+', '.join('{:.3g}'.format(x) if x is not None else 'unmeasured' for x in intervals)+' yr. '
        discrepancy+='Durations and slip are tabulated in the report; exact pulse sequence is not recovered.'
        metadata['slow_slip_diagnostics']=pulses
    elif name=='F6':
        measured=[]
        for case in spec['cases']:
            pair=[]
            for field in ['leading_migration_segments','migration_segments']:
                rates=[s['speed_m_per_year'] for s in analysis(case)[field] if s['r2']>=.8]
                pair.append('{:.3g}'.format(float(np.median(rates))) if rates else 'unmeasured')
            measured.append('/'.join(pair))
        discrepancy='Leading/deepest contour medians at 13--20 km, short through verylong: '+', '.join(measured)+' m/year; published targets 2500, 380, 120, 30. Both contour choices are retained.'
        shallow=[]
        for field in ['leading_shallow_migration_segments','shallow_migration_segments']:
            rates=[s['speed_m_per_year'] for s in analysis('baseline')[field] if s['r2']>=.8]
            shallow.append('{:.3g}'.format(float(np.median(rates))) if rates else 'unmeasured')
        discrepancy+=' Baseline shallow medians: '+ '/'.join(shallow)+' m/year versus 4560.'
    elif name.startswith('F2'):
        case=spec['cases'][0];a=analysis(case);target=50 if case=='reference' else 32
        recurrence=a['median_recurrence_years']
        discrepancy=('No complete large-event recurrence obtained.' if recurrence is None else
            'Median recurrence {:.3g} yr versus approximately {} yr in the displayed reference cycle ({:+.1f}%).'.format(
                recurrence,target,100*(recurrence/target-1)))
        discrepancy+=' Initial state and friction profiles remain uncertain.'
        if case=='baseline' and a['selected_cycle_complete']:
            ruptures=a['phase_diagnostics']['first_partial_ruptures']
            bounds=['{:.2f}--{:.2f}'.format(e['footprint_top_m']/1000,e['footprint_bottom_m']/1000) for e in ruptures]
            discrepancy+=' First partial footprints: '+(', '.join(bounds)+' km.' if bounds else 'none.')
    elif name.startswith('F3'):
        a=analysis('baseline');sample=a['sample_depths']['10000']
        if spec['field']==2:
            discrepancy='At 10 km, stress varies by {:.3g} MPa; article gives an order 10--20 MPa scale. Phase timing differs.'.format(sample['effective_range_mpa'])
        elif spec['field']==3:
            discrepancy='At 10 km, k spans {:.3g}--{:.3g} square metres; exact cycle phases are not recovered.'.format(sample['permeability_min'],sample['permeability_max'])
        else:
            discrepancy='At 10 km, maximum flux is {:.3g} m/s; the article gives an order 1e-7 m/s upper scale. Phase timing differs.'.format(sample['flux_max'])
    elif name.startswith(('S4','S6')):
        sample=analysis(spec['cases'][0])['sample_depths']['10000']
        discrepancy='At 10 km, stress range {:.3g} MPa and permeability ratio {:.3g}. The target is reduced valving with long healing; event sequence and mesh convergence remain uncertain.'.format(
            sample['effective_range_mpa'],sample['permeability_max']/sample['permeability_min'])
    metadata.update(status=status,main_discrepancy=discrepancy,
        model_source_sha256=hashlib.sha256(Path('src/model.cpp').read_bytes()).hexdigest(),
        pipeline_source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path('scripts',name) for name in ['common.py','analyze.py','prepare_panels.py','inventory.py','plot.py']]})
    Path('data/panels').mkdir(exist_ok=True)
    np.savez_compressed(spec['data_file'],**result)
    save_json('data/panels/'+name+'.json',metadata)
    print(name,status,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('panels',nargs='*');args=parser.parse_args()
    for name in args.panels or PANELS:prepare(name)
