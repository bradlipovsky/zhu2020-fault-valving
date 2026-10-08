#!/usr/bin/env python3
"""Measure events and cycle scales without importing published model data."""
from pathlib import Path
import argparse
import json
import numpy as np
from common import YEAR, fields, history, configuration, save_json, MAIN_CASES

def fit_migration(front, depth_min, depth_max, seismic_threshold):
    segments=[];group=[]
    def close_group():
        if len(group)>=5:
            a=np.array(group);span=a[0,1]-a[-1,1];duration=a[-1,0]-a[0,0]
            if span>=1000 and duration>=.1:
                coef=np.polyfit(a[:,0],a[:,1],1);prediction=np.polyval(coef,a[:,0])
                denom=np.sum((a[:,1]-a[:,1].mean())**2)
                r2=1-np.sum((a[:,1]-prediction)**2)/denom if denom else 0
                segments.append(dict(start_year=float(a[0,0]),end_year=float(a[-1,0]),
                    depth_start_m=float(a[0,1]),depth_end_m=float(a[-1,1]),
                    speed_m_per_year=float(-coef[0]),r2=float(r2),samples=len(a)))
    for time,depth,speed in front:
        usable=depth is not None and depth_min<=depth<=depth_max and speed<seismic_threshold
        continuous=usable and (not group or (time-group[-1][0]<=.25 and -250<=depth-group[-1][1]<=30))
        if not usable or not continuous:
            close_group();group=[]
        if usable:group.append([time,depth])
    close_group()
    return segments

def analyze(case):
    z,r=fields(case); h=history(case); t=np.asarray(r['time']); v=np.asarray(r['fields'][:,1,:])
    completed=json.loads(Path('data',case,'completed.json').read_text())
    cfg=configuration(Path('data',case,'config.cfg'))
    seismic=h['vmax_m_s']>cfg['seismic_threshold']
    starts=np.flatnonzero(seismic & ~np.r_[False,seismic[:-1]])
    stops=np.flatnonzero(seismic & ~np.r_[seismic[1:],False])+1
    events=[]
    for i,j in zip(starts,stops):
        t0=float(h['time_s'][max(i-1,0)]); t1=float(h['time_s'][min(j,len(h)-1)])
        # The executable saves each threshold crossing. Earlier snapshots can
        # contain substantial aseismic nucleation slip outside a brief excursion.
        slip_start=float(h['time_s'][i]);slip_end=t1
        k0=int(np.searchsorted(t,slip_start));k1=int(np.searchsorted(t,slip_end))
        if k0>=len(t) or k1>=len(t) or t[k0]!=slip_start or t[k1]!=slip_end:
            raise RuntimeError('Missing threshold-crossing snapshot: '+case)
        slip=np.asarray(r['fields'][k1,0,:]-r['fields'][k0,0,:])
        # A fast rupture tip may pass between the saved 0.5 s snapshots.
        # Integrated event slip is retained even when a local velocity peak is missed.
        rupture=(np.abs(slip)>=.01) & (z<=25000)
        depths=z[rupture]; imax=i+np.argmax(h['vmax_m_s'][i:j])
        shallow=float(depths.min()) if len(depths) else None
        deep=float(depths.max()) if len(depths) else None
        # Do not join disconnected shallow quasi-dynamic motion to a deep rupture.
        edges=np.diff(np.r_[False,rupture,False].astype(int))
        intervals=[[float(z[i]),float(z[j-1])] for i,j in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1))]
        large=any(top<2000 and bottom-top>10000 for top,bottom in intervals)
        events.append(dict(start_s=t0,end_s=t1,slip_start_s=slip_start,slip_end_s=slip_end,
            onset_bracket_s=[t0,slip_start],peak_s=float(h['time_s'][imax]),
            peak_velocity=float(h['vmax_m_s'][imax]),nucleation_depth_m=float(h['z_vmax_m'][i]),
            shallow_m=shallow,deep_m=deep,rupture_intervals_m=intervals,
            max_slip_m=float(np.max(np.abs(slip))),large=large,complete=bool(i>0 and j<len(h))))
    large=[e for e in events if e['large'] and e['complete']]
    cycles=[[large[i]['end_s'],large[i+1]['end_s']] for i in range(len(large)-1)]
    # Fixed, disclosed rule: final complete shallow-reaching large-event cycle.
    selected=cycles[-1] if cycles else [float(t[0]),float(t[-1])]
    mask=(t>=selected[0])&(t<=selected[1]); ids=np.flatnonzero(mask)
    if len(ids)<2: raise RuntimeError('No usable cycle samples: '+case)
    Ne=np.asarray(r['fields'][ids,2,:]); kval=np.asarray(r['fields'][ids,3,:]); flux=np.asarray(r['fields'][ids,5,:])
    selected_events=[e for e in events if e['complete'] and selected[0]<e['peak_s']<=selected[1]]
    sample_depths={}
    for d in [5000,10000,15000,20000]:
        k=int(np.argmin(abs(z-d)))
        sample_depths[str(d)]=dict(sample_depth_m=float(z[k]),
            effective_min_mpa=float(Ne[:,k].min()/1e6),effective_max_mpa=float(Ne[:,k].max()/1e6),
            effective_range_mpa=float(np.ptp(Ne[:,k])/1e6),
            permeability_min=float(kval[:,k].min()),permeability_max=float(kval[:,k].max()),
            flux_min=float(flux[:,k].min()),flux_max=float(flux[:,k].max()))
    # Migration speed: consecutive upward V=Vp crossings connected to deep creep.
    front=[]
    zz=z[z<=25000]
    snapshot_history=np.searchsorted(h['time_s'],t)
    snapshot_maximum=h['vmax_m_s'][snapshot_history]
    for it in ids:
        vv=np.abs(v[it,:len(zz)])
        crossings=np.flatnonzero((vv[:-1]<cfg['Vp'])&(vv[1:]>=cfg['Vp']))
        crossings=crossings[(zz[crossings]>2000)&(zz[crossings]<23000)]
        if len(crossings):
            k=crossings[-1]
            weight=(np.log(cfg['Vp'])-np.log(max(vv[k],1e-35)))/(np.log(max(vv[k+1],1e-35))-np.log(max(vv[k],1e-35)))
            depth=float(zz[k]+weight*(zz[k+1]-zz[k]))
        else: depth=None
        front.append([float((t[it]-selected[0])/YEAR),depth,float(snapshot_maximum[it])])
    segments=fit_migration(front,13000,20000,cfg['seismic_threshold'])
    shallow_segments=fit_migration(front,2000,10000,cfg['seismic_threshold'])
    intervals=np.diff([e['end_s'] for e in large])/YEAR
    slow_slip={}
    for depth in [15000,18000,20000]:
        # The article describes centimetres of slip over about a year, hence
        # speeds near Vp. A 10 Vp cutoff would miss such low-amplitude pulses.
        k=int(np.argmin(abs(z-depth)));speed=np.abs(v[ids,k]);above=speed>=1.1*cfg['Vp']
        changes=np.diff(np.r_[False,above,False].astype(int));episodes=[]
        for begin,end in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)):
            # Require two threshold crossings inside the window. A truncated
            # episode is not assigned a duration or recurrence peak.
            if begin==0 or end==len(ids):continue
            subset=ids[begin:end];duration=(t[subset[-1]]-t[subset[0]])/YEAR
            lo=np.searchsorted(h['time_s'],t[subset[0]])
            hi=np.searchsorted(h['time_s'],t[subset[-1]],side='right')
            if duration<.01 or np.max(h['vmax_m_s'][lo:hi])>=cfg['seismic_threshold']:continue
            peak=subset[np.argmax(speed[begin:end])]
            episodes.append(dict(start_s=float(t[subset[0]]),end_s=float(t[subset[-1]]),
                peak_s=float(t[peak]),duration_years=float(duration),peak_velocity=float(abs(v[peak,k])),
                net_slip_m=float(r['fields'][subset[-1],0,k]-r['fields'][subset[0],0,k])))
        spacings=np.diff([e['peak_s'] for e in episodes])/YEAR
        slow_slip[str(depth)]=dict(sample_depth_m=float(z[k]),threshold_m_s=1.1*cfg['Vp'],episodes=episodes,
            median_interval_years=float(np.median(spacings)) if len(spacings) else None,
            interval_coefficient_of_variation=float(np.std(spacings)/np.mean(spacings)) if len(spacings)>1 else None,
            median_duration_years=float(np.median([e['duration_years'] for e in episodes])) if episodes else None,
            median_net_slip_m=float(np.median([e['net_slip_m'] for e in episodes])) if episodes else None)
    k10=int(np.argmin(abs(z-10000)));peak10=int(np.argmax(Ne[:,k10]))
    partial=[e for e in selected_events if not e['large'] and e['rupture_intervals_m']]
    partial=partial[:2]
    ruptures=[]
    for event in partial:
        top,bottom=max(event['rupture_intervals_m'],key=lambda interval:interval[1]-interval[0])
        ruptures.append(dict(peak_since_window_start_years=(event['peak_s']-selected[0])/YEAR,
            footprint_top_m=top,footprint_bottom_m=bottom,max_slip_m=event['max_slip_m'],
            peak_velocity=event['peak_velocity']))
    phase_diagnostics=dict(sample_depth_m=float(z[k10]),
        maximum_effective_stress_since_window_start_years=float((t[ids[peak10]]-selected[0])/YEAR),
        pressure_drop_from_window_start_mpa=float((Ne[peak10,k10]-Ne[0,k10])/1e6),
        first_partial_ruptures=ruptures,
        rupture_selection='First two complete partial events in chronological order; each depth interval is its longest connected >=1 cm footprint. The full event catalog is retained.')
    result=dict(case=case,configuration=cfg,events=events,complete_cycles_s=cycles,
        resolution={key:completed[key] for key in ['n','dynamic_nodes','dz_m','minimum_Lb_cells','minimum_hstar_cells']},
        event_definition='Catalog every maximum accepted-step speed excursion > 1e-3 m/s; measure slip between saved threshold crossings; resolved rupture requires a connected >= 0.01 m footprint; large if top < 2 km and span > 10 km',
        selected_cycle_s=selected,selection_rule='last complete large-event cycle; entire run if none',
        selected_cycle_complete=bool(cycles),large_event_count=len(large),
        recurrence_years=intervals.tolist(),
        median_recurrence_years=float(np.median(intervals)) if len(intervals) else None,
        selected_event_count=len(selected_events),
        selected_small_ruptures=sum(not e['large'] and bool(e['rupture_intervals_m']) for e in selected_events),
        selected_unresolved_intervals=sum(not e['rupture_intervals_m'] for e in selected_events),
        sample_depths=sample_depths,migration_segments=segments,shallow_migration_segments=shallow_segments,
        slow_slip=slow_slip,phase_diagnostics=phase_diagnostics,
        truncated_event_count=sum(not e['complete'] for e in events),
        minimum_effective_pa=float(h['min_effective_pa'].min()),
        maximum_velocity=float(h['vmax_m_s'].max()),
        negative_velocity_samples=int(np.sum(v<0)),
        minimum_flux=float(np.min(r['fields'][:,5,:])),
        final_time_years=float(t[-1]/YEAR))
    save_json(Path('data',case,'analysis.json'),result)
    np.savetxt(Path('data',case,'front.csv'),np.array([[a,np.nan if b is None else b,c] for a,b,c in front]),
        delimiter=',',header='cycle_time_years,front_depth_m,maximum_velocity_m_s',comments='')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('cases',nargs='*',default=MAIN_CASES);args=parser.parse_args()
    for case in args.cases:
        r=analyze(case)
        print(case,'large events',r['large_event_count'],'recurrence years',r['recurrence_years'],
              'small resolved ruptures in window',r['selected_small_ruptures'],
              'unresolved threshold intervals',r['selected_unresolved_intervals'])
