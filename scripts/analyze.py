#!/usr/bin/env python3
"""Measure events and cycle scales without importing published model data."""
from pathlib import Path
import argparse
import json
import numpy as np
from common import YEAR, fields, history, configuration, save_json, MAIN_CASES

def analyze(case):
    z,r=fields(case); h=history(case); t=np.asarray(r['time']); v=np.asarray(r['fields'][:,1,:])
    cfg=configuration(Path('data',case,'config.cfg'))
    seismic=h['vmax_m_s']>=cfg['seismic_threshold']
    starts=np.flatnonzero(seismic & ~np.r_[False,seismic[:-1]])
    stops=np.flatnonzero(seismic & ~np.r_[seismic[1:],False])+1
    events=[]
    for i,j in zip(starts,stops):
        t0=float(h['time_s'][max(i-1,0)]); t1=float(h['time_s'][min(j,len(h)-1)])
        k0=max(0,np.searchsorted(t,t0)-1); k1=min(len(t)-1,np.searchsorted(t,t1))
        velocity=np.abs(v[k0:k1+1]); peaks=velocity.max(axis=0)
        slip=np.asarray(r['fields'][k1,0,:]-r['fields'][k0,0,:])
        rupture=(peaks>=cfg['seismic_threshold']) & (np.abs(slip)>=.01) & (z<=25000)
        depths=z[rupture]; imax=i+np.argmax(h['vmax_m_s'][i:j])
        shallow=float(depths.min()) if len(depths) else None
        deep=float(depths.max()) if len(depths) else None
        # Do not join disconnected shallow quasi-dynamic motion to a deep rupture.
        edges=np.diff(np.r_[False,rupture,False].astype(int))
        intervals=[[float(z[i]),float(z[j-1])] for i,j in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1))]
        large=any(top<2000 and bottom-top>10000 for top,bottom in intervals)
        events.append(dict(start_s=t0,end_s=t1,peak_s=float(h['time_s'][imax]),
            peak_velocity=float(h['vmax_m_s'][imax]),nucleation_depth_m=float(h['z_vmax_m'][i]),
            shallow_m=shallow,deep_m=deep,rupture_intervals_m=intervals,
            max_slip_m=float(np.max(np.abs(slip))),large=large))
    large=[e for e in events if e['large']]
    cycles=[[large[i]['end_s'],large[i+1]['end_s']] for i in range(len(large)-1)]
    # Fixed, disclosed rule: final complete shallow-reaching large-event cycle.
    selected=cycles[-1] if cycles else [float(t[0]),float(t[-1])]
    mask=(t>=selected[0])&(t<=selected[1]); ids=np.flatnonzero(mask)
    if len(ids)<2: raise RuntimeError('No usable cycle samples: '+case)
    Ne=np.asarray(r['fields'][ids,2,:]); kval=np.asarray(r['fields'][ids,3,:]); flux=np.asarray(r['fields'][ids,5,:])
    selected_events=[e for e in events if selected[0]<e['peak_s']<=selected[1]]
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
    for it in ids:
        vv=np.abs(v[it,:len(zz)])
        crossings=np.flatnonzero((vv[:-1]<cfg['Vp'])&(vv[1:]>=cfg['Vp']))
        crossings=crossings[(zz[crossings]>2000)&(zz[crossings]<23000)]
        if len(crossings):
            k=crossings[-1]
            weight=(np.log(cfg['Vp'])-np.log(max(vv[k],1e-35)))/(np.log(max(vv[k+1],1e-35))-np.log(max(vv[k],1e-35)))
            depth=float(zz[k]+weight*(zz[k+1]-zz[k]))
        else: depth=None
        front.append([float((t[it]-selected[0])/YEAR),depth,float(np.max(vv))])
    segments=[]; group=[]
    def close_group():
        if len(group)>=5:
            a=np.array(group); span=a[0,1]-a[-1,1]; duration=a[-1,0]-a[0,0]
            if span>=1000 and duration>=.1:
                coef=np.polyfit(a[:,0],a[:,1],1); prediction=np.polyval(coef,a[:,0])
                denom=np.sum((a[:,1]-a[:,1].mean())**2)
                r2=1-np.sum((a[:,1]-prediction)**2)/denom if denom else 0
                segments.append(dict(start_year=float(a[0,0]),end_year=float(a[-1,0]),
                    depth_start_m=float(a[0,1]),depth_end_m=float(a[-1,1]),
                    speed_m_per_year=float(-coef[0]),r2=float(r2),samples=len(a)))
    for row in front:
        time,depth,speed=row
        usable=depth is not None and 13000<=depth<=20000 and speed<cfg['seismic_threshold']
        continuous=not group or (time-group[-1][0]<=.25 and -250<=depth-group[-1][1]<=30) if usable else False
        if not usable or not continuous:
            close_group();group=[]
        if usable:group.append([time,depth])
    close_group()
    intervals=np.diff([e['end_s'] for e in large])/YEAR
    slow_slip={}
    maximum_speed=np.max(np.abs(v),axis=1)
    for depth in [15000,18000,20000]:
        k=int(np.argmin(abs(z-depth)));speed=np.abs(v[ids,k]);above=speed>=10*cfg['Vp']
        changes=np.diff(np.r_[False,above,False].astype(int));episodes=[]
        for begin,end in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)):
            # Require two threshold crossings inside the window. A truncated
            # episode is not assigned a duration or recurrence peak.
            if begin==0 or end==len(ids):continue
            subset=ids[begin:end];duration=(t[subset[-1]]-t[subset[0]])/YEAR
            if duration<.01 or np.max(maximum_speed[subset])>=cfg['seismic_threshold']:continue
            peak=subset[np.argmax(speed[begin:end])]
            episodes.append(dict(start_s=float(t[subset[0]]),end_s=float(t[subset[-1]]),
                peak_s=float(t[peak]),duration_years=float(duration),peak_velocity=float(abs(v[peak,k]))))
        spacings=np.diff([e['peak_s'] for e in episodes])/YEAR
        slow_slip[str(depth)]=dict(sample_depth_m=float(z[k]),threshold_m_s=10*cfg['Vp'],episodes=episodes,
            median_interval_years=float(np.median(spacings)) if len(spacings) else None,
            interval_coefficient_of_variation=float(np.std(spacings)/np.mean(spacings)) if len(spacings)>1 else None,
            median_duration_years=float(np.median([e['duration_years'] for e in episodes])) if episodes else None)
    result=dict(case=case,configuration=cfg,events=events,complete_cycles_s=cycles,
        selected_cycle_s=selected,selection_rule='last complete large-event cycle; entire run if none',
        selected_cycle_complete=bool(cycles),large_event_count=len(large),
        recurrence_years=intervals.tolist(),
        median_recurrence_years=float(np.median(intervals)) if len(intervals) else None,
        selected_event_count=len(selected_events),selected_small_events=sum(not e['large'] for e in selected_events),
        sample_depths=sample_depths,migration_segments=segments,slow_slip=slow_slip,
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
              'small events in selected cycle',r['selected_small_events'])
