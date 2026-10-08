"""Compare newly integrated fields with archived fields; retain compact products."""
from pathlib import Path
import json
import re
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_figures import load, mesh, save

ROOT=Path(__file__).resolve().parents[1]
YEAR=31536000

def read_run(path):
    meta=json.loads((path/'metadata.json').read_text())
    width=1+7*meta['nz']
    raw=np.fromfile(path/'fields.bin')
    raw=raw[:len(raw)//width*width].reshape(-1,width)
    t=raw[:,0]/YEAR
    fields=raw[:,1:].reshape(-1,7,meta['nz']).transpose(1,2,0)
    if not np.isfinite(raw).all():raise ValueError('Nonfinite output in '+str(path))
    return dict(t=t,z=np.arange(meta['nz'])*meta['dz']/1000,
        flux_z=(np.arange(meta['nz'])+.5)*meta['dz']/1000,meta=meta,
        slip=fields[0],slip_velocity=fields[1],effective_normal_stress=fields[2]/1e6,
        permeability=fields[3],flux=fields[4],psi=fields[5],kstar=fields[6])

def sample(d,key,z,t):
    points=np.stack(np.meshgrid(z,t,indexing='ij'),axis=-1)
    grid=d.get('flux_z',d['z']) if key=='flux' else d['z']
    return RegularGridInterpolator((grid,d['t']),d[key],bounds_error=True)(points)

def differences(a,b,end,start=0):
    """RMS differences on a common physical-time/depth grid, not output index."""
    z=np.linspace(2,24.5,226);t=np.linspace(start,end,101);result={}
    for key in ['slip','effective_normal_stress','permeability','slip_velocity']:
        x=sample(a,key,z,t);y=sample(b,key,z,t)
        if key in ['permeability','slip_velocity']:
            x=np.log10(np.maximum(x,1e-30));y=np.log10(np.maximum(y,1e-30))
        result[key]={'rms':float(np.sqrt(np.mean((x-y)**2))),
                     'max_abs':float(np.max(np.abs(x-y)))}
    return result

def events(d,start,end):
    """Contiguous saved samples above 1 mm/s within 2--25 km define a burst."""
    use=(d['t']>=start)&(d['t']<=end)
    t=d['t'][use];z=d['z'][(d['z']>=2)&(d['z']<=25)]
    v=d['slip_velocity'][(d['z']>=2)&(d['z']<=25)][:,use]
    peak=v.max(axis=0);active=peak>1e-3
    starts=np.flatnonzero(active&~np.r_[False,active[:-1]])
    stops=np.flatnonzero(active&~np.r_[active[1:],False])
    result=[]
    for first,last in zip(starts,stops):
        j=first+int(np.argmax(peak[first:last+1]))
        result.append({'onset_year':float(t[first]),'end_year':float(t[last]),
            'onset_bracket_years':[float(t[max(0,first-1)]),float(t[first])],
            'peak_year':float(t[j]),'peak_velocity_m_s':float(peak[j]),
            'peak_depth_km':float(z[np.argmax(v[:,j])])})
    return result

def compare(path):
    d=read_run(path);case=path.name.split('_')[1];ref=load(case)
    init=json.loads((ROOT/'inputs'/(case+'_restart.json')).read_text());offset=init['offset_years']
    end=float(d['t'][-1]);comparison_end=min(end,float(ref['t'][-1])-offset)
    z=np.linspace(2,24.5,226);ts=np.linspace(0,comparison_end,101)
    metrics={'case':case,'end_years':end,'spacing_m':d['meta']['dz'],'offset_years':offset,
             'comparison_years':comparison_end,
             'complete':bool(abs(end-d['meta']['years'])<1e-7)}
    for key in ['effective_normal_stress','permeability','slip_velocity','slip']:
        a=sample(d,key,z,ts);b=sample(ref,key,z,np.minimum(ts+offset,ref['t'][-1]))
        if key=='slip':b-=sample(ref,key,z,np.array([offset]))
        if key in ['permeability','slip_velocity']:
            a=np.log10(np.maximum(a,1e-30));b=np.log10(np.maximum(b,1e-30))
        metrics[key]={'rms':float(np.sqrt(np.mean((a-b)**2))),
                      'final_rms':float(np.sqrt(np.mean((a[:,-1]-b[:,-1])**2))),
                      'maximum_absolute':float(np.max(np.abs(a-b)))}
    # Compact original solver output every 0.05 yr, 100 m, retaining units.
    zp=np.linspace(0,25,251)
    regular=np.linspace(0,end,max(2,int(end/.05)+1))
    seismic=d['t'][d['slip_velocity'].max(axis=0)>1e-4]
    tp=np.unique(np.r_[regular,seismic])
    product={key:sample(d,key,zp,tp).astype(np.float32) for key in ['slip','slip_velocity','effective_normal_stress','permeability']}
    flux_z=zp+d['meta']['dz']/2000
    product['flux']=sample(d,'flux',flux_z,tp).astype(np.float32)
    np.savez_compressed(path/'plot_data.npz',z=zp,flux_z=flux_z,t=tp,**product)
    return d,metrics

def main():
    runs={};metrics={}
    for path in sorted((ROOT/'results').glob('restart_*')):
        if not path.is_dir() or not (path/'fields.bin').exists():continue
        runs[path.name],metrics[path.name]=compare(path)
    (ROOT/'results'/'cpp_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    catalog={'definition':'Contiguous saved samples with max(V)>1e-3 m/s in 2--25 km depth; times measured from each original archive start.', 'runs':{}}
    for name,d in runs.items():
        offset=metrics[name]['offset_years'];end=float(d['t'][-1]);ref=load(metrics[name]['case'])
        shifted=dict(d,t=d['t']+offset)
        catalog['runs'][name]={'archive':events(ref,offset,offset+end),
                               'new_solver':events(shifted,offset,offset+end)}
    (ROOT/'results'/'event_catalog.json').write_text(json.dumps(catalog,indent=2)+'\n')
    names=['restart_T1e8_n'+str(n) for n in [8192,16384,32768]]
    timing=''
    if all(name in catalog['runs'] and catalog['runs'][name]['new_solver'] for name in names):
        reference=catalog['runs'][names[0]]['archive'][0]['onset_year']
        timing=r'''\begin{table}[h]\centering\small
\begin{tabular}{rrr}\toprule
$\Delta z$ (m) & First crossing (archive-relative yr) & Difference from archive (days)\\\midrule
'''
        for name in names:
            onset=catalog['runs'][name]['new_solver'][0]['onset_year']
            timing+=f"{metrics[name]['spacing_m']:.2f} & {onset:.6f} & {(onset-reference)*365:.2f} \\\\\n"
        timing+=r'''\bottomrule\end{tabular}
\caption{First sampled crossing of 1 mm\,s$^{-1}$ between 2 and 25 km depth in the featured-case spatial refinements. Times retain the archive's origin; its corresponding crossing occurs at '''+f'{reference:.6f}'+r''' yr. A negative difference denotes an earlier event. The comparison includes uncertainty in the reconstructed initial state and does not measure error relative to an exact solution.}
\end{table}
'''
    (ROOT/'report'/'event_timing.tex').write_text(timing)
    baseline=ROOT/'results'/'seismic_time_baseline'
    refinement=ROOT/'results'/'seismic_time_refinement'
    if (baseline/'fields.bin').exists() and (refinement/'fields.bin').exists():
        a=read_run(baseline);b=read_run(refinement);end=min(a['t'][-1],b['t'][-1])
        comparison={'start_year':12,'end_year':float(end),
                    'complete':bool(abs(end-14.7)<1e-7),
                    'field_differences':differences(a,b,end,12)}
        for name,d in [('baseline',a),('refinement',b)]:
            comparison[name+'_events']=events(d,12,end)
        if comparison['baseline_events'] and comparison['refinement_events']:
            onset=comparison['baseline_events'][0]['onset_year']
            comparison['first_onset_shift_s']=(comparison['refinement_events'][0]['onset_year']-onset)*YEAR
            fig,ax=plt.subplots(figsize=(8,4))
            for d,label in [(a,r'$r_{tol}=10^{-4}$'),(b,r'$r_{tol}=2.5\times10^{-5}$')]:
                t=(d['t']-onset)*YEAR;v=d['slip_velocity'][(d['z']>=2)&(d['z']<=25)].max(axis=0)
                select=(t>=-300)&(t<=600)
                ax.semilogy(t[select],v[select],label=label,lw=1)
            ax.set_xlabel('Time relative to baseline first threshold crossing (s)')
            ax.set_ylabel(r'Maximum velocity at 2--25 km (m s$^{-1}$)')
            ax.set_ylim(1e-5,10);ax.legend()
            save(fig,'cpp_time_accuracy','First seismic event: identical initial checkpoint, two time tolerances')
        (ROOT/'results'/'seismic_time_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    candidate=ROOT/'results'/'time_stepping_candidate'
    if (candidate/'fields.bin').exists() and (baseline/'fields.bin').exists():
        old=read_run(baseline);new=read_run(candidate);refined=read_run(refinement)
        end=min(old['t'][-1],new['t'][-1],refined['t'][-1])
        e1=events(old,12,end)[0];e2=events(new,12,end)[0];e3=events(refined,12,end)[0]
        wall=lambda name:float(re.findall(r'wall_s=([0-9.]+)',(ROOT/'results'/(name+'.log')).read_text())[-1])
        result={'complete':bool(abs(end-14.7)<1e-7),'old_wall_s':wall('seismic_time_baseline'),
            'new_wall_s':wall('time_stepping_candidate'),
            'first_onset_shift_s':(e2['onset_year']-e1['onset_year'])*YEAR,
            'peak_velocity_old_m_s':e1['peak_velocity_m_s'],
            'peak_velocity_new_m_s':e2['peak_velocity_m_s'],
            'peak_velocity_refined_old_m_s':e3['peak_velocity_m_s'],
            'old_new':differences(old,new,end,12),'new_refined_old':differences(new,refined,end,12)}
        (ROOT/'results'/'integrator_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    coarse='restart_T1e8_n16384';fine='restart_T1e8_n32768'
    if coarse in runs and fine in runs:
        offset=metrics[coarse]['offset_years'];common_end=min(12,runs[coarse]['t'][-1],runs[fine]['t'][-1])
        ref=load('T1e8');ref['t']=ref['t']-offset
        fig,axes=plt.subplots(2,3,figsize=(12,7))
        for col,(d,name) in enumerate([(ref,'Archived'),(runs[coarse],'C++: 35.34 m'),(runs[fine],'C++: 17.67 m')]):
            for row,key in enumerate(['slip_velocity','effective_normal_stress']):
                im=mesh(axes[row,col],d,key,(0,common_end),letter='abcdef'[3*row+col])
                if row==1:im.set_clim(0,55)
                axes[row,col].set_title('abcdef'[3*row+col]+'  '+name+': '+('velocity' if row==0 else 'effective stress'),loc='left',fontsize=9)
        save(fig,'cpp_comparison','New C++ integration from approximate archived restart states')
        diff=differences(runs[coarse],runs[fine],common_end)
        spatial={'common_end_years':float(common_end),'medium_fine':diff}
        if 'restart_T1e8_n8192' in runs:
            spatial['coarse_medium']=differences(runs['restart_T1e8_n8192'],runs[coarse],common_end)
        (ROOT/'results'/'spatial_comparison.json').write_text(json.dumps(spatial,indent=2)+'\n')
        refinement=ROOT/'results'/'time_refinement'
        if (refinement/'fields.bin').exists():
            refined=read_run(refinement)
            temporal=differences(runs[coarse],refined,float(refined['t'][-1]))
            (ROOT/'results'/'temporal_comparison.json').write_text(json.dumps(temporal,indent=2)+'\n')
    if coarse in runs:
        ref=load('T1e8');ref['t']-=metrics[coarse]['offset_years'];end=runs[coarse]['t'][-1]
        fig,axes=plt.subplots(2,3,figsize=(12,8))
        for row,(d,name) in enumerate([(ref,'Archived'),(runs[coarse],'New C++')]):
            for col,key in enumerate(['slip_velocity','slip_velocity','effective_normal_stress']):
                im=mesh(axes[row,col],d,key,(0,end),steps=(col==1),letter='abcdef'[3*row+col])
                if col==2:im.set_clim(0,55)
                quantity=['velocity','velocity (sample index)','effective stress'][col]
                axes[row,col].set_title('abcdef'[3*row+col]+'  '+name+': '+quantity,loc='left',fontsize=9)
        save(fig,'cpp_cycle','New C++ earthquake sequence at 35.34 m spacing; approximate archived restart')
        fig,axes=plt.subplots(3,2,figsize=(11,10))
        for col,(d,name) in enumerate([(ref,'Archived'),(runs[coarse],'New C++ integration')]):
            for row,key in enumerate(['effective_normal_stress','permeability','flux']):
                im=mesh(axes[row,col],d,key,(0,end),letter='abcdef'[2*row+col])
                if row==0:im.set_clim(0,55)
                quantity=['effective stress','permeability','upward flux'][row]
                axes[row,col].set_title('abcdef'[2*row+col]+'  '+name+': '+quantity,loc='left',fontsize=9)
        save(fig,'cpp_hydraulics','Independent C++ hydraulic evolution at 35.34 m spacing')
    cases=['T1e7','T1e9','T1e10']
    if all('restart_'+case+'_n16384' in runs for case in cases):
        fig,axes=plt.subplots(3,2,figsize=(11,10))
        for row,case in enumerate(cases):
            name='restart_'+case+'_n16384';d=runs[name];ref=load(case)
            ref['t']-=metrics[name]['offset_years']
            for col,(data,label) in enumerate([(ref,'Archived'),(d,'New C++ integration')]):
                mesh(axes[row,col],data,'slip_velocity',(0,metrics[name]['comparison_years']),letter='abcdef'[2*row+col])
                axes[row,col].set_title('abcdef'[2*row+col]+'  '+label+f": T = {float(case[1:])/YEAR:.3g} yr",loc='left',fontsize=9)
        save(fig,'cpp_healing_cases','Independent coupled integrations at other healing times; 35.34 m spacing')
    rows=[]
    for name,m in sorted(metrics.items(),key=lambda pair:(float(pair[1]['case'][1:]),-pair[1]['spacing_m'])):
        T=float(m['case'][1:])/YEAR
        rows.append(f"{T:.3g} & {m['spacing_m']:.2f} & {m['end_years']:.2f} & {m['effective_normal_stress']['rms']:.3g} & {m['permeability']['rms']:.3g} & {m['slip_velocity']['rms']:.3g} \\\\")
    text=r'''\begin{table}[h]
\centering\small
\begin{tabular}{rrrrrr}\toprule
$T$ (yr) & $\Delta z$ (m) & Duration (yr) & Stress RMS (MPa) & $\log_{10}k$ RMS & $\log_{10}V$ RMS\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}
\caption{New C++ runs compared with the archived fields over 2--24.5 km depth and the available overlap in time. Duration gives the full new integration; the $T=0.317$ yr comparison ends 44.82 yr after its restart because the archive ends there. Logarithmic errors are in decades, with velocities below $10^{-30}$ m\,s$^{-1}$ set to that value for this diagnostic only. Each run starts near four years into its archive. RMS values use 101 uniformly spaced physical times over the overlap and 226 depths at 100~m spacing; this sampling does not resolve every seismic peak.}
\end{table}
'''
    if fine in metrics:
        m=metrics[fine]
        text+=f"The finer featured-case run covers {m['end_years']:.2f} years after the approximate restart. "
        text+=f"Its effective-stress RMS difference from the archive is {m['effective_normal_stress']['rms']:.3g} MPa, "
        text+=f"and its logarithmic slip-rate RMS difference is {m['slip_velocity']['rms']:.3g} decades. "
        text+="These are comparisons with a reconstructed initial condition, rather than errors relative to an exact solution. "
    if coarse in runs and fine in runs:
        text+=f"Over the first {common_end:.2f} years before the earthquake, the two finest featured-case meshes differ by {diff['effective_normal_stress']['rms']:.3g} MPa RMS in effective stress. "
        text+="This smooth-field comparison does not establish convergence of the later seismic sequence.\n"
    (ROOT/'report'/'simulation_results.tex').write_text(text)
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
