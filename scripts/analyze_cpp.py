"""Compare newly integrated fields with archived fields; retain compact products."""
from pathlib import Path
import json
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
    return dict(t=t,z=np.arange(meta['nz'])*meta['dz']/1000,meta=meta,
        slip=fields[0],slip_velocity=fields[1],effective_normal_stress=fields[2]/1e6,
        permeability=fields[3],flux=fields[4],psi=fields[5],kstar=fields[6])

def sample(d,key,z,t):
    points=np.stack(np.meshgrid(z,t,indexing='ij'),axis=-1)
    return RegularGridInterpolator((d['z'],d['t']),d[key],bounds_error=True)(points)

def compare(path):
    d=read_run(path);case=path.name.split('_')[1];ref=load(case)
    init=json.loads((ROOT/'inputs'/(case+'_restart.json')).read_text());offset=init['offset_years']
    end=float(d['t'][-1]);z=np.linspace(2,24.5,226);ts=np.linspace(0,end,101)
    metrics={'case':case,'end_years':end,'spacing_m':d['meta']['dz'],'offset_years':offset,
             'complete':bool(abs(end-d['meta']['years'])<1e-7)}
    for key in ['effective_normal_stress','permeability','slip_velocity','slip']:
        a=sample(d,key,z,ts);b=sample(ref,key,z,ts+offset)
        if key=='slip':b-=sample(ref,key,z,np.array([offset]))
        if key in ['permeability','slip_velocity']:
            a=np.log10(np.maximum(a,1e-30));b=np.log10(np.maximum(b,1e-30))
        metrics[key]={'rms':float(np.sqrt(np.mean((a-b)**2))),
                      'final_rms':float(np.sqrt(np.mean((a[:,-1]-b[:,-1])**2))),
                      'maximum_absolute':float(np.max(np.abs(a-b)))}
    # Compact original solver output every 0.05 yr, 100 m, retaining units.
    zp=np.linspace(0,25,251);tp=np.linspace(0,end,max(2,int(end/.05)+1))
    product={key:sample(d,key,zp,tp).astype(np.float32) for key in ['slip','slip_velocity','effective_normal_stress','permeability','flux']}
    np.savez_compressed(path/'plot_data.npz',z=zp,t=tp,**product)
    return d,metrics

def main():
    runs={};metrics={}
    for path in sorted((ROOT/'results').glob('restart_*')):
        if not path.is_dir() or not (path/'fields.bin').exists():continue
        runs[path.name],metrics[path.name]=compare(path)
    (ROOT/'results'/'cpp_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    coarse='restart_T1e8_n16384';fine='restart_T1e8_n32768'
    if coarse in runs and fine in runs:
        offset=metrics[coarse]['offset_years'];end=min(runs[coarse]['t'][-1],runs[fine]['t'][-1])
        ref=load('T1e8');ref['t']=ref['t']-offset
        fig,axes=plt.subplots(2,3,figsize=(12,7))
        for col,(d,name) in enumerate([(ref,'Archived'),(runs[coarse],'C++: 35.34 m'),(runs[fine],'C++: 17.67 m')]):
            for row,key in enumerate(['slip_velocity','effective_normal_stress']):
                mesh(axes[row,col],d,key,(0,end),letter='abcdef'[3*row+col])
                axes[row,col].set_title(name+': '+('velocity' if row==0 else 'effective stress'),fontsize=9)
        save(fig,'cpp_comparison','New C++ integration from approximate archived restart states')
        z=np.linspace(2,24.5,226);tt=np.linspace(0,end,101);diff={}
        for key in ['slip','effective_normal_stress','permeability','slip_velocity']:
            a=sample(runs[coarse],key,z,tt);b=sample(runs[fine],key,z,tt)
            if key in ['permeability','slip_velocity']:a=np.log10(np.maximum(a,1e-30));b=np.log10(np.maximum(b,1e-30))
            diff[key]=float(np.sqrt(np.mean((a-b)**2)))
        (ROOT/'results'/'spatial_comparison.json').write_text(json.dumps({'common_end_years':float(end),'coarse_fine_rms':diff},indent=2)+'\n')
    rows=[]
    for name,m in metrics.items():
        T=float(m['case'][1:])/YEAR
        rows.append(f"{T:.3g} & {m['spacing_m']:.2f} & {m['end_years']:.2f} & {m['effective_normal_stress']['rms']:.3g} & {m['permeability']['rms']:.3g} & {m['slip_velocity']['rms']:.3g} \\\\")
    text=r'''\begin{table}[h]
\centering\small
\begin{tabular}{rrrrrr}\toprule
$T$ (yr) & $\Delta z$ (m) & Duration (yr) & Stress RMS (MPa) & $\log_{10}k$ RMS & $\log_{10}V$ RMS\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}
\caption{New C++ runs compared with the archived fields over 2--24.5 km depth and each completed time window. Logarithmic errors are in decades. Each run starts near four years into its archive. RMS values include the full sampled time window, not only its endpoint.}
\end{table}
'''
    if fine in metrics:
        m=metrics[fine]
        text+=f"The finer featured-case run covers {m['end_years']:.2f} years after the approximate restart. "
        text+=f"Its effective-stress RMS difference from the archive is {m['effective_normal_stress']['rms']:.3g} MPa, "
        text+=f"and its logarithmic slip-rate RMS difference is {m['slip_velocity']['rms']:.3g} decades. "
        text+="These are comparisons with a reconstructed initial condition, rather than errors relative to an exact solution. "
    if coarse in runs and fine in runs:
        text+=f"Over the common {end:.2f}-year interval, the two finest featured-case meshes differ by {diff['effective_normal_stress']:.3g} MPa RMS in effective stress. "
        text+="Two-grid agreement alone does not establish asymptotic convergence or resolve the later swarm sequence.\n"
    (ROOT/'report'/'simulation_results.tex').write_text(text)
    print(json.dumps(metrics,indent=2))

if __name__=='__main__':main()
