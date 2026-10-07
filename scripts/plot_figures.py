"""Recompute constitutive plots and redraw all published panels from OSF data.

Every archived figure is labelled. These plots do not establish independent
reproduction of the earthquake sequence by the new C++ solver.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
YEAR=31536000
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
    'savefig.dpi':180,'pdf.fonttype':42,'figure.constrained_layout.use':True})
CMAP=LinearSegmentedColormap.from_list('slip',[(0,'#182d54'),(.24,'#274879'),(.374,'#69b7d2'),
    (.4,'#f7fbfc'),(.5,'#f9c3a7'),(.7,'#d94832'),(1,'#76091b')])

def load(name):
    return dict(np.load(ROOT/'build'/'archive'/(name+'_si.npz')))

def save(fig,name,source='Authors’ archived simulation; independently redrawn'):
    fig.suptitle(source,fontsize=9,color='.35')
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/(name+'.png'))
    plt.close(fig)

def label(ax,letter,title):
    ax.set_title(letter+'  '+title,loc='left',fontsize=10)
    ax.set_ylabel('Depth (km)')

def mesh(ax,d,key,window=None,depth=(0,25),steps=False,letter='',bar=True):
    t=d['t']; j=np.ones(len(t),dtype=bool) if window is None else (t>=window[0])&(t<=window[1])
    k=(d['z']>=depth[0])&(d['z']<=depth[1]); z=d['z'][k]
    x=np.arange(j.sum()) if steps else t[j]
    a=d[key][k][:,j]
    opts={}
    if key=='slip_velocity':
        a=np.log10(np.maximum(a,1e-30)); title='Slip velocity';unit=r'$\log_{10}[V/(\mathrm{m\,s^{-1}})]$'
        opts=dict(vmin=-15,vmax=0,cmap=CMAP)
    elif key=='effective_normal_stress':
        title='Effective normal stress';unit='MPa';opts=dict(cmap='viridis')
    elif key=='permeability':
        a=np.log10(a);title='Permeability';unit=r'$\log_{10}[k/\mathrm{m^2}]$';opts=dict(vmin=-19,vmax=-15,cmap='viridis')
    elif key=='flux':
        a=np.ma.masked_less_equal(a,0);a=np.ma.log10(a);title='Upward fluid flux';unit=r'$\log_{10}[q/(\mathrm{m\,s^{-1}})]$'
        opts=dict(vmin=-11,vmax=-7,cmap='viridis')
    im=ax.pcolormesh(x,z,a,shading='auto',rasterized=True,**opts)
    ax.set_ylim(depth[::-1]);ax.set_xlim(x[0],x[-1]);ax.set_xlabel('Saved sample index' if steps else 'Time (yr)')
    label(ax,letter,title)
    if bar:plt.colorbar(im,ax=ax,label=unit,pad=.02,shrink=.9)
    return im

def slip_profiles(ax,d,window,interval,letter):
    t=d['t'];sel=(t>=window[0])&(t<=window[1]);inds=np.flatnonzero(sel)
    slip=d['slip']-d['slip'][:,inds[0]:inds[0]+1]
    for tt in np.arange(window[0],window[1],interval):
        j=np.argmin(abs(t-tt));ax.plot(slip[:,j],d['z'],color='#2468ad',lw=.5,alpha=.8)
    # One-second snapshots, interpolated only within each archived seismic burst.
    seismic=(d['slip_velocity'].max(axis=0)>1e-3)&sel
    starts=np.flatnonzero(seismic & ~np.r_[False,seismic[:-1]])
    ends=np.flatnonzero(seismic & ~np.r_[seismic[1:],False])
    for st,en in zip(starts,ends):
        ts=np.arange(t[st]*YEAR,t[en]*YEAR+1,1)/YEAR
        for tt in ts:
            j=np.searchsorted(t,tt);j=np.clip(j,1,len(t)-1)
            w=(tt-t[j-1])/(t[j]-t[j-1]);ax.plot((1-w)*slip[:,j-1]+w*slip[:,j],d['z'],color='#be302b',lw=.45,alpha=.65)
    ax.set_ylim(25,0);ax.set_xlabel('Slip since window start (m)');label(ax,letter,'Cumulative slip')

def sequence(name,ref,case,end,interval=1.5):
    fig,axs=plt.subplots(2,3,figsize=(12,7),gridspec_kw={'width_ratios':[1.15,1,.8]})
    for row,d in enumerate([load(ref),load(case)]):
        stop=min(end if row else 1e6,d['t'][-1])
        if ref=='T1e8-ref' and row==0:stop=56
        mesh(axs[row,0],d,'slip_velocity',(0,stop),letter='ac'[row])
        mesh(axs[row,1],d,'slip_velocity',(0,stop),steps=True,letter='')
        slip_profiles(axs[row,2],d,(0,stop),interval,'bd'[row])
        axs[row,0].text(.03,.94,'Fixed pressure' if row==0 else 'Fault valving',transform=axs[row,0].transAxes,color='white')
    save(fig,name)

def hydraulic(name,case,end,interval=1.5):
    d=load(case);fig,axs=plt.subplots(3,2,figsize=(10,10),gridspec_kw={'width_ratios':[1.6,1]})
    for row,key in enumerate(['effective_normal_stress','permeability','flux']):
        mesh(axs[row,0],d,key,(0,end),letter='ace'[row])
        ax=axs[row,1]
        times=np.arange(0,end,interval)
        for tt in times:
            j=np.argmin(abs(d['t']-tt));ax.plot(d[key][:,j],d['z'],color=plt.cm.viridis(tt/end),lw=.7)
        z=d['z'];p,k,eff=steady(z*1000,float(case[1:]),3e-9)
        ss={'effective_normal_stress':eff/1e6,'permeability':k,'flux':np.full_like(z,3e-9)}[key]
        ax.plot(ss,z,'k--',lw=1,label='Steady sliding')
        if key!='effective_normal_stress':ax.set_xscale('log')
        ax.set_xlabel({'effective_normal_stress':'Effective normal stress (MPa)','permeability':r'Permeability (m$^2$)','flux':r'Upward flux (m s$^{-1}$)'}[key])
        ax.set_ylim(25,0);label(ax,'bdf'[row],'Depth profiles');ax.legend(fontsize=8)
    save(fig,name)

def four_fields(name,case,window=None,depth=(0,25)):
    d=load(case);fig,axs=plt.subplots(2,2,figsize=(10,7))
    mesh(axs[0,0],d,'slip_velocity',window,depth,letter='a')
    mesh(axs[0,1],d,'slip_velocity',window,depth,steps=True,letter='b')
    mesh(axs[1,0],d,'effective_normal_stress',window,depth,letter='c')
    mesh(axs[1,1],d,'permeability',window,depth,letter='d')
    save(fig,name)

def steady(z,T=1e8,q=3e-9):
    ks=(1e-15*1e-9+1e-19/T)/(1e-9+1/T)
    def fun(zz,p):
        k=1e-19+(ks-1e-19)*np.exp(-(22050*zz-p[0])/30e6)
        return [9800+1e-4*q/k]
    p=solve_ivp(fun,(0,float(max(z))),[0.],t_eval=z,rtol=1e-10,atol=1e-5,max_step=50).y[0]
    eff=22050*z-p;k=1e-19+(ks-1e-19)*np.exp(-eff/30e6)
    return p,k,eff

def figure1():
    fig=plt.figure(figsize=(11,10));gs=fig.add_gridspec(2,2)
    sub=gs[0,0].subgridspec(1,2,width_ratios=[1.2,1]);ax=fig.add_subplot(sub[0,0])
    ax.add_patch(plt.Rectangle((0,0),1,25,color='#eef4f7'))
    ax.plot([0,0],[0,25],lw=6,color='#65acc6');ax.text(.1,4,'Elastic half-space\n'+r'$\mu=32.4$ GPa')
    ax.text(.1,10,'Rate-and-state friction\nPermeability evolution\nUpward Darcy flow')
    ax.annotate('',xy=(0,3),xytext=(0,22),arrowprops=dict(arrowstyle='->',lw=2,color='#2468ad'))
    ax.text(.03,24,r'$q=q_0$ at bottom');ax.text(.03,1,r'$p=0$ at surface');ax.text(.45,19,r'$u=V_pt/2$')
    ax.set_ylim(25,0);ax.set_xlim(-.05,1);ax.set_xticks([]);label(ax,'a','Model geometry')
    ax=fig.add_subplot(sub[0,1]);z=np.linspace(0,30,1000)
    a=np.interp(z,[0,14.9,27.6,60],[.0105,.03,.07,.173]);b=np.interp(z,[0,13.6,14.9,27.6,60],[.02,.0378,.0356,.0375,.0424])
    ax.plot(a,z,label='$a$');ax.plot(a-b,z,label='$a-b$');ax.axvline(0,color='.5',ls='--',lw=.7)
    ax.set_ylim(25,0);ax.set_xlabel('Friction parameters');ax.legend();ax.set_ylabel('Depth (km)')
    ax=fig.add_subplot(gs[0,1]);eff=np.linspace(0,400,500)
    for ks in [1e-15,1e-16,1e-17,1e-18]:ax.semilogy(eff,1e-19+(ks-1e-19)*np.exp(-eff/30),label=r'$k^*=10^{%d}$'%round(np.log10(ks)))
    ax.axhline(1e-19,ls='--',color='.4');ax.set_xlabel('Effective normal stress (MPa)');ax.set_ylabel(r'Permeability (m$^2$)');ax.set_title('b  Stress-dependent permeability',loc='left');ax.legend()
    sub=gs[1,0].subgridspec(2,1);ax=fig.add_subplot(sub[0,0]);t=np.linspace(0,400,1000)
    for T in [1e8,1e9,1e10]:ax.semilogy(t,1e-19+(1e-15-1e-19)*np.exp(-t*YEAR/T),label=f'$T={T/YEAR:.3g}$ yr')
    ax.set_xlabel('Time without slip (yr)');ax.set_ylabel(r'$k^*$ (m$^2$)');ax.set_title('c  Healing and slip enhancement',loc='left');ax.legend(fontsize=8)
    ax=fig.add_subplot(sub[1,0]);sl=np.linspace(0,4,1000);ax.semilogy(sl,1e-15+(1e-19-1e-15)*np.exp(-sl));ax.set_xlabel('Slip (m)');ax.set_ylabel(r'$k^*$ (m$^2$)');ax.text(2,1e-17,'Healing neglected during slip')
    ax=fig.add_subplot(gs[1,1]);z=np.linspace(0,25000,1000);p,k,eff=steady(z)
    for x,lab,ls in [(p/1e6,'Pore pressure','-'),(eff/1e6,'Effective stress','-'),(22050*z/1e6,'Normal stress','--'),(9800*z/1e6,'Hydrostatic','--'),(26460*z/1e6,'Lithostatic','--')]:ax.plot(x,z/1000,ls,label=lab)
    ax.set_ylim(25,0);ax.set_xlim(0,400);ax.set_xlabel('Pressure or stress (MPa)');label(ax,'d','Steady upward flow');ax.legend(fontsize=8,loc='lower right')
    top=ax.twiny();top.plot(k,z/1000,'k');top.set_xscale('log');top.set_xlabel(r'Permeability (m$^2$)')
    save(fig,'figure1','Independent evaluation of published equations and original input profiles')

def front(d):
    # Shallowest V=0.1*Vp crossing below 2 km. This operational definition
    # captures the leading pulse even when the deeper fault relocks.
    threshold=1e-10
    k=np.flatnonzero((d['z']>=2)&(d['z']<=25));v=d['slip_velocity'][k];z=d['z'][k]
    result=np.empty(len(d['t']))
    for j in range(len(result)):
        above=np.flatnonzero(v[:,j]>=threshold)
        if len(above)==0:result[j]=np.nan;continue
        if above[0]==0:result[j]=z[0];continue
        i=above[0]-1
        w=(np.log10(threshold)-np.log10(max(v[i,j],1e-30)))/(np.log10(v[i+1,j])-np.log10(max(v[i,j],1e-30)))
        result[j]=z[i]+w*(z[i+1]-z[i])
    return result

def figure6():
    fig,ax=plt.subplots(figsize=(7,7));stats={}
    for case,color,offset,fit in [('T1e7','#368ed1',12,(12,37)),('T1e8','#eb8634',6,(6,18)),('T1e9','#57bad3',20,(20,45)),('T1e10','#d65cc8',25,(25,50))]:
        d=load(case);z=front(d);t=d['t'];sel=(t>=offset)&(t<=offset+25)
        ax.plot(t[sel]-offset,z[sel],color=color,label=f'$T={float(case[1:])/YEAR:.3g}$ yr',lw=1.1)
        # Declared windows isolate advancing fronts after postseismic retreat.
        # The paper does not specify a front threshold or the Fig. 6 time offsets.
        sel=(t>=fit[0])&(t<=fit[1])&np.isfinite(z)&(d['slip_velocity'].max(axis=0)<1e-5)
        coeff=np.polyfit(t[sel],z[sel],1)
        stats[case]={'front_rate_km_per_year':float(-coeff[0]),'threshold_m_s':1e-10,'plot_offset_years':offset,'fit_window_years':fit}
    ax.set_ylim(25,0);ax.set_xlim(0,25);ax.set_xlabel('Time since declared window start (yr)');ax.set_ylabel('Leading-front depth (km)');ax.legend()
    ax.text(2,3,'Locked');ax.text(2,23,'Creeping');save(fig,'figure6')
    (ROOT/'results'/'front_rates.json').write_text(json.dumps(stats,indent=2)+'\n')

def supplements():
    d=load('T1e8');fig,axs=plt.subplots(4,1,figsize=(10,9),sharex=True)
    for row,key in enumerate(['slip_velocity','effective_normal_stress','permeability','flux']):
        for depth in [5,10,15,20]:
            j=np.argmin(abs(d['z']-depth));axs[row].plot(d['t'],d[key][j],label=f'{depth} km',lw=.7)
        if row!=1:axs[row].set_yscale('log')
        axs[row].set_ylabel([r'$V$ (m s$^{-1}$)','Effective stress\n(MPa)',r'$k$ (m$^2$)',r'$q$ (m s$^{-1}$)'][row]);axs[row].set_xlim(0,35)
        axs[row].set_title('abcd'[row],loc='left')
    axs[0].legend(ncol=4);axs[-1].set_xlabel('Time (yr)');save(fig,'supplement1')
    fig,axs=plt.subplots(2,2,figsize=(10,7))
    for row,case in enumerate(['T1e8-sup1','T1e8-sup2']):
        d=load(case);mesh(axs[row,0],d,'slip_velocity',letter='ac'[row]);mesh(axs[row,1],d,'slip_velocity',steps=True,letter='bd'[row])
    save(fig,'supplement2')
    sequence('supplement3','T1e9-ref','T1e9',140,4)
    hydraulic('supplement4','T1e9',130,4)
    sequence('supplement5','T1e10-ref','T1e10',280,4)
    hydraulic('supplement6','T1e10',280,4)

if __name__=='__main__':
    figure1();print('Figure 1',flush=True)
    sequence('figure2','T1e8-ref','T1e8',35)
    hydraulic('figure3','T1e8',35)
    four_fields('figure4','T1e8',(30.5,31.8),(2,10))
    four_fields('figure5','T1e7',(0,48))
    figure6();print('Main figures complete',flush=True)
    supplements();print('Supplementary figures complete',flush=True)
