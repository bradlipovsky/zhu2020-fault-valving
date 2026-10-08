#!/usr/bin/env python3
"""Inspect an explicitly unfinished calculation without creating report data."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from common import fields,YEAR,configuration
from plot import map_field

parser=argparse.ArgumentParser();parser.add_argument('case');args=parser.parse_args()
z,r=fields(args.case,complete=False)
keep=z<=25000;z=z[keep];t=np.asarray(r['time'])/YEAR
cfg=configuration(Path('data',args.case,'config.cfg'))
fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
for ax,field,title in zip(axs.flat,[1,2,3,5],['Slip speed','Effective stress','Permeability','Fluid flux']):
    map_field(fig,ax,t,z,np.asarray(r['fields'][:,field,:][:,keep]),field,cfg['influx'])
    ax.set_title(title);ax.set_xlabel('Time since independent initialization (yr)')
fig.suptitle('LIVE PREVIEW — {} at {:.2f} yr — incomplete, not report data'.format(args.case,t[-1]))
Path('figures/previews').mkdir(exist_ok=True,parents=True)
fig.savefig('figures/previews/'+args.case+'.png',dpi=130)
plt.close(fig)
