#!/usr/bin/env python3
"""One command generates a named panel's data from the documented inputs."""
import argparse
import subprocess
from inventory import PANELS

parser=argparse.ArgumentParser();parser.add_argument('panel',choices=PANELS);parser.add_argument('--data',action='store_true');args=parser.parse_args()
cases=PANELS[args.panel]['cases']
if cases:
    subprocess.run(['python3','scripts/reproduce.py','--simulate-only','--cases']+cases,check=True)
    subprocess.run(['python3','scripts/analyze.py']+cases,check=True)
else:
    subprocess.run(['bash','scripts/bootstrap.sh'],check=True)
    subprocess.run(['build/valving','--laws','configs/baseline.cfg','data/laws'],check=True)
subprocess.run(['python3','scripts/prepare_panels.py',args.panel],check=True)
if not args.data:subprocess.run(['python3','scripts/plot.py',args.panel],check=True)
