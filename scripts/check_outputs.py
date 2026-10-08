"""Check completion, file integrity, and elementary field constraints."""
from pathlib import Path
import json
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
metrics=json.loads((ROOT/'results'/'cpp_metrics.json').read_text())
assert len(metrics)==6, 'Expected three featured-case meshes and three other healing times'
for name,record in metrics.items():
    assert record['complete'], 'Incomplete run: '+name
    p=ROOT/'results'/name
    m=json.loads((p/'metadata.json').read_text())
    width=8*(1+7*m['nz'])
    assert (p/'fields.bin').stat().st_size%width==0, 'Truncated binary record: '+name
    assert (p/'final_state.bin').stat().st_size==4*m['n']*8, 'Invalid checkpoint: '+name
    with np.load(p/'plot_data.npz') as d:
        assert abs(d['t'][-1]-m['years'])<1e-7
        assert np.all(np.diff(d['t'])>0)
        for key in d.files:assert np.isfinite(d[key]).all(), (name,key)
        assert np.all(d['slip_velocity']>0)
        assert np.all(d['permeability']>0)
        assert np.all(d['effective_normal_stress']>=1-1e-5)
for prefix in ['figure','supplement']:
    for n in range(1,7):
        p=ROOT/'figures'/f'{prefix}{n}.pdf'
        subprocess.run(['pdfinfo',str(p)],check=True,stdout=subprocess.DEVNULL)
for name in ['cpp_cycle','cpp_comparison','cpp_healing_cases','cpp_time_accuracy','cpp_hydraulics']:
    subprocess.run(['pdfinfo',str(ROOT/'figures'/(name+'.pdf'))],check=True,stdout=subprocess.DEVNULL)
subprocess.run(['pdfinfo',str(ROOT/'report'/'reproduction.pdf')],check=True,stdout=subprocess.DEVNULL)
assert 'PASS' in (ROOT/'results'/'verification.txt').read_text()
assert json.loads((ROOT/'results'/'seismic_time_comparison.json').read_text())['complete']
assert json.loads((ROOT/'results'/'integrator_comparison.json').read_text())['complete']
assert json.loads((ROOT/'results'/'threading_comparison.json').read_text())['bitwise_identical']
log=(ROOT/'report'/'reproduction.log').read_text()
assert 'undefined references' not in log
assert 'Overfull' not in log
print('PASS: six completed runs, finite fields, complete checkpoints, twelve reconstructed figures, five independent comparisons, compiled report, and solver verification')
