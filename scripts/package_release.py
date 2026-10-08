"""Package a standalone LaTeX bundle and the newly generated full field records."""
from pathlib import Path
import hashlib
import json
import shutil
import tarfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build'/'release'
OUT.mkdir(parents=True,exist_ok=True)
metrics=json.loads((ROOT/'results'/'cpp_metrics.json').read_text())
assert all(m['complete'] for m in metrics.values()), 'Finish the simulations before packaging'

def archive(name,paths):
    with tarfile.open(str(OUT/name),'w:gz',compresslevel=1) as tar:
        for path in sorted(set(paths)):
            tar.add(str(path),arcname=str(path.relative_to(ROOT)),recursive=False)

latex=list((ROOT/'report').glob('*.tex'))+[ROOT/'README.md',ROOT/'LICENSE']
latex+=list((ROOT/'figures').glob('*.pdf'))
archive('latex-source-with-figures.tar.gz',latex)
raw=[]
for name in list(metrics)+['time_refinement','seismic_time_baseline','seismic_time_refinement','time_stepping_candidate','threading_serial','threading_parallel']:
    path=ROOT/'results'/name
    for filename in ['fields.bin','final_state.bin','phase1_state.bin','phase2_state.bin','phase3_state.bin','checkpoint_state.bin','checkpoint.json','metadata.json','history.csv','plot_data.npz']:
        if (path/filename).exists():raw.append(path/filename)
    raw.append(path.with_suffix('.log'))
raw+=list((ROOT/'inputs').glob('*'))
raw+=list((ROOT/'results').glob('*.json'))
raw+=list((ROOT/'results').glob('*.txt'))
raw+=[ROOT/'scripts'/'analyze_cpp.py',ROOT/'scripts'/'plot_figures.py',ROOT/'scripts'/'check_threading.py',ROOT/'README.md',ROOT/'LICENSE']
archive('new-simulation-fields.tar.gz',raw)
for filename in ['reproduction.tex','reproduction.pdf']:
    shutil.copy2(str(ROOT/'report'/filename),str(OUT/filename))
hashes=[]
for path in sorted(OUT.iterdir()):
    if path.name=='SHA256SUMS' or not path.is_file():continue
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    hashes.append(digest.hexdigest()+'  '+path.name)
(OUT/'SHA256SUMS').write_text('\n'.join(hashes)+'\n')
print('\n'.join(hashes))
