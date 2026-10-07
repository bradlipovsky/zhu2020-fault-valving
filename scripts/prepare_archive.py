"""Extract plotting arrays without changing archived values (50 m near surface)."""
import json
from pathlib import Path
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / 'build' / 'archive'
CACHE.mkdir(parents=True, exist_ok=True)
KEYS = ['slip_velocity', 'slip', 'effective_normal_stress', 'permeability', 'flux', 'pore_pressure']

if __name__ == '__main__':
    summary = {}
    for p in sorted((ROOT / 'reference' / 'data').glob('*.mat')):
        dest = CACHE / (p.stem + '_si.npz')
        with h5py.File(p) as f:
            z = f['depth'][:]
            t = f['time'][:].ravel()
            nz = np.searchsorted(z, 25, side='right')
            summary[p.stem] = dict(samples=len(t), depth_nodes=len(z), bottom_km=float(z[-1]),
                start_seconds=float(t[0]), end_seconds=float(t[-1]), duration_years=float((t[-1]-t[0])/31536000))
            if dest.exists():
                continue
            print('Extracting', p.stem, flush=True)
            data = {key: f[key][:nz:5, :].astype(np.float32) for key in KEYS if key in f}
            # Scycle uses km^2 internally and the OSF archive retains these units.
            if 'permeability' in data:
                data['permeability'] *= 1e6
            data.update(z=z[:nz:5], t=(t-t[0])/31536000, absolute_time_s=t)
            np.savez_compressed(dest, **data)
    (ROOT / 'results' / 'archive_inventory.json').write_text(json.dumps(summary, indent=2)+'\n')
