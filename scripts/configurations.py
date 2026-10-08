#!/usr/bin/env python3
"""Write all declared physical and numerical inputs, with no external data."""
from pathlib import Path

BASE = dict(n=32768, threads=4, Ly=500000, Lz=500000, mu=32.4e9,
    radiation=4.68e6, Vp=1e-9, V0=1e-6, f0=.6, dc=.002,
    rho=1000, gravity=9.8, viscosity=1e-4, storage=1e-11,
    normal_gradient=22000, stress_scale=30e6, kmin=1e-19, kmax=1e-15,
    L=1, T=1e8, influx=3e-9,
    a_surface=.01, a_15km=.03, a_deep_slope=3e-6,
    ab_shallow=-.01, ab_corner=13500, vw_bottom=17000,
    perturbation=1e-4, perturb_depth=12000, perturb_width=1000,
    initial_locking_depth=17000, initial_transition_width=1500, initial_velocity_fraction=1e-4,
    years=200, tolerance=.001, dt_max=1e6,
    output_years=.025, output_seconds=.5, output_depth=30000,
    output_spacing=30, seismic_threshold=.001, fixed_pressure=0)

CASES = {
    'reference': dict(fixed_pressure=1, years=200),
    'baseline': {},
    'short': dict(T=1e7, influx=3.3e-10),
    'long': dict(T=1e9, years=350),
    'long_reference': dict(T=1e9, fixed_pressure=1, years=350),
    'verylong': dict(T=1e10, years=500),
    'verylong_reference': dict(T=1e10, fixed_pressure=1, years=500),
    'smoke': dict(n=2048, threads=1, years=.05),
    'baseline_coarse': dict(n=16384),
    'baseline_fine': dict(n=65536, years=100),
    'baseline_tight': dict(tolerance=.00025),
    'initialization': dict(perturbation=2e-4, years=100),
    'normal_stress': dict(normal_gradient=23000, years=100),
}

if __name__ == '__main__':
    Path('configs').mkdir(exist_ok=True)
    for name, change in CASES.items():
        values = dict(BASE, **change)
        text = ('# Independent model; all values are SI except explicitly named years.\n'
                '# See documentation/independence.md for assumed inputs.\n')
        text += ''.join('{} = {:.16g}\n'.format(k,v) for k,v in values.items())
        Path('configs',name+'.cfg').write_text(text)
