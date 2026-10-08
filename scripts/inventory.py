"""Published panel inventory; labels follow the visible panels, not inferred data."""
# Approximate visual readings of displayed cycles, used only for comparison.
# They are never read by the numerical model or used to select its cycles.
PUBLISHED_CYCLES={
    'reference':(50.,'Main Fig. 2a'),
    'baseline':(32.,'Main Fig. 2c'),
    'long_reference':(160.,'Supplementary Fig. 3a, PDF page 4'),
    'long':(65.,'Supplementary Fig. 3c, PDF page 4'),
    'verylong_reference':(160.,'Supplementary Fig. 5a, PDF page 6'),
    'verylong':(150.,'Supplementary Fig. 5c, PDF page 6')}
PANELS={}
def add(name,quantity,kind,cases,documentation,field=None):
    PANELS[name]=dict(panel=name,quantity=quantity,kind=kind,cases=cases,
        documentation=documentation,field=field,
        data_command='python3 scripts/panel.py '+name+' --data',
        plot_command='python3 scripts/plot.py '+name,
        data_file='data/panels/'+name+'.npz')

add('F1a','Geometry schematic and assumed friction profiles','schematic',[],
    'Main Fig. 1a; Eqs. (7)-(9); Table 1; a(z), b(z) are assumed')
add('F1b','Permeability versus effective normal stress','stress_law',[],'Main Eq. (3), Fig. 1b; Table 1')
add('F1c','Healing and slip-enhancement limiting solutions','evolution_law',[],'Main Eq. (4), Fig. 1c; Table 1')
add('F1d','Steady pressure, stress, and permeability profiles','steady',[],
    'Main Eqs. (2)-(4), Fig. 1d; assumed normal-stress gradient')
for fig,coupled,reference,doc in [('F2','baseline','reference','Main Fig. 2'),
        ('S3','long','long_reference','Supplementary Fig. 3'),
        ('S5','verylong','verylong_reference','Supplementary Fig. 5')]:
    for letter,kind,case in [('a','velocity_pair',reference),('b','slip_profiles',reference),
                             ('c','velocity_pair',coupled),('d','slip_profiles',coupled)]:
        add(fig+letter,'Slip velocity in time and step coordinates' if kind=='velocity_pair' else 'Cumulative slip profiles',
            kind,[case],doc+'; main Eqs. (1)-(9), Table 1',1 if kind=='velocity_pair' else 0)
for fig,case,doc in [('F3','baseline','Main Fig. 3'),('S4','long','Supplementary Fig. 4'),('S6','verylong','Supplementary Fig. 6')]:
    for i,(quantity,field) in enumerate([('Effective normal stress',2),('Permeability',3),('Upward fluid flux',5)]):
        for j,kind in enumerate(['heatmap','profiles']):
            add(fig+chr(97+2*i+j),quantity+(' history' if j==0 else ' depth profiles'),kind,[case],
                doc+'; main Eqs. (1)-(4), Table 1',field)
for fig,case in [('F4','baseline'),('F5','short')]:
    for letter,kind,field,quantity in [('a','heatmap',1,'Slip velocity versus time'),
            ('b','step_heatmap',1,'Slip velocity versus accepted step'),
            ('c','heatmap',2,'Effective normal stress'),('d','heatmap',3,'Permeability')]:
        add(fig+letter,quantity,kind,[case],'Main '+fig.replace('F','Fig. ')+
            '; Eqs. (1)-(9); Table 1'+('; T=1e7 s, q0=3.3e-10 m/s' if fig=='F5' else ''),field)
add('F6','Upward migration of the |V|=Vp contour','fronts',['short','baseline','long','verylong'],
    'Main Fig. 6 and Discussion; reported rates 2.50, 0.38, 0.12, 0.03 km/year')
for letter,field,quantity in [('a',1,'Slip velocity'),('b',3,'Permeability'),('c',5,'Upward flux'),('d',2,'Effective stress')]:
    add('S1'+letter,quantity+' at 5, 10, 15, and 20 km','depth_histories',['baseline'],
        'Supplementary Fig. 1; main Eqs. (1)-(9)',field)
for letter in ['a','b']:
    add('S2'+letter,'Alternative cycle: slip velocity in time and step coordinates','velocity_pair',['baseline'],
        'Visible panels a,b of Supplementary Fig. 2; caption additionally mentions absent c,d',1)

GROUPS={}
for panel in PANELS:
    figure=panel[:2]
    GROUPS.setdefault(figure,[]).append(panel)
