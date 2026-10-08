#!/usr/bin/env python3
"""Write report statements from measured independent outputs and explicit targets."""
from pathlib import Path
from collections import Counter
import json
import math
import numpy as np
from common import YEAR, MAIN_CASES, save_json
from inventory import PANELS
from reproduce import VALIDATION_CASES

def tex(s):
    for a,b in [('&',r'\&'),('%',r'\%'),('_',r'\_'),('#',r'\#')]:s=str(s).replace(a,b)
    return s

def number(x,digits=3):
    return 'not measured' if x is None or not math.isfinite(x) else ('{:.'+str(digits)+'g}').format(x)

def command(s):return r'\texttt{'+tex(s).replace('--','-{}-')+'}'

def write(name,content):Path('report',name+'.tex').write_text(content+'\n')

def front_speed(a,field='migration_segments'):
    values=[s['speed_m_per_year'] for s in a[field] if s['r2']>=.8]
    return float(np.median(values)) if values else None

def main():
    panels={p:json.loads(Path('data/panels',p+'.json').read_text()) for p in PANELS}
    analyses={c:json.loads(Path('data',c,'analysis.json').read_text()) for c in MAIN_CASES+VALIDATION_CASES}
    completed={c:json.loads(Path('data',c,'completed.json').read_text()) for c in analyses}
    baseline=analyses['baseline'];counts=Counter(p['status'] for p in panels.values())
    duration=(baseline['selected_cycle_s'][1]-baseline['selected_cycle_s'][0])/YEAR
    summary=r'\begin{abstract}'+'\n'
    summary+=('We independently implement the coupling between earthquake slip, fluid flow, and permeability '
        'described by Zhu et al. (2020). The new calculation uses the published equations and declared '
        'assumptions where initialization and parameter profiles are incomplete. It generates 49 panel '
        'products: one includes a newly drawn schematic and 48 contain only numerical calculations. '
        'The panel audit classifies {} as independently reproduced, {} as partial, and {} as not reproduced. '
        'The two independently reproduced panels are the constitutive-law limits in Fig. 1b,c. '
        'The baseline calculation produces {} large events over {} years. Its selected comparison window '
        'lasts {} years and contains {} smaller resolved ruptures. Operator tests verify steady Darcy flow, '
        'friction inversion, elasticity, fluid conservation, and fourth-order temporal convergence. '
        'The remaining trajectory comparisons retain explicit uncertainty from friction profiles, initial '
        'conditions, and numerical resolution; we do not claim recovery of the authors\' exact earthquake sequence.'
        ).format(counts['independently reproduced'],counts['partial'],counts['not reproduced'],
                 baseline['large_event_count'],number(completed['baseline']['years']),number(duration),baseline['selected_small_ruptures'])
    summary+='\n'+r'\end{abstract}'
    write('summary',summary)

    verification=r'''\subsection{Verification of equations and limiting cases}
The tests invoke the actual operators in the production source. They check a
steady nonlinear Darcy solution, signed friction inversion over slip speeds
from $10^{-15}$ to $10$ m s$^{-1}$, constant and nonconstant elastic modes,
hydrostatic equilibrium, and constant-speed state/permeability limits. The
nonlinear storage test perturbs both pressure and permeability; its residual
is the storage change minus integrated boundary flux. A separate coupled
calculation refines all time-dependent fields together.

For an independent steady-flow check, define $G=\partial_z\sigma-\rho g$,
$C=k_*-k_{\min}$, $k_{\rm eq}=\eta q_0/G$, and $D=k_{\rm eq}-k_{\min}$.
With constant $k_*$, the documented equations reduce to
$dN/dz=G-\eta q_0/[k_{\min}+C\exp(-N/\sigma_*)]$, where $N=\sigma-p$.
Separating this equation and imposing $N(0)=0$ gives
\begin{equation}
 z=\frac{N}{G}-\frac{\sigma_*k_{\rm eq}}{GD}
 \ln\!\left[\frac{C-D\exp(N/\sigma_*)}{C-D}\right].
\end{equation}
In the additional limit $k_{\min}=0$, the explicit solution is
\begin{equation}
 N(z)=-\sigma_*\ln\!\left[\frac{k_{\rm eq}}{C}
       +\left(1-\frac{k_{\rm eq}}{C}\right)\exp(-Gz/\sigma_*)\right].
\end{equation}
We compare the production initial-pressure solver with these separately
integrated solutions for both the published permeability floor and its zero limit.
Errors decrease by approximately four when grid spacing is halved.
{\small\begin{verbatim}
'''+Path('data/verification.txt').read_text()+Path('data/coupled_verification.txt').read_text()+Path('data/steady_verification.txt').read_text()+Path('data/graded_verification.txt').read_text()+Path('data/analysis_verification.txt').read_text()+r'''\end{verbatim}
}
Catalog and phase-diagnostic checks use synthetic output arrays in a temporary
directory. They verify event boundaries, connected depth intervals, and pressure
extrema; those test arrays never supply the scientific figures.
The diffusion errors decrease by approximately four on each joint space/time
refinement. The coupled test compares fixed steps of 20,000, 10,000, and 5,000 s
with a 2,500 s reference for the second-order comparison integrator. The production
fourth-order method uses approximately 80,000, 40,000, 20,000, and 10,000 s steps,
again against a 2,500 s reference. These are temporal-order tests in a separate
30 km domain, not evidence that the 500 km earthquake calculation is spatially converged.
'''
    write('verification',verification)

    comparison=r'''\begin{center}
\begin{tabular}{lrrrr}
\toprule
Case & $T$ (yr) & Large events & Median recurrence (yr) & Small ruptures in window\\
\midrule
'''
    for case in MAIN_CASES:
        a=analyses[case]
        comparison+=' & '.join([tex(case),number(a['configuration']['T']/YEAR),str(a['large_event_count']),
            number(a['median_recurrence_years']),str(a['selected_small_ruptures'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
The reference cases hold the initial pore pressure fixed. Each median uses all
complete large-event intervals in that run; early intervals may retain startup
effects. The selected figure cycle is specified separately in each panel sidecar.
Small-rupture counts require a resolved 1 cm slip footprint. Every speed-threshold
interval remains in the catalog, including those without such a footprint.
The main reference figure shows an approximately 50-year interval, while the
featured fault-valving cycle ends at approximately 32 years. These are approximate
readings of the displayed time axes, not digitized curves or input data.
'''
    for case,target in [('reference',50.),('baseline',32.)]:
        measured=analyses[case]['median_recurrence_years']
        if measured is not None:
            comparison+='Our {} median is {} yr, a {}\\% difference from that approximate {} yr target.\n'.format(
                tex(case),number(measured),number(100*(measured/target-1)),number(target))
        else:comparison+='No complete large-event recurrence interval was obtained for {}.\n'.format(tex(case))
    comparison+=r'''\begin{center}
\begin{tabular}{lrrr}
\toprule
Case & Published deep-front rate (m/yr) & Calculated median (m/yr) & Difference (\%)\\
\midrule
'''
    target_rates={'short':2500.,'baseline':380.,'long':120.,'verylong':30.}
    metrics={}
    for case,target in target_rates.items():
        rate=front_speed(analyses[case]);diff=None if rate is None else 100*(rate/target-1)
        comparison+=tex(case)+' & '+number(target)+' & '+number(rate)+' & '+number(diff)+r'\\'+'\n'
        metrics[case]=dict(published_rate_m_per_year=target,calculated_rate_m_per_year=rate,relative_difference_percent=diff)
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
Published rates are the numerical annotations in Fig. 6 and the Discussion.
Calculated rates use upward portions of the $|V|=V_p$ contour between 13 and
20 km, spanning at least 1 km and 0.1 yr, with linear-fit $R^2\geq0.8$.
All fitted segments and their fit quality are saved, including rejected low-quality
fits. The rate statistic therefore has an explicit definition and may differ
from the authors' unspecified extraction procedure. A separate calculation uses
the same fit criteria at 2--10 km depth to compare with the published shallow-front
annotation of 4.56 km/yr. The depth interval is an explicit analysis choice.
'''
    shallow=front_speed(baseline,'shallow_migration_segments')
    shallow_difference=None if shallow is None else 100*(shallow/4560-1)
    metrics['baseline']['shallow_migration']=dict(depth_range_m=[2000,10000],published_rate_m_per_year=4560,
        calculated_rate_m_per_year=shallow,relative_difference_percent=shallow_difference)
    if shallow is None:
        comparison+='No shallow segment satisfies the stated duration, distance, and fit-quality criteria.\n'
    else:
        comparison+=('The calculated shallow median is {} m/yr, a {}\\% difference '
            'from the printed 4560 m/yr annotation.\n').format(number(shallow),number(shallow_difference))
    comparison+=r'''\begin{center}
\begin{tabular}{lrrrr}
\toprule
Baseline depth (km) & Minimum $N$ (MPa) & Maximum $N$ (MPa) & Range (MPa) & Maximum $q$ (m/s)\\
\midrule
'''
    for key in ['5000','10000','15000','20000']:
        a=baseline['sample_depths'][key]
        comparison+=' & '.join([number(int(key)/1000),number(a['effective_min_mpa']),number(a['effective_max_mpa']),
            number(a['effective_range_mpa']),number(a['flux_max'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
The article describes effective-stress changes of order 10--20 MPa and an
upper flux scale of order $10^{-7}$ m s$^{-1}$. The table measures those
quantities in the selected independently generated window. It does not imply
that the timing, rupture depths, or four phases of the published cycle match.
'''
    phase=baseline['phase_diagnostics'];metrics['baseline']['phase_diagnostics']=phase
    comparison+=r'''\subsection{Rupture depths and drainage timing}
The published description places two partial ruptures at approximately
13--18 and 8--19 km depth before the final swarm and surface-reaching event.
For a quantitative spatial comparison, we list the first two complete partial
events in chronological order within our comparison window.
Their depth ranges are the longest connected components exceeding 1 cm of event
slip. A later swarm event can have larger local slip than an earlier partial
rupture, so maximum slip does not determine the order. The full event catalog
is retained. This rule does not select the events closest to the published
footprints or establish an automatic correspondence with the published phases.
\begin{center}\begin{tabular}{lrrr}
\toprule
Time in window (yr) & Connected footprint (km) & Maximum slip (m) & Peak speed (m/s)\\
\midrule
'''
    for event in phase['first_partial_ruptures']:
        comparison+=' & '.join([number(event['peak_since_window_start_years']),
            number(event['footprint_top_m']/1000)+'--'+number(event['footprint_bottom_m']/1000),
            number(event['max_slip_m']),number(event['peak_velocity'])])+r'\\'+'\n'
    if not phase['first_partial_ruptures']:
        comparison+=r'\multicolumn{4}{c}{No complete partial rupture with a resolved 1 cm footprint.}\\'+'\n'
    comparison+=r'\bottomrule\end{tabular}\end{center}'+'\n'
    comparison+=('At the sampled depth of {} km, effective stress reaches its window maximum '
        '{} yr after the window starts; pressure has decreased by {} MPa from its initial value.\n').format(
            number(phase['sample_depth_m']/1000),number(phase['maximum_effective_stress_since_window_start_years']),
            number(phase['pressure_drop_from_window_start_mpa']))
    comparison+=r'''The article describes a transition from postseismic drainage to renewed
pressurization after approximately 5--10 years. Our stress-maximum time is an
explicit diagnostic at one depth, rather than the authors' unspecified phase
boundary. It must be read with the pressure maps. When no complete large-event
cycle exists, the comparison window covers the entire run, including startup;
such a window cannot establish the four-phase cycle. Events truncated by an
output boundary cannot close a complete cycle and are excluded from this table.
'''
    comparison+=r'''\subsection{Dependence on healing time}
The article and Supplementary Figs. 3--6 predict reduced pressure cycling when
healing is slow relative to earthquake recurrence. The following quantities
are measured at 10 km in each run's selected window. They test that trend without
equating cycles that begin from different independently generated states.
\begin{center}\begin{tabular}{lrrrr}
\toprule
Case & $N_{\min}$ (MPa) & $N_{\max}$ (MPa) & $k_{\max}/k_{\min}$ & $q_{\max}/q_0$\\
\midrule
'''
    for case in ['short','baseline','long','verylong']:
        sample=analyses[case]['sample_depths']['10000']
        comparison+=' & '.join([tex(case),number(sample['effective_min_mpa']),number(sample['effective_max_mpa']),
            number(sample['permeability_max']/sample['permeability_min']),
            number(sample['flux_max']/analyses[case]['configuration']['influx'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
Here subscripts min and max denote temporal extrema at the sampled depth;
the permeability ratio does not refer to the constitutive bounds. These are
amplitude comparisons. They do not establish agreement in individual event
times or demonstrate grid convergence of the long-healing cases.
Entries from runs with fewer than two complete large events include startup;
they cannot establish the amplitude of a repeating cycle.
'''
    comparison+=r'''\subsection{Aseismic pulse diagnostics}
For the short-healing case, a local episode has $|V|\geq1.1V_p$ for at least
0.01 yr, with the maximum speed over the whole fault remaining below
$10^{-3}$ m/s. Both ends must lie inside the selected comparison window.
This independent measurement convention excludes earthquakes and truncated
episodes; it is not an imposed condition on the equations.
\begin{center}\begin{tabular}{lrrrrr}
\toprule
Depth (km) & Episodes & Duration (yr) & Interval (yr) & Slip (cm) & Interval CV\\
\midrule
'''
    for depth,diagnostic in analyses['short']['slow_slip'].items():
        comparison+=' & '.join([number(float(depth)/1000),str(len(diagnostic['episodes'])),
            number(diagnostic['median_duration_years']),number(diagnostic['median_interval_years']),
            number(100*diagnostic['median_net_slip_m'] if diagnostic['median_net_slip_m'] is not None else None),
            number(diagnostic['interval_coefficient_of_variation'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
Duration, interval, and net slip are medians over complete episodes.
CV is the standard deviation divided by the mean interpeak interval and requires
at least three complete episodes. These measurements assess the paper's
quasi-periodic-pulse behavior without prescribing its pulse times. The Discussion
reports pulses lasting about one year, recurring every few years, and producing
a few centimetres of slip. Those are comparison scales, not precise period or
CV targets. Our threshold defines the measured pulse width; the authors do not
specify an equivalent threshold, so duration comparisons carry that ambiguity.
'''

    validation={};ref_events=[e for e in baseline['events'] if e['large'] and e['complete']]
    comparison+=r'''\subsection{Resolution and initial-condition sensitivity}
\begin{center}
\begin{tabular}{lrrrr}
\toprule
Case & Finest spacing (m) & Tolerance & First large event (yr) & Difference (yr)\\
\midrule
'''
    for case in ['baseline']+VALIDATION_CASES:
        a=analyses[case];events=[e for e in a['events'] if e['complete']];ee=[e for e in events if e['large']]
        first=ee[0]['peak_s']/YEAR if ee else None
        reference=ref_events[0]['peak_s']/YEAR if ref_events else None
        difference=None if first is None or reference is None else first-reference
        validation[case]=dict(first_large_event_year=first,first_event_difference_year=difference,
            final_time_years=a['final_time_years'],large_event_count=a['large_event_count'],
            front_speed_m_per_year=front_speed(a),
            first_seismic_event=events[0] if events else None)
        comparison+=' & '.join([tex(case),number(a['resolution']['dz_m']),number(a['configuration']['tolerance']),number(first),number(difference)])+r'\\'+'\n'
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
The coarse and fine meshes have minimum spacings of 7.63 and 1.91 m, compared
with 3.81 m in the baseline. The maximum deep spacing also halves on refinement.
The tight case
reduces the integration tolerance by four. The initialization case doubles the
small state perturbation; the normal-stress case uses 23 rather than 22 MPa/km.
The fine-grid and input-sensitivity runs cover 100 years; other baseline
comparisons cover 200 years. Consequently their final-cycle statistics need not
sample the same cycle. First-event times compare a common initialization and
are reported without relabeling one event as another to improve agreement.
These tests quantify sensitivity; a small first-event difference alone does
not establish convergence of later swarm sequences.
\begin{center}\begin{tabular}{lrrr}
\toprule
Case & Physical fault nodes & Minimum $L_b/\Delta z$ & Minimum $h^*/\Delta z$\\
\midrule
'''
    for case in MAIN_CASES:
        resolution=analyses[case]['resolution']
        comparison+=' & '.join([tex(case),str(resolution['dynamic_nodes']),number(resolution['minimum_Lb_cells']),
            number(resolution['minimum_hstar_cells'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
Here $L_b=\mu d_c/[N b]$ and $h^*=\mu d_c/[N(b-a)]$, with $N=\sigma-p$.
The minima include every accepted step and every velocity-weakening node;
$\Delta z$ is its actual control-volume width. These diagnostics follow the
length scales discussed by Allison and Dunham\ \cite{allison}, Eqs. (20)--(21).
They quantify local resolution and do not alone establish sequence convergence.
\begin{center}\begin{tabular}{lrrrr}
\toprule
Case & First seismic onset (yr) & Peak speed (m/s) & Maximum slip (m) & Large?\\
\midrule
'''
    for case,metrics_case in validation.items():
        event=metrics_case['first_seismic_event']
        comparison+=' & '.join([tex(case),number(event['start_s']/YEAR if event else None),
            number(event['peak_velocity'] if event else None),number(event['max_slip_m'] if event else None),
            ('yes' if event['large'] else 'no') if event else 'none'])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
The first seismic event can be a small rupture that precedes the first large
event. The table retains this distinction when comparing resolution. Peak speed
comes from every accepted step; slip and the rupture footprint use stored fields.
'''
    comparison+=r'''\begin{center}\begin{tabular}{lrrrr}
\toprule
Case & \shortstack{Onset\\difference (s)} & \shortstack{Duration\\difference (\%)} & \shortstack{Peak speed\\difference (\%)} & \shortstack{Maximum slip\\difference (\%)}\\
\midrule
'''
    reference=validation['baseline']['first_seismic_event']
    for case in ['baseline_coarse','baseline_fine','baseline_tight']:
        event=validation[case]['first_seismic_event'];difference=None
        if event and reference:
            difference=dict(onset_s=event['start_s']-reference['start_s'],
                duration_percent=100*((event['end_s']-event['start_s'])/(reference['end_s']-reference['start_s'])-1),
                peak_speed_percent=100*(event['peak_velocity']/reference['peak_velocity']-1),
                maximum_slip_percent=100*(event['max_slip_m']/reference['max_slip_m']-1)
                    if reference['max_slip_m']>0 else None)
        validation[case]['first_event_difference_from_baseline']=difference
        values=[difference[key] if difference else None for key in
            ['onset_s','duration_percent','peak_speed_percent','maximum_slip_percent']]
        comparison+=' & '.join([tex(case)]+[number(value) for value in values])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
Differences are signed relative to the baseline's first complete seismic event.
Duration is the whole-fault threshold excursion defined above. These comparisons
use a common initialization and include the finite saved-profile spacing in the
slip measurement; they do not establish convergence of subsequent event sequences.
'''
    onset_order=None
    onset_events=[validation[case]['first_seismic_event'] for case in
        ['baseline_coarse','baseline','baseline_fine']]
    if all(onset_events):
        times=[event['start_s'] for event in onset_events]
        coarse_change=abs(times[0]-times[1]);fine_change=abs(times[1]-times[2])
        spacings=[analyses[case]['resolution']['dz_m'] for case in
            ['baseline_coarse','baseline','baseline_fine']]
        ratio=spacings[0]/spacings[1]
        if coarse_change>0 and fine_change>0 and ratio>1 and abs(ratio-spacings[1]/spacings[2])<1e-10:
            onset_order=dict(coarse_to_main_change_s=coarse_change,main_to_fine_change_s=fine_change,
                spacing_ratio=ratio,apparent_order=math.log(coarse_change/fine_change)/math.log(ratio))
    validation['baseline']['first_onset_refinement']=onset_order
    if onset_order:
        comparison+=('Successive mesh refinements change the first seismic onset by {} and {} s, '
            'giving an apparent order of {} from the ratio of these differences. '
            'This estimate includes residual time-integration error and the finite threshold-crossing bracket; '
            'it applies only to this first-onset statistic.\n').format(number(coarse_change),number(fine_change),
                number(onset_order['apparent_order']))
    horizon=min(analyses[case]['final_time_years'] for case in validation)
    large_times={case:[e['end_s']/YEAR for e in analyses[case]['events']
        if e['complete'] and e['large'] and e['end_s']/YEAR<=horizon] for case in validation}
    comparison+=('We also compare every validation case over the first {} years, '
        'the common duration available in all runs. This window includes startup.\n').format(number(horizon))
    comparison+=r'''\begin{center}\begin{tabular}{lrrrr}
\toprule
Case & \shortstack{Fast intervals\\(resolved ruptures)} & \shortstack{Large\\events} & \shortstack{Median recurrence\\(yr)} & \shortstack{Maximum timing\\difference (yr)}\\
\midrule
'''
    for case in validation:
        events=[e for e in analyses[case]['events'] if e['complete'] and e['end_s']/YEAR<=horizon]
        event_count=len(events);resolved_count=sum(bool(e['rupture_intervals_m']) for e in events)
        times=large_times[case];intervals=np.diff(times);matched=min(len(times),len(large_times['baseline']))
        differences=(np.asarray(times[:matched])-np.asarray(large_times['baseline'][:matched])).tolist()
        common=dict(horizon_years=horizon,complete_threshold_intervals=event_count,
            resolved_ruptures=resolved_count,unresolved_threshold_intervals=event_count-resolved_count,large_events=len(times),
            large_event_end_years=times,recurrence_years=intervals.tolist(),
            median_recurrence_years=float(np.median(intervals)) if len(intervals) else None,
            matched_ordinal_large_events=matched,large_event_time_differences_years=differences,
            maximum_absolute_large_event_time_difference_years=max(map(abs,differences)) if differences else None)
        validation[case]['common_horizon_comparison']=common
        comparison+=' & '.join([tex(case),'{} ({})'.format(event_count,resolved_count),str(len(times)),number(common['median_recurrence_years']),
            number(common['maximum_absolute_large_event_time_difference_years'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule\end{tabular}\end{center}
The first count includes all complete excursions above the seismic-speed
threshold; the parenthesized count also requires a resolved 1 cm slip footprint.
Brief threshold recrossings can increase the former without producing an
additional resolved rupture. They are retained separately, without merging
intervals or changing the detection threshold. The timing difference pairs the first classified large event with the first,
the second with the second, and so on, without fitting a time shift or choosing
the closest event. Only pairs present in both catalogs enter that
difference; the event counts expose missing or additional events. The numerical
summary retains every paired difference and recurrence interval. Event counts
and recurrence over this common window test sequence sensitivity, although they
do not alone establish convergence of individual swarm ruptures.

Exact friction profiles, prestress, startup procedure, spatial grid, and the
authors' saved cycle states cannot be recovered from the article and supplement.
The independent calculation therefore cannot establish identity with those
trajectories. Its additional numerical differences are the cosine elasticity
operator, finite-volume flow operator, IMEX integrator, and cell-center surface
treatment. None is concealed by importing an authors' restart or output field.
'''
    save_json('data/quantitative_comparisons.json',metrics);save_json('data/validation_summary.json',validation)
    write('comparisons',comparison.replace(r'\begin{center}',r'\begin{center}\small'))

    provenance=r'''\begin{landscape}
\scriptsize
\setlength{\tabcolsep}{3pt}
\begin{longtable}{P{1.15cm}P{3.5cm}P{3.4cm}P{4.1cm}P{3.4cm}P{2.1cm}P{5.1cm}}
\caption{Figure-by-figure provenance. Input parameters and selections accompany each numerical file in a JSON sidecar.}\label{tab:provenance}\\
\toprule
Original panel & Quantity and output data & Documentation used & Command generating underlying data & Command plotting it & Status & Main discrepancy\\
\midrule\endfirsthead
\toprule
Original panel & Quantity and output data & Documentation used & Command generating underlying data & Command plotting it & Status & Main discrepancy\\
\midrule\endhead
\bottomrule\endfoot
'''
    order=['F1','F2','F3','F4','F5','F6','S1','S2','S3','S4','S5','S6']
    for prefix in order:
        for name,m in panels.items():
            if not name.startswith(prefix):continue
            fields=[name,tex(m['quantity'])+r'\par\panelpath{'+name+'}',tex(m['documentation']),
                command(m['data_command']),command(m['plot_command']),
                tex(m['status']),tex(m['main_discrepancy'])]
            provenance+=' & '.join(fields)+r'\\[5pt]'+'\n'
    provenance+=r'\end{longtable}\end{landscape}'
    write('provenance',provenance)

    captions={
      'F1':r'Independent counterpart to original Fig. 1. (a) Newly drawn geometry schematic and explicitly assumed friction profiles; the diagram is not a simulated result. (b) Stress-dependent permeability evaluated from the published exponential law. (c) Exact healing with zero slip and slip enhancement with negligible healing. (d) Discrete steady Darcy flow and stress profiles for $T=10^8$ s and $q_0=3\times10^{-9}$ m/s. The assumed normal-stress gradient is 22 MPa/km; the illustrative lithostatic line assumes 2700 kg/m$^3$.',
      'F2':r'Independent counterpart to original Fig. 2. Reference fixed-pressure calculation (a,b) and coupled $T=10^8$ s calculation (c,d). Panels a and c retain both physical-time and accepted-step views of slip speed. Panels b and d show accumulated slip relative to the selected window start, with interseismic profiles in blue and seismic profiles in red. The step coordinate is specific to our adaptive integrator and is not a physical comparison metric.',
      'F3':r'Independent counterpart to original Fig. 3 for $T=10^8$ s. Effective stress (a,b), permeability (c,d), and upward fluid flux (e,f) are shown as time--depth maps and profiles. Early profiles are gray, interseismic profiles blue, and earthquake profiles red. Dashed black profiles denote the calculated initial steady state. Time selections and profile intervals are explicit in the panel files; published phase boundaries are not imposed.',
      'F4':r'Independent attempt at original Fig. 4. The fixed two-year window before the selected closing large event shows slip speed in time (a) and accepted-step coordinates (b), effective stress (c), and permeability (d) over 2--10 km. These are generated model fields even when swarm behavior is absent; the provenance table reports whether the target behavior was obtained.',
      'F5':r'Independent attempt at original Fig. 5, with $T=10^7$ s and $q_0=3.3\times10^{-10}$ m/s. Slip speed is shown in time (a) and accepted-step coordinates (b), with effective stress (c) and permeability (d). Any pulses arise from the coupled equations and the declared initial conditions; no pulse timing is prescribed.',
      'F6':r'Independent counterpart to original Fig. 6. Colored contours mark $|V|=V_p$ in the four calculated healing-time cases. Each curve uses its own selected window: the final complete large-event cycle when available, or the entire run otherwise. Time is measured from that window start. This explicit contour definition is an analysis assumption. Quantitative migration fits and the printed published target rates are compared in the text.',
      'S1':r'Independent counterpart to Supplementary Fig. 1. Histories at the cells containing 5, 10, 15, and 20 km show slip speed (a), permeability (b), upward flux (c), and effective stress (d) for the baseline case. Every accepted time sample is retained in the numerical panel data. The window covers the last two complete large-event cycles when available.',
      'S2':r'Independent attempt at the two visible rows of Supplementary Fig. 2. Rows a and b show earlier complete baseline cycles selected in a fixed chronological order, each in physical-time and accepted-step coordinates. They are not selected to resemble the published alternatives. If insufficient cycles exist, the full run is shown and the panel is marked not reproduced.',
      'S3':r'Independent counterpart to Supplementary Fig. 3 for $T=10^9$ s. Fixed-pressure reference (a,b) and coupled flow (c,d) are shown as velocity maps and slip profiles. Interseismic profiles use four-year target intervals; seismic profiles use one-second target intervals sampled from stored output.',
      'S4':r'Independent counterpart to Supplementary Fig. 4 for $T=10^9$ s. Effective stress (a,b), permeability (c,d), and flux (e,f) are shown in time--depth maps and profiles. Gray curves mark the first ten years, blue interseismic profiles use four-year target intervals, and red seismic profiles use one-second target intervals. The initial steady solution is dashed black.',
      'S5':r'Independent counterpart to Supplementary Fig. 5 for $T=10^{10}$ s. Fixed-pressure reference (a,b) and coupled flow (c,d) are shown as velocity maps and slip profiles, using the same independently specified friction and initial conditions as the other experiments.',
      'S6':r'Independent counterpart to Supplementary Fig. 6 for $T=10^{10}$ s. Effective stress (a,b), permeability (c,d), and upward flux (e,f) test the predicted reduction in fault-valving amplitude when healing is slow. All curves and images are generated from the new implementation.'}
    pages=''
    for group in order:
        pages+=r'\clearpage\begin{figure}[p]\centering'+'\n'+r'\includegraphics[width=\textwidth,height=.77\textheight,keepaspectratio]{../figures/'+group+'.pdf}\n'
        pages+=r'\caption{'+captions[group]+r' Comparison target: Zhu et al.~\cite{zhu}.}\label{fig:'+group+r'}\end{figure}'+'\n'
    write('figure_pages',pages)

    execution=r'''The scientific calculation was launched with the complete regeneration
command. It runs the equation checks, numerical cases, analysis, plotting,
report compilation, source-bundle compilation, and artifact audit in sequence:
\begin{verbatim}
python3 scripts/reproduce.py --jobs 12
\end{verbatim}
The initial 200-year fixed-pressure reference completed successfully but
contained only one large earthquake, with no complete recurrence interval.
That duration check is recorded in
{\small\path{documentation/reference_duration.json}} and excluded from the final
figure inputs. The reference configuration was extended to 400 years and
rerun from the same independent initial conditions. The concurrent replacement
run used a separate progress file:
\begin{verbatim}
python3 scripts/reproduce.py --simulate-only --cases reference \
  --jobs 1 --status-file data/reference_extension_status.json
\end{verbatim}
The standard regeneration command now uses the 400-year reference directly.
\begin{center}\small\begin{tabular}{lrrrrr}
\toprule
Case & Years & Fault nodes & FFT points & Accepted steps & Solver time (min)\\\midrule
'''
    for case,c in completed.items():
        execution+=' & '.join([tex(case),number(c['years']),str(c['dynamic_nodes']),str(c['n']),str(c['accepted_steps']),number(c['elapsed_s']/60)])+r'\\'+'\n'
    execution+=r'''\bottomrule\end{tabular}\end{center}
Runs execute concurrently, so the solver wall times must not be added to infer
total project elapsed time. Each run records source, executable, and configuration
hashes, command, process identifier, completion status, and output extent.
Verification and pipeline logs are retained. The tested compiler and library
versions are listed in \code{documentation/environment.json}. Panel data and
figure checksums are recorded in the artifact manifest. Reported computational
times describe this workspace and are not hardware-independent benchmarks.
'''
    write('execution',execution)

if __name__=='__main__':main()
