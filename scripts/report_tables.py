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

def front_speed(a):
    values=[s['speed_m_per_year'] for s in a['migration_segments'] if s['r2']>=.8]
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
        'lasts {} years and contains {} smaller seismic events. Operator tests verify steady Darcy flow, '
        'friction inversion, elasticity, fluid conservation, and fourth-order temporal convergence. '
        'The remaining trajectory comparisons retain explicit uncertainty from friction profiles, initial '
        'conditions, and numerical resolution; we do not claim recovery of the authors\' exact earthquake sequence.'
        ).format(counts['independently reproduced'],counts['partial'],counts['not reproduced'],
                 baseline['large_event_count'],number(completed['baseline']['years']),number(duration),baseline['selected_small_events'])
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
{\small\begin{verbatim}
'''+Path('data/verification.txt').read_text()+Path('data/coupled_verification.txt').read_text()+r'''\end{verbatim}
}
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
Case & $T$ (yr) & Large events & Median recurrence (yr) & Small events in selected cycle\\
\midrule
'''
    for case in MAIN_CASES:
        a=analyses[case]
        comparison+=' & '.join([tex(case),number(a['configuration']['T']/YEAR),str(a['large_event_count']),
            number(a['median_recurrence_years']),str(a['selected_small_events'])])+r'\\'+'\n'
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
The reference cases hold the initial pore pressure fixed. Each median uses all
complete large-event intervals in that run; early intervals may retain startup
effects. The selected figure cycle is specified separately in each panel sidecar.
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
from the authors' unspecified extraction procedure. The additional published
4.56 km/yr shallow-front annotation is not equated to this deep-front statistic.
'''
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
quantities in the selected independently generated cycle. It does not imply
that the timing, rupture depths, or four phases of the published cycle match.
'''
    comparison+=r'''\subsection{Dependence on healing time}
The article and Supplementary Figs. 3--6 predict reduced pressure cycling when
healing is slow relative to earthquake recurrence. The following quantities
are measured at 10 km in each run's selected cycle. They test that trend without
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

    validation={};ref_events=[e for e in baseline['events'] if e['large']]
    comparison+=r'''\subsection{Resolution and initial-condition sensitivity}
\begin{center}
\begin{tabular}{lrrrr}
\toprule
Case & Cells & Tolerance & First large event (yr) & Difference from baseline (yr)\\
\midrule
'''
    for case in ['baseline']+VALIDATION_CASES:
        a=analyses[case];ee=[e for e in a['events'] if e['large']]
        first=ee[0]['peak_s']/YEAR if ee else None
        reference=ref_events[0]['peak_s']/YEAR if ref_events else None
        difference=None if first is None or reference is None else first-reference
        validation[case]=dict(first_large_event_year=first,first_event_difference_year=difference,
            final_time_years=a['final_time_years'],large_event_count=a['large_event_count'],
            front_speed_m_per_year=front_speed(a),
            first_seismic_event=a['events'][0] if a['events'] else None)
        comparison+=' & '.join([tex(case),str(int(a['configuration']['n'])),number(a['configuration']['tolerance']),number(first),number(difference)])+r'\\'+'\n'
    comparison+=r'''\bottomrule
\end{tabular}
\end{center}
The coarse and fine meshes have spacings of 30.52 and 7.63 m. The tight case
reduces the integration tolerance by four. The initialization case doubles the
small state perturbation; the normal-stress case uses 23 rather than 22 MPa/km.
The fine-grid and input-sensitivity runs cover 100 years; other baseline
comparisons cover 200 years. Consequently their final-cycle statistics need not
sample the same cycle. First-event times compare a common initialization and
are reported without relabeling one event as another to improve agreement.
These tests quantify sensitivity; a small first-event difference alone does
not establish convergence of later swarm sequences.
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
      'F2':r'Independent counterpart to original Fig. 2. Reference fixed-pressure calculation (a,b) and coupled $T=10^8$ s calculation (c,d). Panels a and c retain both physical-time and accepted-step views of slip speed. Panels b and d show accumulated slip relative to cycle onset, with interseismic profiles in blue and seismic profiles in red. The step coordinate is specific to our adaptive integrator and is not a physical comparison metric.',
      'F3':r'Independent counterpart to original Fig. 3 for $T=10^8$ s. Effective stress (a,b), permeability (c,d), and upward fluid flux (e,f) are shown as time--depth maps and profiles. Early profiles are gray, interseismic profiles blue, and earthquake profiles red. Dashed black profiles denote the calculated initial steady state. Time selections and profile intervals are explicit in the panel files; published phase boundaries are not imposed.',
      'F4':r'Independent attempt at original Fig. 4. The fixed two-year window before the selected closing large event shows slip speed in time (a) and accepted-step coordinates (b), effective stress (c), and permeability (d) over 2--10 km. These are generated model fields even when swarm behavior is absent; the provenance table reports whether the target behavior was obtained.',
      'F5':r'Independent attempt at original Fig. 5, with $T=10^7$ s and $q_0=3.3\times10^{-10}$ m/s. Slip speed is shown in time (a) and accepted-step coordinates (b), with effective stress (c) and permeability (d). Any pulses arise from the coupled equations and the declared initial conditions; no pulse timing is prescribed.',
      'F6':r'Independent counterpart to original Fig. 6. Colored contours mark $|V|=V_p$ in the four calculated healing-time cases. Each curve uses the selected cycle of its own run, with time measured from that cycle onset. This explicit contour definition is an analysis assumption. Quantitative migration fits and the printed published target rates are compared in the text.',
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

    execution=r'''The calculation was executed in two stages while the report and plotting
code were developed. Both stages are the same stages invoked by the single
command above:
\begin{verbatim}
python3 scripts/reproduce.py --simulate-only --jobs 12
python3 scripts/reproduce.py --postprocess-only
\end{verbatim}
\begin{center}\begin{tabular}{lrrrr}
\toprule
Case & Years & Cells & Accepted steps & Solver wall time (min)\\\midrule
'''
    for case,c in completed.items():
        execution+=' & '.join([tex(case),number(c['years']),str(c['n']),str(c['accepted_steps']),number(c['elapsed_s']/60)])+r'\\'+'\n'
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
