#!/usr/bin/env python3
"""Compile a dated scientific report using only completed, audited case outputs."""
from pathlib import Path
from collections import Counter
import datetime
import hashlib
import json
import os
import subprocess
import tempfile
import time
import zipfile
from common import MAIN_CASES, YEAR
from inventory import PANELS, GROUPS, PUBLISHED_CYCLES
from reproduce import VALIDATION_CASES
from report_tables import tex, number, command, front_speed, baseline_supplementary

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    os.chdir(str(ROOT))
    scratch = ROOT / '.tmp'
    scratch.mkdir(exist_ok=True)
    os.environ['TMPDIR'] = str(scratch)
    stamp = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
    inputs = {}

    def read(path):
        inputs[str(path)] = digest(path)
        return json.loads(Path(path).read_text())

    completed, analyses, audits = {}, {}, {}
    source_hash = digest('src/model.cpp')
    for case in MAIN_CASES + VALIDATION_CASES:
        directory = Path('data', case)
        if not (directory / 'completed.json').exists():
            continue
        done = read(directory / 'completed.json')
        audit = read(directory / 'output_verification.json')
        provenance = read(directory / 'provenance.json')
        assert done['completed'] and audit['checks_passed'] and provenance['exit_code'] == 0
        assert audit['source_sha256'] == provenance['source_sha256'] == source_hash
        assert audit['config_sha256'] == provenance['config_sha256'] == digest('configs/' + case + '.cfg')
        analyses[case] = read(directory / 'analysis.json')
        completed[case], audits[case] = done, audit

    panels = {}
    for name in sorted(PANELS, key=lambda p: (p[0] == 'S', p)):
        entry = PANELS[name]
        sidecar = Path('data/panels', name + '.json')
        if set(entry['cases']).issubset(completed) and sidecar.exists():
            panel = read(sidecar)
            assert panel['model_source_sha256'] == source_hash
            for case in entry['cases'] if len(entry['cases']) == 1 else []:
                assert panel['source_provenance']['config_sha256'] == audits[case]['config_sha256']
            for path in [Path(entry['data_file']), Path('figures/panels', name + '.pdf')]:
                inputs[str(path)] = digest(path)
            panels[name] = panel
    counts = Counter(p['status'] for p in panels.values())
    pending = [p for p in PANELS if p not in panels]
    missing_cases = [c for c in MAIN_CASES + VALIDATION_CASES if c not in completed]

    for path in ['report/reproduction.tex', 'scripts/inventory.py', 'scripts/report_tables.py']:
        inputs[path] = digest(path)
    source = Path('report/reproduction.tex').read_text()
    preamble = source.split(r'\title{', 1)[0]
    bibliography = source[source.index(r'\begin{thebibliography}'):source.index(r'\end{document}')]
    model = source[source.index(r'\section{Physical model'):source.index(r'\input{verification.tex}')]
    # The unchanged full-report methods are shared; all interim results are separate.
    model = model.replace('reported below.', 'recorded in the completed-case audit files.')
    out = preamble + r'''
\pagestyle{myheadings}
\markright{INTERIM REPORT: incomplete reproduction study}
\title{Independent calculations of fault valving: interim results}
\author{Reproducibility study of Zhu, Allison, Dunham, and Yang (2020)}
\date{''' + stamp + r'''}
\begin{document}
\maketitle
\textbf{This is an interim scientific report. The full study remains incomplete.}
'''
    out += ('At this snapshot, {} of 12 configured cases have completed, and {} of 49 panel packages '
            'are available. Of these, {} are independently reproduced constitutive-law panels and {} '
            'are partial comparisons, including one newly drawn schematic. ').format(
                len(completed), len(panels), counts['independently reproduced'], counts['partial'])
    if counts['not reproduced']:
        out += '{} generated attempts are classified as not reproduced. '.format(counts['not reproduced'])
    if pending:
        out += ('The {} pending panels are labeled as not reproduced (pending) in the inventory; '
                'this describes their current availability, not a failed physical behavior. ').format(len(pending))
    out += ('Pending runs are {}. Their configured durations remain unchanged.\n').format(
        ', '.join(command(c) for c in missing_cases) if missing_cases else 'none')
    out += r'''
\section{Scope and independence}
Fault slip increases permeability, while healing between earthquakes restricts
drainage. The resulting pore-pressure changes alter frictional strength. We
independently implemented the coupled earthquake-cycle and fluid-flow model of
Zhu et al.~\cite{zhu} in C++. The scientific inputs are the published article,
its supplement, and comparable written model documentation. The earlier project
attempt used authors' archived outputs; that attempt and all its simulation
inputs are excluded here. No authors' code, simulation output, source dataset,
digitized output curve, or published image supplies a numerical result or figure
in this report. The schematic is newly drawn and is not a reproduced simulation.

Every displayed numerical panel comes from independent arrays in
\code{data/panels}. Its JSON sidecar records parameters, source hashes, selection
window, and provenance. A constitutive-law check establishes a narrower result
than recovery of an earthquake sequence. Displayed cycle calculations retain
the status in each inventory row because timing, event patterns, assumed inputs,
or spatial convergence remain unresolved. Published cycle durations below are approximate
visual readings used only as comparison targets, never computational inputs.
'''
    out += model
    out += r'''
\section{Completed runs and quantitative comparisons}
Every completed case passed an audit of all accepted history rows and saved
field profiles, including samples outside the figure windows. The audit checks
finite values, positive effective stress, permeability bounds, accepted error
estimates, consecutive steps, field/history consistency, and source and parameter
hashes. Fixed-pressure cases also retain pressure throughout their saved profiles.
These checks establish internal consistency; they do not establish sequence
convergence or agreement with the paper. The audit files are
\code{data/<case>/output_verification.json}.

\begin{center}\small
\begin{tabular}{lrrrr}\toprule
Case & Duration (yr) & Accepted steps & Profiles & Large events\\\midrule
'''
    for case in completed:
        a, audit = analyses[case], audits[case]
        out += ' & '.join([tex(case), number(audit['years']), str(audit['accepted_steps']),
            str(audit['saved_profiles']), str(a['large_event_count'])]) + r'\\' + '\n'
    out += r'''\bottomrule\end{tabular}\end{center}
A seismic interval is a whole-fault maximum speed excursion above
$10^{-3}$ m s$^{-1}$. A resolved rupture accumulates at least 1 cm of slip in a
connected depth interval. A large event has a connected footprint reaching above
2 km and spanning more than 10 km. Slip is measured between the saved first-above
and first-below threshold profiles. Incomplete events cannot define cycles.
Velocity maps use the last complete large-event cycle; cumulative-slip profiles
use the last two complete cycles. All available recurrence intervals are retained
in each analysis file. Neither event selection nor time windows are fitted to
the published images.

\begin{center}\small
\begin{tabular}{lrrrr}\toprule
Case & Last cycle (yr) & All-cycle median (yr) & Paper (yr) & Difference (\%)\\\midrule
'''
    cycle_metrics = {}
    for case, (target, document) in PUBLISHED_CYCLES.items():
        if case not in analyses:
            continue
        a = analyses[case]
        last = (a['selected_cycle_s'][1] - a['selected_cycle_s'][0]) / YEAR
        difference = 100 * (last / target - 1)
        cycle_metrics[case] = dict(last_years=last, median_years=a['median_recurrence_years'],
            approximate_target_years=target, approximate_difference_percent=difference, source=document)
        out += ' & '.join([tex(case), number(last), number(a['median_recurrence_years']),
                          number(target), number(difference)]) + r'\\' + '\n'
    out += r'''\bottomrule\end{tabular}\end{center}
The reference and baseline targets are main Figs. 2a and 2c; the long-healing targets are supplementary
Figs. 3a,c and 5a,c. The differences describe the selected final cycles and do not
imply stationary recurrence. In particular, the two long-healing fixed-pressure
references differ substantially from their approximately 160-year targets.
Uncertain friction profiles and initialization, finite simulated duration, and
unestablished spatial convergence prevent assigning these discrepancies to one
cause. No high-healing-time mesh study has been completed.

\begin{center}\small
\begin{tabular}{lrrrr}\toprule
Case & $T$ (s) & $\Delta(\sigma-p)$ at 10 km (MPa) & Leading front & Deepest front\\
 & & & (m/yr) & (m/yr)\\\midrule
'''
    for case in ['short', 'baseline', 'long', 'verylong']:
        if case not in analyses:
            continue
        a = analyses[case]
        out += ' & '.join([tex(case), number(a['configuration']['T']),
            number(a['sample_depths']['10000']['effective_range_mpa']),
            number(front_speed(a, 'leading_migration_segments')),
            number(front_speed(a, 'migration_segments'))]) + r'\\' + '\n'
    out += r'''\bottomrule\end{tabular}\end{center}
The stress range uses the nearest saved depth, about 10.004 km. Front rates are
medians of all retained upward segments of the $|V|=V_p$ contour, with $R^2\geq0.8$,
at 13--20 km depth. The leading and deepest branches are contour extrema, not
tracked identities of individual pulses. Their definitions need not coincide
with the paper's reported migration speeds of 2500, 380, 120, and 30 m/yr for the
short, baseline, long, and very-long healing cases. The discrepancies are retained rather than selecting the branch
closest to the published rate.
'''
    if 'baseline' in analyses:
        a = analyses['baseline']
        sample = a['sample_depths']['10000']
        phase = a['phase_diagnostics']
        closing = [e for e in a['events'] if e['large'] and e['complete']][-1]
        top, bottom = max(closing['rupture_intervals_m'], key=lambda pair: pair[1]-pair[0])
        small = panels['F4a']['middepth_small_events'] if 'F4a' in panels else []
        out += ('\nThe baseline recurrence intervals range from {} to {} yr; the close agreement '
                'of its final interval with the displayed 32-year example does not establish '
                'a stationary cycle. At the saved depth nearest 10 km, effective stress spans '
                '{}--{} MPa, a variation of {} MPa compared with the article\'s 10--20 MPa scale. '
                'The maximum occurs {} yr after cycle opening. This maximum is a diagnostic, '
                'not an identification of the paper\'s phase boundary. Maximum upward flux '
                'there is {} m/s, compared with the article\'s order $10^{{-7}}$ m/s upper scale.\n').format(
            number(min(a['recurrence_years'])), number(max(a['recurrence_years'])),
            number(sample['effective_min_mpa']), number(sample['effective_max_mpa']),
            number(sample['effective_range_mpa']),
            number(phase['maximum_effective_stress_since_window_start_years']), number(sample['flux_max']))
        first = phase['first_partial_ruptures']
        out += '\nThe first two resolved partial ruptures are retained in chronological order: '
        out += '; '.join('{} yr after opening, with a connected footprint at {}--{} km'.format(
            number(e['peak_since_window_start_years']), number(e['footprint_top_m']/1000),
            number(e['footprint_bottom_m']/1000)) for e in first) + '. '
        out += ('The extra early shallow rupture differs from the published sequence. '
                'The final two-year window contains {} smaller ruptures intersecting 2--10 km. '
                'This recovers swarm-like behavior, while its timing and depth range differ. '
                'The closing large event first crosses the speed threshold at {} km, '
                'and its connected footprint spans {}--{} km. The paper gives a 12 km '
                'nucleation depth for the preceding 8--19 km rupture; it does not specify '
                'a numerical nucleation depth for the final surface-breaking event in this description. '
                'Our threshold definition does not establish identical published event boundaries.\n').format(
            len(small), number(closing['nucleation_depth_m']/1000),
            number(top/1000), number(bottom/1000))
        supplementary, _ = baseline_supplementary(a, panels)
        out += '\n' + supplementary
    if 'short' in analyses:
        out += r'''
The short-healing case produces aseismic pulses. The following medians compare
with the article's approximate one-year duration, few-year recurrence, and
few-centimeter slip scales. Pulse widths depend on our declared
threshold, $|V|\geq1.1V_p$. Each counted episode lasts at least 0.01 yr, lies
within the selected window, and remains below the seismic threshold at every
accepted whole-fault history sample. All qualifying episodes are retained in
the case analysis. The depth dependence limits any claim of agreement.
\begin{center}\small
\begin{tabular}{rrrrr}\toprule
Depth (km) & Episodes & Duration (yr) & Recurrence (yr) & Slip (cm)\\\midrule
'''
        for key in ['15000', '18000', '20000']:
            s = analyses['short']['slow_slip'][key]
            out += ' & '.join([number(s['sample_depth_m'] / 1000), str(len(s['episodes'])),
                number(s['median_duration_years']), number(s['median_interval_years']),
                number(100 * s['median_net_slip_m'])]) + r'\\' + '\n'
        out += r'\bottomrule\end{tabular}\end{center}' + '\n'
    out += r'''
The final model assessment requires the baseline, mesh-refined case, and tighter
time-tolerance case together. This interim report does not treat completed
coarse-grid, normal-stress, or initialization runs as a completed sensitivity
comparison. Software fixtures and exploratory short calculations are excluded
from all scientific figures here.

\section{Panel inventory and provenance}
Each row specifies the original label, quantity, documentation, independent-data
command, plotting command, status, and discrepancy. Data paths follow
\code{data/panels/<label>.npz}, with parameters and provenance in the adjacent
JSON file. Commands are intended for a clean checkout or completed cases; do not
start a duplicate simulation in a directory containing a live run.
Supplementary Fig. 2 visibly labels rows a and b; its caption also mentions c
and d, which are not separately visible. The inventory therefore includes the
two displayed rows and records that inconsistency.
\begin{landscape}\scriptsize
\setlength{\tabcolsep}{3pt}
\begin{longtable}{P{10mm}P{30mm}P{34mm}P{39mm}P{32mm}P{23mm}P{65mm}}
\toprule Panel & Quantity & Documentation & Data command & Plot command & Status & Discrepancy\\\midrule\endhead
'''
    for name in sorted(PANELS, key=lambda p: (p[0] == 'S', p)):
        entry = PANELS[name]
        p = panels.get(name, entry)
        fields = [name, tex(p['quantity']), tex(p['documentation']), command(p['data_command']),
                  command(p['plot_command']), tex(p.get('status', 'not reproduced (pending)')),
                  tex(p.get('main_discrepancy', 'Pending completed baseline calculation and final extraction; no result supplied in this snapshot.'))]
        out += ' & '.join(fields) + r'\\[5pt]' + '\n'
    out += r'\bottomrule\end{longtable}\end{landscape}'
    out += '\n'
    figures = []
    for group in sorted(GROUPS, key=lambda g: (g[0] == 'S', g)):
        names = GROUPS[group]
        available = [p for p in names if p in panels]
        if not available:
            continue
        if len(available) == len(names):
            paths = [Path('figures', group + '.pdf')]
        else:
            paths = [Path('figures/panels', p + '.pdf') for p in available]
        description = '; '.join(p + ': ' + PANELS[p]['quantity'] for p in available)
        missing = [p for p in names if p not in panels]
        caption = ('Independent outputs corresponding to ' + group + '. ' + description +
                   '. Units and color scales are shown on the panels. Parameter files and selection rules are specified in the provenance sidecars.')
        if missing:
            caption += ' Pending and absent here: ' + ', '.join(missing) + '.'
        if group == 'F1':
            caption += ' Panel F1a is a new schematic, not a simulated reproduction.'
        else:
            cases = sorted(set(c for p in available for c in PANELS[p]['cases']))
            settings = []
            for case in cases:
                cfg = analyses[case]['configuration']
                settings.append(case + ' (T=' + number(cfg['T']) + ' s; ' +
                    ('fixed pressure' if cfg['fixed_pressure'] else 'coupled pressure') + ')')
            caption += ' Cases: ' + ', '.join(settings) + '. These panels are partial comparisons.'
        out += r'\clearpage\begin{figure}[p]\centering' + '\n'
        for path in paths:
            figures.append(path)
            inputs[str(path)] = digest(path)
            out += r'\includegraphics[width=\textwidth,height=' + ('.34' if len(paths) > 1 else '.72') + r'\textheight,keepaspectratio]{../' + str(path) + '}\n'
        out += r'\caption{' + tex(caption) + r' Comparison target: Zhu et al.~\cite{zhu}.}\end{figure}' + '\n'
    out += r'''
\clearpage
\section{Reproduction and remaining work}
The full project command builds and tests the C++ solver, runs all twelve cases,
extracts every panel, compiles the final report, and verifies its source bundle:
\begin{verbatim}
python3 -m pip install -r requirements.txt
python3 scripts/reproduce.py
\end{verbatim}
That full-duration run is ongoing at this snapshot. This interim report is
regenerated from completed, audited case metadata and generated panel packages:
\begin{verbatim}
python3 scripts/interim_report.py
\end{verbatim}
The command compiles \code{report/interim.pdf}, packages its standalone LaTeX
source and all included figures, and verifies compilation of the extracted bundle.
It records commands, exit codes, input hashes, and output hashes in
\code{data/interim_report_verification.json}. The source archive is
\code{report/interim-latex-source.zip}; the unexpanded source is already a single
file, \code{report/interim.tex}. All scientific arrays and their independent
generation commands are available in the repository.

Completion still requires the remaining calculations, any pending
panels listed in the inventory, common-duration validation comparisons, the full scientific report,
visual and numerical artifact review, execution of the final regeneration
procedure, and publication of those final deliverables. This snapshot makes
the completed scientific results reviewable without declaring that work finished.
'''
    out += bibliography + '\n' + r'\end{document}' + '\n'
    Path('report/interim.tex').write_text(out)
    commands = []
    with Path('data/interim_report_build.log').open('w') as log:
        def run(args, cwd):
            start = time.time()
            result = subprocess.run(args, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT)
            commands.append(dict(command=args, cwd=str(cwd), exit_code=result.returncode,
                                 elapsed_s=time.time() - start))
            assert result.returncode == 0, args
        args = ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', 'interim.tex']
        run(args, ROOT / 'report')
        with zipfile.ZipFile('report/interim-latex-source.zip', 'w', zipfile.ZIP_DEFLATED) as bundle:
            for path in [Path('report/interim.tex')] + figures:
                bundle.write(str(path), str(path))
        with tempfile.TemporaryDirectory(prefix='interim-compile-', dir=str(scratch)) as directory:
            with zipfile.ZipFile('report/interim-latex-source.zip') as bundle:
                bundle.extractall(directory)
            run(args, Path(directory) / 'report')
        run(['pdftotext', 'interim.pdf', str(scratch / 'interim-report.txt')], ROOT / 'report')
    rendered = (scratch / 'interim-report.txt').read_text()
    assert 'interim' in rendered.lower() and 'incomplete' in rendered.lower()
    assert all(name in rendered for name in PANELS)
    assert 'Overfull' not in Path('report/interim.log').read_text()
    outputs = {str(p): digest(p) for p in [Path('report/interim.tex'), Path('report/interim.pdf'),
        Path('report/interim-latex-source.zip'), Path('data/interim_report_build.log')]}
    record = dict(scope='Dated interim report only; not completion of the independent reproduction goal.',
        recorded_utc=stamp, checks_passed=True, completed_cases=list(completed), pending_cases=missing_cases,
        included_panels=list(panels), pending_panels=pending, included_status_counts=dict(counts),
        cycle_comparisons=cycle_metrics, commands=commands, input_sha256=inputs,
        script_sha256=digest(__file__), output_sha256=outputs,
        checks=['Completed cases have successful exit and raw-output audits with matching source and configuration hashes.',
                'All 49 panel labels and explicit incompleteness appear in the extracted PDF text.',
                'Both working-tree and isolated source-bundle compilations succeed.',
                'LaTeX reports no overfull boxes.'])
    Path('data/interim_report_verification.json').write_text(json.dumps(record, indent=2) + '\n')
    print('Compiled interim report: {} cases, {} panels; {} panels pending.'.format(
        len(completed), len(panels), len(pending)))


if __name__ == '__main__':
    main()
