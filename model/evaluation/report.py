"""Render the same requirement records for the candidate index, comparison and detailed ledger."""
from __future__ import annotations

import os
from pathlib import Path

from .assessment import ROOT, SHORT, STATUS, DECISIONS, comparison_studies

HERE = Path(__file__).resolve().parent
CHECKS = {
    'native_cell': 'Native cell has the fixed 12-port, 30-seam incidence',
    'hermitian': 'Cell generator is Hermitian',
    'locality': 'Every cell coupling joins incident registers',
    'rotations': 'All 60 proper rotations preserve the generator',
    'rank_three_response': 'Lowest nonconstant cell response has multiplicity three',
    'federation': 'Registered generator covers native 20-, 80- and 320-cell boards',
    'parameters': 'Declared parameters pass provenance and forbidden-input checks',
}
GRADES = {**STATUS, 'pass': 'Pass', 'fail': 'Fail', 'fulfilled': 'Demonstrated', 'fails': 'Counterexample'}


def local_link(source, directory=HERE):
    if source.startswith('https://'):
        return source
    path, marker, fragment = source.partition('#')
    return Path(os.path.relpath(ROOT / path, directory)).as_posix() + marker + fragment


def links(sources, directory=HERE):
    out = []
    for source in dict.fromkeys(sources):
        path, _, fragment = source.partition('#')
        name = Path(path).name
        label = name + (' ' + fragment if fragment.startswith('L') else '')
        out.append(f'[{label}]({local_link(source, directory)})')
    return ' · '.join(out)


def table_row(*cells):
    return '| ' + ' | '.join(str(c).replace('|', r'\|').replace('\n', ' ') for c in cells) + ' |'


def table(*columns):
    return [table_row(*columns), table_row(*['---' for _ in columns])]


def criterion_grade(case):
    grade = STATUS[case['status']]
    measurement = case.get('demonstrated_measurement')
    if measurement is not None:
        source = local_link(measurement['source'])
        linked = f"**Demonstrated** · [{measurement['label']}]({source})"
        grade = linked if case['status'] == 'demonstrated' else f"{linked}. Full criterion · {grade}"
    filing = case.get('filing')
    if filing:
        grade += (f". OpenWave · [{filing['official_status']}]({local_link(filing['decision'])})"
                  f". Resubmission · [{filing['resubmission']}]({local_link(filing['review'])})")
    return grade


def comparison_link(comparison):
    law = comparison['law']
    label = f"{law} · {comparison['label']}" if 'label' in comparison else law
    return f"[{label}]({local_link(comparison['source'])})"


def short(result):
    return SHORT.get(result['candidate'], result['title'])


def measured(result, directory):
    m = result['measurement']
    if m is None:
        return 'Awaiting a registered readout.'
    number = format(m['value'], '.6g')
    if m['standard_error'] is not None:
        number += ' ± ' + format(m['standard_error'], '.2g') + ' (SE)'
    return f"{m['label']} **{number}**. [Data]({local_link(m['source'], directory)})"


def candidate_table(results, directory):
    lines = table('Candidate', 'Role', 'Measured example', 'Reason for its role')
    for r in results:
        doc = f"model/candidate_laws/{r['candidate']}/README.md"
        name = f"[{short(r)}]({local_link(doc, directory)})"
        lines.append(table_row(name, r['role'], measured(r, directory), DECISIONS.get(r['candidate'], r['summary'])))
    return '\n'.join(lines)


def grouped(results):
    for i, group in enumerate(results[0]['requirements']):
        yield group, list(zip(*(r['requirements'][i]['criteria'] for r in results), strict=True))


def scorecard(results):
    names = [short(r) for r in results]
    lines = ['# Candidate law scorecard', '',
             'Each row asks one question. Select its ID to see the acceptance criterion, the evidence for each law and the scope of earlier results.', '',
             candidate_table(results, HERE), '',
             '**Demonstrated** means the evidence meets the stated criterion on its declared domain. **Partial** identifies a measured or proved part and a remaining condition. **Open** means the available evidence leaves the question undecided. **Counterexample** identifies a tested construction or claim contradicted in its stated domain.', '',
             'Roles reflect the current investigation of patch clocks and bound matter on the spherical federation. Several rows use the same experiment at different scopes. Adding their grades would count the same evidence several times.', '',
             'The signed-defect experiment tests whether a transported orbital qubit acquires spinorial rotation and exchange phases. The gravity experiment tests how a second probe and an equal-energy seam source respond to a calibration on ordinary matter. The semiclassical comparison uses supplied Einstein evolution. [Matter and its controls](EVIDENCE.md#ft-matter) · [Universal coupling and its controls](EVIDENCE.md#g6) · [Null balance](EVIDENCE.md#gr-null)', '',
             '## Common carrier and parameter checks', '',
             'These checks execute on each registered generator and its declared parameters at every evaluation. The cell spectral test concerns the finite response operator. Observer reconstruction is scored separately below.', '',
             *table('What is checked', *names)]
    for key, title in CHECKS.items():
        lines.append(table_row(title, *(GRADES.get(r['checks'].get(key, {}).get('status'), 'Unmeasured') for r in results)))
    for group, rows in grouped(results):
        lines.extend(['', '## ' + group['title'], '', group['purpose']])
        if group.get('reference'):
            lines.extend(['', group['reference']])
        if group.get('sources'):
            lines.extend(['', 'Sources · ' + links(group['sources'])])
        if group['id'] in {'sm-inventory', 'openwave'}:
            for i, result in enumerate(results):
                counts = {key: sum(cases[i]['status'] == key for cases in rows) for key in STATUS}
                nonzero = ', '.join(f"{count} {STATUS[key].lower()}" for key, count in counts.items() if count)
                lines.extend(['', f"**{short(result)}** · {nonzero} of {len(rows)}."])
        domain = None
        for i, cases in enumerate(rows):
            row = cases[0]
            if row.get('domain') != domain or i == 0:
                domain = row.get('domain')
                if domain:
                    lines.extend(['', '### ' + domain])
                headings = ['Requirement', *names]
                if group['id'] == 'observer_benchmarks':
                    headings.append('Historical apparatus and result')
                lines.extend(['', *table(*headings)])
            key = row['id']
            regime = ' (' + row['regime'] + ')' if row.get('regime') else ''
            cells = [f"[{key}](EVIDENCE.md#{key.lower()}) {row['title']}{regime}"]
            for case in cases:
                grade = criterion_grade(case)
                for study in comparison_studies(case):
                    grade += f". Comparison · {comparison_link(study)}"
                cells.append(grade)
            if group['id'] == 'observer_benchmarks':
                cells.append(row['context'])
            lines.append(table_row(*cells))
    lines.extend(['', '## Agreement across theorems and experiments', '',
                  'These requirements demand that the same physical system supply the quantities used in multiple arguments. Each [candidate receipt](../candidate_laws/CANDIDATE.md) records the declared scope and any counterexample.', '',
                  *table('Joint requirement', *names)])
    for key, row in results[0]['joint'].items():
        lines.append(table_row(f"[{row['requirement']}](EVIDENCE.md#joint-{key})", *(GRADES[r['joint'][key]['status']] for r in results)))
    lines.extend(['', 'The [full evidence ledger](EVIDENCE.md#full-premise-register) also reports every numbered premise, its owner and its status for all candidates.', ''])
    return '\n'.join(lines)


def evidence_ledger(results):
    lines = ['# Candidate evidence ledger', '', 'Use the [scorecard](SCORECARD.md) for the cross-candidate comparison. Each entry gives the observation or construction required for full credit. The scope in each verdict applies to that candidate.', '']
    for group, rows in grouped(results):
        lines.extend(['## ' + group['title'], '', group['purpose'], ''])
        for cases in rows:
            row = cases[0]
            lines.extend([f"### {row['id']}", '', f"**{row['title']}.** {row['requires']}", ''])
            if row.get('premises'):
                lines.extend(['Register crosswalk · ' + ', '.join(row['premises']), ''])
            if row.get('context'):
                lines.extend([row['context'], ''])
            lines.extend(['Criterion · ' + links(row['sources']), ''])
            if all(c['status'] == 'open' and c['note'] == row['note'] and c['evidence'] == row['evidence']
                   and not comparison_studies(c) for c in cases):
                detail = links(row['evidence'])
                lines.extend([f"**All three candidates · Open.** {row['note']} {detail}".rstrip() if len(results) == 3
                              else f"**All listed candidates · Open.** {row['note']} {detail}".rstrip(), ''])
            else:
                for r, c in zip(results, cases, strict=True):
                    heading = f"**{short(r)} · {STATUS[c['status']]}.**"
                    if 'demonstrated_measurement' in c:
                        heading = f"**{short(r)}** · {criterion_grade(c)}."
                    lines.append(f"- {heading} {c['note']} {links(c['evidence'])}".rstrip())
                    for comparison in comparison_studies(c):
                        lines.append(f"  - **{comparison['law']}** · {comparison['finding']} "
                                     f"[Comparison study]({local_link(comparison['source'])})")
                lines.append('')
    lines.extend(['## Agreement across theorems and experiments', ''])
    for key, definition in results[0]['joint'].items():
        lines.extend(['### JOINT-' + key, '', definition['requirement'], ''])
        for r in results:
            case = r['joint'][key]
            sources = [os.path.normpath(f"model/candidate_laws/{r['candidate']}/{s}") for s in case['evidence']]
            note = case['note'] or 'No joint realization is registered.'
            lines.append(f"- **{short(r)} · {GRADES[case['status']]}.** {note} {links(sources)}".rstrip())
        lines.append('')
    lines.extend(['## Full premise register', '',
                  'The register records each theorem assumption and its formal consumers. The role identifies who must supply it. Candidate declarations, including supplied apparatus and imports, are preserved in the individual JSON receipts. An open entry has no full discharge registered for that row.', '',
                  '[Formal premise register](../../docs/registers/PREMISE_REGISTER_V3.md) · [Role assignments](roles.json)', '',
                  *table('Premise', 'Role', *(short(r) for r in results))])
    for cases in zip(*(r['premises'] for r in results), strict=True):
        row = cases[0]
        lines.append(table_row(f"{row['id']} {row['name']}", row['role'], *(GRADES[c['status']] for c in cases)))
    return '\n'.join(lines) + '\n'
