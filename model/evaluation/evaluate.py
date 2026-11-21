"""Evaluate candidate laws against the fixed carrier and the premises of the finite theorems.

    python -m model evaluate                 # print the evaluation of every candidate
    python -m model evaluate --law NAME      # one candidate
    python -m model evaluate --write         # write EVALUATION.json, SCORECARD.md and the table of CANDIDATE.md
    python -m model evaluate --check         # fail when a written file differs from a fresh evaluation

The carrier checks run on each candidate's own generator: it must act on the fixed echosahedral cell, couple only
incident registers, commute with the sixty rotations of the cell, build on the federations, and keep the
rank-three slowest response of the carrier. The parameter check reads the declared parameters and rejects any that
would hardcode an observed outcome. The premise matrix reads every row of tracking/premise_register.json, assigns
it a role from roles.json and records the candidate's declared status with its evidence. Undeclared premises are
open.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REGISTER = ROOT / 'tracking' / 'premise_register.json'
ROLES = HERE / 'roles.json'
SCORECARD = HERE / 'SCORECARD.md'
INDEX = ROOT / 'model' / 'candidate_laws' / 'CANDIDATE.md'
BEGIN, END = '<!-- CANDIDATES:BEGIN -->', '<!-- CANDIDATES:END -->'
STATUSES = ('fulfilled', 'partial', 'open', 'fails')
PROVENANCE = ('designed', 'calibrated', 'derived')
FORBIDDEN = ('particle mass', 'mass ratio', 'mixing angle', 'spatial dimension', 'number of dimensions',
             'spacetime metric', 'lattice', 'gauge group')
JOINT = {
    'one-law': 'Every premise assigned to the law is fulfilled by one unchanged law on one process.',
    'einstein-together': 'The stress that observers infer from modular flow is the stress carried by the excitations '
                         'that propagate and interact under the law.',
    'matter-together': "Gauge currents act on those excitations through the law's operations; detector, generator, "
                       'clock and currents refer to the same object.',
    'patch-measurement': 'Finite patches measure both the current response and the stress of one bound excitation '
                         'with their own clocks and records.',
    'relativistic-matter': 'The rest energy and the inertia of every bound excitation follow one dictionary.',
}


def fail(message):
    raise SystemExit(f'evaluation error: {message}')


def dense(H):
    return H.toarray() if hasattr(H, 'toarray') else np.asarray(H, dtype=float)


def carrier_checks(cand):
    from model.executable_model import carrier as fixed
    from model.executable_model import incident, rotations
    seams = [tuple(sorted(map(int, e))) for e in cand.LABELLING]
    iso = fixed.isomorphism(fixed.CELL_SEAMS, seams) if len(seams) == 30 else None
    out = {'native_cell': {'status': 'pass' if iso is not None else 'fail',
                           'detail': 'the labelling is the fixed cell' if iso is not None else
                                     'the labelling is not the 12-port, 30-seam echosahedron'}}
    if iso is None:
        return out
    cell = fixed.cell()
    H = dense(cand.generator(cell))
    labels = cand.basis(cell)
    if H.shape != (len(labels), len(labels)):
        fail(f'{cand.NAME}: generator and basis disagree in size')
    out['hermitian'] = {'status': 'pass' if np.abs(H - H.T).max() < 1e-12 else 'fail',
                        'detail': f'largest asymmetry {np.abs(H - H.T).max():.1e}'}
    i, j = np.nonzero(np.abs(H) > 1e-12)
    far = sum(1 for a, b in zip(i, j) if not incident(labels[a], labels[b], seams))
    out['locality'] = {'status': 'pass' if far == 0 else 'fail',
                       'detail': f'{far} couplings between registers that share no seam' if far else
                                 'every coupling joins incident registers of one cell'}
    inverse = {v: k for k, v in enumerate(iso)}
    worst = 0.0
    for g in rotations():
        moved = [iso[g[inverse[p]]] for p in range(12)]
        P = cand.represent(moved, cell)
        worst = max(worst, float(np.abs(P @ H @ P.T - H).max()))
    out['rotations'] = {'status': 'pass' if worst < 1e-10 else 'fail',
                        'detail': f'commutes with the sixty rotations of the cell (largest defect {worst:.1e})'}
    levels = np.linalg.eigvalsh(dense(cand.response(cell)))
    groups = []
    for x in np.sort(levels):
        if groups and abs(x - groups[-1][0]) < 1e-8:
            groups[-1][1] += 1
        else:
            groups.append([float(x), 1])
    slow = groups[1][1] if len(groups) > 1 else 0
    out['rank_three_response'] = {'status': 'pass' if slow == 3 else 'fail',
                                  'detail': 'levels above the lowest: ' + ', '.join(
                                      f'{x - groups[0][0]:.6f} x{m}' for x, m in groups[1:5])}
    built = ['cell']
    note = ''
    for level in (0, 1, 2):
        try:
            Hf = cand.generator(fixed.federation(level))
        except NotImplementedError as exc:
            note = str(exc)
            break
        asym = abs(Hf - Hf.T).max() if hasattr(Hf, 'toarray') else np.abs(Hf - Hf.T).max()
        if asym > 1e-12:
            fail(f'{cand.NAME}: generator on federation level {level} is not Hermitian')
        built.append(fixed.NAMES[level])
    out['federation'] = {'status': 'pass' if len(built) == 4 else 'partial',
                         'detail': 'builds on ' + ', '.join(built) + (f' ({note})' if note else '')}
    return out


def parameter_check(cand):
    findings = []
    for name, p in cand.PARAMETERS.items():
        if p.get('provenance') not in PROVENANCE:
            findings.append(f'{name}: provenance must be one of {PROVENANCE}')
        if p.get('provenance') == 'calibrated' and not p.get('observation'):
            findings.append(f'{name}: a calibrated parameter names the observation it is calibrated on')
        text = (p.get('quantity', '') + ' ' + p.get('note', '')).lower()
        findings += [f'{name}: declares a {word}' for word in FORBIDDEN if word in text]
    return {'status': 'pass' if not findings else 'fail',
            'detail': 'no parameter hardcodes an observed outcome' if not findings else '; '.join(findings)}


def evidence_ok(cand, paths):
    base = ROOT / 'model' / 'candidate_laws' / cand.NAME
    missing = [p for p in paths if not (base / p).resolve().exists()]
    if missing:
        fail(f'{cand.NAME}: missing evidence {missing}')


def premise_rows(cand, register, roles):
    ids = {r['id'] for r in register['rows']}
    unknown = sorted(set(cand.PREMISES) - ids)
    if unknown:
        fail(f'{cand.NAME}: premises not in the register: {unknown}')
    rows = []
    for r in register['rows']:
        role = roles['overrides'].get(r['id'], {}).get('role') or roles['defaults'][r['type']]
        decl = cand.PREMISES.get(r['id'], {})
        status = decl.get('status', 'open')
        if status not in STATUSES:
            fail(f'{cand.NAME} {r["id"]}: unknown status {status!r}')
        if status == 'fails' and not decl.get('counterexample'):
            fail(f'{cand.NAME} {r["id"]}: a failed premise names its counterexample')
        evidence_ok(cand, decl.get('evidence', []))
        rows.append({'id': r['id'], 'name': r['name'], 'type': r['type'], 'role': role, 'status': status,
                     'evidence': decl.get('evidence', []), 'note': decl.get('note', ''),
                     **({'counterexample': decl['counterexample']} if 'counterexample' in decl else {})})
    return rows


def joint_rows(cand):
    unknown = sorted(set(cand.JOINT) - set(JOINT))
    if unknown:
        fail(f'{cand.NAME}: unknown joint requirements {unknown}')
    out = {}
    for key, text in JOINT.items():
        decl = cand.JOINT.get(key, {})
        status = decl.get('status', 'open')
        if status not in STATUSES or (status == 'fails' and not decl.get('counterexample')):
            fail(f'{cand.NAME} {key}: invalid declaration')
        evidence_ok(cand, decl.get('evidence', []))
        out[key] = {'requirement': text, 'status': status, 'evidence': decl.get('evidence', []),
                    'note': decl.get('note', ''),
                    **({'counterexample': decl['counterexample']} if 'counterexample' in decl else {})}
    return out


def evaluate(name):
    from model.executable_model import laws
    cand = laws.load(name)
    register = json.loads(REGISTER.read_text(encoding='utf-8'))
    roles = json.loads(ROLES.read_text(encoding='utf-8'))
    checks = carrier_checks(cand)
    checks['parameters'] = parameter_check(cand)
    premises = premise_rows(cand, register, roles)
    summary = {role: {s: sum(1 for r in premises if r['role'] == role and r['status'] == s) for s in STATUSES}
               for role in ('law', 'designer', 'import')}
    return {'schema': 'oph.law_evaluation.v1', 'candidate': cand.NAME, 'title': cand.TITLE, 'role': cand.ROLE,
            'summary': cand.SUMMARY, 'parameters': cand.PARAMETERS, 'checks': checks, 'joint': joint_rows(cand),
            'premise_counts': summary, 'premises': premises}


def dumps(value):
    return json.dumps(value, indent=1, ensure_ascii=False) + '\n'


def cell(status):
    return {'pass': 'pass', 'fail': 'FAIL', 'partial': 'partial'}.get(status, status)


def scorecard(results):
    names = [r['candidate'] for r in results]
    lines = ['# Scorecard of the candidate laws', '',
             'Generated by `python -m model evaluate --write`; do not edit by hand. The criteria are in',
             '[REQUIREMENTS.md](../REQUIREMENTS.md): the law runs on the fixed spherical echosahedral federation and',
             'fulfils the premises of the finite theorems together. Calibration of anything a designer can implement',
             'is allowed; hardcoding an observed outcome is not.', '',
             '## Carrier and parameter checks', '',
             '| Check | ' + ' | '.join(names) + ' |', '| --- |' + ' --- |' * len(names)]
    keys = ['native_cell', 'hermitian', 'locality', 'rotations', 'rank_three_response', 'federation', 'parameters']
    for k in keys:
        lines.append(f'| {k.replace("_", " ")} | ' + ' | '.join(
            f"{cell(r['checks'][k]['status'])}: {r['checks'][k]['detail']}" if k in r['checks'] else 'not run'
            for r in results) + ' |')
    lines += ['', '## Premises of the finite theorems', '',
              'Counts over the rows of the premise register. A premise with role "law" is supplied by the law, one '
              'with role "designer" may be calibrated, and an import is supplied from outside.', '',
              '| Role | ' + ' | '.join(names) + ' |', '| --- |' + ' --- |' * len(names)]
    for role in ('law', 'designer', 'import'):
        lines.append(f'| {role} | ' + ' | '.join(
            ', '.join(f"{r['premise_counts'][role][s]} {s}" for s in STATUSES) for r in results) + ' |')
    lines += ['', '## Joint requirements', '', '| Requirement | ' + ' | '.join(names) + ' |',
              '| --- |' + ' --- |' * len(names)]
    for key, text in JOINT.items():
        lines.append(f'| {text} | ' + ' | '.join(r['joint'][key]['status'] for r in results) + ' |')
    lines += ['', '## Premises with evidence', '', '| Premise | Role | ' + ' | '.join(names) + ' |',
              '| --- | --- |' + ' --- |' * len(names)]
    for i, row in enumerate(results[0]['premises']):
        cells = [r['premises'][i]['status'] for r in results]
        if any(c != 'open' for c in cells):
            lines.append(f"| {row['id']} {row['name']} | {row['role']} | " + ' | '.join(cells) + ' |')
    lines += ['', 'Every other premise of the register is open for every candidate. The evaluation of each candidate',
              'with its notes and evidence is its `EVALUATION.json`.', '']
    return '\n'.join(lines)


def index_table(results):
    lines = [BEGIN, '| Candidate | Role | Carrier checks | Law premises | Joint requirements |',
             '| --- | --- | --- | --- | --- |']
    for r in results:
        checks = sum(1 for c in r['checks'].values() if c['status'] == 'pass')
        law = r['premise_counts']['law']
        joint = {s: sum(1 for j in r['joint'].values() if j['status'] == s) for s in STATUSES}
        lines.append(f"| [{r['title']}]({r['candidate']}/README.md) | {r['role']} | {checks} of {len(r['checks'])} "
                     f"pass | {law['fulfilled']} fulfilled, {law['partial']} partial, {law['open']} open, "
                     f"{law['fails']} failed | {joint['fulfilled']} fulfilled, {joint['partial']} partial, "
                     f"{joint['open']} open, {joint['fails']} failed |")
    lines.append(END)
    return '\n'.join(lines)


def outputs(results):
    files = {ROOT / 'model' / 'candidate_laws' / r['candidate'] / 'EVALUATION.json': dumps(r) for r in results}
    files[SCORECARD] = scorecard(results)
    text = INDEX.read_text(encoding='utf-8')
    if BEGIN not in text or END not in text:
        fail('CANDIDATE.md has no generated table markers')
    head, rest = text.split(BEGIN, 1)
    files[INDEX] = head + index_table(results) + rest.split(END, 1)[1]
    return files


def main(argv=None):
    from model.executable_model import laws
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--law', choices=laws.names())
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = ap.parse_args(argv)
    if (args.write or args.check) and args.law:
        ap.error('--write and --check evaluate every candidate')
    order = sorted(laws.names(), key=lambda n: (laws.load(n).ROLE != 'main', n))
    results = [evaluate(n) for n in ([args.law] if args.law else order)]
    if args.write:
        for path, text in outputs(results).items():
            path.write_text(text, encoding='utf-8')
            print('wrote', path.relative_to(ROOT))
        return 0
    if args.check:
        stale = [str(p.relative_to(ROOT)) for p, t in outputs(results).items()
                 if not p.exists() or p.read_text(encoding='utf-8') != t]
        if stale:
            print('stale evaluation files: ' + ', '.join(stale) + '; run python -m model evaluate --write')
            return 1
        print('evaluation files are current')
        return 0
    for r in results:
        law = r['premise_counts']['law']
        print(f"{r['title']} ({r['role']} candidate)")
        for k, c in r['checks'].items():
            print(f"  {k:22s} {c['status']:8s} {c['detail']}")
        print(f"  law premises: {law['fulfilled']} fulfilled, {law['partial']} partial, {law['open']} open, "
              f"{law['fails']} failed")
        for key, j in r['joint'].items():
            print(f"  joint {key:18s} {j['status']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
