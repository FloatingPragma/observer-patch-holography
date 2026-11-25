"""Compare registered candidate laws on the spherical carrier and the scientific requirements.

Carrier and parameter checks run live. Evidence-backed assessments retain the domain of each theorem and
experiment. Use --criterion ID for one question, --write to generate the comparison, or --check to verify it.
"""
from __future__ import annotations

import argparse
import json
import operator
import sys
from pathlib import Path

import numpy as np

from . import assessment, report

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REGISTER = ROOT / 'docs' / 'registers' / 'OBSERVABLES.json'
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
    return H.toarray() if hasattr(H, 'toarray') else np.asarray(H)


def register_at(label, port):
    """Map the ports of a register. A wedge has an unoriented support and a fixed middle port."""
    kind, where = label
    if kind == 'port':
        return kind, port(where)
    if kind == 'seam':
        a, b = map(port, where)
        return kind, tuple(sorted((a, b)))
    if kind == 'wedge':
        a, b, c = map(port, where)
        return kind, (min(a, c), b, max(a, c))
    raise ValueError(f'unknown register kind {kind!r}')


def native_basis(labels, port_count, seamset):
    try:
        keys = {register_at(label, operator.index) for label in labels}
    except (TypeError, ValueError):
        return set()
    if not keys or len(keys) != len(labels):
        return set()
    for kind, where in keys:
        if kind == 'port' and not 0 <= where < port_count:
            return set()
        if kind == 'seam' and where not in seamset:
            return set()
        if kind == 'wedge':
            a, b, c = where
            if a == c or any(tuple(sorted(edge)) not in seamset for edge in ((a, b), (b, c))):
                return set()
    return keys


def generator_defects(H, labels, seamset):
    from scipy import sparse
    from model.executable_model import incident
    H = H.tocsr() if sparse.issparse(H) else np.asarray(H)
    data = H.data if sparse.issparse(H) else H
    if not np.isfinite(data).all():
        return float('inf'), None
    asym = float(abs(H - H.conj().T).max())
    i, j, values = sparse.find(H)
    far = sum(1 for a, b, v in zip(i, j, values)
              if abs(v) > 1e-12 and not incident(labels[a], labels[b], seamset))
    return asym, far


def carrier_checks(cand):
    from model.executable_model import carrier as fixed
    from model.executable_model import rotations
    seams = {tuple(sorted(map(int, e))) for e in cand.LABELLING}
    iso = fixed.isomorphism(fixed.CELL_SEAMS, seams) if len(cand.LABELLING) == len(seams) == 30 and all(
        0 <= a < b < 12 for a, b in seams) else None
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
    cell_keys = native_basis(labels, cell.n, seams)
    if not cell_keys:
        out['locality'] = {'status': 'fail', 'detail': 'the basis contains missing, duplicate or nonnative registers'}
        return out
    asym, far = generator_defects(H, labels, seams)
    out['hermitian'] = {'status': 'pass' if asym < 1e-12 else 'fail',
                        'detail': f'largest asymmetry {asym:.1e}'}
    if far is None:
        return out
    out['locality'] = {'status': 'pass' if far == 0 else 'fail',
                       'detail': f'{far} couplings between registers that share no seam' if far else
                                 'every coupling joins incident registers of one cell'}
    inverse = {v: k for k, v in enumerate(iso)}
    index = {register_at(label, int): k for k, label in enumerate(labels)}
    worst = 0.0
    for g in rotations():
        moved = [iso[g[inverse[p]]] for p in range(12)]
        P = dense(cand.represent(moved, cell))
        targets = [index.get(register_at(label, moved.__getitem__)) for label in labels]
        if None in targets or P.shape != H.shape or not np.isfinite(P).all():
            worst = float('inf')
            break
        support = np.zeros(P.shape)
        support[targets, np.arange(len(labels))] = 1
        worst = max(worst, float(np.abs(abs(P) - support).max()),
                    float(np.abs(P @ H @ P.conj().T - H).max()))
    out['rotations'] = {'status': 'pass' if worst < 1e-10 else 'fail',
                        'detail': f'commutes with the sixty rotations of the cell (largest defect {worst:.1e})'}
    out['rank_three_response'] = rank_check(cand.response(cell))
    out['federation'] = federation_check(cand, cell_keys, iso, fixed)
    return out


def rank_check(response):
    matrix = dense(response)
    if (matrix.ndim != 2 or not matrix.shape[0] or matrix.shape[0] != matrix.shape[1]
            or not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.conj().T, atol=1e-12, rtol=0)):
        return {'status': 'fail', 'detail': 'the declared cell response must be finite and Hermitian'}
    levels = np.linalg.eigvalsh(matrix)
    groups = []
    for x in np.sort(levels):
        if groups and abs(x - groups[-1][0]) < 1e-8:
            groups[-1][1] += 1
        else:
            groups.append([float(x), 1])
    slow = groups[1][1] if len(groups) > 1 else 0
    return {'status': 'pass' if slow == 3 else 'fail',
            'detail': 'levels above the lowest: ' + ', '.join(
                f'{x - groups[0][0]:.6f} x{m}' for x, m in groups[1:5])}


def federation_check(cand, cell_keys, iso, fixed):
    def relabel(p):
        return 12 * (p // 12) + iso[p % 12]

    built = ['cell']
    status, note = 'pass', ''
    for level in (0, 1, 2):
        fed = fixed.federation(level)
        try:
            Hf = cand.generator(fed)
        except NotImplementedError as exc:
            status, note = 'partial', str(exc)
            break
        fed_labels = cand.basis(fed)
        fed_seams = {tuple(sorted((relabel(a), relabel(b)))) for a, b in fed.seams}
        keys = native_basis(fed_labels, fed.n, fed_seams)
        expected = {register_at(label, lambda p: offset + p) for offset in range(0, fed.n, 12)
                    for label in cell_keys}
        if Hf.shape != (len(fed_labels), len(fed_labels)) or not expected <= keys:
            note = f'level {level}: the generator and native basis must cover every cell'
            status = 'fail'
            break
        asym, far = generator_defects(Hf, fed_labels, fed_seams)
        if asym > 1e-12 or far != 0:
            note = f'level {level}: Hermitian defect {asym:.1e}, nonlocal couplings {far}'
            status = 'fail'
            break
        built.append(fixed.NAMES[level])
    return {'status': status, 'detail': 'builds on ' + ', '.join(built) + (f' ({note})' if note else '')}


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
            'detail': 'declared parameter descriptions pass provenance and forbidden-input checks' if not findings else '; '.join(findings)}


def evidence_ok(cand, paths):
    base = ROOT / 'model' / 'candidate_laws' / cand.NAME
    missing = [p for p in paths if not (base / p).resolve().exists()]
    if missing:
        fail(f'{cand.NAME}: missing evidence {missing}')


def premise_register():
    """The premise rows of the observables register, in register order."""
    return [r for r in json.loads(REGISTER.read_text(encoding='utf-8'))['rows'] if r.get('kind') == 'premise']


def premise_rows(cand, premises, roles):
    ids = {r['id'] for r in premises}
    unknown = sorted(set(cand.PREMISES) - ids)
    if unknown:
        fail(f'{cand.NAME}: premises not in the register: {unknown}')
    rows = []
    for r in premises:
        kind = roles['types'].get(r['id']) or r['premise_kind']
        role = roles['overrides'].get(r['id'], {}).get('role') or roles['defaults'][kind]
        decl = cand.PREMISES.get(r['id'], {})
        status = decl.get('status', 'open')
        if status not in STATUSES:
            fail(f'{cand.NAME} {r["id"]}: unknown status {status!r}')
        if status == 'fails' and not decl.get('counterexample'):
            fail(f'{cand.NAME} {r["id"]}: a failed premise names its counterexample')
        evidence_ok(cand, decl.get('evidence', []))
        rows.append({'id': r['id'], 'name': r['name'], 'statement': r['statement'],
                     'type': kind, 'role': role, 'status': status,
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


def evaluate(name, shared=None):
    from model.executable_model import laws
    cand = laws.load(name)
    roles = json.loads(ROLES.read_text(encoding='utf-8'))
    checks = carrier_checks(cand)
    checks['parameters'] = parameter_check(cand)
    premises = premise_rows(cand, premise_register(), roles)
    summary = {role: {s: sum(1 for r in premises if r['role'] == role and r['status'] == s) for s in STATUSES}
               for role in ('law', 'designer', 'import')}
    return {'schema': 'oph.law_evaluation.v2', 'candidate': cand.NAME, 'title': cand.TITLE, 'role': cand.ROLE,
            'summary': cand.SUMMARY, 'parameters': cand.PARAMETERS, 'checks': checks, 'joint': joint_rows(cand),
            'premise_counts': summary, 'premises': premises, 'measurement': assessment.readout(name),
            'requirements': assessment.for_candidate(cand, {r['id'] for r in premises}, laws.names(), shared)}


def dumps(value):
    return json.dumps(value, indent=1, ensure_ascii=False) + '\n'


def outputs(results):
    files = {ROOT / 'model' / 'candidate_laws' / r['candidate'] / 'EVALUATION.json': dumps(r) for r in results}
    files[SCORECARD] = report.scorecard(results)
    files[HERE / 'EVIDENCE.md'] = report.evidence_ledger(results)
    text = INDEX.read_text(encoding='utf-8')
    if BEGIN not in text or END not in text:
        fail('CANDIDATE.md has no generated table markers')
    head, rest = text.split(BEGIN, 1)
    files[INDEX] = head + BEGIN + '\n' + report.candidate_table(results, INDEX.parent) + '\n' + END + rest.split(END, 1)[1]
    return files


def criterion(result, identifier):
    key = identifier.casefold()
    for group in result['requirements']:
        for row in group['criteria']:
            if row['id'].casefold() == key:
                return row
    for row in result['premises']:
        if row['id'].casefold() == key:
            return {**row, 'title': row['name'], 'requires': row['statement']}
    if key.startswith('joint-') and key[6:] in result['joint']:
        row = result['joint'][key[6:]]
        return {**row, 'id': key.upper(), 'title': 'Agreement across tests', 'requires': row['requirement']}
    return None


def main(argv=None):
    from model.executable_model import laws
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--law', choices=laws.names())
    ap.add_argument('--criterion', metavar='ID', help='show an individual requirement and candidate evidence')
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = ap.parse_args(argv)
    if (args.write or args.check) and (args.law or args.criterion):
        ap.error('--write and --check evaluate every candidate')
    order = sorted(laws.names(), key=lambda n: (laws.load(n).ROLE != 'main', n))
    shared = assessment.definitions()
    results = [evaluate(n, shared) for n in ([args.law] if args.law else order)]
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
    if args.criterion:
        for r in results:
            row = criterion(r, args.criterion)
            if row is None:
                ap.error(f'unknown criterion {args.criterion!r}')
            if r is results[0]:
                print(row['id'], row['title'], row['requires'], sep='\n')
            note = row['note'] or 'No candidate-specific discharge is registered.'
            print(f"{r['title']} · {report.GRADES[row['status']]} · {note}")
            if 'demonstrated_measurement' in row:
                measurement = row['demonstrated_measurement']
                print(f"  Demonstrated · {measurement['label']} ({measurement['source']})")
            for source in row['evidence']:
                if row.get('role') or args.criterion.casefold().startswith('joint-'):
                    source = (INDEX.parent / r['candidate'] / source).resolve().relative_to(ROOT)
                print('  ', source)
            for study in assessment.comparison_studies(row):
                print(f"  Comparison · {study['law']} · {study['finding']} ({study['source']})")
        return 0
    print(report.scorecard(results), end='')
    return 0


if __name__ == '__main__':
    sys.exit(main())
