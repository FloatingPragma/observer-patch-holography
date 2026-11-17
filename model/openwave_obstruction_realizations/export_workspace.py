#!/usr/bin/env python3
"""Export the declared evidence closure, preserving its workspace-relative paths."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
BASE = WORKSPACE / 'openwave-tuning/next/OBSTRUCTION_REPAIR'
SHA = re.compile(r'^[0-9a-f]{64}$')


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    ledger = json.loads((BASE / 'ISSUES.json').read_text())
    selected, pending, edges, missing, changed, leaves = {}, [], [], [], [], []
    resolved = {}
    basename_index = {}
    for directory, dirs, names in os.walk(BASE):
        dirs[:] = [x for x in dirs if x not in {'.git', '.venv', 'venv', '__pycache__'}]
        for name in names:
            basename_index.setdefault(name, []).append(Path(directory)/name)
    digest_cache = {}

    def cached_digest(p):
        if p not in digest_cache:
            digest_cache[p] = digest(p)
        return digest_cache[p]

    def resolve(raw, origin):
        if not isinstance(raw, str) or '\n' in raw:
            return None
        cache_key = (str(origin.parent), raw)
        if cache_key in resolved:
            return resolved[cache_key]
        candidates = [origin.parent / raw, BASE / raw, WORKSPACE / raw,
                      WORKSPACE / 'reverse-engineering-reality' / raw]
        for anchor in ('openwave-tuning/', 'openwave-push/', 'reverse-engineering-reality/'):
            if anchor in raw:
                candidates.append(WORKSPACE / raw[raw.index(anchor):])
        for p in candidates:
            p = Path(os.path.abspath(p))
            if p.is_relative_to(WORKSPACE) and p.is_file():
                resolved[cache_key] = p
                return p
        resolved[cache_key] = None
        return None

    def add(p, origin=None, expected=None):
        key = str(p.relative_to(WORKSPACE))
        actual = selected[key]['sha256'] if key in selected else cached_digest(p)
        if expected and actual != expected:
            changed.append({'origin': str(origin.relative_to(WORKSPACE)),
                            'path': key, 'expected': expected, 'actual': actual})
        if origin:
            edges.append({'origin': str(origin.relative_to(WORKSPACE)),
                          'path': key, 'sha256': expected})
        if key not in selected:
            selected[key] = {'sha256': actual, 'bytes': p.stat().st_size}
            pending.append(p)

    for name in ('ISSUES.json', 'README.md', 'SOURCE_JOIN_REGISTER.md',
                 'CONTINUUM_GEOMETRY.md', 'LATEST_MANIFEST_AUDIT.json'):
        add(BASE / name)
    for issue in ledger['issues']:
        for raw in issue.get('evidence', []):
            p = resolve(raw, BASE / 'ISSUES.json')
            if p:
                add(p)
            else:
                missing.append({'origin': 'ISSUES.json', 'path': raw})
    for subdir in ('gravity', 'species_leptons', 'hadrons'):
        for p in (BASE / subdir).rglob('*MANIFEST*'):
            if (p.suffix in ('.json', '.sha256') and p.is_file()
                    and not {'inputs', 'provenance', 'consumed'}.intersection(p.relative_to(BASE).parts)):
                add(p)

    def pin(raw, expected, origin):
        p = resolve(raw, origin)
        if p is None or cached_digest(p) != expected:
            matches = [candidate for candidate in basename_index.get(Path(raw).name, [])
                       if cached_digest(candidate) == expected]
            if matches:
                p = sorted(matches, key=lambda candidate: (len(candidate.parts), str(candidate)))[0]
        if p:
            add(p, origin, expected)
        else:
            missing.append({'origin': str(origin.relative_to(WORKSPACE)),
                            'path': raw, 'sha256': expected})

    def walk(value, origin):
        if isinstance(value, dict):
            if isinstance(value.get('path'), str) and SHA.fullmatch(str(value.get('sha256', ''))):
                pin(value['path'], value['sha256'], origin)
            for key, child in value.items():
                if isinstance(child, str) and SHA.fullmatch(child) and ('.' in key or '/' in key):
                    pin(key, child, origin)
                elif isinstance(child, (dict, list)):
                    walk(child, origin)
        elif isinstance(value, list):
            for child in value:
                walk(child, origin)

    while pending:
        p = pending.pop()
        if len(selected) % 100 == 0:
            print(f'Closure: {len(selected)} files; reading {p.relative_to(WORKSPACE)}', flush=True)
        if p.is_relative_to(BASE):
            leaf = bool({'inputs', 'provenance', 'consumed'}.intersection(p.relative_to(BASE).parts))
        else:
            leaf = True
        external_own_manifest = (not p.is_relative_to(BASE)
            and p.name == 'MANIFEST.json'
            and ('openwave_row1_integrity' in p.parts or 'reviews' in p.parts)
            and not {'inputs', 'provenance', 'consumed'}.intersection(p.parts))
        if leaf and not external_own_manifest:
            if p.suffix == '.json':
                leaves.append(str(p.relative_to(WORKSPACE)))
            continue
        if p.suffix == '.json':
            data = json.loads(p.read_text())
            if external_own_manifest:
                data = {k: data[k] for k in ('files', 'files_sha256') if k in data}
            walk(data, p)
        elif p.suffix == '.sha256':
            for line in p.read_text().splitlines():
                fields = line.split(maxsplit=1)
                if len(fields) == 2 and SHA.fullmatch(fields[0]):
                    pin(fields[1].lstrip('*'), fields[0], p)

    manifests = {k: v['sha256'] for k, v in selected.items()
                 if 'MANIFEST' in Path(k).name and k.endswith(('.json', '.sha256'))}
    record = {'schema': 'oph.obstruction_realizations.export.v1',
              'source_layout': 'Workspace-relative paths below inputs/.',
              'files': dict(sorted(selected.items())),
              'manifests': dict(sorted(manifests.items())),
              'dependency_edges': edges, 'missing_pins': missing,
              'changed_pins': changed, 'failed_construction_count': len(ledger['failed_constructions']),
              'consumed_leaf_documents': sorted(leaves),
              'dependency_boundary': 'Own construction manifests and declared input pins are checked. Copied input records and external receipts are consumed leaves; their upstream scientific claims are not independently rederived by this export. The native row-one implementation and external review own-file manifests are expanded for portability.',
              'physical_rows_solved': False, 'official_rows_changed': False}
    report = {k: record[k] for k in ('missing_pins', 'changed_pins', 'failed_construction_count')}
    report.update(files=len(selected), manifests=len(manifests),
                  bytes=sum(v['bytes'] for v in selected.values()))
    (HERE / 'EXPORT_DRY_RUN.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if args.write:
        if missing or changed:
            raise SystemExit('Refusing export with unresolved or changed pins.')
        for key in selected:
            target = HERE / 'inputs' / key
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(WORKSPACE / key, target)
        (HERE / 'EXPORT_MANIFEST.json').write_text(json.dumps(record, indent=2) + '\n')
        keys = ('id', 'rows', 'problem', 'observed_failure', 'possible_solution',
                'completion_test', 'status', 'evidence', 'latest_result')
        register = {
            'schema': 'oph.obstruction_realizations.sources.v1',
            'source_ledger_sha256': record['files'][str((BASE/'ISSUES.json').relative_to(WORKSPACE))]['sha256'],
            'source_layout': 'Evidence paths are relative to inputs/openwave-tuning/next/OBSTRUCTION_REPAIR/.',
            'components': [{k: issue[k] for k in keys} for issue in ledger['issues']],
            'failed_constructions': ledger['failed_constructions'],
            'source_join_statement': 'Separate constructions require explicit state, operator, reference and clock transfer before joint physical claims.',
            'physical_rows_solved': False, 'official_rows_changed': False,
            'owning_papers': [
                'paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex',
                'paper/deriving_the_particle_zoo_from_observer_consistency.tex',
                'paper/executable_observer_construction.tex']}
        (HERE/'SOURCE_REGISTER.json').write_text(json.dumps(register, indent=2)+'\n')


if __name__ == '__main__':
    main()
