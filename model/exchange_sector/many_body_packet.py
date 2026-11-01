"""Fresh replay of the finite occupation catalogue and complete source readout."""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'many_body_results.json'
CAMPAIGN = 'native-many-particle-statistics-20261031-v1'
CONTRACT_SHA256 = 'ca959087e9578610ec993623926f3e526cdec984ecc236923a94416d6266c49c'
RESULT_KEYS = frozenset(('schema', 'campaign_id', 'contract_sha256', 'scope', 'domain',
    'deformations', 'topology', 'compatibility', 'complete_instruments', 'source_observer',
    'physical_particle_identification', 'physical_comparison_count_added',
    'localized_protected_defect_constructed', 'nonabelian_many_particle_transport_classified',
    'refinement_or_continuum_equivalence_proved', 'integrated_completion_gate_accepted', 'source_sha256'))


def source_bindings():
    """Static repository import closure plus the source's archived data roots."""
    from .statistics_packet import source_bindings as inherited
    names = set(inherited())
    names.update(p.relative_to(ROOT).as_posix() for pattern in
                 ('many_body*.py', 'test_many_body*.py') for p in HERE.glob(pattern))
    names.update(('model/exchange_sector/MANY_BODY.md',
                  'model/exchange_sector/many_body_contract.json',
                  'model/requirements.txt', 'requirements.txt'))
    pending, visited = list(names), set()

    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        path = module.replace('.', '/')
        parts = path.split('/')
        candidates = [path+'.py', path+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        path = ROOT/name
        module = name.removesuffix('.py').replace('/', '.')
        package = (module.removesuffix('.__init__') if module.endswith('.__init__')
                   else module.rpartition('.')[0])
        add(module)
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _report():
    """Owned child entry point; public replay always starts a fresh interpreter."""
    from . import many_body, many_body_compatibility, many_body_service, many_body_topology
    before = source_bindings()
    if before['model/exchange_sector/many_body_contract.json'] != CONTRACT_SHA256:
        raise ValueError('many-particle contract changed')
    result = {'schema': 'oph.native_many_body_statistics.receipt.v1',
              'campaign_id': CAMPAIGN, 'contract_sha256': CONTRACT_SHA256,
              'scope': 'FINITE_COMPONENT_FAMILY',
              'domain': 'all 4096 binary occupations of the fixed native twelve-port carrier',
              'deformations': 'one-particle native face triangles and disjoint-move squares',
              'topology': many_body_topology.certificate(),
              'compatibility': many_body_compatibility.certificate(),
              'complete_instruments': many_body.certificate(),
              'source_observer': many_body_service.certificate(),
              'physical_particle_identification': False,
              'physical_comparison_count_added': 0,
              'localized_protected_defect_constructed': False,
              'nonabelian_many_particle_transport_classified': False,
              'refinement_or_continuum_equivalence_proved': False,
              'integrated_completion_gate_accepted': False,
              'source_sha256': before}
    if source_bindings() != before:
        raise RuntimeError('many-particle sources changed during replay')
    return result


def report():
    before = source_bindings()
    command = ('import sys,json; sys.path.insert(0,'+repr(str(ROOT))+'); '
               'from model.exchange_sector.many_body_packet import _report; '
               'print(json.dumps(_report(),sort_keys=True,allow_nan=False))')
    environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    # -B disables writes but still permits reading stale timestamp bytecode.
    # An empty private cache prefix prevents both warm modules and stale pyc
    # from being rebound to a newly computed source manifest.
    with tempfile.TemporaryDirectory(prefix='oph-many-particle-') as cache:
        run = subprocess.run([sys.executable, '-I', '-B', '-W', 'error', '-X', 'utf8',
                              '-X', 'pycache_prefix='+cache, '-c', command],
                             cwd=ROOT, env=environment, check=True, capture_output=True,
                             text=True, encoding='utf-8', timeout=600)
    result = strict_json.loads(run.stdout)
    if (type(result) is not dict or set(result) != RESULT_KEYS
            or result.get('source_sha256') != before
            or source_bindings() != before):
        raise RuntimeError('fresh replay and current source custody disagree')
    if (result['schema'] != 'oph.native_many_body_statistics.receipt.v1'
            or result['campaign_id'] != CAMPAIGN or result['contract_sha256'] != CONTRACT_SHA256
            or result['scope'] != 'FINITE_COMPONENT_FAMILY'
            or any(result[key] is not False for key in RESULT_KEYS if key in (
                'physical_particle_identification', 'localized_protected_defect_constructed',
                'nonabelian_many_particle_transport_classified', 'refinement_or_continuum_equivalence_proved',
                'integrated_completion_gate_accepted'))
            or type(result['physical_comparison_count_added']) is not int
            or result['physical_comparison_count_added'] != 0):
        raise ValueError('fresh replay changed the finite research contract')
    # Exact JSON comparison rejects NaN, type changes and promoted scope too.
    compare(result, result)
    return result


def compare(actual, expected):
    def encode(value):
        return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
    if encode(actual) != encode(expected):
        raise ValueError('many-particle receipt differs from complete fresh replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    mode.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = report()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n',
                              encoding='utf-8', newline='\n')
    else:
        compare(strict_json.load(args.check), result)
    print('Verified all occupation sectors, compatible exchange lifts and complete native-source readouts.')


if __name__ == '__main__':
    main()
