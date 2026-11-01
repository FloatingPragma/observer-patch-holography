"""Cold, complete source replay of the bounded local exchange-charge result."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..exchange_sector import statistics_packet

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENTS = ('model/defect_field_bridge/full_domain_results.json',)
CAMPAIGN = 'OPH-LOCAL-EXCHANGE-CHARGE-2026-11-01-v1'


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/local_exchange_charge/README.md',
        'source_receipts': list(PARENTS), 'independent_checks': [],
    })
    # The archive geometry is imported through its own search root. Bind its
    # current static closure as well. The ordered-cover result is recomputed
    # on every column in core.py; no saved protection receipt is consumed.
    paths.update(statistics_packet.source_bindings())
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def _produce():
    # Imports belong inside the cold worker. The public report always runs a
    # fresh interpreter with a private cache; parent module/cache mutations
    # cannot manufacture a successful source replay.
    from . import core, optimality, service
    pins = source_bindings()
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        inherited = parent.get('source_bindings', parent.get('source_sha256'))
        if not inherited or any(pins.get(path) != digest for path, digest in inherited.items()):
            raise ValueError('predecessor source custody changed: '+name)
    value = {
        'schema': 'oph.local_exchange_charge.v1', 'campaign_id': CAMPAIGN,
        'scope': 'FINITE_COMPONENT_FAMILY', 'whole_program_frozen': False,
        'physical_comparison_count_added': 0,
        'local_instrument_and_transport': core.certificate(),
        'locality_optimality': optimality.certificate(),
        'source_execution': service.certificate(),
        'remaining_obligations': {
            'physical_particle_or_electron_identified': False,
            'sharp_local_charge_constructed': False,
            'bare_binary_law_supplies_coherent_labels': False,
            'nontrivial_modular_clock_or_mass_derived': False,
            'protection_against_label_sensitive_services': False,
            'spatial_locality_of_processor_hardware_derived': False,
            'refinement_limit_proved': False,
            'canonical_physical_premises_discharged': [],
        },
        'source_bindings': pins,
    }
    if source_bindings() != pins:
        raise ValueError('source changed during local charge replay')
    # Normalize tuples to their public JSON representation before returning.
    return strict_json.loads(canonical(value))


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-local-charge-') as directory:
        folder = Path(directory)
        path = folder/'result.json'
        env = dict(os.environ)
        env['OPENBLAS_NUM_THREADS'] = env['OMP_NUM_THREADS'] = '1'
        env['PYTHONUTF8'] = '1'
        env['PYTHONPYCACHEPREFIX'] = str(folder/'pycache')
        completed = subprocess.run([
            sys.executable, '-W', 'error', '-X', 'pycache_prefix='+env['PYTHONPYCACHEPREFIX'],
            '-m', 'model.local_exchange_charge.packet', '--worker', str(path),
        ], cwd=ROOT, env=env, capture_output=True, text=True)
        if completed.returncode:
            raise RuntimeError('cold local-charge replay failed:\n'+completed.stderr)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold local charge replay')
    return value


def compare(retained, fresh):
    # Exact JSON types reject bool/int substitutions and invented resource,
    # output, scope and custody values. Numerical admission happens in the
    # fresh scientific producer before this comparison.
    if canonical(retained) != canonical(fresh):
        raise ValueError('local exchange-charge receipt differs from complete source replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    modes.add_argument('--worker', type=Path)
    args = parser.parse_args()
    if args.worker:
        args.worker.write_text(canonical(_produce()), encoding='utf-8', newline='\n')
        return
    value = report()
    if args.write:
        args.write.write_text(canonical(value), encoding='utf-8', newline='\n')
    else:
        compare(strict_json.load(args.check), value)
    print('Verified local exchange-charge realization, optimality and complete source history.')


if __name__ == '__main__':
    main()
