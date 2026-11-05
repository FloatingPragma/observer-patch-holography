"""Cold exact replay of the composite spatial and source-rate certificate."""
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
from ..charged_composite_search import packet as inherited
from ..centered_particle_spectrum.packet import inherited_inputs
from . import certificate, source

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/composite_spatial_response/PROOF.md',
        'source_receipts': ['model/charged_composite_search/results.json',
                            'model/centered_particle_spectrum/results.json'],
        'independent_checks': []})
    paths.update({'model/composite_spatial_response/README.md',
                  'model/charged_composite_search/witness.npz',
                  'model/centered_particle_spectrum/witness.npz'})
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def produce():
    pins = source_bindings()
    block = source.previous.block(tuple(map(int, inherited.load_witness()['columns'])))
    result = {'schema': 'oph.composite_spatial_response.v1', 'all_pass': True,
              'issues': [1141, 1143], 'gates': ['CP-05', 'CP-07'],
              'interface': {
                  'state': 'unique centered N8 hidden-spin-zero eta-zero ground of the eight-site cube',
                  'algebra': 'complete sixteen-mode CAR source; cubic probes on all eight source sites',
                  'operation': 'signed spatial transports, exact spin/eta ladders and original-source kinetic/contact commutators',
                  'clock': 'K=H-N/2, c=1/8 U=1 and hbar_model=1; no SI scale',
                  'refinement': 'fixed open 2x2x2 cube; no continuum or infinite-volume identification'},
              'source_covariance': source.covariance_certificate(),
              'source_spatial_commutant': source.spatial_commutant_certificate(),
              'source_kinetic_anchor': source.kinetic_source_anchor(),
              'source_rate_algebra': source.commutator_certificate(),
              'spin_eta_ladders': source.multiplet_certificate(block),
              'spectrum': certificate.verify(),
              'inherited_ground_inputs': inherited_inputs(),
              'canonical_physical_premises_discharged': [],
              'integrated_gate_accepted': False, 'source_law_changed': False,
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('composite spatial sources changed during exact replay')
    return result


def validate(value, expected):
    if canonical(value) != canonical(expected):
        raise ValueError('composite spatial receipt differs from exact reconstruction')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-composite-spatial-') as temp:
        destination = Path(temp)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        process = subprocess.run([sys.executable, '-W', 'error', '-X',
                                  'pycache_prefix='+str(Path(temp)/'pycache'),
                                  '-m', 'model.composite_spatial_response.packet',
                                  '--worker', str(destination)], cwd=ROOT, env=env,
                                 capture_output=True, text=True)
        if process.returncode:
            raise RuntimeError('cold composite spatial replay failed:\n'+process.stderr)
        result = strict_json.load(destination)
    if pins != source_bindings() or result['source_bindings'] != pins:
        raise ValueError('composite spatial sources changed across cold replay')
    return result


def verify(path=OUTPUT):
    validate(strict_json.load(Path(path)), report())
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    result = produce() if args.worker else report()
    if args.check:
        validate(strict_json.load(args.check), result)
    else:
        (args.worker or args.write).write_text(canonical(result), encoding='utf-8')
    if not args.worker:
        print('Composite spatial residues and source-rate certificate OK')


if __name__ == '__main__':
    main()
