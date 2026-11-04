"""Cold replay of the complete radial local-current classification."""
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
from ..common_source import native_matter_attachment
from ..radial_interacting_core import packet as core_packet
from . import algebra, controls, representation

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/radial_current_map'
OUTPUT = Path(__file__).with_name('results.json')
PARENTS = ('model/radial_interacting_core/results.json',
           'model/common_source/native_matter_attachment_results.json')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': list(PARENTS), 'independent_checks': [PACKAGE+'/PROOF.md']})
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def custody():
    for path, producer, schema in zip(PARENTS, (core_packet, native_matter_attachment),
            ('oph.radial_interacting_core.v1', 'oph.native_matter_attachment.v1'), strict=True):
        value = strict_json.load(ROOT/path)
        if value.get('schema') != schema or value.get('source_bindings') != producer.source_bindings():
            raise ValueError('the current complete source parent is required: '+path)
        if path == PARENTS[0] and value.get('all_pass') is not True:
            raise ValueError('the complete interacting isolation receipt must pass')
        for source, digest in value['source_bindings'].items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != digest:
                raise ValueError('an inherited source changed: '+source)
    return {'parents': list(PARENTS), 'all_current_source_bindings_verified': True,
            'expensive_native_detector_execution_repeated': False}


def scope():
    return {'all_onsite_even_Hermitian_CAR_currents_classified': True,
            'quartic_parity_and_pairing_directions_retained': True,
            'same_interacting_source_eta_current_intertwiner_exact': True,
            'actual_charged_pair_is_eta_doublet_and_ground_pair_eta_singlets': True,
            'center_spin_leakage_outside_complete_low_quartet_certified': True,
            'faithful_full_native_SM_current_in_classified_local_symmetry_algebra': False,
            'auxiliary_matter_register_added_to_radial_source': False,
            'eta_identified_with_electroweak_or_Lorentz_spin': False,
            'all_nonlocal_or_broken_gauge_currents_classified': False,
            'native_repair_preserves_eta_proved': False,
            'physical_defect_to_particle_map_complete': False,
            'source_texture_vacuum_or_couplings_selected': False,
            'interacting_continuum_or_Einstein_limit_proved': False,
            'completion_gate_accepted': False}


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after current-map producer import')
    value = {'schema': 'oph.radial_current_map.v1', 'all_pass': True,
             'issues': [1140, 1141, 1142, 1143], 'gates': ['CP-04', 'CP-05', 'CP-06', 'CP-07'],
             'physical_comparison_count_added': 0, 'scope': scope(),
             'inherited_custody': custody(), 'source_bindings': pins,
             'exact_source_current_classification': algebra.certificate(),
             'actual_low_current_map_and_leakage': representation.certificate(),
             'independent_controls': controls.certificate()}
    if source_bindings() != pins:
        raise ValueError('source changed during current-map production')
    return value


def validate(value, expected):
    if (type(value) is not dict or value.get('schema') != 'oph.radial_current_map.v1'
            or value.get('all_pass') is not True or value.get('scope') != scope()
            or type(value.get('physical_comparison_count_added')) is not int
            or value.get('physical_comparison_count_added') != 0):
        raise ValueError('invalid current-map scope')
    if canonical(value) != canonical(expected):
        raise ValueError('current-map receipt differs from cold replay')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-current-') as directory:
        output = Path(directory)/'receipt.json'
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(directory)/'pycache'), '-m',
                        'model.radial_current_map.packet', '--worker', str(output)],
                       cwd=ROOT, check=True, capture_output=True, text=True, timeout=600,
                       env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
        value = strict_json.load(output)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold current-map replay')
    validate(value, value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = _produce() if args.worker else report()
    if args.check: validate(strict_json.load(args.check), value)
    else: (args.worker or args.write).write_text(canonical(value), encoding='utf-8', newline='\n')
    if not args.worker: print('radial current map: PASS; exact eta multiplets and complete local classification')


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
