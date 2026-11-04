"""Cold source-bound replay for one-preparation repeated radial readout."""
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
from ..radial_geometry_join import packet as parent
from ..radial_interacting_core.core import encode
from . import controls, reuse

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/radial_reused_readout'
OUTPUT = Path(__file__).with_name('results.json')
PARENT = 'model/radial_geometry_join/results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [PARENT], 'independent_checks': [PACKAGE+'/PROOF.md'],
    })
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def custody():
    value = strict_json.load(ROOT/PARENT)
    if (type(value) is not dict or value.get('schema') != 'oph.radial_geometry_join.v1'
            or value.get('all_pass') is not True
            or value.get('source_bindings') != parent.source_bindings()):
        raise ValueError('the complete current same-field native instrument receipt is required')
    for path, digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('an inherited source binding changed: '+path)
    return {'receipt': PARENT, 'all_current_parent_source_bindings_verified': True,
            'parent_expensive_controls_reexecuted': False}


def scope():
    return {'one_prepared_actual_radial_field_reused': True,
            'complete_conditional_field_reference_coherences_retained': True,
            'complete_native_work_history_mean_derived': True,
            'ground_versus_either_charged_error_below_one_in_2500': True,
            'complete_control_failures_and_halted_suffix_scored': True,
            'full_Fock_measurement_energy_injection_accounted': True,
            'fresh_probe_registers_and_independent_coins_required': True,
            'readout_only_block_without_uncompensated_Hamiltonian_or_repair_required': True,
            'independent_fresh_field_preparation_per_round_required': False,
            'initial_parity_dephasing_inserted': False,
            'reused_mean_replaced_by_power_of_averaged_matrix': False,
            'ensemble_parity_probability_estimated_consistently_from_one_field': False,
            'postmeasurement_field_stays_in_initial_spectral_sector': False,
            'initial_spectral_state_preparation_executed': False,
            'full_32_round_physical_readout_executed': False,
            'charge_sign_resolved': False,
            'autonomous_Hamiltonian_and_repair_detector_derived': False,
            'continuum_interacting_QFT_or_Einstein_dynamics_derived': False,
            'physical_species_identified': False, 'completion_gate_accepted': False}


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after reused readout producer import')
    value = {'schema': 'oph.radial_reused_readout.v1', 'all_pass': True,
             'issues': [1143, 1144], 'gates': ['CP-07', 'CP-08'],
             'physical_comparison_count_added': 0,
             'scope': scope(), 'inherited_custody': custody(), 'source_bindings': pins,
             'history': encode({'rounds': 32, 'ordered_binary_histories': 2**32,
                               'unread': reuse.unread_coefficients(32),
                               'one_ordered_balanced_history': reuse.history_coefficients(32, 16)}),
             'native_reused_mean': encode(reuse.native_transition_control()),
             'single_preparation_classifier': encode(reuse.classification_bound()),
             'resources': reuse.resources(), 'independent_controls': controls.certificate()}
    if source_bindings() != pins:
        raise ValueError('source changed during reused readout production')
    return value


def validate(value, expected):
    if (type(value) is not dict or value.get('schema') != 'oph.radial_reused_readout.v1'
            or value.get('all_pass') is not True or value.get('scope') != scope()
            or type(value.get('physical_comparison_count_added')) is not int):
        raise ValueError('invalid complete reused-readout receipt or scope')
    if canonical(value) != canonical(expected):
        raise ValueError('reused-readout receipt differs from its cold replay')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-reuse-') as directory:
        output = Path(directory)/'receipt.json'
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(directory)/'pycache'), '-m',
                        'model.radial_reused_readout.packet', '--worker', str(output)],
                       cwd=ROOT, check=True, capture_output=True, text=True, timeout=300,
                       env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
        value = strict_json.load(output)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold reused-readout replay')
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
    if args.check:
        validate(strict_json.load(args.check), value)
    else:
        (args.worker or args.write).write_text(canonical(value), encoding='utf-8', newline='\n')
    if not args.worker:
        print('reused radial readout: PASS; one preparation, error <1/2500')


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
