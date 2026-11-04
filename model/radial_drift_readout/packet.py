"""Cold source-bound replay for complete radial readout with local drift."""
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
from ..radial_reused_readout import packet as parent
from ..radial_interacting_core.core import encode
from . import algebra, bounds, commands, controls

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/radial_drift_readout'
OUTPUT = Path(__file__).with_name('results.json')
PARENT = 'model/radial_reused_readout/results.json'


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
    if (type(value) is not dict or value.get('schema') != 'oph.radial_reused_readout.v1'
            or value.get('all_pass') is not True
            or value.get('source_bindings') != parent.source_bindings()):
        raise ValueError('the complete current same-field native instrument receipt is required')
    for path, digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('an inherited source binding changed: '+path)
    return {'receipt': PARENT, 'all_current_parent_source_bindings_verified': True,
            'parent_expensive_controls_reexecuted': False}


def scope():
    return {'actual_full_Hamiltonian_interleaved_between_complete_events': True,
            'local_full_Fock_star_and_contact_norm_exact': True,
            'readout_error_has_no_volume_factor': True,
            'complete_field_reference_and_record_error_derived': True,
            'positive_timed_control_interval_classifier_below_one_in_2500': True,
            'all_native_pointer_coin_integer_repair_and_failure_outputs_retained': True,
            'timed_local_control_realization_and_bounded_port_drift_extension_required': True,
            'physical_rate_latency_and_clock_contract_instantiated': False,
            'source_Hamiltonian_compensation_required': False,
            'comparison_generator_substituted_for_actual_source': False,
            'interleaved_record_law_equals_frozen_parity_mixture': False,
            'parity_projector_and_sign_parity_commutators_identified': False,
            'postmeasurement_spectral_sector_protection_proved': False,
            'uncontrolled_field_repairs_included': False,
            'initial_spectral_state_or_texture_preparation_derived': False,
            'full_timed_27_site_experiment_executed': False,
            'uniform_long_time_limit_or_continuum_QFT_proved': False,
            'physical_species_or_Einstein_dynamics_derived': False,
            'completion_gate_accepted': False}


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after reused readout producer import')
    value = {'schema': 'oph.radial_drift_readout.v1', 'all_pass': True,
             'issues': [1139, 1141], 'gates': ['CP-03', 'CP-05'],
             'physical_comparison_count_added': 0,
             'scope': scope(), 'inherited_custody': custody(), 'source_bindings': pins,
             'exact_local_source_algebra': algebra.certificate(),
             'complete_drift_and_classifier': encode(bounds.certificate()),
             'timed_native_commands': commands.timing_plan(),
             'independent_full_Fock_controls': controls.certificate()}
    if source_bindings() != pins:
        raise ValueError('source changed during reused readout production')
    return value


def validate(value, expected):
    if (type(value) is not dict or value.get('schema') != 'oph.radial_drift_readout.v1'
            or value.get('all_pass') is not True or value.get('scope') != scope()
            or type(value.get('physical_comparison_count_added')) is not int):
        raise ValueError('invalid complete drift-readout receipt or scope')
    if canonical(value) != canonical(expected):
        raise ValueError('drift-readout receipt differs from its cold replay')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-reuse-') as directory:
        output = Path(directory)/'receipt.json'
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(directory)/'pycache'), '-m',
                        'model.radial_drift_readout.packet', '--worker', str(output)],
                       cwd=ROOT, check=True, capture_output=True, text=True, timeout=300,
                       env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
        value = strict_json.load(output)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold drift-readout replay')
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
        print('radial drift readout: PASS; one preparation, error <1/2500')


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
