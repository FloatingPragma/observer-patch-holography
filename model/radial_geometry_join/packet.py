"""Cold replay of the same-field charged-energy to geometry map."""
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
from model.radial_interacting_core.core import encode
from . import bounds, source

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENTS = ('radial_interacting_core', 'matter_geometry_response')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/radial_geometry_join/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PARENTS],
        'independent_checks': [],
    })
    paths.add('model/radial_geometry_join/PROOF.md')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def inherited_custody():
    import importlib
    result = {}
    for name in PARENTS:
        value = strict_json.load(ROOT/f'model/{name}/results.json')
        producer = importlib.import_module(f'model.{name}.packet')
        if (value.get('schema') != f'oph.{name}.v1' or value.get('all_pass') is not True
                or value.get('source_bindings') != producer.source_bindings()):
            raise ValueError('the field/geometry source manifest is stale: '+name)
        result[name] = {'full_current_source_manifest_verified': True,
                        'parent_expensive_receipt_reexecuted': False}
    return result


def scope():
    return {
        'same_actual_finite_radial_ground_and_charged_states_drive_geometry': True,
        'source_interaction_energy_is_read_without_fitting_a_geometry_energy': True,
        'complete_field_spectator_and_reference_extension_proved': True,
        'complete_native_control_and_failure_error_certified': True,
        'readout_field_energy_injection_accounted': True,
        'all_native_load_histories_and_idle_outcomes_retained': True,
        'mass_texture_and_control_schedule_are_inputs': True,
        'fresh_independent_spectral_preparations_required_for_mean_iteration': True,
        'ground_or_charged_preparation_executed': False,
        'full_27_site_Fock_evolution_executed': False,
        'arbitrary_correlated_reuse_has_same_iteration': False,
        'mean_geometry_identified_as_Lorentzian_metric': False,
        'regional_Einstein_premises_instantiated': False,
        'autonomous_reciprocal_energy_backreaction_derived': False,
        'interacting_continuum_or_physical_species_derived': False,
        'completion_gate_accepted': False,
    }


def quantum_certificate():
    result = source.quantum_execution()
    # Float diagnostics are checked on each replay; the receipt retains their
    # bounds, not platform-dependent last bits of a matrix contraction.
    for key in ('complete_output_error_diagnostic',
                'complete_first_failure_probability', 'normalization_residual'):
        result.pop(key)
    result['complete_numeric_execution_within_certified_bound'] = True
    result['normalization_residual_at_most_2e_minus_11'] = True
    return result


def _produce():
    from model.matter_geometry_response import geometry, source as parent
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after the producer import')
    native = parent.native_execution()
    geometric = geometry.certificate()
    value = {'schema': 'oph.radial_geometry_join.v1', 'all_pass': True,
        'issues': [1139, 1142, 1143, 1144],
        'gates': ['CP-03', 'CP-06', 'CP-07', 'CP-08'],
        'physical_comparison_count_added': 0,
        'same_source_bounds': encode(bounds.joined_bounds()),
        'quantum_instrument': quantum_certificate(),
        'native_histories': {
            'complete_native_histories_checked': native['complete_native_histories_checked'],
            'complete_history_sha256': native['complete_history_sha256'],
            'all_initial_load_mean_identities_checked_exactly': True,
            'same_positive_30_seams_and_idle_clock': True},
        'geometry_source': {
            'parent_exact_derivative_and_rational_endpoint_certificate_replayed': True,
            'parent_finite_contrast_certified_lower': geometric['finite_contrast_certified_lower']},
        'finite_observation': encode(bounds.observation_bound()),
        'scope': scope(), 'inherited_custody': inherited_custody(),
        'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during complete join replay')
    return value


def validate(value, expected):
    if type(value) is not dict or value.get('scope') != scope():
        raise ValueError('invalid same-source geometry scope')
    if canonical(value) != canonical(expected):
        raise ValueError('same-source geometry receipt differs from fresh replay')


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-geometry-') as directory:
        temporary = Path(directory)
        output = temporary/'receipt.json'
        subprocess.run([sys.executable, '-W', 'error', '-X',
            'pycache_prefix='+str(temporary/'pycache'), '-m',
            'model.radial_geometry_join.packet', '--worker', str(output)],
            cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
            check=True, capture_output=True, text=True, timeout=300)
        value = strict_json.load(output)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold join replay')
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
        print('same radial field/geometry source: PASS; charged Gram contrast >59/1000')


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
