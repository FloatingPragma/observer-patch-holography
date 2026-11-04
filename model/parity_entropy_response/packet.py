"""Cold source-bound replay of regional parity entropy and modular response."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..matter_geometry_response import packet as geometry_packet
from ..radial_interacting_core import packet as core_packet
from . import modular, records, regional

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/parity_entropy_response'
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = {'matter_geometry_response': geometry_packet, 'radial_interacting_core': core_packet}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PRIORS],
        'independent_checks': [PACKAGE+'/PROOF.md'],
    })
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def custody():
    result = {}
    for name, producer in PRIORS.items():
        path = f'model/{name}/results.json'
        value = strict_json.load(ROOT/path)
        if (type(value) is not dict or value.get('schema') != 'oph.'+name+'.v1'
                or value.get('all_pass') is not True
                or value.get('source_bindings') != producer.source_bindings()):
            raise ValueError('complete current parent source receipt required: '+name)
        for source, digest in value['source_bindings'].items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != digest:
                raise ValueError('inherited source binding changed: '+source)
        result[name] = {'receipt': path, 'complete_parent_source_manifest_verified': True,
                        'every_parent_source_digest_verified': True,
                        'expensive_parent_controls_reexecuted': False}
    return result


def scope():
    return {'same_actual_bounded_radial_interacting_Hamiltonian': True,
            'actual_ground_central_M4_faithfulness_derived': True,
            'actual_ground_full_Fock_faithfulness': False,
            'stationary_faithful_affine_source_family_constructed': True,
            'same_family_regional_entropy_full_energy_and_native_bias_related': True,
            'complete_site_modular_word_with_native_failures_executed': True,
            'regional_contact_clock_and_full_Gibbs_clock_distinguished': True,
            'entropy_maximum_at_fixed_native_kinetic_bias_proved': True,
            'complete_committed_record_entropy_identity_proved': True,
            'same_H_Gibbs_measurement_work_balance_proved': True,
            'spectral_state_preparation_or_actual_ground_q_calibration_executed': False,
            'physical_precision_drift_and_fresh_probe_service_are_inputs': True,
            'parity_copy_instrument_implemented_in_this_package': False,
            'regional_entropy_identity_substituted_for_global_entropy': False,
            'signed_kinetic_bias_identified_with_area': False,
            'full_Hamiltonian_flow_preserves_complete_site_algebra': False,
            'finite_G3_or_G4_Einstein_premise_discharged': False,
            'Lorentzian_stress_or_fixed_volume_area_law_derived': False,
            'continuum_interacting_QFT_or_Einstein_dynamics_constructed': False,
            'physical_species_identified': False,
            'completion_gate_accepted': False}


def _same(actual, expected, path='receipt'):
    if type(actual) is not type(expected):
        raise ValueError('receipt type changed: '+path)
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError('receipt fields changed: '+path)
        for key in expected:
            _same(actual[key], expected[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError('receipt array length changed: '+path)
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            _same(left, right, path+'.'+str(index))
    elif isinstance(expected, float):
        # Floats are independently gated matrix-execution diagnostics. All
        # theorem coefficients, source bounds and scope flags are exact.
        if not math.isfinite(actual) or abs(actual-expected) > 2e-11*max(1., abs(expected)):
            raise ValueError('numeric control differs from cold replay: '+path)
    elif actual != expected:
        raise ValueError('exact receipt value changed: '+path)


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'regional', 'modular',
                'records', 'scope', 'inherited_custody', 'source_bindings',
                'physical_comparison_count_added'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.parity_entropy_response.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid regional entropy/modular source receipt')
    for key, expected in [('scope', scope()), ('regional', regional.certificate()),
                          ('inherited_custody', custody()), ('source_bindings', source_bindings())]:
        _same(value[key], expected, key)
    if (value['issues'] != [1144, 1145] or value['gates'] != ['CP-08', 'CP-09']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('regional entropy cannot accept a physical comparison')
    if replay is not None:
        _same(value, replay)


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after entropy producer import')
    value = {'schema': 'oph.parity_entropy_response.v1', 'all_pass': True,
             'issues': [1144, 1145], 'gates': ['CP-08', 'CP-09'],
             'regional': regional.certificate(), 'modular': modular.certificate(),
             'records': records.certificate(), 'scope': scope(),
             'inherited_custody': custody(), 'source_bindings': pins,
             'physical_comparison_count_added': 0}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during regional entropy replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-parity-entropy-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(folder)/'pycache'),
                        '-m', 'model.parity_entropy_response.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True, timeout=300)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold regional entropy replay')
    validate(value)
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
        print(json.dumps({'all_pass': True, 'source_bindings': len(value['source_bindings']),
                          'actual_ground_local_faithfulness': True,
                          'same_source_regional_entropy_response': True,
                          'Einstein_premises_discharged': False}))


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
