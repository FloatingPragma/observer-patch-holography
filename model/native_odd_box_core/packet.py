"""Cold source-bound replay of the unmodified finite charged core and decoder."""
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
from ..centered_source_ground import packet as centered
from ..native_spin_observer import packet as observer
from . import core, decoder


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/native_odd_box_core'
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = {'centered_source_ground': centered, 'native_spin_observer': observer}


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PRIORS],
        'independent_checks': [],
    })
    paths.update({PACKAGE+'/PROOF.md', 'model/native_spectral_detector/contract.cpp'})
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def custody():
    result = {}
    for name, producer in PRIORS.items():
        row = strict_json.load(ROOT/f'model/{name}/results.json')
        if (row.get('schema') != 'oph.'+name+'.v1'
                or type(row.get('all_pass')) is not bool
                or row.get('source_bindings') != producer.source_bindings()):
            raise ValueError('inherited source manifest is incomplete or stale: '+name)
        for path, digest in row['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('inherited source custody changed: '+path)
        result[name] = {'historical_all_pass': row['all_pass'],
                        'complete_inherited_source_manifest_verified': True,
                        'inherited_expensive_receipt_reexecuted': False}
    if result['centered_source_ground']['historical_all_pass'] is not True:
        raise ValueError('the centered source clock requires its accepted receipt')
    # The observer's broader protected-spin dependency has a false aggregate
    # status. Local native identities, errors and this readout are reexecuted;
    # this theorem does not import protected-spin convergence or repair.
    if result['native_spin_observer']['historical_all_pass'] is not False:
        raise ValueError('the broader observer historical status changed')
    return result


def scope():
    return {'actual_unmodified_centered_source_ground_proved': True,
            'source_contact_and_hopping_retained': True,
            'all_particle_number_sectors_and_virtual_modes_retained': True,
            'source_selected_onsite_mass_or_texture_assumed': False,
            'positive_finite_charged_gap_proved': True,
            'odd_open_boundary_required': True,
            'native_complete_four_outcome_occupation_decoder_constructed': True,
            'chain_full_carrier_native_readout_executed': True,
            'cube_full_Fock_matrix_or_native_readout_executed': False,
            'filled_sea_native_preparation_executed': False,
            'exact_interacting_spectral_projector_readout_proved': False,
            'detector_protected_against_interleaved_repair_or_Hamiltonian_drift': False,
            'volume_uniform_nonzero_charged_gap_proved': False,
            'topological_particle_localization_proved': False,
            'physical_particle_species_identified': False,
            'interacting_relativistic_continuum_constructed': False,
            'Einstein_dynamics_constructed': False,
            'completion_gate_accepted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


DIAGNOSTICS = {'complete_output_error_diagnostic',
               'complete_first_failure_probability', 'normalization_residual'}


def compare(left, right, path='receipt'):
    if type(left) is not type(right):
        raise ValueError(path+': JSON type changed')
    if type(right) is dict:
        if left.keys() != right.keys():
            raise ValueError(path+': keys changed')
        for key in left:
            compare(left[key], right[key], path+'.'+key)
    elif type(right) is list:
        if len(left) != len(right):
            raise ValueError(path+': list length changed')
        for index, (a, b) in enumerate(zip(left, right, strict=True)):
            compare(a, b, path+'.'+str(index))
    elif type(right) is float:
        if (path.rsplit('.', 1)[-1] not in DIAGNOSTICS or not math.isfinite(right)
                or not math.isfinite(left) or min(left, right) < 0
                or abs(left-right) > 2e-11):
            raise ValueError(path+': invalid or changed host diagnostic')
    elif type(right) not in (str, int, bool, type(None)) or left != right:
        raise ValueError(path+': value changed')


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'scope', 'core', 'decoder',
                'native_execution', 'inherited_custody', 'source_bindings',
                'physical_comparison_count_added'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.native_odd_box_core.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid odd-box source core receipt')
    for key, expected in [('scope', scope()), ('core', core.certificate()),
                          ('decoder', decoder.certificate()), ('inherited_custody', custody())]:
        if canonical(value[key]) != canonical(expected):
            raise ValueError('exact source theorem, type or boundary changed: '+key)
    if (value['issues'] != [1140, 1142, 1143]
            or value['gates'] != ['CP-04', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('finite native spectrum cannot accept a physical comparison')
    execution = value['native_execution']
    fixed = {'field_modes': 6, 'full_field_rows': 64, 'quantum_reference_columns': 64,
             'work_rows_including_two_pointer_bits': 256, 'linear_reference_column_batches': 2,
             'native_pulses': 278, 'complete_logical_pointer_outcomes': 4,
             'physical_two_pointer_outcomes': 36, 'ideal_zero_reserved_pointer_outcomes': 32,
             'native_first_failure_slots': 980,
             'reference_columns_dephased': False,
             'success_postselected_or_renormalized': False,
             'all_first_failure_weights_retained_during_execution': True}
    if (type(execution) is not dict or set(execution) != set(fixed)|DIAGNOSTICS
            or canonical({key: execution[key] for key in fixed}) != canonical(fixed)):
        raise ValueError('complete native carrier or output scope changed')
    for key in DIAGNOSTICS:
        if type(execution[key]) is not float or not 0 <= execution[key] < 2e-11:
            raise ValueError('invalid complete native host diagnostic')
    if value['source_bindings'] != source_bindings():
        raise ValueError('odd-box source manifest is incomplete or stale')
    if replay is not None:
        compare(value, replay)


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after odd-box producer import')
    value = {'schema': 'oph.native_odd_box_core.v1', 'all_pass': True,
             'issues': [1140, 1142, 1143], 'gates': ['CP-04', 'CP-06', 'CP-07'],
             'scope': scope(), 'core': core.certificate(), 'decoder': decoder.certificate(),
             'native_execution': decoder.native_execution(), 'inherited_custody': custody(),
             'physical_comparison_count_added': 0, 'source_bindings': pins}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('sources changed during odd-box replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-odd-core-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.native_odd_box_core.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold odd-box replay')
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
                          'actual_source_finite_charged_gap_positive': True,
                          'physical_species_identified': False, 'completion_gate_accepted': False}))


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
