"""Cold replay of complete finite charged spectral measures and currents."""
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
from ..topological_defect_spectrum import packet as textured
from ..native_odd_box_core import packet as original
from . import algebra, response


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/defect_core_response'
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = {'topological_defect_spectrum': textured, 'native_odd_box_core': original}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PRIORS],
        'independent_checks': [],
    })
    paths.add(PACKAGE+'/PROOF.md')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def inherited_custody():
    result = {}
    for name, producer in PRIORS.items():
        parent = strict_json.load(ROOT/f'model/{name}/results.json')
        if (parent.get('schema') != 'oph.'+name+'.v1'
                or parent.get('all_pass') is not True
                or parent.get('source_bindings') != producer.source_bindings()
                or parent.get('scope', {}).get('completion_gate_accepted') is not False):
            raise ValueError('missing, stale or promoted finite source: '+name)
        for path, digest in parent['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('inherited source digest changed: '+path)
        result[name] = {'complete_current_source_manifest_verified': True,
                        'every_inherited_digest_verified': True,
                        'parent_expensive_receipt_reexecuted': False}
    return result


def scope():
    return {
        'two_distinct_finite_source_domains_retained': True,
        'textured_branch_requires_supplied_mass_field': True,
        'original_odd_box_requires_no_added_mass_field': True,
        'actual_interacting_ground_and_charged_projectors_used': True,
        'full_Fock_complement_retained': True,
        'positive_native_CAR_pole_residues_proved': True,
        'complete_number_and_energy_currents_retained': True,
        'all_native_coin_raw_and_event_outputs_retained': True,
        'no_event_conditioning_success_probability_retained': True,
        'explicit_positive_rate_continuing_window_proved': True,
        'continuing_window_rate_is_chosen_model_parameter': True,
        'two_to_54_matrix_executed': False,
        'long_time_spectral_detector_executed': False,
        'ground_sea_preparation_executed': False,
        'native_repair_protects_the_spectral_poles': False,
        'pole_matrix_rank_is_physical_species_count': False,
        'thermodynamic_particle_mass_proved': False,
        'physical_species_and_energy_units_identified': False,
        'interacting_relativistic_continuum_constructed': False,
        'Einstein_stress_or_regional_entropy_premises_supplied': False,
        'completion_gate_accepted': False,
    }


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'response', 'currents',
                'scope', 'physical_comparison_count_added', 'inherited_custody',
                'source_bindings'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.defect_core_response.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid finite core response receipt')
    if canonical(value['scope']) != canonical(scope()):
        raise ValueError('spectral-response domain or physical scope changed')
    if (value['issues'] != [1140, 1142, 1143]
            or value['gates'] != ['CP-04', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('finite pole evidence cannot accept physical completion')
    if value['source_bindings'] != source_bindings():
        raise ValueError('response source manifest is incomplete or stale')
    for key, expected in [('response', response.certificate()),
                          ('currents', algebra.certificate()),
                          ('inherited_custody', inherited_custody())]:
        if canonical(value[key]) != canonical(expected):
            raise ValueError('finite source result or custody changed: '+key)
    if replay is not None and canonical(value) != canonical(replay):
        raise ValueError('response receipt differs from fresh replay')


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after response producer import')
    value = {'schema': 'oph.defect_core_response.v1', 'all_pass': True,
             'issues': [1140, 1142, 1143], 'gates': ['CP-04', 'CP-06', 'CP-07'],
             'response': response.certificate(), 'currents': algebra.certificate(),
             'scope': scope(), 'physical_comparison_count_added': 0,
             'inherited_custody': inherited_custody(), 'source_bindings': pins}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during response replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-defect-core-response-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.defect_core_response.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold response replay')
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
        print('finite charged response: PASS; source domains and complete outputs retained')


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
