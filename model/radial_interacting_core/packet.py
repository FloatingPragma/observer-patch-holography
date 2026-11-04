"""Cold custody-bound replay of the same bounded radial interacting source."""
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
from ..radial_defect_gap import packet as gap_packet
from ..defect_core_response import packet as pole_packet
from . import algebra, core


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/radial_interacting_core'
OUTPUT = Path(__file__).with_name('results.json')
PRIORS = {'radial_defect_gap': gap_packet, 'defect_core_response': pole_packet}


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [f'model/{name}/results.json' for name in PRIORS],
        'independent_checks': [],
    })
    paths.add(PACKAGE+'/PROOF.md')
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def custody():
    result = {}
    for name, producer in PRIORS.items():
        path = f'model/{name}/results.json'
        value = strict_json.load(ROOT/path)
        if (type(value) is not dict or value.get('schema') != 'oph.'+name+'.v1'
                or value.get('all_pass') is not True
                or value.get('source_bindings') != producer.source_bindings()):
            raise ValueError('complete current parent proof receipt required: '+name)
        for source, digest in value['source_bindings'].items():
            if hashlib.sha256((ROOT/source).read_bytes()).hexdigest() != digest:
                raise ValueError('inherited source binding changed: '+source)
        result[name] = {'receipt': path, 'complete_parent_source_manifest_verified': True,
                        'every_parent_source_digest_verified': True,
                        'expensive_parent_controls_reexecuted': False}
    return result


def scope():
    return {'same_bounded_radial_one_body_and_interacting_Hamiltonian': True,
            'supplied_bounded_mass_texture_required': True,
            'strong_mass_domain_required': True,
            'complete_Fock_complement_retained': True,
            'exact_actual_ground_doublet_and_positive_charged_gap': True,
            'same_source_positive_Lehmann_pole_weights': True,
            'actual_local_center_parity_contrast_proved': True,
            'actual_center_parity_is_exact_zero_orbital_occupation': False,
            'faithful_finite_interacting_thermal_reference_constructed': True,
            'finite_modular_clock_equals_negative_beta_Hamiltonian_time': True,
            'source_selects_or_prepares_radial_mass_texture': False,
            'Gibbs_preparation_or_repair_equilibration_executed': False,
            'two_to_54_Fock_matrix_executed': False,
            'volume_uniform_interacting_gap_proved': False,
            'strong_mass_bound_is_fixed_mass_spatial_continuum_limit': False,
            'native_repair_protection_of_poles_proved': False,
            'physical_species_identified': False,
            'interacting_relativistic_continuum_constructed': False,
            'responsive_Einstein_geometry_instantiated': False,
            'completion_gate_accepted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'algebra', 'core', 'scope',
                'inherited_custody', 'source_bindings', 'physical_comparison_count_added'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.radial_interacting_core.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid same-source radial interacting receipt')
    for key, expected in [('scope', scope()), ('algebra', algebra.certificate()),
                          ('core', core.certificate()), ('inherited_custody', custody())]:
        if canonical(value[key]) != canonical(expected):
            raise ValueError('source theorem, type or boundary changed: '+key)
    if (value['issues'] != [1139, 1142, 1143]
            or value['gates'] != ['CP-03', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('finite radial source cannot accept a physical comparison')
    if value['source_bindings'] != source_bindings():
        raise ValueError('radial source manifest is incomplete or stale')
    if replay is not None and canonical(value) != canonical(replay):
        raise ValueError('retained radial source receipt differs from the cold replay')


def _produce():
    pins = source_bindings()
    if pins != _IMPORT_BINDINGS:
        raise ValueError('source changed after radial interacting producer import')
    value = {'schema': 'oph.radial_interacting_core.v1', 'all_pass': True,
             'issues': [1139, 1142, 1143], 'gates': ['CP-03', 'CP-06', 'CP-07'],
             'algebra': algebra.certificate(), 'core': core.certificate(), 'scope': scope(),
             'inherited_custody': custody(), 'source_bindings': pins,
             'physical_comparison_count_added': 0}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during radial interacting replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-interaction-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(folder)/'pycache'),
                        '-m', 'model.radial_interacting_core.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across cold radial interacting replay')
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
                          'same_bounded_radial_interacting_gap': True,
                          'physical_identification': False, 'completion_gate_accepted': False}))


_IMPORT_BINDINGS = source_bindings()

if __name__ == '__main__':
    main()
