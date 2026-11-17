"""Cold source-bound exact inversion index and localization receipt."""
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
from ..topological_defect_spectrum import packet as parent_packet
from . import algebra, bounds


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/bounded_defect_index'
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/topological_defect_spectrum/results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [PRIOR], 'independent_checks': [],
    })
    paths.add(PACKAGE+'/PROOF.md')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in sorted(paths)}


def inherited_custody():
    value = strict_json.load(ROOT/PRIOR)
    if (type(value) is not dict or value.get('schema') != 'oph.topological_defect_spectrum.v1'
            or value.get('all_pass') is not True
            or value.get('source_bindings') != parent_packet.source_bindings()):
        raise ValueError('the complete unchanged topological spectrum source is required')
    for path, digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('inherited source custody changed: '+path)
    return {'source_receipt': PRIOR, 'complete_parent_manifest_verified': True,
            'every_inherited_source_digest_verified': True,
            'parent_spectrum_proof_reexecuted': False}


def scope():
    return {'same_supplied_bounded_radial_texture_required': True,
            'site_centered_nonlocal_spatial_inversion_grading_required': True,
            'all_source_valleys_retained': True,
            'exact_lattice_zero_mode_lower_bound_proved': True,
            'small_spacing_domain_required_for_uniform_physical_localization': True,
            'reflection_index_identified_with_internal_chiral_index_or_degree': False,
            'exact_kernel_exhaustion_proved': False,
            'native_repair_protection_proved': False,
            'source_selected_texture_proved': False,
            'interacting_many_body_zero_modes_constructed': False,
            'arbitrary_nonlocal_compact_perturbation_localization_proved': False,
            'physical_species_identified': False,
            'physical_energy_units_calibrated': False,
            'interacting_relativistic_continuum_constructed': False,
            'Einstein_dynamics_constructed': False,
            'completion_gate_accepted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'index', 'localization',
                'scope', 'physical_comparison_count_added', 'inherited_custody', 'source_bindings'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.bounded_defect_index.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid bounded defect index receipt')
    if canonical(value['scope']) != canonical(scope()):
        raise ValueError('the supplied texture or symmetry domain was promoted')
    if (value['issues'] != [1139, 1142, 1143] or value['gates'] != ['CP-03', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('an operator index cannot accept a physical comparison')
    if (value['index']['exact_Fredholm_index'] != 2
            or value['index']['exact_normalizable_zero_mode_lower_bound'] != 2
            or value['index']['atomic']['reflection_index_antihedgehog'] != 2
            or value['index']['exact_kernel_exhaustion_proved'] is not False
            or value['index']['native_repair_preserves_required_operator_symmetry_proved'] is not False
            or value['localization']['localization_extended_to_arbitrary_nonlocal_compact_perturbations'] is not False):
        raise ValueError('the exact index, localization or admitted perturbation class changed')
    expected_custody = {'source_receipt': PRIOR, 'complete_parent_manifest_verified': True,
                        'every_inherited_source_digest_verified': True,
                        'parent_spectrum_proof_reexecuted': False}
    if canonical(value['inherited_custody']) != canonical(expected_custody):
        raise ValueError('the inherited proof verification boundary changed')
    if replay is not None and canonical(value) != canonical(replay):
        raise ValueError('the bounded index receipt differs from the fresh replay')
    if value['source_bindings'] != source_bindings():
        raise ValueError('the bounded index source manifest is incomplete or stale')


def _produce():
    pins = source_bindings()
    custody = inherited_custody()
    value = {'schema': 'oph.bounded_defect_index.v1', 'all_pass': True,
             'issues': [1139, 1142, 1143], 'gates': ['CP-03', 'CP-06', 'CP-07'],
             'index': algebra.certificate(), 'localization': bounds.certificate(),
             'scope': scope(), 'physical_comparison_count_added': 0,
             'inherited_custody': custody, 'source_bindings': pins}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during the index replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-bounded-defect-index-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.bounded_defect_index.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold index replay')
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
                          'reflection_Fredholm_index': 2, 'exact_zero_mode_lower_bound': 2,
                          'native_repair_protection': False, 'physical_species_identified': False}))


if __name__ == '__main__':
    main()
