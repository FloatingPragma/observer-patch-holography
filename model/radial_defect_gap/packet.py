"""Cold source-bound exact radial kernel and spectral-gap receipt."""
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
from ..bounded_defect_index import packet as parent_packet
from . import algebra, bounds


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/radial_defect_gap'
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/bounded_defect_index/results.json'


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
    if (type(value) is not dict or value.get('schema') != 'oph.bounded_defect_index.v1'
            or value.get('all_pass') is not True
            or value.get('source_bindings') != parent_packet.source_bindings()):
        raise ValueError('the complete unchanged bounded radial source is required')
    for path, digest in value['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('inherited source custody changed: '+path)
    return {'source_receipt': PRIOR, 'complete_parent_manifest_verified': True,
            'every_inherited_source_digest_verified': True,
            'parent_index_proof_reexecuted': False}


def scope():
    return {'same_supplied_bounded_radial_texture_required': True,
            'strong_mass_domain_required': True,
            'all_source_valleys_retained': True,
            'exact_kernel_dimension_two_proved': True,
            'volume_uniform_one_body_gap_proved': True,
            'finite_volume_kernel_projector_convergence_quantified': True,
            'local_zero_frame_and_origin_observable_join_proved': True,
            'origin_observable_identity_is_a_free_reference_identity': True,
            'fixed_physical_mass_continuum_refinement_proved': False,
            'full_bound_spectrum_outside_the_central_gap_classified': False,
            'interacting_many_body_gap_proved_in_this_package': False,
            'source_selected_texture_proved': False,
            'native_repair_protection_proved': False,
            'physical_species_identified': False,
            'physical_energy_units_calibrated': False,
            'interacting_relativistic_continuum_constructed': False,
            'Einstein_dynamics_constructed': False,
            'completion_gate_accepted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'algebra', 'bounds',
                'scope', 'physical_comparison_count_added', 'inherited_custody', 'source_bindings'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.radial_defect_gap.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid radial defect gap receipt')
    if canonical(value['scope']) != canonical(scope()):
        raise ValueError('the supplied texture or strong mass domain was promoted')
    if (value['issues'] != [1139, 1142, 1143] or value['gates'] != ['CP-03', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('the one-body theorem cannot accept a physical comparison')
    if (canonical(value['algebra']) != canonical(algebra.certificate())
            or canonical(value['bounds']) != canonical(bounds.certificate())):
        raise ValueError('the exact source algebra or admitted spectral bounds changed')
    expected_custody = {'source_receipt': PRIOR, 'complete_parent_manifest_verified': True,
                        'every_inherited_source_digest_verified': True,
                        'parent_index_proof_reexecuted': False}
    if canonical(value['inherited_custody']) != canonical(expected_custody):
        raise ValueError('the inherited proof verification boundary changed')
    if replay is not None and canonical(value) != canonical(replay):
        raise ValueError('the radial gap receipt differs from the fresh replay')
    if value['source_bindings'] != source_bindings():
        raise ValueError('the radial gap source manifest is incomplete or stale')


def _produce():
    pins = source_bindings()
    custody = inherited_custody()
    value = {'schema': 'oph.radial_defect_gap.v1', 'all_pass': True,
             'issues': [1139, 1142, 1143], 'gates': ['CP-03', 'CP-06', 'CP-07'],
             'algebra': algebra.certificate(), 'bounds': bounds.certificate(),
             'scope': scope(), 'physical_comparison_count_added': 0,
             'inherited_custody': custody, 'source_bindings': pins}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during the radial gap replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-radial-defect-gap-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.radial_defect_gap.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold radial gap replay')
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
                          'exact_kernel_dimension': 2, 'volume_uniform_one_body_gap': True,
                          'fixed_physical_mass_continuum_refinement': False,
                          'physical_species_identified': False}))


if __name__ == '__main__':
    main()
