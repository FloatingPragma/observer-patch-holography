"""Cold source-bound replay of the conditional topological defect spectra."""
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
from ..native_mass_channels import packet as parent_packet
from . import algebra, bulk, core, lattice, lattice_bounds, sampling


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = 'model/topological_defect_spectrum'
OUTPUT = Path(__file__).with_name('results.json')
PRIOR = 'model/native_mass_channels/results.json'
SECTIONS = ('algebra', 'bulk', 'lattice', 'lattice_bounds', 'core', 'sampling')


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': PACKAGE+'/README.md',
        'source_receipts': [PRIOR], 'independent_checks': [],
    })
    paths.update(PACKAGE+'/'+name for name in
                 ('PROOF.md', 'CORE_PROOF.md', 'spectral_witness.json'))
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in sorted(paths)}


def inherited_custody():
    prior = strict_json.load(ROOT/PRIOR)
    if (type(prior) is not dict
            or prior.get('schema') != 'oph.native_mass_channels.receipt.v1'
            or prior.get('all_pass') is not True
            or prior.get('integrated_gate_accepted') is not False
            or type(prior.get('physical_comparison_count_added')) is not int
            or prior['physical_comparison_count_added'] != 0):
        raise ValueError('the complete native mass source receipt is required')
    if prior.get('source_bindings') != parent_packet.source_bindings():
        raise ValueError('native mass source manifest is incomplete or stale')
    for path, digest in prior['source_bindings'].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('native mass source custody changed: '+path)
    return {'source_receipt': PRIOR,
            'complete_inherited_source_manifest_verified': True,
            'every_inherited_source_digest_verified': True,
            'inherited_centered_ground_proof_reexecuted': False}


def scope():
    return {'supplied_mass_texture_required': True,
            'all_eight_source_valleys_retained': True,
            'constant_magnitude_radial_profile_required_for_Coulomb_levels': True,
            'full_Fock_virtual_states_retained_in_finite_core_bound': True,
            'bounded_hedgehog_lattice_spectral_inclusion_proved': True,
            'additional_bounded_lattice_gap_levels_excluded': False,
            'bounded_lattice_exact_zero_modes_proved': False,
            'physical_species_identified': False,
            'source_selected_texture_proved': False,
            'physical_energy_units_calibrated': False,
            'autonomous_defect_rest_energy_derived': False,
            'common_massless_bulk_and_gapped_defect_regime_constructed': False,
            'interacting_relativistic_continuum_constructed': False,
            'native_repair_protection_proved': False,
            'Einstein_dynamics_constructed': False,
            'two_to_54_Fock_matrix_executed': False,
            'completion_gate_accepted': False}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value, replay=None):
    required = {'schema', 'all_pass', 'issues', 'gates', 'scope',
                'physical_comparison_count_added', 'inherited_custody',
                'source_bindings', *SECTIONS}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.topological_defect_spectrum.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid topological defect spectrum receipt')
    if canonical(value['scope']) != canonical(scope()):
        raise ValueError('conditional texture premises or physical scope changed')
    if (value['issues'] != [1139, 1142, 1143]
            or value['gates'] != ['CP-03', 'CP-06', 'CP-07']
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0):
        raise ValueError('the conditional spectrum cannot accept a physical comparison')
    expected_custody = {'source_receipt': PRIOR,
                        'complete_inherited_source_manifest_verified': True,
                        'every_inherited_source_digest_verified': True,
                        'inherited_centered_ground_proof_reexecuted': False}
    if canonical(value['inherited_custody']) != canonical(expected_custody):
        raise ValueError('the inherited source verification boundary changed')
    if (value['algebra']['clifford']['source_folded_dimension'] != 16
            or value['algebra']['clifford']['weyl_valleys_retained'] != 8
            or value['bulk']['Coulomb_zero_multiplicity'] != 2
            or len(value['lattice_bounds']['positive_levels']) != 6
            or value['lattice']['full_exact_lattice_zero_mode_dimension'] != 2
            or value['sampling']['source_valleys_retained'] != 8
            or value['sampling']['all_continuum_subgap_levels_approached_by_lattice_discrete_spectrum'] is not True
            or value['sampling']['additional_lattice_gap_levels_excluded'] is not False
            or value['core']['virtual_nonzero_modes_retained_in_error_bound'] is not True
            or value['core']['controlled_27_site_cube']['charged_gap_certified_positive'] is not True):
        raise ValueError('complete source spectrum or finite interacting core was lost')
    if replay is not None and canonical(value) != canonical(replay):
        raise ValueError('defect spectrum receipt differs from the fresh replay')
    if value['source_bindings'] != source_bindings():
        raise ValueError('defect spectrum source manifest is incomplete or stale')


def _produce():
    pins = source_bindings()
    custody = inherited_custody()
    value = {'schema': 'oph.topological_defect_spectrum.v1', 'all_pass': True,
             'issues': [1139, 1142, 1143], 'gates': ['CP-03', 'CP-06', 'CP-07'],
             'algebra': algebra.certificate(), 'bulk': bulk.certificate(),
             'lattice': lattice.certificate(), 'lattice_bounds': lattice_bounds.certificate(),
             'core': core.certificate(), 'sampling': sampling.certificate(), 'scope': scope(),
             'physical_comparison_count_added': 0,
             'inherited_custody': custody, 'source_bindings': pins}
    validate(value)
    if source_bindings() != pins:
        raise ValueError('source changed during the defect spectrum replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-topological-spectrum-') as temporary:
        path = Path(temporary)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(temporary)/'pycache'),
                        '-m', 'model.topological_defect_spectrum.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source changed across the cold defect spectrum replay')
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
                          'topological_zero_modes': 2, 'finite_interacting_core_gap_positive': True,
                          'physical_species_identified': False, 'completion_gate_accepted': False}))


if __name__ == '__main__':
    main()
