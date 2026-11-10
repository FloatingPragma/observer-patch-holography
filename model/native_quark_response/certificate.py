"""Rebuild native quark spectral images, residues and full-source records."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

from model.particle_common_contract.provenance import import_closure
from .quark import copy_alignment, target_module
from .local import CHANNELS
from .response import certificate as local_certificate, positive_response
from .full import windows


PARENTS = {
    'model/coherent_defect_wave/chiral_joint/local_reference/results.json':
        '0f8db24b8e5ed7312f39484247a482be216bbd6944ab3b56cf0929ba8b4e138c',
    'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json':
        'e135f7b15460427068dd18541b0f7990a3d259d298467604283841c68f75ab34',
    'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json':
        '7335d80394ab25ab949ca0e359e70a8eca55636256255b684c68086be2a8c7a6',
    'model/native_charge_resolution/results.json':
        'f22c9a266a43de8aaae776b1293894081c3dc20b7d8e5678e30ef6af993dd19e',
    'model/coherent_defect_wave/chiral_joint/interacting_charge_conjugation/results.json':
        'f5852be27204dde1bc46df2517d476e737506fc36758c3497536f37ff00a3f64',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(pins):
    for name, expected in pins.items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('repository-relative source pin required')
        if digest(ROOT/name) != expected:
            raise ValueError('frozen source changed: '+name)


def calculate():
    from fractions import Fraction
    frames = copy_alignment()
    modules = [target_module(create) for create in (False, True)]
    local = {family: local_certificate(family) for family in CHANNELS}
    uc, dc = positive_response('uc'), positive_response('dc')
    assert uc['numerator'] == dc['numerator'] and uc['denominator'] == dc['denominator']
    native = {family: windows(family) for family in CHANNELS}
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import operator_certificate
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.model import separate_spin
    spin = separate_spin()
    return {
        'schema': 'oph.native_quark_response.v1',
        'same_source': {
            'local': 'full native thirty-CAR plus five-state scalar, all Yukawa coefficients one',
            'local_reference': 'actual physical ground projector divided by two',
            'interacting': 'h_active+h_partner+4N_ref, h=H0+2^-34 T on the full native eight-cell carrier',
            'full_source_dimension': 'D^2, D=2^240 5^8 119^24',
            'static_reference_dimension': 2**30,
            'original_physical_preparation_retained': True,
            'additional_exact_ground_projection': False,
            'all_failure_and_no_operation_records_retained': True,
        },
        'actual_quark_CAR_images': {
            'six_native_components': [list(rows) for rows, _ in frames],
            'signed_unitary_alignment_matrices': [[[str(x) for x in row] for row in u.tolist()]
                                                   for _, u in frames],
            'source_harmonic_dimension': 140,
            'source_Casimir': '54',
            'target_ambient_dimension': 504,
            'target_image_dimensions': [m['image_dimension'] for m in modules],
            'target_Casimir': '95/2',
            'transformed_target_numbers': [m['number'] for m in modules],
            'target_metric_split_factor': '1/6',
            'exact_full_image_rank_checks_by_every_weak_hypercharge_weight': True,
            'ten_generator_CAR_covariance_identities_per_mode_and_branch': True,
            'sum_four_elementary_effects': '2 I_140, for each transformed creation/removal branch',
            'actual_quark_hypercharge_offsets': {'creation': '1/6', 'removal': '-1/6'},
            'full_native_separate_quark_lepton_spin_H_commutators':
                spin['independent_quark_and_lepton_spin_H_commutators'],
        },
        'full_probe_reference_operator': operator_certificate(),
        'exact_local_response': local,
        'computed_equal_coupling_uc_dc_resolvent_equality': True,
        'full_interacting_prepared_records': native,
        'apparatus': {
            'pump_duration_source_units': 1,
            'energy_duration_source_units': 1,
            'autonomous_H_on_every_branch': 'h_active+h_partner+4N_ref',
            'energy_service': 'exact full-source spectral bins, complete complementary pointer outcome',
            'unitary_service_completion': 'V(t)=exp(-itH)Lexp(itH), exp(-iL)=U; full evolution exp(-iH)U',
            'work_upper_per_completed_service': '2 ||H_total||',
            'full_H_total_norm_upper': str(Fraction(1080)+Fraction(720, 2**34)),
            'bare_source_work_upper_per_completed_service': str(Fraction(2160)+Fraction(1440, 2**34)),
            'retained_storage_upper': 'D^2 2^30 2^32 Mprep^2 (degree(q)+1) 2^3',
            'Mprep': '1197721483502453176152084669',
            'controller_synthesis_workspace_included': False,
            'efficient_native_gate_synthesis_claimed': False,
        },
        'scope': {
            'all_local_probe_roots_and_exact_zero_residues_retained': True,
            'full_compensated_charged_ground_identified': False,
            'quark_propagation_or_asymptotic_color_states': False,
            'rest_mass_or_GeV_calibration': False,
            'computed_degeneracy_is_physical_up_down_mass_equality': False,
            'Standard_Model_matching': False,
        },
    }


def raw_certificate():
    pins = {}
    saved = {}
    for name, expected in PARENTS.items():
        if digest(ROOT/name) != expected:
            raise ValueError('frozen parent changed: '+name)
        value = json.loads((ROOT/name).read_text())
        validate(value['input_sha256'])
        saved[name] = value
        for path, sha in value['input_sha256'].items():
            if path in pins and pins[path] != sha:
                raise ValueError('inconsistent frozen import closure')
            pins[path] = sha
        pins[name] = expected
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.certificate import certificate as charged
    from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.certificate import certificate as prep
    from model.native_charge_resolution.bounds import calculate as refined
    from model.coherent_defect_wave.chiral_joint.interacting_charge_conjugation.certificate import certificate as conjugation
    assert charged() == saved[list(PARENTS)[1]]
    assert prep() == saved[list(PARENTS)[2]]
    conjugate = conjugation()
    assert conjugate == saved[list(PARENTS)[4]]
    assert conjugate['same_source_spectral_consequences']['all_five_ground_family_pump_failure_probabilities'] == '1/2'
    resolution = refined()
    for key, value in resolution.items():
        if key == 'schema':
            continue  # The public parent renames only its provenance schema.
        assert saved[list(PARENTS)[3]][key] == value
    value = calculate()
    value['prerequisite_receipts'] = PARENTS
    seeds = set(pins) | {
        'model/coherent_defect_wave/chiral_joint/local_reference/algebra.py',
        'model/coherent_defect_wave/chiral_joint/local_reference/reduction.py',
        'model/coherent_defect_wave/chiral_joint/local_reference/physical.py',
        'model/coherent_defect_wave/chiral_joint/ground_composite_response/model.py',
        'model/native_charge_resolution/bounds.py',
    }
    for name in import_closure(seeds):
        pins[name] = digest(ROOT/name)
    value['input_sha256'] = dict(sorted(pins.items()))
    value['owned_source_sha256'] = {p.name: digest(p) for p in sorted(HERE.iterdir())
                                    if p.suffix in ('.py', '.md')}
    value['fresh_replay'] = {
        'actual_source_reference_and_physical_preparation_recomputed': True,
        'all_full_CAR_images_and_Hamiltonian_inclusions_recomputed': True,
        'all_target_roots_and_number_field_residues_recomputed': True,
        'full_interacting_probability_bounds_recomputed': True,
        'numerical_eigensolver_used': False,
        'large_full_source_pointer_matrix_executed': False,
    }
    return value


REFERENCE_SHA = '48027932e62481a7c33d6ed7ba34afffc1497aa83cef4ca9e5011f9ff19d4a17'


def certificate():
    value = raw_certificate()
    path = HERE/'research_reference.json'
    if digest(path) != REFERENCE_SHA:
        raise ValueError('immutable research receipt changed')
    reference = json.loads(path.read_text())
    validate(reference['input_sha256'])
    ignored = {'input_sha256', 'owned_source_sha256'}
    if {k: v for k, v in value.items() if k not in ignored} != {
            k: v for k, v in reference.items() if k not in ignored}:
        raise ValueError('fresh exact quark science differs from research receipt')
    value.pop('owned_source_sha256')
    seeds = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in ('.py', '.md')}
    seeds.add(path.relative_to(ROOT).as_posix())
    for name in import_closure(seeds):
        value['input_sha256'][name] = digest(ROOT/name)
    value['input_sha256'] = dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison'] = {
        'reference_sha256': REFERENCE_SHA,
        'every_scientific_field_equal': True,
        'ignored_only_path_and_owned_file_hash_metadata': sorted(ignored),
        'fresh_full_CAR_images_residues_and_interacting_records_recomputed': True,
    }
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    path = HERE/'results.json'
    if args.check:
        saved = json.loads(path.read_text())
        validate(saved['input_sha256'])
        for name, expected in saved.get('owned_source_sha256', {}).items():
            if digest(HERE/name) != expected:
                raise ValueError('owned source changed: '+name)
    value = certificate()
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    elif value != saved:
        raise ValueError('fresh complete native quark response differs')
    print('PASS: exact native quark CAR images, local measures and complete prepared eight-cell bins')


if __name__ == '__main__':
    main()
