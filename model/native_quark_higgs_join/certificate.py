"""Replay complete native quark/Higgs operators and finite interference records."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
from model.native_quark_response.certificate import certificate as parent_certificate
from model.particle_common_contract.provenance import import_closure
from .operators import operator_controls, reference_controls, reduced_EOM
from .bounds import calculate

PARENT = 'model/native_quark_response/results.json'
PARENT_SHA = '06c336e1f691db866896de6123fea34937dcd462d92765f594600744c54fd40f'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(pins):
    for name, expected in pins.items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts or digest(ROOT/path) != expected:
            raise ValueError('frozen source pin changed: '+name)


def raw_certificate():
    if digest(ROOT/PARENT) != PARENT_SHA:
        raise ValueError('native quark parent changed')
    parent = json.loads((ROOT/PARENT).read_text())
    validate(parent['input_sha256'])
    assert parent_certificate() == parent
    result = dict(schema='oph.native_quark_higgs_join.v1',
                  source_parent={PARENT: PARENT_SHA},
                  exact_original_CAR_operators=operator_controls(),
                  full_static_reference_controls=reference_controls(),
                  complete_target_Higgs_images=[reduced_EOM(f) for f in ('uc', 'dc', 'Q')],
                  full_interacting_interference_records=calculate(),
                  same_source=dict(H='h_active+h_partner+4N_ref', native_h='H0+2^-34 T',
                                   native_D='2^240 5^8 119^24', source_dimension='D^2',
                                   full_reference_dimension=2**30, reference_cost_per_occupation=4,
                                   original_physical_preparation=True, additional_P0=False,
                                   all_complements_and_failure_records_retained=True),
                  apparatus=dict(pump_source_duration=1, energy_source_duration=1,
                                 full_source_clock_on_no_operation_branches=True,
                                 service='finite programmable Halmos pump and exact spectral-bin readout',
                                 interaction_picture='V(t)=exp(-itH)Lexp(itH), exp(-iL)=U',
                                 full_service_evolution='exp(-iH)U',
                                 norm_H_upper='1080+720/2^34', work_per_service_upper='2160+1440/2^34',
                                 retained_storage_upper='D^2 2^30 2^32 Mprep^2 2^4',
                                 Mprep=1197721483502453176152084669,
                                 full_setting_pump_energy_labels=8,
                                 potentially_nonzero_successful_preparation_labels=6,
                                 controller_synthesis_workspace_included=False,
                                 efficient_native_gate_synthesis=False),
                  scope=dict(actual_shared_positive_spectral_bin_contribution=True,
                             all_atomic_roots_and_zero_residues_retained=True,
                             classical_scalar_background=False,
                             canonical_composite_CAR=False,
                             full_compensated_charged_ground_identified=False,
                             propagation_or_asymptotic_colored_particle=False,
                             physical_Lorentz_chirality_or_Standard_Model_matching=False,
                             source_energies_are_calibrated_rest_masses=False,
                             local_equal_coupling_uc_dc_resolvent_identity_is_physical_mass_equality=False),
                  fresh_replay=dict(parent_actual_producer_recomputed=True,
                                    all_twenty_four_original_CAR_EOM_and_currents_recomputed=True,
                                    all_sixty_CAR_reference_commutators_recomputed=True,
                                    literal_full_target_image_and_H_EOM_recomputed=True,
                                    all_atomic_weights_and_actual_source_errors_recomputed=True,
                                    numerical_eigensolver_used=False,
                                    enormous_full_source_pointer_matrix_executed=False))
    pins = dict(parent['input_sha256'])
    pins[PARENT] = PARENT_SHA
    seeds = set(pins) | {'model/coherent_defect_wave/chiral_joint/higgs_composites/model.py'}
    for name in import_closure(seeds):
        pins[name] = digest(ROOT/name)
    result['input_sha256'] = dict(sorted(pins.items()))
    result['owned_source_sha256'] = {p.name: digest(p) for p in sorted(HERE.iterdir())
                                     if p.suffix in ('.py', '.md')}
    return result


REFERENCE_SHA = 'd724a495ce3e114f05fe0832a94f7b868aacbcea48b4f962b2863260eec402d8'


def certificate():
    value = raw_certificate()
    reference_path = HERE/'research_reference.json'
    if digest(reference_path) != REFERENCE_SHA:
        raise ValueError('immutable research receipt changed')
    reference = json.loads(reference_path.read_text())
    validate(reference['input_sha256'])
    ignored = {'input_sha256', 'owned_source_sha256'}
    if {k: v for k, v in value.items() if k not in ignored} != {
            k: v for k, v in reference.items() if k not in ignored}:
        raise ValueError('fresh exact quark/Higgs science differs from research receipt')
    from .alignment import certificate as admitted_states
    value['normalized_admitted_states_and_histories'] = admitted_states()
    value.pop('owned_source_sha256')
    seeds = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in ('.py', '.md')}
    seeds.add(reference_path.relative_to(ROOT).as_posix())
    seeds.update(value['normalized_admitted_states_and_histories']['pure_inequality_inputs'])
    for name in import_closure(seeds):
        value['input_sha256'][name] = digest(ROOT/name)
    value['input_sha256'] = dict(sorted(value['input_sha256'].items()))
    value['independent_research_comparison'] = {
        'reference_sha256': REFERENCE_SHA,
        'every_research_interference_scientific_field_equal': True,
        'additional_public_theorem': 'normalized_admitted_states_and_histories',
        'ignored_only_path_and_owned_file_hash_metadata': sorted(ignored),
        'fresh_complete_operator_image_and_probability_replay': True,
    }
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif json.loads(path.read_text()) != result:
        raise SystemExit('native quark/Higgs exact replay mismatch')
    print('PASS native quark/Higgs complete same-source records', digest(path), flush=True)


if __name__ == '__main__':
    main()
