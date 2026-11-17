"""Bind the official CCF-LT-2 finite evidence without physical promotion."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
PACKAGE = ROOT / 'model/ccf_lt2'
OUTPUT = ROOT / 'claims/ccf_lt2.json'
PARENT_DIGEST = '591949feb5120a050c49a8b4edf52db9a920bfaf0db920d8c7937089b4613013'
COMPONENTS = {'parent', 'dirac', 'LL', 'quark_family', 'QQ_nonleptonic'}
MATTER_DIGESTS = {
    'source_edge': '80ea27a13e1f31629d6a9867352816334a80938e7ffc3740aede4803c721a2f2',
    'source_triangle': 'c824e41e0d7f3a98f2755f1202cf9e92067c7c479d317336b79d9a4453b23669',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_receipt(receipt, *, check_bindings=True):
    require(receipt.get('model') == 'CCF-LT-2', 'receipt belongs to another model')
    require(receipt.get('all_passed') is True and receipt.get('full_source_included') is True,
            'complete native source verification is required')
    if check_bindings:
        pins = receipt.get('scientific_inputs_sha256', {})
        require(bool(pins), 'scientific input bindings are required')
        from model.ccf_lt2.verify import scientific_inputs
        require(set(pins) == set(scientific_inputs()), 'scientific input inventory omitted')
        for name, digest in pins.items():
            path = (PACKAGE / name).resolve()
            require((path.is_relative_to(PACKAGE) or name in {'../SPEC.md', '../AUTHORITY.json'}) and path.is_file(), 'invalid scientific input path')
            require(sha256(path.read_bytes()).hexdigest() == digest, 'stale scientific input: '+name)
    for key, modes, terms in [('source_edge', 192, 63589), ('source_triangle', 288, 109675)]:
        source = receipt[key]
        require(source['all_passed'] and source['native_CAR_modes'] == modes, 'native carrier changed')
        require(source['parent_export_digest'] == PARENT_DIGEST and source['parent_factors'] == 114,
                'parent snapshot changed')
        require(source['joint'] == {'all_local_gauss_identities': True, 'hermitian': True, 'terms': terms},
                'combined native polynomial failed')
        require(set(source['components']) == COMPONENTS, 'native source component omitted')
        require(source['combined_source']['matter_digest'] == MATTER_DIGESTS[key], 'native matter polynomial changed')
        require(all(c['hermitian'] and c['all_local_gauss_identities'] for c in source['components'].values()),
                'component Hermiticity or Gauss failed')
        require(source['remove_all_additions_recovers_replicated_parent'], 'parent recovery failed')
        require(source['local_matter_shards_sum_to_full_polynomial'], 'energy allocation does not sum')
        require(source['wrong_link_mutation_rejected'], 'wrong-link control was accepted')
        require(source['phase_current_and_contact']['all_passed'], 'complete source current/contact failed')
        require(source['combined_source']['full_source_state_vector_executed'] is False,
                'basis columns cannot be promoted to arbitrary state propagation')
        require(source['physical_rows_passed_by_this_check'] == [], 'finite checks cannot award physical rows')
        mass = source['actual_parent_mass_kinetic_join']
        require(mass['mass_coefficients_unchanged'] and mass['charge_preserving_spin_summed_density_control_rejected'],
                'actual mass/kinetic bridge failed')
        require(all(k['all_three_mass_dirac_anticommutators_exact_zero'] and k['complete_nambu_closure']
                    for k in mass['kinds'].values()), 'charged or neutral mass compatibility failed')
    require(receipt['source_edge']['exact_two_cell_parent_recovery'], 'two-cell parent recovery failed')
    local = receipt['actual_local_capacity_join']
    require(receipt.get('specification') == 'A-ccf2-local-capacity-v2 finite matter, field and response construction',
            'receipt belongs to another capacity law')
    require(local['all_passed'] and local['native_modes'] == 288 and
            local['local_response_dimension'] == 27 and local['capacity_auxiliary_dimension_per_vertex'] == 16 and
            local['capacity_auxiliary_count'] == 3, 'local-capacity source carrier changed')
    require(local['matter_digest'] == MATTER_DIGESTS['source_triangle'], 'local-capacity native source changed')
    columns = local['columns']
    require([c['preparation'] for c in columns] == ['native_empty', 'native_neutral'],
            'complete native column preparations changed')
    require([c['complete_nonzero_outputs'] for c in columns] == [4567, 4564],
            'complete native output count changed')
    require(all(c['all_outputs_gauss_consistent'] and c['independent_diagonal_matches'] and
                c['no_native_occupation_projection'] and c['local_phase_current_and_contact_exact'] and
                c['exact_total_charge_continuity'] for c in columns),
            'complete local-capacity column identity failed')
    require(columns[0]['diagonal_exact'] != columns[1]['diagonal_exact'],
            'neutral energy decoupled from total response')
    require(set(local['all_source_components_retained']) == COMPONENTS | {'electric', 'magnetic'},
            'coupled source component omitted')
    for name in ['field', 'per_vertex_capacity', 'local_geometry', 'local_capacity']:
        require(receipt[name]['all_passed'], 'finite diagnostic failed: '+name)
    capacity = receipt['local_capacity']
    require(capacity['matrix_columns']['complete_columns_checked'] == 256 and
            capacity['matrix_columns']['max_column_error'] == 0, 'complete capacity matrix failed')
    require(capacity['disconnected_factorization']['exact_component_sum'] and
            capacity['disconnected_factorization']['disjoint_even_terms_commute'], 'disconnected factorization failed')
    require(capacity['reference_convention']['one_common_F'] and
            capacity['reference_convention']['nonuniform_shift_not_absorbed_by_common_F'], 'energy reference convention changed')
    require(all(row['mode'] == 'localized' and row['response'] < 1e-10 for row in capacity['clock']['records']),
            'disconnected source clock depends on remote label')
    require(receipt['local_geometry']['joint_dimension'] == 72, 'local diagnostic domain changed')
    require(abs(receipt['local_geometry']['neutral_distant_clock_response']) > 1e-4,
            'local neutral response absent')
    require(receipt['field']['full_flux_states'] == 343, 'complete field carrier changed')
    interactions = receipt['interactions']
    require(interactions['complete_native_canonical_LL_identity'] and all(interactions['controls'].values()),
            'contact interaction identity or controls failed')
    inherited = receipt['inherited_charge_and_preparation']
    require(set(inherited) == {'charge','compiler','full_fock','quantum','topology_and_repair'},
            'inherited charge or preparation service omitted')
    require(inherited['compiler']['nonneutral_input_rejected'] and
            inherited['full_fock']['charge_transform_exact'] and
            inherited['topology_and_repair']['unbalanced_sector_rejected'],
            'inherited charge/preparation integrity failed')


def validate_bridges(bridges):
    required={'model','candidate','kinematics','phase_instrument','matter','matter_columns','reference',
              'gravity','physical_rows_awarded','all_passed','scientific_inputs_sha256','elapsed_seconds'}
    require(set(bridges)==required, 'finite bridge group omitted or added')
    require(bridges['all_passed'] is True and bridges['physical_rows_awarded']==[], 'bridge physical promotion')
    kin=bridges['kinematics']
    require(kin['all_finite_algebraic_gates_passed'] is True and kin['physical_rows_claimed_solved']==[],
            'kinematics scope changed')
    require(len(kin['continuity'])==2 and all(c['all_vertex_continuity_coefficients_exact_zero'] is True and
            c['actual_currents_hermitian'] is True and c['omitted_dirac_current_control_nonzero_residual_terms']>0
            for c in kin['continuity']), 'current continuity or omission control failed')
    matter=bridges['matter']
    require(matter['all_passed'] is True and matter['source_cluster_additivity'] is True and
            matter['native_colour_singlet_checks']==8 and matter['complete_disconnected_empty_column_identity'] is True,
            'complete matter certificate failed')
    for label in ('edge','disconnected'):
        require(matter[label]['conserved']==['Q3','B3','L','fermion_parity'] and matter[label]['hermitian'] is True,
                'native sector certificate changed')
    ref=bridges['reference']
    require(ref['complete_fixture']['all_passed'] is True and ref['complete_fixture']['wrong_gauss_mutation_rejected'] is True,
            'regional reference support fixture failed')
    members=ref['actual_source']['members']
    require(len(members)==2, 'complete reference member omitted')
    for v,row in zip((2,3),members):
        require(row['beta_exact']=='1/1024' and row['geometry_dimension']==48**v and
                int(row['complete_physical_rank'])==int(row['matter_field_gauss_rank'])*48**v and
                row['physical_gauss_support_nonempty'] is True, 'full Gauss reference changed')
    gravity=bridges['gravity'];coframe=gravity['coframe'];mesh=gravity['refinement'];ward=gravity['energy_ward']
    require(coframe['identity_member_exact'] is True and len(coframe['strain_derivatives'])==6 and
            len(coframe['shift_derivatives'])==3 and coframe['omitted_wilson_factor_detected'] is True,
            'same-source mechanical coframe failed')
    require([r['level'] for r in mesh['members']]==[0,1,2,3] and
            all(r['exact_boundary_squared_zero'] is True and r['reversed_face_detected'] is True for r in mesh['members']),
            'shrinking mesh certificate failed')
    require(ward['all_passed'] is True and ward['complete_local_capacity_columns_checked']==256 and
            ward['continuity_residual']==ward['antisymmetry_residual']==0 and ward['omitted_current_residual']>0,
            'local-capacity energy Ward certificate failed')


def generate():
    receipt_path = PACKAGE / 'RESULTS.json'
    receipt = json.loads(receipt_path.read_text())
    validate_receipt(receipt)
    from model.ccf_lt2.bridges.verify import validate_saved
    bridges = validate_saved()
    validate_bridges(bridges)
    phase = bridges['phase_instrument']
    require(phase['local_capacity_entry_exact'] == '9/32' and
            phase['initial_probability_derivatives_exact'] == ['-9/32', '9/32', '0'] and
            phase['complete_PVM_exact'] and phase['dirac_deleted_same_entry_exact'] == '0',
            'complete current phase instrument failed')
    require(bridges['reference']['stationarity_boundary']['entropy_plus_area_first_variation_exact'] == '9/16384',
            'stationarity countervariation changed')
    evidence = ['model/ccf_lt2/bridges/RESULTS.json', 'model/ccf_lt2/bridges/PROOF.md', 'model/SPEC.md', 'model/ccf_lt2/RESULTS.json',
                'model/ccf_lt2/EINSTEIN_BRIDGE.md',
                'model/particle_joint_reconstruction/dependencies.json',
                'Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean',
                'model/matched_einstein_limit/PROOF.md',
                'model/openwave_gravity_capacity/Transfer.lean']
    return {
        'schema': 'oph.ccf_lt2_claims.v1',
        'official_model': 'CCF-LT-2',
        'specification': 'model/SPEC.md',
        'runtime': 'model/ccf_lt2',
        'candidate': 'A-ccf2-local-capacity-v2',
        'law': 'H_v2 = sum_v [A_v + M_v (g_v + N_v (h_v-e_v I) + e_v I - F I) + F I]',
        'parent_polynomial_sha256': PARENT_DIGEST,
        'native_car_modes_per_vertex': 96,
        'edge': {'modes': 192, 'matter_monomials': 63589},
        'triangle': {'modes': 288, 'matter_monomials': 109675},
        'local_capacity_column_nonzero_outputs': [c['complete_nonzero_outputs'] for c in receipt['actual_local_capacity_join']['columns']],
        'reference': {'beta_exact': '1/1024', 'Gauss_reference_charge': 0, 'flux_cap': 6, 'geometry_dimension_per_vertex': 48},
        'phase_instrument_initial_derivatives': ['-9/32', '9/32', '0'],
        'stationarity_countervariation_exact': '9/16384',
        'verified_scope': [
            'Exact native parent regeneration, signed Dirac and current-square/family polynomial joins.',
            'Full-polynomial Hermiticity, Gauss identities, complete energy allocation and literal local-capacity native columns.',
            'Complete finite field, Dirac and local-capacity diagnostics, exact phase instrument, source sectors, faithful reference, coframe derivatives, shell and energy Ward checks.',
        ],
        'physical_boundary': {
            'arbitrary_full_native_state_propagation': False,
            'propagated_master_source_remote_detector_records': False,
            'independent_neutral_spectator_factorization': False,
            'relativistic_interacting_continuum': False,
            'Einstein_premises_instantiated': False,
            'physical_Newton_G_or_cosmological_Lambda_derived': False,
            'independent_physical_particle_masses_GeV': [],
            'official_physical_row_awards': [],
        },
        'calibration': 'Parent empirical inputs and binary coefficients are retained; contact, family, kinetic, field, geometry and capacity parameters are declared construction inputs.',
        'Einstein_route': {
            'theorem': 'OPH.EinsteinBranch.Quantitative.small_ball_tensor_bound',
            'decoder_norm': '17/8',
            'tensor_error': '255*(e_bulk+eta*e_area+e_stationarity)/(32*pi*eta*ell^4)',
            'required': 'Actual common-source ten-probe regional data; errors o(ell^4), eta -> eta_*>0, common tensor/metric and reference limit, Ward/Bianchi and a shared physical scale.',
            'Newton_transfer': 'The same tensor residual gives source error <=2 epsilon_E and exterior force error <=2 epsilon_E R^3/(3 r^2) under the static weak-field hypotheses.',
        },
        'evidence_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in evidence},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    payload = json.dumps(generate(), sort_keys=True, indent=2)+'\n'
    if args.write:
        OUTPUT.write_text(payload)
    elif OUTPUT.read_text() != payload:
        raise SystemExit('CCF-LT-2 claim register differs from verified current inputs')
    print('PASS_CCF_LT2_FINITE_CLAIMS_AND_PHYSICAL_BOUNDARY')


if __name__ == '__main__':
    main()
