"""Conditional chiral extension and inherited calibrated spectrum controls.

The selected primary source is bound as an architectural ancestor. The
extension has no supplied source embedding, vacuum or interacting limit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as s

from model.common_source import calibrated_matter
from model.conditional_standard_model import contract as sm_contract
from model.conditional_standard_model import matching, packet, spectrum
from model.particle_candidate_atlas import conditional_standard_model as atlas
from model.particle_common_contract import contract as primary, provenance
from model.particle_yang_mills_gap_contract import transfer as gap
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('candidate_results.json')
SCHEMA = 'oph.joint_reconstruction.candidate.v1'
EXTENSION_ID = 'oph.conditional_chiral_gauge_higgs_extension.v1'
DIRECT_INPUTS = (
    'model/particle_joint_reconstruction/candidate.py',
    'model/particle_joint_reconstruction/test_candidate.py',
    'model/particle_joint_reconstruction/CANDIDATE.md',
    'model/conditional_standard_model/PROOF.md',
    'model/common_source/calibrated_matter_parameters.json',
    'model/common_source/calibrated_neutrino_parameters.json',
    'model/particle_yang_mills_gap_contract/PROOF.md',
)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n'


def assumption_ledger():
    """Specific construction premises; none is discharged on the primary."""
    rows = (
        ('extension_carrier_and_action',
         ['full_fock_spin_chiral_multiplicity', 'effective_action_and_normalization', 'physical_global_group'],
         'Declare three copies of the all-left M15 representation, two-component Weyl CAR, '
         'a complex Higgs doublet and canonically normalized gauge and family kinetic terms. '
         'Choose the global gauge group and its charge lattice. Supply a regulated extension '
         'and an embedding or controlled record map from the primary; preserve all original '
         'source/state/reference/clock data on the claimed common domain.',
         ['model/conditional_standard_model/PROOF.md', 'model/common_source/chiral_candidate/action.py']),
        ('chiral_measure_and_mirror_control',
         ['same_source_continuum_matching', 'full_fock_spin_chiral_multiplicity'],
         'Construct a gauge-invariant chiral measure, reflection-positive physical correlations '
         'and controlled removal of mirror species. Local anomaly sums alone supply neither '
         'this measure nor a regulator limit.',
         ['model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md']),
        ('common_vacuum_reference_clock',
         ['vacuum_and_phase', 'parameter_and_unit_provenance'],
         'Supply compatible vacuum-limit preparation and reference maps for the extension, '
         'with a single source-energy-to-GeV conversion and shared physical time on every '
         'spectral and record claim. Control the ordered Euclidean correlator error at each '
         'fixed physical time; low energy excess alone is insufficient.',
         ['model/particle_common_contract/PROOF.md', 'model/particle_yang_mills_gap_contract/PROOF.md']),
        ('renormalized_four_dimensional_correlations',
         ['same_source_continuum_matching', 'effective_action_and_normalization'],
         'On a common local test domain prove joint infinite-volume, cutoff and spacing limits '
         'of renormalized physical correlations with positivity, Euclidean covariance, '
         'regularity, locality and clustering sufficient for relativistic reconstruction '
         'on R4. Identify the limiting currents and action coefficients with this extension.',
         ['model/conditional_standard_model/PROOF.md', 'extra/yang_mills_gap_clay_problem.tex']),
        ('nonzero_connected_interaction',
         ['same_source_continuum_matching', 'interacting_poles_and_observables'],
         'Exhibit a nonvanishing renormalized interaction witness in a fixed generating-field '
         'and current dictionary that is incompatible with a free/Wick theory under '
         'the admitted identifications, or an appropriate nonfree scattering criterion. '
         'Connected correlations of nonlinear composites in a free theory do not suffice; '
         'nonzero bare couplings alone do not prove this premise.',
         ['model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md']),
        ('physical_current_spin_detector_dictionary',
         ['defect_sector_intertwiner', 'interacting_poles_and_observables'],
         'Identify conserved charge, Lorentz spin and the physical spectral or inclusive '
         'detector channels on the same Hilbert domain. A claimed source-sector map must '
         'control the full rectangular generator residual, including image leakage; '
         'spectrum exhaustion is a separate premise.',
         ['model/conditional_standard_model/matching.py', 'model/conditional_standard_model/PROOF.md']),
        ('radiative_infrared_resonance_matching',
         ['parameter_and_unit_provenance', 'interacting_poles_and_observables'],
         'Supply scheme- and scale-consistent radiative matching, uncertainties, inclusive '
         'infrared prescriptions and unstable-resonance widths. Relate quark running '
         'parameters to confined observables. Tree coefficients at MZ are not measured poles.',
         ['model/common_source/calibrated_matter_parameters.json', 'model/conditional_standard_model/PROOF.md']),
        ('neutrino_EFT_completion_and_running',
         ['effective_action_and_normalization', 'parameter_and_unit_provenance'],
         'For the separate calibrated Weinberg arm, declare a dimension-five effective '
         'theory window, transport its low-energy parameters to the common scale with '
         'controlled errors, and supply a UV completion for an all-energy claim. Normal '
         'ordering, the lightest mass and Majorana phases retain their declared ancestry.',
         ['model/common_source/calibrated_neutrino_parameters.json', 'model/common_source/calibrated_matter.py']),
        ('pure_YM_source_and_sector_attachment',
         ['same_source_continuum_matching', 'physical_global_group'],
         'Declare a separate pure compact-simple gauge source and prove its four-dimensional '
         'Yang-Mills identification and nontrivial local curvature observables. Any transfer '
         'back to the matter-containing product-group theory requires an invariant sector '
         'or decoupling theorem, including a clock map.',
         ['model/particle_yang_mills_gap_contract/PROOF.md']),
        ('pure_YM_dense_decay_and_noncollapse',
         ['vacuum_and_phase', 'interacting_poles_and_observables'],
         'For that pure source supply a reconstructed positive Hamiltonian, dense centered '
         'local vacuum vectors and common-clock connected Euclidean decay C_A(kh)<=B_A q^k '
         'with 0<q<1. Supply a surviving vector with positive variance and finite energy '
         'to bound the nonzero spectral threshold above. These are separate from finite gaps.',
         ['model/particle_yang_mills_gap_contract/PROOF.md']),
    )
    return [{'id': key, 'conditional_standard_model_interfaces': interfaces,
             'required_statement': statement, 'references': refs,
             'assumed_in_conditional_scenario': True,
             'proved_on_primary': False, 'witness_supplied_by_candidate': False}
            for key, interfaces, statement, refs in rows]


def _bindings(theorem):
    inherited = packet.source_bindings()
    paths = set(DIRECT_INPUTS) | set(inherited) | set(atlas.SOURCE_INPUTS) | set(gap.SOURCE_INPUTS)
    for row in assumption_ledger():
        paths.update(row['references'])
    for section in ('theorem_chain', 'conditional_interfaces'):
        for row in theorem[section]:
            paths.update(ref['path'] for ref in row['references'])
    return provenance.definition_hashes(provenance.definition_inputs(paths))


def _matrix(record):
    return s.ImmutableSparseMatrix(*record['shape'],
                                  {(i, j): s.sympify(v) for i, j, v in record['entries']})


def quadratic_consistency(catalogue, symbolic):
    """Cross-check two existing action producers at the same supplied inputs."""
    inputs = spectrum.calibrated_inputs()
    q = atlas.quadratic_action()
    substitutions = dict(zip(q['symbols'],
                             (inputs.gW, inputs.gY, inputs.vEW, inputs.lambda_SM), strict=True))
    gauge = q['gauge_mass_squared'].subs(substitutions)
    neutral = _matrix(catalogue['bosons']['neutral_real_vector_mass_squared_matrix'])
    if gauge[2:4, 2:4] != neutral:
        raise ArithmeticError('the inherited neutral gauge Hessians disagree')
    bosons = {row['name']: row for row in catalogue['bosons']['rows']}
    for name, key in (('W+', 'm_W_squared'), ('Z', 'm_Z_squared'), ('Higgs_radial', 'm_H_squared')):
        if s.simplify(q[key].subs(substitutions)-s.sympify(bosons[name]['tree_mass_squared_exact'])) != 0:
            raise ArithmeticError('the inherited boson action coefficients disagree')
    count = sum(row['multiplicity']*row['free_polarizations_per_species'] for row in bosons.values())
    census = next(row for row in symbolic['field_censuses'] if row['neutrino_completion'] == 'majorana_eft')
    if count != census['boson_on_shell_degrees_broken'] or count != 28:
        raise ArithmeticError('quadratic boson census changed')
    mass2 = s.sympify(bosons['photon']['tree_mass_squared_exact'])
    if mass2 != 0:
        raise ArithmeticError('unbroken photon null direction changed')
    soft = [{'momentum_magnitude_GeV': str(s.Rational(1, n)),
             'quadratic_energy_GeV': str(s.sqrt(s.Rational(1, n)**2+mass2))}
            for n in (1, 10, 100, 1000)]
    return {'independent_inherited_gauge_Hessians_agree': True,
            'independent_inherited_boson_mass_coefficients_agree': True,
            'boson_on_shell_degrees': count, 'massless_photon_exact': True,
            'soft_photon_sequence_control': soft,
            'arbitrarily_soft_photon_follows_under_continuum_dispersion': True,
            'positive_global_gap_of_photon_containing_continuum_scenario': False,
            'this_is_a_gap_statement_about_primary': False}


def build_candidate():
    """Recompute inherited action controls and bind every calibration ancestor."""
    theorem = sm_contract.theorem_contract()
    bindings = _bindings(theorem)
    member = primary.member(0)
    primary.verify_member(member)
    catalogue = spectrum.build_catalogue()
    symbolic = atlas.certificate()
    neutrinos = calibrated_matter.neutrino_report()
    ledger = assumption_ledger()
    known_interfaces = {row['id'] for row in theorem['conditional_interfaces']}
    if any(set(row['conditional_standard_model_interfaces'])-known_interfaces for row in ledger):
        raise ValueError('candidate assumption has an unknown inherited interface')
    checks = quadratic_consistency(catalogue, symbolic)
    if (member['source']['CAR_modes_per_cell'] != 60
            or catalogue['counts']['static_spin_times_internal_CAR_annihilation_modes'] != 90):
        raise ArithmeticError('source/extension carrier comparison changed')
    bosons = {row['name']: row for row in catalogue['bosons']['rows']}
    result = {
        'schema': SCHEMA, 'candidate_id': EXTENSION_ID,
        'pass_scope': 'conditional specification, exact finite action consistency, calibration custody and replay',
        'finite_sanity_checks_pass': True,
        'primary_member': member,
        'extension_boundary': {
            'primary_CAR_modes_per_cell': 60, 'extension_static_spin_internal_CAR_modes': 90,
            'primary_representation': member['source']['gauge_representation'],
            'extension_representation': 'three all-left M15 families times two-component Weyl spin, plus Higgs',
            'same_carrier_or_basis_change': False,
            'common_native_internal_current_tensors': True,
            'shared_internal_tensors_imply_same_spatial_action': False,
            'source_embedding_supplied': False, 'state_reference_clock_transfer_supplied': False,
            'primary_GeV_conversion': member['clock']['GeV_or_proper_time_conversion'],
            'extension_energy_unit': catalogue['inputs']['energy_unit'],
            'extension_is_executed_primary_program': False,
            'extension_finite_preparation_and_locality_inherited_automatically': False,
        },
        'assumptions': ledger,
        'inherited_conditional_theorem_contract': theorem,
        'inherited_symbolic_chiral_action': symbolic,
        'inherited_running_coefficient_arm': {
            'catalogue': catalogue,
            'role': 'empirical MSbar common-scale tree-action consistency, not measured-pole recovery',
            'minimal_neutrinos_fit_nonzero_oscillation_splittings': False,
            'tree_boson_coefficients_are_measured_poles': False,
            'boson_tree_coefficient_GeV_display': {name: bosons[name]['tree_mass_decimal_display']
                                                  for name in ('W+', 'Z', 'Higgs_radial')},
            'direct_full_physical_spectrum_assignment_pass': False,
        },
        'conditional_physical_completion_arm': {
            'charged_and_boson_input_source': 'inherited_running_coefficient_arm.catalogue',
            'neutrino_choice': 'separate calibrated dimension-five Weinberg EFT',
            'neutrino_producer': neutrinos,
            'minimal_catalogue_neutrino_entries_overwritten': False,
            'neutrino_mass_eigenstates_identified_with_flavor_labels': False,
            'oscillation_input_recovery_is_independent_prediction': False,
            'all_reference_targets_numerically_matched': False,
            'compatibility_with_full_reference_data': 'conditional feasibility only; no complete matching witness or error budget',
            'all_energy_UV_completion_supplied': False,
            'gauge_theta_angles_vacuum_energy_and_global_charge_lattice_selected_by_primary': False,
            'physical_poles_inclusive_rates_and_resonance_widths_supplied': False,
        },
        'quadratic_cross_checks': checks,
        'full_rectangular_matching_control': matching.example_certificate(),
        'pure_Yang_Mills_branch': {
            'criterion': 'model/particle_yang_mills_gap_contract/PROOF.md',
            'finite_transfer_controls': gap.certificate(),
            'target': 'separate pure compact-simple gauge theory on R4',
            'physical_lower_bound': '-log(q)/h under the dense-vacuum-correlation hypotheses',
            'finite_upper_bound': 'M/v0 under a surviving variance v0 and energy bound M',
            'required_assumptions': ['pure_YM_source_and_sector_attachment', 'pure_YM_dense_decay_and_noncollapse'],
            'same_as_matter_containing_product_group_primary': False,
            'full_photon_containing_theory_has_positive_global_gap': False,
            'strong_sector_gap_transferred_to_full_theory': False,
            'invariant_sector_or_decoupling_proved': False,
            'Clay_problem_solved': False,
        },
        'scope': {
            'interacting_four_dimensional_QFT_constructed': False,
            'conditional_interacting_reconstruction_statement': True,
            'full_elementary_quadratic_comparison_channels': True,
            'complete_hadron_resonance_or_physical_spectrum': False,
            'source_selected_physical_masses': False,
            'independent_predicted_reference_mass_count': 0,
            'canonical_physical_premises_discharged': [],
            'integrated_gate_accepted': False,
            'fresh_full_spatial_preparation_executed': False,
        },
        'source_sha256': bindings,
    }
    if _bindings(theorem) != bindings:
        raise ValueError('candidate inputs changed during calculation')
    return result


def validate_candidate(value):
    if canonical(value) != canonical(build_candidate()):
        raise ValueError('conditional candidate differs from the fresh inherited calculations')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(build_candidate()), encoding='utf-8')
    else:
        validate_candidate(strict_json.load(OUTPUT))
    print('PASS: conditional candidate and inherited action controls; no physical continuum discharge')


if __name__ == '__main__':
    main()
