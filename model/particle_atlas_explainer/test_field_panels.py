"""Reject physical promotions and omitted source sectors in the field panel."""
import copy
import re
from fractions import Fraction
import pytest

from . import field_panels
from .canvas import Canvas


@pytest.fixture(scope='module')
def data():
    return field_panels.load()


def test_source_receipts_and_domains(data):
    field_panels.validate(data)
    assert len(data['pins']) > 350


def test_printed_energy_and_fidelity_bounds_are_outward(data):
    canvas = Canvas()
    field_panels.charged_field_panel(canvas, data, 0)
    text = ''.join(canvas.items)
    match = re.search(r'Total charged energy: \[([\d.]+), ([\d.]+)\]', text)
    lo, hi = map(Fraction, match.groups())
    actual_lo, actual_hi = map(Fraction, data['charged_field']['actual_charged_ground_energy_interval'])
    assert lo <= actual_lo <= actual_hi <= hi
    shown = Fraction(re.search(r'fidelity above ([\d.]+)\.', text).group(1))
    assert shown <= Fraction(data['charged_clock']['finite_rational_gate_record_bounds']['accepted_source_ground_fidelity_lower'])


def test_native_field_panels_retain_their_distinct_domains(data):
    canvas=Canvas()
    end=field_panels.native_energy_panel(canvas,data,0)
    field_panels.dirac_panel(canvas,data,end+35)
    text=''.join(canvas.items)
    assert 'all 16 occupations' in text
    assert 'admitted-bin energy' in text
    assert 'no electron weak chirality' in text
    assert 'full spatial range' in text
    assert 'Same-H arrival then energy records' in text


def test_complete_current_and_field_panels(data):
    canvas=Canvas()
    end=field_panels.higgs_current_panel(canvas,data,0)
    field_panels.coupled_field_panel(canvas,data,end+35)
    text=''.join(canvas.items)
    assert 'all 80 apparatus labels' in text
    assert '1,161 outcomes' in text
    assert 'not the coupled ground' in text
    for key, label in [('actual_unconditional_target_probability_lower','Pair + field probability'),
                       ('actual_unconditional_binary_connected_response_lower','Connected matter–field response')]:
        shown=Fraction(re.search(re.escape(label)+r' ≥ ([\d.e-]+)',text).group(1))
        assert shown <= Fraction(data['coupled_field']['complete_interaction_records']['exact'][key])


def test_photon_cutoff_chart_keeps_outward_complete_record_bounds(data):
    canvas=Canvas();field_panels.photon_cutoff_panel(canvas,data,0)
    text=''.join(canvas.items)
    shown=[Fraction(x)for x in re.findall(r'error ≤ ([\d.e-]+)',text)]
    rows=[r for r in data['photon_cutoff']['cutoff_table']if r['cap']in(1,2,3,4,8)]
    assert len(shown)==len(rows)==5
    assert all(a>=Fraction(b['state_vector_error_upper'])for a,b in zip(shown,rows))
    assert '1,451 labels' in text and 'fixed graph' in text


def test_quark_energy_bins_and_original_weights_are_outward(data):
    canvas=Canvas();field_panels.quark_response_panel(canvas,data,0)
    text=''.join(canvas.items)
    shown=re.findall(r'\[([\d.]+), ([\d.]+)\]',text)
    weights=re.findall(r'original attempts ≥ ([\d.]+)%',text)
    assert len(shown)==len(weights)==3
    for key,bounds,weight in zip(('uc','dc','Q'),shown,weights):
        row=max(data['quark_response']['full_interacting_prepared_records'][key]['records'],
                key=lambda r:Fraction(r['original_attempt_probability_lower']))
        lo,hi=map(Fraction,row['source_excitation_relative_actual_ground'])
        assert Fraction(bounds[0])<=lo<=hi<=Fraction(bounds[1])
        assert Fraction(weight)<=100*Fraction(row['original_attempt_probability_lower'])
    assert '12 exact local zeros' in text and 'does not count particle species' in text


def test_complete_ground_and_local_control_domains(data):
    canvas=Canvas();end=field_panels.coupled_ground_panel(canvas,data,0)
    field_panels.local_current_panel(canvas,data,end+35)
    text=' '.join(re.sub('<[^>]+>',' ',''.join(canvas.items)).split())
    assert '96 active lowest states' in text
    assert 'not radiation from a vacuum' in text
    assert 'source Hamiltonian continues' in text
    assert 'Earlier state preparation' in text
    p=Fraction(re.search(r'success ≥ ([\d.]+)%',text).group(1))/100
    assert p<=Fraction(data['coupled_ground']['complete_repaired_ground_preparation']['complete_tree_success_lower'])


def test_ground_cap_and_quark_interference_display_bounds(data):
    ground=data['ground_refinement']['uniform_ground_bounds_exact']
    assert Fraction('3.094e-9')>=Fraction(ground['ideal_ground_two_field_record_error_at_K1024_upper'])
    assert Fraction('9.549e-9')>=max(Fraction(ground[k])for k in (
        'fixed_pointer_complete_prep_and_count_error_upper',
        'fixed_pointer_complete_prep_and_two_field_tests_error_upper'))
    records=data['quark_higgs']['full_interacting_interference_records']
    assert all(Fraction('0.096796971')<=Fraction(records[k]['original_attempt_joint_contrast_lower'])for k in ('uc','dc'))
    assert Fraction('0.096907802')<=Fraction(records['Q']['original_attempt_joint_contrast_lower'])
    states=data['quark_higgs']['normalized_admitted_states_and_histories']['rows']
    assert all(Fraction(r['actual_same_history_total_variation_upper'])<=Fraction('5.193e-6')for r in states.values())


@pytest.mark.parametrize('name,change', [
    ('interacting limit', lambda d: d['field']['theorems'].__setitem__('nontrivial_interacting_relativistic_limit_proved', True)),
    ('vector identification', lambda d: d['field']['theorems'].__setitem__('Lorentz_vector_or_Proca_identification_proved', True)),
    ('hidden empirical input', lambda d: d['field']['construction']['measured_inputs'].append('electron mass')),
    ('dropped screening', lambda d: d['charge']['full_source_charged_bins']['ec']['energy_filtered_localization'].__setitem__('origin_probability_lower_inferred_by_discarding_screened_complement', True)),
    ('cycle promoted to continuum', lambda d: d['collision']['interval_certified_five_block_collision'].__setitem__('three_dimensional_continuum_claim', True)),
    ('charged GeV promotion', lambda d: d['charged_field']['scope'].__setitem__('mass_GeV', 0.511)),
    ('missing charged complement', lambda d: d['charged_field'].__setitem__('additional_state_truncation', True)),
    ('paused clock', lambda d: d['charged_clock']['phase_apparatus'].__setitem__('common_source_evolution_retained', False)),
    ('exact projection claim', lambda d: d['charged_clock']['scope'].__setitem__('exact_ground_projection', True)),
    ('named vector boson', lambda d: d['vector']['construction'].__setitem__('Standard_Model_particle_assignment', 'Z')),
    ('silent reference conjugation', lambda d: d['conjugation']['reference_boundary'].__setitem__('same_reference_total_energy_equality_claimed', True)),
    ('prepared ground substituted', lambda d: d['observed_charge'].__setitem__('exact_ground_projection_or_stationary_post_clock_state_claimed', True)),
    ('native energy promoted', lambda d: d['native_energy'].__setitem__('physical_spin_Dirac_or_mass_identification', True)),
    ('Dirac interaction transferred', lambda d: d['dirac']['field_comparison'].__setitem__('interacting_continuum_or_electromagnetic_Gauss_law_claimed', True)),
    ('Dirac species assigned', lambda d: d['dirac']['field_comparison'].__setitem__('SM_weak_chirality_or_observed_mass_identity_claimed', True)),
    ('Higgs selectors identified', lambda d: d['higgs_join']['same_H_cross_residue_and_records']['positive_cross_selector'].__setitem__('maximum_elementary_E_star_common_composite_overlap_proved', True)),
    ('nested window relabelled', lambda d: d['resolution'].__setitem__('original_maximum_residue_eigenvalue_relocalized', True)),
    ('centroid ground projection inserted', lambda d: d['centroid']['actual_physical_preparation'].__setitem__('extra_P0_ground_projection_used', True)),
    ('unknown ground phase inserted', lambda d: d['centroid']['clock'].__setitem__('exact_native_Eg_used_as_free_pointer_phase', True)),
    ('filtered preparation failures dropped', lambda d: d['higgs_state']['complete_instruments'].__setitem__('every_fine_preparation_record_and_partner_retained', False)),
    ('filtered state promoted to SM', lambda d: d['higgs_state']['typed_free_Dirac_corollary'].__setitem__('interacting_native_Lorentz_or_SM_chirality_transfer', True)),
    ('weak response promoted to W particle', lambda d: d['weak_response']['scope'].__setitem__('W_boson_or_Lorentz_vector_identified', True)),
    ('weak current loses original weight', lambda d: d['weak_response']['finite_inelastic_experiment'].__setitem__('original_attempt_nonzero_phase_probability_lower', '0')),
    ('current conjugation erased', lambda d: d['higgs_current']['conjugation'].__setitem__('W_invariant', True)),
    ('current labels discarded', lambda d: d['higgs_current']['apparatus'].__setitem__('complete_setting_pump_energy_current_alphabet_size', 5)),
    ('finite interaction promoted', lambda d: d['coupled_field']['boundaries'].__setitem__('interacting_Lorentz_or_four_dimensional_QED_limit_proved', True)),
    ('coupled vacuum substituted', lambda d: d['coupled_field']['complete_interaction_records'].__setitem__('source_state_is_actual_coupled_ground', True)),
    ('Gauss failure removed', lambda d: d['coupled_field']['count_instrument'].__setitem__('full_pointer_dimension', 1160)),
    ('photon cutoff promoted to spatial continuum', lambda d: d['photon_cutoff']['boundaries'].__setitem__('growing_volume_or_UV_limit', True)),
    ('photon overflow omitted', lambda d: d['photon_cutoff']['cap3_pair_photon_transfer'].__setitem__('common_overflow_padded_labels', 1161)),
    ('cutoff changes apparatus', lambda d: d['photon_cutoff']['cap3_pair_photon_transfer'].__setitem__('padding_changes_parent_executed_apparatus', True)),
    ('quark pole assigned', lambda d: d['quark_response']['scope'].__setitem__('rest_mass_or_GeV_calibration', True)),
    ('local equality promoted', lambda d: d['quark_response']['scope'].__setitem__('computed_degeneracy_is_physical_up_down_mass_equality', True)),
    ('ground complement discarded', lambda d: d['coupled_ground']['unchanged_source'].__setitem__('ground_in_code', True)),
    ('ground energy reset', lambda d: d['coupled_ground']['declared_onsite_completion'].__setitem__('global_scalar_energy_subtraction', 'fitted offset')),
    ('ground projector assumed', lambda d: d['coupled_ground']['complete_repaired_ground_preparation'].__setitem__('ground_energy_or_exact_ground_projector_used', True)),
    ('local detector stops source', lambda d: d['local_current']['source'].__setitem__('source_hamiltonian_continues_through_final_pulse', False)),
    ('local detector closes all controls', lambda d: d['local_current']['control_resources'].__setitem__('full_matrix_preparation_and_exact_spectral_bin_premises_removed', True)),
    ('local current promoted', lambda d: d['local_current']['scope'].__setitem__('Standard_Model_species_or_GeV_mass_identified', True)),
    ('quark composite made canonical', lambda d: d['quark_higgs']['scope'].__setitem__('canonical_composite_CAR', True)),
    ('quark reference omitted', lambda d: d['quark_higgs']['full_static_reference_controls'].__setitem__('full_joint_CAR_modes', 30)),
    ('quark mass assigned', lambda d: d['quark_higgs']['scope'].__setitem__('source_energies_are_calibrated_rest_masses', True)),
    ('quark setting labels silently equated', lambda d: d['quark_higgs']['normalized_admitted_states_and_histories']['rows']['uc'].__setitem__('same_history_policy_may_read_pump_setting', True)),
    ('ground heat clock substituted', lambda d: d['ground_refinement']['fixed_complete_energy_clock'].__setitem__('imaginary_heat_time_is_physical_clock', True)),
    ('ground complete filter alphabet lost', lambda d: d['ground_refinement']['fixed_complete_energy_clock'].__setitem__('same_pointer_alphabet_for_all_K_and_mathematical_comparison', False)),
    ('ground photon cap promoted to space', lambda d: d['ground_refinement']['same_source_family'].__setitem__('spatial_volume_or_spacing_refinement', True)),
])
def test_invalid_interpretation_is_rejected(data, name, change):
    altered = copy.deepcopy(data)
    change(altered)
    with pytest.raises(ValueError):
        field_panels.validate(altered)
