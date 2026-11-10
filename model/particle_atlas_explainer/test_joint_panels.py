"""The joined interaction graphic uses exact receipts without promoting scope."""
from copy import deepcopy
import xml.etree.ElementTree as ET

import pytest

from . import joint_panels
from .canvas import Canvas


@pytest.fixture(scope='module')
def joint():
    return joint_panels.load()


def visible(value):
    c = Canvas()
    joint_panels.draw(c, {'joint': value}, 0)
    svg = '<svg xmlns="http://www.w3.org/2000/svg">' + ''.join(c.items) + '</svg>'
    return ' '.join(' '.join(ET.fromstring(svg).itertext()).split()), c


def test_both_reciprocal_readouts_and_preparation_domains_are_visible(joint):
    text, c = visible(joint)
    assert not c.overflow
    for phrase in ('72,000 basis states', 'Click probability: 1/120 versus 0',
                   'Face detector: 1 versus 0', 'Return: 1 unread; 1/2 after readout',
                   '281/1080', '36,000 basis states', '281/2160',
                   'the retained link frame to compare amplitudes',
                   'Every failure is retained.',
                   'relativistic masses, generations and GeV values are unproved.'):
        assert phrase in text


def test_frequency_and_species_boundaries_are_visible(joint):
    text, c = visible(joint)
    assert not c.overflow
    for phrase in ('Fresh preparation success: 1/10', '32 records; 256 on the larger cube',
                   '1 : 2.17082 : 2.618034', 'source energies, not lepton masses',
                   'Maximum double occupation: π²/(π² + 16)',
                   'Five internal face locations describe one protected type',
                   'ω² = |p|² + 4; |ω(0)| = 2 (declared units)',
                   'No GeV or Standard Model identity follows.',
                   'Its spinor memory and coupling law are declared architecture.',
                   'Fermionic hopping and tensor-wave exchange are distinct spatial laws'):
        assert phrase.lower() in text.lower()


def test_all_construction_receipts_are_loaded_and_pinned(joint):
    assert set(joint['extensions'])==set(joint_panels.EXTENSIONS)
    pins=joint_panels.input_pins(joint)
    for name in joint_panels.EXTENSIONS:
        assert 'model/coherent_defect_wave/'+name+'/results.json' in pins


def test_rest_parameter_and_collision_number_follow_receipts(joint):
    changed=deepcopy(joint)
    changed['extensions']['mobile_flux']['protected_object_route']['positive_branch_rest_parameter']=3
    changed['extensions']['energy']['finite_dark_collision']['maximum_double_occupation_probability']='1/7'
    text,_=visible(changed)
    assert 'ω² = |p|² + 9; |ω(0)| = 3 (declared units)' in text
    assert 'Maximum double occupation: 1/7' in text


@pytest.mark.parametrize('claim', ['SM_identity_or_GeV_mass_derived','GeV_prediction',
                                 'energies_are_physical_rest_masses', 'interacting_Lorentz_limit',
                                 'GeV_mass_derived','chiral_SM_particle_identification',
                                 'interacting_relativistic_limit_or_SM_or_GeV'])
def test_nested_extension_promotion_is_rejected(claim):
    with pytest.raises(ValueError,match='physical promotion'):
        joint_panels.check_extension_scope({'nested':{'proof':[{'scope':{claim:True}}]}},'test')


def test_shared_register_and_complete_join_bounds_are_required(joint):
    changed=deepcopy(joint['extensions']['native_routing_join'])
    changed['shared_tuple']['native_U_reads_the_same_link_register_as_routing']=False
    with pytest.raises(ValueError,match='shared link'):
        joint_panels.check_native_join(changed)
    changed=deepcopy(joint['extensions']['native_routing_join'])
    changed['finite_join']['Taylor_output_renormalized']=True
    with pytest.raises(ValueError,match='approximate readout'):
        joint_panels.check_native_join(changed)
    changed=deepcopy(joint['extensions']['native_routing_join'])
    changed['finite_join']['exact_unitary_probability_enclosures'].pop()
    with pytest.raises(ValueError,match='detector outcome'):
        joint_panels.check_native_join(changed)


def test_numeric_preparation_label_is_receipt_driven(joint):
    changed = deepcopy(joint)
    changed['reference']['success_probability_from_uniform_basis_ensemble'] = '2/7'
    assert 'probability 2/7' in visible(changed)[0]


def test_unconstrained_return_cannot_be_advertised_as_invariant(monkeypatch, joint):
    changed = deepcopy(joint)
    changed['interaction']['constrained_coherence']['return_projector_commutes_with_Gamma'] = False
    monkeypatch.setattr(joint_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='admitted gauge sector'):
        joint_panels.load()


@pytest.mark.parametrize('claim', ['physical_particle_species', 'physical_energy_or_GeV_mass',
                                 'interacting_3plus1_continuum'])
def test_physical_promotion_rejected(monkeypatch, joint, claim):
    changed = deepcopy(joint)
    changed['scope'][claim] = True
    monkeypatch.setattr(joint_panels.json, 'loads', lambda _: changed)
    with pytest.raises(ValueError, match='scope promotion'):
        joint_panels.load()


def test_positive_free_mass_current_binding_and_complete_apparatus_are_visible(joint):
    text,c=visible(joint)
    assert not c.overflow
    for phrase in ('20 graded flux slots per cell', 'Positive particles and holes',
                   'native gauge currents act trivially on these bare flux objects',
                   'Free rest energy 2; neutral-pair threshold 4',
                   '180 actual current identities; native K = 0',
                   'Reservoir + battery', '36 outcomes, including failures',
                   'E² = |p|² + 4; rest energy = 2 source units in this free limit',
                   'A chiral field dictionary from the actual source currents',
                   'No electron/muon/tau generation is selected',
                   '3,424 neutral physical states and one certified ground',
                   '0.97009613 source units',
                   'An actual Higgs composite excites this same line',
                   'Weak currents act only on left components',
                   'one flux and one native occupation, bounded coupling gives vanishing contact response',
                   'fixed-volume, fixed native-momentum window'):
        assert phrase in text


def test_periodic_story_keeps_preparation_failure_and_composite_scope(joint):
    periodic = deepcopy(joint['extensions']['chiral_joint/periodic_spectrum'])
    higgs = deepcopy(joint['extensions']['chiral_joint/higgs_composites'])
    periodic['compensated_charged_spectroscopy']['original_ground_preparation_failure'] = '0'
    with pytest.raises(ValueError, match='complete reference records'):
        joint_panels.check_periodic(periodic, higgs)
    periodic = deepcopy(joint['extensions']['chiral_joint/periodic_spectrum'])
    higgs['same_periodic_H_spectral_bridge']['residues_are_instrument_probabilities'] = True
    with pytest.raises(ValueError, match='residue interpretation'):
        joint_panels.check_periodic(periodic, higgs)


def test_local_reference_retains_dark_sector_and_charge_boundaries(joint):
    original = joint['extensions']['chiral_joint/neutral_local_spectrum']
    changed = deepcopy(original)
    changed['complete_physical_spectrum']['full_component_excitation_gap_inferred_from_visible_line_spacing'] = True
    with pytest.raises(ValueError, match='reference scope'):
        joint_panels.check_local_reference(changed)
    changed = deepcopy(original)
    changed['physical_symmetries']['raw_L_is_identified_with_measured_lepton_charge'] = True
    with pytest.raises(ValueError, match='charge or spin'):
        joint_panels.check_local_reference(changed)
    text, c = visible(joint)
    assert not c.overflow
    assert 'empty-register cyclic ground cannot be the universal vacuum' in text
    assert 'not a negative particle mass' in text


def test_neutral_detectors_preserve_failures_and_shared_source(joint):
    neutral = deepcopy(joint['extensions']['chiral_joint/neutral_composite'])
    bosons = deepcopy(joint['extensions']['chiral_joint/periodic_boson_records'])
    neutral['complete_readout']['success_postselection_renormalization'] = True
    with pytest.raises(ValueError, match='original response weights'):
        joint_panels.check_neutral_detectors(neutral, bosons)
    neutral = deepcopy(joint['extensions']['chiral_joint/neutral_composite'])
    bosons['complete_readout']['settings'][0]['lowest_certified_positive_inelastic_record']['energy_transfer_interval'] = ['2', '3']
    with pytest.raises(ValueError, match='shared response'):
        joint_panels.check_neutral_detectors(neutral, bosons)


def test_apparatus_cannot_drop_failures_or_switch_the_clock(joint):
    free=joint['extensions']['many_body_route']
    original=joint['extensions']['mass_instrument']
    changed=deepcopy(original);changed['experiments'][0]['complete_records'].pop()
    with pytest.raises(ValueError,match='complete reservoir'):
        joint_panels.check_mass_instrument(changed,free)
    changed_free=deepcopy(free);changed_free['reference']['unique_positive_excitation_ground']=False
    with pytest.raises(ValueError,match='positive free reference'):
        joint_panels.check_mass_instrument(original,changed_free)
    changed=deepcopy(original);changed['source']='different clock'
    with pytest.raises(ValueError,match='different source clocks'):
        joint_panels.check_mass_instrument(changed,free)
    changed=deepcopy(original);changed['unit_map']['source_energy_unit_in_GeV']='1'
    with pytest.raises(ValueError,match='unowned GeV scale'):
        joint_panels.check_mass_instrument(changed,free)


def test_dark_current_attachment_cannot_be_presented_as_native_binding(joint):
    changed=deepcopy(joint['extensions']['field_identity_audit'])
    changed['native_dark_current_module']['native_K_on_module']='1'
    with pytest.raises(ValueError,match='native current module'):
        joint_panels.check_identity(changed)
    changed=deepcopy(joint['extensions']['field_identity_audit'])
    changed['mass_normalization_and_chirality']['Dirac_R15_contains_both_chiralities_in_same_internal_representation']=False
    with pytest.raises(ValueError,match='chirality boundary'):
        joint_panels.check_identity(changed)


def test_conditional_dictionary_cannot_be_promoted_to_measured_species_or_masses(joint):
    value=deepcopy(joint['extensions']['chiral_particle_map'])
    value['conditional_field_rows'][0]['empirically_identified_species']='top quark'
    with pytest.raises(ValueError,match='physical promotion'):
        joint_panels.check_extension_scope(value,'chiral_particle_map')
    value=deepcopy(joint['extensions']['chiral_particle_map'])
    value['identification']['measured_GeV_predictions']={'electron':'0.000511'}
    with pytest.raises(ValueError,match='physical promotion'):
        joint_panels.check_extension_scope(value,'chiral_particle_map')


def test_two_cell_energy_story_preserves_full_ground_and_preparation_failures(joint):
    value=deepcopy(joint['extensions']['chiral_joint/two_cell'])
    value['original_attempt_readout']['records'].pop(0)
    with pytest.raises(ValueError,match='preparation failures'):
        joint_panels.check_two_cell(value)
    value=deepcopy(joint['extensions']['chiral_joint/two_cell'])
    value['full_spectra'][0]['ground_is_certifiably_simple']=False
    with pytest.raises(ValueError,match='physical ground'):
        joint_panels.check_two_cell(value)
    value=deepcopy(joint['extensions']['chiral_joint/two_cell'])
    value['original_attempt_readout']['records'][1]['probability_interval']=['-1','2']
    with pytest.raises(ValueError,match='invalid probability interval'):
        joint_panels.check_two_cell(value)


def test_full_quantum_transport_and_local_ground_cannot_drop_their_domains(joint):
    import copy
    cases=[('chiral_joint/peter_weyl_links', ('complete_neutral_transfer_experiment','extra_neutral_hopping_term'), True),
           ('chiral_joint/family_transport', ('same_clock_full_H_probability_records','source_or_reference_states_discarded'), True),
           ('chiral_joint/local_reference', ('global_physical_ground_and_gap','dark_energy_one_rank'), 0),
           ('chiral_joint/local_reference', ('global_physical_ground_and_gap','interacting_many_cell_vacuum_proved'), True)]
    for package,path,value in cases:
        altered=copy.deepcopy(joint['extensions'])
        altered[package][path[0]][path[1]]=value
        with pytest.raises(ValueError):
            joint_panels.check_PW_transport(altered)


def test_interacting_line_cannot_drop_complement_or_change_mass_domain(joint):
    import copy
    cases=[('chiral_joint/atomic_spin_geometry','full_interacting_transfer',
            'unit_hopping_ground_statement_transferred',True),
           ('chiral_joint/atomic_spin_geometry','full_interacting_transfer',
            'lowest_gap_is_a_zero_momentum_particle_rest_mass',True),
           ('chiral_joint/atomic_spin_geometry','full_interacting_transfer',
            'all_full_Fock_scalar_link_and_physical_Gauss_sectors_retained',False)]
    for package,section,key,value in cases:
        changed=copy.deepcopy(joint['extensions'])
        changed[package][section][key]=value
        with pytest.raises(ValueError):
            joint_panels.check_interacting_reference(changed)
    changed=copy.deepcopy(joint['extensions'])
    line=changed['chiral_joint/collective_spin_response']['same_source_interacting_records']['actual_interacting_energy_line']
    line['outside_line_probability_conditioned_on_ground_preparation']=['0','0']
    with pytest.raises(ValueError):
        joint_panels.check_interacting_reference(changed)


def test_charged_panel_rejects_free_reference_and_unproved_pole(joint):
    import copy
    for section,key,value in (
        ('source','reference_particle_cost',0),
        ('scope','physical_rest_mass_or_free_charge_dispersion',True),
        ('full_interacting_energy_readout','isolated_charged_eigenvalue_or_unchanged_multiplicity_claimed',True),
        ('constructive_approximate_preparation_records','all_preparation_success_and_failure_branches_obey_source_Gauss',False),
    ):
        changed=copy.deepcopy(joint['extensions'])
        changed['chiral_joint/interacting_charged_response']['same_H_interacting_channels']['ec'][section][key]=value
        with pytest.raises(ValueError):
            joint_panels.check_charged_reference(changed)
