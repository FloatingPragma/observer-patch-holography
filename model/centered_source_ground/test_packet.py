from copy import deepcopy
import pytest
from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet._produce()


def test_cold_replay_matches_direct_source(receipt):
    assert packet.report()==receipt


def test_centered_reference_cannot_be_reported_as_the_old_one(receipt):
    value=receipt;packet.validate(value)
    for key in ('old_ground_reused_without_check','interacting_Lorentz_continuum_proved','physical_mass_or_species_identified'):
        wrong=deepcopy(value);wrong['bounds'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)
    wrong=deepcopy(value);wrong['cube']['bounds']['neutral_gap_from_charge_bounds']=True
    with pytest.raises(ValueError):packet.validate(wrong)
    for key, changed in (('full_centered_gap_strict_lower', '1/64'),
                         ('balanced_N8_dimension', 588),
                         ('thermodynamic_gap_proved', True),
                         ('relativistic_mass_identified', True)):
        wrong=deepcopy(value);wrong['finite_gap'][key]=changed
        with pytest.raises(ValueError):packet.validate(wrong)
    for key,changed in (('ground_multiplicity',5),('bipartition_sizes',[3,5]),
                        ('proof_imported_not_reproved',False)):
        wrong=deepcopy(value);wrong['cube']['imported_Lieb_theorem'][key]=changed
        with pytest.raises(ValueError):packet.validate(wrong)
    for key in ('infinite_ground_local_faithfulness_inferred','geometric_modular_charge_identified'):
        wrong=deepcopy(value);wrong['cube_local_state'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)


@pytest.mark.parametrize('section,key,value', [
    ('selection','joint_singlet_isotypic_dimension',1),
    ('selection','internal_spin_is_physical_Lorentz_spin',True),
    ('selection','every_isotypic_copy_reachable_claimed',True),
    ('complete_repair','native_raw_slots_retained',119),
    ('complete_repair','every_outside_map_identity_on_entire_domain',False),
    ('complete_repair','nontrivial_native_defect_record_obtained',True),
    ('complete_repair','new_native_preparation_or_detector_executed',True),
    ('response_and_current','physical_electromagnetic_normalization_or_CPT_claimed',True),
    ('response_and_current','oriented_bond_currents_checked',1),
    ('onsite_mass_order','mass_order_fluctuations_claimed_zero',True),
    ('onsite_mass_order','interacting_gap_or_thermodynamic_symmetry_breaking_excluded',True),
])
def test_charged_response_cannot_be_promoted_or_truncated(receipt,section,key,value):
    wrong=deepcopy(receipt);wrong['charged_particles'][section][key]=value
    with pytest.raises(ValueError):packet.validate(wrong)


def test_empty_custody_cannot_authenticate_a_centered_reference(receipt):
    wrong=deepcopy(receipt);wrong['source_bindings']={}
    with pytest.raises(ValueError,match='custody'):packet.validate(wrong)
