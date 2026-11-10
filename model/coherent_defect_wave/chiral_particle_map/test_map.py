"""Physical-chirality, complete canonical and parameter ancestry checks."""
import sympy as s
from .model import dictionary,physical_currents,spectrum_identity,fock_pair_check,certificate


def test_bijective_bogoliubov_map_and_actual_currents():
    rows,nu,A,B,S=dictionary()
    assert len(rows)==7 and nu==10
    assert A.rank()==16 and B.rank()==14
    assert S*S.H==s.eye(30)
    _,currents=physical_currents()
    assert currents['weak_generators_supported_only_on_left_components']
    assert currents['color_generators_equal_on_physical_left_and_right']
    assert len(currents['charged_current_nonzero_spin_resolved_entries'])==8


def test_massive_dirac_assembly_and_no_independent_nambu_mirrors():
    result=spectrum_identity()
    assert result['Dirac_fields_including_color']==7
    assert result['massless_left_Weyl_fields']==1
    assert not result['Nambu_rows_are_independent_mirror_species']
    assert fock_pair_check()['both_source_pair_signs_match_Dirac_mass_operator']


def test_conditional_roles_do_not_claim_generation_or_mass_selection():
    result=certificate()
    assert len(result['conditional_field_rows'])==8
    assert all(row['empirically_identified_species'] is None for row in result['conditional_field_rows'])
    assert result['identification']['measured_GeV_predictions']=={}
    assert result['identification']['identity_probability_assigned'] is None
    assert result['parameter_controls']['v1_yu1_yd1_ye1']!=result['parameter_controls']['v1_yu2_yd3_ye5']
