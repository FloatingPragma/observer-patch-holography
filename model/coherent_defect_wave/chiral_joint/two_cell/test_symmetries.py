import sympy as s
from .symmetries import (open_edge_obstruction,spin_obstruction,periodic_checks,
                        momentum_projector,periodic_translation,apply_column)


def test_single_quantum_link_cannot_reverse_all_raising_channels():
    result=open_edge_obstruction()
    assert result['link_common_domain_rank']==1
    assert result['link_total_range_rank']==5
    assert not result['single_link_unitary_or_antiunitary_orientation_reversal']


def test_axial_spin_conservation_is_not_a_full_rotation_group():
    result=spin_obstruction()
    assert result['coefficient']=='-I'
    assert not result['physical_little_group_spin_identified']


def test_two_quantum_links_restore_complete_AP_momentum():
    result=periodic_checks()
    assert result['declared_periodic_full_dimension']==3686400
    assert result['complete_four_phase_readout_on_odd_one_particle_input']==['0','1/2','0','1/2']
    assert not result['fixed_charged_boundary_preserved_by_translation']
    assert result['translated_boundary_restores_Gauss_intertwiner']


def test_phase_readout_retains_even_sector_zero_outcomes():
    seed={(0,0,0,0,0):s.Integer(1)}
    assert apply_column(seed,periodic_translation)==seed
    assert momentum_projector(seed,0)==seed
    assert all(not momentum_projector(seed,k) for k in (1,2,3))


def test_actual_antiunitary_maps_boundary_and_probe_together():
    from .symmetries import antiunitary_checks
    result=antiunitary_checks()
    assert result['complete_CAR_basis_columns_checked']==49152
    assert result['same_open_two_cell_H_antiunitary_symmetry']
    assert result['probe_map']=='right addition -> left removal, up to phase -I'
    assert not result['same_side_addition_removal_equality_implied']


def test_fixed_boundary_odd_sectors_have_exact_Kramers_pairs():
    from .symmetries import spin_reversal_checks
    result=spin_reversal_checks()
    assert result['complete_CAR_columns_checked']==49152
    assert result['odd_L_energy_degeneracy_at_least_two_exact']
    assert result['sector_map']['boundary']=='unchanged'


def test_physical_periodic_pair_preparation_has_complete_momentum_records():
    from .symmetries import physical_periodic_momentum_control
    result=physical_periodic_momentum_control()
    assert result['all_24_local_Gauss_actions_annihilate_state']
    assert result['complete_phase_probabilities']==['1/2','0','1/2','0']
