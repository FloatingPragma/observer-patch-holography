"""Exact current, full spectrum and complete instrument regression controls."""
import sympy as s
import pytest
from .model import (tensor_checks,background_checks,quantum_response,
                    quantum_ground_and_readout,neutral_matrix,bdg_symbol)


def test_quantum_scalar_and_Weyl_vertex_are_gauge_covariant():
    value=tensor_checks()
    assert value['full_twelve_current_pair_Ward_checks']==72
    assert value['quantum_scalar_commutator_checks']==24
    assert value['gauge_invariant_reference_Higgs_expectation']=='0'
    assert value['canonical_scalar_CCR_at_cutoff'] is False


def test_background_mass_is_chiral_pairing_with_Nambu_redundancy():
    value=background_checks()
    assert value['canonical_Weyl_modes_per_site']==30
    assert value['Nambu_matrix_dimension']==60
    assert value['analytic_zero_momentum_symbol_eigenvalues']=={'-1':12,'1':12,'-2':12,'2':12,'-3':4,'3':4,'0':4}
    assert neutral_matrix(1,1,1,1)!=neutral_matrix(1,2,3,1)
    with pytest.raises(ValueError):neutral_matrix(1.0,1,1,1)
    assert bdg_symbol((1,2,3),1,2,3,1).is_hermitian


def test_quantum_conversion_retains_entire_active_source():
    value=quantum_response()
    assert value['full_dimension']==320 and value['cyclic_block_dimension']==5
    assert value['all_gauge_current_commutators_checked']==12


def test_full_quantum_ground_and_complete_charged_readout():
    value=quantum_ground_and_readout()
    assert value['full_ground_energy']=='-2' and value['full_ground_rank']==2
    assert value['full_ground_fermion_parity']=='odd'
    assert value['ground_preparation_success_probability']=='1/160'
    assert len(value['complete_reference_preparation_records'])==14
    assert sum(s.sympify(r['probability']) for r in value['complete_reference_preparation_records'])==1
    assert value['neutral_gauge_singlet_reference']['ground_is_full_H_ground'] is False
    records=value['complete_driven_reservoir_energy_records']
    assert len(records)==56 and s.simplify(sum(s.sympify(r['probability']) for r in records))==1
    assert sum(value['full_energy_multiplicities'].values())==320
    conversion=value['vacuum_input_interaction']['complete_basis_probabilities']
    assert len(conversion)==320 and sum(map(s.sympify,conversion))==1
    assert value['nonzero_mean_not_classical_field_control']['field_variance']=='3/4'
    assert value['onsite_gap_is_relativistic_particle_mass'] is False


def test_charged_reservoir_retains_its_Gauss_compensator():
    from .apparatus import certificate
    value=certificate()
    assert value['with_scalar_dimension']==1280
    assert value['physical_reference_preparation']['actual_vectorized_Gauss_checks']==12
    assert value['physical_reference_preparation']['ground_success_probability']=='1/160'
    assert value['with_retained_dual_preparation_reference_dimension']==409600
    assert value['actual_total_current_commutators']==12
    assert value['complete_records_preserve_all_56_source_energy_weights']
    assert value['complete_preparation_and_readout_record_count']==69
    assert sum(s.sympify(r['probability']) for r in value['complete_preparation_and_readout_records'])==1
    assert not value['uncompensated_charged_probe_is_source_only_Gauss_preserving']
