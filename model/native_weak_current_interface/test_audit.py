import importlib.util
from pathlib import Path
from fractions import Fraction as F
import pytest
import sympy as s

spec=importlib.util.spec_from_file_location('native_weak_current_public_control',Path(__file__).with_name('model.py'))
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


@pytest.fixture(scope='module')
def result():return audit.certificate()


def test_full_current_map_needs_no_scalar_state(result):
    interface=result['exact_source_current_interface']
    assert interface['same_full_CAR_automorphism_independent_of_scalar_state']
    assert interface['actual_current_map_controls']['weak_generators_supported_only_on_left_components']
    assert not interface['state_or_Hamiltonian_replaced_by_field_variable_change']
    assert not interface['constant_background_BdG_intertwiner_applies_to_quantum_scalar_without_error']


def test_literal_nilpotent_scalars_and_strict_background_residual(result):
    data=result['fixed_scalar_carrier_obstruction']
    matrices=[s.Matrix(x) for x in data['actual_phi_matrices']]
    assert all(x**3==s.zeros(5) and x**2!=s.zeros(5) for x in matrices)
    assert data['common_phi_and_adjoint_kernel_dimension']==0
    assert F(data['native_unit_Yukawa_shared_energy_as_background_residual_lower_exact'])>F(15199,100000)


def test_obstruction_does_not_exclude_collective_expectations(result):
    control=result['fixed_scalar_carrier_obstruction']['explicit_expectation_control']
    assert control['means_phi0_phi1_and_adjoints']==['sqrt(2)/2','0','sqrt(2)/2','0']
    assert control['one_copy_variances']==['1/4','1/2','1/4','1/2']
    phi=audit.higgs_memory()[0][0];eta=s.Matrix([s.sqrt(2),1,0,1,0])/2
    two=s.kronecker_product(eta,eta)
    average=(s.kronecker_product(phi,s.eye(5))+s.kronecker_product(s.eye(5),phi))/2
    residual=(average-s.eye(25)/s.sqrt(2))*two
    assert s.simplify((residual.H*residual)[0])==s.Rational(1,8)


def test_actual_dressed_current_operator_norm_and_scope(result):
    data=result['exact_Higgs_dressed_weak_current_interface']
    assert data['all_twelve_Gauss_current_commutators_zero'] and data['spin_scalar_and_even_parity']
    assert data['actual_full_operator_norm']=='2'
    assert data['top_eigenprojector_rank_on_active_carrier']==2
    assert not data['actual_eight_cell_ground_expectation_or_form_factor_calculated']


def test_identical_family_copies_do_not_split_into_three_masses(result):
    data=result['family_copy_symmetry']
    assert data['three_copy_permutation_commutant_dimension']==2
    assert data['one_particle_family_eigenvalues']==['d+2b','d-b','d-b']
    assert not data['family_breaking_operator_selected_by_copy_count']


def test_public_provider_and_frozen_research_comparison(result):
    W,D,S,Ka,nL=audit.dressed_weak_polynomial()
    assert W and len(Ka)==3 and S.rows==5
    assert audit.dressed_weak_matrix().shape==(80,80)
    comparison=result['independent_research_comparison']
    assert comparison['all_scientific_fields_equal'] and comparison['fresh_calculation_was_replayed']
    assert comparison['reference_sha256']=='63bc7d02fbd19066ed6afb09e79713d6ba330073838a6b874568f39172a1bb81'
