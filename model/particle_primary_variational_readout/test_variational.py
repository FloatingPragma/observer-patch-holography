"""Independent characters, full small matrices and actual endpoint authentication."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import variational as v


@pytest.fixture(scope='module')
def member():
    return v.contract.member(0)


@pytest.fixture(scope='module')
def result(member):
    return v.selected_bounds(member)


def test_two_independent_native_character_calculations():
    value = v.character_control()
    assert value['real_character_sum_moments_1_2_3'] == ['0', '10', '54']
    matrix = np.array(value['cubic_loop_singlet_matrix'])
    assert np.array_equal(matrix, matrix.T)
    assert matrix.sum() == 54
    assert matrix.tolist() == [[0, 2, 2, 4, 2], [2, 2, 4, 2, 2], [2, 4, 2, 4, 0],
                               [4, 2, 4, 2, 2], [2, 2, 0, 2, 0]]
    assert sum(x['ordered_multiplicity'] for x in value['singlet_charge_multisets']) == 54


def matrix_and_trial(value):
    data = value['dimensionless_trial']
    norms = np.array([float(F(x)) for x in data['hopping_packet_norms_squared']]
                     +[float(F(data['real_loop_packet_norm_squared']))]*5)
    matrix = np.zeros((11, 11))
    matrix[0, 1:] = matrix[1:, 0] = np.sqrt(norms)
    matrix[1:6, 1:6] = np.diag([float(F(x)) for x in data['hopping_diagonal']])
    matrix[6:, 6:] = np.array([[float(F(x)) for x in row] for row in data['loop_packet_matrix']])
    coefficients = np.array([float(F(x)) for x in data['rational_hopping_coefficients']+data['rational_loop_coefficients']])
    trial = np.r_[1, -np.sqrt(norms)*coefficients]
    return matrix, trial/np.linalg.norm(trial)


def test_exact_schur_root_matches_independent_small_matrix(result):
    matrix, trial = matrix_and_trial(result)
    values = np.linalg.eigvalsh(matrix)
    data = result['dimensionless_trial']
    low, high = map(F, data['lower_Ritz_shift_interval'])
    assert high-low <= F(1, 2**60)
    assert abs(values[0]-float(low)) < 1e-12
    assert sum(values < 0) == 1
    assert abs(trial@matrix@trial-float(F(data['trial_energy_shift']))) < 1e-13
    assert np.linalg.norm(matrix@trial-values[0]*trial) < 1e-12


def test_an_unchanged_compression_does_not_certify_an_invariant_space(result):
    matrix, trial = matrix_and_trial(result)
    full = np.zeros((12, 12)); full[:11, :11] = matrix
    full[-1, -1] = 7
    full[1, -1] = full[-1, 1] = 1
    vector = np.r_[trial, 0]
    rayleigh = vector@full@vector
    residual = np.linalg.norm(full@vector-rayleigh*vector)
    eigenvalues = np.linalg.eigvalsh(full)
    assert eigenvalues[0] < rayleigh
    assert residual > 0.01
    assert min(abs(eigenvalues-rayleigh)) <= residual
    assert not result['dimensionless_trial']['trial_subspace_claimed_invariant']


def test_centered_noncommuting_local_fourth_moment_control():
    X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1); I=s.eye(2)
    terms=[s.kronecker_product(X,I,I),s.kronecker_product(I,X,I),
           s.kronecker_product(Y,Y,I),s.kronecker_product(I,Z,X)]
    vector=s.zeros(8,1); vector[0]=1
    assert all((vector.H*t*vector)[0]==0 for t in terms)
    assert terms[0]*terms[2] != terms[2]*terms[0]
    total=sum(terms,s.zeros(8))
    fourth=(vector.H*total**4*vector)[0]
    assert fourth > 0 and fourth <= v.local_fourth_bound(4,4)


def test_disconnected_covariances_use_variance_without_a_Wick_assumption():
    A=s.Matrix([[0,s.Rational(1,10)],[s.Rational(1,10),1]])
    left=s.kronecker_product(A,s.eye(2)); right=s.kronecker_product(s.eye(2),A)
    omega=s.Matrix([1,0,0,0])
    assert (omega.H*left*right*omega)[0]==0
    exact=(omega.H*(left+right)**4*omega)[0]
    variance=F(1,50)
    estimate=v.orthogonal_packet_fourth(4,2,variance)
    assert exact==s.Rational(13,625)
    assert F(exact)<=estimate['fourth_moment_upper']<v.local_fourth_bound(4,2)


def test_same_source_cubic_and_two_vector_moments(result):
    b=result['actual_bonds']; p=result['actual_plaquettes']
    assert (b,p)==(12,6)
    assert F(result['variance'])==15*b+F(p,4410)
    assert F(result['H0_packet_first_moment'])==720*b+F(11*p,350)
    assert F(result['cubic_perturbation_moment'])==-F(p,171500)
    assert F(result['dimensionless_trial']['trial_energy_shift']) < F(result['dimensionless_two_vector_control']['trial_energy_shift'])


def test_actual_endpoint_ground_bound_and_broad_residual_have_distinct_domains(result):
    upper=F(result['ground_energy_upper'])
    assert F(-724)<upper<F(-723)
    assert upper<F(result['product_trial_ground_energy_upper'])==-720
    assert result['coupling_evaluated']==1
    assert result['first_excitation_gap_or_rest_mass'] is None
    assert not result['trial_is_selected_polynomial_preparation']
    assert not result['residual_interval_identified_as_ground_or_excitation']
    assert not result['original_primary_probe_residue_computed']
    assert F(result['trial_residual_upper'])>100
    assert not result['source_reference_state_and_clock_replaced']


def test_exact_bindings_reject_energy_and_source_mutations(member):
    bound=v.bind_member(member)
    altered=deepcopy(bound); altered['result']['ground_energy_upper']='-1000'
    with pytest.raises(ValueError): v.verify_bound(member,altered)
    altered_member=deepcopy(member); altered_member['state']['other_preparation']=True
    with pytest.raises(ValueError): v.selected_bounds(altered_member)


def test_resolvent_solver_uses_exact_positive_pivots():
    matrix=[[F(2),F(-1)],[F(-1),F(2)]]
    assert v.solve(matrix,[1,1])==[1,1]
    with pytest.raises(ArithmeticError):v.solve([[-1]],[1])
    with pytest.raises(ValueError):v.local_fourth_bound(-1,1)


def test_generic_two_vector_trial_keeps_the_omitted_space():
    # H0=diag(0,3,5), V couples0↔1 and1↔2; the first trial space is not invariant.
    H=s.Matrix([[0,1,0],[1,3,1],[0,1,5]])
    result=v.ritz_trial(1,3,9,2,0)
    t=result['rational_trial_t']; norm=1+t*t
    psi=s.Matrix([1,-s.Rational(t.numerator,t.denominator),0])/s.sqrt(s.Rational(norm.numerator,norm.denominator))
    e=s.Rational(result['trial_energy_shift'].numerator,result['trial_energy_shift'].denominator)
    residual=(H*psi-e*psi).norm()**2
    assert F(residual)<=result['trial_residual_squared_upper']
    assert residual>0
