"""Exact source spectral hierarchy and original-probability error controls."""
from fractions import Fraction as F
from math import comb
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import gamma
import sympy as s
from . import branch


def test_block_krylov_complex_moments_and_dependent_columns():
    k = s.diag(0, 1, 2, 4, 7)
    seed = s.Matrix([[1, 0], [s.I, 1], [1, -s.I], [0, 1], [1, 1]])/3
    row = branch.krylov(k, seed, 2, F(1))
    assert row['rank'] == 4 and row['moments_matched_through'] == 3
    assert row['Gram_floor'] > 0
    assert row['residual_Gram'] != s.zeros(4)
    complete = branch.krylov(k, seed, 3, F(1))
    assert complete['rank'] == 5 < complete['input_column_count'] == 6
    assert complete['residual_Gram'] == s.zeros(5)
    assert complete['moments_matched_through'] == 5


def test_compression_is_quadrature_not_an_exact_energy_assignment():
    k, seed = s.diag(1, 3), s.Matrix([1, 1])/s.sqrt(2)
    one = branch.krylov(k, seed, 1)
    assert one['compression_in_Gram_coordinates'][0] == s.Rational(3, 8)
    assert s.Rational(3, 8) not in one['S'].eigenvals()
    two = branch.krylov(k, seed, 2)
    assert two['residual_Gram'] == s.zeros(2)
    assert set(two['compression_in_Gram_coordinates'].eigenvals()) == {s.Rational(1, 2),s.Rational(1, 4)}


def test_gram_floor_controls_actual_coherent_preparation_probability():
    h, omega, _, x, z = branch.closure._native_dimer()
    k = h+s.eye(4)*9/8
    bounded = (k+s.eye(4)/8).inv()/8
    psi = (s.eye(4)-bounded)*(x+z)*omega/2
    row = branch.krylov(k, psi, 2, F(1, 8))
    for coeff in (s.Matrix([1, 1]),s.Matrix([1, -2]),s.Matrix([s.I, 2])):
        alpha = sum(abs(c) for c in coeff)
        probability = s.simplify((coeff.H*row['Gram']*coeff)[0]/alpha**2)
        assert probability >= row['polynomial_preparation_original_success_lower'] > 0
        assert probability < 1
    floor = row['Gram_floor']
    assert branch.gram_floor(row['Gram'],floor/4) == floor/2
    with pytest.raises(ValueError):
        branch.gram_floor(row['Gram'],floor/2)


@pytest.mark.parametrize('matrix', [s.zeros(2),s.diag(1,-1),s.Matrix([[1,2],[0,1]]),s.Matrix([[1,1],[1,1]])])
def test_gram_rank_or_positivity_cannot_be_guessed(matrix):
    with pytest.raises(ValueError):
        branch.gram_floor(matrix)


def test_bernstein_power_conversion_against_independent_binomial_sum():
    values = (F(0), F(1,4), F(3,4), F(1), F(1,2))
    coefficients = branch.bernstein_coefficients(values)
    for u in (F(0),F(1,7),F(1,2),F(1)):
        direct = sum(comb(4,j)*u**j*(1-u)**(4-j)*value for j,value in enumerate(values))
        power = sum(c*u**j for j,c in enumerate(coefficients))
        assert direct == power
        assert 0 <= power <= 1


def test_original_window_brackets_with_boundary_atom_and_moment_errors():
    support, weights = (F(1,4), F(1,2), F(9,10)), (F(1,5), F(2,5), F(1,10))
    degree=64
    moments = tuple(sum(w*u**j for u,w in zip(support,weights,strict=True)) for j in range(degree+1))
    exact = weights[0]+weights[1]
    row = branch.window_from_moments(moments,F(1,4),F(3,4),F(1,4))
    assert row['original_interval_lower'] <= exact <= row['original_interval_upper'] <= sum(weights)
    assert row['original_interval_lower'] > 0
    eps=F(1,10**70)
    errors=(F(0),)+(eps,)*degree
    perturbed = (moments[0],)+tuple(v+(-1)**j*eps for j,v in enumerate(moments[1:],1))
    measured = branch.window_from_moments(perturbed,F(1,4),F(3,4),F(1,4),errors)
    assert measured['original_interval_lower'] <= exact <= measured['original_interval_upper']
    assert measured['lower_score_error'] > eps
    assert row['boundary_weight_proved_zero'] is row['spectral_atom_inferred'] is False
    # At the interval endpoint the limiting ramp gap retains the entire atom.
    values_lower = branch.ramp_values(4,F(1,4),F(3,4),F(1,8))
    values_upper = branch.ramp_values(4,F(1,4),F(3,4),F(1,8),upper=True)
    assert values_lower[1] == 0 and values_upper[1] == 1


def test_positive_fourier_mixture_and_tail_against_gamma_integrals():
    for alpha,z,energy,T in ((.5,1.,.7,20.),(1.5,.5,2.,40.),(3.,2.,1.,10.)):
        density=lambda u:z**alpha*u**(alpha-1)*np.exp(-z*u)/gamma(alpha)
        laplace=quad(lambda u:density(u)*np.exp(-energy*u),0,np.inf,epsabs=1e-11)[0]
        assert abs(laplace-(z/(z+energy))**alpha)<1e-10
        tail=quad(lambda u:density(u)*(2/np.pi)*np.arctan(u/T),0,np.inf,epsabs=1e-11)[0]
        assert 0 < tail <= 2*alpha/(np.pi*z*T)
        row=branch.fourier_preparation_budget(F(str(alpha)),F(str(z)),F(str(T)),F(1,10**5),2)
        assert tail <= float(row['Fourier_tail_upper'])


def test_local_time_discretization_selected_column_bound_on_actual_source():
    # Alpha=1/2 implements a moment as an original success norm, without a projection.
    h,omega,_,x,z=branch.closure._native_dimer()
    op=(x+z)/2
    hn=np.array(h,float); an=np.array(op,float); ground=np.array(omega,float).reshape(-1)
    energies, vectors=np.linalg.eigh(hn)
    comm=np.linalg.norm(hn@an-an@hn,2)
    assert comm <= 2
    # Direct Gamma Laplace quadrature independently recovers the positive spectral filter.
    gaps=energies-energies[0]
    scalar=np.sqrt(1/(1+gaps))
    exact=vectors@(scalar*(vectors.T@(an@ground)))
    integral=[]
    for gap in gaps:
        value=quad(lambda u:np.exp(-u*(1+max(0.,gap)))/np.sqrt(np.pi*u),0,np.inf,epsabs=1e-11)[0]
        integral.append(value)
    reconstructed=vectors@(np.array(integral)*(vectors.T@(an@ground)))
    assert np.linalg.norm(exact-reconstructed)<1e-10
    budget=branch.fourier_preparation_budget(F(1,2),1,10000,F(1,10**6),2)
    assert budget['selected_creation_column_error'] == F(1,30000)+F(1,10**6)
    assert budget['arbitrary_complement_channel_error_inferred'] is False


def test_positive_gap_filter_does_not_supply_a_uniform_joint_limit_gap():
    fixed=branch.positive_gap_filter(F(1,32),1,64)
    assert 0 < fixed['ground_projection_operator_error'] < F(1,5)
    errors=[]
    for n in (2,4,8):
        row=branch.positive_gap_filter(F(1,8*n**3),1,64)
        errors.append(row['ground_projection_operator_error'])
        assert row['continuum_gap_inferred'] is False
    assert errors[0]<errors[1]<errors[2]<1


def test_joint_dimer_labels_preserve_failure_and_zero_weight_records():
    row=branch.dimer_joint_response()
    weights=[F(r['original_weight']) for r in row['complete_joint_success_rows']]
    assert weights==[F(81,4000),F(729,4000),F(16,121),0]
    assert sum(weights)+F(row['preparation_failure_probability'])+F(row['filter_failure_probability'])==1
    assert row['joint_block_rank']==3 and row['joint_block_input_columns']==4
    assert F(row['joint_polynomial_preparation_original_success_lower'])>0
    assert row['joint_translation_compression_commutes_with_clock'] is True
    assert row['translation_is_only_a_finite_graph_symmetry'] is True
    assert row['continuum_pole_certified'] is False
    assert row['physical_mass_assignment'] is None


def test_source_mutation_invalidates_joint_spectrum_anchor(monkeypatch):
    original=branch.closure.source.matrix
    def altered(*args,**kwargs):
        raw,scale=original(*args,**kwargs)
        return raw+np.eye(len(raw),dtype=raw.dtype),scale
    monkeypatch.setattr(branch.closure.source,'matrix',altered)
    with pytest.raises(AssertionError):
        branch.dimer_joint_response()


def test_invariant_interval_is_conditional_and_does_not_fit_speed_or_mass():
    row=branch.invariant_interval(3,4,1,2,1)
    assert row['energy_momentum_invariant_interval']==[5,15]
    assert row['conditional_squared_mass_interval']==[5,15]
    assert row['physical_particle_mass_assigned'] is False
    assert row['Lorentz_action_supplied_by_interval'] is False
    mixed=branch.invariant_interval(1,2,1,2,1)
    assert mixed['energy_momentum_invariant_interval']==[-3,3]
    assert mixed['conditional_squared_mass_interval'] is None


@pytest.mark.parametrize('args', [(0,1,0,1,0),(-1,1,0,1,1),(2,1,0,1,1),(0,1,2,1,1)])
def test_invalid_energy_momentum_intervals_rejected(args):
    with pytest.raises(ValueError):
        branch.invariant_interval(*args)


def test_complete_creator_probability_and_direct_difference_normalization():
    with pytest.raises(ValueError):
        branch.krylov(s.eye(2),s.Matrix([2,0]),1)
    row=branch.positive_gap_filter(F(1,32),1,64)
    assert row['direct_difference_LCU_normalization']==2
    assert row['direct_difference_LCU_retention_lower']*4==row['excited_original_weight_retention_factor_lower']
    assert row['unit_normalized_filter_native_dilation_executed'] is False
