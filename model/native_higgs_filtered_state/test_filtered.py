import importlib.util
from pathlib import Path
from fractions import Fraction as F
import pytest
import sympy as s

spec=importlib.util.spec_from_file_location('filtered_native_probe_public',Path(__file__).with_name('model.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


@pytest.fixture(scope='module')
def result():return m.certificate()


def test_geometric_bound_is_exact_and_sharp():
    a=s.Matrix([1,0]); b=s.Matrix([s.Rational(16,25),s.Rational(12,25)])
    residual=b-a
    assert (residual.H*residual)[0]==s.Rational(9,25)
    overlap_squared=(a.H*b)[0]**2/(b.H*b)[0]
    assert 1-overlap_squared==s.Rational(9,25)


def test_complete_instrument_retains_inside_outside_and_failure():
    C=s.diag(1,s.Rational(1,2),s.Rational(1,3));P=s.diag(1,1,0)
    K=(P*C,(s.eye(3)-P)*C,s.diag(0,s.sqrt(3)/2,2*s.sqrt(2)/3))
    assert sum((x.H*x for x in K),s.zeros(3))==s.eye(3)
    psi=s.ones(3,1)/s.sqrt(3)
    probabilities=[s.simplify((psi.H*x.H*x*psi)[0]) for x in K]
    assert all(p>0 for p in probabilities) and sum(probabilities)==1


def test_imperfect_entangled_preparation_is_normalized_after_filter():
    g=s.Matrix([1,0,0]);env=s.Matrix([1,0]);exc=s.Matrix([0,0,1])
    cos=s.Rational(999999,1000001);sin=s.Rational(2000,1000001)
    psi=cos*s.kronecker_product(g,env)+sin*s.kronecker_product(exc,s.Matrix([0,1]))
    assert (psi.H*psi)[0]==1
    eps=F(1,500)
    ideal=s.kronecker_product(g,env)
    assert (psi-ideal).dot(psi-ideal)<s.Rational(eps.numerator,eps.denominator)**2
    A=s.zeros(3);A[1,0]=1;A[2,2]=1
    B=s.zeros(3);B[1,0]=1;B[2,0]=s.Rational(1,10);B[1,2]=-1
    a=s.kronecker_product(A,s.eye(2))*psi;b=s.kronecker_product(B,s.eye(2))*psi
    infidelity=s.simplify(1-abs((a.H*b)[0])**2/((a.H*a)[0]*(b.H*b)[0]))
    bound=m.filtered_bounds((F(9,10),F(11,10)),(1,1),eps,source_remainder=F(0))
    assert infidelity<=s.Rational(bound['actual_normalized_state_infidelity_upper'])
    assert (a.H*a)[0]>=s.Rational(bound['actual_elementary_in_bin_probability_lower'])
    assert (b.H*b)[0]/4>=s.Rational(bound['actual_composite_B_over_2_in_bin_probability_lower'])


def test_all_history_outcomes_contract_without_length_factor():
    a=s.Matrix([1,0]);b=s.Matrix([s.Rational(4,5),s.Rational(3,5)])
    H=s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    P=(s.diag(1,0),s.diag(0,1));laws=[]
    for vector in (a,b):
        law=[]
        for first in P:
            for second in P:
                for third in P:
                    v=third*H*second*H*first*vector
                    law.append(s.simplify((v.H*v)[0]))
        assert sum(law)==1; laws.append(law)
    tv=sum(abs(x-y) for x,y in zip(*laws))/2
    assert tv==s.Rational(9,25)<s.Rational(3,5)


def test_depth_zero_has_success_one_and_nonzero_state_error():
    p=m.preparation(0)
    assert p['complete_preparation_success_lower']=='1'
    assert p['agreement_filter_count']==0
    assert F(p['retained_purification_error_upper'])>0
    assert p['deterministic_initial_trial'] and p['filter_depth_counts_all_rounds']


def test_actual_two_bins_and_all_preparation_flags(result):
    narrow,wide=result['filtered_state_settings'];p=result['physical_preparation']
    assert p['native_trial_copies']==8 and p['agreement_filter_count']==7
    assert not p['exact_interacting_ground_projection_assumed']
    assert F(narrow['actual_normalized_state_infidelity_upper'])<F(3,10**14)
    assert F(wide['actual_normalized_state_infidelity_upper'])<F(2,10**11)
    assert F(narrow['unconditional_elementary_in_bin_probability_lower'])>F(3147,10000)
    assert F(wide['unconditional_elementary_in_bin_probability_lower'])>F(4197,10000)
    assert not narrow['retains_original_certified_elementary_probability_lower']
    assert wide['retains_original_certified_elementary_probability_lower']
    for row in (narrow,wide):
        assert F(row['positive_cross_residue_lower'])>0
        assert row['spectral_projector_uses_absolute_bin_not_excitation_enclosure']
        assert not row['unconditional_A_and_B_over_2_complete_laws_are_close_claimed']
        assert not row['implemented_bin_QPE_claimed']
    assert result['complete_instruments']['full_setting_pump_energy_alphabet_size']==16


def test_typed_mass_transfer_keeps_selectors_and_hardware_separate(result):
    transfer=result['typed_free_Dirac_corollary']
    assert transfer['actual_positive_cross_G1_eigen_doublet_exists']
    assert transfer['parametric_four_mode_CAR_and_rotation_algebra_replayed']['complete_four_mode_Fock_code_dimension']==16
    assert transfer['spectral_matrix_units_and_full_range_kinetic_operator_must_be_redeclared']
    assert not transfer['old_Dirac_receipt_or_old_spectral_couplers_modified']
    assert not transfer['source_composite_equals_code_field_on_all_occupations']
    assert not transfer['interacting_native_Lorentz_or_SM_chirality_transfer']
    assert m.dirac.four_mode_code()['mass_term_on_code']=='E_star times total canonical particle-plus-antiparticle number'


def test_bounds_reject_lost_normalization_margin():
    with pytest.raises(ValueError):m.filtered_bounds((1,2),(F(1,100),F(1,50)),F(1,2))
    with pytest.raises(ValueError):m.filtered_bounds((0,1),(1,1),0)


def test_public_science_matches_immutable_research_reference(result):
    comparison=result['independent_research_comparison']
    assert comparison['all_scientific_fields_equal']
    assert comparison['fresh_source_and_record_bounds_replayed']
    assert comparison['reference_sha256']=='13200ce38d74e97c6fee01f985032f6d5a3a3c18a0beb8d6aa8aa399d34561a8'
    m.higgs.validate(result['input_sha256'])
