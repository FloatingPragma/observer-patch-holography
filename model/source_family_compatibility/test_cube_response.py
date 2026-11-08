"""Exact source, fidelity, complement and matrix-error controls for the cube."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy import sparse

from model.interacting_vacuum import global_ground
from . import composite
from . import cube_response as cube


@pytest.fixture(scope='module')
def data():
    ground=cube.certified_ground()
    response=cube.response_on_trial(ground['full_trial'],ground['norm_squared'])
    bonds=cube.bond_response_on_trial(ground['full_trial'],ground['norm_squared'])
    return ground,response,bonds,cube.enclosure(ground,response,bonds)


def test_actual_default_full_ground_and_exact_fidelity(data):
    ground,_,_,_=data
    assert ground['full_trial'].shape==(65536,)
    assert len(ground['balanced_indices'])==4900
    assert ground['infidelity_upper']<F(1,1920)**2
    assert ground['half_trace_distance_upper']**2>=ground['infidelity_upper']
    assert ground['full_gap_premise']['full_centered_gap_strict_lower']=='1/128'
    assert ground['ground_energy_interval'][0]<ground['ground_energy_interval'][1]<F(-1055,256)
    assert np.array_equal(ground['full_trial'][ground['balanced_indices']],ground['trial'])


def test_source_clock_restricted_on_exact_ordered_addresses(data):
    ground,_,_,_=data;indices=ground['balanced_indices'];h=global_ground.source_block(8)
    assert not (ground['centered_numerator'][indices][:,indices]-h+
                64*sparse.eye(4900,dtype=np.int64)).nnz
    assert (ground['centered_numerator']-ground['centered_numerator'].T).nnz==0


def test_literal_pi_flux_source_and_force_norms():
    ops=cube._operators();one=ops['one_numerator']
    assert np.array_equal(one@one,3*np.eye(16,dtype=np.int64))
    assert np.count_nonzero(one)==48
    assert ops['Y_norm_upper']==(F(9,16),)*16
    unsigned=np.abs(one)
    assert not np.array_equal(unsigned@unsigned,one@one)


def test_bond_norm_is_exact_full_two_site_algebra():
    small=composite.source((2,1,1),F(1,8),F(1))
    operator=composite.bond_observable(small,0,1)
    identity=np.eye(16,dtype=np.int64)
    assert composite.psd_rank(3*identity-operator)>=0
    assert composite.psd_rank(3*identity+operator)>=0
    assert np.any(operator@operator-identity)


def test_source_bond_and_force_readouts_enclose_same_C(data):
    ground,response,bonds,enclosed=data
    assert len(bonds['four_spin_eta_bond_expectations'])==12
    for value in bonds['four_spin_eta_bond_expectations'].values():
        assert F(-157,100)<value-6*F(1,1920)
        assert value+6*F(1,1920)<F(-155,100)
    assert enclosed['C_operator_error_upper']==F(3,5120)
    assert enclosed['Y_Gram_operator_error_upper']==F(27,2560)
    assert cube.enclosure(ground,response,bonds)['residual_spectrum_interval']==(F(3,80),F(61,1000))


def test_full_residual_positive_in_every_mode_and_all_walsh_probes(data):
    _,_,_,enclosed=data;R=enclosed['trial_residual_Gram'];error=enclosed['residual_operator_error_upper']
    points=global_ground.cube.GEOMETRY.sites
    for k,spin in product(product((0,1),repeat=3),(0,1)):
        f=[0]*16
        for x,point in enumerate(points):
            f[2*x+spin]=(-1)**sum(a*b for a,b in zip(k,point,strict=True))
        norm=sum(x*x for x in f)
        value=sum((f[i]*R[i][j]*f[j]for i in range(16)for j in range(16)),F())/norm
        assert F(3,80)<value-error<value+error<F(61,1000)


def test_both_original_charge_probabilities_are_retained(data):
    ground,_,_,_=data;p=ground['full_trial'];norm=ground['norm_squared']
    for a in cube._operators()['a']:
        plus=F(cube._norm2(cube._action(a.T,p)),norm)
        minus=F(cube._norm2(cube._action(a,p)),norm)
        assert plus+minus==1
        assert abs(plus-F(1,2))<F(1,960)
        assert abs(minus-F(1,2))<F(1,960)


def test_removed_residual_is_rejected(data):
    ground,response,bonds,_=data;mutated=deepcopy(response)
    mutated['Y_Gram']=cube._square(bonds['C'])
    with pytest.raises(ArithmeticError,match='residual'):
        cube.enclosure(ground,mutated,bonds)


def test_halved_complete_force_score_is_rejected(data):
    ground,response,bonds,_=data;mutated=deepcopy(response)
    mutated['Y_Gram']=tuple(tuple(x/2 for x in row)for row in response['Y_Gram'])
    with pytest.raises(ArithmeticError,match='residual'):
        cube.enclosure(ground,mutated,bonds)


def test_wrong_independent_C_score_is_rejected(data):
    ground,response,bonds,_=data;mutated=deepcopy(response)
    C=[list(row)for row in response['C']];C[0][2]+=1
    mutated['C']=tuple(map(tuple,C))
    with pytest.raises(ArithmeticError,match='readouts disagree'):
        cube.enclosure(ground,mutated,bonds)


def test_loose_fidelity_cannot_certify_fixed_response_interval(data):
    ground,response,bonds,_=data;mutated=deepcopy(ground)
    mutated['infidelity_upper']=F(1,1000)**2
    with pytest.raises(ArithmeticError,match='fidelity'):
        cube.enclosure(mutated,response,bonds)


def test_trial_and_clock_mutations_are_rejected(data):
    ground,_,_,_=data;h=global_ground.source_block(8)
    with pytest.raises(ValueError,match='literal'):
        cube.verify_trial(h+sparse.eye(4900,dtype=np.int64),ground['trial'])
    bad=np.zeros(4900,dtype=np.int64);bad[0]=1
    with pytest.raises(ArithmeticError,match='separation'):
        cube.verify_trial(h,bad)
    with pytest.raises(ValueError,match='norm'):
        cube.response_on_trial(ground['full_trial'],ground['norm_squared']+1)


def test_reference_arrays_are_detached(data):
    ground,_,_,_=data;copy=cube.certified_ground()
    copy['full_trial'][0]=7
    copy['trial'][0]+=11
    copy['full_gap_premise']['full_ground_multiplicity']=9
    fresh=cube.certified_ground()
    assert np.array_equal(fresh['trial'],ground['trial'])
    assert np.array_equal(fresh['full_trial'],ground['full_trial'])
    assert fresh['full_gap_premise']['full_ground_multiplicity']==1


def test_exact_sqrt_and_integer_domain_controls():
    assert cube.upper_sqrt(F(9,16),100)==F(3,4)
    assert cube.upper_sqrt(F(2),100)**2>=2
    with pytest.raises(ValueError):cube.upper_sqrt(-1)
    huge=np.array([2**62,2**62],dtype=np.int64)
    with pytest.raises(OverflowError):cube._dot(huge,huge)


def test_receipt_scope_and_controlled_fourth_moment(data):
    result=cube.certificate()
    assert result['cyclic_response_dimension_lower']==48
    assert result['fourth_original_moment_defect_operator_interval']==(F(3,320),F(61,4000))
    assert result['original_response_unit_time_error']['original_response_vector_error_squared_upper']==F(61,16000)
    assert result['each_local_original_addition_removal_weight']==F(1,2)
    assert result['normalized_smearing_total_signed_weight']==1
    assert not result['all_smearings_have_equal_charge_weights']
    assert result['physical_mass'] is None
    assert not result['retarded_compression_is_arbitrary_reference_physical_channel']
    assert not result['native_ground_preparation_executed']
    assert not result['nonzero_momentum_pole_identified']
    assert result['physical_species_added']==result['extra_independent_CAR_flavours']==0
