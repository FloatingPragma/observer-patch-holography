from fractions import Fraction as F

import numpy as np
import pytest

from . import algebra as a


def test_complete_local_and_disjoint_CAR_identities():
    value=a.contact_identities()
    assert value['onsite_Fock_dimension']==4
    assert value['all_spin_pairs_checked']==4
    assert value['two_site_Fock_dimension']==16
    assert value['different_site_mixed_pairs_checked']==4


def test_spin_coherence_cannot_be_dropped_from_the_memory_formula():
    # A pure coherent one-fermion state has rho^2=rho and therefore R^2=R.
    # Keeping only the diagonal would incorrectly predict nonzero memory.
    rho=np.array([[1,1],[1,1]],dtype=np.int64)
    r=2*np.eye(2,dtype=np.int64)-rho  # numerator with denominator two
    assert np.array_equal(r@r,2*r)
    diagonal=np.diag(np.diag(r))
    assert np.any(2*diagonal-diagonal@diagonal)


def test_exhaustive_period_four_frame_translation_normalizes_SU2():
    value=a.frame_translation_certificate()
    assert value['exact_translation_cases']==4096
    assert value['central_sign_counts']=={'+':2560,'-':1536}
    assert value['ground_state_functional_faithfulness_inferred'] is False


def test_selected_density_interval_gives_strict_memory_and_keeps_boundaries():
    for j in range(8,33):
        factor=a.memory_factor(F(j,64))
        assert F(7,64)<=factor<=F(1,4)
    assert a.memory_factor(F(1,8))==F(7,64)
    assert a.memory_factor(F(1,2))==F(1,4)
    for bad in (0,1,F(1,16),F(3,4),0.25,True):
        with pytest.raises(ValueError):a.memory_factor(bad)


def test_nonzero_response_memory_is_not_a_scattering_or_norm_limit_claim():
    value=a.certificate()
    assert value['strictly_nonzero_infinite_ground_response_memory'] is True
    assert value['same_CAR_number_conserving_static_quadratic_response_excluded'] is True
    assert value['nontrivial_scattering_proved'] is False
    assert value['operator_norm_large_frequency_limit_claimed'] is False


def test_pairing_quadratic_normal_subblock_is_outside_the_excluded_class():
    # A Nambu matrix [[0,Delta],[Delta,0]] is quadratic but not number
    # conserving. Its normal resolvent has nonzero Delta^2/z memory.
    z,delta=F(2),F(1)
    normal=z/(z*z-delta*delta)
    assert z-1/normal==delta*delta/z!=0
    assert 'same_CAR_static_quadratic_response_excluded' not in a.certificate()
