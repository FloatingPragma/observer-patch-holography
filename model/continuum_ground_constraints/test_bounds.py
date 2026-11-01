from fractions import Fraction as F

import pytest

from . import bounds as b


def test_optimized_threshold_identity_and_density_control():
    for delta in (F(1,2),F(1,8),F(1,64)):
        epsilon=24*delta**4
        assert b.occupation_defect(epsilon,delta)==min(2,8*delta**3)
    result=b.best_dyadic(F(1,8192))
    assert result['threshold']==F(1,16)
    assert result['occupation_defect_upper']==F(1,1024)
    assert result['density_lower']==F(1023,1024)
    assert result['per_spin_density_error_upper']==F(1,2048)


def test_rational_square_root_bound_is_outward_and_tight():
    for value in (F(0),F(1,4),F(1,3),F(2),F(999,71)):
        upper=b.sqrt_upper(value,20)
        assert upper*upper>=value
        assert upper==0 or (upper-F(1,2**20))**2<value


def test_Wick_recursion_counts_two_boundary_errors_and_all_contractions():
    assert [b.wick_error_coefficient(n) for n in (0,2,4,6,8)]==[0,2,8,42,296]
    for n in range(2,16,2):
        assert b.wick_error_coefficient(n)==2+(n-1)*b.wick_error_coefficient(n-2)


def test_positive_measure_loses_second_moment_despite_uniform_bound():
    for n in (2,3,8,100):
        row=b.moment_escape(n);p=row['remote_weight']
        assert 0<p<1 and row['first_moment']==p*n
        assert row['second_moment']==p*n*n==1
        assert row['memory_Gram']==row['second_moment']-row['first_moment']**2
        assert row['weak_limit_second_moment']==row['weak_limit_memory_Gram']==0
        # Direct resolvent calculation at z=i: squared error is rational.
        error_squared=p*p*n*n/F(1+n*n)
        assert error_squared<=row['resolvent_error_at_i_upper']**2
    assert b.moment_escape(100)['resolvent_error_at_i_upper']<F(1,1000)


def test_small_integrated_defect_does_not_bound_a_scaling_packet():
    row=b.momentum_concentration(100)
    assert row['integrated_occupation_defect']==F(1,10**6)
    assert row['packet_defect_expectation']==1


@pytest.mark.parametrize('call',[
    lambda:b.occupation_defect(1,0),lambda:b.occupation_defect(1,2),
    lambda:b.occupation_defect(-1,F(1,2)),lambda:b.occupation_defect(0.1,F(1,2)),
    lambda:b.best_dyadic(0),lambda:b.best_dyadic(1,True),
    lambda:b.wick_error_coefficient(3),lambda:b.moment_escape(1),
    lambda:b.momentum_concentration(True),lambda:b.sqrt_upper(1,0)])
def test_domain_errors_are_rejected(call):
    with pytest.raises(ValueError):call()
