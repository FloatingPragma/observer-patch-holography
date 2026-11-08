"""Independent moment, complete-character and optimal-defect controls."""
from fractions import Fraction as F
import numpy as np
import pytest
from . import dispersion


def test_conditional_variance_is_exact_scalar_fit_minimum():
    weights=(F(1,7),F(2,7),F(4,7));energies=(F(1),F(2),F(5))
    mean=sum(w*e for w,e in zip(weights,energies))
    variance=sum(w*(e-mean)**2 for w,e in zip(weights,energies))
    for candidate in (F(-1),F(),mean,F(10)):
        loss=sum(w*(e-candidate)**2 for w,e in zip(weights,energies))
        assert loss==variance+(candidate-mean)**2
    assert variance>0


def test_outward_variance_intervals_and_degenerate_weight_rejection():
    n=F(1,3);first=F(-2,5);second=F(1)
    exact=second-first*first/n
    interval=dispersion.variance_interval((n,n),(first,first),(second,second))
    assert interval==(exact,exact)
    epsilon=F(1,1000)
    lo,hi=dispersion.variance_interval((n-epsilon,n+epsilon),(first-epsilon,first+epsilon),(second-epsilon,second+epsilon))
    assert lo<exact<hi
    with pytest.raises(ValueError):dispersion.variance_interval((0,n),(first,first),(second,second))


def test_actual_source_complete_eight_record_instrument():
    op,projectors=dispersion.source_projectors()
    summed=np.zeros((256,256),dtype=np.int64)
    for a in op['twice_creation'].values():
        for p in projectors:
            k=p@a;summed+=k.T@k
    assert np.array_equal(summed,64*np.eye(256,dtype=np.int64))
    omitted=projectors[0]@op['twice_creation'][1]
    assert not np.array_equal(summed-omitted.T@omitted,summed)
    # An odd translation preserves H but reverses Q, so it cannot be
    # silently added as a charge-preserving character resolution.
    perm=np.array([sum(((word>>j)&1)<<(j^1)for j in range(8))for word in range(256)])
    assert np.array_equal(op['H'][np.ix_(perm,perm)],op['H'])
    assert np.array_equal(op['twice_charge'][np.ix_(perm,perm)],-op['twice_charge'])


def test_exact_source_moment_intervals_against_independent_eigensolver():
    value=dispersion.build();op,projectors=dispersion.source_projectors()
    h=op['H']/op['scale'];energies,vectors=np.linalg.eigh(h);ground=vectors[:,0]
    for row in value['channels']:
        q=row['internal_charge'];created=op['twice_creation'][q]@ground/2
        residual=0
        for detail,p in zip(row['characters'],projectors):
            v=p@created/4;n=v@v;first=v@h@v;second=(h@v)@(h@v)
            variance=second/n-(first/n)**2;residual+=n*variance
            lo,hi=map(float,detail['conditional_energy_variance_interval'])
            assert lo-1e-12<=variance<=hi+1e-12
            lo,hi=map(float,detail['conditional_excitation_mean_interval'])
            assert lo-1e-12<=first/n-energies[0]<=hi+1e-12
        lo,hi=map(float,row['optimal_unnormalized_squared_dispersion_defect_interval'])
        assert lo-1e-12<=residual<=hi+1e-12
        assert residual>float(F(row['strict_simple_original_weight_defect_lower']))
    assert not value['source_spectral_subbranch_or_distinct_interpolator_excluded']
    assert not value['physical_dispersion_or_mass_fit_supplied']
