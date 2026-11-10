"""Full residual admission and independent finite-exponential controls."""
from fractions import Fraction as F
import math
import numpy as np
import mpmath as mp
import pytest

from . import certificate as cert
from . import operator as op


@pytest.fixture(scope='module')
def receipt():return cert.certificate()


def test_full_source_ground_and_complete_records(receipt):
    lo,hi=map(F,receipt['actual_charged_ground_energy_interval'])
    assert F('1.873809754965')<lo<hi<F('1.873809754968')
    assert hi-lo<=F(21,10**13)
    assert receipt['physical_one_charge_dimension']==13*5456
    assert receipt['additional_state_truncation'] is False
    assert F(receipt['energy_lowering_from_full_transverse_dressing_interval'][0])>F('.000175')
    p,q=map(F,receipt['actual_ground_nonvacuum_transverse_probability_interval'])
    assert F('.00007224')<p<q<F('.00007225')
    positions=receipt['actual_ground_position_records']
    assert len(positions)==13
    assert sum(F(x['actual_ground_position_probability'][0]) for x in positions)<=1
    assert sum(F(x['actual_ground_position_probability'][1]) for x in positions)>=1
    prep=receipt['complete_ground_preparation']
    assert F(prep['unconditional_ground_window_probability_interval'][0])>F('.24998')
    assert prep['pump_failure_probability']=='3/4'
    assert receipt['scope']['mass_GeV'] is None


def test_trial_range_rejects_tiny_or_incomplete_vectors():
    with pytest.raises(ValueError):cert.norm_admission(np.zeros(70928))
    with pytest.raises(ValueError):cert.norm_admission(np.ones(13))
    bad=np.zeros(70928);bad[0]=math.nan
    with pytest.raises(ValueError):cert.norm_admission(bad)


def test_every_spectator_block_matches_literal_exponential():
    mp.mp.dps=65;mp.iv.dps=80
    for edge in (0,12):
        variance=op.ivalue(op.parent.exact_peierls_data(edge)[2])
        coefficients=op.peierls_coefficients(variance)
        cc=[mp.mpf(str((op.parent.endpoint(c)+op.parent.endpoint(c,True))/2)) for c in coefficients]
        sigma=mp.sqrt(mp.mpf(str(op.parent.exact_peierls_data(edge)[2].evalf(70))))
        for size in range(1,5):
            a=mp.zeros(size)
            for n in range(1,size):a[n-1,n]=mp.sqrt(n)
            D=sigma*(a.T-a)
            polynomial=sum((c*D**k for k,c in enumerate(cc)),mp.zeros(size))
            exact=mp.expm(D)
            assert mp.norm(polynomial-exact)<mp.mpf('1e-60')


def test_published_radius_requires_actual_full_replay(receipt):
    diagnostics=cert.replay()[3]
    assert 0<diagnostics['computed_actual_residual_upper']<float(cert.RADIUS)
    assert diagnostics['computed_H_evaluation_budget']>0
    assert receipt['arithmetic_admission']['numerical_eigensolver_used_in_verifier'] is False
