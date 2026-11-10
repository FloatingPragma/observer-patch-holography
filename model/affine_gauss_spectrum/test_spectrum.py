"""Independent finite exponential and full-complement spectral controls."""
from fractions import Fraction as F
import mpmath as mp
import numpy as np
import pytest

from . import model


@pytest.fixture(scope='module')
def result():
    return model.certificate()


def test_all_low_levels_and_preparation_have_a_full_complement(result):
    data=result['spectral'];rows=data['low_charged_levels']
    assert len(rows)==13
    lower,upper=map(F,rows[0]['actual_H_energy_interval'])
    assert F('1.86')<lower<upper<F('1.88')
    assert F(data['one_charge_ground_gap_lower'])>F('.12')
    assert F(data['strict_ground_energy_reduction_by_photon_dressing_lower'])>F('0.00000026')
    assert 0<F(data['ground_nonvacuum_transverse_probability_lower'])<F(data['ground_nonvacuum_transverse_probability_upper'])<F('.01')
    assert F(rows[-1]['actual_H_energy_interval'][1])<F(data['fourteenth_and_higher_energy_lower'])
    prep=data['ground_preparation']
    hit=list(map(F,prep['ground_window_hit_probability_interval']))
    other=list(map(F,prep['energy_complement_probability_interval']))
    assert hit[0]>F('.23') and hit[1]==F(1,4)
    assert hit[0]+other[1]+F(prep['pump_failure_probability'])==1
    assert hit[1]+other[0]+F(prep['pump_failure_probability'])==1
    assert result['one_charge_dimension']==13*5456
    assert result['scope']['mass_predictions_GeV']==[]


def test_peierls_compression_against_independent_literal_exponential(result):
    mp.mp.dps=75
    data=result['spectral']['operator_bounds']
    for j,edge in enumerate((0,12)):
        variance=mp.mpf(str(model.mobile.exact_peierls_data(edge)[2].evalf(80)))
        A=mp.zeros(4)
        for k in range(3):A[k,k+1]=A[k+1,k]=mp.sqrt(variance*(k+1))
        W=mp.expm(1j*A)
        for value,interval in ((W[0,0].real,data['radial_and_surface_vacuum_overlaps'][j]),
                               (W[2,0].real,data['radial_and_surface_two_photon_amplitudes'][j])):
            assert F(interval[0])<F(str(value))<F(interval[1])
        assert abs(sum(abs(W[k,0])**2 for k in range(4))-1)<mp.mpf('1e-70')


def test_inertia_refuses_an_unresolved_pivot():
    iv=mp.iv
    A=[[iv.mpf(1),iv.mpf(0)],[iv.mpf(0),iv.mpf(2)]]
    assert model.inertia(A,F(3,2))==1
    with pytest.raises(ArithmeticError):model.inertia(A,F(1))


def test_schur_bounds_with_a_noncommuting_finite_complement():
    # Independent literal finite matrices exercise the full-complement bound,
    # not only the thirteen-dimensional source compression.
    A=np.array([[1.2,.1],[.1,1.9]])
    C=np.array([[4.,.3],[.3,5.]])
    R=np.array([[.3,-.2],[.1,.4]])
    H=np.block([[A,R.T],[R,C]])
    values=np.linalg.eigvalsh(H);a=np.linalg.eigvalsh(A)
    c=F(str(np.linalg.eigvalsh(C)[0]-1e-12))
    r=F(str(np.linalg.norm(R,2)+1e-12))
    for j in range(2):
        low=float(model.schur_lower(F(str(a[j]-1e-12)),c,r))
        assert low<=values[j]<=a[j]+1e-12
    assert values[2]>=float(c)
