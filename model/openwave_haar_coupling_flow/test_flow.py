from fractions import Fraction as F
from .flow import run,quotient
from .certificate import calculate
from model.openwave_coupled_haar_spectrum.certificate import calculate as parent
import pytest


def test_complete_family_and_positive_derivative():
    r=run();assert len(r['records'])==7 and len(r['finite_coupling_steps'])==6
    for row in r['records']:
        assert F(row['gap_derivative_interval_exact'][0])>0
        assert not row['slope_is_a_renormalization_group_beta_function']
    assert not r['spatial_refinement_performed']


def test_reused_certificate_has_identical_source_results():
    for key in ('shared_g2_2_cap8','shared_g2_4_cap6','shared_g2_8_cap6'):
        a,b=calculate(key)[0],parent(key)[0]
        # Independently generated fixed coordinates may differ in irrelevant
        # rounding; all complete full-source intervals must overlap.
        for field in ('ground_energy_interval_exact','full_Haar_gap_interval_exact',
                      'connected_Wilson_covariance_interval_exact'):
            x,y=tuple(map(F,a[field])),tuple(map(F,b[field]))
            assert max(x[0],y[0])<=min(x[1],y[1])


def test_interval_ratio_rejects_unproved_positive_denominator():
    with pytest.raises(AssertionError):quotient((F(1),F(2)),(F(0),F(1)))
    assert quotient((F(2),F(4)),(F(1),F(2)))==(F(1),F(4))
