from fractions import Fraction as F
import numpy as np
import pytest
from .operators import polynomial_member
from .certificates import complete_eigensystem,actual_vector_certificate
from .records import execute


def test_complete_instrument_and_work_retain_field_vacuum():
    r=execute(12,-3,F(1),(2,))
    x=r['records'][0]
    assert abs(sum(x['conditional_position_flux_probabilities'])-1)<1e-12
    assert abs(sum(x['native_probe_probabilities_96'])-1)<1e-12
    assert x['preparation_failure_probability_exact']=='1/2'
    assert abs(x['complete_work_balance_residual_diagnostic'])<3e-12
    assert x['initial_and_failure_field_energy_exact']=='1/8'
    assert len(x['complete_position_flux_labels'])==100


def test_exact_complete_spin_intertwiner_and_charge_reference():
    for q in(-3,3):
        m=polynomial_member(5,q,1);n=len(m['basis'])//2;H=m['H_minus_mass']
        assert np.array_equal(H[n:,n:],H[:n,:n].conj())
        assert np.count_nonzero(H[:n,n:])==0
        assert m['relative_dipole_charge_unit']=='native third-charge q=3Q'


def test_actual_amplitude_corruption_is_rejected():
    m=polynomial_member(5,-3,1);n=len(m['basis'])//2
    c=complete_eigensystem(m['H_minus_mass'][:n,:n],1,6)
    i=m['index'][(0,0,0)];Q=c['Q'];psi=Q@(np.exp(-2j*c['eigenvalues'])*Q[i].conj())
    actual_vector_certificate(c,psi,i,2)
    psi[0]+=.01
    with pytest.raises(AssertionError):actual_vector_certificate(c,psi,i,2)
