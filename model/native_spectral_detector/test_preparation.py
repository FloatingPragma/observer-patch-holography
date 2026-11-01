from fractions import Fraction as F
import numpy as np
import pytest

from model.ground_spectral_source import preparation as old
from model.protected_particle_readout import execute
from . import preparation as p,experiment as x


def test_exact_finite_twirl_and_all_selector_outcomes():
    proof=p.quadrature_certificate()
    assert proof['all_five_diagonal_ground_weights']==['1/5']*5
    assert sum(p.SELECTOR_WEIGHTS)==90 and p.SELECTOR_WEIGHTS[-1]==0
    run=p.selector_run()
    assert np.allclose(abs(run.success[:,0])**2,np.asarray(p.SELECTOR_WEIGHTS)/90,atol=2e-14)
    assert run.failure_weights.shape==(147,4)
    assert F(p.selector_certificate()['complete_selector_error'])<F(1,10**11)


def test_every_source_rotation_including_vacuum_and_doublon():
    for record in range(15):
        instrument=execute.compile_block(execute.key(p.rotation_word(record)),2)
        expected=np.eye(4,dtype=complex)
        azimuth,polar=divmod(record,3);v=(-np.sqrt(3/5),0.,np.sqrt(3/5))[polar]
        co,si=np.sqrt((1+v)/2),np.sqrt((1-v)/2);phi=2*np.pi*azimuth/5
        expected[1:3,1:3]=np.diag([np.exp(-.5j*phi),np.exp(.5j*phi)])@np.array([[co,-si],[si,co]])
        assert np.linalg.norm(instrument.success-expected)<2e-13
        assert F(p.rotation_certificate(record)['all_eight_sites_error'])<F(1,10**10)
    with pytest.raises(ValueError):p.rotation_word(15)


def test_executed_full_ground_mixture_without_projection():
    result=x.preparation_diagnostics();columns=x.prepared()[3]
    eigenvalues=np.linalg.eigvalsh(columns.conj().T@columns)
    assert np.allclose(eigenvalues[-5:],np.full(5,.2),atol=2e-12)
    assert np.max(abs(eigenvalues[:-5]))<2e-12
    assert result['all15_conditional_ground_words_executed']
    assert result['complete_preparation_error_diagnostic']<1e-10
    assert not result['native_ground_rank_projection_used']
    # A single highest-weight preparation is a different density operator.
    highest=old.ideal_ground()
    overlap=np.linalg.norm(highest.conj()@columns)**2
    assert abs(overlap-.2)<2e-12
