from fractions import Fraction as F
import numpy as np
import pytest

from model.protected_particle_readout import execute
from model.common_source.chiral_candidate import gaussian_source as gaussian
from . import circuits as c,execution as e,bounds,experiment as x


def test_exact_charge_instrument_sign_and_complete_branches():
    instrument=c.charge_blocks()[0][0]
    for mask in range(4):
        actual=instrument.success[:,mask]
        expected=np.zeros(8,complex)
        target=(mask^2)+4*((mask>>1)&1)
        expected[target]=(-1)**(mask&1)
        assert np.linalg.norm(actual-expected)<1e-13
    # Omitting the parity sign changes the occupied-up amplitude.
    assert instrument.success[3,1].real<-.99
    assert np.linalg.norm(instrument.success.conj().T@instrument.success+
                          instrument.effects.sum(axis=(0,1))-np.eye(8))<1e-12


def test_controlled_hopping_and_contact_on_all_columns():
    for z in (1/16,-1/16,1j/16,-1j/16):
        word=c.emit(c.controlled_hopping_template(z,F(1,7)))
        actual=execute.compile_block(execute.key(word),3).success
        target=np.eye(8,dtype=complex)
        target[4:,4:]=gaussian.local_gate('hopping',z,1/7).fock_unitary
        assert np.linalg.norm(actual-target)<1e-12
    word=c.emit(c.controlled_contact_template(F(1,7)))
    actual=execute.compile_block(execute.key(word),3).success
    target=np.eye(8,dtype=complex);target[7,7]=np.exp(-1j/7)
    assert np.linalg.norm(actual-target)<1e-12


def test_positive_fourier_all_clock_columns():
    for bits in (1,2,3):
        blocks=tuple((instrument,tuple(j-16 for j in targets)) for instrument,targets in c.fourier_blocks(bits))
        d=1<<bits;actual=e.complete(blocks,bits,np.eye(d)).success
        target=np.exp(2j*np.pi*np.arange(d)[:,None]*np.arange(d)[None,:]/d)/np.sqrt(d)
        assert np.linalg.norm(actual-target)<2e-13


def test_complete_probability_and_resource_certificate():
    proof=bounds.certificate()
    assert F(proof['complete_state_and_total_variation_error'])<F(3,100)
    assert proof['largest_coherent_history_time']=='3'
    assert proof['pulse_counts']['worst_case_complete_history']==612666
    rows=bounds.ideal_probabilities()
    for tag,mass in ((0,F(3,4)),(1,F(1,4))):
        selected=[r for r in rows if r['charge_record']==tag]
        assert sum(F(r['ideal_joint_probability_interval'][0]) for r in selected)<=mass
        assert sum(F(r['ideal_joint_probability_interval'][1]) for r in selected)>=mass
    assert bounds.finite_shot_certificate()['shots_or_random_measurement_records_generated'] is False
