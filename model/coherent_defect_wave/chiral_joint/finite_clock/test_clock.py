"""Physical completeness, approximation limits and hostile clock controls."""
import pytest
import sympy as s
from . import model as m


def test_exact_error_identity_and_refinement():
    result = m.polynomial_certificate()
    assert result['exact_unitary_on_every_finite_self_adjoint_H']
    assert m.clock_error(s.Rational(-1,10),3,2) == m.clock_error(s.Rational(1,10),3)/64
    assert m.clock_error(1,0) == 0
    for bad in (0,-1,True,1.5):
        with pytest.raises(ValueError): m.clock_error(1,1,bad)
    for time,norm in ((0.1,1),(1,-1),(1,s.sqrt(2))):
        with pytest.raises(ValueError): m.clock_error(time,norm)


def test_clock_preserves_rare_Peter_Weyl_event_at_amplitude_level():
    r = m.peter_weyl_transfer()
    t, B = s.Rational(r['source_time']),s.Rational(r['full_H_norm_bound'])
    leading = s.Rational(r['source_target_amplitude_leading_magnitude'])
    error = s.Rational(r['finite_R33_operator_error'])
    lower,upper = map(s.Rational,r['exclusive_arrival_probability_interval'])
    assert 0 < error < leading/4 and 0 < lower < upper < 1
    # A low-order clock cannot be substituted into this proof: its available
    # error bound overwhelms the fifth-order transition.
    assert (B*t)**3/12 > leading
    # Event total-variation error is too coarse for this rare-event certificate.
    assert error > upper
    assert lower == leading**2/16 and upper == 49*lower


def test_actual_complete_ring_clock_and_independent_numerical_exponential():
    import numpy as np
    from scipy.linalg import eigh
    r = m.rational_ring_records()
    rows = r['outcomes']
    assert sum(s.Rational(row['probability']) for row in rows) == 1
    assert [row['effect_rank'] for row in rows] == [1,1,70]
    assert all(0 < s.Rational(row['probability']) < 1 for row in rows)
    block = m.physical_split()
    g = np.array(block['g'],dtype=float)
    k = np.array(block['g']*block['h'],dtype=float)
    e,v = eigh(k,g)
    initial = np.zeros(72); initial[block['p'][0]]=1
    exact = v @ (np.exp(-1j*float(s.Rational(r['source_time']))*e)*(v.T@g@initial))
    error = float(s.Rational(r['exponential_comparison_error']))
    for row,j in zip(rows[:2],(2,4)):
        target = np.zeros(72)
        target[block['p'][j]]=float(block['normalization'][j,j])
        amplitude = complex(float(s.Rational(row['amplitude_real'])),
                            float(s.Rational(row['amplitude_imaginary'])))
        assert abs(amplitude-target@g@exact) < error+1e-13
    assert not r['Peter_Weyl_spectrum_transfer']


def test_common_family_clock_retains_preparation_off_control():
    r = m.family_transfer()
    time = s.Rational(r['source_time'])
    prepared = s.Rational(r['prepared_local_arrival_probability_interval'][0])
    background = s.Rational(r['preparation_off_local_arrival_probability_interval'][1])
    assert prepared-background == time*time/48 > 0
    assert s.Rational(r['full_clock_error']) < time/96
    assert m.clock_error(time,1896,1) > time/96


def test_family_clock_rejects_changed_source_descriptor(monkeypatch):
    original = m.json.loads
    def tampered(text):
        receipt=original(text)
        if 'source_and_static_reference' in receipt:
            receipt['source_and_static_reference']['source']={'different':'source'}
        return receipt
    monkeypatch.setattr(m.json,'loads',tampered)
    with pytest.raises(ArithmeticError,match='different source descriptors'):
        m.family_transfer()


def test_custody_rejects_changed_consumed_proof(tmp_path,monkeypatch):
    import hashlib
    from ..higgs_composites import certificate as custody
    name='source.md'
    (tmp_path/name).write_text('original')
    pins={name:hashlib.sha256(b'original').hexdigest()}
    monkeypatch.setattr(custody,'ROOT',tmp_path)
    custody.validate_pins(pins)
    (tmp_path/name).write_text('different')
    with pytest.raises(ValueError): custody.validate_pins(pins)
