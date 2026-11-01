from fractions import Fraction as F

import numpy as np
import pytest

from model.interacting_vacuum import ground, charged
from . import bounds, circuits as c


def test_complete_source_bound_counts_every_segment():
    packet = bounds.certificate()
    assert max(map(F, packet['complete_native_state_and_probability_errors'])) < F(3, 5000)
    assert F(packet['full_Fock_Strang_error']) > 0
    assert packet['source_pulse_counts']['preparation'] == 5563
    assert packet['source_pulse_counts']['field_per_segment'] > 100000
    assert packet['all_native_export_and_repair_outcomes_retained']
    assert not packet['hardware_or_host_floating_accuracy_proved']


def test_exact_positive_gap_sum_and_quadrature_separation():
    rows = bounds.spectral_rows()
    assert len(rows) == 78
    assert all(F(row['gap_above_ground_interval'][0]) > 0 for row in rows)
    actual = bounds.probability_intervals()
    assert F(actual[1]['literal_complete_pointer_zero_interval'][0]) > F(3, 5)
    atzero = bounds.characteristic_interval(F())
    assert F(atzero['real'][0]) <= 1 <= F(atzero['real'][1])
    assert F(atzero['imag'][0]) <= 0 <= F(atzero['imag'][1])


def test_source_free_ground_probe_has_real_cosine_control():
    a, *_ = ground.band_data()
    assert np.array_equal(a@a, 3*np.eye(8, dtype=np.int64))
    assert a[0, 0] == 0
    # For U0 the same G is a kinetic eigenstate, and the added empty down
    # orbital has h=A/16. Thus C0(t)=cos(sqrt3*t/16), with zero imaginary part.
    moments = charged.spectral_certificate()['down_addition']
    assert moments['gap_mean'] == '1/2' and moments['gap_second_moment'] == '131/256'


@pytest.mark.parametrize('time', [1., -F(1), F(3), True])
def test_characteristic_domain_rejected(time):
    with pytest.raises(ValueError): bounds.characteristic_interval(time)


def test_all_literal_tape_templates_are_consumed():
    word = c.pointer_word(); template = c.pointer_template()
    assert bounds.local(word, template) < F(1, 10**11)
    with pytest.raises(ArithmeticError): bounds.local(word+word[:1], template)
