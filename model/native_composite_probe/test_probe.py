"""Complete-domain detector checks and controls for coherent CAR signs."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import circuits, detection, packet


@pytest.fixture(scope='module')
def retained():
    return packet.read_pair(packet.OUTPUT)


def test_integer_selector_is_full_cubic_car_operator():
    probe = detection.exact_selector_identity()
    assert probe.dtype == np.int64 and probe.shape == (256, 256)
    assert np.count_nonzero(probe) == 96
    # Removing the source occupation sign changes the actual odd CAR probe.
    wrong = (1-2*(np.arange(256) & 1))[:, None]*probe
    assert not np.array_equal(probe, wrong)


def test_coherent_preparation_and_success_kraus():
    p = detection.selector_preparation()
    assert np.allclose(p.conj().T @ p, np.eye(32), atol=1e-15, rtol=0)
    assert np.allclose(abs(p[:24, 0])**2, 1/24, atol=1e-16, rtol=0)
    assert np.array_equal(p[24:, 0], np.zeros(8))
    rng = np.random.default_rng(931)
    field = rng.normal(size=65536)+1j*rng.normal(size=65536)
    field /= np.linalg.norm(field)
    output = detection.ideal_output(field)
    expected = detection.exact_selector_identity() @ detection.as_reference_columns(field)/3
    assert np.allclose(output[0], expected, atol=2e-17, rtol=0)
    assert abs(np.linalg.norm(output)-1) < 1e-14
    # A classical uniform choice of unitary words has unit total weight;
    # it cannot replace the interference-selected cubic filter.
    assert np.linalg.norm(expected)**2 < .1


def test_literal_complete_control_bound():
    certificate = circuits.certificate()
    assert certificate['native_pulses'] == 1008
    assert F(certificate['complete_operator_and_untouched_reference_error']) < F(16, 10**14)
    assert certificate['declared_native_pointer_outcomes'] == 6**5


def test_retained_complete_probability_and_spectral_attachment(retained):
    value, arrays = retained
    run = value['execution']
    lo, hi = map(F, run['source_success_probability_interval'])
    error = F(run['complete_state_and_probability_error_upper'])
    assert float(lo-error) <= run['success_pointer_probability'] <= float(hi+error)
    assert abs(arrays['native_pointer_probabilities'].sum()+
               arrays['first_export_failure_weights'].sum()+
               run['inherited_tree_abort_probability']-1) < 2e-11
    assert arrays['complete_pointer_field_reference_amplitudes'].shape == (32, 256, 256)


def test_all_six_level_pointer_records_including_zeros():
    logical = np.arange(1, 33, dtype=float)/528
    raw = detection.native_pointer_probabilities(logical)
    assert len(raw) == 7776 and np.count_nonzero(raw) == 32
    for row, label in enumerate(np.ndindex(*(6,)*5)):
        if all(v < 2 for v in label):
            assert raw[row] == logical[sum(v << j for j, v in enumerate(label))]
        else:
            assert raw[row] == 0


@pytest.mark.parametrize('mutation', ['truncate', 'probability', 'raw_outcome', 'scope', 'hash'])
def test_receipt_rejects_incomplete_or_promoted_evidence(retained, mutation):
    value, arrays = deepcopy(retained[0]), dict(retained[1])
    if mutation == 'truncate':
        arrays['complete_pointer_field_reference_amplitudes'] = arrays['complete_pointer_field_reference_amplitudes'][:, :, :255]
    elif mutation == 'probability':
        arrays['pointer_probabilities'] = arrays['pointer_probabilities'].copy()
        arrays['pointer_probabilities'][0] += .01
    elif mutation == 'raw_outcome':
        arrays['native_pointer_probabilities'] = arrays['native_pointer_probabilities'].copy()
        arrays['native_pointer_probabilities'][2] = .01
    elif mutation == 'scope':
        value['execution']['native_energy_clock_or_resolved_138_lines_executed'] = True
    else:
        first = next(iter(value['source_bindings']))
        value['source_bindings'][first] = '0'*64
    with pytest.raises(ValueError):
        packet.validate(value, arrays)


@pytest.mark.parametrize('value', [np.zeros(256), np.zeros((65536, 2)), np.full(65536, np.nan)])
def test_full_field_required(value):
    with pytest.raises(ValueError):
        detection.as_reference_columns(value)


def test_entangled_input_native_and_abort_reference():
    field = np.zeros(65536, complex)
    # Opposite local and reference occupations test retained off-diagonal
    # reference coherences, including complex phase.
    field[6], field[6+(1 << 6)+(1 << 0)] = 1/np.sqrt(2), 1j/np.sqrt(2)
    run, arrays = detection.complete_native(field)
    assert np.linalg.norm(arrays['complete_pointer_field_reference_amplitudes']-
                          detection.ideal_output(field)) < 2e-12
    failure = detection.failure_reference_amplitudes(field, 0, 0, 0)
    assert failure.shape[1] == 256
    expected = arrays['first_export_failure_weights'][0, 0]
    assert abs(np.linalg.norm(failure)**2-expected) < 1e-38
    assert run['success_postselection_or_renormalization'] is False
    with pytest.raises(ValueError, match='padding'):
        detection.failure_reference_amplitudes(field, 0, 0, 99)
