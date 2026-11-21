import numpy as np
import pytest

from . import formed_clock_readout as clock


@pytest.mark.parametrize("size", clock.SIZES)
def test_ordinary_source_and_competing_birth(size):
    for f in clock.FAR_POWER:
        for competing in (False, True):
            refused, state = clock.ordinary_birth(size, f, competing=competing)
            expected = (1-f)**2/2
            other = f*f/2 if competing else 0.
            np.testing.assert_allclose(refused, [0, other, 1-expected-other], atol=1e-14)
            np.testing.assert_allclose(np.vdot(state,state).real, expected, atol=1e-14)
            np.testing.assert_allclose(state[0,0,0,0], np.sqrt(expected)*clock.EQUAL, atol=1e-14)
        refused, state = clock.ordinary_birth(size, f, reverse=True)
        np.testing.assert_allclose(refused, [0,0,1], atol=1e-14)
        assert not np.any(state)


@pytest.mark.parametrize("size", clock.SIZES)
def test_spectral_mass_and_spatial_wash(size):
    e, a, lines = clock.poles(size)
    times = clock.AGES-clock.OPENING
    frequencies = np.array([e[1]-e[0], sum(e), 2*e[0], 2*e[1]])
    assert len(set(np.round(frequencies, 9))) == 4
    reconstructed = (abs(a)**2).sum() + 2*(np.exp(-1j*np.outer(times, frequencies))@lines).real
    np.testing.assert_allclose(reconstructed, clock.predicted_intensity(size), atol=1e-13)
    assert 0 < abs(a[0,1]*a[1,1]) < 1e-6*abs(a[0,0]*a[1,0])
    q = (1-.25)**2/2
    native = clock.native_words(size, .25, "ordinary")
    np.testing.assert_allclose(native[:,-1], q*reconstructed, atol=3e-13)
    np.testing.assert_allclose(native, clock.native_words(size,.25,"reflected"), atol=1e-13)
    washed = clock.native_words(size,.25,"translated")
    assert np.ptp(washed[:,-1])/q < .00025
    # A birth from another wire has its own refusal; it cannot supply this
    # herald's later point clicks even if its pair reaches the same receiver.
    competitor = clock.native_words(size,.25,"competing_birth")
    assert np.all(competitor[:,1] > .03)
    np.testing.assert_allclose(competitor[:,3:], native[:,3:], atol=1e-13)


def test_alias_reader_from_original_start_intensity():
    t = np.arange(128)
    counts = np.zeros((128, 6), np.int64)
    counts[:, -1] = np.rint(clock.SHOTS*(.012+.006*np.cos(.4*t+.6)))
    frequency, prediction = clock.prefix_fit(counts)
    assert abs(frequency-.4) < 5e-5
    assert np.max(abs(prediction-counts[:,-1]/clock.SHOTS)) < 1e-7
    bands = clock.frequency_intervals(counts, 0., .5)
    assert len(bands) == 1 and bands[0][0] <= .4 <= bands[0][1]
    counts[:, -1] = 0
    assert clock.frequency_intervals(counts, 0., .5) == [[0., np.pi]]


def test_empirical_spread_requires_both_ordered_controls():
    ordinary = np.zeros((128, 6), np.int64)
    ordinary[1::2,-1] = clock.SHOTS//100
    washed = np.zeros_like(ordinary)
    assert clock.empirical_wash(ordinary, washed, .5)["positive_lower_bound"]
    assert not clock.empirical_wash(washed, ordinary, .5)["positive_lower_bound"]
    assert not clock.empirical_wash(ordinary, ordinary, .5)["positive_lower_bound"]
    with pytest.raises(ValueError, match="both complete"):
        clock.empirical_wash(ordinary, washed[:-1], .5)
