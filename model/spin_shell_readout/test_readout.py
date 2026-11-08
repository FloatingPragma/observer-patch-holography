from copy import deepcopy
from fractions import Fraction as F
import math
import numpy as np
import pytest

from . import bounds as b, finite, packet


def fejer(size, phase):
    turns = phase-np.arange(size)/size
    values = abs(np.mean(np.exp(2j*np.pi*np.outer(turns, np.arange(size))), axis=1))**2
    assert abs(values.sum()-1) < 1e-12
    return values


@pytest.mark.parametrize('size', [4, 16, 64, 128])
def test_summed_circular_tail_controls_every_distant_bin(size):
    for phase in (F(), F(1, 17), F(1, 2), F(999, 1000)):
        weights = fejer(size, float(phase))
        for margin in (F(1, 100), F(1, 8), F(1, 3), F(1, 2)):
            distances = [min(abs(phase-F(j, size)), 1-abs(phase-F(j, size))) for j in range(size)]
            tail = sum(w for w, d in zip(weights, distances, strict=True) if d >= margin)
            assert tail <= float(b.tail_bound(size, margin))+1e-12


def test_clock_choice_and_huge_counter_do_not_allocate_a_source():
    plan = b.clock_plan(F(1, 10**122), F(1, 100))
    assert plan['clock_size'] & (plan['clock_size']-1) == 0
    assert plan['complete_clock_tail_upper'] <= F(1, 100)
    assert plan['controlled_binary_duration_sum_over_tau'] == plan['clock_size']-1


def test_inner_outer_bounds_hold_for_coherent_clock_laws_and_uncertain_center():
    size = 256
    centers, width, margin = (F(39, 100), F(41, 100)), F(1, 10), F(1, 50)
    plan = b.bin_ranges(size, centers, width, margin)
    for center in (centers[0], F(2, 5), centers[1]):
        phases = [F(1, 10), center-width, center, center+width, F(4, 5)]
        masses = [F(1, 20), F(1, 10), F(1, 5), F(1, 10), F(1, 20)]
        law = sum(float(w)*fejer(size, float(p)) for p, w in zip(phases, masses, strict=True))
        inner = sum(law[plan['inner'][0]:plan['inner'][1]+1])
        outer = sum(law[plan['outer'][0]:plan['outer'][1]+1])
        actual = sum(w for p, w in zip(phases, masses, strict=True) if abs(p-center) <= width)
        def enclose(x):
            return F(str(max(0, x-1e-12))), F(str(min(1, x+1e-12)))
        report = b.tube_bounds((sum(masses), sum(masses)), enclose(inner), enclose(outer),
                               plan['complete_clock_tail_upper'])
        assert report['original_tube_weight_lower'] <= actual <= report['original_tube_weight_upper']
        assert report['fraction_of_original_success_lower'] <= actual/sum(masses)


def test_absolute_error_can_destroy_a_small_success_claim():
    result = b.tube_bounds((F(1, 10**9),)*2, (F(9, 10**10),)*2,
                          (F(1, 10**9),)*2, F(1, 100), F(1, 10**8))
    assert result['original_tube_weight_lower'] == 0
    assert result['fraction_of_original_success_lower'] == 0


def test_zero_success_and_all_failures_are_retained():
    result = b.tube_bounds((F(),)*2, (F(),)*2, (F(),)*2, F(1, 10))
    assert result['original_tube_weight_upper'] == 0
    assert result['fraction_of_original_success_lower'] == 0


def test_alias_rejection_uses_full_spectrum_and_periodic_control():
    finite.energy_chart(F(8), F(1, 2))
    with pytest.raises(ValueError, match='aliases'):
        finite.energy_chart(F(8), F(1))
    # Energy shifts by a clock period leave every output unchanged.
    assert np.max(abs(fejer(64, .17)-fejer(64, 1.17))) < 1e-12


def test_scan_confidence_charges_the_predeclared_candidate_count():
    ordinary = b.sampled_intervals(750000, 10000, 8000, 9000, F(1, 100), F(1, 100))
    assert ordinary['all_other_outcome_count'] == 740000
    assert ordinary['required_trials'] == 750000
    with pytest.raises(ValueError, match='too few trials'):
        b.sampled_intervals(750000, 10000, 8000, 9000, F(1, 100), F(1, 100), searches=2)


@pytest.mark.parametrize('call', [
    lambda: b.tail_bound(True, F(1, 10)),
    lambda: b.tail_bound(16, .1),
    lambda: b.tail_bound(16, F(3, 4)),
    lambda: b.bin_ranges(16, (F(1, 20),)*2, F(1, 10), F(1, 100)),
    lambda: b.bin_ranges(16, (F(19, 20),)*2, F(1, 10), F(1, 100)),
    lambda: b.bin_ranges(16, (F(1, 2),)*2, F(1, 10), F(1, 10)),
    lambda: b.tube_bounds((F(1, 5),)*2, (F(1, 2),)*2, (F(3, 5),)*2, F(1, 10)),
    lambda: b.sampled_intervals(1000, 10, 9, 8, F(1, 10), F(1, 10)),
])
def test_invalid_clock_alias_and_incomplete_count_inputs_fail(call):
    with pytest.raises(ValueError):
        call()


@pytest.mark.parametrize('energy,j', [(F(), 0), (F(), 3), (F(1, 2), 1), (F(11, 10), 11)])
def test_fast_kernel_encloses_independent_clock_dilation(energy, j):
    lo, hi = finite.kernel_interval(energy, j, 128, F(1, 2))
    value = fejer(128, float(energy)/(4*math.pi))[j]
    assert float(lo)-1e-13 <= value <= float(hi)+1e-13


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_complete_selected_source_has_positive_unnormalized_tube_weight(result):
    row = result['finite_source_readout']
    assert F(row['original_band_weight_interval'][0]) > F(1, 5)
    assert row['all_8_times_128_joint_bins_retained'] is True
    assert row['source_spectral_cluster_indices']
    assert row['character_is_translation_momentum'] is False
    assert row['additional_native_128_bin_tape_executed'] is False
    assert row['sharp_scaling_mass_shell_or_physical_mass_proved'] is False
    assert all(not name.startswith('docs/') for name in result['source_bindings'])


@pytest.mark.parametrize('mutation', ['mass', 'weight', 'omit_failure', 'numeric_alias', 'source'])
def test_scientific_promotions_and_source_mutations_fail(result, mutation):
    damaged = deepcopy(result)
    if mutation == 'mass':
        damaged['scope']['sharp_continuum_mass_shell_or_measured_mass_derived'] = True
    elif mutation == 'weight':
        damaged['finite_source_readout']['band_weight_strict_lower'] = '1/2'
    elif mutation == 'omit_failure':
        del damaged['finite_source_readout']['all_8_times_128_joint_bins_retained']
    elif mutation == 'numeric_alias':
        damaged['scope']['same_selected_finite_spin_state_and_clock'] = 1
    else:
        damaged['source_bindings'][next(iter(damaged['source_bindings']))] = '0'*64
    with pytest.raises(ArithmeticError):
        packet.validate(damaged, result)
