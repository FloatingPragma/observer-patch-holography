from fractions import Fraction as F
import math
import numpy as np
import pytest
from . import algebra, spectral, circuits, experiment, bounds
from model.native_ordered_ground import source
from model.native_spin_observer import circuits as native
from model.native_spin_observer import execute as native_execute
from model.protected_particle_readout import execute as local


def test_exact_complete_cube_group_and_character_resolution():
    h, scale = source.matrix()
    assert algebra.audit(h, int(scale))['exact_group_and_source_checks']
    vector = np.arange(256, dtype=object) - 127
    parts = [algebra.character_projection(vector, k) for k in range(8)]
    assert np.array_equal(sum(parts), 8 * vector)
    for k, part in enumerate(parts):
        for j in range(8):
            assert np.array_equal(algebra.character_projection(part, j), 8 * part if j == k else 0 * part)


@pytest.mark.parametrize('energy', [F(), F(1, 3), F(3), F(-7, 2)])
def test_fejer_complete_probability_and_cosine_intervals(energy):
    intervals = [spectral.fejer_interval(energy, j) for j in range(4)]
    assert sum(x[0] for x in intervals) <= 1 <= sum(x[1] for x in intervals)
    for j, (lo, hi) in enumerate(intervals):
        direct = abs(sum(np.exp(-1j * n * (float(energy) - 2 * math.pi * j / 4)) for n in range(4)) / 4) ** 2
        assert float(lo) - 2e-15 <= direct <= float(hi) + 2e-15


def test_actual_source_common_kernel_and_nontransverse_countercontrol():
    checked = algebra.dimer_common_kernel()
    assert checked['exchange_annihilates'] and checked['selected_field_annihilates']
    r = math.sqrt(5) - 2
    ground = np.array([r, 1, -1, -r])
    exchange = (np.eye(4)[[0, 2, 1, 3]] - np.eye(4)) / 2
    xprobe = ground[[1, 0, 3, 2]]
    assert np.linalg.norm(exchange @ xprobe) > .1


def test_controlled_reflection_literal_source_includes_relative_phase():
    template = circuits.combined_dimer_reflection(control=4, bits=5)
    events = tuple(native.emit(template))
    initial = np.zeros((32, 8), complex)
    target = np.zeros((32, 8), complex)
    for logical in range(4):
        encoded = sum(1 << (2 * j + (0 if logical >> j & 1 else 1)) for j in range(2))
        for control in (0, 1):
            column = 2 * logical + control
            initial[encoded + (control << 4), column] = 1 / math.sqrt(8)
            swapped = ((logical & 1) << 1) | ((logical >> 1) & 1)
            out = sum(1 << (2 * j + (0 if swapped >> j & 1 else 1)) for j in range(2)) if control else encoded
            target[out + (control << 4), column] = ((-1) ** logical.bit_count() if control else 1) / math.sqrt(8)
    result = native_execute.complete(events, 5, initial)
    assert np.linalg.norm(result.success - target) < 1e-11


def test_positive_fourier_reads_positive_excitation_energy():
    template = circuits.positive_fourier((0, 1))
    instrument = local.compile_block(local.key(tuple(native.emit(template))), 2)
    for j in range(4):
        state = np.exp(-2j * math.pi * j * np.arange(4) / 4) / 2
        actual = instrument.success @ state
        assert abs(actual[j] - 1) < 1e-12
        assert np.linalg.norm(np.delete(actual, j)) < 1e-12


def test_dimer_normalization_and_exact_reflection_weights():
    target, probabilities = experiment.ideal_target()
    assert np.linalg.norm(target) == pytest.approx(1)
    assert probabilities.sum() == pytest.approx(1)
    assert probabilities.sum(axis=1) == pytest.approx([(1 - 2 / math.sqrt(5)) / 2, (1 + 2 / math.sqrt(5)) / 2])


def test_operator_bound_retains_frame_and_state_errors():
    result = spectral.operator_error(4, F(1), F(1, 100), F(1, 200), F(1, 1000), F(1, 10000))
    assert result == F(5, 4) * F(3, 200) + F(1, 1000) * F(2001, 1000) + F(1, 5000)
    with pytest.raises(ValueError):
        spectral.operator_error(True, F(1), F(), F(), F(), F())
    with pytest.raises(ValueError):
        spectral.operator_error(4, F(1), F(), F(), F(1), F())


def test_conditional_bound_requires_nonzero_approximate_success():
    value = bounds.conditional_bound(F(1, 100), F(1, 10))
    assert value['approximate_success_lower'] == F(19, 200)
    assert value['conditional_unhalved_trace_error'] == F(1, 5)
    for delta, success in ((F(), F()), (F(1, 10), F(1, 10)), (F(1), F(1, 2))):
        with pytest.raises(ValueError):
            bounds.conditional_bound(delta, success)


def test_clock_leakage_and_control_duration_budgets():
    assert bounds.leakage_upper(64, 1, F(1)) == F(1, 256)
    result = bounds.finite_program_error(4, F(1, 2), F(1), F(1, 4), 12, 24, F(1, 8192), F(1, 10**12), 100000, 64)
    assert result['clock_duration_sum'] == F(3, 2)
    assert result['operator_splitting_error'] == F(81, 16384)


def test_warmed_caches_do_not_admit_numeric_type_aliases():
    spectral.cosine_interval(F(), F())
    algebra.reflections(1)
    for args in ((False, F()), (F(), 0)):
        with pytest.raises(ValueError):
            spectral.cosine_interval(*args)
    with pytest.raises(ValueError):
        algebra.reflections(True)
