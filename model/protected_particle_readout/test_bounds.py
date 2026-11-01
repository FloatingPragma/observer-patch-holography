from fractions import Fraction as F

import numpy as np

from . import bounds, experiment


def test_exact_native_source_and_full_Fock_field_error_budget():
    value = bounds.source_certificate()
    assert value['full_Fock_Strang_error'] == '629/1048576'
    assert all(F(x) < F(601, 10**6) for x in value['probability_error_bounds'])
    assert value['floating_block_residuals_used_as_theorem_premises'] is False


def test_exact_characteristic_intervals_enclose_independent_full_N6_evolution():
    value = bounds.characteristic_interval(F(1))
    for index, phase in enumerate((F(), F(1, 2))):
        actual = experiment.ideal_analyzed(F(1), phase)
        probability = float(np.linalg.norm(actual[:65536])**2)
        low, high = map(F, value['real' if index == 0 else 'imag'])
        endpoints = ((1+float(low))/2, (1+float(high))/2) if index == 0 else ((1-float(high))/2, (1-float(low))/2)
        assert endpoints[0]-3e-14 <= probability <= endpoints[1]+3e-14


def test_common_scalar_energy_offset_cannot_be_measured_by_this_ramsey():
    for phase in (F(), F(1, 2)):
        original = experiment.ideal_analyzed(F(1), phase)
        shifted = experiment.ideal_analyzed(F(1), phase, scalar_shift=137.5)
        np.testing.assert_allclose(np.abs(original)**2, np.abs(shifted)**2, atol=1e-16, rtol=0)
