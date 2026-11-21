"""Check the single-parity obstruction and its collapsed-axis exception on the native route."""
import math

import numpy as np
import pytest
from scipy.linalg import null_space

from .. import law


@pytest.mark.parametrize("parity", (-1, 1))
def test_all_dark_edge_fields_in_a_four_dimensional_cell(parity):
    # Every edge coefficient in the 16-vertex cube is free before the
    # zero local brightness conditions. Keep both endpoint slot labels.
    sites = list(np.ndindex((2,) * 4))
    edges = []
    for start in sites:
        for axis in range(4):
            if start[axis] == 0:
                end = list(start)
                end[axis] = 1
                edges.append((start, tuple(end), axis))
    incidence = np.zeros((len(sites), len(edges)))
    for j, (start, end, _) in enumerate(edges):
        incidence[sites.index(start), j] = 1
        incidence[sites.index(end), j] = parity
    dark = null_space(incidence)
    assert dark.shape == (32, 17)
    # The Fourier component on each separate pair of slots vanishes at
    # zero momentum for F=-1, and at the all-pi momentum for F=+1.
    for coefficients in dark.T:
        field = np.zeros((4, 4, 4, 4, 8))
        for amplitude, (start, end, axis) in zip(coefficients, edges):
            field[(*start, 2*axis)] = amplitude
            field[(*end, 2*axis+1)] = parity * amplitude
        np.testing.assert_allclose(law.route(field), parity*field, atol=1e-14, rtol=0)
        np.testing.assert_allclose(field.sum(axis=-1), 0, atol=1e-14, rtol=0)
        k = ((0 if parity == -1 else 2),) * 4
        np.testing.assert_allclose(np.fft.fftn(field, axes=(0,1,2,3))[k], 0, atol=1e-14, rtol=0)
    for age in (0, 127, 128):
        theta = math.pi*(1-math.cos(math.pi*(age+1)/128))/4 if age < 128 else math.pi/2
        for r in (0, 1):
            phase = -np.exp(1j*(128+r)*math.pi/4096)
            if parity == -1:
                phase *= np.exp(-2j*theta)
            np.testing.assert_allclose(law.step_maturing(field, law.LAWS["torus"], age, r),
                                       phase*field, atol=1e-14, rtol=0)


def test_private_period_one_has_orthogonal_onsite_negative_parity_translates():
    shape = (3, 3, 3, 1)
    values = np.arange(1, 1+math.prod(shape)).reshape(shape)
    field = np.zeros((*shape, 8))
    field[..., 6] = values/math.sqrt(2)
    field[..., 7] = -field[..., 6]
    np.testing.assert_array_equal(law.route(field), -field)
    np.testing.assert_array_equal(field.sum(axis=-1), np.zeros(shape))
    # Each origin is encoded into its own pair of slots, with unit gain.
    np.testing.assert_allclose(np.sum(abs(field)**2, axis=-1), abs(values)**2, rtol=2e-15)
    negative_mature = np.exp(1j*math.pi/32)
    np.testing.assert_allclose(law.step_maturing(field, law.LAWS["torus"], 128, 0),
                               negative_mature*field, atol=3e-14, rtol=0)
