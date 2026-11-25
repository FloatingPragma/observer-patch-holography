from fractions import Fraction
from itertools import product
from math import log, sqrt

import numpy as np
import pytest

from . import measure


@pytest.fixture(scope="module")
def acquisition():
    return measure.run()


def test_independent_and_persistent_tags_on_complete_native_record(acquisition):
    np.testing.assert_allclose(acquisition["native_forward"], [7/12, 1/2, 5/12], rtol=0, atol=1e-13)
    np.testing.assert_allclose(acquisition["no_reversal_forward"], [1/3]*3, rtol=0, atol=1e-13)
    assert acquisition["energy_multiplicities_per_36"] == {0: 4, 32: 8, 33: 8, 64: 4, 65: 8, 66: 4}
    two = acquisition["two_writes"]
    assert two["two_write_trace_distance"] == "1/72"
    assert two["forward_covariance"] == "1/144"
    assert abs(two["signed_modular_slopes"]["0.5"] + log((36.5*4.5) / 11.5**2)/144) < 1e-15
    assert two["signed_modular_slopes"]["-1.0"] > 0
    assert abs(two["signed_modular_slopes"]["0.0"]) < 1e-15
    assert two["folded_amplitude_slope"] == pytest.approx(log(3)/144)


def test_finite_convolution_against_all_binary_words():
    n, k = 2, 3
    expected = Fraction(0)
    for tags in product((Fraction(5, 12), Fraction(7, 12)), repeat=n):
        for clicks in product((0, 1), repeat=n*k):
            probability = Fraction(1, 2**n)
            for i, value in enumerate(clicks):
                p = tags[i//k]
                probability *= p if value else 1-p
            expected += probability * abs(Fraction(sum(clicks), k) - Fraction(n, 2))
    folded = measure.finite_aggregate_mean(n, k, 5/12, 7/12)
    assert abs(folded - log(3)*float(expected)) < 1e-14
    assert measure.mean_absolute_signs(1) == 1
    assert measure.mean_absolute_signs(3) == Fraction(3, 2)


def test_bank_controls_and_noise_limit(acquisition):
    for row in acquisition["comparisons"]:
        n, k = row["sources"], row["reads_per_source"]
        if k == 1:
            assert row["paired"] == pytest.approx(row["fresh"], abs=1e-11)
        assert row["fresh"] == pytest.approx(float(measure.mean_absolute_signs(n*k)/(2*k))*log(3), abs=4e-11)
        assert row["no_reversal_exchange"] == pytest.approx(row["no_reversal_fresh"], abs=4e-11)
        assert row["paired"] - row["infinite_reads"] <= log(3)*sqrt(35*n/k)/12 + 1e-11


def test_restricted_material_modular_record():
    from . import material_algebra

    result = material_algebra.run()
    assert result["hidden_material_slope"] == 0 and result["hidden_hamiltonian_energy"] == "1/4"
    for row in result["visible_returns"]:
        if row["read_contrast"] == 0:
            np.testing.assert_allclose(row["formed"], row["independent"], atol=1e-13)
        elif row["exposure"] == 2.:
            assert row["modular_slope"] == pytest.approx(-0.0645423913551)
            assert row["independent_reference_subtracted"] == pytest.approx(-0.0058383417769)
