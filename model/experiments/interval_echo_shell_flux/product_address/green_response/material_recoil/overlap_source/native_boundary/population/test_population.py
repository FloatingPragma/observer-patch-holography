from fractions import Fraction

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.population import (
    measure,
)


def test_integer_operator_and_saturating_dilute_inputs():
    ops = measure.certificate(8)
    vacuum = np.eye(9)[0]
    symmetric = (np.eye(9)[1] + np.eye(9)[3]) / np.sqrt(2)
    for p in (1.0, 0.25, 1 / 1024):
        state = (1 - p) * np.outer(vacuum, vacuum) + p * np.outer(symmetric, symmetric)
        incident, _, result = measure.compare(state, ops)
        np.testing.assert_allclose(result["coherent"]["gain"] / incident, 1 - p / 4)
        np.testing.assert_allclose(
            result["sharp"]["gain"] / incident, (17 / 32) * (1 - p / 2)
        )
    # Removing correlations twice on already independent factors changes nothing.
    joint = ops[0] @ state @ ops[0].T / 2
    independent = measure.scramble(joint)
    np.testing.assert_allclose(measure.scramble(independent), independent)


def test_gate_enumeration_and_witness():
    result = measure.geometry()
    assert result["words_checked"] == 3**12
    witness = result["witness"]
    assert Fraction(
        witness["degree"] + witness["internal"], witness["gate_square"]
    ) == Fraction(3, 2)
