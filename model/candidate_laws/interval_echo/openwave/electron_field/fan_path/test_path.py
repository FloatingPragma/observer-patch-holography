"""Check distinct terminal, remote and energy-filter instrument records."""

from fractions import Fraction

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
    FanPath,
    archived_fixture,
)
from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.measure import (
    terminal_law,
)
from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.spectral import (
    verified_inertia,
)


def test_remote_record_survives_tracing_and_changes_the_accepted_channel():
    with archived_fixture() as archived:
        path = FanPath(archived, 3, 1, asymmetric=True)
    end = path.blocks[2, 0]
    first_vacuum = np.ones(end.shape[0]) / np.sqrt(end.shape[0])
    first_excited = path.code[:, 0]
    middle_vacuum = np.ones(end.shape[1]) / np.sqrt(end.shape[1])
    columns = np.zeros((len(path.basis), 2))
    for i, far in enumerate((first_vacuum, first_excited)):
        columns[end.rows, i] = np.kron(np.kron(far, middle_vacuum), path.code[:, i])
    bank = terminal_law(path, columns)
    np.testing.assert_allclose(bank["joint_decoder_probability_range"], [1, 1])
    np.testing.assert_allclose(bank["conditional_reference_negativity"], 0, atol=1e-14)
    quiet = bank["all_other_fans_vacuum_and_decoded"]
    np.testing.assert_allclose(quiet["probability_range"], [0, 1], atol=1e-14)
    np.testing.assert_allclose(
        quiet["receiver_accepted_but_spectator_vacuum_refused_per_input"],
        [0, 1],
        atol=1e-14,
    )
    np.testing.assert_allclose(
        bank["remote_positive_ring_and_decoded"][-1][
            "conditional_on_received_charge_and_code"
        ],
        0.5,
    )


def test_dyadic_shift_applied_before_inertia_count():
    mat = sparse.csr_matrix([[4.5, -0.5], [-0.5, 4.5]])
    assert verified_inertia(mat, Fraction(15, 4))["count_below"] == 0
    assert verified_inertia(mat, Fraction(9, 2) - Fraction(1, 16))["count_below"] == 1
    assert verified_inertia(mat, Fraction(81, 16))["count_below"] == 2
