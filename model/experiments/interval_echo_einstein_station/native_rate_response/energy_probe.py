"""A supplied dynamical energy probe and two controls on the native four-word source."""

import numpy as np
from scipy.linalg import eigh

from model.experiments.interval_echo_einstein_station import backreaction as seam
from .measure import integer_native_parts


def source_states():
    law = seam.law
    legs = np.eye(law.N_SRC, dtype=int)
    blank = np.eye(law.N_LAB + 1, dtype=int)[0]
    vacuum, left, right = [
        legs[law.SRC.index(pair)] for pair in ((0, 0), (1, 0), (0, 1))
    ]
    ordinary = [np.kron(blank, load) for load in (vacuum, left + right, left - right)]
    protected = []
    for mag, parity, phase in ((0, 1, -1), (0, 1, 1), (2, -1, -1), (2, 1, -1)):
        memory = np.r_[0, seam.coin_basis(mag, parity) @ [1, phase]]
        protected.append(np.kron(memory, vacuum))
    return np.column_stack(ordinary + protected)


def measured_probe():
    exchange, contrast = integer_native_parts()
    native = exchange + int(seam.KAPPA) * contrast
    states = source_states()
    exchange_values = np.array([0, 0, 1, 0, 0, 1, 0])
    contrast_values = np.array([0, 1, 1, 0, 0, 4, 4])
    np.testing.assert_array_equal(exchange @ states, 8 * states * exchange_values)
    np.testing.assert_array_equal(contrast @ states, 8 * states * contrast_values)
    high, low = np.array([[1, 3, 3, 1], [3, 1, 1, 3]]) / 8
    for axis in (0, 1):
        np.testing.assert_array_equal(
            high.reshape(2, 2).sum(axis), low.reshape(2, 2).sum(axis)
        )

    frequency, coupling = 1.0, 1 / 8
    qx, qz = np.array([[0, 1], [1, 0]]), np.diag([1, -1])
    state_columns = np.kron(
        np.array([[1.0], [1.0]]) / np.sqrt(2), states / np.linalg.norm(states, axis=0)
    )
    records = {}
    for name, interaction in (
        ("whole", native),
        ("contrast", int(seam.KAPPA) * contrast),
        ("idle", np.zeros_like(native)),
    ):
        source = np.sum(states * (interaction @ states), axis=0) / (
            8 * np.sum(states**2, axis=0)
        )
        total = (
            np.kron(np.eye(2), native / 8)
            + frequency / 2 * np.kron(qx, np.eye(native.shape[0]))
            + coupling * np.kron(qz, interaction / 8)
        )
        eigenvalues, basis = eigh(total)
        records[name] = []
        for duration in (1 / 64, 1 / 8, 1 / 2):
            final = basis @ (
                np.exp(-1j * duration * eigenvalues)[:, None]
                * (basis.T @ state_columns)
            )
            upper, lower = np.vsplit(final, 2)
            yes = np.sum(abs((upper - 1j * lower) / np.sqrt(2)) ** 2, axis=0)
            angular = np.hypot(frequency, 2 * coupling * source)
            np.testing.assert_allclose(
                yes,
                0.5 + coupling * source * np.sin(duration * angular) / angular,
                atol=2e-13,
                rtol=0,
            )
            records[name].append(
                dict(
                    duration=duration,
                    ordinary_yes=yes[:3].tolist(),
                    high_yes=float(high @ yes[3:]),
                    low_yes=float(low @ yes[3:]),
                )
            )
    return records


if __name__ == "__main__":
    for law, rows in measured_probe().items():
        print(law, [row["high_yes"] - row["low_yes"] for row in rows])
