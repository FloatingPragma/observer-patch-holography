"""Independent spin-space audit and checks of the archived full readout."""

from itertools import combinations
import json

import numpy as np
from scipy.linalg import expm

from .measure import HERE, HAND, SHOTS, SOURCE, TIMES, material, read


def physical_pair(n, a, b):
    # Enumerate actual dimer orientations; use the virtual end signs only in
    # computing the physical wall word. Act i gamma_a gamma_b from the right.
    states = np.arange(1 << n)
    spin = 1 - 2 * ((states[:, None] >> np.arange(n)) & 1)
    z = np.column_stack([np.ones(len(spin), int), spin, -np.ones(len(spin), int)])
    wall_z = -z[:, :-1] * z[:, 1:]
    amp = np.full(len(spin), 1j)
    for index in (b, a):
        j, odd = divmod(index, 2)
        amp *= np.prod(wall_z[:, :j], axis=1)
        if odd:
            amp *= 1j * wall_z[:, j]
        wall_z[:, j] *= -1
    dest_sign = np.cumprod(-wall_z, axis=1)[:, :n]
    dest = ((1 - dest_sign) // 2) @ (1 << np.arange(n))
    op = np.zeros((len(spin), len(spin)), complex)
    op[dest, states] = amp
    return op


def test_every_small_chain_even_pointer_is_a_physical_pauli_and_two_time_trace():
    n = 4
    states = np.arange(1 << n)
    spin = 1 - 2 * ((states[:, None] >> np.arange(n)) & 1)
    pairs = list(combinations(range(2 * (n + 1)), 2))
    op = np.stack([physical_pair(n, a, b) for a, b in pairs])
    np.testing.assert_array_equal(op, op.conj().transpose(0, 2, 1))
    np.testing.assert_array_equal(op @ op, np.broadcast_to(np.eye(1 << n), op.shape))
    # Inspect the actual n-site Pauli factorization of every W (including sign).
    for mat in op:
        dest = np.argmax(abs(mat), axis=0)
        assert np.all((states ^ dest) == dest[0])
        signs = mat[dest, states] / mat[dest[0], 0]
        for bit in range(n):
            np.testing.assert_array_equal(
                signs[states ^ (1 << bit)], signs * signs[1 << bit]
            )
    for j in (8, 4):
        h = 0.5
        ham = np.diag(
            j * (np.sum(spin[:, :-1] * spin[:, 1:], axis=1) + spin[:, 0] - spin[:, -1])
        ).astype(float)
        for k in range(n):
            ham[states, states ^ (1 << k)] -= h
        u = expm(-0.732j * ham)
        # Tr[U†W_ab U W_ik]/2^n uses no single-wall sector, JW energy,
        # dispersion, prepared wave number or forecast.
        evolved = u.conj().T @ op @ u
        traced = np.einsum("aij,bji->ab", evolved, op) / (1 << n)
        r = expm(0.732 * material(n, j, h))
        a, b = np.array(pairs).T
        expected = (
            r[a[:, None], a] * r[b[:, None], b] - r[a[:, None], b] * r[b[:, None], a]
        )
        np.testing.assert_allclose(traced, expected, atol=3e-13)


def test_projector_recovery_with_orthogonal_auxiliaries():
    rng = np.random.default_rng(89)
    r, _ = np.linalg.qr(rng.normal(size=(24, 24)))
    i, j, k = SOURCE
    mats = []
    for x, y in ((i, j), (i, k), (j, k)):
        cross = np.outer(r[:, x], r[:, y])
        m = cross - cross.T
        mats.append(m @ m.T)
    np.testing.assert_allclose(
        (mats[0] + mats[1] - mats[2]) / 2, np.outer(r[:, i], r[:, i]), atol=4e-15
    )


def test_record_coverage_waveforms_and_falsifiers():
    data = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "counts.npz") as counts:
        banks, hand = counts["banks"], counts["hand"]
        assert banks.shape == (2, 61, 3, 435, 4) and hand.shape == (2, 61, 2)
        assert np.all(banks.sum(axis=-1) == SHOTS) and np.all(hand <= HAND)
        for bank, clock, row in zip(banks, hand, data["rows"], strict=True):
            vec, times, overlap, gap = read(bank, clock)
            np.testing.assert_allclose(vec[:, ::2], row["reconstructed_A"], atol=1e-11)
            assert max(abs(times - TIMES)) < data["time_error_bound"]
            assert overlap > 0.5 and gap > 0.99
            assert row["own_law_max"] < row["kg_count_and_clock"]
            assert row["remote_A9_max"] > 10 * row["column_error_bound"]
        native, changed = data["rows"]
        assert native["native_kg_max"] < native["kg_count_and_clock"]
        assert native["hopping_max"] > 2 * native["hopping_total_allowance"]
        assert native["immobile_max"] > 1
        assert (
            changed["native_kg_max"]
            > changed["calibration_box_kg"] + changed["kg_count_and_clock"]
        )
        long_wave = native["ordinary_long_wave"]
        assert len(long_wave["k"]) == 2
        assert 0 < min(long_wave["k"]) < max(long_wave["k"]) < 0.5
        assert (
            long_wave["max_error"]
            < native["kg_count_and_clock"] + long_wave["lattice_to_continuum_bound"]
        )
