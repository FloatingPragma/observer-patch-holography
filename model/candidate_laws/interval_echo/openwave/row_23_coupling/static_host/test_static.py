"""Independent microscopic operators and exact static-source controls."""

from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation

from . import measure as m


def physical_seam(mix):
    exchange, contrast = law.local_term_parts()
    h = (exchange + candidate.KAPPA * contrast)[:, :9].copy()
    assert np.max(np.abs(h[9:, :])) < 1e-13
    h = h[:9]
    # Independent index order for the two-load, two-port channel.
    mix_state = np.zeros(9)
    mix_state[[2, 6]] = 1 / np.sqrt(2)
    middle = np.eye(9)[4]
    return h - mix * (np.outer(mix_state, middle) + np.outer(middle, mix_state))


def test_reduced_columns_and_held_energy_clock_from_full_four_site_operator():
    for setting, mix in (("native", 0), ("conversion", 8), ("cut_conversion", 8)):
        h = physical_seam(mix)
        full = np.kron(h, np.eye(9)) + np.kron(np.eye(9), h)
        if setting != "cut_conversion":
            full += np.kron(np.eye(3), np.kron(h, np.eye(3)))
        for left, right in product(range(3), repeat=2):
            fixed = [27 * left + 3 * j + right for j in range(9)]
            block = full[np.ix_(fixed, fixed)]
            selected = [4, 6, 2]
            assert (
                np.max(
                    np.abs(
                        block[np.ix_(sorted(set(range(9)) - set(selected)), selected)]
                    )
                )
                < 1e-13
            )
            small = m.source_block(setting, left, right)
            np.testing.assert_allclose(
                block[np.ix_(selected, selected)], small, atol=1e-13
            )
            # Time amplitude of the unresolved word 11 checks the sign of the
            # Fourier energy label, independently of recorded count values.
            pointer = m.energy_pointer(small).sum(axis=0)
            native_step = expm(-2j * np.pi * block / m.PERIOD)[4, 4]
            # The finite pointer wraps its last lag to its first.
            last_step = expm(2j * np.pi * (m.BINS - 1) * block / m.PERIOD)[4, 4]
            expected = ((m.BINS - 1) * native_step + last_step) / m.BINS
            np.testing.assert_allclose(
                np.fft.fft(pointer)[1], expected, atol=m.CLOCK_ERR
            )


def test_all_open_nonadjacent_source_spectra_and_local_contact_on_full_board():
    board = federation(0)
    adjacent = (1, 82)
    a, b = m.PORTS[0], m.PORTS[-1]
    assert b not in board.neighbours[a]
    assert not set(board.neighbours[a]) & set(board.neighbours[b])
    assert adjacent[1] in board.neighbours[adjacent[0]]
    diag = np.diag(physical_seam(0)).reshape(3, 3)
    edges = np.array(board.seams)
    for endpoints, expected_odd, expected_even in (
        (adjacent, -16.25, -0.75),
        ((a, b), 0, 0),
    ):
        energies = np.zeros((3, 3))
        for left, right in product(range(3), repeat=2):
            word = np.ones(board.n, dtype=int)
            word[list(endpoints)] = (left, right)
            x, y = edges.T
            energies[left, right] = diag[word[x], word[y]].sum() + 8 * np.sum(
                (6 - board.degree) * word**2
            )
            # A possible off-diagonal swap must touch a guarded source.
            unequal = edges[word[x] != word[y]]
            assert all(set(edge) & set(endpoints) for edge in unequal)
        connected = energies + energies[1, 1] - energies[:, [1]] - energies[[1], :]
        corners = connected[::2, ::2]
        np.testing.assert_allclose(corners.mean(), expected_even, atol=1e-11)
        np.testing.assert_allclose(
            (np.outer([-1, 1], [-1, 1]) * corners).mean(), expected_odd, atol=1e-11
        )


def test_every_launch_and_sign_channel_resolve_the_registered_comparison():
    result, archives = m.run()
    assert len(archives) == 36 and all(h.sum() == m.SHOTS for h in archives.values())
    with np.load(m.HERE / "counts.npz", allow_pickle=False) as saved:
        for name, counts in archives.items():
            np.testing.assert_array_equal(saved[name], counts)
    tables = result["tables"]
    for setting, row in tables.items():
        certificate = row["signed_error"]
        assert certificate is not None
        assert all(s["modal_energy_error"] is not None for s in row["sources"].values())
        if setting in ("native", "cut_conversion"):
            assert abs(row["signed_interaction"]) <= certificate
            assert all(
                s["central_counts"] == [m.SHOTS, 0, 0] for s in row["sources"].values()
            )
        elif setting == "direct":
            assert abs(row["signed_interaction"] + 16.25) <= certificate
        else:
            assert row["signed_interaction"] > certificate
            sources = row["sources"]
            assert sources["1,-1"]["polarization"] > result["polarization_error"]
            assert sources["-1,1"]["polarization"] < -result["polarization_error"]
