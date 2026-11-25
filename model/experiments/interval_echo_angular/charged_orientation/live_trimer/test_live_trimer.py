"""Independently evolve a native connected component; check the phase and null controls."""

import json
from itertools import product

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from . import detectors, measure


def records():
    return json.loads((measure.HERE / "results.json").read_text())


def test_native_active_sectors_and_pinned_parity():
    blocks = measure.blocks()
    assert sum(a.shape[0] for a in blocks.values()) == 99
    for label, expected in [
        ("blank", {-4}),
        ((0, 1), set()),
        ((1, 1), {-32}),
        ((1, -1), {-36}),
        ((2, 1), {-128}),
        ((2, -1), {-132}),
    ]:
        h = blocks[label]
        off = h[np.triu_indices(len(h), 1)]
        assert set(off[off != 0]) == expected

    # In the original independent sign basis cross terms cancel when the
    # launch head is traced, even when the two parity branches evolve differently.
    a, ix = measure.sector(1, True, True)
    b, _ = measure.sector(-1, True, True)
    even = measure.branch(measure.spectrum(a, ix), ix, 131072)[:, :, 0]
    odd = measure.branch(measure.spectrum(b, ix), ix, 131072)[:, :, 0]
    rho = (even.T @ even.conj() + odd.T @ odd.conj()) / 2
    for phase in [1, -1, 1j, -1j]:
        plus, minus = (even + phase * odd) / 2, (even - phase * odd) / 2
        np.testing.assert_allclose(
            plus.T @ plus.conj() + minus.T @ minus.conj(), rho, atol=1e-13
        )
    assert np.trace(rho).real > 0.4


def test_long_flight_direct_exponential_and_receiver_transpose():
    h, ix = measure.sector(1, True, True)
    sources = [ix[w, (0, 0, 0)] for w in [(0, 2, 0, 1, 1, 1), (0, 0, 2, 1, 1, 1)]]
    _, labels = connected_components(h - np.diag(h.diagonal()), directed=False)
    # The tested orbitals are in the same dynamical component.
    assert labels[sources[0]] == labels[sources[1]]
    component = np.flatnonzero(labels == labels[sources[0]])
    work = expm((-1j / 8) * h[np.ix_(component, component)])
    for _ in range(17):
        work = work @ work
    full = np.zeros((len(ix), 2), complex)
    full[component] = work[:, np.searchsorted(component, sources)]
    maps = []
    for coins in product(range(2), repeat=3):
        maps.append(
            full[[ix[w, coins] for w in [(1, 1, 1, 0, 2, 0), (1, 1, 1, 0, 0, 2)]]]
        )
    # Direct partial transpose from |in> <in'| times the transposed
    # receiver operator on every independent final memory word.
    basis = np.eye(2)
    partial = np.zeros((4, 4), complex)
    for mat in maps:
        for i, j in product(range(2), repeat=2):
            partial += (
                np.kron(
                    np.outer(basis[i], basis[j]),
                    np.outer(mat[:, i], mat[:, j].conj()).T,
                )
                / 2
            )
    d = records()
    row = d["settings"]["matched"]["reads"]["131072"]
    np.testing.assert_allclose(
        np.linalg.eigvalsh(partial),
        row["choi_partial_transpose_eigenvalues"],
        atol=1e-8,
    )
    assert np.isclose(np.trace(partial).real, row["mixed_code_acceptance"], atol=1e-8)
    estimated = detectors.measured_partial(
        row["pauli_counts"], d["starts_per_pauli_read"]
    )[0]
    assert estimated + d["simultaneous_pt_radius"] < -0.01
    assert (
        row["body_full_turn_mixed_counts"][1] == 0
        and row["body_full_turn_mixed_counts"][0] > 28000000
    )
    for settings in d["settings"].values():
        for obs in settings["reads"].values():
            assert np.all(
                np.sum(obs["pauli_counts"], axis=-1) == d["starts_per_pauli_read"]
            )
            assert (
                abs(
                    obs["measured_pt_eigenvalues"][0]
                    - obs["choi_partial_transpose_eigenvalues"][0]
                )
                <= d["simultaneous_pt_radius"]
            )


def test_destination_mirror_and_closed_gate():
    h, ix = measure.sector(1, True)
    reorder = []
    for loads, coins in ix:
        mirror = list(loads)
        mirror[4], mirror[5] = mirror[5], mirror[4]
        reorder.append(ix[tuple(mirror), coins])
    np.testing.assert_array_equal(h, h[np.ix_(reorder, reorder)])
    for data in records()["settings"]["even"]["reads"].values():
        assert (
            np.linalg.norm(
                np.asarray(data["sharp_load_outputs"])[:, 1]
                - np.asarray(data["sharp_load_outputs"])[:, 2]
            )
            < 1e-12
        )
        assert data["choi_partial_transpose_eigenvalues"][0] > -1e-16
    closed, ix = measure.sector(1, False, True)
    _, labels = connected_components(
        closed - np.diag(closed.diagonal()), directed=False
    )
    start = labels[ix[(0, 2, 0, 1, 1, 1), (0, 0, 0)]]
    assert all(labels[idx] != start for (w, c), idx in ix.items() if w[:3] == (1, 1, 1))
    for v in records()["settings"]["closed"]["reads"].values():
        assert v["mixed_code_acceptance"] == 0
    # One proper face turn moves the three physical positions; three restore all.
    turn = detectors.body_turn()
    assert np.trace(turn) == 0
    np.testing.assert_array_equal(np.linalg.matrix_power(turn, 3), np.eye(3))
