"""Checks independent of the Chebyshev and four-phase record producer."""

import json
import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.row_31_orbitals.test_orbitals import (
    assembled_cell,
)
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index
from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair import (
    measure as m,
)


def test_seam_gap_and_full_cell():
    law = candidate.native_law
    k, v = law.local_term_parts()
    ports = [law.SRC.index(x) for x in ((1, 2), (2, 1))]
    local = (k + candidate.KAPPA * v)[np.ix_(ports, ports)]
    np.testing.assert_allclose(np.linalg.eigvalsh(local), [8, 9], atol=1e-12)
    np.testing.assert_allclose(expm(-1j * np.pi * local), [[0, 1], [1, 0]], atol=1e-12)
    board = cell()
    pairs, direct = assembled_cell(0, (1, 2))
    _, _, adj, kinetic, _, _ = generators(board)
    np.testing.assert_array_equal(
        [ordered_index(*p, board.n) for p in pairs], np.arange(len(pairs))
    )
    expected = kinetic.toarray() + np.diag(240 - 32 * adj)
    np.testing.assert_allclose(direct, expected, atol=1e-12)


def test_crossing_with_independent_full_exponential():
    board = federation(0)
    x, y, adjacent, kinetic, _, _ = generators(board)
    h = kinetic - sparse.diags(32 * adjacent, dtype=float)
    a, b = m.fixtures(board)[0]
    source = np.zeros(len(x))
    source[ordered_index(*a, board.n)] = 1
    reverse = ordered_index(y, x, board.n)
    assert (h[reverse][:, reverse] != h).nnz == 0
    t = m.TIMES[16]
    state = expm_multiply(-1j * t * h, source)
    f, r = state[[ordered_index(*b, board.n), ordered_index(*b[::-1], board.n)]]
    assert np.linalg.norm(state[adjacent]) ** 2 >= 27 / 32
    reference = (f + r).conj() * (f - r)
    archive = json.loads((m.HERE / "records.json").read_text())
    scan = archive["settings"]["registered"]["scans"][0]
    _, z, dz = m.decode(scan["orders"])
    assert 0 < abs(z[16] - reference) < dz[16] < abs(z[16])
    assert h.diagonal().max() <= 12
    assert h[ordered_index(*a, board.n), ordered_index(*a, board.n)] == -27


def test_decode_independent_followup_and_erasure():
    archive = json.loads((m.HERE / "records.json").read_text())
    for arm in archive["settings"].values():
        for route in arm["scans"]:
            assert route["measured"]["best_resolved_hand"]["excludes_zero"]
            assert not route["erased"]["best_resolved_hand"]["excludes_zero"]
    for scan in archive["independent_followup"]["scans"]:

        def four_phase(counts):
            outcomes = np.array(counts)
            assert np.all(outcomes.sum(axis=-1) == 2**38)
            differences = (outcomes[:, :, 0] - outcomes[:, :, 1]) / 2**38
            return (differences @ np.array([1, 1j, -1, -1j])) / 2

        _, transmitted, radius = m.decode(scan["orders"], 2**38)
        # The two quarter-turn source phases have identical sharp-load populations.
        assert np.all(abs(transmitted.imag) > 4 * radius)
        corrected = four_phase(scan["orders"]) * four_phase(scan["bench"]).conj()
        witness = scan["witness"]
        value = abs((corrected[0] * corrected[1].conj()).imag)
        np.testing.assert_allclose(
            value, witness["imaginary_product"], atol=2e-25, rtol=0
        )
        assert value > 5 * witness["simultaneous_allowance"]
    for row in archive["initial_energy_readings"]:
        for read in row:
            assert read["even"] + read["odd"] == 2**20
            assert abs(read["mean_energy"] - 213) < 0.02
