"""Compare the pair-sector reduction directly to native 99-dimensional seam gates."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations

import numpy as np
from scipy.linalg import eigh

from openwave.interval_echo.row_19.selected_collar import parked_pair_flight as flight


def test_two_load_reduction_and_exponential_against_native_seams():
    routes = flight.source.reference.local.survey.routes
    volume = routes.volume
    law = volume.candidate.native_law
    board = routes.cell()
    home = (0, 1)
    pairs = list(combinations(range(board.n), 2))
    idx = {pair: i for i, pair in enumerate(pairs)}
    endpoint = ((0, 0), (1, 0), (0, 1), (1, 1))
    one = np.eye(law.N_SRC)[:, [law.SRC.index(x) for x in endpoint]]
    exchange, contrast = law.local_term_parts()
    h = exchange + volume.candidate.KAPPA * contrast
    words = routes.codebook()[0]
    marks = flight.marks_for(board, (0,), (5,), words)
    memories = {
        None: np.eye(11)[0],
        0: volume.label(0, 1),
        32: volume.label(2, 1),
        33: volume.label(2, -1),
    }
    small = {}
    for tag, mem in memories.items():
        sharp = np.kron(mem[:, None], one)
        small[tag] = sharp.conj().T @ h @ sharp
        np.testing.assert_allclose(h @ sharp, sharp @ small[tag], atol=1e-12, rtol=0)
    H = np.zeros((len(pairs), len(pairs)))
    for k, pair in enumerate(pairs):
        H[k, k] = sum(8 * (6 - board.degree[a]) for a in pair)
        for e in board.seams:
            a, b = e
            block = small[marks.get(e)]
            j = endpoint.index((int(a in pair), int(b in pair)))
            for outcome, (left, right) in enumerate(endpoint):
                if abs(block[outcome, j]) < 1e-13:
                    continue
                new = set(pair) - {a, b}
                if left:
                    new.add(a)
                if right:
                    new.add(b)
                H[idx[tuple(sorted(new))], k] += block[outcome, j]
    sector = flight.FreePair(board, home)
    hop, diag, _ = sector.operator(marks)
    reduced = (hop.toarray() + np.diag(diag)) / 2 + np.eye(len(pairs)) * sum(
        marks.values()
    )
    np.testing.assert_allclose(H, reduced, rtol=0, atol=5e-13)
    evals, basis = eigh(H)
    expected = abs(basis[idx[home]] ** 2 @ np.exp(-1j * evals)) ** 2
    weights, error = flight.flight_weights()
    assert abs(sector.chebyshev((hop, diag), 138, weights) - expected) < error
    series = flight.return_series(board, marks, home, 2)
    variance = (H @ H)[idx[home], idx[home]] - H[idx[home], idx[home]] ** 2
    assert abs(float(series[1]) + variance) < 1e-10


def test_free_response_and_ordinary_controls():
    z = flight.source.reference.local.inputs("parked_pair_flight")
    for name, digest in z["inputs"].items():
        assert (
            sha256((flight.source.reference.local.HERE / name).read_bytes()).hexdigest()
            == digest
        )
    native = z["change_in_joint_return"]["native"]
    changed = z["change_in_joint_return"]["remove_marked_pair_binding"]
    assert z["release_and_return_ports"] == [10, 20]
    assert F(native["shared"][4]) == F(-95, 13934592) and F(native["plus"][4]) > 0
    at_one = z["unit_time_joint_return"]
    observed = at_one["native"]["changes"]["shared"]
    control = at_one["remove_marked_pair_binding"]["changes"]["shared"]
    assert -5.008e-8 < observed[0] < observed[1] < -5.007e-8
    assert 1.883e-10 < observed[0] - control[1] <= observed[1] - control[0] < 1.910e-10
    assert z["ordinary_unit_time_joint_return"][0] > 0.976
    for tag in ("plus", "minus", "fourth_order", "shared"):
        assert native[tag][:4] == ["0"] * 4
        lead = 5 if tag == "shared" else 6
        assert (
            native[tag][:lead] == changed[tag][:lead]
            and native[tag][lead] != changed[tag][lead]
        )
        for mode in z["unit_time_joint_return"].values():
            assert len(mode["by_source_word_pair"]) == 81
            assert all(0 <= x <= 1 for x in mode["by_source_word_pair"])


def test_equal_positive_charge_preserves_class_records_and_changes_transport():
    local = flight.source.reference.local
    z = local.inputs("matched_charge_transport")
    measured = local.inputs("parked_energy_charge")
    for name, digest in z["inputs"].items():
        assert sha256((local.HERE / name).read_bytes()).hexdigest() == digest
    u, v = np.asarray(z["u"]), np.asarray(z["v"])
    diff = np.outer(u, v) + np.outer(v, u)
    words = measured["classes_in_calibration_order"]
    flags = np.eye(3, dtype=int)[words]
    assert not np.any(flags.T @ diff @ flags)
    assert sum(sum(r) for r in z["joint_class_counts_on_both"]) == z["denominator"]
    for pair, bins in z["nonzero_relay_histogram_count_differences"].items():
        assert not sum(bins.values())
        mean = F(sum(int(t) * count for t, count in bins.items()), z["denominator"])
        assert mean == F(z["nonzero_mean_relay_changes"].get(pair, "0"))
    assert z["nonzero_mean_relay_changes"]["1-16"] == "-2453/221184"
    assert len(z["nonzero_relay_histogram_count_differences"]) == 7
    assert z["common_charge_relative_to_fresh"][0] > 0
    assert z["unit_time_free_pair_return_change"]["native"][1] < -4.36e-7
