import hashlib
import json
from itertools import combinations, product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import Carrier, federation
from . import free_exchange, measure


def test_blank_sector_is_an_invariant_sector_of_the_native_contact():
    exchange, contrast = law.local_term_parts()
    h = exchange + 8 * contrast
    np.testing.assert_allclose(h[:9, :9], law.seam_term(8), atol=2e-14)
    np.testing.assert_allclose(h[9:, :9], 0, atol=2e-14)
    swap = np.eye(9)[[law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(expm(-1j * np.pi * h[:9, :9]), swap, atol=4e-14)


def test_both_blank_generators_against_port_load_basis():
    graph = Carrier("four-port branch", 4, [(0, 1), (1, 2), (2, 3), (1, 3)])
    for kappa in (0, 8):
        states = list(product(range(3), repeat=4))
        full = np.zeros((len(states), len(states)))
        for col, n in enumerate(states):
            full[col, col] = kappa * np.dot(6 - graph.degree, np.square(n))
            for i, j in graph.seams:
                seam = law.seam_term(kappa)[:, law.SRC.index((n[i], n[j]))]
                for out in np.flatnonzero(seam):
                    next_n = list(n)
                    next_n[i], next_n[j] = law.SRC[out]
                    full[states.index(tuple(next_n)), col] += seam[out]
        for sign in (-1, 1):
            for count in (1, 2):
                cols = []
                for occupied in combinations(range(4), count):
                    n = np.ones(4, dtype=int)
                    n[list(occupied)] += sign
                    cols.append(states.index(tuple(n)))
                other = list(set(range(len(states))) - set(cols))
                np.testing.assert_array_equal(full[np.ix_(other, cols)], 0)
                h = full[np.ix_(cols, cols)] - kappa * (
                    sum(6 - graph.degree) + 6 * count
                ) * np.eye(len(cols))
                if count == 1:
                    reduced = measure.generator(graph, sign, kappa)
                else:
                    reduced = free_exchange.generator(
                        graph, sign, False, free_exchange.pair_basis(4), kappa
                    )
                np.testing.assert_allclose(reduced.toarray(), h, atol=1e-13)


def test_single_defect_counts_and_rotations():
    result = json.loads((measure.HERE / "measure.json").read_text())
    assert (
        result["protocol_sha256"]
        == hashlib.sha256((measure.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    passed = 0
    for board in result["native"]:
        for arm in board["arms"]:
            for appointment in arm["appointments"]:
                for receiver in appointment["receivers"]:
                    groups = np.asarray(receiver["counts"])
                    assert groups.shape == (6, 3, 4)
                    np.testing.assert_array_equal(groups.sum(axis=-1), measure.SHOTS)
                    assert (
                        abs(
                            receiver["measured_negative_pt"]
                            - receiver["exact_negative_pt"]
                        )
                        < 5 * measure.EPS
                    )
                    passed += bool(
                        arm["kappa"]
                        and appointment["exposure"]
                        and receiver["resolved"]
                    )
                    for name, turns in receiver["turns"].items():
                        for counts in turns.values():
                            np.testing.assert_array_equal(
                                np.sum(counts, axis=1), measure.SHOTS
                            )
                        if name in ("idle", "360", "720", "xy_vs_yx"):
                            assert all(row[1] == 0 for row in turns["orbital"])
                            dark = name in ("360", "xy_vs_yx")
                            assert all(
                                row[0 if dark else 1] == 0
                                for row in turns["imposed_spin_half"]
                            )
    assert passed == 28


def test_compiled_exchange_reaches_the_other_reference():
    result = json.loads((measure.HERE / "measure.json").read_text())
    code = measure.code_columns([0, 1, 2], 3)
    doublet = np.kron(code, code)
    for level, record in enumerate(result["native"]):
        board = federation(level)
        src, dests, eligible = measure.design(board)
        encoded_src, encoded_dst = eligible[src][1], eligible[dests[0]][1]
        compiled = measure.exchange_word(board, encoded_src, encoded_dst)
        assert compiled == [tuple(e) for e in record["exchange_word"]]
        pos, neg = measure.apply_pair_word(encoded_src, encoded_dst, compiled)
        swap = np.eye(4)[[0, 2, 1, 3]]
        np.testing.assert_allclose(doublet.conj().T @ pos @ doublet, swap, atol=2e-14)
        np.testing.assert_allclose(doublet.conj().T @ neg @ doublet, -swap, atol=2e-14)
        for arm in record["exchange"].values():
            refs = arm["registered"]["references"]
            np.testing.assert_allclose(
                [refs["original_pairing"], refs["crossed_pairing"]], [0.25, 1]
            )
            assert arm["registered"]["parity"]["symmetric"][1] == 0
            assert arm["odd_worldlines"]["parity"]["symmetric"][0] == 0
    for field in result["written_extension"]:
        refs = field["readout"]["written_head"]["references"]
        np.testing.assert_allclose(
            [refs["original_pairing"], refs["crossed_pairing"]], [1, 0.25]
        )


def test_free_scattering_has_complements_and_a_distinguishable_exchange_lift():
    panel = json.loads((measure.HERE / "free_exchange.json").read_text())
    assert (
        panel["protocol_sha256"]
        == hashlib.sha256((measure.HERE / "FREE_EXCHANGE.md").read_bytes()).hexdigest()
    )
    assert panel["basis_size"] == 28680
    for arm in panel["arms"]:
        for row in arm["appointments"]:
            assert row["separated"] == bool(row["exposure"])
            for law_name in ("registered", "odd_worldlines"):
                for read_name in ("character_counts", "reference_counts"):
                    np.testing.assert_array_equal(
                        np.sum(row[law_name][read_name], axis=1), measure.SHOTS
                    )
            if row["exposure"] == 0:
                np.testing.assert_allclose(row["registered"]["accepted"], 1)
                assert row["odd_return_weight"] == 0
    last = panel["arms"][1]["appointments"][-1]
    np.testing.assert_allclose(
        [last["registered"]["accepted"], last["odd_worldlines"]["accepted"]],
        [0.1887418465, 0.1837129453],
        atol=2e-10,
    )


def test_tomography_separates_a_full_qubit_from_a_single_transmitted_mode():
    analyzer = np.eye(2)
    rng = np.random.default_rng(327)
    coherent = measure.tomography(np.eye(2) / 2, analyzer, rng)
    single_mode = measure.tomography(np.diag([1 / np.sqrt(2), 0]), analyzer, rng)
    np.testing.assert_allclose(coherent["present"], single_mode["present"])
    np.testing.assert_allclose(coherent["exact_negative_pt"], 0.125)
    assert coherent["resolved"]
    assert single_mode["exact_negative_pt"] == 0
    assert not single_mode["resolved"]
