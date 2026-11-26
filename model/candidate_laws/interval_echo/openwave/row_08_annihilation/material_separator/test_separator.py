"""Held word banks, native matrices, the sharpened counterfactual and neutral relay."""

import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm

from . import apparatus, measure, neutral_field, separator


def test_local_exact_columns_and_physical_separators():
    rows = separator.audit()
    assert len(rows) == 6
    assert all(length == 33 and forbidden == [0, 0] for _, _, length, forbidden in rows)
    left, right = apparatus.geometry()
    assert len(set(left) & set(right)) == 0
    # Each successive field seam is physically adjacent even through the pin.
    for path in (left, right):
        assert len(set(path)) == len(path) == 7


def test_native_matrix_and_unilateral_factorization():
    h, initial, _, _ = apparatus.setup()
    # In a fixed left-middle sector the last left dimer and receiver are a
    # separate system. Test the full matrix split using bit permutations.
    idx = np.arange(h.shape[0])
    far = (1 << 5) | 2
    for pin in (0, 1):
        distant = idx[((idx & ~far) == 0)]
        source = idx[((idx & far) == 0) & (((idx & apparatus.MIDDLE[0]) != 0) == pin)]
        order = (source[:, None] | distant).ravel()
        physical = h[order, :][:, order]
        left = physical[::4, ::4]
        right = physical[:4, :4].toarray() - np.eye(4) * physical[0, 0]
        split = sparse.kron(left, sparse.eye(4)) + sparse.kron(
            sparse.eye(left.shape[0]), right
        )
        assert abs(physical - split).max() < 1e-12
    # A middle flip has been deleted; the adjacent diagonal contact energies
    # in both the tested physical blocks above remain in the pinned operator.
    pinned = apparatus.frozen(h, (0,))
    assert abs(pinned[idx, idx ^ apparatus.MIDDLE[0]]).max() == 0
    removed = (pinned - h).tocoo()
    np.testing.assert_array_equal(removed.row ^ removed.col, apparatus.MIDDLE[0])
    np.testing.assert_allclose(removed.data, 0.5, atol=2e-15, rtol=0)
    assert removed.nnz == h.shape[0]
    assert np.allclose(initial.T @ initial, np.eye(2))


def test_field_operator_and_direct_read():
    xhop = 3.1
    paths = neutral_field.field(xhop)
    k = np.arange(1, 7)
    single = [1 << j for j in range(6)]
    np.testing.assert_allclose(
        np.linalg.eigvalsh(paths[np.ix_(single, single)]),
        np.sort(2 * xhop * np.cos(k * np.pi / 7)),
        atol=1e-13,
    )
    cut = neutral_field.field(xhop, True)
    left_bits = np.array([int(i >> 3).bit_count() for i in range(64)])
    assert np.max(abs(cut * (left_bits[:, None] - left_bits))) == 0
    assert np.max(abs(paths * (left_bits[:, None] - left_bits))) > 0
    # Independently exponentiate the physical far dimer and meter with the
    # field contact absent. The pin does not inhibit local far/meter exchange.
    gap = np.hypot(32.5, 1)
    native = np.array([[0, -0.5], [-0.5, 32.5]])
    q = np.array([[0, 1.0], [1.0, 0]])
    h = (
        np.kron(native, np.eye(2))
        + np.kron(np.eye(2), np.diag([0, gap]))
        + np.kron(q, q) / 16
    )
    excited = eigh(native)[1][:, -1]
    state = expm(-25j * h) @ np.kron(excited, [1, 0])
    assert np.linalg.norm(state[1::2]) ** 2 > 0.95


def test_held_words_and_detector_gates():
    with np.load(measure.HERE / "counts.npz") as raw:
        result = measure.decode(raw)
        for key in raw.files:
            if not key.endswith("_counts"):
                continue
            bank = raw[key]
            law = raw[key.removesuffix("_counts") + "_probabilities"]
            np.testing.assert_array_equal(
                bank.sum(axis=-1), np.full(bank.shape[:-1], measure.SHOTS)
            )
            np.testing.assert_allclose(law.sum(axis=-1), 1, atol=5e-11, rtol=0)
            se = np.sqrt(np.maximum(0, law * (1 - law)) / measure.SHOTS)
            assert np.all(abs(bank / measure.SHOTS - law) < 8 * se + 3 / measure.SHOTS)
        saved = json.loads((measure.HERE / "results.json").read_text())
        for arm, sides in apparatus.ARMS.items():
            sample = result["native"][arm]["observed"]
            exact = result["native"][arm]["born"]
            for side, name in enumerate(("left", "right")):
                metric = name + "_excess"
                assert (
                    abs(
                        sample[-1][metric]
                        - saved["native"][arm]["observed"][-1][metric]
                    )
                    < 1e-12
                )
                if side in sides:
                    assert max(abs(row[metric]) for row in sample) < 0.001
                    assert max(abs(row[metric]) for row in exact) < 1e-9
                elif not sides:
                    assert sample[-1][metric] > 0.002
            if len(sides) == 2:
                assert max(row["largest_bin_difference"] for row in sample) < 0.001
                assert max(row["largest_bin_difference"] for row in exact) < 1e-9
            for pins in ("00", "01", "10", "11"):
                rows = result["field"][pins][arm]["observed"]
                if not sides:
                    assert any(
                        min(r["left_excess"], r["right_excess"]) > 0.01
                        and r["both_excess"] > 0.0001
                        for r in rows
                        if r["time"] in (64, 128)
                    )
                for side in sides:
                    metric = ("left", "right")[side] + "_excess"
                    assert max(abs(row[metric]) for row in rows) < 0.001


def test_known_open_receivers_and_independent_field_propagation():
    from functools import reduce
    from scipy.sparse.linalg import expm_multiply

    with np.load(measure.HERE / "counts.npz") as raw:
        # Independent tensor-product construction of a calibration arm. The
        # far dimer at middle orientation zero has diagonal costs 32.5 and 0
        # in some order. Global material and meter X sign conventions give
        # the same energy probabilities for the reversed pin.
        omega = float(raw["receiver_gap"])
        x = np.array([[0, 1], [1, 0]])
        y = np.array([[0, -1j], [1j, 0]])
        i = np.eye(2)

        def local(ops):
            return reduce(np.kron, [ops.get(k, i) for k in range(8)])

        def xx(j):
            return local({j: x, j + 1: x})

        def yy(j):
            return local({j: y, j + 1: y})

        fields = [float(raw["field_hopping"]) * (xx(j) + yy(j)) / 2 for j in range(5)]
        chain = sum(fields)
        # fields are tensored with the far and meter identities.
        values, vec = eigh(chain[::4, ::4])
        ground = vec[:, 0]
        assert values[0] < values[1]
        far = np.array([[0, -0.5], [-0.5, 32.5]])
        far_ground = eigh(far)[1][:, 0]
        start = np.kron(ground, np.kron(far_ground, [1, 0]))
        states = np.column_stack([start, local({0: x}) @ start])
        full = (
            chain
            + local({6: far})
            + local({7: np.diag([0, omega])})
            + (xx(5) + xx(6)) / 16
        )
        for name, cut in [("open", False), ("both", True)]:
            h = sparse.csr_matrix(full - (fields[2] if cut else 0))
            times = 64
            output = expm_multiply(
                -1j * times * h, states, traceA=-1j * times * h.diagonal().sum()
            )
            clicks = (abs(output[1::2]) ** 2).sum(axis=0)
            recorded = raw["field_00_" + name + "_probabilities"][5]
            np.testing.assert_allclose(
                clicks, recorded[:, 2:].sum(axis=1), atol=2e-10, rtol=0
            )
        original = json.loads((apparatus.continuous.HERE / "results.json").read_text())
        panel = original["arms"]["connected"]["observations"]["64"]
        prob = raw["open_probabilities"][-1].reshape(2, 16, 4, 4).sum(axis=2)
        for src, key in enumerate(apparatus.continuous.BASES):
            np.testing.assert_allclose(
                prob[src], panel[key]["ring_word_by_receiver"], atol=2e-11, rtol=0
            )
