import json
from itertools import product

import numpy as np
from scipy.linalg import eigh, expm

from . import clock_check, measure as m


def labelled_wing(k, ports):
    """Compress five literal two-symbol swap matrices on the 64-word tensor product."""
    words = np.array(list(product((0, 2), repeat=6)))
    code = np.where(words.sum(axis=1) == 6)[0]
    h = np.zeros((64, 64))
    for a in range(5):
        swap = np.eye(4)[[0, 2, 1, 3]]
        contact = (np.eye(4) - swap) / 2 + k * np.diag([0.0, 4.0, 4.0, 0.0])
        h += np.kron(np.kron(np.eye(2**a), contact), np.eye(2 ** (4 - a)))
    h += np.diag(k * (words**2) @ (6 - m.federation(0).degree[ports]))
    return words[code], h[np.ix_(code, code)]


def test_released_exact_floor_and_native_compression():
    for ports in m.wall_flight.extended_geometry()[2]:
        words, h = labelled_wing(8, ports)
        _, exchange = labelled_wing(0, ports)
        actual, labels = m.full_wing(8, ports)
        indices = [np.nonzero((labels == w).all(1))[0].item() for w in words]
        np.testing.assert_allclose(actual[np.ix_(indices, indices)], h, atol=1e-14)
        cost = h - exchange
        assert np.array_equal(cost, np.diag(np.diag(cost)))
        values, counts = np.unique(np.diag(cost), return_counts=True)
        np.testing.assert_array_equal(values, [64, 96, 128, 160, 192, 224])
        np.testing.assert_array_equal(counts, [1, 1, 5, 5, 6, 2])
        assert np.linalg.eigvalsh(exchange).min() >= -1e-13
        arg = np.diag(cost).argmin()
        assert tuple(words[arg]) == (2, 2, 2, 0, 0, 0) and h[arg, arg] == 64.5
        observed = np.diff(np.linalg.eigvalsh(h)[:2])[0]
        assert observed >= 96 - 64.5 and abs(observed - 32.49357657675) < 1e-10


def test_soft_edge_doublet_and_uniform_gap():
    # Independent free-chain Majorana singular-value route and variational trial.
    for n in range(1, 9):
        coupled = np.diag(np.full(n, 0.5)) + np.diag(np.full(n - 1, 16.25), 1)
        d = 2 * np.linalg.svd(coupled, compute_uv=False)[::-1]
        assert d[0] <= 1 * (0.5 / 16.25) ** (n - 1) * (1 + 2e-12)
        if n >= 2:
            assert d[1] >= 31.5
        if n <= 3:
            e = np.linalg.eigvalsh(m.walls.wing_ref(8, 0, True, n)[0])
            assert abs(d[0] - (e[1] - e[0])) < 1e-12


def test_all_records_count_charge_and_counted_clock_transfer():
    r = json.loads((m.HERE / "records.json").read_text())
    clocks = json.loads((m.HERE / "clocks.json").read_text())
    for arm in r["arms"]:
        for switch, condition in arm["conditions"].items():
            for name, v in condition["sources"].items():
                assert sum(v["joint_energy_counts"]) == m.N
                for side, readings in enumerate(v["local"]):
                    labels = np.array(condition["calibration"][0][side]["word_labels"])
                    for record in readings.values():
                        n = np.array(record["counts"])
                        assert n.sum() == m.N
                        byword = n[:-1].reshape(2, -1)[:, :-1].sum(axis=0)
                        q = labels[:, -2:].sum(axis=1) - 2
                        assert [
                            int(byword[q == a].sum()) for a in (-2, 0, 2)
                        ] == record["outer_charge_counts"]
                        if switch != "release":
                            assert not (byword[q != 0]).any()
            if arm["k"] == 8 and not arm["cut_in_flight"]:
                c = clocks["settings"][switch]
                for phase, read in enumerate(c["clocks"]):
                    lo, hi = clock_check.read_clock(read["quadrature_counts"])
                    assert (
                        lo
                        <= condition["calibration"][phase][0]["excitation_costs"][1]
                        <= hi
                    )
                for i, name in enumerate(["sham", "wall"]):
                    other = c["held"][i]
                    assert sum(other["complete_counts"]) == m.N
                    a, b = condition["sources"][name]["local"][0]["ideal"][
                        "line_given_ordered_interval"
                    ]
                    lo, hi = other["conditional_rate"]
                    assert max(a, lo) <= min(b, hi)
    frequency = clocks["settings"]["detach"]["clocks"][0]["gap_interval"]
    assert frequency[1] - frequency[0] < 4.1e-10
    control = r["arms"][2]["conditions"]["release"]["sources"]["wall"]["local"][0][
        "ideal"
    ]
    assert control["split_given_line_interval"][0] > 0.79


def test_independent_source_and_dyadic_receiver():
    ring, anchors, wings = m.wall_flight.extended_geometry()
    states = np.array(list(product((0, 1), repeat=10)))
    loads = np.stack([2 * states, 2 - 2 * states], axis=-1).reshape(-1, 20)
    contacts = [
        (1, 2),
        (3, 4),
        (5, 6),
        (7, 0),
        (anchors[0], 8),
        (anchors[1], 14),
        (9, 10),
        (11, 12),
        (15, 16),
        (17, 18),
    ]
    mismatch = sum((loads[:, a] != loads[:, b]).astype(int) for a, b in contacts)
    h = np.diag(32.5 * mismatch)
    for j in range(10):
        h[np.arange(1024), np.arange(1024) ^ (1 << j)] -= 0.5
    e, v = eigh(h)
    _, starts, order = m.reference()
    wave = ((v * np.exp(-32j * e)) @ (v.T @ starts)).reshape(16, 8, 8, 2)
    file = json.loads((m.HERE / "clocks.json").read_text())
    base = [tuple(m.walls.loads(b)) for b in product((0, 1), repeat=3)]
    words, local = labelled_wing(8, wings[0])
    ix = [next(i for i, w in enumerate(words) if tuple(w) == target) for target in base]
    for name in m.SETTINGS:
        bank = file["settings"][name]
        for src in (0, 1):
            total = np.zeros(len(words) if name == "release" else 8)
            found = wave[order, :, :, src]
            for phase, state in enumerate(found):
                if name == "release":
                    target = local
                    psi = np.zeros((20, 8), complex)
                    psi[ix] = state
                else:
                    target = local[np.ix_(ix, ix)]
                    if name == "keep":
                        boundary = loads[order[phase] * 64, anchors[0]]
                        target = target + np.diag(32.5 * (words[ix, 0] != boundary))
                    psi = state
                levels = np.linalg.eigvalsh(target)
                aim = np.mean(bank["clocks"][phase]["gap_interval"])
                step = expm(
                    -1j * m.DT * (target - (levels[0] + aim) * np.eye(len(target)))
                )
                # Coherent binary holds are implemented by direct 8/20-state unitaries.
                for _ in range(28):
                    psi = (psi + step @ psi) / 2
                    step = step @ step
                total += (abs(psi) ** 2).sum(axis=1)
            n = np.asarray(bank["held"][src]["complete_counts"])
            weights = n[:-1].reshape(2, -1)[:, :-1].sum(axis=0) / m.N
            # Per-word and absolute rate tolerance for these independent simulated draws.
            np.testing.assert_allclose(weights, total, atol=9e-5, rtol=0)
            assert abs(n[-1] / m.N - (1 - np.linalg.norm(found) ** 2)) < 9e-5
