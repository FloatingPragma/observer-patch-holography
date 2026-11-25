"""Phase transport of three net-negative defects and the complete local parity basis."""

import argparse
import json
from functools import cache
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import detectors

LAW = candidate.native_law
HERE = Path(__file__).resolve().parent
TIMES = (0, 0.25, 0.5, 1, 2, 4, 8, 16, 32)


@cache
def blocks(kappa=candidate.KAPPA):
    """Integer 8H on all 99 states. Compress in exact orthogonal sector coordinates."""
    swap, cost = LAW._load_ops()
    reverse, d2 = LAW._label_ops()
    w = LAW._write_pattern()
    q = 2 * np.eye(90, dtype=int) - w @ w.T
    k = int(kappa)
    assert k == kappa
    h8 = np.zeros((99, 99), dtype=int)
    h8[:9, :9] = 4 * (np.eye(9, dtype=int) - swap) + 8 * k * cost
    h8[9:, 9:] = (
        q
        @ np.kron(np.eye(10, dtype=int) - reverse + 2 * k * d2, np.eye(9, dtype=int))
        @ q
    )
    ex, con = LAW.local_term_parts()
    np.testing.assert_allclose(h8, 8 * (ex + k * con), atol=1e-12)
    orientation = list(range(99))
    for word, (i, j) in enumerate(LAW.SRC):
        other = LAW.SRC.index((j, i))
        orientation[word] = other
        for label, (coin, d) in enumerate(LAW.LABELS):
            orientation[9 + 9 * label + word] = (
                9 + 9 * LAW.LABELS.index((1 - coin, -d)) + other
            )
    assert np.array_equal(h8[np.ix_(orientation, orientation)], h8)
    assert np.array_equal(np.array(orientation)[orientation], np.arange(99))
    bases = {"blank": np.eye(99, 9, dtype=int)}
    for m in range(3):
        for sign in [1] if m == 0 else [-1, 1]:
            b = np.zeros((99, 18), dtype=int)
            for coin, word in product(range(2), range(9)):
                j = 9 * coin + word
                b[9 + 9 * LAW.LABELS.index((coin, m)) + word, j] = 1
                if m:
                    b[9 + 9 * LAW.LABELS.index((coin, -m)) + word, j] = sign
            bases[m, sign] = b
    basis = np.column_stack(list(bases.values()))
    gram = basis.T @ basis
    assert basis.shape == (99, 99) and np.all(gram == np.diag(gram.diagonal()))
    reduced = {}
    for label, b in bases.items():
        norm = int(b[:, 0] @ b[:, 0])
        matrix = b.T @ h8 @ b
        assert not np.any(matrix % norm)
        matrix //= norm
        assert np.array_equal(h8 @ b, b @ matrix)
        off = matrix - np.diag(matrix.diagonal())
        assert np.max(off) == 0
        reduced[label] = matrix
    return reduced


@cache
def geometry():
    board = federation(0)

    def triangle(g):
        for i, j in combinations(board.neighbours[g], 2):
            if board.degree[i] == board.degree[j] == 5 and j in board.neighbours[i]:
                return [g, i, j]
        raise ValueError("The gateway has no triangle with two degree-five tips")

    a, b = board.external[0]
    return board, triangle(a) + triangle(b)


def sector(source_parity, bridge_live, matched=False):
    board, ports = geometry()
    raw = blocks()
    # Coin and load conventions follow the displayed edge orientations.
    edges = [
        (0, 1, (1, source_parity), 0),
        (0, 2, "blank", None),
        (1, 2, "blank", None),
        (3, 4, (1, 1) if matched else "blank", 2 if matched else None),
        (3, 5, "blank", None),
        (4, 5, "blank", None),
        (0, 3, (1 if bridge_live else 0, 1), 1),
    ]
    assert all(ports[j] in board.neighbours[ports[i]] for i, j, _, _ in edges)
    words = []
    for holes in combinations(range(6), 2):
        for excess in set(range(6)) - set(holes):
            word = tuple(0 if i in holes else 2 if i == excess else 1 for i in range(6))
            words.extend((word, c) for c in product(range(2), repeat=2 + int(matched)))
    ix = {w: i for i, w in enumerate(words)}
    matrix = np.zeros((len(words), len(words)), dtype=int)
    for col, (loads, coins) in enumerate(words):
        matrix[col, col] = 64 * sum(
            (6 - int(board.degree[p])) * (n * n - 1)
            for p, n in zip(ports, loads, strict=True)
        )
        for i, j, label, mark in edges:
            a = raw[label]
            start = LAW.SRC.index((loads[i], loads[j])) + (
                9 * coins[mark] if mark is not None else 0
            )
            for terminal in np.flatnonzero(a[:, start]):
                new_loads, new_coins = list(loads), list(coins)
                new_loads[i], new_loads[j] = LAW.SRC[terminal % 9]
                if mark is not None:
                    new_coins[mark] = terminal // 9
                row = ix[tuple(new_loads), tuple(new_coins)]
                matrix[row, col] += a[terminal, start]
    assert np.array_equal(matrix, matrix.T)
    return matrix, ix


def spectrum(h8, ix):
    coins = min(c for _, c in ix)
    src = [ix[(0, 2, 0, 1, 1, 1), coins], ix[(0, 0, 2, 1, 1, 1), coins]]
    energies, vec = eigh(h8 / 8)
    return energies, vec, src


def branch(spec, ix, time):
    energies, vec, src = spec
    if time == 0:
        evolved = np.eye(len(ix), dtype=complex)[:, src]
    else:
        evolved = (vec * np.exp(-1j * time * energies)) @ vec[src, :].T
    # Axes are the retained head word, output excess address, and input code rail.
    pieces = []
    for coins in sorted({c for _, c in ix}):
        rows = []
        for tip in range(3):
            out = tuple([1] * 3 + [2 if i == tip else 0 for i in range(3)])
            rows.append(ix[out, coins])
        pieces.append(evolved[rows])
    return np.array(pieces)


def records(kraus, rng):
    outputs = np.einsum("cai,cai->ia", kraus.conj(), kraus).real
    assert np.all(outputs.sum(axis=1) <= 1 + 1e-12)
    material = kraus[:, 1:, :]
    bell = material.transpose(0, 2, 1).reshape(len(material), 4) / np.sqrt(2)
    choi = bell.T @ bell.conj()
    partial = choi.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    channel = np.einsum("cai,cbj->abij", material, material.conj())
    projectors = detectors.spin_readers()
    shots = 2**26
    banks = []
    for inp in projectors:
        out = np.einsum("abij,ij->ab", channel, inp)
        settings = []
        for a, b in [(0, 1), (2, 3), (4, 5)]:
            pair = np.array([np.trace(projectors[p] @ out).real for p in (a, b)])
            all_prob = np.maximum(np.r_[pair, 1 - pair.sum()], 0)
            settings.append(rng.multinomial(shots, all_prob / all_prob.sum()).tolist())
        banks.append(settings)
    sharp = np.maximum(np.c_[outputs, 1 - outputs.sum(axis=1)], 0)
    face = detectors.body_turn()
    relative = np.linalg.matrix_power(face, 3)
    dark = (
        sum(np.linalg.norm((np.eye(3) - relative) @ piece) ** 2 for piece in kraus) / 8
    )
    accepted = outputs.sum() / 2
    turn_counts = rng.multinomial(shots, [accepted - dark, dark, 1 - accepted])
    return dict(
        measured_pt_eigenvalues=detectors.measured_partial(banks, shots).tolist(),
        body_full_turn_mixed_counts=turn_counts.tolist(),
        sharp_load_counts=[
            rng.multinomial(shots, row / row.sum()).tolist() for row in sharp
        ],
        sharp_load_outputs=np.c_[outputs, 1 - outputs.sum(axis=1)].tolist(),
        mixed_code_acceptance=float(np.trace(choi).real),
        choi_partial_transpose_eigenvalues=np.linalg.eigvalsh(partial).tolist(),
        pauli_counts=banks,
    )


def run():
    out = {}
    for idx, (name, sign, live, matched) in enumerate(
        [
            ("even", 1, True, False),
            ("odd", -1, True, False),
            ("closed", 1, False, False),
            ("matched", 1, True, True),
            ("matched_odd", -1, True, True),
        ]
    ):
        h, ix = sector(sign, live, matched)
        n, labels = connected_components(h - np.diag(h.diagonal()), directed=False)
        spec = spectrum(h, ix)
        component = labels[spec[2][0]]
        sub = h[np.ix_(labels == component, labels == component)] / 8
        lowest, v = eigh(sub, subset_by_index=(0, min(1, len(sub) - 1)))
        assert min(v[:, 0] * np.sign(v[0, 0])) > 0
        times = TIMES + tuple(2**k for k in range(6, 19)) if matched else TIMES
        out[name] = dict(
            components=int(n),
            occupied_component_size=len(sub),
            occupied_component_lowest=lowest.tolist(),
            reads={
                str(t): records(
                    branch(spec, ix, t), np.random.default_rng([1967, idx, j])
                )
                for j, t in enumerate(times)
            },
        )
    appointments = sum(len(case["reads"]) for case in out.values())
    comparison = {}
    for name in ("even", "matched"):
        odd = out["odd" if name == "even" else "matched_odd"]["reads"]
        comparison[name] = {}
        for j, (t, row) in enumerate(out[name]["reads"].items()):
            # Trace over the orthogonal source parity labels. The full target
            # instrument has zero off-diagonal operator even using its coins.
            p = (
                np.array(row["sharp_load_outputs"][0]) + odd[t]["sharp_load_outputs"][0]
            ) / 2
            p = np.maximum(p, 0)
            rng = np.random.default_rng([1973, j, int(name == "matched")])
            comparison[name][t] = dict(
                opposite_phase_sharp_banks=[
                    rng.multinomial(2**26, p / p.sum()).tolist() for _ in range(2)
                ]
            )
    return dict(
        ports=geometry()[1],
        single_head_coin_herald=detectors.written_source(),
        parity_phase_controls=comparison,
        simultaneous_pt_radius=detectors.simultaneous_bound(2**26, appointments),
        times=TIMES,
        starts_per_pauli_read=2**26,
        compared_source_parity_after_departure="opposite X phases have identical destination records by local parity conservation",
        settings=out,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    r = run()
    path = HERE / "results.json"
    if write:
        path.write_text(json.dumps(r, indent=2) + "\n")
    else:
        old = json.loads(path.read_text())
        for case in r["settings"]:
            for t in r["settings"][case]["reads"]:
                for key in ("sharp_load_outputs", "choi_partial_transpose_eigenvalues"):
                    np.testing.assert_allclose(
                        r["settings"][case]["reads"][str(t)][key],
                        old["settings"][case]["reads"][str(t)][key],
                        atol=1e-11,
                    )
    for t, v in r["settings"]["matched"]["reads"].items():
        print(t, v["mixed_code_acceptance"], v["choi_partial_transpose_eigenvalues"][0])


if __name__ == "__main__":
    main()
