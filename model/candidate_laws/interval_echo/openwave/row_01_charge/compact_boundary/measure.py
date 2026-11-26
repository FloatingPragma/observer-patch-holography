"""Remote scattering from a compact signed load in the registered two-owner aperture."""

import argparse
import json
import math
from functools import cache
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from .. import measure as original

HERE = Path(__file__).resolve().parent
TIMES = (0, 2, 4, 6, 8)
STARTS = 2**22
SEED = 26112701
SIGNS = (-1, 0, 1)
PORTS = tuple(range(12)) + tuple(range(72, 84))
BRIDGE = (1, 82)


@cache
def geometry():
    board = federation(0)
    ext = nx.Graph(board.seams).subgraph(PORTS[12:])
    distances = nx.single_source_shortest_path_length(ext, 82)
    detector = min(
        (p for p in ext if board.degree[p] == 5), key=lambda p: (-distances[p], p)
    )
    shared = [
        tuple(sorted(z))
        for z in board.seams
        if set(z) <= set(PORTS) and z[0] // 12 != z[1] // 12
    ]
    if shared != [BRIDGE]:
        raise ValueError("unexpected two-owner aperture")
    _, spec = original.source_and_reader(board, 0)
    mode = np.array(spec["integer_amplitudes_in_cell"], dtype=int)
    if np.any(mode[spec["gateway_ports"]]):
        raise ValueError("source reaches the aperture")
    return board, detector, mode


def word(stock):
    return tuple(sorted((p, s) for p, s in stock if s))


@cache
def sector(signs, strength, opened):
    """Use literal native diagonal terms on all incident seams and private ports."""
    board, _, _ = geometry()
    if len(signs) == 0:
        words = [()]
    elif len(signs) == 1:
        words = [((p, signs[0]),) for p in PORTS]
    else:
        pairs = (
            combinations(PORTS, 2) if signs[0] == signs[1] else permutations(PORTS, 2)
        )
        words = [word(zip(pair, signs, strict=True)) for pair in pairs]
    idx = {w: i for i, w in enumerate(words)}
    if len(idx) != len(words):
        raise ValueError("duplicate defect word")
    block = blocks(strength)["blank"]
    host = LAW.SRC.index((1, 1))
    operator = sparse.dok_matrix((len(words), len(words)), dtype=np.int64)
    for col, w in enumerate(words):
        stock = dict(w)
        diagonal = (
            8
            * strength
            * sum((6 - int(board.degree[p])) * ((1 + s) ** 2 - 1) for p, s in w)
        )
        seams = {tuple(sorted((p, q))) for p, s in w for q in board.neighbours[p]}
        for a, b in seams:
            loads = (1 + stock.get(a, 0), 1 + stock.get(b, 0))
            j = LAW.SRC.index(loads)
            diagonal += block[j, j] - block[host, host]
            if (
                a in PORTS
                and b in PORTS
                and loads[0] != loads[1]
                and (opened or (a, b) != BRIDGE)
            ):
                switched = {a: b, b: a}
                moved = word((switched.get(p, p), s) for p, s in w)
                operator[idx[moved], col] += block[LAW.SRC.index(loads[::-1]), j]
        operator[col, col] += diagonal
    h8 = operator.tocsr()
    if (h8 != h8.T).nnz:
        raise ValueError("native restriction must be symmetric")
    return idx, h8


def compact_energy(sign, strength):
    board, _, local = geometry()
    source = np.zeros(board.n, dtype=int)
    source[:12] = local
    full = defect_generator(board, sign, strength, 0)
    # Half-integer generator with integer orbitals makes this check exact.
    op = full * 8
    np.testing.assert_array_equal(op.data, op.data.round())
    evolved = op.astype(int) @ source
    first = np.flatnonzero(source)[0]
    energy8 = int(evolved[first] // source[first])
    np.testing.assert_array_equal(evolved, energy8 * source)
    return energy8


def closed_certificate(source_sign, probe_sign, strength):
    _, _, source = geometry()
    singles, one = sector((probe_sign,), strength, False)
    order = [singles[((p, probe_sign),)] for p in PORTS[12:]]
    indices, h8 = sector(tuple(sorted((source_sign, probe_sign))), strength, False)
    embedding = np.zeros((len(indices), 12), dtype=int)
    for j, p in enumerate(PORTS[12:]):
        for x in np.flatnonzero(source):
            embedding[indices[word(((x, source_sign), (p, probe_sign)))], j] = source[x]
    exterior = one[order, :][:, order].toarray()
    rhs = embedding @ (
        exterior + compact_energy(source_sign, strength) * np.eye(12, dtype=int)
    )
    np.testing.assert_array_equal(h8 @ embedding, rhs)
    np.testing.assert_array_equal(
        embedding.T @ embedding, np.dot(source, source) * np.eye(12, dtype=int)
    )


def initial(idx, source_sign, probe_sign):
    mode = geometry()[2]
    first = np.zeros(len(idx), dtype=float)
    for x in np.flatnonzero(mode) if source_sign else [0]:
        config = word(((x, source_sign), (82, probe_sign)))
        first[idx[config]] = mode[x] / np.linalg.norm(mode) if source_sign else 1
    return first


def outcome_bins(idx):
    detector = geometry()[1]
    bins = []
    for w in idx:
        loads = dict(w)
        count = sum(p < 12 for p, s in w)
        bins.append(9 * count + 3 * (1 + loads.get(82, 0)) + 1 + loads.get(detector, 0))
    return np.asarray(bins)


@cache
def probabilities(source_sign, probe_sign, strength, opened):
    idx, h8 = sector(
        tuple(sorted(s for s in (source_sign, probe_sign) if s)), strength, opened
    )
    first = initial(idx, source_sign, probe_sign)
    if h8.nnz:
        op = -0.125j * h8
        trajectory = expm_multiply(
            op, first, start=0, stop=8, num=5, traceA=op.diagonal().sum()
        )
    else:
        trajectory = np.tile(first, (5, 1))
    bins = outcome_bins(idx)
    panel = []
    for psi in trajectory:
        rates = np.bincount(bins, weights=abs(psi) ** 2, minlength=27)
        np.testing.assert_allclose(rates.sum(), 1, rtol=0, atol=1e-10)
        panel.append(rates / rates.sum())
    return panel


def enrollment():
    return product((int(candidate.KAPPA), 0), (False, True), SIGNS, SIGNS, TIMES)


def acquire():
    for strength, ss, ps in product((int(candidate.KAPPA), 0), (-1, 1), (-1, 1)):
        closed_certificate(ss, ps, strength)
    for n in range(3):
        assert original.sharp_contrast(n, 0) == n
    rows = []
    for strength, opened, ss, ps, t in enrollment():
        rates = probabilities(ss, ps, strength, opened)[TIMES.index(t)]
        key = [strength, opened, ss, ps, t]
        digest = sha256(json.dumps([SEED, *key]).encode()).digest()
        gen = np.random.default_rng(int.from_bytes(digest[:16], "big"))
        rows.append([*key, gen.multinomial(STARTS, rates).tolist()])
    return dict(seed=SEED, starts=STARTS, ports=[82, int(geometry()[1])], banks=rows)


def decode(raw):
    if (
        raw["seed"] != SEED
        or raw["starts"] != STARTS
        or raw["ports"] != [82, geometry()[1]]
    ):
        raise ValueError("wrong enrollment")
    data = {}
    for *key, counts in raw["banks"]:
        if (
            len(counts) != 27
            or any(type(c) is not int or c < 0 for c in counts)
            or sum(counts) != STARTS
        ):
            raise ValueError("invalid terminal bank")
        if tuple(key) in data:
            raise ValueError("duplicate terminal bank")
        data[tuple(key)] = np.array(counts).reshape(3, 3, 3) / STARTS
    if set(data) != set(enrollment()):
        raise ValueError("incomplete terminal history")
    eps = math.sqrt(math.log(2 * 32768 / 0.01) / (2 * STARTS))
    studies = []
    for strength, opened in product((int(candidate.KAPPA), 0), (False, True)):
        examples = []
        for a, b in combinations(SIGNS, 2):
            diffs = []
            for ps, t in product((-1, 1), TIMES):
                delta = (
                    data[strength, opened, a, ps, t] - data[strength, opened, b, ps, t]
                ).sum(axis=0)
                i, j = np.unravel_index(abs(delta).argmax(), delta.shape)
                diffs.append(
                    (float(abs(delta[i, j])), ps, t, int(i), int(j), float(delta[i, j]))
                )
            value, ps, t, i, j, signed = max(diffs)
            examples.append(
                dict(
                    sources=[a, b],
                    probe=ps,
                    time=t,
                    word=[i, j],
                    difference=signed,
                    resolved=value > 2 * eps,
                )
            )
        odd = []
        for t in TIMES:

            def returned(ss, ps):
                return data[strength, opened, ss, ps, t][:, 1 + ps, :].sum()

            odd.append(
                float(
                    (returned(1, 1) - returned(-1, 1))
                    - (returned(1, -1) - returned(-1, -1))
                )
            )
        studies.append(
            dict(
                strength=strength,
                bridge="open" if opened else "parked",
                exterior=examples,
                signed_return_difference=odd,
            )
        )
    return dict(half_width=eps, studies=studies)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire", type=Path)
    args = cli.parse_args()
    if args.acquire:
        raw = acquire()
        args.acquire.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
    else:
        raw = json.loads((HERE / "counts.json").read_text())
    print(json.dumps(decode(raw), indent=2))


if __name__ == "__main__":
    main()
