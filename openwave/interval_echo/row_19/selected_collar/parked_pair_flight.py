"""Free motion of two ordinary loads on the full parked collar board."""

from fractions import Fraction as F
from math import factorial, fsum
from os.path import relpath

import numpy as np
from scipy.sparse import coo_matrix, diags

from model.executable_model import carrier
from . import parked_energy_charge as source


class FreePair:
    """Native two-unit sector with all 28,680 terminal port pairs."""

    def __init__(self, board, home):
        a, b = np.triu_indices(board.n, 1)
        self.a, self.b = a, b
        index = np.zeros((board.n, board.n), dtype=int)
        index[a, b] = index[b, a] = np.arange(len(a))
        self.home = index[home]
        self.board = board
        first, second, contact = [], [], []
        self.contact_ids = np.full(len(a), len(board.seams), int)
        for i, (u, v) in enumerate(board.seams):
            other = np.arange(board.n)
            other = other[(other != u) & (other != v)]
            first.extend(index[u, other])
            second.extend(index[v, other])
            contact.extend([i] * len(other))
            self.contact_ids[index[u, v]] = i
        self.first, self.second, self.contact = map(np.array, (first, second, contact))

    def operator(self, marks):
        board = self.board
        marked = np.array([e in marks for e in board.seams])
        edges = np.asarray(board.seams)
        degree = np.bincount(edges[~marked].ravel(), minlength=board.n)
        site = 17 * degree + 16 * (6 - board.degree)
        interactions = np.array(
            [-2 * marks[e] if e in marks else -34 for e in board.seams] + [0]
        )
        binding = np.array([-2 * marks.get(e, 0) for e in board.seams] + [0])
        select = ~marked[self.contact]
        x, y = self.first[select], self.second[select]
        size = len(self.a)
        hop = coo_matrix(
            (-np.ones(2 * len(x)), (np.r_[x, y], np.r_[y, x])), shape=(size, size)
        ).tocsr()
        diagonal = site[self.a] + site[self.b] + interactions[self.contact_ids]
        return hop, diagonal, binding[self.contact_ids]

    def chebyshev(self, operator, shift, weights):
        hop, diagonal = operator
        # Every matrix entry is an exactly representable binary fraction.
        x = (hop + diags(diagonal - shift, format="csr", dtype=float)) / 64
        assert abs(x).sum(axis=1).max() <= 1 and np.diff(x.indptr).max() <= 13
        one = np.zeros(x.shape[0])
        one[self.home] = 1
        prev, curr = one, x @ one
        moments = [1.0, curr[self.home]]
        for _ in range(2, len(weights)):
            prev, curr = curr, 2 * (x @ curr) - prev
            moments.append(curr[self.home])
        # exp(-i t J/2) up to a phase. Separate the real and imaginary terms.
        real = fsum(
            (-1) ** (k // 2) * weights[k] * moments[k]
            for k in range(0, len(weights), 2)
        )
        imag = fsum(
            (-1) ** ((k - 1) // 2) * weights[k] * moments[k]
            for k in range(1, len(weights), 2)
        )
        return real * real + imag * imag


def flight_weights():
    # At t=1 the Chebyshev argument is 32. Bound the tail by the
    # Poisson integral |J_k(32)| <= 16^k/k!, for integer k >= 0.
    order = 84
    tail = 2 * F(16 ** (order + 1), factorial(order + 1)) / (1 - F(16, order + 2))
    assert tail < F(1, 10**25)
    weights = []
    for k in range(order + 1):
        term = F(16**k, factorial(k))
        total = term
        j = 0
        while True:
            j += 1
            term *= -F(256, j * (k + j))
            if abs(term) < F(1, 10**30) and j * (k + j) > 256:
                break
            total += term
        weights.append(float((1 if k == 0 else 2) * total))
    # In the real three-term recurrence the local norm error is below 6e-15
    # (13-term binary64 sparse dot, ||abs(X)||_2 <= 1). The inhomogeneous
    # recurrence propagates one error by U_m(X), whose norm is at most m+1.
    # Thus errors in T_k(X)e are bounded by 6e-15 k(k+1)/2. Include
    # conversion of the exact rational Bessel bounds and final summation.
    amplitude_error = 1e-13 + sum(
        abs(w) * 6e-15 * k * (k + 1) / 2 for k, w in enumerate(weights)
    )
    return weights, 2 * amplitude_error + amplitude_error**2 + 1e-14


def marks_for(board, stations, words, code):
    marks = {}
    for cell in range(board.n // 12):
        if cell in stations:
            word = words[stations.index(cell)]
            matching = code[word + 1]
            digits = divmod(word, 3)
        else:
            matching = code[0]
            digits = (0, 0)
        for k, (a, b) in enumerate(matching):
            marks[12 * cell + a, 12 * cell + b] = (0, 32, 33)[digits[k]] if k < 2 else 0
    return marks


def return_series(board, marks, home, half_order, binding=True):
    # J is twice the sector Hamiltonian, with the initial diagonal removed.
    # Constant written energies and this shift cancel from probabilities.
    blank = [[] for _ in range(board.n)]
    for a, b in board.seams:
        if (a, b) not in marks:
            blank[a].append(b)
            blank[b].append(a)
    homes = [
        17 * len(row) + 16 * (6 - len(board.neighbours[a]))
        for a, row in enumerate(blank)
    ]

    def diagonal(a, b):
        e = marks.get((a, b))
        if e is not None:
            return homes[a] + homes[b] - (2 * e if binding else 0)
        return homes[a] + homes[b] - 34 * int(b in blank[a])

    shift = diagonal(*home)
    even = {home: 1}
    moments = [1]
    columns = {}
    for _ in range(half_order):
        odd = {}
        for (a, b), value in even.items():
            if (a, b) not in columns:
                columns[a, b] = (
                    diagonal(a, b) - shift,
                    [(min(t, b), max(t, b)) for t in blank[a] if t != b]
                    + [(min(t, a), max(t, a)) for t in blank[b] if t != a],
                )
            center, neighbors = columns[a, b]
            odd[a, b] = odd.get((a, b), 0) + center * value
            for key in neighbors:
                odd[key] = odd.get(key, 0) - value
        moments.extend(
            [
                sum(v * even.get(key, 0) for key, v in odd.items()),
                sum(v * v for v in odd.values()),
            ]
        )
        even = odd
    # |<home|exp(-itJ/2)|home>|^2. All odd coefficients vanish by reality.
    return [
        sum(
            F(
                (-1) ** (degree // 2 + k) * moments[k] * moments[degree - k],
                2**degree * factorial(k) * factorial(degree - k),
            )
            for k in range(degree + 1)
        )
        for degree in range(0, len(moments), 2)
    ]


def run(order=6):
    survey = source.reference.local.survey
    board = survey.bulk.federation(0)
    code = survey.routes.codebook()[0]
    # The lexicographically first pair of neighbouring cells, then their
    # first direct bridge. Release one unit load on each end of that bridge.
    left, right, a, b = min((x // 12, y // 12, x, y) for x, y in board.external)
    home = (a, b)
    _, fresh, shared = survey.standards(board)[0][:3]
    _, _, stocks = source.controls(fresh)
    source_record = source.reference.local.inputs("parked_energy_charge")
    np.testing.assert_array_equal(stocks["fresh"], source_record["fresh_counts"])
    w = np.array(source_record["w"])
    stocks["fourth_order"] = stocks["fresh"] + 32 * np.outer(w, w)
    stocks["shared"] = 8 * shared
    weights = {
        k: (v - stocks["fresh"]).ravel() for k, v in stocks.items() if k != "fresh"
    }
    modes = {}
    sector = FreePair(board, home)
    cheb, error = flight_weights()
    held = {}
    for name, binding in (("native", True), ("remove_marked_pair_binding", False)):
        contrasts = {tag: [F(0)] * (order + 1) for tag in weights}
        probabilities = []
        for row, (i, j) in enumerate(np.ndindex(9, 9)):
            marks = marks_for(board, (left, right), (i, j), code)
            coef = return_series(board, marks, home, order, binding)
            hop, diag, bind = sector.operator(marks)
            probabilities.append(
                sector.chebyshev((hop, diag if binding else diag - bind), 138, cheb)
            )
            for tag, vals in weights.items():
                if vals[row]:
                    for n, c in enumerate(coef):
                        contrasts[tag][n] += int(vals[row]) * c
        modes[name] = {
            tag: [str(c / int(stocks["fresh"].sum())) for c in cc]
            for tag, cc in contrasts.items()
        }
        held[name] = {
            "by_source_word_pair": probabilities,
            "changes": {},
        }
        for tag, v in weights.items():
            total = int(stocks["fresh"].sum())
            delta = (
                fsum(int(h) * p for h, p in zip(v, probabilities, strict=True)) / total
            )
            half = error * sum(abs(int(h)) for h in v) / total + 2e-15
            held[name]["changes"][tag] = [delta - half, delta + half]
    hop, diag, _ = sector.operator({})
    ordinary = sector.chebyshev((hop, diag), 186, cheb)
    return {
        "ports": board.n,
        "source_stations": [left, right],
        "release_and_return_ports": [a, b],
        "native_time_even_powers": list(range(0, 2 * order + 1, 2)),
        "ordinary_blank_calibration": [
            str(c) for c in return_series(board, {}, home, order)
        ],
        "change_in_joint_return": modes,
        "flight_time": 1,
        "single_return_probability_error": error,
        "ordinary_unit_time_joint_return": [ordinary - error, ordinary + error],
        "unit_time_joint_return": held,
        "inputs": source.reference.coherent.steady.exp.save_inputs(
            "parked_energy_charge.json",
            relpath(
                carrier.HERE / "carriers" / "FEDERATION_LEVEL0.npz",
                source.reference.local.HERE,
            ),
        ),
    }


if __name__ == "__main__":
    result = run()
    source.reference.coherent.steady.exp.emit("parked_pair_flight", result)
    print({mode: z["changes"] for mode, z in result["unit_time_joint_return"].items()})
