"""Build each connector from a sharp native write and a two-port return test."""

from collections import Counter
from hashlib import sha256

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from .. import experiment, formation

SEED = 723091
LIMIT = 24
MATERIAL_SEED = 312141


def instrument():
    """Terminal-load Kraus maps, directly from the native seam on both arms."""
    law = experiment.candidate.native_law
    exchange, contrast = law.local_term_parts()
    full = exchange + experiment.K * contrast
    echo = law.echo()
    base = np.eye(law.DIM)
    loads = np.eye(law.N_SRC)

    def word(a, b):
        return loads[:, law.SRC.index((a, b))]

    # The borrowed hole and probe carry one unit; blank one-unit pulses split,
    # recombine and, if necessary, swap it back to the connector.
    dim = base[:, [law.SRC.index((1, 0)), law.SRC.index((0, 1))]]
    split = dim.T @ expm(-1j * np.pi / 2 * full) @ dim
    undo = dim.T @ expm(-1j * 3 * np.pi / 2 * full) @ dim
    swap = dim.T @ expm(-1j * np.pi * full) @ dim
    np.testing.assert_allclose(up_to_phase(undo @ split), np.eye(2), atol=1e-12)
    np.testing.assert_allclose(up_to_phase(swap), [[0, 1], [1, 0]], atol=1e-12)
    wait = expm(-1j * np.pi * full)
    donors = [base[:, law.SRC.index(w)] for w in ((2, 0), (0, 2))]
    marks = []
    energies = []
    for parity in (1, -1):
        expected = np.kron(label(2, parity), word(1, 1))
        prepared = echo @ ((donors[0] + parity * donors[1]) / np.sqrt(2))
        np.testing.assert_allclose(
            up_to_phase(prepared), up_to_phase(expected), atol=1e-13
        )
        mark = prepared.reshape(-1, law.N_SRC) @ word(1, 1)
        marks.append(mark)
        columns = np.kron(mark[:, None], loads)
        block = columns.T @ full @ columns
        np.testing.assert_allclose(full @ columns, columns @ block, atol=3e-13)
        gap = float(block[law.SRC.index((0, 1)), law.SRC.index((0, 1))])
        np.testing.assert_allclose(
            block, np.diag([0 if w == (1, 1) else gap for w in law.SRC]), atol=3e-13
        )
        energies.append(gap)

    def read(mark, unitary):
        # The complementary auxiliary load labels these two paths.
        branches = []
        for i, a in enumerate((1, 0)):
            inc = np.kron(mark, word(a, 1))
            after = unitary @ inc
            col = after.reshape(-1, law.N_SRC) @ word(a, 1)
            np.testing.assert_allclose(after, np.kron(col, word(a, 1)), atol=2e-13)
            branches.append(split[i, 0] * col)
        return undo @ np.asarray(branches)

    sharp = [(echo @ q).reshape(-1, law.N_SRC) @ word(1, 1) for q in donors]
    transfers = []
    for initial in sharp:
        output = read(initial, wait)
        pmf = [float(np.vdot(x, x).real) for x in output]
        np.testing.assert_allclose(pmf, [0.5, 0.5], atol=2e-13)
        parked = read(initial, base)
        np.testing.assert_allclose(abs(parked[1]), 0, atol=2e-13)
        for bit, state in enumerate(output):
            normalized = state / np.sqrt(pmf[bit])
            np.testing.assert_allclose(
                up_to_phase(normalized), up_to_phase(marks[bit]), atol=2e-13
            )
            second = read(normalized, wait)
            np.testing.assert_allclose(np.linalg.norm(second[1 - bit]), 0, atol=3e-13)
            back = echo @ np.kron(normalized, word(1, 1))
            weights = np.abs(np.array([q @ back for q in donors])) ** 2
            np.testing.assert_allclose(weights, [0.5, 0.5], atol=2e-13)
        transfers.append([round(x, 12) for x in pmf])
    return dict(
        sharp_pointer_laws=transfers,
        repeat_pointer_error=0,
        parked_moved_probability=0,
        contact_gaps=[round(x, 12) for x in energies],
    )


def up_to_phase(col):
    flat = col.ravel()
    idx = np.flatnonzero(abs(flat) > 1e-9)[0]
    return col * np.conjugate(flat[idx]) / abs(flat[idx])


def commission(level, standard):
    board = federation(level)
    components = experiment.grow(board, MATERIAL_SEED)
    route = max((s for cycle, s in components if not cycle), key=len)
    links = sorted(
        {
            formation.edge(left[1], right[0])
            for cycle, strand in components
            for left, right in zip(strand, strand[1:] + (strand[:1] if cycle else []))
        }
    )
    blocked = set(links)
    allowed = [
        [
            int(v)
            for v in sorted(board.neighbours[i])
            if formation.edge(i, v) not in blocked
        ]
        for i in range(board.n)
    ]
    used = {v for link in links for v in link}
    holes = sorted(set(range(board.n)) - used)[: len(route)]
    assert len(holes) == len(route)
    permanent = holes[0]
    stock = np.ones(board.n, dtype=int)
    stock[holes] = 0
    for a, b in links:
        stock[a], stock[b] = 2, 0
    rng = np.random.default_rng([SEED, level, standard, 1])
    trace, routes, max_swaps = [], sha256(), 0
    for index, (a, b) in enumerate(links):
        # Search for a blank path ending on a neighbour of a; avoid disturbing either donor load.
        safe = [
            neighbours if i not in (a, b) else []
            for i, neighbours in enumerate(allowed)
        ]
        targets = set(allowed[a]) - {b}
        path = formation.path_to(permanent, safe, lambda port: port in targets)
        v = path[-1]
        walk = list(zip(path, path[1:]))

        def exchange(x, y):
            stock[x], stock[y] = stock[y], stock[x]

        for x, y in walk:
            exchange(x, y)
        assert stock[v] == 0 and (stock[a], stock[b]) == (2, 0)
        routes.update(f"{a},{b},{path}\n".encode())
        max_swaps = max(max_swaps, len(walk))
        # Exact conditional law of the first returned pointer defines this cohort.
        attempts, resets = [], []
        for attempt in range(LIMIT):
            bit = standard if index == 0 else int(rng.integers(2))
            attempts.append(bit)
            if bit == standard:
                break
            # Inverse echo on restored 11, project onto 20 or 02, then reuse that word.
            resets.append(int(rng.integers(2)))
        for x, y in reversed(walk):
            exchange(x, y)
        if attempts[-1] != standard:
            return dict(
                complete=False,
                failed_link=index,
                reads=trace + [dict(pointer=attempts, reset=resets)],
            )
        stock[a], stock[b] = 1, 1
        trace.append(dict(pointer=attempts, reset=resets))
    assert set(np.flatnonzero(stock == 0)) == set(holes)
    assert np.all(stock[stock != 0] == 1) and sum(stock) == board.n - len(route)
    layouts = formation.prepare(level, MATERIAL_SEED)
    assert layouts["load_total"] == int(stock.sum())
    n = len(links)
    return dict(
        complete=True,
        cells=board.cells,
        enrolled_route=route,
        connectors=n,
        vacancy_route_sha256=routes.hexdigest(),
        max_one_way_swaps=max_swaps,
        stock_total=int(stock.sum()),
        attempt_histogram=dict(
            sorted(Counter(len(r["pointer"]) for r in trace).items())
        ),
        reads=trace,
        completion_probability=(1 - 2**-LIMIT) ** (n - 1),
        uniform_without_feedback=2.0 ** (2 - len(route)),
    )
