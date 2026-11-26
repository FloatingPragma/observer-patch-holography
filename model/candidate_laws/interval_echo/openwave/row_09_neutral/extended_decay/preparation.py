"""Commission a native real-amplitude cold stock on the five crossed dimers."""

from collections import Counter, deque
import json
import argparse
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from ..passive_decay.readers import rotate, Z, X
from ..passive_decay.material import edge
from ..passive_decay.measure import compare_records, hamiltonians
from ..joined_decay.measure import core_pulse
from ..passive_decay.readers import circuit
from .material import route
from .measure import HERE, rebuild


def mount(board, dimers, heads):
    # Select a single real interacting contact for each edge. Primary word
    # order is always 00,01,10,11, independent of contact orientation.
    sites = {p: (i, k) for i, d in enumerate(dimers) for k, p in enumerate(d)}
    links = {}
    for (p, q), cost in sorted(heads.items()):
        if p not in sites or q not in sites:
            continue
        a, ia = sites[p]
        b, ib = sites[q]
        if a == b:
            continue
        h = np.diag(
            [cost * (1 - ((s ^ ia) * (t ^ ib))) for s in (0, 1) for t in (0, 1)]
        )
        cz = -np.diag(np.exp(-1j * np.pi * np.diag(h) / cost))
        # Flip either complementary rail into the logical primary convention.
        orient = np.kron(X if ia else np.eye(2), X if ib else np.eye(2))
        cz = np.real_if_close(orient @ cz @ orient)
        np.testing.assert_allclose(cz, np.diag([1, 1, 1, -1]), atol=1e-12)
        links.setdefault(edge(a, b), ((p, q), cz))
    private_or_pin = []
    for d in dimers:
        if any(board.degree[p] == 5 for p in d):
            private_or_pin.append("private")
        else:
            assert any(
                (p in d and q not in sites) or (q in d and p not in sites)
                for p, q in heads
            )
            private_or_pin.append("exterior unit pin")
    return links, private_or_pin


def synthesize(target, links):
    n = int(np.log2(len(target)))
    state = np.zeros(2**n)
    state[0] = 1
    tally = Counter()
    had = rotate(np.pi / 4) @ Z

    def gate(u, sites):
        a = np.moveaxis(state.reshape((2,) * n), sites, range(len(sites)))
        shape = a.shape
        b = (u @ a.reshape(2 ** len(sites), -1)).reshape(shape)
        state[:] = np.moveaxis(b, range(len(sites)), sites).ravel()

    def local_cn(a, b):
        ports, cz = links[edge(a, b)]
        frame = np.kron(np.eye(2), had)
        gate(frame @ cz @ frame, (a, b))
        tally[f"conditional contact {ports}"] += 1

    def swap(a, b):
        local_cn(a, b)
        local_cn(b, a)
        local_cn(a, b)

    def routed(a, b):
        routes = deque([[a]])
        while routes:
            path = routes.popleft()
            if path[-1] == b:
                break
            for q in range(n):
                if q not in path and edge(path[-1], q) in links:
                    routes.append([*path, q])
        else:
            raise ValueError((a, b))
        moves = list(zip(path[:-2], path[1:-1]))
        for x, y in moves:
            swap(x, y)
        local_cn(path[-2], path[-1])
        for x, y in reversed(moves):
            swap(x, y)

    for depth in range(n):
        size = 2**depth
        block = target.reshape(size, 2, -1)
        weights = block[:, :, 0] if depth == n - 1 else np.linalg.norm(block, axis=2)
        angles = np.arctan2(weights[:, 1], weights[:, 0])
        gray = [j ^ (j >> 1) for j in range(size)]
        for j, code in enumerate(gray):
            angle = (
                sum(
                    ((-1) ** ((word & code).bit_count())) * angles[word]
                    for word in range(size)
                )
                / size
            )
            gate(rotate(angle), (depth,))
            tally["single-dimer Y rotation"] += 1
            change = code ^ gray[(j + 1) % size]
            if change:
                routed(depth - change.bit_length(), depth)
    error = float(np.linalg.norm(state - target))
    assert error < 2e-12
    return error, dict(tally)


def run():
    b, d, h, _, _ = route((3,), True)
    s, *_ = hamiltonians(4)
    _, core = eigh(s)
    flip = core_pulse(core)
    np.testing.assert_allclose(circuit(flip, 32), flip, atol=1e-12)
    runs = []
    for cut in (False, True):
        material = dict(h)
        if cut:
            material[edge(d[-2][0], d[-1][0])] = 32
        links, phase = mount(b, d, material)
        _, v = eigh(rebuild(b, d, h, cut, True))
        err, tally = synthesize(v[:, 0], links)
        runs.append(
            dict(
                interior_cut=cut,
                rows=len(v),
                error_below_2e_12=bool(err < 2e-12),
                native_gate_counts=tally,
                diagonal_phase=phase,
            )
        )
    return dict(
        boards=runs, local_pulse="real two-source-dimer reader, outer-level reflection"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    dest = parser.parse_args().output
    r = run()
    if dest:
        dest.write_text(json.dumps(r, indent=1) + "\n")
    else:
        compare_records(json.loads((HERE / "preparation.json").read_text()), r)
    print(
        "Verified preparation of both complete 32-word cold stocks and the local core pulse."
    )
