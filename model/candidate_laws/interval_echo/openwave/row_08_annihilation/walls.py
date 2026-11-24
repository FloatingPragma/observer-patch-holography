"""An unmeasured collision of opposite material kinks on a complete native ring."""

import argparse
from itertools import product
import json
import numpy as np
from scipy.linalg import eigh
from model.executable_model import federation
from . import experiment as ex

TIMES = (0, 1, 2, 4, 8, 16, 32, 64, 128)
SHOTS = 2**28


def geometry():
    b = federation(0)
    five = {i for i in range(12) if b.degree[i] == 5}

    def paths(trace):
        if len(trace) == 8:
            if trace[0] in b.neighbours[trace[-1]] and trace[1] < trace[-1]:
                yield trace
            return
        for p in sorted(set(b.neighbours[trace[-1]]) & five - set(trace)):
            if p > trace[0]:
                yield from paths(trace + [p])

    def wings(port):
        for p in sorted(
            q for q in b.neighbours[port] if q // 12 == 0 and b.degree[q] == 6
        ):
            other = next(q for q in b.neighbours[p] if q // 12 != 0)
            for a in sorted(
                q
                for q in b.neighbours[other]
                if q // 12 == other // 12 and b.degree[q] == 5
            ):
                for z in sorted(
                    q
                    for q in b.neighbours[a]
                    if q // 12 == other // 12 and b.degree[q] == 5
                ):
                    yield [p, other, a, z]

    for p in sorted(five):
        for ring in paths([p]):
            for ai, aj in product(range(2), range(2, 4)):
                for left, right in product(wings(ring[ai]), wings(ring[aj])):
                    if set(left).isdisjoint(right):
                        return ring, [ai, aj], [left, right]
    raise ValueError("No native ring has two independent receiving wings")


def reference(inverted):
    return tuple((j + inverted) % 2 for j in range(4))


def loads(bits):
    return [y for i in bits for y in (2 * i, 2 - 2 * i)]


def operators(k, closed, n=2):
    ring, anchors, wings = geometry()
    dim = 2 ** (4 + 2 * n)
    h = np.zeros((dim, dim))
    states = list(product((0, 1), repeat=4 + 2 * n))
    lookup = {s: i for i, s in enumerate(states)}
    t = ex.term(k)
    edges = [(1, 2), (3, 4), (5, 6), (7, 0)]
    for offset in (8, 8 + 2 * n):
        edges.extend((offset + 2 * i + 1, offset + 2 * i + 2) for i in range(n - 1))
    if not closed:
        edges.extend(zip(anchors, (8, 8 + 2 * n)))
    for i, bits in enumerate(states):
        w = loads(bits)
        for j in range(4 + 2 * n):
            v = list(bits)
            v[j] = 1 - v[j]
            a, b = w[2 * j : 2 * j + 2]
            src = ex.LAW.SRC.index((a, b))
            h[lookup[tuple(v)], i] += t[ex.LAW.SRC.index((b, a)), src]
            h[i, i] += t[src, src]
        for a, b in edges:
            src = ex.LAW.SRC.index((w[a], w[b]))
            h[i, i] += t[src, src]
    # Every dimer spans equal-degree ports: its native degree term is scalar in this sector.
    board = federation(0)
    active = ring + wings[0] + wings[1]
    assert all(
        board.degree[a] == board.degree[b] for a, b in zip(active[::2], active[1::2])
    )
    return h


def sample(prob, rng):
    p = np.maximum(np.array(prob), 0).round(12)
    order = np.argsort(p, kind="stable")
    counts = np.empty(len(p), np.int64)
    q = p[order]
    q[-1] = 1 - q[:-1].sum()
    counts[order] = rng.multinomial(SHOTS, q)
    return dict(exact=p.tolist(), counts=counts.tolist())


def wing_ref(k, q, cut, n=2):
    h = ex.wing(k, n, q)
    if cut:
        # Remove the source-link diagonal exactly, keeping all intra-wing terms.
        t = ex.term(k)
        adj = [2 * b[0] for b in product((0, 1), repeat=n)]
        h -= np.diag(
            [t[ex.LAW.SRC.index((q, x)), ex.LAW.SRC.index((q, x))] for x in adj]
        )
    e, v = eigh(h)
    return h, e[0], v[:, 0]


def experiment(k, cut, rng):
    _, ports, _ = geometry()
    h = operators(k, cut)
    e, v = eigh(h)
    rows = []
    for inv in (0, 1):
        ref = reference(inv)
        bits = list(product((0, 1), repeat=4))
        order = [bits.index(ref), bits.index(tuple(1 - b for b in ref))]
        q = np.array(loads(ref))[ports]
        init = [wing_ref(k, b, cut) for b in q]
        flat = [ex.wing(k, 2, 1) for _ in q]
        for wall in (False, True):
            seq = tuple(z ^ int(wall and i >= 2) for i, z in enumerate(ref))
            x = np.zeros((16, 4, 4))
            x[bits.index(seq)] = np.outer(init[0][2], init[1][2])
            psi = x.ravel()
            energy0 = float(psi @ h @ psi)
            eigen = v.T @ psi
            clocks = []
            for time in TIMES:
                w = (v @ (np.exp(-1j * time * (e - e[0])) * eigen)).reshape(16, 4, 4)
                assert abs(np.vdot(w.ravel(), h @ w.ravel()).real - energy0) < 1e-8
                found = w[order]
                p = np.zeros(4)
                fixed = np.zeros(4)
                for i, a in enumerate(found):
                    local = [wing_ref(k, b if i == 0 else 2 - b, cut) for b in q]
                    bright = [np.eye(4) - np.outer(t[2], t[2]) for t in local]
                    fixed_bright = [np.eye(4) - np.outer(t[2], t[2]) for t in init]
                    p += ex.read_joint(a, *bright)
                    fixed += ex.read_joint(a, *fixed_bright)
                miss = 1 - p.sum()
                end = ex.read_joint(found, ex.effects(2), ex.effects(2))

                def means(arr, ops):
                    return [
                        np.vdot(arr, ops[0] @ arr).real,
                        np.vdot(arr, arr @ ops[1].T).real,
                    ]

                wing_delta = np.array(means(w, flat)) - means(x, flat)
                conditional = []
                # Energy of both wings measured in the original calibration, including all heralds.
                for a in found:
                    conditional.append(
                        means(a, [t[0] - t[1] * np.eye(4) for t in init])
                    )
                clocks.append(
                    dict(
                        time=time,
                        wing_energy_change=wing_delta.tolist(),
                        ordered_wing_energy_per_launch=conditional,
                        wing_energies=sample(np.r_[p, miss], rng),
                        same_reference=sample(np.r_[fixed, miss], rng),
                        end_parities=sample(np.r_[end, miss], rng),
                        ring_wall_count=sample(
                            [
                                sum(
                                    abs(w[j]).ravel() @ abs(w[j]).ravel()
                                    for j, b in enumerate(bits)
                                    if sum(b[i] == b[(i + 1) % 4] for i in range(4))
                                    == n
                                )
                                for n in (0, 2, 4)
                            ],
                            rng,
                        ),
                    )
                )
            rows.append(
                dict(invert=inv, wall_pair=wall, initial_energy=energy0, clocks=clocks)
            )
    return dict(contrast=k, cut_wings=cut, rows=rows)


def run():
    rng = np.random.default_rng(81128)
    ring, ports, wings = geometry()
    return dict(
        shots=SHOTS,
        seed=81128,
        ring=ring,
        attachments=ports,
        wings=wings,
        cases=[
            experiment(k, cut, rng) for k, cut in ((8, False), (8, True), (0, False))
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    d = json.loads(json.dumps(run()))
    p = ex.HERE / "walls.json"
    if args.write:
        p.write_text(json.dumps(d, indent=1, allow_nan=False) + "\n")
    else:
        from .verify import compare

        compare(d, json.loads(p.read_text()), "walls")
        print(
            "Verified closed-ring disappearance, two wing spectra, remote parity and controls."
        )


if __name__ == "__main__":
    main()
