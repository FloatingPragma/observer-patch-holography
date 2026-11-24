"""Stationary exterior Ramsey probes with all source-to-gateway exchanges open."""

import argparse
from functools import lru_cache
from itertools import combinations
import json

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1 import law

from . import experiment as base


@lru_cache(maxsize=None)
def chamber(n, strength, reference=1):
    geometry = base.geometry()
    inner = list(base.belt.ring(10).seams)
    open_edges = inner + [(i, 10) for i in geometry["near_indices"]]
    # Each belt site sees one polar contact and one supplied host-1 reference.
    # The five contacts to the held opposite pole and the far-to-gateway contact
    # contribute their exact diagonals; each uses an actual native seam.
    fixed = [i for i in range(10) if i not in geometry["near_indices"]] + [10]
    words = list(combinations(range(11), n))
    number = {word: j for j, word in enumerate(words)}
    seam = base.free_seam(strength)
    entries, rows, cols = [], [], []
    second = np.empty(len(words), dtype=int)
    belt_labels = base.detector(n)
    for j, word in enumerate(words):
        occupied = set(word)
        diag = (
            strength * (3 - 2 * reference) * len(occupied - {10})
        )  # unused sixth reference of belt ports
        for a, b in open_edges:
            states = (1 + (a in occupied), 1 + (b in occupied))
            loc = law.SRC.index(states)
            diag += seam[loc, loc]
            if states[0] != states[1]:
                target = tuple(sorted(occupied ^ {a, b}))
                exchange = seam[law.SRC.index(states[::-1]), loc]
                entries.append(exchange)
                rows.append(number[target])
                cols.append(j)
        for a in fixed:
            loc = law.SRC.index((1 + (a in occupied), 1))
            diag += seam[loc, loc]
        entries.append(diag)
        rows.append(j)
        cols.append(j)
        second[j] = 11 if 10 in occupied else belt_labels[base.sector(n).index[word]]
    h = sparse.csr_matrix((entries, (rows, cols)), shape=(len(words), len(words)))
    # Numerical common offset does not alter Ramsey or position observables.
    return h - sparse.eye(len(words)) * h.diagonal().mean(), number, second


def outcomes(ref, sample, labels, rng):
    frames = []
    for angle in (1, -1j):
        prob = [
            np.bincount(
                labels, weights=abs((ref + s * angle * sample) / 2) ** 2, minlength=12
            )
            for s in (1, -1)
        ]
        frames.append(base.record(prob, rng))
    v = np.vdot(ref, sample)
    return dict(
        visibility=base.rounded([v.real, v.imag]),
        quadratures=frames,
        gateway_probabilities=base.rounded(
            [np.sum(abs(s[labels == 11]) ** 2) for s in (ref, sample)]
        ),
    )


def run(reference=1):
    seed = 1211250 - reference
    rng = np.random.default_rng(seed)
    rows = []
    ref_bulk = base.sector(1)
    triangle = (4, 5, 6)
    assert set(combinations(triangle, 2)).issubset(set(ref_bulk.car.seams))
    triangle_h = np.array([[1, -0.5, -0.5], [-0.5, 1, -0.5], [-0.5, -0.5, 1]])
    bulk = base.record(
        abs(base.evolve(sparse.csr_matrix(triangle_h), np.eye(3)[:, 0], 1.0)) ** 2, rng
    )
    for strength in (8, 0):
        for n, kind, shift, sense in [*base.inputs(), (1, "gateway_standard", 0, 1)]:
            h, number, labels = chamber(n, strength, reference)
            gateway = (labels == 11).astype(float)
            state = np.zeros(len(number), complex)
            if kind == "gateway_standard":
                state[number[(10,)]] = 1
            else:
                source_h = sparse.csr_matrix(base.sector(n).hamiltonian(strength))
                prep = base.preparation(kind, n, source_h, shift, sense)
                for k, word in enumerate(base.sector(n).states):
                    state[number[word]] = prep[k]
            # Guard only the entrance exchanges; retain their native diagonal.
            g = sparse.diags(gateway)
            background = sparse.eye(len(number)) - g
            held = (g @ h @ g + background @ h @ background).tocsr()
            row = dict(
                kappa=strength,
                loads=n,
                kind=kind,
                shift=shift,
                sense=sense,
                settings=[],
            )
            for t in base.LONG:
                for opened, free in [(False, held), (True, h)]:
                    ref = base.evolve(free, state, t)
                    for sign in (-1, 1):
                        _, gain, _ = base.contact(strength, sign)
                        sample = base.evolve(free + gain * g, state, t)
                        row["settings"].append(
                            dict(
                                time=t,
                                open=opened,
                                sign=sign,
                                **outcomes(ref, sample, labels, rng),
                            )
                        )
            rows.append(row)
    return dict(
        seed=seed,
        spare_slot_reference=reference,
        shots=base.SHOTS,
        belt_reference_indices=list(triangle),
        bulk_pointer_words=[[2, 1, 1], [1, 2, 1], [1, 1, 2]],
        bulk=bulk,
        rows=rows,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    for reference, name in [(1, "exterior.json"), (0, "exterior_native.json")]:
        path = base.HERE / name
        data = run(reference)
        if args.write:
            path.write_text(
                json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"
            )
        elif data != json.loads(path.read_text()):
            raise SystemExit(f"Exterior probe replay differs from {name}")
        print(
            f"{len(data['rows'])} exterior preparations verified for reference {reference}"
        )


if __name__ == "__main__":
    main()
