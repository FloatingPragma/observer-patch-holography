"""Total cell-translation spectrum and the physical cost of restoring winding heads."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.linalg import eigh, eigvalsh
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components

from . import geometry, reservoir

HERE = reservoir.HERE
PROTOCOL = "647593dbf1e96fc613e9517fc1173f8036250cd2e51de729c31cdff449e7456e"
SHOTS = 2**22


def rotate(v, shift, n):
    return ((v << shift) | (v >> (n - shift))) & ((1 << n) - 1)


def charges(words, n):
    L = n // 2
    return np.column_stack(
        [
            ((words >> (2 * i)) & 1)
            + ((words >> (2 * i + 1)) & 1)
            + ((words >> (n + i)) & 1)
            - ((words >> (n + (i - 1) % L)) & 1)
            for i in range(L)
        ]
    )


def blocks(h, shift, cells):
    rows, columns = h.nonzero()
    assert not (h[shift, :][:, shift] - h).nnz
    cycles = []
    seen = np.zeros(len(shift), dtype=bool)
    for start in range(len(shift)):
        if seen[start]:
            continue
        cycle = [start]
        x = shift[start]
        while x != start:
            cycle.append(x)
            x = shift[x]
        seen[cycle] = True
        cycles.append(cycle)
    spectra = []
    for j in range(cells):
        r, c, z = [], [], []
        for cycle in cycles:
            length = len(cycle)
            if j * length % cells:
                continue
            col = 0 if not c else c[-1] + 1
            r.extend(cycle)
            c.extend([col] * length)
            z.extend(
                np.exp(-2j * np.pi * j * np.arange(length) / cells) / np.sqrt(length)
            )
        W = coo_matrix((z, (r, c)), shape=(len(shift), max(c) + 1)).tocsr()
        K = (W.getH() @ h @ W).tocsr()
        K.data[abs(K.data) < 1e-11] = 0
        K.eliminate_zeros()
        nr, tags = connected_components(K, directed=False)
        pieces = []
        for group in range(nr):
            ix = np.flatnonzero(tags == group)
            pieces.extend(eigvalsh(K[np.ix_(ix, ix)].toarray() / 8))
        spectra.append(np.sort(pieces))
    assert sum(map(len, spectra)) == len(shift)
    return spectra


def spectral_result(spectra):
    compar = [max(abs(s - spectra[0])) for s in spectra if len(s) == len(spectra[0])]
    return dict(
        dimensions=[len(s) for s in spectra],
        ground_by_momentum=[s[0] for s in spectra],
        maximum_by_momentum=[s[-1] for s in spectra],
        all_same_size=len(set(map(len, spectra))) == 1,
        max_difference_at_equal_dimension=max(compar),
    )


def blank(n):
    native = reservoir.native_contacts()["blank"]
    # Every single vacancy touches two unlike blank bonds.
    h = np.zeros((n, n), dtype=int)
    for i in range(n):
        h[i, i] = 2 * native[2, 2]
        h[i, (i + 1) % n] = h[(i + 1) % n, i] = native[1, 2]
    return h


def mode_comparison(n, tables):
    L = n // 2
    rows = []
    for N in sorted({1, L - 1, L, L + 1}):
        words, h = reservoir.integer_sector(n, N, tables)
        G = charges(words, n)
        r, c = h.nonzero()
        assert np.array_equal(G[r], G[c])
        assert np.all(G.sum(axis=1) == N)
        shifted = rotate(words & ((1 << n) - 1), 2, n) + (rotate(words >> n, 1, L) << n)
        shift = np.searchsorted(words, shifted)
        np.testing.assert_array_equal(words[shift], shifted)
        spectra = blocks(h, shift, L)
        rows.append(
            dict(
                number=N,
                distinct_local_charge_vectors=len(np.unique(G, axis=0)),
                integer_nonzero_entries=len(r),
                **spectral_result(spectra),
            )
        )
    h = blank(n)
    alt = spectral_result(blocks(csr_matrix(h), np.roll(np.arange(n), -2), L))
    return dict(n=n, sectors=rows, all_blank_one_vacancy=alt)


def propagator():
    flights = json.loads((HERE / "geometry.json").read_text())["flights"]
    route = max((r["route"] for r in flights if r["level"] == 2), key=len)
    n = len(route)
    assert n == 82
    h, ports = geometry.microscopic(route)
    # Replay each native written hop on the accessible finite path. A blank
    # changes the vacancy position while leaving the complete history unchanged.
    head = 0
    words = [head]
    for step in range(n + 1):
        if step % 2:
            head ^= 1 << ((step % n) // 2)
        words.append(head)
    runs = []
    rng = np.random.default_rng(317144)
    for name, op, where, history in [
        ("winding", h, ports, words),
        ("blank", blank(n), np.arange(n), [0] * n),
    ]:
        start = 40
        e, v = eigh(op / 8)
        t = np.array([0, 2, 8, 20, 32])
        p = abs((np.exp(-1j * np.outer(t, e)) * v[start]) @ v.T) ** 2
        counts = np.array([rng.multinomial(SHOTS, x / x.sum()) for x in p])
        same = where // 2 == where[start] // 2
        same_heads = np.array(history) == history[start]
        categories = np.array(
            [
                same & same_heads,
                same & ~same_heads,
                ~same & same_heads,
                ~same & ~same_heads,
            ]
        )
        assert np.all(categories.sum(axis=0) == 1)
        runs.append(
            dict(
                name=name,
                times=t.tolist(),
                initial_port=int(route[where[start]]),
                probability_by_category=(p @ categories.T).tolist(),
                counts_by_category=(counts @ categories.T).tolist(),
                full_word_counts=counts.tolist(),
            )
        )
    return runs


def run():
    assert sha256((HERE / "MOMENTUM.md").read_bytes()).hexdigest() == PROTOCOL
    tables = reservoir.native_contacts()
    rows = [mode_comparison(n, tables) for n in [4, 6, 10]]
    return dict(
        protocol_sha256=PROTOCOL,
        finite_momentum=rows,
        charge_position_and_restored_heads=propagator(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    out = json.dumps(run(), indent=2, sort_keys=True, default=lambda v: v.item()) + "\n"
    target = HERE / "momentum.json"
    if write:
        target.write_text(out)
    elif target.read_text() != out:
        raise AssertionError("Translation-spectrum replay mismatch")
    print("Verified all registered translation blocks and sharp-vacancy readouts.")
