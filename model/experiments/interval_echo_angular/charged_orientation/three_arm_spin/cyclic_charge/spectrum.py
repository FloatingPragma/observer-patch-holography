"""Exact character symmetries and source energy moments of the native chart."""

from collections import Counter
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix

from ..apparatus import permutations_of_chart
from ..measure import CY, SWAP, ONE, reduce


def symmetry(h):
    perms = permutations_of_chart()
    c, s = perms[CY], perms[SWAP]
    assert np.array_equal(c[c[c]], np.arange(len(c)))
    assert np.array_equal(s[c[s]], c[c])
    for perm in (c, s):
        assert (h[perm, :][:, perm] != h).nnz == 0
    fixed = int(sum(c == np.arange(len(c))))
    rank = (len(c) - fixed) // 3
    return (
        dict(
            native_states=len(c),
            fixed_words=fixed,
            cyclic_ranks=[rank + fixed, rank, rank],
        ),
        perms,
    )


def exact_moments(h):
    twice = (2 * h).tocsc()
    np.testing.assert_array_equal(twice.data, np.rint(twice.data))
    # By arm permutation symmetry the three sharp-source diagonal moments
    # are equal and all six off-diagonals coincide exactly.
    x = {ONE[0]: 1}
    for step in range(32):
        cross = x.get(ONE[1], 0)
        assert cross == x.get(ONE[2], 0)
        if cross:
            diagonal = x.get(ONE[0], 0)
            return dict(
                power_of_2H=step,
                sharp_diagonal=diagonal,
                sharp_cross=cross,
                symmetric=diagonal + 2 * cross,
                each_winding=diagonal - cross,
                difference=3 * cross,
            )
        y = Counter()
        for j, amp in x.items():
            for k in range(twice.indptr[j], twice.indptr[j + 1]):
                y[int(twice.indices[k])] += amp * int(twice.data[k])
        x = y
    raise AssertionError("No distinction in explored moment range")


def symmetric_basis(h, perms):
    row = []
    col = []
    values = []
    seen = set()
    n = 0
    for i in range(h.shape[0]):
        if i in seen:
            continue
        orbit = {int(p[i]) for p in perms.values()}
        seen.update(orbit)
        for j in orbit:
            row.append(j)
            col.append(n)
            values.append(1 / np.sqrt(len(orbit)))
        n += 1
    return coo_matrix((values, (row, col)), shape=(h.shape[0], n)).tocsr()


def source_spectra(h, perms):
    one = symmetric_basis(h, perms)
    standard, small = reduce(h, perms)
    source = np.zeros(h.shape[0])
    source[ONE] = 1 / np.sqrt(3)
    odd = np.zeros_like(source)
    odd[ONE] = [0, 1 / np.sqrt(2), -1 / np.sqrt(2)]
    output = {}
    for name, mat, vec in [
        ("symmetric", (one.T @ h @ one).toarray(), one.T @ source),
        ("each_winding", small, standard.T @ odd),
    ]:
        np.testing.assert_allclose(np.linalg.norm(vec), 1, atol=2e-12)
        ev, V = eigh(mat, driver="evd")
        weights = (vec @ V) ** 2
        breaks = np.r_[0, np.flatnonzero(np.diff(ev) > 1e-8) + 1, len(ev)]
        lines = [(ev[i], float(sum(weights[i:j]))) for i, j in zip(breaks, breaks[1:])]
        visible = [(e, w) for e, w in lines if w > 1e-8]
        output[name] = dict(
            representation_rows=len(ev),
            threshold_weight=1e-8,
            number_above_threshold=len(visible),
            lowest_visible=visible[0],
            strongest=max(visible, key=lambda z: z[1]),
            spectral_mean=float(ev @ weights),
            spectral_variance=float((ev**2) @ weights - (ev @ weights) ** 2),
        )
    return output
