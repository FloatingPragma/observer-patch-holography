"""A posteriori enclosure for the measured positive pole on the exact symmetry row."""

import numpy as np


def isolate(h, groups, q, vec, energies, index, *, dyadic_bits=43):
    # Quantize to a dyadic basis *in the exact standard representation*. All
    # additions and integer permutations in 6 times the projector are exact.
    scale = 2**dyadic_bits
    rotate = groups[1, 2, 0]
    mirror = groups[0, 2, 1]
    source = np.asarray(q @ vec)
    z = np.rint(source * (scale / 6)).astype(np.int64)
    z = 2 * z - z[rotate] - z[rotate[rotate]]
    z = z - z[mirror]
    w = z.astype(float) / scale
    n, m = w.shape
    indices = np.arange(n)

    # trace[(I-S)(2I-C-C²)] / 6 is the exact dimension of this row.
    def fixed(perm):
        return int(np.count_nonzero(perm == indices))

    numerator = 2 * n - fixed(rotate) - fixed(rotate[rotate]) - 2 * fixed(mirror)
    numerator += fixed(rotate[mirror]) + fixed(rotate[rotate[mirror]])
    assert numerator == 6 * m
    eps = 2.0**-53
    # Row sums of 2H are integers. Dot products against z are exact within
    # 53 bits after allowing for each sparse row's full absolute coefficient sum.
    row_max = float(np.max(abs(h).sum(axis=1)))
    assert np.max(abs(z)) * 2 * row_max < 2**53
    gram = w.T @ w - np.eye(m)
    # Gram dot error in Frobenius norm ≤ γ_n ||W||_F², with
    # ||W||_F² < 2m. Leave a factor of two also on computed norm evaluations.
    assert np.max(abs(np.diag(gram))) < 1e-7
    delta = 2 * np.linalg.norm(gram, "fro") + 4 * n * eps * m
    residual = h @ w - w * energies
    frequency = float(max(abs(energies)))
    rho = 2 * np.linalg.norm(residual, "fro") + 8 * eps * (
        row_max + frequency
    ) * np.sqrt(m)
    # Polar-correct W to an isometry and apply the self-adjoint eigenvalue
    # perturbation bound to the resulting approximate diagonalization.
    radius = rho / np.sqrt(1 - delta) + 4 * frequency * delta / (1 - delta)
    lo = energies[index] - energies[index - 1]
    hi = energies[index + 1] - energies[index]
    gap = min(lo, hi) - 2 * radius
    assert gap > 0
    v = w[:, index]
    specific = 2 * np.linalg.norm(residual[:, index]) + 8 * eps * (row_max + frequency)
    norm = np.linalg.norm(v)
    error = (
        2 * specific / (norm * gap)
        + abs(norm - 1)
        + 2 * np.linalg.norm(v - source[:, index])
        + 4 * n * eps
    )
    # Each projected source column differs by at most 2*eigenvector error;
    # the second symmetry row comes from (C-C²)/√3. A factor four covers both.
    return dict(
        standard_row_dimension=m,
        gram_error=float(delta),
        residual=float(rho),
        all_mode_radius=float(radius),
        interval=[
            float(energies[index] - specific / norm),
            float(energies[index] + specific / norm),
        ],
        nearest_other_standard_distance=float(gap),
        eigenvector_error=float(error),
        bipartite_amplitude_error=float(4 * error),
    )
