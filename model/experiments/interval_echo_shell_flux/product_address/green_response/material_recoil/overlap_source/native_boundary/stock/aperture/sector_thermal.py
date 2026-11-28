"""Ordinary energy spectra and modular observables on the measured four-home support."""

from functools import lru_cache

import numpy as np
from scipy.special import logsumexp

from . import closed_pairing as pair
from . import thermal_kernel

BETA = thermal_kernel.BETA


@lru_cache(maxsize=1)
def apparatus():
    _, words, cells, _, _, _, groups, count = pair.setup()
    edges = pair.stock.carrier()[3]
    _, _, native, base = pair.stock.kernels()
    index = {tuple(w): i for i, w in enumerate(words)}
    operations = np.zeros((2, 3, len(words), len(words)), dtype=np.int16)
    for j, (left, right) in enumerate(edges):
        for i, w in enumerate(words):
            col = 9 * w[4 + j] + 3 * w[left] + w[right]
            for layer, h in zip(operations, (native, native - base), strict=True):
                for k in np.flatnonzero(h[:, col]):
                    to = w.copy()
                    to[4 + j], ports = divmod(k, 9)
                    to[left], to[right] = divmod(ports, 3)
                    end = index[tuple(to)]
                    assert groups[end] == groups[i]
                    layer[j, end, i] = h[k, col]
    assert np.array_equal(operations, operations.transpose(0, 1, 3, 2))
    blocks = [np.flatnonzero(groups == j) for j in range(count)]
    return words, cells, groups, blocks, operations


def bath(beta=BETA, stiffness=1.0):
    w, _, _, blocks, ops = apparatus()
    H, W = ops.sum(1)
    H1 = (H + (stiffness - 1) * W) / 4
    raw = np.zeros_like(H1)
    kept = np.zeros_like(H1)
    for group in blocks:
        ix = np.ix_(group, group)
        eig, basis = np.linalg.eigh(H1[ix])
        log_z = logsumexp(-beta * eig)
        rho = (basis * np.exp(-beta * eig - log_z)) @ basis.T
        kept[ix] = rho * (len(group) / len(w))
        raw[ix] = rho * np.exp(log_z)
    return {"retained": kept, "exchanged": raw / np.trace(raw)}


@lru_cache(maxsize=None)
def region_indices(start, stop):
    w = apparatus()[0]
    columns = [*range(start, stop + 1), *range(4 + start, 4 + stop)]
    outside = [i for i in range(7) if i not in columns]
    _, inside = np.unique(w[:, columns], axis=0, return_inverse=True)
    if outside:
        _, complement = np.unique(w[:, outside], axis=0, return_inverse=True)
    else:
        complement = np.zeros(len(w), int)
    ii, jj = np.nonzero(complement[:, None] == complement[None, :])
    size = inside.max() + 1
    return inside, ii, jj, inside[ii] * size + inside[jj], size


def restrict(rho, region):
    _, i, j, target, n = region
    return np.bincount(target, weights=rho[i, j], minlength=n * n).reshape(n, n)


def operators(reference, region):
    loc = restrict(reference, region)
    eig, basis = np.linalg.eigh(loc)
    assert min(eig) > 0
    K = (basis * (-np.log(eig))) @ basis.T
    sharp = -np.log(np.diag(loc))
    return K, sharp


def panel(source, ordinary, therm, shifted, start, stop):
    _, _, _, _, ops = apparatus()
    region = region_indices(start, stop)
    inside = region[0]
    # Calibrate a common ruler for this region. An ordinary heating step
    # and a stiffness step give independent source directions where possible.
    K, sharp = operators(ordinary, region)
    joint_scores = K.diagonal()[inside]
    sharp_scores = sharp[inside]
    eta = np.array(
        [
            0.5 * (int(start <= i <= stop) + int(start <= i + 1 <= stop))
            for i in range(3)
        ]
    )
    E, W = np.einsum("j,ljab->lab", eta, ops) / 4
    measurements = []
    for d in (therm, shifted):
        # Both controls are read at the original instrument settings.
        measurements.append(
            [np.sum(d * E), np.sum(d * W), np.sum(restrict(d, region) * K)]
        )
    A = np.array(measurements)
    y = [source @ E.diagonal(), source @ W.diagonal(), source @ joint_scores]
    singular = np.linalg.svd(A[:, :2], compute_uv=False)
    fit = None
    if min(singular) > 1e-9 * max(singular):
        coeff = np.linalg.solve(A[:, :2], A[:, 2])
        predicted = float(coeff @ y[:2])
        fit = dict(
            energy_and_work_coefficients=coeff.tolist(),
            prediction=predicted,
            held_residual=float(y[2] - predicted),
            condition_number=float(max(singular) / min(singular)),
        )
    return (
        dict(
            dimension=int(region[-1]),
            quantum_modular_nats=float(y[2]),
            sharp_modular_nats=float(source @ sharp_scores),
            energy=float(y[0]),
            fractional_stiffness_work=float(y[1]),
            ordinary_control_rows=A.tolist(),
            calibrated=fit,
        ),
        joint_scores,
        sharp_scores,
    )
