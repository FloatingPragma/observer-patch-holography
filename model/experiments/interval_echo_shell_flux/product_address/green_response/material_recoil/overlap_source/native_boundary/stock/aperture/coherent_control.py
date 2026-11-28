"""Local discarded readings with matched ordinary heat on the occupied chart."""

import numpy as np

from . import regional_tomography as meter

AGES = (0, 1, 4, 16, 64)


def ports(live, survey, volume, clamp, edges):
    masks = [mask for (n, _, _, _), (mask, _) in survey.items() if n == volume]
    interior = np.logical_and.reduce(masks)
    exterior = ~np.logical_or.reduce(masks)
    reserved = set([clamp, *np.ravel(edges)])

    def first(side):
        sites = [i for i in np.flatnonzero(side) if i not in reserved]
        deep = [
            i
            for i in sites
            if side[live.indices[live.indptr[i] : live.indptr[i + 1]]].all()
        ]
        return int((deep or sites)[0]), bool(deep)

    x, deep_x = first(interior)
    y, deep_y = first(exterior)
    return (x, y), [deep_x, deep_y]


def work(reconstructed, errors, h, x):
    """Energy inserted when the defect-presence reading at x is discarded."""
    row = -2 * h[x].copy()
    row[x] = 0
    return np.array([row @ reconstructed[x].real, abs(row) @ errors[x]])


def propagation(e, v, rho, x):
    phase = np.exp(-1j * e[:, None] * AGES)
    off = rho[:, x].copy()
    off[x] = 0
    return v @ (v[x, :, None] * phase), v @ ((v.T @ off)[:, None] * phase)


def readouts(mask, cut, h, rho, x, propagated, fringe_error, rec):
    inside = np.flatnonzero(mask)
    rb = rho[np.ix_(inside, inside)]
    vals, basis = np.linalg.eigh(rb)
    log = (basis * np.log(vals)) @ basis.T
    log.flat[:: len(log) + 1] -= np.log(1 - np.trace(rb))
    u, v = propagated
    # The pinching tangent at age t is -u v* - v u*.
    entropy = 2 * np.sum(v[inside].conj() * (log @ u[inside]), axis=0).real
    i, j = cut
    fringe = -2 * np.mean((u[i] * v[j].conj() + v[i] * u[j].conj()).real, axis=0)
    baseline = 2 * np.mean(rho[i, j])
    # Jeffreys divergence after a balanced interior port read. Use a single
    # measured real quadrature, with its lower contrast and upper mean hit rate.
    bound = 0.0
    if mask[x]:
        candidates = np.flatnonzero(mask & (h[x] < 0))
        for z in candidates:
            s = max(0, abs(rec[x, z].real) - fringe_error[x, z])
            p = (
                rec[x, x].real
                + rec[z, z].real
                + fringe_error[x, x]
                + fringe_error[z, z]
            ) / 2
            if s:
                bound = max(bound, s * np.log1p(2 * s / (p - s)))
    assert entropy[0] >= bound - 2e-12
    return entropy, fringe, baseline, bound


def measured_difference(baseline, first, second, rng):
    # Select a uniformly random crossing link; split idle trials equally
    # between the two recorded signs. Mean sign is twice its real fringe.
    probs = (1 + baseline + np.stack((first, second))) / 2
    counts = rng.binomial(meter.N, probs)
    middle = 2 * (counts[0] - counts[1]) / meter.N
    width = 4 * meter.HOEFFDING
    return counts.tolist(), np.stack((middle - width, middle + width), axis=1)
