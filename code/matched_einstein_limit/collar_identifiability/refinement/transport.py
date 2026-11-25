"""Visit density from a single compact injection on a finite three-dimensional comparison chart."""

import json
from math import pi

import numpy as np
from scipy.fft import dstn, idstn
from scipy.ndimage import convolve

from model.candidate_laws.interval_echo import candidate


def native_fair_bit():
    from scipy.linalg import expm

    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    indices = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    start = np.eye(len(h))[:, indices[0]]
    probability = abs(expm(-1j * (pi/2) * h) @ start)**2
    np.testing.assert_allclose(probability[indices], .5, rtol=0, atol=2e-14)
    np.testing.assert_allclose(probability[indices].sum(), 1, rtol=0, atol=2e-14)


def chart(radius, survival):
    # Outside the cube a mark is absorbed. The clock begins at injection.
    width = 2*radius-1
    c = radius-1
    inject = np.zeros((width,) * 3)
    inject[c, c, c] = 1.
    cos = np.cos(pi*np.arange(1, width+1) / (width+1))
    eigen = 1 - survival * (1 + cos[:, None, None] + cos[None, :, None] + cos[None, None, :])/4
    visits = idstn(dstn(inject, type=1, norm="ortho") / eigen, type=1, norm="ortho")
    stencil = np.zeros((3, 3, 3))
    stencil[1, 1, 1] = .25
    for a in range(3):
        for step in (-1, 1):
            index = [1, 1, 1]
            index[a] += step
            stencil[tuple(index)] = .125
    residual = visits - survival*convolve(visits, stencil, mode="constant", cval=0.) - inject
    if abs(residual).max() > 1e-12:
        raise ArithmeticError("the reconstructed visits fail the local injection equation")
    escape = sum(np.take(visits, j, axis=axis).sum() for j in (0, -1) for axis in (0, 1, 2))
    sink = survival*escape/8 + (1-survival)*visits.sum()
    np.testing.assert_allclose(sink, 1., rtol=1e-11)
    rows = []
    for r in (2, 4, 8, 16):
        first, second = [visits[c + s, c, c] for s in (r, 2*r)]
        rows.append(dict(r=r, visits=float(first), doubled_visits=float(second),
                         local_density_ratio=float((first/second)**2)))
    return dict(radius=radius, survival=survival, mean_opportunities=float(visits.sum()),
                total_sink=float(sink), maximum_cell_residual=float(abs(residual).max()), rows=rows)


def axial_visits(radius, offsets, survival=1.):
    """Exact finite-box Green function on the central axis, with two sine transforms summed."""
    if not 0 < survival <= 1 or any(s < 0 or s >= radius for s in offsets):
        raise ValueError("need a positive survival and offsets inside the cube")
    odd = pi*np.arange(1, 2*radius, 2) / (2*radius)
    cosine = np.cos(odd)
    angle = np.arccosh(4/survival - 1 - cosine[:, None] - cosine[None, :])
    denominator = ((survival/4) * np.sinh(angle) * (1 + np.exp(-2*radius*angle)))
    return [float(np.sum(np.exp(-r*angle) * -np.expm1(-2*(radius-r)*angle) / denominator)/radius**2)
            for r in offsets]


def run():
    native_fair_bit()
    boxes = [chart(64, survival) for survival in (1., 7/8)]
    for box in boxes:
        axial = axial_visits(64, (2, 4, 8, 16), box['survival'])
        np.testing.assert_allclose(axial, [row['visits'] for row in box['rows']], atol=2e-14)
    larger = []
    for radius in (128, 256, 512, 1024):
        r = int(np.sqrt(radius) / 2)
        near, far = axial_visits(radius, (r, 2*r))
        larger.append(dict(radius=radius, r=r, density_ratio=(near/far)**2))
    return dict(complete_box=boxes, larger_boxes=larger)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
