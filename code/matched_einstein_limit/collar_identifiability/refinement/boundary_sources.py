"""Sign-resolved collar response for a patch of registered, shared native seam outcomes."""

import json
from math import atanh, comb

import numpy as np
from scipy.stats import binom

from model.executable_model import federation
from . import admission, transport


def boundary_coefficients(board, enrolled, placement):
    """Coefficients of each independent signed gateway word at the source."""
    owners = np.array(board.external, dtype=int) // 12
    weight = np.zeros(board.cells)
    weight[enrolled] = placement
    return weight[owners[:, 0]] - weight[owners[:, 1]]


def separated(board, count):
    neighbors = [set() for _ in range(board.cells)]
    for p, q in board.external:
        a, b = p//12, q//12
        neighbors[a].add(b)
        neighbors[b].add(a)
    enrolled, excluded = [], set()
    for cell in range(board.cells):
        if cell not in excluded:
            enrolled.append(cell)
            excluded.update(neighbors[cell])
            if len(enrolled) == count:
                return enrolled
    raise ValueError("Insufficient separated cells on this board")


def conditional_charge(edges, gain, reference):
    """Exact binomial outcome law. The returned value is the weakening slope per unit pairing."""
    k = np.arange(edges + 1)
    y = gain*(2*k-edges)
    if np.max(abs(y)) >= 1 or not 0 < reference < 1:
        raise ValueError("the two-sided pairing variation requires faithful binary states")
    return float(binom.pmf(k, edges, .5) @ (y * np.arctanh(reference*y)))


def absolute_flux(edges):
    """Mean absolute sum of independent fair signs, from the exact binomial identity."""
    if edges < 0:
        raise ValueError("the number of boundary edges must be nonnegative")
    if edges == 0:
        return 0.
    return edges * comb(edges - 1, (edges - 1)//2) / 2**(edges - 1)


def run():
    transport.native_fair_bit()
    reading = admission.material_panel(1, 1, 8)["records"][-1]
    z = reading["d"][-1]
    reference = .5
    rows = []
    for level in (0, 1, 2):
        board = federation(level)
        ports = np.asarray(board.external)
        assert np.unique(ports).size == ports.size
        n = 4**level
        for name, sites in (("triangular ancestor", list(range(n))),
                            ("isolated cells", separated(board, n))):
            weights = boundary_coefficients(board, sites, np.ones(n))
            count = int(weights @ weights)
            assert count == (3*2**level if name == "triangular ancestor" else 3*n)
            u = 1/4096
            hit = admission.any_hit(u, n)
            gain = z*hit/(3*n)
            slope = conditional_charge(count, gain, reference)
            expected = reference*z*z*count/9
            mean_flux = absolute_flux(count)
            aligned = gain*mean_flux
            variance = gain**2*count
            rows.append(dict(cells=board.cells, enrolled=n, source_shape=name, boundary=count,
                             source_variance=count, coefficient=expected, observed_per_u_squared=slope/u**2,
                             mean_absolute_flux=mean_flux, one_bit_coefficient=reference*z*z*mean_flux**2/9,
                             one_bit_observed_per_u_squared=aligned*atanh(reference*aligned)/u**2,
                             two_unaligned_returns=2*reference*variance*atanh(reference**2*variance),
                             all_flags_hidden=0))
    return dict(material_return_contrast=z, reference_pairing=reference, rows=rows)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
