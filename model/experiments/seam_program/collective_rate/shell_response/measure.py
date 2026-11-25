"""Four-port parity inside the eight-port thermal region, at fixed pair and shell records."""

import json
from fractions import Fraction
from itertools import combinations
from math import comb, factorial

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import cell, rotations
from ..measure import clock_seam, reduced_gibbs_states, spectrum


# Addresses chosen from the four-port parity scan at beta=1/16.
PARITY_PORTS = (0, 1, 7, 10)


def region():
    board = cell()
    a, b = clock_seam()
    shell = sorted((set(board.neighbours[a]) | set(board.neighbours[b])) - {a, b})
    ports = [a, b, *shell]
    return ports


def parities(ports):
    orbit = sorted({tuple(sorted(r[i] for i in PARITY_PORTS)) for r in rotations()})
    probes = {"single": [PARITY_PORTS], "rotated": orbit, "stationary": list(combinations(range(12), 4))}
    numerators = []
    for name, supports in probes.items():
        visible = [s for s in supports if set(s) <= set(ports)]
        if name == "rotated" and ports == region():
            assert len(supports) == 30 and len(visible) == 9
            assert all(set(s) & set(ports[:2]) for s in visible)
        masks = [sum(1 << ports.index(i) for i in support) for support in visible]
        numerators.append([sum((-1)**((word & m).bit_count()) for m in masks)
                           for word in range(2**len(ports))])
    return list(probes), np.array(numerators), [len(s) for s in probes.values()]


def moments(ports, order, exchange):
    """Unnormalized partial traces of the integer powers of -2(H-mu N)."""
    board = cell()
    hidden = 2**(board.n - len(ports))
    dim = 2**len(ports)
    total = [np.zeros((dim, dim), dtype=np.int64) for _ in range(order + 1)]
    for number in range(board.n + 1):
        sector = BinarySector(board, number)
        h = sector.hamiltonian(KAPPA, shifted=False)
        np.fill_diagonal(h, h.diagonal() - KAPPA * number)
        if not exchange:
            h -= sector.K
        integral = (-2*h).astype(np.int64)
        assert np.array_equal(integral, -2*h)
        # Bounds every intermediate dot product and the sum over exterior words.
        row_bound = int(np.max(abs(integral).sum(axis=1)))
        assert hidden * row_bound**order < 2**63
        generator = csr_matrix(integral)
        power = np.eye(len(integral), dtype=np.int64)
        masks = np.array([sum(1 << p for p in state) for state in sector.states])
        inside = sum(((masks >> p) & 1) << j for j, p in enumerate(ports))
        outside = masks & ~sum(1 << p for p in ports)
        groups = [np.flatnonzero(outside == x) for x in np.unique(outside)]
        for k in range(1, order + 1):
            power = generator @ power
            for selected in groups:
                total[k][np.ix_(inside[selected], inside[selected])] += power[np.ix_(selected, selected)]
    return total


def log_series(order=5, exchange=True):
    """Exact Taylor coefficients of Tr[(Z/256)(-log rho_region)] in inverse temperature."""
    ports = region()
    names, parity, denominators = parities(ports)
    dim = 2**len(ports)
    hidden = 2**(cell().n - len(ports))
    raw = moments(ports, order, exchange)
    coefficients = {name: [Fraction(0) for _ in range(order + 1)] for name in names}
    for occupation in range(len(ports) + 1):
        ids = [i for i in range(dim) if i.bit_count() == occupation]
        # F = E_exterior[exp(-beta H')]. Its j-th derivative in beta/(2*hidden)
        # is hidden^(j-1) times the unnormalized partial trace of (-2H')^j.
        one = {j: raw[j][np.ix_(ids, ids)].astype(object) * hidden**(j-1)
               for j in range(1, order + 1)}
        powers = one
        for factors in range(1, order + 1):
            for j, derivative in powers.items():
                traces = parity[:, ids] @ derivative.diagonal()
                for name, norm, trace in zip(names, denominators, traces, strict=True):
                    coefficients[name][j] += Fraction((-1)**factors * int(trace),
                        norm * dim * factors * factorial(j) * (2*hidden)**j)
            powers = {n: sum(comb(n, j) * (one[j] @ powers[n-j])
                             for j in range(1, n - factors + 1))
                      for n in range(factors + 1, order + 1)}
    return coefficients


def finite_measurements():
    ports = region()
    names, parity, denominators = parities(ports)
    blocks = spectrum(KAPPA)
    rows = []
    for beta in (1/32, 1/16, 1/8):
        state, = reduced_gibbs_states(blocks, beta, [ports])
        weights, vectors = eigh(state)
        diagonal = (vectors**2) @ -np.log(weights)
        slopes = parity @ diagonal / len(diagonal) / denominators
        rows.append(dict(beta=beta, modular_slopes=dict(zip(names, slopes.tolist(), strict=True)),
                         minimum_region_probability=float(weights.min())))
    return rows


if __name__ == "__main__":
    print(json.dumps(dict(ports=region(), parity_ports=PARITY_PORTS,
                          series={k: list(map(str, v)) for k,v in log_series().items()},
                          finite=finite_measurements()), indent=2))
