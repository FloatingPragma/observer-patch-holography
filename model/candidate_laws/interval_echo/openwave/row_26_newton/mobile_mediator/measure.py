"""Connected thermal and Ramsey responses of one mobile between two clamped dipoles."""

import argparse
from collections import deque
from fractions import Fraction
from itertools import product
import json
from math import factorial, log, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.material_inertia.pair_recoil.force import cell_paths


HERE = Path(__file__).parent
SIGNS = np.array([1, -1, -1, 1])
SETTINGS = ((0, 0), (1, 0), (0, 1), (1, 1))
TIMES = (4, 8)
SHOTS = 2**28


def contact(k):
    """Twice the exact connected diagonal of the native nine-load seam."""
    law = candidate.native_law
    h = law.seam_term(k)
    exchange, contrast = law.local_term_parts()
    np.testing.assert_allclose((exchange + k*contrast)[:, :9], np.vstack(
        (h, np.zeros((90, 9)))), atol=4e-14, rtol=0)
    diagonal = np.array([[h[law.SRC.index((i, j)), law.SRC.index((i, j))]
                          for j in range(3)] for i in range(3)])
    twice = 2*(diagonal - diagonal[:, 1:2] - diagonal[1:2, :] + diagonal[1, 1])
    assert np.array_equal(twice, twice.astype(int))
    return twice.astype(int)


def sources(board):
    seams, paths = cell_paths(board)
    for distance, path in sorted(paths.items()):
        a = seams[path[0], path[1]][0]
        b = seams[path[-2], path[-1]][1]
        def companion(p):
            return min(j for j in board.neighbours[p] if board.degree[j] == 5)
        yield distance, path, (a, companion(a), b, companion(b))


def operators(board, pins, k, q, reverse):
    adjacency = board.adjacency().toarray().astype(int)
    allowed = np.setdiff1d(np.arange(board.n), pins)
    a, a2, b, b2 = pins
    d = board.degree[allowed]
    graph = adjacency[np.ix_(allowed, allowed)]
    # Changing one occupied site from 1 to 1+q, including its counting
    # term on the original board and every diagonal touch of a clamped site.
    diagonal = (2*k + 1)*d + 2*k*(6-d)*((1+q)**2-1)
    h2 = np.diag(diagonal) - graph
    v = contact(k)[1+q]
    va = v[0]*adjacency[allowed, a] + v[2]*adjacency[allowed, a2]
    vb = v[2*reverse]*adjacency[allowed, b] + v[2-2*reverse]*adjacency[allowed, b2]
    return graph, h2, va, vb


def first_mixed_moment(graph, va, vb):
    """The first possible two-insertion trace coefficient, in exact integers."""
    left, right = np.flatnonzero(va), np.flatnonzero(vb)
    neighbours = [np.flatnonzero(row) for row in graph]
    near = len(graph)
    numerator = 0
    for a in left:
        distances = {a: 0}
        counts = {a: 1}
        queue = deque([a])
        while queue:
            vertex = queue.popleft()
            for neighbour in neighbours[vertex]:
                step = distances[vertex] + 1
                if neighbour not in distances:
                    distances[neighbour] = step
                    counts[neighbour] = 0
                    queue.append(neighbour)
                if distances[neighbour] == step:
                    counts[neighbour] += counts[vertex]
        for b in right:
            d = distances[b]
            if d < near:
                near, numerator = d, 0
            if d == near:
                numerator += int(va[a])*int(vb[b])*counts[b]**2
    order = 2*near + 2
    moment = Fraction(order*numerator, 2**order)
    return dict(neighbourhood_distance=near, order=order,
                signed_path_weight=numerator, trace_moment=str(moment))


def thermal_trace(h2, va, vb):
    """Enclose the mixed heat trace at β=1/8 without subtracting four traces."""
    beta = 1/8
    shift = 128
    n = len(va)
    # M_ij = β(shift I − H_ij). The four sectors below propagate the
    # baseline, each single change, and the mixed change on the same bank.
    mats = [sparse.csr_matrix((2*shift*np.eye(n) - h2 - np.diag(i*va+j*vb))*beta/2)
            for i, j in SETTINGS]
    a, b = sparse.diags(-beta*va/2), sparse.diags(-beta*vb/2)
    generator = sparse.bmat([[mats[0], None, None, None],
                             [a, mats[1], None, None],
                             [b, None, mats[2], None],
                             [None, b, a, mats[3]]], format="csr")
    majorant = abs(generator)
    limit = 180
    maxsum = int(np.ceil(majorant.sum(axis=1).max()))
    assert maxsum <= 24 and max(np.diff(generator.indptr)) <= 9

    def series(matrix):
        term = np.zeros((4*n, n))
        term[:n] = np.eye(n)
        total = 0.
        for k in range(1, limit+1):
            term = (matrix @ term) / k
            total += np.trace(term[3*n:])
        return float(total)

    center = series(generator)
    absolute_sum = series(majorant)
    # Each product has <=9 dyadic entries per row. Through 180 powers,
    # division by k, traces of length n and summation, gamma < 2^-40.
    # The same estimate applies to the computed positive majorant.
    assert (20*limit+n+limit)*2.0**-53 < 0.99*2.0**-40
    # Tail bounded by first omitted norm power and the geometric ratio.
    tail = Fraction(n*maxsum**(limit+1), factorial(limit+1)) / (
        1 - Fraction(maxsum,limit+2))
    error = 2**-39*absolute_sum + 2*float(tail) + 1e-240
    # exp(16) is between these two decimals, verified below using rational
    # Taylor terms and a geometric bound. Dividing an interval of fixed sign
    # by both endpoints gives an outward enclosure for the unshifted trace.
    e0, e1 = 8886110, 8886111
    s = sum(Fraction(16**k, factorial(k)) for k in range(81))
    assert e0 < s < s + Fraction(16**81,factorial(81))/(1-Fraction(16,82)) < e1
    bounds = [(center + sign*error)/e for sign in (-1,1) for e in (e0,e1)]
    return [float(np.nextafter(min(bounds), -np.inf)), float(np.nextafter(max(bounds), np.inf))]


def spectra(h):
    """Full spectrum and a residual enclosure for its lowest energy."""
    eigen, vectors = eigh(h, driver="evd")
    residual = np.linalg.norm(h @ vectors - vectors*eigen, "fro")
    gram_error = np.linalg.norm(vectors.T @ vectors - np.eye(len(h)), "fro")
    # Classical n-term dot product bounds. Use a large dyadic margin for
    # both computed Frobenius norms, their products and the norm reductions.
    n = len(h)
    margin = 16 * n * np.finfo(float).eps * sqrt(n)
    magnitude = max(abs(h).sum(axis=1).max(), max(abs(eigen)), 1.)
    assert gram_error + margin*n < .01
    error = 2*(residual + margin*magnitude)
    return eigen, error


def run():
    rng = np.random.default_rng(260027)
    board = federation(0)
    out = []
    # 40 panels, four settings, two quadratures and two appointments.
    banks = 40*4*2*2
    eps = sqrt(log(2*banks/.01)/(2*SHOTS))
    for k, (distance,path,pins), q, reverse in product(
            (candidate.KAPPA, 0), tuple(sources(board)), (1,-1), (False,True)):
        graph, h2, va, vb = operators(board, pins, k, q, reverse)
        assert not np.any(va * vb)
        moment = first_mixed_moment(graph, va, vb)
        readings = [spectra((h2 + np.diag(i*va+j*vb))/2) for i, j in SETTINGS]
        eigen = np.array([x for x, _ in readings])
        lower = [float(x[0]) for x in eigen]
        connected = np.dot(SIGNS, lower)
        err = sum(e for _, e in readings)
        frozen = [float((h2.diagonal()+i*va+j*vb).min()/2) for i,j in SETTINGS]
        meters = []
        for t in TIMES:
            r = np.exp(-1j*t*eigen).mean(axis=1)
            for name, amplitude in (("real",r.real), ("imag",r.imag)):
                clicks = rng.binomial(SHOTS, (1+amplitude)/2)
                meters.append(dict(wait=t, component=name, clicks=clicks.tolist(),
                                   predicted=float(SIGNS @ amplitude),
                                   measured=float(2*int(SIGNS @ clicks)/SHOTS)))
        heat = thermal_trace(h2,va,vb)
        finite_potential = None
        if k == candidate.KAPPA:
            # Positive mobile energies in all four settings bound Z by n.
            row_hops = graph.sum(axis=1)
            assert all(np.min(h2.diagonal()+i*va+j*vb-row_hops) >= 0 for i,j in SETTINGS)
            remainder = 16*len(va)**2*1e-55
            finite_potential = [float(np.nextafter(-8*heat[1]-remainder, -np.inf)),
                                float(np.nextafter(-8*heat[0]+remainder, np.inf))]
        out.append(dict(kappa=k, distance=distance, path=path, pins=pins, mobile_charge=q,
                        reversed_receiver=reverse, first_mixed=moment,
                        ground_energies=lower, connected_ground=float(connected),
                        ground_error_radius=float(err),
                        frozen_connected_ground=float(np.dot(SIGNS,frozen)),
                        heat_trace_beta_eighth=heat, ramsey=meters,
                        potential_per_activity_beta_eighth_z1e_minus55=finite_potential))
    return dict(source="20-cell 240-port blank-memory Interval Echo; 236 mobile positions",
                setting_order=["00","10","01","11"], shots_per_bank=SHOTS,
                joint_99pct_connected_ramsey_radius=8*eps, panels=out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE/"results.json").write_text(result)
    else:
        print(result, end="")
