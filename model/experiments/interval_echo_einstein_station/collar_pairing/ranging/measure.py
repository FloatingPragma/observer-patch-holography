"""Censored two-way clock ranges and propagation fronts on the freely formed collar."""

import json
from itertools import combinations, permutations
from math import factorial

import numpy as np
from scipy import sparse
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import expm_multiply

from ...mobile_collar_source import measure as formation
from ...shared_geometry_probe.flight import prepared_states
from ..measure import charge_slope


def parity_blocks(q, k, parity):
    """Check all three load words and both coins within the 99-state contact."""
    emission = formation.emission
    basis = np.zeros((emission.LAW.DIM, 6), int)
    for i, word in enumerate(((1, 1), (1+q, 1), (1, 1+q))):
        for c in range(2):
            for d, s in ((1, 1), (-1, parity)):
                address = emission.columns(word)[1 + emission.LAW.LABELS.index((c, d))]
                basis[address, 2*i+c] = s
    native = emission.four_h(k)
    projected = basis.T @ native @ basis
    np.testing.assert_array_equal(2 * native @ basis, basis @ projected)
    return (projected / 8).reshape(3, 2, 3, 2).transpose(0, 2, 1, 3)


def generator(board, route, q, k, left_odd):
    """Change only the prepared reversal sector of the left head."""
    native = formation.hamiltonian(board, route, q, k)
    if not left_odd:
        return native
    change = parity_blocks(q, k, -1) - parity_blocks(q, k, 1)
    np.testing.assert_array_equal(change[0, 0], np.eye(2))
    ends = sorted(route[:2])
    for i in range(2):
        for j in range(2):
            contact = change[i+1, j+1] - int(i == j)*change[0, 0]
            dyad = sparse.csr_matrix(([1], ([ends[i]], [ends[j]])), shape=(board.n, board.n))
            native += sparse.kron(dyad, np.kron(contact, np.eye(2)), format="csr")
    return native - native.diagonal().mean()*sparse.eye(len(native.diagonal()), format="csr")


def ordinary_write():
    """Verify the ordinary exchange gap and its native write into reversal parity."""
    law = formation.emission.LAW
    native = formation.emission.four_h(8)
    signs = np.array([[1, 1], [1, -1]])
    even, odd = (np.outer(x, x) for x in signs)
    np.testing.assert_array_equal(2*np.diag(even.diagonal()), even+odd)
    for parity in (1, -1):
        source = np.zeros(law.DIM, int)
        source[[law.SRC.index((2, 1)), law.SRC.index((1, 2))]] = (1, parity)
        energy = 8 + (1 - parity)//2
        np.testing.assert_array_equal(native @ source, 4*energy*source)
        written = law._write_pattern() @ source[:law.N_SRC]
        target = np.zeros_like(written).reshape(law.N_LAB, law.N_SRC)
        for c, load in enumerate(((1, 2), (2, 1))):
            for d, s in ((1, 1), (-1, parity)):
                target[law.LABELS.index((c, d)), law.SRC.index(load)] = s
        np.testing.assert_array_equal(written, target.ravel())
        np.testing.assert_array_equal(parity_blocks(1, 8, parity)[0, 0], energy*np.eye(2))


def kernels(h, route, step, steps):
    launch = np.zeros((h.shape[0], 16), complex)
    stations = np.array([4*p+v for p in route for v in range(4)])
    launch[stations] = np.eye(16)
    sampled = expm_multiply(-1j*h, launch, start=0, stop=step*steps,
                            num=steps+1, traceA=0)
    return sampled[:, stations].reshape(steps+1, 4, 4, 4, 4).transpose(0, 1, 3, 2, 4)


def first_hits(u, p, q):
    """First visit to q when q alone is interrogated once per appointment."""
    result = [np.zeros((4, 4), complex)]
    for n in range(1, len(u)):
        prior = sum((u[n-j, q, q] @ result[j] for j in range(1, n)), start=result[0])
        result.append(u[n, q, p] - prior)
    return result


def clocks(step, steps):
    t = step*np.arange(steps+1)
    populations = np.stack((np.sin(t/2)**2, (1+np.sin(t))/2), axis=1)
    # Both pointer settings are supplied on independent copies of the same isolated dimer.
    phase = np.unwrap(np.arctan2(2*populations[:, 1]-1, 1-2*populations[:, 0]))
    return populations, phase


def range_panel(u, step, states):
    steps = len(u)-1
    clock, elapsed = clocks(step, steps)
    rows = []
    for p, q in permutations(range(4), 2):
        outward, back = first_hits(u, p, q), first_hits(u, q, p)
        effects = np.zeros((steps+1, 4, 4), complex)
        for n in range(2, steps+1):
            for j in range(1, n):
                amplitude = back[n-j] @ outward[j]
                effects[n] += amplitude.conj().T @ amplitude
        outbound = np.array([a.conj().T @ a for a in outward])
        arm = {}
        for name, rho in states.items():
            returns = np.einsum("nij,ji->n",effects,rho).real
            reaches = np.einsum("nij,ji->n",outbound,rho).real
            ranges = []
            for horizon in (2, 4, 6):
                end = round(horizon/step)
                if end > steps:
                    continue
                hit = returns[:end+1].sum()
                out = reaches[:end+1].sum()
                if min(hit, out-hit, 1-out) < -1e-10:
                    raise ArithmeticError("inconsistent first-hit branches")
                distance = (elapsed[:end+1] @ returns[:end+1] + elapsed[end]*(1-hit))/2
                ranges.append(dict(horizon=horizon, range=float(distance),
                                   complete=float(hit), waiting_for_return=float(out-hit),
                                   missed_outward=float(1-out)))
            arm[name] = dict(first_return=returns.tolist(), first_relay=reaches.tolist(), horizons=ranges)
        rows.append(dict(p=p, q=q, arms=arm))
    return dict(step=step, clock=clock.tolist(), pairs=rows)


def fronts(h, route, states):
    """Leading one-appointment-per-leg return probability from the shortest paths."""
    rows = []
    for i, j in combinations(range(4), 2):
        d = j-i
        injection = np.zeros((h.shape[0], 4), complex)
        injection[4*route[i]:4*route[i]+4] = np.eye(4)
        for order in range(d):
            z = injection[4*route[j]:4*route[j]+4]
            if np.any(z):
                raise ValueError("a shorter path joins the terminals")
            injection = h @ injection
        block = injection[4*route[j]:4*route[j]+4] / factorial(d)
        turn = block.conj().T @ block
        effect = turn @ turn
        coefficients = {name: float(np.trace(effect @ state).real) for name, state in states.items()}
        rows.append(dict(p=i, q=j, graph_separation=d, tested_return_power=4*d,
                         leading_order_certified_on_faithful_heads=bool(np.any(block)),
                         coefficients=coefficients))
    return rows


def run():
    ordinary_write()
    board = formation.federation(1)
    route = formation.region(board)
    joint, product, _ = prepared_states(board, route, 1)
    joint = joint/np.trace(joint).real
    product = product/np.trace(product).real
    states = dict(formed=joint, independent=product)
    comparison = eigvalsh(joint,product)
    report = dict(cells=board.cells, ports=list(route), forming_load=1,
                  probability_ratio_bounds=comparison[[0,-1]].tolist(),
                  improving_pairing_modular_slope=charge_slope(product, joint-product, 1),
                  ordinary_exchange_gap=1, runs=[])
    for sign, k, odd in ((1,8,False),(1,8,True),(-1,8,False),(-1,8,True),(1,0,False)):
        h = generator(board,route,sign,k,odd)
        u = kernels(h,route, .25,24)
        panels = [range_panel(u[::stride],.25*stride,states) for stride in (1,2)]
        report["runs"].append(dict(probe=sign, contrast=k,left_odd=odd,
                                   fronts=fronts(h,route,states),panels=panels))
    report["clock_comparisons"] = clock_comparisons(report["runs"])
    return report


def clock_comparisons(runs):
    """Ratios that eliminate a free positive clock normalization at each station."""
    rows = []
    for run in runs:
        for panel in run["panels"]:
            for h_index, horizon in enumerate((2,4,6)):
                arms = {}
                for arm in ("formed","independent"):
                    ranges = {(p["p"],p["q"]):p["arms"][arm]["horizons"][h_index]["range"]
                              for p in panel["pairs"]}
                    ratios = []
                    for a,b,c in combinations(range(4),3):
                        ratios.append((ranges[a,b]*ranges[b,c]*ranges[c,a])
                                      /(ranges[b,a]*ranges[c,b]*ranges[a,c]))
                    arms[arm] = ratios
                rows.append(dict(probe=run["probe"], contrast=run["contrast"],
                                 left_odd=run["left_odd"], step=panel["step"],
                                 horizon=horizon, triangles_in_lexical_order=arms))
    return rows


if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
