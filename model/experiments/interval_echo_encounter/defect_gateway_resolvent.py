"""Exact gateway resolvent for a single signed defect on the blank, unit-filled host.

From the OPH repository root, run
python -m model.experiments.interval_echo_encounter.defect_gateway_resolvent
"""
from collections import Counter
from functools import cache
import json

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply, spsolve
import sympy as sp

from model.executable_model.carrier import CELL_SEAMS, automorphisms, cell, federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator


CELL = cell().adjacency().toarray().astype(int)


@cache
def gateway_chart(gates):
    """The sorted representative and order permutation for a physical triple of gateway labels."""
    def image(g):
        moved = tuple(g[p] for p in gates)
        return tuple(sorted(moved)), moved

    representative, moved = min(image(g) for g in cell_symmetries())
    return representative, tuple(representative.index(p) for p in moved)


@cache
def cell_symmetries():
    return automorphisms(CELL_SEAMS)


def local_blocks(gates):
    interior = sorted(set(range(12)) - set(gates))
    b = -CELL[np.ix_(gates, interior)] / 2
    kii = (5 * np.eye(9) - CELL[np.ix_(interior, interior)]) / 2
    kgg = 3 * np.eye(3) - CELL[np.ix_(gates, gates)] / 2
    return kgg, b, kii


@cache
def gateway_block(gates, sign, kappa, z):
    kgg, b, kii = local_blocks(gates)
    interior_energy = (6 + 2 * sign) * kappa
    correction = b @ np.linalg.solve((z - interior_energy) * np.eye(9) - kii, b.T)
    return 6 * kappa * np.eye(3) + kgg + correction


def gateway_operator(board, sign, kappa, z):
    """Return gateway addresses and H_eff(z); (z I - H_eff)^-1 is the exact gateway resolvent."""
    if sign not in (-1, 1):
        raise ValueError("The one-defect sector has sign -1 or +1")
    ports = np.flatnonzero(board.degree == 6)
    triples = np.split(ports, np.flatnonzero(np.diff(ports // 12)) + 1)
    if len(triples) != board.cells or any(len(x) != 3 for x in triples):
        raise ValueError("Every owner must have three gateways")
    blocks = []
    for triple in triples:
        chart, order = gateway_chart(tuple(int(p % 12) for p in triple))
        blocks.append(sparse.coo_matrix(gateway_block(chart, sign, kappa, z)[np.ix_(order, order)]))
    left, right = np.searchsorted(ports, np.asarray(board.external)).T
    ends = np.r_[left, right], np.r_[right, left]
    links = sparse.coo_matrix((np.full(len(ends[0]), -0.5), ends), shape=(len(ports), len(ports)))
    return ports, (sparse.block_diag(blocks) + links).tocsc()


def exact_mixing(sign):
    """Exact off-diagonal virtual coupling on the most frequent gateway chart at E=51."""
    gates = (0, 2, 8)
    interior = sorted(set(range(12)) - set(gates))
    a = sp.Matrix(CELL.tolist())
    b = -a.extract(gates, interior) / 2
    kii = (5 * sp.eye(9) - a.extract(interior, interior)) / 2
    sigma = b * ((51 - (6 + 2 * sign) * 8) * sp.eye(9) - kii).inv() * b.T
    assert sigma[0, 1] == sigma[0, 2] == sigma[1, 2]
    return sigma[0, 1]


def check_board(level):
    board = federation(level)
    ports = np.flatnonzero(board.degree == 6)
    direct = board.adjacency()[ports][:, ports]
    count, labels = connected_components(direct, directed=False)
    sizes = np.bincount(labels)
    assert direct.nnz == 2 * (len(ports) - count) and max(direct.sum(axis=1).A1) <= 2
    owners = [len(set(ports[labels == component] // 12)) for component in range(count)]
    home_path_owners = np.unique(ports[labels == labels[0]] // 12)
    outside = ~np.isin(np.arange(board.n) // 12, home_path_owners)
    launch = np.zeros(board.n)
    launch[ports[0]] = 1
    rng = np.random.default_rng(127 + level)
    sources = rng.normal(size=(len(ports), 3)) + 1j * rng.normal(size=(len(ports), 3))
    rhs = np.zeros((board.n, sources.shape[1]), dtype=complex)
    rhs[ports] = sources
    worst = 0.
    outside_by_sign = {}
    for sign in (-1, 1):
        for kappa, z in ((8., 51 + .2j), (1., 9 + .3j)):
            gates, effective = gateway_operator(board, sign, kappa, z)
            np.testing.assert_array_equal(gates, ports)
            native = defect_generator(board, sign, kappa, reference=0)
            full = spsolve((z * sparse.eye(board.n) - native).tocsc(), rhs)[ports]
            reduced = spsolve((z * sparse.eye(len(ports)) - effective).tocsc(), sources)
            worst = max(worst, float(np.linalg.norm(full - reduced) / np.linalg.norm(full)))
            if kappa == 8:
                flight = expm_multiply(-32j * native, launch)
                assert abs(np.vdot(flight, flight).real - 1) < 1e-10
                outside_by_sign["vacancy" if sign == -1 else "excess"] = float(np.linalg.norm(flight[outside])**2)
    assert worst < 1e-12
    return dict(cells=board.cells, full=board.n, reduced=len(ports),
                gateway_path_lengths=dict(sorted(Counter(map(int, sizes)).items())),
                maximum_owners_on_one_path=max(owners), relative_resolvent_error=worst,
                launch_port=int(ports[0]), launch_path_owners=home_path_owners.tolist(),
                outside_launch_path_owners_at_32=outside_by_sign)


if __name__ == "__main__":
    vacancy, excess = exact_mixing(-1), exact_mixing(1)
    assert (vacancy, excess) == (sp.Rational(286176, 10605875), -sp.Rational(237088, 6416685))
    boards = [check_board(level) for level in (0, 1, 2)]
    assert [b["maximum_owners_on_one_path"] for b in boards] == [4, 4, 7]
    np.testing.assert_allclose(
        [[b["outside_launch_path_owners_at_32"][s] for s in ("vacancy", "excess")] for b in boards],
        [[0.05584569, 0.08541939], [0.28738970, 0.44646610], [0.28906122, 0.45149586]],
        atol=1e-8, rtol=0)
    print(json.dumps(dict(boards=boards, vacancy_coupling=str(vacancy), excess_coupling=str(excess)), indent=2))
