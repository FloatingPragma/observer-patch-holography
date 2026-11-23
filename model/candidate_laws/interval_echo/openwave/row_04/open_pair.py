"""Two-load interference at an open exit on the scored 20-cell federation."""
from functools import lru_cache
from math import factorial, isfinite
from numbers import Real

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import Carrier
from .corridor_pair import certify_native_collar, seal_path


PATH = (157, 156, 166)
BRIDGE = 163
BASIS = ((156, 157), (157, 166), (156, 166))


@lru_cache(maxsize=1)
def geometry():
    """Native neutral seams and the unmodified geometric receiver addresses."""
    board = federation(0)
    a, b, c = PATH
    # Start on the M10.3 sharp seam. Its degree-five endpoint owns the
    # reference. Choose the other nonadjacent tip with the gateway's degree.
    peers = [v for v in board.neighbours[b] if v != a and v not in board.neighbours[a]
             and board.degree[v] == board.degree[a]]
    common = set(board.neighbours[a]) & set(board.neighbours[c])
    if (peers != [c] or common != {b, BRIDGE}
            or tuple(board.degree[list(PATH)]) != (6, 5, 6)):
        raise ValueError("Unexpected source path, exit or native reference at its center")
    protected = set(PATH)
    openings = {tuple(sorted((v, BRIDGE))) for v in (a, c)}
    cut = tuple(edge for edge in board.seams if len(protected.intersection(edge)) == 1 and edge not in openings)
    home = tuple(edge for edge in board.seams if edge[0]//12 == edge[1]//12 == 13)
    if len(cut) != 11 or len(home) != 30:
        raise ValueError("Incomplete neutral boundary or scored pair detector")
    return cut, home


def seal_from_sharp():
    return seal_path(PATH, geometry()[0])


@lru_cache(maxsize=1)
def source_parts():
    """Two addressed preparation seams and the native empty reference at b."""
    board = federation(0)
    exchange, contrast = np.zeros((2, 3, 3))
    one = candidate.native_law.seam_term(1)
    kinetic = candidate.native_law.seam_term(0)
    for i, pair in enumerate(BASIS):
        for edge in (BASIS[0], BASIS[2]):
            word = tuple(int(v in pair) for v in edge)
            j = candidate.native_law.SRC.index(word)
            exchange[i, i] += kinetic[j, j]
            contrast[i, i] += one[j, j]-kinetic[j, j]
            if sum(word) == 1:
                after = tuple(sorted(set(pair)^set(edge)))
                dest = candidate.native_law.SRC.index(word[::-1])
                exchange[BASIS.index(after), i] += kinetic[dest, j]
        contrast[i, i] += (6-board.degree[PATH[1]])*int(PATH[1] in pair)
    if not np.array_equal(contrast, 2*np.eye(3)):
        raise ArithmeticError("The prepared three-state source has a relative contrast phase")
    exchange.setflags(write=False)
    contrast.setflags(write=False)
    return exchange, contrast


@lru_cache(maxsize=1)
def full_parts():
    """Full distinct-port sector, exchange, native contrast, detector and source indices."""
    certify_native_collar()
    board = federation(0)
    cut, home = geometry()
    sealed = set(cut)
    # A dirty d=0 contact contributes zero for every endpoint state; the
    # blank contacts act natively. Reference terms belong to the original
    # ports, including when a material mark changes their active degree.
    live = Carrier("ports connected by blank contacts", board.n, [e for e in board.seams if e not in sealed])
    if connected_components(live.adjacency(), directed=False, return_labels=False) != 1:
        raise ArithmeticError("The prepared contacts disconnect the federation")
    sector = PairSector(live)
    exchange = sector.exchange()
    single_cost = live.degree + (6-board.degree)
    contrast = single_cost[sector.aa] + single_cost[sector.bb] - 2*sector.adjacent
    detector = np.array([sector.index(*e) for e in home])
    source = np.array([sector.index(*e) for e in BASIS])
    return sector, exchange, contrast, detector, source


def dark_certificate():
    """Verify exact cancellation including the common exit and its entire exterior."""
    sector, exchange, contrast, detector, source = full_parts()
    if not np.array_equal(2*exchange.data, np.rint(2*exchange.data)):
        raise ArithmeticError("Non-half-integral binary exchange")
    kernel = (2*exchange).astype(np.int64)
    dark = np.zeros(sector.P, dtype=np.int64)
    dark[source[[0, 2]]] = (1, -1)
    if (not np.array_equal(kernel @ dark, 2*dark) or not np.array_equal(contrast*dark, 3*dark)
            or dark[detector] @ dark[detector] != 2):
        raise ArithmeticError("Native dark state leaks, changes phase nonuniformly or misses the receiver")
    bright = abs(dark)
    first = kernel @ bright
    if first[source[1]] != -2 or first[sector.index(PATH[1], BRIDGE)] != -2:
        raise ArithmeticError("The symmetric pair does not reach both exit configurations")
    return {"dimension": sector.P, "twice_exchange_eigenvalue": 2, "contrast_eigenvalue": 3,
            "dark_bins": (dark[detector]**2 / 2).tolist()+[0.]}


def bright_escape():
    """First strictly positive outside-cell probability coefficient of the symmetric pair."""
    sector, exchange, _, _, source = full_parts()
    # At the first distance reaching outside, no diagonal (including κ)
    # can contribute. Every shortest path has the same exchange sign.
    steps = (-2*exchange).astype(np.int64)
    steps.setdiag(0)
    steps.eliminate_zeros()
    if np.any(steps.data != 1):
        raise ArithmeticError("Unexpected binary exchange hopping weight")
    signs = np.zeros(sector.P, dtype=np.int64)
    signs[source[[0, 2]]] = 1
    outside = (sector.aa//12 != 13) | (sector.bb//12 != 13)
    for length in range(1, 12):
        if int(signs.max()) * int(np.diff(steps.indptr).max()) > np.iinfo(np.int64).max:
            raise ArithmeticError("Shortest-path count exceeds the integer accumulator")
        signs = steps @ signs
        if np.any(signs[outside]):
            numerator = sum(int(v)**2 for v in signs[outside] if v)
            return length, numerator, 2*(2**length * factorial(length))**2
    raise ArithmeticError("The symmetric pair has no certified route out of its original home")


def read_release(kappa, time, source_density, *, flight="coherent"):
    """Thirty geometric seam bins and quiet for an unnormalized three-state release."""
    if (any(not isinstance(v, Real) or isinstance(v, (bool, np.bool_)) or not isfinite(v)
            for v in (kappa, time)) or time < 0):
        raise ValueError("Expected finite real strength and nonnegative exposure")
    rho = np.asarray(source_density, dtype=complex)
    if (rho.shape != (3, 3) or not np.isfinite(rho).all()
            or not np.allclose(rho, rho.conj().T, atol=1e-13, rtol=0)
            or np.linalg.eigvalsh(rho)[0] < -2e-13 or rho.trace().real > 1+2e-12):
        raise ValueError("Expected the unnormalized accepted density on the source words")
    sector, exchange, contrast, detector, source = full_parts()
    if flight == "stochastic_swap":
        diag = np.zeros(sector.P)
        diag[source] = rho.diagonal().real
        populations = expm_multiply((-time)*exchange, diag)
    elif flight in ("coherent", "dephased"):
        if flight == "dephased":
            rho = np.diag(rho.diagonal())
        populated = np.flatnonzero(np.any(rho != 0, axis=0) | np.any(rho != 0, axis=1))
        if not len(populated):
            return np.zeros(31)
        initial = np.zeros((sector.P, len(populated)))
        initial[source[populated], np.arange(len(populated))] = 1
        # Removing the uniform energy translates only a common phase.
        h = exchange + sparse.diags(kappa*(contrast-3), format="csr", dtype=np.float64)
        columns = expm_multiply((-1j*time)*h, initial)
        reduced = rho[np.ix_(populated, populated)]
        populations = np.sum((columns @ reduced)*columns.conj(), axis=1).real
    else:
        raise ValueError("Unknown free-flight law")
    quiet = np.ones(sector.P, dtype=bool)
    quiet[detector] = False
    bins = np.r_[populations[detector], populations[quiet].sum()]
    if not np.isfinite(bins).all() or bins.min() < -2e-12 or abs(bins.sum()-rho.trace().real) > 2e-10:
        raise ArithmeticError("Incomplete or negative terminal law")
    return np.maximum(bins, 0.)


if __name__ == "__main__":
    import json
    print(json.dumps({"setup": seal_from_sharp(), "dark": dark_certificate(), "bright_escape": bright_escape()},
                     indent=2, allow_nan=False))
