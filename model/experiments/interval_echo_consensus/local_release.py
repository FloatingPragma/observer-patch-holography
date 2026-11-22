"""A four-port release using blank Interval Echo contacts and signed vacancy reads."""
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import cell, federation

PROTOCOL = Path(__file__).with_name("local_release_protocol.json")
SUPPORT = (4, 5, 6, 7)
VACANCY = ((0, 1), (2, 3), (8, 9), (10, 11))


def layers(matchings=None):
    if matchings is None:
        matchings = json.loads(PROTOCOL.read_text(encoding="utf-8"))["matchings"]
    seams = set(cell().seams)
    result = []
    for matching in matchings:
        if sorted(v for edge in matching for v in edge) != list(range(12)):
            raise ValueError("The layer must cover each port exactly once")
        perm = np.arange(12)
        for a, b in matching:
            if (a, b) not in seams:
                raise ValueError("The preparation addresses a missing native seam")
            perm[a], perm[b] = b, a
        result.append(perm)
    return result


def native_instruments():
    """Check the local reduction, and all three signed echo instrument maps on 00,10,01."""
    basis = np.eye(law.DIM)[:, [law.SRC.index(word) for word in ((0, 0), (1, 0), (0, 1))]]
    echo = law.echo()
    written = echo @ basis
    contrast = np.repeat([3] + [d for _, d in law.LABELS], law.N_SRC)
    for sign in (-2, -1, 0, 1, 2, 3):
        out = echo @ ((contrast == sign)[:, None] * written)
        np.testing.assert_allclose(out, basis * (np.array([0, 1, -1]) == sign), atol=1e-14, rtol=0)
    h0, h1 = law.local_term_parts()
    exchange = np.zeros((3, 3))
    exchange[1:, 1:] = np.array([[1, -1], [-1, 1]]) / 2
    np.testing.assert_allclose(h0 @ basis, basis @ exchange, atol=1e-14, rtol=0)
    np.testing.assert_allclose(h1 @ basis, basis @ np.diag([0, 1, 1]), atol=1e-14, rtol=0)


def prepare():
    """Gaussian integer numerators. Every half layer divides by 1+i up to a common phase."""
    if {v for edge in VACANCY for v in edge} != set(range(12)).difference(SUPPORT):
        raise ValueError("The vacancy readers must cover the complement")
    if not set(VACANCY).issubset(cell().seams):
        raise ValueError("Missing vacancy reader seam")
    re, im = np.eye(12, dtype=int)[5], np.zeros(12, dtype=int)
    plan = layers()
    for p in plan:
        re, im = re - im[p], im + re[p]
    denominator = 2**len(plan)
    weights = re**2 + im**2
    inside = list(SUPPORT)
    signed = np.zeros(12, dtype=int)
    signed[inside] = [1, -1, -1, 1]
    if (weights.sum() != denominator or weights[inside].sum()*4 != denominator
            or np.any(im[inside]) or not np.array_equal(re[inside], -signed[inside])):
        raise ArithmeticError("The native preparation failed its signed herald")
    outcomes = {str(v): int(weights[v]) for edge in VACANCY for v in edge}
    return signed, {"vacancy": int(weights[inside].sum()), "occupied": outcomes, "denominator": denominator}


def retry_budget(herald, attempts):
    """Route each observed failure back to the source on ordinary sharp pi exchanges."""
    if type(attempts) is not int or attempts < 1:
        raise ValueError("Expected a positive attempt budget")
    board = cell()
    routes = {5: [5]}
    queue = [5]
    for port in queue:
        for neighbor in board.neighbours[port]:
            if neighbor not in routes:
                routes[neighbor] = [neighbor, *routes[port]]
                queue.append(neighbor)
    paths = {v: routes[int(v)] for v, n in herald["occupied"].items() if n}
    failure = Fraction(sum(herald["occupied"].values()), herald["denominator"])
    if failure + Fraction(herald["vacancy"], herald["denominator"]) != 1:
        raise ArithmeticError("Lost a herald branch")
    return {"attempts": attempts, "all_failed": str(failure**attempts),
            "prepared": str(1 - failure**attempts), "failure_routes": paths,
            "maximum_pi_swaps": (attempts-1)*max(len(p)-1 for p in paths.values()),
            "maximum_half_contacts": 24*attempts, "maximum_signed_queries": 4*attempts}


def moving_profile():
    """Two native half layers and a zero-contrast herald create an internally evolving state."""
    matchings = (
        ((0, 1), (2, 4), (3, 6), (5, 9), (7, 8), (10, 11)),
        ((0, 1), (2, 6), (3, 4), (5, 9), (7, 8), (10, 11)),
    )
    local = np.eye(12, dtype=complex)[4]
    for perm in layers(matchings):
        local = (local + 1j*local[perm])/(1+1j)
    signed = np.eye(12, dtype=int)[4] - np.eye(12, dtype=int)[6]
    # The rejection addresses 2 and 3 share a seam, so one signed query suffices.
    failed = {str(v): abs(local[v])**2 for v in (2, 3)}
    local[[2, 3]] = 0
    np.testing.assert_allclose(local, -1j*signed/2, rtol=0, atol=1e-14)
    result = {}
    for level, home in ((1, 5), (2, 8)):
        b = federation(level)
        boundary = sorted(v % 12 for edge in b.external for v in edge if v//12 == home)
        lap = cell().laplacian().toarray().astype(int)
        # The minimal polynomial gives an exact all-time test from four boundary moments.
        eye = np.eye(12, dtype=int)
        if np.any(lap @ (lap-6*eye) @ (lap@lap-10*lap+20*eye)):
            raise ArithmeticError("Unexpected primitive spectrum")
        word = signed
        for _ in range(4):
            if np.any(word[boundary]):
                raise ArithmeticError("A free flight can reach the gateways")
            word = lap @ word
        start = np.zeros(b.n, complex)
        start[12*home:12*(home+1)] = signed/np.sqrt(2)
        rows = {}
        for t in (0, 1, 7, 37):
            state = expm_multiply((-.5j*t)*b.laplacian(), start)
            rows[t] = record(state, home).tolist()
        if rows[37][4]+rows[37][6] >= .01 or max(v[-1] for v in rows.values()) > 1e-20:
            raise ArithmeticError("The internal movement or exterior-quiet control failed")
        result[b.cells] = {"home": home, "conditional_ports_and_quiet": rows}
    return {"success": .5, "failed_addresses": failed, "releases": result}


def all_cells(signed):
    if (signed.shape != (12,) or not np.issubdtype(signed.dtype, np.integer)
            or sorted(signed.tolist()) != [-1, -1, *([0]*8), 1, 1]):
        raise ValueError("Expected two positive and two negative unit amplitudes")
    inventory = {}
    for level in range(3):
        b = federation(level)
        # Column c is supported only in cell c. Check every external row as well.
        columns = sparse.kron(sparse.eye(b.cells, dtype=int), signed[:, None], format="csr")
        native = b.laplacian()
        if not np.array_equal(native.data, np.rint(native.data)):
            raise ValueError("Expected the integer native Laplacian")
        lap = native.astype(np.int64)
        if (lap @ columns - 6*columns).nnz:
            raise ArithmeticError("A port on the full federation leaks or changes phase")
        inventory[b.cells] = columns.shape[1]
    return inventory


def record(state, home, reverse=False):
    local = state[12*home:12*(home+1)].copy()
    if reverse:
        for p in reversed(layers()):
            local = (local - 1j*local[p])/(1-1j)
    quiet = np.vdot(state[:12*home], state[:12*home]).real
    quiet += np.vdot(state[12*(home+1):], state[12*(home+1):]).real
    out = np.r_[abs(local)**2, quiet]
    np.testing.assert_allclose(out.sum(), 1, rtol=0, atol=2e-12)
    return out


def appointments():
    native_instruments()
    signed, herald = prepare()
    coverage = all_cells(signed)
    home = 13
    b = federation(0)
    support = 12*home + np.array(SUPPORT)
    start = np.zeros((b.n, 1), complex)
    start[support, 0] = signed[list(SUPPORT)] / 2
    # Each resolved preparation address is an orthogonal retained outcome.
    resolved = np.eye(b.n, dtype=complex)[:, support]
    h = b.laplacian()/2
    a, d = 12*home + 4, 12*home + 9
    if (a, d) not in b.seams:
        raise ValueError("Missing rate-control seam")
    edge = sparse.csr_matrix(([1., -1.], ([a, d], [0, 0])), shape=(b.n, 1))
    arms = (("coherent", start, h), ("resolved", resolved, h),
            ("one_exchange_disabled", start, h - edge @ edge.T/2))
    report = {}
    for name, initial, generator in arms:
        table = {}
        for t in (0, 1, 7, 37):
            states = expm_multiply(-1j*t*generator, initial) if t else initial
            direct, inverse = (
                np.array([record(col, home, rev) for col in states.T]) for rev in (False, True))
            prob = direct.mean(axis=0)[list(SUPPORT)].sum()
            if name == "coherent":
                np.testing.assert_allclose([prob, inverse[0, 5]], [1, .25], atol=2e-12, rtol=0)
            elif t == 1 and prob >= .99:
                raise ArithmeticError("Control lost less than one percent at exposure one")
            table[t] = {"addresses_and_quiet": direct.tolist(), "inverse_addresses_and_quiet": inverse.tolist()}
        report[name] = table
    return {"prepared": herald, "bounded_preparation": retry_budget(herald, 16),
            "moving_profile": moving_profile(),
            "cells_checked": coverage, "conditional_on_vacancy": report}


if __name__ == "__main__":
    print(json.dumps(appointments(), indent=2, allow_nan=False))
