"""Independent raw-board replay of the connected two-load experiment."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import acosh, ceil, exp, log, log1p, sinh
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import jv

from tools.strict_json import loads
from ..submission.row_04_check import BOARD_HASH


BOARD = Path(__file__).resolve().parents[4]/"executable_model/carriers/FEDERATION_LEVEL0.npz"
ARCHIVE = Path(__file__).with_name("connected_records.json")
STRENGTHS, TIMES = (0, .5, 8), (1, 7, 37)
ARMS = {"coherent", "dephased", "unfiltered", "stochastic_swap"}
FAIL = Fraction(364, 729)
RETRY = sum(FAIL**attempt for attempt in range(8))
VECTOR = np.array([-13, 14])/27  # signed amplitudes after three clear readers
TOLERANCE = 2e-10


def board_generator():
    """Assemble all binary swaps seam by seam directly from the frozen carrier."""
    if sha256(BOARD.read_bytes()).hexdigest() != BOARD_HASH:
        raise ValueError("Unexpected source board")
    with np.load(BOARD, allow_pickle=False) as board:
        edges = sorted(tuple(sorted(map(int, edge))) for edge in board["edges"])
    n = 240
    pairs = np.array(list(combinations(range(n), 2)))
    lookup = np.full((n, n), -1)
    lookup[tuple(pairs.T)] = lookup[tuple(pairs[:, ::-1].T)] = np.arange(len(pairs))
    a, b, c, exit_port = 157, 156, 166, 163
    openings = {tuple(sorted((tip, exit_port))) for tip in (a, c)}
    cut = {edge for edge in edges if len({a,b,c}.intersection(edge)) == 1} - openings
    home = [edge for edge in edges if edge[0]//12 == edge[1]//12 == 13]
    if len(cut) != 11 or len(home) != 30:
        raise ArithmeticError("Changed neutral boundary or receiver")
    sources, targets = [], []
    for x, y in edges:
        if (x,y) in cut:
            continue
        stationary = np.delete(np.arange(n), [x,y])
        sources.append(lookup[x, stationary])
        targets.append(lookup[y, stationary])
    x, y = np.concatenate(sources), np.concatenate(targets)
    rows, columns = np.r_[x,y], np.r_[y,x]
    moves = sparse.csr_matrix((np.ones(len(rows)), (rows,columns)), shape=(len(pairs),)*2)
    counts = np.asarray(moves.sum(axis=1)).ravel()
    exchange = ((sparse.diags(counts)-moves)/2).tocsr()
    reference = 6-np.bincount(np.ravel(edges), minlength=n)
    contrast = counts + reference[pairs[:, 0]] + reference[pairs[:, 1]]
    launch = np.zeros((len(pairs), 2))
    launch[[lookup[a,b], lookup[c,b]], [0,1]] = 1
    # Visit the three vacancies in order and undo both exchanges. Only the
    # first visit to each vacancy writes its incident neutral records.
    setup = {"neutral_edges": [edge for hole in (c,b,a) for edge in sorted(cut) if hole in edge],
             "live_swaps": [tuple(sorted(edge)) for edge in ((b,c), (a,b), (a,b), (b,c))],
             "final_pair": tuple(sorted((a,b)))}
    return exchange, contrast, launch, home, lookup[tuple(np.array(home).T)], setup


def flight_amplitudes(exchange, contrast, kappa, times, launch, detector):
    """Chebyshev–Bessel propagation with only the detector coordinates accumulated."""
    offset = kappa*contrast
    low, high = min(offset), max(2*exchange.diagonal()+offset)
    radius, center = (high-low)/2, (high+low)/2
    scaled = (exchange+sparse.diags(offset-center, format="csr"))/radius
    argument = max(times)*radius
    degree = ceil(argument)+1
    # The generating function on |z|=exp(s) gives the absolute tail bound
    # 2 exp(x sinh(s) - (degree+1)s)/(1-exp(-s)).
    while True:
        s = acosh((degree+1)/argument)
        if log(2)+argument*sinh(s)-(degree+1)*s-log1p(-exp(-s)) < log(1e-17):
            break
        degree += 1
    orders = np.arange(degree+1)
    weights = 2*(-1j)**(orders % 4)*jv(orders, radius*np.array(times)[:, None])
    weights[:, 0] /= 2
    samples = np.empty((degree+1, len(detector), 2))
    prev, current = launch, scaled @ launch
    samples[:2] = prev[detector], current[detector]
    for i in range(2, degree+1):
        prev, current = current, 2*(scaled @ current)-prev
        samples[i] = current[detector]
    return np.einsum("tn,nij->tij", weights, samples)


def stochastic_population(exchange, times, launch, detector):
    """Poisson uniformization of the continuous-time rate-half seam swaps."""
    rate = 6
    step = sparse.eye(exchange.shape[0], format="csr")-exchange/rate
    if step.data.min() < 0 or not np.allclose(np.asarray(step.sum(axis=0)), 1, atol=1e-15, rtol=0):
        raise ArithmeticError("Invalid uniformized stochastic matrix")
    state = launch @ (float(RETRY)*VECTOR**2)
    means = rate*np.array(times)
    weights = np.exp(-means)
    result = weights[:, None]*state[detector]
    n = 0
    while True:
        n += 1
        weights *= means/n
        state = step @ state
        result += weights[:, None]*state[detector]
        # Bound all omitted Poisson weights by a geometric series.
        if n > max(means) and max(weights*means/(n+1-means)) < 1e-17:
            return result


def independently_replay(strengths, times):
    exchange, contrast, launch, home, detector, setup = board_generator()
    stochastic = stochastic_population(exchange, times, launch, detector)
    records = {}
    for kappa in strengths:
        amplitudes = flight_amplitudes(exchange, contrast, kappa, times, launch, detector)
        for i, time in enumerate(times):
            a = amplitudes[i]
            bins = {"coherent": float(RETRY)*abs(a @ VECTOR)**2,
                    "dephased": (abs(a)**2) @ (float(RETRY)*VECTOR**2),
                    "unfiltered": abs(a[:,0])**2, "stochastic_swap": stochastic[i]}
            # The complement includes every missed detector and exhausted preparation.
            records[kappa, time] = {arm: np.r_[v, 1-v.sum()] for arm,v in bins.items()}
    return setup, home, records


def appointments(entries, expected_keys):
    rows = {}
    for entry in entries:
        key = entry["kappa"], entry["exposure"]
        if (type(key[0]) not in (int, float) or type(key[1]) not in (int, float)
                or key in rows or set(entry["bins"]) != ARMS):
            raise ValueError("Invalid appointment, duplicate or wrong flight arms")
        rows[key] = entry["bins"]
    if set(rows) != expected_keys:
        raise ValueError("Missing or unexpected appointment")
    return rows


def archive_bins(archive, home):
    unit = archive["units_per_probability"]
    if (type(unit) is not int or unit != 10**12
            or json.dumps(archive["outcome_order"]) != json.dumps([*home, "quiet"])):
        raise ValueError("Wrong archive unit or receiver")
    rows = appointments(archive["appointments"], {(k,t) for k in STRENGTHS for t in TIMES})
    for record in rows.values():
        for values in record.values():
            if (not isinstance(values, list) or len(values) != 31
                    or any(type(n) is not int or not 0 <= n <= unit for n in values)):
                raise ValueError("Invalid archived probability units")
    return {key: {arm: np.asarray(v)/unit for arm,v in record.items()} for key,record in rows.items()}


def max_bin_error(actual, expected):
    maximum = 0.
    for key, arms in expected.items():
        for arm, probabilities in arms.items():
            values = actual[key][arm]
            if isinstance(values, (list, tuple)) and any(type(v) not in (float, int) for v in values):
                raise ValueError("Expected real numeric terminal probabilities")
            received = np.asarray(values)
            if (received.shape != (31,) or received.dtype.kind not in "fiu"
                    or not np.isfinite(received).all() or received.min() < 0
                    or abs(received.sum()-1) > TOLERANCE):
                raise ValueError("Invalid terminal distribution")
            delta = float(np.max(np.abs(received-probabilities)))
            if not np.isfinite(delta) or delta > TOLERANCE:
                raise ArithmeticError(f"Independent replay disagrees at {key}, {arm} by {delta}")
            maximum = max(maximum, delta)
    return maximum


def verify(report, strengths, times, archive):
    setup, home, expected = independently_replay(strengths, times)
    # JSON comparison identifies tuples and lists while distinguishing real
    # integer addresses from booleans and floats that compare equal in Python.
    if json.dumps(report["setup"], sort_keys=True) != json.dumps(setup, sort_keys=True):
        raise ValueError("Wrong contact installation")
    if (json.dumps(report["home_seams_in_order"]) != json.dumps(home)
            or report["guaranteed_joint_home"] != str(RETRY/2)):
        raise ValueError("Wrong receiver or protected weight")
    for field, values in (("first_vacancy_in_one_attempt", [4/9**i for i in range(1,4)]),
                          ("attempt_reach", [float(FAIL**i) for i in range(8)]),
                          ("exhaustion", float(FAIL**8))):
        supplied = report[field]
        elements = supplied if isinstance(supplied, (tuple, list)) else (supplied,)
        if any(type(v) not in (int, float) for v in elements) or np.shape(supplied) != np.shape(values):
            raise ValueError("Invalid native preparation probabilities")
        np.testing.assert_allclose(supplied, values, rtol=0, atol=5e-14, equal_nan=False)
    raw = appointments(report["appointments"], set(expected))
    return {"checked_bins": len(expected)*4*31, "max_producer_error": max_bin_error(raw, expected),
            "max_archive_error": max_bin_error(archive_bins(archive, home), expected)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="check all nine appointments")
    parser.add_argument("--receipt", type=Path, help="verify saved producer output of matching scope")
    args = parser.parse_args()
    strengths, times = (STRENGTHS, TIMES) if args.full else ((.5,), (1,))
    if args.receipt:
        report = loads(args.receipt.read_text(encoding="utf-8-sig"))
    else:
        from . import open_pair, open_pair_release
        rho, report = open_pair_release.source()
        report["home_seams_in_order"] = open_pair.geometry()[1]
        report["appointments"] = [dict(kappa=k, exposure=t, bins=open_pair_release.terminal(
            rho, report["exhaustion"], k, t)) for k in strengths for t in times]
    print(json.dumps(verify(report, strengths, times, loads(ARCHIVE.read_text())), allow_nan=False))


if __name__ == "__main__":
    main()
