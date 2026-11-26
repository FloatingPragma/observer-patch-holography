"""Which clock quadratures survive a free flight between disjoint blank dimers?"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from mpmath import iv

from ...full_algebra import native_precision as native
from ..measure import density_center

HERE = Path(__file__).resolve().parent
PORTS = (5, 9, 40, 45)
TIMES = (2, 3, 4)
N = 2**48
GRID = 2**51
SEED = 116403
Z = np.diag([1, -1])
Y = np.array([[0, -1j], [1j, 0]])
READ = (Z, Y)
I2 = np.eye(2)
PHASES = {"Z+": (1, 0), "Z-": (-1, 0), "Y+": (0, 1), "Y-": (0, -1),
          "held": (math.sqrt(0.5), math.sqrt(0.5)), "erased+": (0, 0), "erased-": (0, 0)}
# Upper bound for simultaneous event frequencies (all ternary effects, rest and free).
LOG = math.log(2 * 1024 / 0.01)


def outputs(kernel, state):
    """Trace out the two parked coins, retain both coherent receiver modes."""
    return np.einsum("ijab,bc,klac->ikjl", kernel, state, kernel.conj())


def native_flows():
    source, _, _, _ = density_center()
    board = native.iso.pairing.formation.federation(1)
    marks = native.iso.pairing.formation.region(board)
    assert all(b in board.neighbours[a] for a, b in (PORTS[:2], PORTS[2:]))
    ids = np.array([4 * p + j for p in PORTS for j in range(4)])
    coefficients = native.prior.weights()
    channels = {}
    for q in (-1, 1):
        for t in (0, *TIMES):
            channels[q, t] = np.zeros((2, 2, 2, 2), complex)
        for parity, rho in enumerate(source):
            h = native.generator(board, marks, q, bool(parity))
            blocks = native.kernel_panels(h, ids, coefficients)
            for t in TIMES:
                a, b = blocks[:, 2*t-1, 2:, :2]
                kernel = (np.asarray(a, float) + 1j * np.asarray(b, float)) / native.prior.FRAC
                if np.max(np.abs(kernel)) >= 0.05:
                    raise ArithmeticError("native block exceeds the contraction error enclosure")
                channels[q, t] += outputs(kernel, rho)
    return channels


def rest_flow(time):
    c, s = math.cos(time/2), math.sin(time/2)
    u = np.array([[c, 1j*s], [1j*s, c]])
    return np.einsum("ij,kl->ikjl", u, u.conj())


def score_probabilities(flow, phase):
    x, y = PHASES[phase]
    incident = (I2 + x*Z + y*Y) / 2
    rho = np.einsum("ikjl,jl->ik", flow, incident)
    return [[float(np.trace((I2 + sign*b) @ rho).real/2) for sign in (1, -1)] for b in READ]


def acquire(flow, rng):
    settings = {}
    for label in PHASES:
        counts = []
        for pair in score_probabilities(flow, label):
            grid = [round(p * GRID) for p in pair]
            grid.append(GRID - sum(grid))
            if min(grid) < 0:
                raise ArithmeticError("invalid native counting law")
            counts.append(rng.multinomial(N, np.array(grid) / GRID).tolist())
        settings[label] = counts
    return settings


def probability_interval(count):
    y = count / N
    a = LOG / N
    radius = math.sqrt(2*y*a) + 4*a
    # 2e-7 source HS uncertainty and a certified reference lower bound 0.02 I;
    # integer block enclosure, conversion, trace contraction and sampling grid.
    radius += 1.00002e-5 * (y + radius) + 5e-13
    return iv.mpf([max(0., math.nextafter(y-radius, -math.inf)), min(1., math.nextafter(y+radius, math.inf))])


def scores(bank):
    answer = []
    for plus, minus, missed in bank:
        if min(plus, minus, missed) < 0 or plus+minus+missed != N:
            raise ValueError("incomplete signed hit record")
        answer.append(probability_interval(plus) - probability_interval(minus))
    return answer


def limits(x):
    return [math.nextafter(float(x.a), -math.inf), math.nextafter(float(x.b), math.inf)]


def matrix_record(settings):
    obs = {k: scores(rows) for k, rows in settings.items()}
    m = [[(obs[k+"+"][i]-obs[k+"-"][i])/2 for k in ("Z", "Y")] for i in (0, 1)]
    det = m[0][0]*m[1][1]-m[0][1]*m[1][0]
    differences = []
    for i in (0, 1):
        b = (obs["erased+"][i]+obs["erased-"][i])/2
        predicted = b+(m[i][0]+m[i][1])/iv.sqrt(2)
        differences.append(limits(obs["held"][i]-predicted))
    # Invert the four calibrated columns against a separate unknown-input arm.
    # The zero-flight calibration has no inverse.
    phase = None
    if float(det.b) < 0 or float(det.a) > 0:
        residual = [obs["held"][i] - (obs["erased+"][i]+obs["erased-"][i])/2 for i in (0, 1)]
        x = (m[1][1]*residual[0] - m[0][1]*residual[1]) / det
        y = (m[0][0]*residual[1] - m[1][0]*residual[0]) / det
        phase = limits(iv.atan2(y, x))
    return dict(entries=[[limits(x) for x in row] for row in m], det=limits(det), held_phase=phase,
                erased_label_difference=[limits(a-b) for a, b in zip(obs["erased+"], obs["erased-"], strict=True)],
                held_observed_minus_forecast=differences)


def conditional_matrix(settings):
    """Remove the measured, phase-dependent efficiency of each receiving aperture."""
    conditional = {}
    for label in ("Z+", "Z-", "Y+", "Y-"):
        conditional[label] = []
        for plus, minus, _ in settings[label]:
            p, n = probability_interval(plus), probability_interval(minus)
            conditional[label].append((p-n)/(p+n))
    return [(conditional[k+"+"][i]-conditional[k+"-"][i])/2 for i in (0, 1) for k in ("Z", "Y")]


def rotation_invariant(first, second):
    a, b, c, d = first
    u, v, w, x = second
    # trace(second adj(first)) squared divided by the determinant product.
    return limits((u*d-v*c-w*b+x*a)**2 / ((a*d-b*c)*(u*x-v*w)))


def comparisons(records):
    result = {}
    for q in (-1, 1):
        observed = [conditional_matrix(row["counts"]) for row in records
                    if row["charge"] == q and row["mode"] == "free" and row["time"]]
        result[str(q)] = [rotation_invariant(a, b) for a, b in zip(observed, observed[1:])]
    return result


def make():
    flows = native_flows()
    rng = np.random.default_rng(SEED)
    records = []
    for q in (-1, 1):
        for mode in ("rest", "free"):
            for t in (0, *TIMES):
                settings = acquire(rest_flow(t) if mode == "rest" else flows[q, t], rng)
                records.append(dict(charge=q, mode=mode, time=t, counts=settings,
                                    measured=matrix_record(settings)))
    return dict(port_pairs=[list(PORTS[:2]), list(PORTS[2:])], shots=N, seed=SEED, rows=records, conditional_transfer=comparisons(records))


def replay():
    saved = json.loads((HERE / "results.json").read_text())
    flows = native_flows()
    for charge, bounds in comparisons(saved["rows"]).items():
        np.testing.assert_allclose(bounds, saved["conditional_transfer"][charge], atol=1e-14, rtol=1e-13)
    for row in saved["rows"]:
        mode, t, q = row["mode"], row["time"], row["charge"]
        measured = matrix_record(row["counts"])
        np.testing.assert_allclose(measured["det"], row["measured"]["det"], rtol=1e-13, atol=0)
        for label, matrix in row["counts"].items():
            predicted = score_probabilities(rest_flow(t) if mode == "rest" else flows[q, t], label)
            for j, counts in enumerate(matrix):
                for i, phat in enumerate(predicted[j]):
                    lo, hi = limits(probability_interval(counts[i]))
                    if not lo <= phat <= hi:
                        raise AssertionError((q, mode, t, label, i, j, phat, (lo, hi)))
    return saved


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE / "results.json").write_text(json.dumps(make(), indent=2) + "\n")
    else:
        replay()
        print("Verified the disjoint two-dimer clock records.")
