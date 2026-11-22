"""Fresh-word cost of a coupling that agrees with the scheduled one-load contact."""

from collections import Counter
from fractions import Fraction
import json
from math import log, pi, sin, sqrt

import numpy as np
from scipy.linalg import expm
import sympy as sp

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import offer_service
from .modular_standard import stationary_operator

SHOTS = 2**24
# Per-category probability allowance for dyadic cumulative cuts and numerical pulse columns.
SAMPLER_ERROR = 1e-10


def sampling_weights(prob):
    """Round cumulative cuts to reproducible dyadic probabilities before the three-way draw."""
    if np.any(prob < 0):
        raise ValueError("A terminal probability is negative")
    np.testing.assert_allclose(prob.sum(), 1, atol=4e-12, rtol=0)
    scale = 2**36
    cuts = np.rint(scale*np.cumsum(prob/prob.sum())[:-1]).astype(np.int64)
    weights = np.diff(np.r_[0, cuts, scale])/scale
    np.testing.assert_allclose(weights, prob, atol=SAMPLER_ERROR, rtol=0)
    return weights


def portable_counts(cuts, arm, menu):
    """PCG64 high-36-bit categorical sampler, version 1; no floating point draws."""
    scale = 2**36
    if (len(cuts) != 2 or any(type(v) is not int for v in cuts)
            or not 0 <= cuts[0] <= cuts[1] <= scale):
        raise ValueError("Invalid three-port integer cumulative cuts")
    rng = np.random.PCG64([20261123, arm, menu])
    first = below_second = 0
    for start in range(0, SHOTS, 2**18):
        ranks = rng.random_raw(min(2**18, SHOTS - start)) >> np.uint64(28)
        first += int(np.count_nonzero(ranks < cuts[0]))
        below_second += int(np.count_nonzero(ranks < cuts[1]))
    return [first, below_second - first, SHOTS - below_second]


def native_triangle(hamiltonian, duration):
    """Two terminal measurements, each on a fresh sharp input at the first port."""
    law = candidate.native_law
    single = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    empty = law.SRC.index((0, 0))
    pulses = {}
    for t in (pi / 2, duration, 3 * pi / 2):
        unitary = expm(-1j * t * hamiltonian)
        code = np.eye(law.DIM)[:, single]
        block = code.T @ unitary @ code
        np.testing.assert_allclose(unitary @ code, code @ block, rtol=0, atol=2e-12)
        np.testing.assert_allclose(unitary[:, empty], np.eye(law.DIM)[:, empty], rtol=0, atol=2e-12)
        pulses[t] = block
    results = []
    for partner in (1, 2):
        state = np.array([1, 0, 0], complex)
        for pair, t in (([0, 1], pi / 2), ([1, 2], duration), ([0, partner], 3 * pi / 2)):
            state[pair] = pulses[t] @ state[pair]
        results.append(abs(state)**2)
    return np.array(results)


def triangle_forecast(strength, duration):
    z = np.exp(-1j * strength * duration)
    swap_phase = np.exp(-1j * duration)
    amplitudes = [(1 + swap_phase) / 2, (1 - swap_phase) / 2]
    rows = []
    for j, f in enumerate(amplitudes):
        row = np.zeros(3)
        row[0] = abs(1 + z * f)**2 / 4
        row[j + 1] = abs(1 - z * f)**2 / 4
        row[2 - j] = abs(amplitudes[1 - j])**2 / 2
        rows.append(row)
    return np.array(rows)


def word_masses(weights, witnesses, reader):
    zero, unit, double = weights
    even, odd = {"common": (zero, unit), "fair": (1, sp.Rational(1, 2))}[reader]
    return (zero**witnesses * even, unit**witnesses * even, unit**witnesses * odd,
            double**witnesses * even, double**witnesses * odd)


def mixture_exponent(masses, t):
    """Logarithmic derivative on a positive mixture interval; test the integration constant."""
    zero, even, odd = masses[:3]
    r, s = sp.cancel(even / odd), sp.cancel(zero / even)
    if sp.diff(r, t) == 0:
        return None
    k = sp.cancel(sp.diff(s, t) * r / (s * sp.diff(r, t)))
    if not k.is_Rational or k <= 0:
        return None
    # The derivative fixes the exponent. One interior value fixes the constant.
    middle = sp.Rational(1, 2)
    if s.subs(t, middle)**k.q != r.subs(t, middle)**k.p:
        return None
    return k


def source_check():
    t = sp.Symbol("mixture", real=True)
    starts = ((1, 1, 1), (7, 1, 1))
    words = [offer_service.signed_probabilities(1, counts=s) for s in starts]
    mix = [(1-t)*a + t*b for a, b in zip(*words, strict=True)]
    out = {}
    for m, reader in ((7, "common"), (8, "common"), (8, "fair"), (16, "common")):
        masses = word_masses(mix, m, reader)
        for arm, start in enumerate(starts):
            native = offer_service.finite_service(1, m, counts=start, reader=reader)
            admitted = native["conditional_on_filling_pointer_flip"]
            if tuple(x.subs(t, arm) for x in masses) != tuple(admitted*x for x in native["modes"]):
                raise ArithmeticError("The five joint probabilities differ from the finite offer service")
        degree = max(sp.degree(x, t) for x in masses)
        k = mixture_exponent(masses, t)
        if k != (m if reader == "common" else None) or degree != m + int(reader == "common"):
            raise ArithmeticError("Wrong source degree or mixture response")
        zero, even, _, two_even, _ = masses
        held = sp.cancel((two_even * zero**3 / even**4).subs(t, 1))
        if held == 1:
            raise ArithmeticError("The perturbed preparation lost its higher-mode failure")
        out[f"{m}_{reader}"] = dict(degree=int(degree), exponent=None if k is None else int(k),
                                    held_two_unit_ratio=float(held))
    return out


def finite_resolution(budget):
    candidates = {Fraction(p, q) for q in range(1, budget) for p in range(1, budget - q + 1)}
    aligned = sorted(k for k in candidates if k % 8 == 0)
    # sin²(pi k/8) depends only on distance to the nearest multiple of eight.
    distances = {k: min(k % 8, 8 - k % 8) for k in candidates if k not in aligned}
    distance = min(distances.values())
    return aligned, distance, sorted(k for k, d in distances.items() if d == distance), len(candidates)


def reflecting_bank(strength, up, down):
    """Check detailed balance on levels 0 through strength+1, including both reflecting ends."""
    n = strength + 2
    transition = sp.zeros(n)
    for level in range(n):
        transition[level, min(level+1, n-1)] += up
        transition[level, max(level-1, 0)] += down
        transition[level, level] += 1-up-down
    weights = sp.Matrix([(up / down)**j for j in range(n)])
    if weights.T * transition != weights.T:
        raise ArithmeticError("The proposed bank law is not stationary")
    return weights[-1]/weights[1], weights[1]/weights[0]


def run():
    board = federation(0)
    face = sorted(board.faces())[0]
    budgets = {}
    for level in range(3):
        carrier = federation(level)
        counts = Counter(len(offer_service.incident_witnesses(carrier, e)) + 1 for e in carrier.seams)
        if not set(counts) <= {9, 10, 11}:
            raise ValueError("The available one-pass record budgets changed")
        budgets[carrier.cells] = dict(sorted(counts.items()))
    exchange, contrast = candidate.native_law.local_term_parts()
    derived, _ = stationary_operator()
    n = max(max(counts) for counts in budgets.values())
    aligned, distance, rivals, size = finite_resolution(n)
    if aligned != [8] or distance != Fraction(1, 10) or rivals != [Fraction(1, 10)]:
        raise ArithmeticError("The finite positive alias class changed")
    strength, = aligned
    np.testing.assert_allclose(derived, exchange + float(strength)*contrast, rtol=0, atol=2e-13)
    rows = {}
    for k in (-8, 0, float(distance), 4, 7, int(strength), 9, 10, 16, 128):
        data = []
        for t in (pi/4, pi/(4*sqrt(2))):
            native = native_triangle(exchange + k*contrast, t)
            np.testing.assert_allclose(native, triangle_forecast(k, t), rtol=0, atol=4e-12)
            np.testing.assert_allclose(1.25-native[:,0].sum(), sin(k*t/2)**2, rtol=0, atol=4e-12)
            data.append(native.tolist())
        rows[str(k)] = data
    # Both materials can host the protected eight-witness standard. Their
    # scheduled pulses coincide. A separate one-load read resolves the materials.
    np.testing.assert_allclose(rows["8"][0], rows["16"][0], rtol=0, atol=4e-12)
    tv = abs(np.array(rows["8"][1]) - np.array(rows["16"][1])).sum(axis=1)/2
    if max(tv) < 0.1:
        raise ArithmeticError("The second pulse failed to distinguish the material counterexample")
    # Sample the material derived operator, with every terminal address retained.
    acquired = native_triangle(derived, pi/4)
    counts = [np.random.default_rng(seed).multinomial(SHOTS, sampling_weights(prob)).tolist()
              for seed, prob in zip((74504, 74505), acquired, strict=True)]
    radius = sqrt(log(100) / SHOTS) + 2*SAMPLER_ERROR
    deficit = 1.25-sum(row[0] for row in counts)/SHOTS
    threshold = sin(float(distance)*pi/8)**2
    if deficit + radius >= threshold:
        raise ValueError("The two readings do not resolve the nearest bounded-budget rival")
    for up, down in [(sp.Rational(1,3), sp.Rational(1,4)), (sp.Rational(1,2), sp.Rational(1,5))]:
        gap, parity = reflecting_bank(16, up, down)
        if gap != parity**16:
            raise ArithmeticError("The persistent bank lost the alternative exponent")
    return dict(triangle=list(map(int, face)), word_budgets=budgets, word_sources=source_check(),
                positive_candidates=size, selected=int(strength), nearest_rival=str(rivals[0]),
                probabilities=rows, shots_per_arm=SHOTS, terminal_counts=counts,
                alternate_material_second_pulse_tv=tv.tolist(),
                observed_deficit=deficit, one_sided_99pct_radius=radius, nearest_competitor_deficit=threshold)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
