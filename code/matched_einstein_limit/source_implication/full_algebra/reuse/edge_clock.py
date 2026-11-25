"""Certify common timestamp information on one native edge for different charges."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import factorial, isqrt
from pathlib import Path

import numpy as np

from . import clock_history as clock, full_clock, measure

DEN = 8 * clock.F
ORDER = 104
SCALE = 2**76
MIDDLE = (1, 2)


def coordinates(effect):
    """Four times the Pauli/2 coefficients, plus the unsigned trace coordinate."""
    paulis = 2 * full_clock.station.iso.basis(False)
    moments = (paulis.real.astype(int)*effect[0] + paulis.imag.astype(int)*effect[1]).sum(axis=(1, 2))
    return 2*moments, int(np.trace(effect[0]))


def native_rows():
    rows = np.zeros((2, 12, 10, 31), dtype=object)
    for q, parity, route, tick, effect in clock.native_events():
        moments, trace = coordinates(effect)
        rows[q, route, tick-2, 15*parity:15*(parity+1)] = 2*moments
        rows[q, route, tick-2, -1] += (2 if parity else -2)*trace
    return rows


def exact_half_step(h, source, relay):
    """Integer Taylor numerator for the relay/source block of exp(-ih/4)."""
    assert h.dtype == np.int64 and not (h-h.T).nnz
    assert int(abs(h).sum(axis=1).max()) <= 43
    h = h.tocsr()
    degree = np.diff(h.indptr)
    panels = []
    for slot in range(int(max(degree))):
        addresses = np.flatnonzero(degree > slot)
        entries = h.indptr[addresses] + slot
        panels.append((addresses, h.indices[entries], h.data[entries, None]))
    power = np.zeros((h.shape[0], 4), dtype=object)
    power[4*source:4*source+4] = np.eye(4, dtype=int)
    denominator = 4**ORDER * factorial(ORDER)
    coefficient = denominator
    value = np.zeros((2, 4, 4), dtype=object)
    for k in range(ORDER+1):
        sign = 1 if k % 4 in (0, 3) else -1
        value[k % 2] += sign*coefficient*power[4*relay:4*relay+4]
        if k < ORDER:
            coefficient //= 4*(k+1)
            following = np.zeros_like(power)
            for addresses, neighbors, weights in panels:
                following[addresses] += weights*power[neighbors]
            power = following
    return value, denominator


def product(a, b):
    return np.stack((a[0]@b[0] - a[1]@b[1], a[1]@b[0] + a[0]@b[1]))


def precise_middle_contrasts():
    # For a real symmetric generator the reverse block is EXACTLY the transpose.
    # The half-step scalar phase cancels in each two-leg probability.
    station = full_clock.station.iso.pairing.formation
    board = station.federation(1)
    ports = station.region(board)
    common = np.zeros((2, 31), dtype=object)
    for qi, q in enumerate((-1, 1)):
        for parity in (0, 1):
            h = clock.native.generator(board, ports, q, bool(parity))
            k, d = exact_half_step(h, ports[1], ports[2])
            effects = []
            for leg in (k, k.transpose(0, 2, 1)):
                amp = product(leg.transpose(0, 2, 1), leg)
                effects.append(product(np.stack((amp[0].T, -amp[1].T)), amp))
            moments, trace = coordinates(effects[1] - effects[0])
            common[qi, 15*parity:15*(parity+1)] = moments
            common[qi, -1] += trace if parity else -trace
    # Each block is an orthogonal compression of a unitary. Its Taylor error is e.
    # Telescoping the four factors of one probability costs (1+e)^4-1.
    error = Fraction(43, 4)**(ORDER+1) / factorial(ORDER+1)
    coord_error = 4*SCALE*((1+error)**4-1) + Fraction(1, DEN)
    assert coord_error < Fraction(1, 10**28)
    return common*(SCALE*DEN)//(4*d**4), coord_error


def upper_sqrt(x):
    """Strict dyadic upper bound for the square root of a nonnegative fraction."""
    bits = 120
    return Fraction(isqrt((x.numerator << (2*bits)) // x.denominator)+1, 2**bits)


def upper(x):
    return float(np.nextafter(float(x), np.inf))


def lower(x):
    return float(np.nextafter(float(x), -np.inf))


def rank_certificate(rows, entry_error):
    """Verify a numerical right inverse by exact integer matrix multiplication."""
    a = rows if len(rows) <= rows.shape[1] else rows.T
    m, n = a.shape
    approx = np.asarray(a, dtype=float) / DEN
    proposal = np.linalg.pinv(approx)
    unit = 2**40
    right = np.array([[round(float(v)*unit) for v in row] for row in proposal], dtype=object)
    residual = np.eye(m, dtype=object)*(DEN*unit) - a@right
    rnorm = upper_sqrt(Fraction(int((residual*residual).sum()), (DEN*unit)**2))
    pnorm = upper_sqrt(Fraction(int((right*right).sum()), unit**2))
    defect = rnorm + upper_sqrt(Fraction(m*n))*entry_error*pnorm
    assert defect < 1, (m, n, float(defect))
    return dict(rank=m, inverse_defect_upper=upper(defect),
                smallest_nonzero_singular_value_lower=lower((1-defect)/pnorm))


def run():
    rows = native_rows()
    contrasts, narrow = precise_middle_contrasts()
    epsilon = 2*clock.timing_error()[0]
    assert narrow < epsilon
    evidence = []
    for edge in combinations(range(4), 2):
        forward = measure.ROUTES.index(edge)
        reverse = measure.ROUTES.index(edge[::-1])
        blocks = [rows[q, route] for q in (0, 1) for route in (forward, reverse)]
        # Indices 0,1,2,3 denote vacancy-forward/reverse, excess-forward/reverse.
        def stacked(indices):
            entries = []
            for i in indices:
                group = blocks[i].copy()
                if edge == MIDDLE and i % 2 and i-1 in indices:
                    group[0] = contrasts[i//2]
                entries.append(group)
            return np.vstack(entries)

        def cert_groups(indices):
            return rank_certificate(stacked(indices), epsilon)

        pairs = {}
        for name, negative, positive in (("forward", (0,), (2,)), ("reverse", (1,), (3,)),
                ("both_tagged", (0, 1), (2, 3)), ("opposite_assignments", (0, 3), (2, 1))):
            neg, pos, joint = (cert_groups(g) for g in (negative, positive, negative+positive))
            pairs[name] = dict(negative=neg, positive=pos, joint=joint,
                               common_dimension=neg["rank"] + pos["rank"] - joint["rank"])
        averages = []
        for left, right in ((0, 1), (2, 3)):
            total = blocks[left]+blocks[right]
            assert np.all(total % 2 == 0)
            averages.append(total//2)
        neg, pos, joint = (rank_certificate(a, epsilon)
                           for a in (averages[0], averages[1], np.vstack(averages)))
        pairs["orientations_pooled"] = dict(negative=neg, positive=pos, joint=joint,
                                           common_dimension=neg["rank"]+pos["rank"]-joint["rank"])
        direct_sums = {name: cert_groups(indices) for name, indices in
                       (("vacancy_forward_and_both_reverse", (0, 1, 3)),
                        ("excess_forward_and_both_reverse", (2, 3, 1)))}
        assert all(r["rank"] == 30 for r in direct_sums.values())
        evidence.append(dict(edge=list(edge), menus=pairs, direct_sums=direct_sums,
                             common_dimension_all_four_charge_assignments=0))
    contrast_floors = [lower(Fraction(int(max(abs(row))), SCALE*DEN)-narrow/SCALE)
                       for row in contrasts]
    contrast_norms = [upper((upper_sqrt(Fraction(int((row*row).sum()), DEN**2))
                            +upper_sqrt(Fraction(31))*narrow)/SCALE) for row in contrasts]
    return dict(ports=list(map(int, full_clock.station.iso.pairing.formation.region(
                    full_clock.station.iso.pairing.formation.federation(1)))),
                source_coordinates=31, ticks=list(range(2, 12)), ordinary_entry_error_upper=upper(epsilon),
                middle_precision=dict(Taylor_order=ORDER, tested_row_scale=str(SCALE),
                    rescaled_entry_error_upper=upper(narrow),
                    raw_coordinate_l_infinity_lower=contrast_floors, raw_row_l2_upper=contrast_norms),
                edges=evidence)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False)+"\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result, end="")
