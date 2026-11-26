import json
import math

import numpy as np
import sympy as sp
from mpmath import iv

from ...full_algebra.global_survey import positive_definite
from ...full_geometry.measure import leg_kernels, pairing, ranging
from ..measure import HERE as CLOCKS, blank_clock, density_center
from . import measure as m


def test_reference_state_lower_bound_and_rest_gap():
    raw = json.loads((CLOCKS.parent / "full_algebra/results.json").read_text())
    blocks = [m.native.prior.sym_matrix(a) for a in np.array(raw["x_real"]) + 1j*np.array(raw["x_imag"])]
    blocks = [(a+a.H)/2 for a in blocks]
    trace = sum(a.trace() for a in blocks)
    for a in blocks:
        assert positive_definite(a/trace - sp.eye(4)/50)
    for charge in (-1, 1):
        np.testing.assert_array_equal(2*blank_clock(charge, 8), [[17, -1], [-1, 17]])


def test_held_phase_controls_and_determinant_reversal():
    record = json.loads((m.HERE / "results.json").read_text())
    comparisons = m.comparisons(record["rows"])
    assert comparisons["-1"][0][0] > 8.58
    assert comparisons["1"][0][0] > comparisons["1"][1][1]
    for row in record["rows"]:
        r = m.matrix_record(row["counts"])
        np.testing.assert_array_equal(r["det"], row["measured"]["det"])
        for field in ("erased_label_difference", "held_observed_minus_forecast"):
            assert all(lo < 0 < hi for lo, hi in r[field])
        if row["mode"] == "rest":
            assert r["det"][0] < 1 < r["det"][1]
        elif row["time"] == 0:
            assert all(counts[2] == m.N for bank in row["counts"].values() for counts in bank)
            assert r["held_phase"] is None
        elif row["time"] == 2:
            assert r["det"][1] < 0
        else:
            assert r["det"][0] > 0
        if row["time"]:
            lo, hi = r["held_phase"]
            assert lo < math.pi/4 < hi


def test_population_signs_and_independent_flight():
    states = density_center()[0]
    board = pairing.formation.federation(1)
    source = pairing.formation.region(board)
    flows = m.native_flows()
    for charge in (-1, 1):
        independent = np.zeros((3, 2, 2, 2, 2), complex)
        for parity, state in enumerate(states):
            h = ranging.generator(board, source, charge, 8, bool(parity))
            u = leg_kernels(h, m.PORTS)
            for i, t in enumerate(m.TIMES):
                independent[i] += m.outputs(u[2*t, 2:, :2], state)
        conditional = []
        for i, t in enumerate(m.TIMES):
            flow = flows[charge, t]
            np.testing.assert_allclose(independent[i], flow, rtol=0, atol=1e-14)
            data = {}
            normalized = {}
            for label in ("Z+", "Z-", "Y+", "Y-"):
                data[label] = []
                normalized[label] = []
                for plus, minus in m.score_probabilities(flow, label):
                    e = 1.00002e-5*(plus+minus)+1e-12
                    signed = iv.mpf([np.nextafter(plus-minus-e, -np.inf),
                                     np.nextafter(plus-minus+e, np.inf)])
                    total = iv.mpf([np.nextafter(plus+minus-e, -np.inf),
                                    np.nextafter(plus+minus+e, np.inf)])
                    data[label].append(signed)
                    normalized[label].append(signed / total)
            a, b = ((data[k+"+"][0]-data[k+"-"][0])/2 for k in ("Z", "Y"))
            c, d = ((data[k+"+"][1]-data[k+"-"][1])/2 for k in ("Z", "Y"))
            det = a*d-b*c
            assert det.b < 0 if t == 2 else det.a > 0
            conditional.append([(normalized[k+"+"][j]-normalized[k+"-"][j])/2
                                for j in (0, 1) for k in ("Z", "Y")])
        lower, upper = m.rotation_invariant(conditional[0], conditional[1])
        lo, hi = (9.2758, 9.4121) if charge == -1 else (4.6654, 4.7131)
        assert lo < lower <= upper < hi
