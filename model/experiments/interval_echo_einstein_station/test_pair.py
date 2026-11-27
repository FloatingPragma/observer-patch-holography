import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import cap
from .pair import measure


def test_one_scalar_normalization_cannot_repair_the_subcap_mismatch():
    caps = json.loads(Path(__file__).with_name("cap.json").read_text())["nested_caps"]["regions"]
    errors = {}
    for region in caps:
        m = np.array([t["modular_energy"] for t in region["tangents"]])
        s = np.array([t["seam_energy"] for t in region["tangents"]])
        gain = (s @ m) / (s @ s)
        errors[region["pole"], len(region["ports"])] = np.linalg.norm(m - gain*s)/np.linalg.norm(m)
    assert errors[1, 12] < 0.000023
    assert errors[1, 6] > 0.378 and errors[1, 11] > 0.344 and errors[0, 6] > 0.509


def test_restriction_keeps_a_load_on_an_enrolled_port_for_any_port_iterable():
    state = np.zeros((len(cap.WORDS), len(cap.WORDS)))
    i = cap.WORD_INDEX[(0,)]
    state[i, i] = 1
    for ports in [(0, 8), iter((0, 8)), (0, 0, 8)]:
        region = cap.CapRestriction(cap.WORDS, ports)
        assert region.ports == (0, 8)
        np.testing.assert_array_equal(region.reduce(state), np.diag([0, 1, 0, 0]))


def test_entire_single_seam_observable_span_misses_the_six_port_modular_response():
    receipt = json.loads(cap.HERE.joinpath("cap.json").read_text())
    witness = receipt["nested_caps"]["off_seam_witness"]
    a, b = witness["nonadjacent_ports"]
    edges = cap.cell().seams
    assert all({a, b} != set(edge) for edge in edges)

    # Difference of the one-load states (|a> ± |b>)/sqrt(2).
    tangent = np.zeros((len(cap.WORDS), len(cap.WORDS)))
    i, j = (cap.WORD_INDEX[(p,)] for p in (a, b))
    tangent[i, j] = tangent[j, i] = 1
    for edge in edges:
        np.testing.assert_array_equal(cap.CapRestriction(cap.WORDS, edge).reduce(tangent), 0)

    weights, _ = cap.cap_geometry(receipt["pole"])
    h, _, _ = cap.native_energy([cap.seam_indices(edge) for edge in edges], weights)
    gibbs = expm(-receipt["exchange_beta"] * h)
    region = cap.CapRestriction(cap.WORDS, receipt["nested_caps"]["regions"][0]["ports"])
    modular, _ = region.charge(gibbs / np.trace(gibbs))
    signed = np.trace(modular @ region.reduce(tangent))
    np.testing.assert_allclose(signed, witness["exact_gibbs_limit_nats"], atol=1e-11, rtol=0)
    assert signed > 0.00459


def test_modular_energy_comparison_retains_quiet_and_null_support():
    report, readings = measure()
    assert sum(report['cutoffs'][1]['ranks']) == 17
    for fit in report['cutoffs']:
        assert fit['max_outside_support_frobenius'] < 1e-11
        assert fit['port9_prediction'] < -0.01 and fit['energy_contrasts'][1] > 0.09
        assert fit['matches_held'] is False
    budget = report['reference_survey']['shots_per_setting']
    for name, counts in readings.items():
        if name == 'diagonal':
            assert counts.sum() == budget
        else:
            np.testing.assert_equal(counts.sum(axis=1), budget)
