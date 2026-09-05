"""False-green controls for serial feedback, causal consumption and field joins."""
import copy
from fractions import Fraction as Q
from pathlib import Path
import sys

import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import serial_maxwell_readout as producer
import verify_serial_maxwell_readout as verifier


@pytest.fixture(scope="module")
def receipt():
    return verifier.load()


def test_pinned_receipt_and_producer_replay(receipt):
    assert producer.build() == receipt
    assert verifier.verify(receipt) == {
        "action": sp.Rational(-304804229, 21123216), "events": 585,
        "cycles": 180, "path_variations": 8, "field_variations": 54}


def test_cp1252_default_cannot_change_source_decoding(receipt, monkeypatch):
    original = Path.read_text
    def windows_default(path, encoding=None, errors=None):
        return original(path, encoding=encoding or "cp1252", errors=errors)
    monkeypatch.setattr(Path, "read_text", windows_default)
    assert producer.build() == receipt
    assert verifier.verify(verifier.load())["events"] == 585


@pytest.mark.parametrize("mutation", ["missing_feedback", "wrong_feedback", "wrong_writer",
    "forged_parent", "synthetic_read", "wrong_response", "mixed_frame", "omitted_probe",
    "wrong_current", "external_reset", "clock_substitution", "direct_reference_advance",
    "promoted_scope", "omitted_gauge", "wrong_pin", "omitted_assumption",
    "boolean_writer", "boolean_parent"])
def test_corrupt_execution_is_rejected(receipt, mutation):
    value = copy.deepcopy(receipt)
    events = value["executions"][0]["events"]
    def event(op, frame=None):
        return next(e for e in events if e["op"] == op and (frame is None or e["args"][0] == frame))
    if mutation == "missing_feedback":
        events.remove(event("feedback"))
    elif mutation == "wrong_feedback":
        event("feedback")["writes"]["x/12"] = "1"
    elif mutation == "wrong_writer":
        event("feedback")["reads"]["r/0/0/0"]["writer"] = 0
    elif mutation == "forged_parent":
        event("feedback")["parents"].append(0)
    elif mutation == "synthetic_read":
        # A genuine value and writer still cannot fabricate opcode support.
        e = event("feedback")
        e["reads"]["q/0"] = {"value": "1", "writer": 0}
        e["parents"] = sorted(set(e["parents"]+[0]))
    elif mutation == "wrong_response":
        event("response")["writes"]["r/0/0/0"] = "1"
    elif mutation == "mixed_frame":
        e = event("feedback", 1)
        e["reads"]["b/0/0"] = e["reads"].pop("b/1/0")
    elif mutation == "omitted_probe":
        events.remove(event("probe"))
    elif mutation == "wrong_current":
        event("sources")["writes"]["J/0/0"] = "0"
    elif mutation == "external_reset":
        event("feedback")["op"] = "restore_checkpoint"
    elif mutation == "clock_substitution":
        event("public")["writes"]["completed_probe_cycles"] = "585"
    elif mutation == "direct_reference_advance":
        event("advance")["reads"] = {}
        event("advance")["parents"] = []
    elif mutation == "promoted_scope":
        value["scope"] = "PHYSICAL_MAXWELL_DERIVED"
    elif mutation == "omitted_gauge":
        value["executions"].pop()
    elif mutation == "wrong_pin":
        value["pins"]["Lean/Screen/SerialMaxwellReadout.lean"] = "0"*64
    elif mutation == "boolean_writer":
        event("sources")["reads"]["h"]["writer"] = False
    elif mutation == "boolean_parent":
        event("sources")["parents"] = [False]
    else:
        value["assumptions"].pop()
    with pytest.raises(ValueError):
        verifier.verify(value)


def test_arbitrary_overlapping_word_and_omitted_feedback_counterexample():
    initial = [Q(i*i-3*i+1, 2*i+1) for i in range(42)]
    state, baseline = initial[:], initial[:12]
    edges = verifier.carrier()[0]
    word = [(u, 12+e) for e, pair in enumerate(edges) for u in pair]
    word = word+word[::-1]+word[::3]
    for u, m in word:
        response = (state[u]+state[m])/2
        state[u] = state[m] = response
        state[u], state[m] = baseline[u], 2*response-baseline[u]
        assert state == initial
    state = initial[:]
    for u, m in word[:2]:
        state[u] = state[m] = (state[u]+state[m])/2
    u, m = word[1]
    assert 2*state[u]-baseline[u] != initial[m]


def test_local_record_error_law_does_not_hide_accumulation():
    b, m, eb, er = Q(2, 3), Q(-5, 7), Q(1, 101), Q(-1, 89)
    decoded = 2*((b+m)/2+er)-(b+eb)
    assert decoded-m == 2*er-eb
    assert abs(decoded-m) <= 2*abs(er)+abs(eb)


def test_source_discontinuity_passes_ampere_but_breaks_gauss():
    edges, d, c = verifier.carrier()
    a = verifier.reference()
    h = sp.Rational(1, 2)
    rho = [d.T*(-(a[n+1]-a[n])/h) for n in range(2)]
    current = [sp.zeros(30, 1), sp.zeros(30, 1)]
    a[2] = 2*a[1]-a[0]-h*h*c.T*c*a[1]
    e, b, _ = verifier.fields_and_action(d, c, h, a, [sp.zeros(12, 1)]*2, rho, current)
    assert e[1]-e[0] == h*c.T*b[1]
    assert d.T*e[1] != rho[1]
    with pytest.raises(ValueError, match="Gauss"):
        verifier.physics(d, c, edges, h, a, [sp.zeros(12, 1)]*2, rho, current)


def test_endpoint_gauge_restriction_is_material():
    _, d, c = verifier.carrier()
    a = verifier.reference()
    h = sp.Rational(1, 2)
    phi = [sp.zeros(12, 1)]*2
    rho = [d.T*(-(a[n+1]-a[n])/h) for n in range(2)]
    j = [sp.Matrix([-2 if e == 0 else 2 if e == 29 else 0 for e in range(30)]), sp.zeros(30, 1)]
    e, b, action = verifier.fields_and_action(d, c, h, a, phi, rho, j)
    chi0 = sp.eye(12)[:, 0]
    aa, pp = [a[0]+d*chi0, a[1], a[2]], [phi[0]+chi0/h, phi[1]]
    ee, bb, changed = verifier.fields_and_action(d, c, h, aa, pp, rho, j)
    assert (ee, bb) == (e, b)
    assert changed-action == rho[0].dot(chi0) == 1


@pytest.mark.parametrize("raw", ['{"id":0,"id":1}', '{"x":NaN}'])
def test_noncanonical_json_fails(tmp_path, raw):
    path = tmp_path/"bad.json"
    path.write_text(raw)
    with pytest.raises(ValueError):
        verifier.load(path)
