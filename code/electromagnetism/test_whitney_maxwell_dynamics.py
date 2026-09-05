"""Adversarial controls for the same-metric stationary field history."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_whitney_maxwell_dynamics as verifier


def test_independent_stationary_volume_replay():
    result = verifier.verify(verifier.load())
    assert result["field_variations_per_history"] == 68
    assert result["events_per_history"] == 805
    assert result["exact_global_stiffness_bound"] == 24


def test_fresh_producer_passes_independent_quadrature_and_event_replay():
    import whitney_maxwell_dynamics as producer
    fresh = producer.build()
    result = verifier.verify(fresh)
    assert result["instrumented_slices_per_history"] == 3
    assert fresh["stability_certificate"] == verifier.load()["stability_certificate"]


def replace(path, value):
    def mutate(packet):
        row = packet
        for key in path[:-1]:
            row = row[key]
        row[path[-1]] = value
    return mutate


def event_mutation(op, mutate):
    def apply(packet):
        event = next(e for e in packet["executions"][0]["events"] if e["op"] == op)
        mutate(event)
    return apply


@pytest.mark.parametrize("mutation", [
    replace(["scope"], "PHYSICAL_CONTINUUM_CONFIRMED"),
    lambda p: p["pins"].pop("Lean/Screen/WhitneyMaxwellDynamics.lean"),
    replace(["pins", "code/electromagnetism/verify_cone_whitney_bridge.py"], "0"*64),
    replace(["numeric_policy", "atol"], 1.0),
    replace(["dynamics", "source"], "dynamical charged matter from the OPH source"),
    replace(["dynamics", "energy"], "raw velocity energy in arbitrary gauge"),
    replace(["stability_certificate", "mass_over_volume", 0, 0], ["1", "0"]),
    replace(["stability_certificate", "ldl", "24", "diagonal", 0], ["-1", "0"]),
    replace(["stability_certificate", "ldl", "24", "lower", 2, 1], ["0", "0"]),
    replace(["stability_certificate", "bindings", 0, "edge_signs", 0], -1),
    replace(["stability_certificate", "bindings", 0, "vertices", 0], False),
    lambda p: p["stability_certificate"]["bindings"].pop(),
    replace(["executions", 0, "instrumented_slices"], 65),
    replace(["executions", 0, "writable_slots"], 42),
    replace(["executions", 0, "gauge"], 0),
    replace(["executions", 0, "decoded", 2, 13], "100"),
    replace(["executions", 0, "decode_event_ids", 0], 0),
    replace(["executions", 0, "projection", "projected_E0", 0], "0"),
    lambda p: p["executions"][0]["metrics"]["ampere_hat"].__delitem__(slice(0, 12)),
    event_mutation("inputs", lambda e: e["writes"].__setitem__("h", "1")),
    event_mutation("gauss_project", lambda e: e["reads"].__delitem__("rho/0/0")),
    event_mutation("advance", lambda e: e["reads"]["d/1/13"].__setitem__("writer", 0)),
    event_mutation("advance", lambda e: e["reads"]["d/1/13"].__setitem__("value", "100")),
    event_mutation("advance", lambda e: e["parents"].append(999)),
    event_mutation("advance", lambda e: e["writes"].__setitem__("x/13", "0")),
    event_mutation("probe", lambda e: e["writes"].__setitem__("x/0", "100")),
    event_mutation("feedback", lambda e: e["reads"].__delitem__("clock")),
    event_mutation("decode", lambda e: e["writes"].__setitem__("d/0/13", "100")),
    event_mutation("public", lambda e: e["writes"].__setitem__("public/E/0/0", "100")),
    replace(["continuation", "scope"], "all time steps authenticated by serial instrument"),
    replace(["continuation", "A", 20, 0], 100.0),
    replace(["continuation", "rho_load", 20, 0], 1.0),
    replace(["continuation", "J_load", 20, 0], 1.0),
    replace(["continuation", "energy", 20], 100.0),
    replace(["continuation", "source_work", 0], 0.0),
    replace(["scalar_controls", "histories", 1, "A", 16], "0"),
    lambda p: p["nonclaims"].pop(),
])
def test_false_green_mutations_fail(mutation):
    original = verifier.load()
    packet = deepcopy(original)
    mutation(packet)
    assert json.dumps(packet, sort_keys=True) != json.dumps(original, sort_keys=True)
    with pytest.raises(ValueError):
        verifier.verify(packet)


def test_invisible_nonbinary_solver_write_rejected():
    packet = verifier.load()
    ex = packet["executions"][0]
    event = next(e for e in ex["events"] if e["op"] == "gauss_project")
    old = Q(event["writes"]["projection_z/0"])
    forged = old+Q(1, 10**1000)
    assert float(old) == float(forged) and forged != old
    event["writes"]["projection_z/0"] = str(forged)
    ex["projection"]["z"][0] = str(forged)
    with pytest.raises(ValueError, match="exact float64 output encoding"):
        verifier.verify(packet)


@pytest.mark.parametrize("raw", ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_ambiguous_json_rejected(tmp_path, raw):
    path = tmp_path / "bad.json"
    path.write_text(raw, encoding="utf-8")
    with pytest.raises(ValueError):
        verifier.load(path)


def test_fresh_transitive_custody_after_prior_success(monkeypatch):
    packet = verifier.load()
    verifier.verify(packet)
    original = Path.read_bytes
    upstream = verifier.ROOT / "Lean/Screen/SerialMaxwellReadout.lean"
    def tampered(path):
        content = original(path)
        return content+b"\n-- altered transitive provider\n" if path == upstream else content
    monkeypatch.setattr(Path, "read_bytes", tampered)
    # Neither the new receipt nor the unchanged immediate parent's hash is
    # enough to authenticate all of the immediate parent's current providers.
    with pytest.raises(ValueError, match="source pin"):
        verifier.verify(packet)


def test_exact_variational_map_and_modified_hamiltonian():
    q, next_q, p, h, lam = sp.symbols("q next_q p h lam", positive=True)
    discrete_l = (next_q-q)**2/(2*h)-h*lam*(q*q+q*next_q+next_q*next_q)/6
    solution = sp.solve(sp.Eq(p, -sp.diff(discrete_l, q)), next_q)[0]
    next_p = sp.diff(discrete_l, next_q).subs(next_q, solution)
    f = 1+h*h*lam/6
    a, b, cc = (1-h*h*lam/3)/f, h/f, -h*lam*(1-h*h*lam/12)/f
    assert sp.simplify(solution-a*q-b*p) == 0
    assert sp.simplify(next_p-cc*q-a*p) == 0
    assert sp.simplify(a*a-b*cc) == 1
    # Consistent canonical initial data give second-order time convergence.
    errors = []
    exact = np.array([[np.cos(1), np.sin(1)], [-np.sin(1), np.cos(1)]])
    for count in (16, 32, 64):
        step = 1/count
        z = step*step
        factor = 1+z/6
        matrix = np.array([[1-z/3, step], [-step*(1-z/12), 1-z/3]])/factor
        errors.append(np.linalg.norm(np.linalg.matrix_power(matrix, count)-exact))
    assert 3.99 < errors[0]/errors[1] < 4.01
    assert 3.99 < errors[1]/errors[2] < 4.01


def test_stable_initial_value_problem_can_have_singular_endpoint_problem():
    h, lam = Q(1, 2), Q(12)  # z=3, strictly inside the IVP stability window.
    f, middle = 1+h*h*lam/6, 2-2*h*h*lam/3
    assert f > 0 and middle == 0 and h*h*lam < 12
    # With zero endpoints, every midpoint solves J=0; none solves J=1.
    for q1 in (Q(-3), Q(0), Q(7, 2)):
        assert f*0-middle*q1+f*0 == h*h*0
        assert f*0-middle*q1+f*0 != h*h*1


def test_explicit_utf8_file_reads(monkeypatch):
    original = Path.read_text
    def windows_default(path, *args, **kwargs):
        if not args and "encoding" not in kwargs:
            kwargs["encoding"] = "cp1252"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", windows_default)
    assert verifier.verify(verifier.load())["field_variations_per_history"] == 68
