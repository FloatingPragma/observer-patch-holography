"""Independent full-catalog replay and hostile evidence controls."""

import copy
import fnmatch
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from m1_vacuum_fidelity import check, model, verify


@pytest.fixture(scope="module")
def packet():
    return verify.strict_load(verify.HERE/"receipt.json")


def save(tmp_path, value):
    path = tmp_path/"candidate.json"
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_committed_receipt_and_independent_producer(packet):
    verify.verify()
    check.same(packet["evidence"], model.candidate(), "producer differs")


@pytest.mark.parametrize("scale", [0, 2])
def test_stability_certificate_is_bound_to_executed_polynomials(monkeypatch, scale):
    original = check.polynomial_step
    # Both replacements remain symmetric and symplectic; those tests alone
    # cannot justify applying the intended step's stability coefficients.
    monkeypatch.setattr(check, "polynomial_step", lambda weight: original(scale*weight))
    with pytest.raises(ValueError, match="stability coefficients"):
        check.exact_algebra.__wrapped__()


@pytest.mark.parametrize("order", [2, 4])
@pytest.mark.parametrize("z_text", ["0.001", "0.01", "0.05", "0.1"])
def test_explicit_coherent_band_bound(order, z_text):
    with mp.workdps(65):
        z = mp.mpf(z_text)
        step = model.mode_step(z, order)
        for j in (1, 11, 100, 1001):
            exact = mp.matrix([[mp.cos(j*z), mp.sin(j*z)], [-mp.sin(j*z), mp.cos(j*z)]])
            delta = step**j-exact
            square = delta.T*delta
            largest = (square[0, 0]+square[1, 1]
                       +mp.sqrt((square[0, 0]-square[1, 1])**2+4*square[0, 1]**2))/2
            assert mp.sqrt(largest) <= z**order*(1+j*z)


@pytest.mark.parametrize("path,value", [
    (("algebra", "4", "b", 5, 0), "0"),
    (("algebra", "4", "c", 7, 1), "1"),
    (("graphs", "axis3", "sites"), True),
    (("graphs", "axis4", "degree"), 0),
    (("graphs", "ball6", "moments", "5"), "0"),
    (("graphs", "ball6", "spectrum", 0, 1), 0),
    (("graphs", "two_scale6", "vectors", 0, 0), 99),
    (("graphs", "axis4", "quantum", "4:1/256", "observations", 9, 1), "0.0"),
    (("graphs", "axis3", "quantum", "2:1/512", "mean", 2), "NaN"),
    (("graphs", "axis3", "quantum", "4:1/512", "times", 7), 7),
    (("programs", "4", "nonlocal_reads"), 0),
    (("programs", "4", "two_mode_phase_gates"), 0),
    (("programs", "4", "signed_word", 3, 1), "1.0"),
    (("programs", "2", "trace_sha256"), "0"*64),
    (("programs", "4", "inverse_writes"), 0),
    (("programs", "4", "final_p", 3), "0.0"),
    (("first_sparse_level", "second_moment"), 1),
    (("first_sparse_level", "stability_square_bound"), "0"),
    (("exponents", "sparse", "4", "energy_tick_power"), "33/64"),
    (("exponents", "explicit_repair", "second_order_energy_power"), "67/32"),
    (("extreme_signals", "4:3", "outermost_signal"), "0.0"),
    (("extreme_signals", "4:3", "max_p_distance"), 3),
    (("extreme_signals", "2:2", "p", 10), "1.0"),
])
def test_changed_evidence_is_rejected(packet, path, value):
    changed = copy.deepcopy(packet["evidence"])
    target = changed
    for component in path[:-1]:
        target = target[component]
    target[path[-1]] = value
    with pytest.raises(ValueError, match="independent reconstruction"):
        check.verify_evidence(changed)


@pytest.mark.parametrize("section", ["algebra", "graphs", "programs", "exponents", "extreme_signals"])
def test_omitted_catalog_rejected(packet, section):
    changed = copy.deepcopy(packet["evidence"])
    changed[section].pop(next(iter(changed[section])))
    with pytest.raises(ValueError):
        check.verify_evidence(changed)


def forged_stationary_vacuum(packet):
    changed = copy.deepcopy(packet)
    for graph in changed["evidence"]["graphs"].values():
        for case in graph["quantum"].values():
            case["observations"] = [["0.0"]*3 for _ in case["times"]]
            case["mean"] = ["0.0"]*3
    return changed


def test_coordinated_modified_vacuum_substitution_rejected(packet):
    with pytest.raises(ValueError):
        check.verify_evidence(forged_stationary_vacuum(packet)["evidence"])


def test_producer_cannot_filter_observation_times(packet, monkeypatch):
    monkeypatch.setattr(model, "TIMES", (1, 2, 3))
    with pytest.raises(ValueError, match="32 declared observation times"):
        model.quantum_case(model.spectrum(3, model.stencil("axis")), model.TICKS[0], 2)
    changed = copy.deepcopy(packet["evidence"])
    case = changed["graphs"]["axis3"]["quantum"]["2:1/256"]
    case["times"] = case["times"][:3]
    case["observations"] = case["observations"][:3]
    with pytest.raises(ValueError):
        check.verify_evidence(changed)


def test_inverse_sign_change_is_not_hidden_by_recomputed_outputs(packet, monkeypatch):
    original = model.gate_word
    monkeypatch.setattr(model, "gate_word", lambda order, steps:
                        [(k, abs(c)) for k, c in original(order, steps)])
    changed = copy.deepcopy(packet["evidence"])
    changed["programs"]["4"] = model.program(4)
    with pytest.raises(ValueError):
        check.verify_evidence(changed)


def test_stale_writers_rejected_even_with_unchanged_outputs(packet, monkeypatch):
    original = model.digest

    def stale_digest(records):
        changed = copy.deepcopy(records)
        for record in changed:
            for read in record[3]:
                read[1] = 0
        return original(changed)

    monkeypatch.setattr(model, "digest", stale_digest)
    changed = copy.deepcopy(packet["evidence"])
    changed["programs"]["4"] = model.program(4)
    honest = packet["evidence"]["programs"]["4"]
    forged = changed["programs"]["4"]
    assert forged["final_q"] == honest["final_q"]
    assert forged["final_p"] == honest["final_p"]
    assert forged["trace_sha256"] != honest["trace_sha256"]
    with pytest.raises(ValueError, match="independent reconstruction"):
        check.verify_evidence(changed)


@pytest.mark.parametrize("bad", ["{}", "null", "[]", '{"x":NaN}', '{"x":Infinity}', '{"x":1,"x":2}'])
def test_malformed_receipts_fail(tmp_path, bad):
    path = tmp_path/"bad.json"
    path.write_text(bad, encoding="utf-8")
    with pytest.raises((ValueError, TypeError)):
        verify.verify(path)


def test_stale_pin_and_extra_envelope_are_rejected(packet, tmp_path):
    changed = copy.deepcopy(packet)
    changed["sources"][next(iter(changed["sources"]))] = "0"*64
    with pytest.raises(ValueError, match="source pins"):
        verify.verify(save(tmp_path, changed))
    changed = copy.deepcopy(packet)
    changed["approved"] = True
    with pytest.raises(ValueError, match="envelope"):
        verify.verify(save(tmp_path, changed))


def test_parent_rows_must_be_unique(tmp_path):
    registry = verify.strict_load(verify.ROOT/"claims/claim_registry.yaml")
    row = next(x for x in registry["claims"] if x["claim_id"] == verify.PARENTS[0])
    registry["claims"].append(copy.deepcopy(row))
    with pytest.raises(ValueError, match="exactly once"):
        verify.parent_claims(save(tmp_path, registry))


@pytest.mark.parametrize("forgery", ["vacuum", "spectrum", "empty", "clock"])
def test_real_optimized_interpreter_rejects_forgery(packet, tmp_path, forgery):
    changed = forged_stationary_vacuum(packet) if forgery == "vacuum" else copy.deepcopy(packet)
    if forgery == "spectrum":
        changed["evidence"]["graphs"]["ball6"]["spectrum"].pop()
    if forgery == "empty":
        changed = {}
    if forgery == "clock":
        for row in changed["evidence"]["extreme_signals"].values():
            row["p"] = ["0.0"]*row["sites"]
            row["outermost_signal"] = "0.0"
            row["max_p_distance"] = 0
    result = subprocess.run([sys.executable, "-O", str(verify.HERE/"verify.py"), str(save(tmp_path, changed))],
                            capture_output=True, text=True, cwd=verify.ROOT)
    assert result.returncode != 0
    assert "ValueError" in result.stderr


def test_verifier_runs_with_producer_imports_blocked():
    script = """
import importlib.abc, runpy, sys
class RejectProducer(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in ('m1_vacuum_fidelity.model', 'm1_vacuum_fidelity.build'):
            raise RuntimeError('producer import forbidden')
sys.meta_path.insert(0, RejectProducer())
sys.argv = [sys.argv[1]]
runpy.run_path(sys.argv[0], run_name='__main__')
"""
    result = subprocess.run([sys.executable, "-c", script, str(verify.HERE/"verify.py")],
                            capture_output=True, text=True, cwd=verify.ROOT)
    assert result.returncode == 0, result.stderr


def test_every_pinned_source_triggers_replay():
    workflow = yaml.load((verify.ROOT/".github/workflows/m1-vacuum-fidelity.yml").read_text(encoding="utf-8"),
                         Loader=yaml.BaseLoader)
    for event in ("push", "pull_request"):
        patterns = workflow["on"][event]["paths"]
        for path in verify.SOURCE_FILES+["claims/claim_registry.yaml"]:
            assert any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns), (event, path)
    lean = (verify.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    assert '"Geometry.M1VacuumFidelityAxiomAudit"' in lean
