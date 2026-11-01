"""Replay, scientific-boundary and false-green controls for the captured bus."""
import ast
from copy import deepcopy
from fractions import Fraction as F
from fnmatch import fnmatchcase
import inspect
import re
import shlex
import shutil
import subprocess
import sys

import pytest
import yaml

from . import build, codec, verify


@pytest.fixture(scope="module")
def controls():
    return codec.load(codec.HERE/"controls.json")


def test_independent_replay_and_retained_receipt(controls):
    codec.equal(verify.verify(controls), codec.load(codec.HERE/"receipt.json"), "receipt")


def test_strict_cross_platform_artifact_bytes():
    assert codec.canonical(build.build()) == (codec.HERE/"controls.json").read_bytes()
    assert codec.canonical(verify.verify()) == (codec.HERE/"receipt.json").read_bytes()


def test_verifier_has_no_producer_import():
    tree = ast.parse((codec.HERE/"verify.py").read_text(encoding="utf-8"))
    imports = [node for node in ast.walk(tree) if isinstance(node, (ast.Import,ast.ImportFrom))]
    assert not any("build" in (getattr(node,"module","") or "") or
                   any("build" in name.name for name in node.names) for node in imports)


@pytest.mark.parametrize("mutation", [
    "precision", "initial_error", "per_mean_error", "readout_error", "pin", "port_alias",
    "remote_port", "requests", "round_count", "payload", "prepared_bus", "delete_mean",
    "extra_mean", "changed_edge", "changed_input", "stale_writer", "assigned_reset",
    "wrong_round", "noninteger", "boolean", "decoded", "read_version", "read_use",
    "read_position", "read_writer", "read_units", "final_load", "final_writer",
    "register_cost", "preparation_cost", "mean_cost", "read_cost", "sample_cost",
    "write_cost", "state_bits", "workspace_bits", "program_bits", "port_bits",
    "counter_bits", "version_bits", "request_bits", "receiver_workspace", "read_schedule",
    "receiver_address", "cross_carrier_cost", "missing_history", "extra_field",
])
def test_rejects_corrupt_or_optimistic_artifacts(controls, mutation):
    data = deepcopy(controls)
    case = data["executions"][0]
    if mutation == "precision": data["precision"] = 19
    elif mutation == "initial_error": data["assumed_error_bounds"]["initial"] = "0"
    elif mutation == "per_mean_error": data["assumed_error_bounds"]["extra_per_mean"] = "0"
    elif mutation == "readout_error": data["assumed_error_bounds"]["terminal_per_rail"] = "0"
    elif mutation == "pin": data["source_sha256"][codec.PIN_PATHS[0]] = "0"*64
    elif mutation == "port_alias": data["global_ports"][6] = data["global_ports"][4]
    elif mutation == "remote_port": data["global_ports"][10] = 44
    elif mutation == "requests": data["requests"] = [0,1,1]
    elif mutation == "round_count": data["cleanup_sweeps"] = 79
    elif mutation == "payload": case["payloads"][0] = "-1/2"
    elif mutation == "prepared_bus": case["initial_units"][4] += 1
    elif mutation == "delete_mean": case["tape"].pop(8)
    elif mutation == "extra_mean": case["tape"].append(case["tape"][-1])
    elif mutation == "changed_edge": case["tape"][2][0] = 0
    elif mutation == "changed_input": case["tape"][0][2] += 1
    elif mutation == "stale_writer": case["tape"][9][4] = -5
    elif mutation == "assigned_reset": case["tape"][9][-1] = 2*2**20
    elif mutation == "wrong_round": case["tape"][30][-1] += 1
    elif mutation == "noninteger": case["tape"][0][0] = 0.0
    elif mutation == "boolean": case["tape"][0][0] = False
    elif mutation == "decoded": case["reads"][2]["decoded"] = "3/2"
    elif mutation == "read_version": case["reads"][2]["record"] = 1
    elif mutation == "read_use": case["reads"][2]["use"] = 1
    elif mutation == "read_position": case["reads"][0]["after_mean"] = 6
    elif mutation == "read_writer": case["reads"][1]["writers"][0] = 7
    elif mutation == "read_units": case["reads"][0]["local_units"][0] += 1
    elif mutation == "final_load": case["final_units"][0] += 1
    elif mutation == "final_writer": case["final_writers"][0] = -1
    elif mutation == "missing_history": data["executions"].pop()
    elif mutation == "extra_field": case["remote_payload_side_channel"] = "-3/2"
    else:
        keys = {"register_cost":"scalar_registers", "preparation_cost":"preparation_writes",
                "mean_cost":"means", "read_cost":"mean_reads", "sample_cost":"receiver_scalar_samples",
                "write_cost":"writes", "state_bits":"scalar_storage_bits_bound",
                "workspace_bits":"mean_adder_bits_bound", "program_bits":"expanded_local_schedule_bits",
                "port_bits":"global_port_map_bits", "counter_bits":"program_counter_bits",
                "version_bits":"version_use_counter_bits", "request_bits":"request_record_bits",
                "receiver_workspace":"receiver_arithmetic_bits_bound",
                "read_schedule":"read_schedule_bits", "receiver_address":"receiver_address_bits",
                "cross_carrier_cost":"cross_carrier_means"}
        case["resources"][keys[mutation]] -= 1
    with pytest.raises(ValueError):
        verify.verify(data)


def test_deleted_physical_seam_is_rejected(monkeypatch, controls):
    edges, count = verify.support_edges()
    assert (23,45) in edges
    edges.remove((23,45))
    monkeypatch.setattr(verify,"support_edges",lambda: (edges,count))
    with pytest.raises(ValueError, match="non-native seam"):
        verify.layout(controls["global_ports"])


def test_uniform_coefficient_error_is_not_only_a_sample_check():
    rows, _ = verify.ideal_maps()
    assert rows[2][1] != 0
    corrupt = (rows[2][0],rows[2][1]+F(1,4))
    with pytest.raises(ValueError,match="universal ideal transfer error"):
        verify.read_bounds(2,corrupt)


def test_old_record_after_other_record_with_all_interventions(controls):
    observations = {}
    for case in controls["executions"]:
        p0,p1 = case["payloads"]
        assert [r["decoded"] for r in case["reads"]] == [p0,p1,p0]
        observations[(p0,p1)] = [r["decoded"] for r in case["reads"]]
    levels = ("-3/2","-1/2","1/2","3/2")
    for first in levels:
        assert {observations[(first,other)][2] for other in levels} == {first}
    receipt = codec.load(codec.HERE/"receipt.json")
    for case in receipt["histories"]:
        assert set(range(4)) <= set(case["reads"][2]["initial_ancestor_ports"])


def test_rounded_resource_bounds_and_identity_steps_are_retained(controls):
    for case in controls["executions"]:
        assert all(0 <= row[6] < 2**23 and row[2]+row[3] < 2**24 for row in case["tape"])
        assert any(row[2] == row[3] == row[6] for row in case["tape"])
        mass_changes = [2*row[6]-row[2]-row[3] for row in case["tape"]]
        assert max(map(abs,mass_changes)) == 1  # Rounding is not exact conservation.
    altered = deepcopy(controls)
    tape = altered["executions"][0]["tape"]
    index = next(i for i,row in enumerate(tape) if row[2] == row[3] == row[6])
    tape.pop(index)
    with pytest.raises(ValueError,match="missing/extra mean"):
        verify.verify(altered)


@pytest.mark.parametrize("u", range(-10,11))
def test_integer_rounding_bound_and_ties_match_theorem(u):
    for v in range(-10,11):
        result = build.round_mean(u,v)
        assert result == (u+v)//2 + int((u+v)%4 == 3)
        assert abs(2*result-u-v) <= 1
        assert min(u,v) <= result <= max(u,v)
        if (u+v)%2:
            assert result%2 == 0


def test_low_precision_collapses_the_receiver_record():
    case = build.execute((F(1,2),F(-3,2)),precision=2)
    first = case["reads"][0]
    assert first["local_units"][0] == first["local_units"][1]
    assert first["decoded"] is None
    assert verify.decode_receiver(*first["local_units"],1,4) is None
    assert F(8,2*4) > F(1,32)


def test_receiver_uses_local_samples_and_declared_counter_only():
    assert list(inspect.signature(build.decode_local).parameters) == ["positive","negative","use","grid"]
    assert list(inspect.signature(verify.decode_receiver).parameters) == ["positive","negative","use","grid"]
    rows,_ = verify.ideal_maps()
    for cycle,row in enumerate(rows):
        bound = verify.read_bounds(cycle,row)
        radius = F(bound["total_error_bound"])
        use = bound["use"]
        for a in (F(-3,2),F(-1,2),F(1,2),F(3,2)):
            for sign in (-1,1):
                contrast = a/2**(use+3)+sign*radius
                p,m = (2+contrast)*2**20,(2-contrast)*2**20
                assert verify.decode_receiver(p,m,use,2**20) == str(a)


def test_exact_midpoint_is_ambiguous():
    for midpoint in (-1,0,1):
        contrast = F(midpoint,16)
        p,m = (2+contrast)*2**20,(2-contrast)*2**20
        assert build.decode_local(p,m,1,2**20) is None
        assert verify.decode_receiver(p,m,1,2**20) is None


def test_strict_json_rejects_duplicate_keys(tmp_path):
    p=tmp_path/"duplicate.json"
    p.write_text('{"x":1,"x":2}',encoding="utf-8")
    with pytest.raises(ValueError,match="duplicate"):
        codec.load(p)


def test_every_public_theorem_is_transitively_audited():
    proof=(codec.ROOT/"Lean/Geometry/SourceReusableBus.lean").read_text(encoding="utf-8")
    audit=(codec.ROOT/"Lean/Geometry/SourceReusableBusAxiomAudit.lean").read_text(encoding="utf-8")
    names=re.findall(r"^theorem (\w+)",proof,re.M)
    assert len(names)==25
    covered=re.findall(r"^audit_reusable_bus_axioms OPH.SourceReusableBus\.(\w+)$",audit,re.M)
    assert covered==names
    assert "Lean.collectAxioms" in audit and "Lean.ofReduceBool" in audit and "sorryAx" in audit


def require_ci_audit_target(workflow):
    parsed = yaml.load(workflow, Loader=yaml.BaseLoader)
    selector = next(s["run"] for s in parsed["jobs"]["build"]["steps"]
                    if s.get("name") == "Detect changed Lean modules")
    block = re.search(r"(?ms)^targets=\((.*?)^\)", selector)
    assert block is not None
    words = shlex.split(block.group(1), comments=True)
    assert "Geometry.SourceReusableBusAxiomAudit" in words


def test_ci_enforces_audit_on_dependency_changes():
    require_ci_audit_target((codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8"))
    dedicated = yaml.load((codec.ROOT/".github/workflows/source-reusable-bus.yml")
                          .read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    for event in ("push", "pull_request"):
        for path in codec.PIN_PATHS + (".github/workflows/lean-ci.yml", ".gitattributes"):
            assert any(fnmatchcase(path, pattern) for pattern in dedicated["on"][event]["paths"])
    job = dedicated["jobs"]["controls"]
    assert job["strategy"]["matrix"]["os"] == ["ubuntu-latest"]
    assert "python -m pytest -q code/source_reusable_bus" in [s.get("run") for s in job["steps"]]


@pytest.mark.parametrize("replacement", [
    '            # "Geometry.SourceReusableBusAxiomAudit"',
    '            "Geometry.SourceReusableBus"', '',
])
def test_ci_guard_rejects_disabled_audit_target(replacement):
    workflow = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    target = '            "Geometry.SourceReusableBusAxiomAudit"'
    assert workflow.count(target) == 1
    with pytest.raises(AssertionError):
        require_ci_audit_target(workflow.replace(target, replacement))


@pytest.mark.parametrize("cycle", range(3))
@pytest.mark.parametrize("direction", [-1, 1])
def test_adversarial_signed_errors_on_every_register(cycle, direction):
    # Backward adjoints give a maximizing disturbance for this read. This
    # tests the abstract one-step allowance, not a physical noise mechanism
    # or the production integer rounding rule. Idle rails also receive error.
    count = 568*cycle+8
    e0, eta, rho = F(1,2**18), F(1,2**28)+F(1,2**21), F(1,2**16)
    edges = [verify.operation(k) for k in range(count)]
    weights = [F(0)]*10+[F(1,2), F(-1,2)]
    influence = [None]*count
    for k in reversed(range(count)):
        influence[k] = weights.copy()
        u,v = edges[k]
        weights[u] = weights[v] = (weights[u]+weights[v])/2
    def signed(w):
        return direction*(-1 if w < 0 else 1)
    exact = [F(3,2), F(5,2), F(5,2), F(3,2)]+[F(2)]*8
    perturbed = [x+e0*signed(w) for x,w in zip(exact, weights)]
    idle_drift = False
    for k,(u,v) in enumerate(edges):
        exact[u] = exact[v] = (exact[u]+exact[v])/2
        perturbed[u] = perturbed[v] = (perturbed[u]+perturbed[v])/2
        perturbed = [x+eta*signed(w) for x,w in zip(perturbed, influence[k])]
        assert max(abs(x-y) for x,y in zip(perturbed, exact)) <= e0+(k+1)*eta
        idle_drift |= any(i not in (u,v) and abs(perturbed[i]-exact[i]) > e0
                          for i in range(4))
    ideal_contrast = (exact[10]-exact[11])/2
    p, m = perturbed[10]+direction*rho, perturbed[11]-direction*rho
    error = (p-m)/2-ideal_contrast
    dual_bound = e0*sum(map(abs,weights))+eta*sum(sum(map(abs,w)) for w in influence)+rho
    assert error == direction*dual_bound
    assert dual_bound <= e0+count*eta+rho
    assert idle_drift
    use = (1,1,2)[cycle]
    assert verify.decode_receiver(p,m,use,1) == ("-1/2","1/2","-1/2")[cycle]


def test_cleanup_does_not_remove_common_mode_error():
    # Contraction concerns balanced amplitudes, not arbitrary raw offsets.
    offset = F(1,2**18)
    raw = [F(2)+offset]*12
    for k in range(8,568):
        u,v = verify.operation(k)
        raw[u] = raw[v] = (raw[u]+raw[v])/2
    assert raw[4:] == [2+offset]*8
    assert offset > offset*F(7,8)**80


def test_resource_ledger_counts_actual_edges_and_encoding_limits(controls):
    case = deepcopy(controls["executions"][0])
    ports = controls["global_ports"]
    original,_ = verify.resource_accounting(case,ports,15360,2**20)
    k = next(k for k,e in enumerate(case["tape"]) if ports[e[0]]//12 != ports[e[1]]//12)
    case["tape"].pop(k)
    fewer,_ = verify.resource_accounting(case,ports,15360,2**20)
    assert fewer["cross_carrier_means"] == original["cross_carrier_means"]-1
    assert fewer["means"] == original["means"]-1
    assert fewer["mean_reads"] == original["mean_reads"]-2
    assert fewer["writes"] == original["writes"]-2
    case["reads"][0]["local_units"] = [2**30,0]
    with pytest.raises(ValueError,match="receiver arithmetic overflow"):
        verify.resource_accounting(case,ports,15360,2**20)


@pytest.mark.parametrize("mutation,diagnostic", [
    ("promoted_receipt", "retained receipt"),
    ("control_bytes", "noncanonical artifact bytes: controls.json"),
    ("receipt_bytes", "noncanonical artifact bytes: receipt.json"),
    ("changed_proof", "source pins"),
    ("changed_tape", "native mean, rounding, value or consumed writer"),
])
def test_real_cli_rejects_tampered_evidence(tmp_path, mutation, diagnostic):
    for relative in codec.PIN_PATHS + ("code/source_reusable_bus/controls.json",
                                       "code/source_reusable_bus/receipt.json"):
        destination = tmp_path/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(codec.ROOT/relative, destination)
    command = [sys.executable, str(tmp_path/"code/source_reusable_bus/verify.py")]
    clean = subprocess.run(command,cwd=tmp_path,capture_output=True,text=True,encoding="utf-8")
    assert clean.returncode == 0, clean.stderr
    packet = tmp_path/"code/source_reusable_bus"
    if mutation == "promoted_receipt":
        p = packet/"receipt.json"
        data = codec.load(p)
        data["physical_M1_derived"] = True
        p.write_bytes(codec.canonical(data))
    elif mutation in ("control_bytes", "receipt_bytes"):
        p = packet/("controls.json" if mutation == "control_bytes" else "receipt.json")
        p.write_bytes(b" "+p.read_bytes())
    elif mutation == "changed_proof":
        p = tmp_path/"Lean/Geometry/SourceReusableBus.lean"
        p.write_bytes(p.read_bytes()+b"\n-- changed proof input\n")
    else:
        p = packet/"controls.json"
        data = codec.load(p)
        data["executions"][0]["tape"][0][-1] += 1
        p.write_bytes(codec.canonical(data))
    failed = subprocess.run(command,cwd=tmp_path,capture_output=True,text=True,encoding="utf-8")
    assert failed.returncode != 0
    assert diagnostic in failed.stderr
