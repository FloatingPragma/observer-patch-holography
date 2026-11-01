"""Scientific replay, mutation gates and encoded-interface boundary controls."""
from fractions import Fraction as F
import importlib
import importlib.util
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import pytest
import yaml

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "oph_encoded_memory_tests", HERE/"__init__.py", submodule_search_locations=[str(HERE)])
PACKAGE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PACKAGE
SPEC.loader.exec_module(PACKAGE)
codec = importlib.import_module(SPEC.name+".codec")
build = importlib.import_module(SPEC.name+".build")
verify = importlib.import_module(SPEC.name+".verify")


@pytest.fixture
def packet():
    return codec.load(HERE/"controls.json")


def test_retained_evidence_replays(packet):
    receipt = verify.verify(packet)
    assert codec.canonical(receipt) == (HERE/"receipt.json").read_bytes()
    assert (receipt["executed_means"], receipt["checked_copy_states"]) == (626, 209)


def test_builder_reproduces_exact_bytes():
    assert codec.canonical(build.build()) == (HERE/"controls.json").read_bytes()


def test_independent_replay_does_not_call_producer(packet, monkeypatch):
    def broken(*args, **kwargs):
        raise AssertionError("producer called during independent replay")
    for name in ("build", "execute", "specifications", "margins"):
        monkeypatch.setattr(build, name, broken)
    verify.verify(packet)


@pytest.mark.parametrize("key,value", [
    ("index", True), ("index", 1), ("op", "copy"), ("op", "reset"),
    ("edge", [0, 3]), ("edge", [0, 500]), ("inputs", ["2", "2"]),
    ("writers", ["init:0", "init:0"]), ("output", "0"),
    ("quadratic_loss", "0"),
])
def test_forged_native_events_rejected(packet, key, value):
    packet["executions"][1]["events"][0][key] = value
    with pytest.raises(ValueError, match="native event"):
        verify.verify(packet)


@pytest.mark.parametrize("key", ["scalar_registers", "preparation_writes", "means",
                                  "scalar_reads", "scalar_writes", "peak_scalar_value_bits"])
def test_hidden_resource_cost_rejected(packet, key):
    packet["executions"][1]["resources"][key] = 0
    with pytest.raises(ValueError, match="resource account"):
        verify.verify(packet)


@pytest.mark.parametrize("mutation", ["missing_clear", "invented_prepare", "stale_final_writer",
                                     "hide_amplitude_loss", "false_checkpoint", "fake_embedding",
                                     "wrong_blank", "extra_field", "source_pin"])
def test_interface_and_custody_forgeries_rejected(packet, mutation):
    c = packet["executions"][2]
    if mutation == "missing_clear":
        c["events"].pop(2)
    elif mutation == "invented_prepare":
        c["events"].insert(3, {"op": "prepare", "value": "2"})
    elif mutation == "stale_final_writer":
        c["final_writers"][0] = "init:0"
    elif mutation == "hide_amplitude_loss":
        c["final"][:2] = c["initial"][:2]
    elif mutation == "false_checkpoint":
        c["copy_checkpoints"][0]["rails"][0] = "7/2"
    elif mutation == "fake_embedding":
        c["global_ports"] = [12, 13, 17, 21]
    elif mutation == "wrong_blank":
        c["initial"][2] = "7/2"
    elif mutation == "extra_field":
        c["physical_M1_derived"] = True
    else:
        packet["source_sha256"][codec.SUPPORT] = "0"*64
    with pytest.raises(ValueError):
        verify.verify(packet)


def test_coherent_different_experiment_is_rejected(packet):
    # A valid negative preparation cannot replace the specified positive run.
    c = build.execute("reread_positive", "read", 16, F(-3, 2), "captured_square")
    packet["executions"][2] = c
    with pytest.raises(ValueError, match="prepared input"):
        verify.verify(packet)


@pytest.mark.parametrize("field,value", [("sufficient_sign_margin", True),
                                         ("total_bound", "0"),
                                         ("opposite_bits_terminal_ambiguity_radius", "1")])
def test_false_noise_certificates_rejected(packet, field, value):
    packet["analytic_margin_controls_not_noisy_executions"][-1][field] = value
    with pytest.raises(ValueError, match="analytic margins"):
        verify.verify(packet)


@pytest.mark.parametrize("value", [True, 1, 0.5, "2/2", "1.0", "1/0", "NaN"])
def test_noncanonical_rational_rejected(value):
    with pytest.raises(ValueError):
        codec.rational(value)


@pytest.mark.parametrize("data", ['{"x":1,"x":2}', '{"x":0.5}', '{"x":NaN}'])
def test_ambiguous_json_rejected(tmp_path, data):
    p = tmp_path/"bad.json"
    p.write_text(data, encoding="utf-8")
    with pytest.raises(ValueError):
        codec.load(p)


def mean(x, u, v):
    """A third, test-only arithmetic model; no tape or verifier imports."""
    y = list(x)
    y[u] = y[v] = (x[u]+x[v])/2
    return y


def test_unknown_payload_copy_and_cleanup_for_a_rational_family():
    for a in (F(n, 7) for n in range(-14, 15)):
        x = [2+a, 2-a, F(2), F(2)]
        for _ in range(8):
            x = mean(mean(x, 0, 2), 1, 3)
            a /= 2
            assert x == [2+a, 2-a, 2+a, 2-a]
            x = mean(x, 2, 3)
            assert x == [2+a, 2-a, 2, 2]


def test_preconditions_cannot_be_removed():
    # An occupied opposite target destroys both signals.
    assert mean(mean([F(3), F(1), F(1), F(3)], 0, 2), 1, 3) == [2]*4
    # Arbitrary means on the connected menu are not allowed during storage.
    assert mean([F(3), F(1), F(2), F(2)], 0, 1) == [2]*4
    # Common-mode error survives encoded reset; it is not a precision refresh.
    assert mean([F(31, 10), F(11, 10)], 0, 1) == [F(21, 10)]*2


def test_bounded_noisy_reuse_includes_preparation_and_idle_drift():
    delta, e0, rho = F(1, 2**18), F(1, 2**19), F(1, 2**20)
    for sign in (-1, 1):
        exact = [2+sign*F(3, 2), 2-sign*F(3, 2), F(2), F(2)]
        noisy = [x+(-1)**i*e0 for i, x in enumerate(exact)]
        for k in range(36):
            # Reuse the same four rails; includes error on untouched rails.
            source = (k//3) % 2
            u, v = ((2*source, 2*(1-source)), (2*source+1, 2*(1-source)+1),
                    (2*source, 2*source+1))[k % 3]
            exact = mean(exact, u, v)
            noisy = [x+(-1)**(k+i)*delta*F(i+1, 4)
                     for i, x in enumerate(mean(noisy, u, v))]
            assert max(abs(x-y) for x, y in zip(noisy, exact)) <= e0+(k+1)*delta
            if k % 3 == 2:
                h = (k+1)//3
                t = h % 2
                amplitude = F(3, 2**(h+1))
                if amplitude > e0+(k+1)*delta+rho:
                    yp, ym = noisy[2*t]+rho, noisy[2*t+1]-rho
                    assert (yp > ym) == (sign > 0)


def test_sharp_terminal_ambiguity_and_analog_gain():
    for h in (0, 1, 12, 64):
        amp = F(3, 2**(h+1))
        plus, minus = [2+amp, 2-amp], [2-amp, 2+amp]
        assert max(abs(x-2) for x in plus) == max(abs(x-2) for x in minus) == amp
        e = F(1, 10**6)
        assert 2**h*(amp+e)-F(3, 2) == 2**h*e


def test_captured_locality_boundary_with_positive_finder_control():
    support = codec.load(codec.ROOT/codec.SUPPORT)
    rungs = {tuple(sorted(e)) for e in support["intra_carrier_seams"]}

    def square_count(glued):
        by_carriers = {}
        for a, u, b, v in glued:
            if a > b:
                a, u, b, v = b, v, a, u
            by_carriers.setdefault((a, b), set()).add((u, v))
        count = 0
        for pairs in by_carriers.values():
            for u, v in pairs:
                for w, z in pairs:
                    if u < w and v != z and (u, w) in rungs and tuple(sorted((v, z))) in rungs:
                        count += 1
        return count

    assert {(0, 1), (5, 9), (0, 5), (1, 9)} <= rungs
    assert square_count(support["glued_pairs"]) == 0
    assert square_count([[0, 0, 1, 5], [0, 1, 1, 9]]) == 1


def require_downstream_audit(workflow):
    # Read the executable selector, not a matching string elsewhere in YAML.
    doc = yaml.load(workflow, Loader=yaml.BaseLoader)
    selector = next(step["run"] for step in doc["jobs"]["build"]["steps"]
                    if step.get("name") == "Detect changed Lean modules")
    block = re.search(r"(?ms)^targets=\((.*?)^\)", selector)
    assert block is not None, "unconditional default target list is missing"
    # Shell comments are not targets. The old raw-substring guard could be
    # fooled by commenting out the actual downstream audit line.
    active = shlex.split(block.group(1), comments=True)
    assert "Geometry.SourceEncodedMemoryAxiomAudit" in active


@pytest.mark.parametrize("replacement", [
    '            # "Geometry.SourceEncodedMemoryAxiomAudit"',
    '            "Geometry.SourceEncodedMemory"',
    '',
])
def test_removing_or_commenting_the_ci_audit_target_fails(replacement):
    workflow = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    original = '            "Geometry.SourceEncodedMemoryAxiomAudit"'
    assert workflow.count(original) == 1
    with pytest.raises(AssertionError):
        require_downstream_audit(workflow.replace(original, replacement))


@pytest.mark.parametrize("mutation,diagnostic", [
    ("promoted_receipt", "retained receipt mismatch"),
    ("noncanonical_bytes", "noncanonical controls bytes"),
    ("changed_proof", "source pins mismatch"),
])
def test_real_cli_rejects_tampered_artifacts_in_isolation(tmp_path, mutation, diagnostic):
    # The real entry point must fail, not only an in-memory helper. Copy all
    # pinned inputs byte for byte; never alter the actual retained artifacts.
    for relative in codec.PINS + ("code/source_encoded_memory/controls.json",
                                  "code/source_encoded_memory/receipt.json"):
        destination = tmp_path/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(codec.ROOT/relative, destination)
    command = [sys.executable, str(tmp_path/"code/source_encoded_memory/verify.py")]
    clean = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, encoding="utf-8")
    assert clean.returncode == 0, clean.stderr
    if mutation == "promoted_receipt":
        p = tmp_path/"code/source_encoded_memory/receipt.json"
        data = codec.load(p)
        data["physical_M1_derived"] = True
        p.write_bytes(codec.canonical(data))
    elif mutation == "noncanonical_bytes":
        p = tmp_path/"code/source_encoded_memory/controls.json"
        p.write_bytes(b" " + p.read_bytes())
    else:
        p = tmp_path/"Lean/Geometry/SourceEncodedMemory.lean"
        p.write_bytes(p.read_bytes() + b"\n-- changed proof input\n")
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, encoding="utf-8")
    assert result.returncode != 0
    assert diagnostic in result.stderr


def test_no_hidden_untested_theorems_and_downstream_audit():
    lean = codec.ROOT/"Lean"
    proof = (lean/"Geometry/SourceEncodedMemory.lean").read_text(encoding="utf-8")
    audit = (lean/"Geometry/SourceEncodedMemoryAxiomAudit.lean").read_text(encoding="utf-8")
    names = set(re.findall(r"^theorem (\w+)", proof, re.M))
    covered = set(re.findall(r"^audit_encoded_axioms OPH\.SourceEncodedMemory\.(\w+)", audit, re.M))
    assert names == covered and len(names) >= 25
    declarations = re.sub(r"/\-.*?\-/", "", proof, flags=re.S)
    declarations = re.sub(r"--[^\n]*", "", declarations)
    assert not re.search(r"\b(sorry|admit|axiom|native_decide)\b", declarations)
    assert "Lean.collectAxioms" in audit and "unexpected.isEmpty" in audit
    assert "audit_encoded_axioms sorryAx" in audit
    assert "audit_encoded_axioms Lean.ofReduceBool" in audit
    assert "import Geometry.SourceEncodedMemoryAxiomAudit" in (lean/"Geometry.lean").read_text(encoding="utf-8")
    workflow = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    # Dependency-only changes must rerun this importer; changed-file builds do not suffice.
    require_downstream_audit(workflow)


def test_dedicated_ci_covers_receipts_and_every_input():
    path = codec.ROOT/".github/workflows/source-encoded-memory.yml"
    workflow = yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    assert workflow["jobs"]["controls"]["strategy"]["matrix"]["os"] == ["ubuntu-latest"]
    commands = [step.get("run", "") for step in workflow["jobs"]["controls"]["steps"]]
    assert "python -m pytest -q code/source_encoded_memory" in commands
    from fnmatch import fnmatchcase
    for event in ("push", "pull_request"):
        paths = workflow["on"][event]["paths"]
        for source in codec.PINS + (".github/workflows/lean-ci.yml", ".gitattributes"):
            assert any(fnmatchcase(source, pattern) for pattern in paths), source
