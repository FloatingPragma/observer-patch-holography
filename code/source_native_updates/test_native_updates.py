from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from fnmatch import fnmatchcase
from pathlib import Path
import ast
import re
import shlex
import shutil
import subprocess
import sys

import pytest
import yaml

from . import build, codec, verify


@pytest.fixture(scope="module")
def artifact():
    return codec.load_artifact(codec.HERE/"controls.json")


def test_full_replay_and_reproducible_bytes(artifact):
    codec.equal(build.build(), artifact, "producer drift")
    codec.equal(verify.verify(artifact), codec.load_artifact(codec.HERE/"receipt.json"), "receipt drift")
    assert len(artifact["representative_tape"]) == 2279
    assert (codec.HERE/"controls.json").stat().st_size < 250_000


def test_verifier_does_not_call_producer(artifact, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("verifier called producer")
    monkeypatch.setattr(build, "program", forbidden)
    monkeypatch.setattr(build, "history", forbidden)
    monkeypatch.setattr(build, "build", forbidden)
    verify.verify(artifact)
    tree = ast.parse((codec.HERE/"verify.py").read_text(encoding="utf-8"))
    imports = [node for node in ast.walk(tree) if isinstance(node, (ast.Import,ast.ImportFrom))]
    assert all("build" not in ast.unparse(node) for node in imports)


@pytest.mark.parametrize("mutation", [
    lambda a: a.update(schema="wrong"),
    lambda a: a.update(cases=[]),
    lambda a: a["cases"].pop(),
    lambda a: a["cases"].append(deepcopy(a["cases"][0])),
    lambda a: a["cases"].__setitem__(1,deepcopy(a["cases"][0])),
    lambda a: a.update(representative_tape=None),
    lambda a: a.update(representative_tape=[]),
    lambda a: a["representative_tape"].pop(),
    lambda a: a["representative_tape"].append(a["representative_tape"][-1]),
    lambda a: a.update(representative_payloads=["-3/2","1/2"]),
    lambda a: a["ports"].__setitem__(0,True),
    lambda a: a["ports"].__setitem__(1,a["ports"][0]),
    lambda a: a["ports"].__setitem__(0,12),
    lambda a: a["cases"][0].update(event_count=True),
    lambda a: a["cases"][0].update(event_count=2279.0),
    lambda a: a["cases"][0].update(payloads=["-6/4","-3/2"]),
    lambda a: a["cases"][0]["initial_units"].__setitem__(4,2**20),
    lambda a: a["cases"][0]["commit_units"].__setitem__(4,0),
    lambda a: a["cases"][0]["final_units"].__setitem__(4,0),
    lambda a: a["cases"][0]["reads"][0].update(record=0),
    lambda a: a["cases"][0]["reads"][0].update(use=0),
    lambda a: a["cases"][0]["reads"][0].update(exponent=5),
    lambda a: a["cases"][0]["reads"][0].update(decoded="3"),
    lambda a: a["cases"][0]["reads"][0].update(after_mean=7),
    lambda a: a["cases"][0]["reads"][0]["writers"].__setitem__(0,0),
    lambda a: a["cases"][0].update(tape_sha256="0"*64),
    lambda a: a["pins"].__setitem__(codec.SUPPORT,"0"*64),
    lambda a: a.update(extra=True),
])
def test_semantic_mutations(artifact, mutation):
    changed = deepcopy(artifact)
    mutation(changed)
    with pytest.raises((ValueError,TypeError,KeyError)):
        verify.verify(changed)


@pytest.mark.parametrize("event,field", list(product((0,2,4,6,7,14,15,574,1150,2278),range(7))))
def test_resealed_event_mutations(artifact,event,field):
    changed = deepcopy(artifact)
    changed["representative_tape"][event][field] += 1
    # The altered producer's digest is coherent with its altered evidence.
    changed["cases"][8]["tape_sha256"] = codec.digest(changed["representative_tape"])
    with pytest.raises(ValueError):
        verify.verify(changed)


def test_valid_alternative_program_rejected(monkeypatch):
    original = build.program()
    # Both independent exports still compute the same sum. This is a valid
    # native program, but is not the exact controller specified by the packet.
    alternate = original[2:4]+original[:2]+original[4:]
    monkeypatch.setattr(build,"program",lambda: alternate)
    changed = build.build()
    with pytest.raises(ValueError,match="independently replayed"):
        verify.verify(changed)


@pytest.mark.parametrize("bad", ['{"x":1,"x":2}', '{"x":1.0}', '{"x":NaN}', '{"x":Infinity}'])
def test_strict_json(tmp_path,bad):
    path = tmp_path/"bad.json"
    path.write_text(bad,encoding="ascii")
    with pytest.raises(ValueError):
        codec.load(path)


@pytest.mark.parametrize("bad", [True,1,1.0,"0/1","2/4","+1","1.0"])
def test_strict_rationals(bad):
    with pytest.raises(ValueError):
        codec.rational(bad)


def cli(controls,receipt):
    return subprocess.run([sys.executable,str(codec.HERE/"verify.py"),
                           "--controls",str(controls),"--receipt",str(receipt)],
                          capture_output=True,text=True,timeout=30)


def test_cli_success():
    result = cli(codec.HERE/"controls.json",codec.HERE/"receipt.json")
    assert result.returncode == 0, result.stderr


def local_proof_import_closure():
    pending = ["Geometry.SourceNativeUpdatesAxiomAudit"]
    seen = set()
    while pending:
        module = pending.pop()
        path = "Lean/"+module.replace(".","/")+".lean"
        if path in seen or not (codec.ROOT/path).is_file():
            continue
        seen.add(path)
        source = (codec.ROOT/path).read_text(encoding="utf-8")
        pending.extend(word for imports in re.findall(r"^import (.+)$",source,re.M)
                       for word in imports.split())
    return seen


def test_all_local_proof_dependencies_are_pinned():
    assert local_proof_import_closure() <= set(codec.PINS)


@pytest.mark.parametrize("path", [
    "Lean/Geometry/SourceReusableBusAxiomAudit.lean",
    "Lean/Geometry/SourceNativeUpdatesAxiomAudit.lean",
    "Lean/Geometry/SourceBusScaling.lean",
    "Lean/Geometry/SourceNativeRecords.lean",
    "Lean/ObserverPatchHolography/ScalarSeamRepair.lean",
    "Lean/lean-toolchain",
    codec.SUPPORT,
])
def test_real_cli_rejects_source_file_tampering(tmp_path,path):
    # Use an isolated checkout surface and a fresh interpreter. This changes
    # actual source bytes, not just a claimed hash in an artifact.
    paths = set(codec.PINS) | local_proof_import_closure() | {
        "code/source_native_updates/controls.json", "code/source_native_updates/receipt.json"}
    for relative in paths:
        target = tmp_path/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(codec.ROOT/relative,target)
    command = [sys.executable,str(tmp_path/"code/source_native_updates/verify.py")]
    baseline = subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert baseline.returncode == 0, baseline.stderr
    target = tmp_path/path
    target.write_bytes(target.read_bytes()+b"\n")
    changed = subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert changed.returncode != 0
    assert "source pins" in changed.stderr


@pytest.mark.parametrize("kind", ["bytes","optimistic_margin","omitted_history","optimistic_cost",
                                      "lost_ancestry","fake_commit","missing_trace","bool_event"])
def test_cli_rejects_bad_evidence(artifact,tmp_path,kind):
    controls = deepcopy(artifact)
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    if kind == "optimistic_margin":
        receipt["read_bounds"][2]["margin"] = "1"
    elif kind == "omitted_history":
        receipt["histories"].pop()
    elif kind == "optimistic_cost":
        receipt["resources"]["means_per_history"] = 1
    elif kind == "lost_ancestry":
        receipt["histories"][0]["reads"][0]["initial_ancestors"] = [0]
    elif kind == "fake_commit":
        controls["cases"][0]["commit_units"][4] += 1
    elif kind == "missing_trace":
        controls["representative_tape"] = None
    elif kind == "bool_event":
        controls["representative_tape"][0][0] = False
    receipt["controls_sha256"] = codec.digest(controls)
    cp,rp = tmp_path/"controls.json",tmp_path/"receipt.json"
    cp.write_bytes(codec.canonical(controls)+(b" " if kind == "bytes" else b""))
    rp.write_bytes(codec.canonical(receipt))
    result = cli(cp,rp)
    assert result.returncode != 0
    assert "Verified native commit" not in result.stdout


@pytest.mark.parametrize("payloads", [(F(-3,2),F(-1,2)),(F(1,2),F(-3,2)),(F(3,2),F(3,2))])
def test_adversarial_noise_in_full_native_program(payloads):
    q = 2**20
    state = [F(2)]*14
    for j,amplitude in enumerate(payloads):
        state[2*j] += amplitude
        state[2*j+1] -= amplitude
    state = [x+(F(1,2**18) if i%2 else -F(1,2**18)) for i,x in enumerate(state)]
    a,b = payloads
    for k in range(2279):
        u,v = verify.operation(k)
        output = F(round((state[u]+state[v])*q/2),q)
        state[u] = state[v] = output
        state = [x+(F(1,2**28) if (k+3*i)%5 < 2 else -F(1,2**28))
                 for i,x in enumerate(state)]
        if k >= 7 and (k-7)%568 == 7:
            cycle = (k-7)//568
            report = verify.bounds(cycle,verify.ideal_reads()[cycle])
            expected = (a,b,a+b)[report["record"]]/2**report["exponent"]
            sign = 1 if (state[12]-state[13])/2 >= expected else -1
            measured = (state[12]-state[13])/2+sign*F(1,2**16)
            assert abs(measured-expected) <= F(report["total_error"])
            assert abs(measured-expected) < F(report["half_spacing"])


def test_low_precision_fails_and_actual_decoding_is_wrong():
    q = 2**8
    wrong = ambiguous = 0
    for a,b in product((F(-3,2),F(-1,2),F(1,2),F(3,2)),repeat=2):
        state = [2*q]*14
        state[:4] = [int((2+a)*q),int((2-a)*q),int((2+b)*q),int((2-b)*q)]
        for k in range(2279):
            u,v = verify.operation(k)
            state[u] = state[v] = round(F(state[u]+state[v],2))
            if k >= 7 and (k-7)%568 == 7:
                cycle = (k-7)//568
                report = verify.bounds(cycle,verify.ideal_reads()[cycle],q)
                assert F(report["margin"]) < 0
                expected = (a,b,a+b)[report["record"]]
                observation = F(state[12]-state[13],2*q)*2**report["exponent"]
                alphabet = list(map(F,range(-3,4))) if report["record"] == 2 else [F(x,2) for x in (-3,-1,1,3)]
                candidates = sorted((abs(observation-value),value) for value in alphabet)
                if candidates[0][0] == candidates[1][0]:
                    ambiguous += 1
                else:
                    wrong += candidates[0][1] != expected
    # Exceeding half spacing is not itself a decoding failure: saturation
    # can still recover the correct extreme codeword. Check actual decisions.
    assert wrong == 12
    assert ambiguous == 16


@pytest.mark.parametrize("n",[0,1,2,3,8,16])
def test_chain_basis_and_signed_weight_bound(n):
    weights = [F((n+1)**2-(n-i)**2) for i in range(n+1)]
    rate = 1-F(1,(n+1)**2)
    columns = []
    for j in range(n+1):
        values = [F(i == j) for i in range(n+1)]
        values[0] = 0
        for i in range(n):
            values[i] = values[i+1] = (values[i]+values[i+1])/2
        columns.append(values)
    for i in range(n+1):
        assert all(column[i] >= 0 for column in columns)
        assert sum(column[i]*w for column,w in zip(columns,weights)) <= rate*weights[i]
    if n == 0:
        assert columns == [[0]]


def test_every_theorem_transitively_audited():
    lean = codec.ROOT/"Lean/Geometry"
    audit = (lean/"SourceNativeUpdatesAxiomAudit.lean").read_text(encoding="utf-8")
    expected = []
    for module in ("SourceBusScaling","SourceNativeRecords"):
        source = (lean/(module+".lean")).read_text(encoding="utf-8")
        assert not re.search(r"\b(sorry|axiom|native_decide)\b",source)
        expected.extend("OPH."+module+"."+name for name in re.findall(r"^theorem (\w+)",source,re.M))
    actual = re.findall(r"^audit_reusable_bus_axioms (\S+)",audit,re.M)
    assert sorted(actual) == sorted(expected)
    assert "import Geometry.SourceNativeUpdatesAxiomAudit" in (lean.parent/"Geometry.lean").read_text(encoding="utf-8")


def require_ci_audit(workflow):
    parsed = yaml.load(workflow,Loader=yaml.BaseLoader)
    step = next(s for s in parsed["jobs"]["build"]["steps"] if s.get("name") == "Detect changed Lean modules")
    assert "if" not in step and "continue-on-error" not in step
    block = re.search(r"(?ms)^targets=\((.*?)^\)",step["run"])
    assert block is not None
    assert "Geometry.SourceNativeUpdatesAxiomAudit" in shlex.split(block.group(1),comments=True)


def test_ci_dependency_coverage():
    require_ci_audit((codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8"))
    workflow = yaml.load((codec.ROOT/".github/workflows/source-native-updates.yml")
                         .read_text(encoding="utf-8"),Loader=yaml.BaseLoader)
    for event in ("push","pull_request"):
        for path in codec.PINS+(".github/workflows/lean-ci.yml",".gitattributes"):
            assert any(fnmatchcase(path,pattern) for pattern in workflow["on"][event]["paths"])
    job = workflow["jobs"]["controls"]
    assert job["strategy"]["matrix"]["os"] == ["ubuntu-latest"]
    assert "if" not in job and "continue-on-error" not in job
    commands = [s for s in job["steps"] if s.get("run") == "python -m pytest -q code/source_native_updates"]
    assert len(commands) == 1 and "if" not in commands[0] and "continue-on-error" not in commands[0]


@pytest.mark.parametrize("replacement", ["", '# "Geometry.SourceNativeUpdatesAxiomAudit"',
                                        '"Geometry.SourceNativeRecords"'])
def test_ci_audit_guard_negative_controls(replacement):
    workflow = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    target = '"Geometry.SourceNativeUpdatesAxiomAudit"'
    assert workflow.count(target) == 1
    with pytest.raises(AssertionError):
        require_ci_audit(workflow.replace(target,replacement))
