from copy import deepcopy
from fractions import Fraction as F
from fnmatch import fnmatchcase
from pathlib import Path
import ast
import re
import shutil
import shlex
import subprocess
import sys

import pytest
import yaml

from . import build,codec,verify


@pytest.fixture(scope="module")
def artifact():
    return codec.load_artifact(codec.HERE/"controls.json")


def cli(controls,receipt,script=None):
    return subprocess.run([sys.executable,str(script or codec.HERE/"verify.py"),
                           "--controls",str(controls),"--receipt",str(receipt)],
                          capture_output=True,text=True,timeout=30)


def test_complete_independent_replay(artifact):
    codec.equal(build.build(),artifact,"producer drift")
    codec.equal(verify.verify(artifact),codec.load_artifact(codec.HERE/"receipt.json"),"receipt drift")
    result = cli(codec.HERE/"controls.json",codec.HERE/"receipt.json")
    assert result.returncode==0,result.stderr
    assert sum((codec.HERE/name).stat().st_size for name in ("controls.json","receipt.json"))<40_000


def test_no_producer_dependency(artifact,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("verifier called producer")
    for name in ("edges","execute","toy","probability_row","captured_controls","build"):
        monkeypatch.setattr(build,name,forbidden)
    verify.verify(artifact)
    tree = ast.parse((codec.HERE/"verify.py").read_text(encoding="utf-8"))
    assert all("build" not in ast.unparse(node) for node in ast.walk(tree)
               if isinstance(node,(ast.Import,ast.ImportFrom)))


MUTATIONS = [
    lambda a:a.update(scope="full_A1_A2_A3_countermodel"),
    lambda a:a.update(schema="wrong"),
    lambda a:a.update(extra="ignored"),
    lambda a:a.update(source_ports=[0]),
    lambda a:a.update(receiver_ports=[0,9]),
    lambda a:a.update(cut=[]),
    lambda a:a["cut"].pop(),
    lambda a:a["cut"].append(a["cut"][0]),
    lambda a:a["alphabet"].update(size=46080),
    lambda a:a["alphabet"].update(identity="directed_integer_completions"),
    lambda a:a["alphabet"].update(sha256="0"*64),
    lambda a:a["probabilities"].pop(),
    lambda a:a["probabilities"][0].update(horizon=True),
    lambda a:a["probabilities"][0].update(horizon=0.0),
    lambda a:a["probabilities"][2].update(lower_units=2**32),
    lambda a:a["probabilities"][2].update(upper_units=0),
    lambda a:a["probabilities"][2].update(denominator=2**33),
    lambda a:a["probabilities"][2].update(balanced_error_lower="1"),
    lambda a:a["probabilities"][2]["avoid_mass"].update(base="0"),
    lambda a:a["probabilities"][2]["avoid_mass"].update(exponent=3),
    lambda a:a["toy"].update(decoder="know_the_source"),
    lambda a:a["toy"].update(payloads=[1,1]),
    lambda a:a["toy"]["alphabet"].pop(),
    lambda a:a["toy"]["rows"].pop(),
    lambda a:a["toy"]["rows"][3].update(words=1),
    lambda a:a["toy"]["rows"][3].update(avoid_words=0),
    lambda a:a["toy"]["rows"][3].update(failure_counts=[0,0]),
    lambda a:a["toy"]["rows"][3].update(scalar_means=0),
    lambda a:a["toy"]["rows"][3].update(scalar_reads=0),
    lambda a:a["toy"]["rows"][3].update(scalar_writes=0),
    lambda a:a["toy"]["rows"][3].update(histories_sha256="0"*64),
    lambda a:a["captured_controls"].pop(),
    lambda a:a["captured_controls"][0].update(word=[]),
    lambda a:a["captured_controls"][0].update(means_per_case=1),
    lambda a:a["captured_controls"][0]["cases"].pop(),
    lambda a:a["captured_controls"][0]["cases"][0].update(tape=None),
    lambda a:a["captured_controls"][0]["cases"][0].update(tape=[]),
    lambda a:a["captured_controls"][0]["cases"][0].update(payload=1),
    lambda a:a["captured_controls"][1]["cases"][0].update(decoded=-1),
    lambda a:a["captured_controls"][0]["cases"][0]["initial"][0].__setitem__(1,"3"),
    lambda a:a["captured_controls"][0]["cases"][0].update(local=["3","1"]),
    lambda a:a["pins"].pop(codec.SUPPORT),
]


def reseal(a):
    for control in a["captured_controls"]:
        for case in control["cases"]:
            case["tape_sha256"] = codec.digest(case["tape"])


@pytest.mark.parametrize("mutation",MUTATIONS)
def test_resealed_semantic_mutations(artifact,mutation,tmp_path):
    changed = deepcopy(artifact)
    mutation(changed)
    reseal(changed)
    with pytest.raises((ValueError,TypeError,KeyError)):
        verify.verify(changed)
    controls,receipt = tmp_path/"controls.json",tmp_path/"receipt.json"
    controls.write_bytes(codec.canonical(changed))
    forged_receipt = codec.load_artifact(codec.HERE/"receipt.json")
    forged_receipt["controls_sha256"] = codec.digest(changed)
    receipt.write_bytes(codec.canonical(forged_receipt))
    result = cli(controls,receipt)
    assert result.returncode!=0
    assert "verification failed" in result.stderr


@pytest.mark.parametrize("control,event",[(0,0),(0,3),(0,7),(1,0),(1,7),(1,15)])
@pytest.mark.parametrize("field",range(7))
def test_actual_event_mutations(artifact,control,event,field):
    changed = deepcopy(artifact)
    row = changed["captured_controls"][control]["cases"][0]["tape"][event]
    row[field] = str(F(row[field])+1) if isinstance(row[field],str) else row[field]+1
    reseal(changed)
    with pytest.raises(ValueError,match="reconstructed"):
        verify.verify(changed)


@pytest.mark.parametrize("raw",["null","[]","{}","null\n","[]\n","{}\n",'{"x":1,"x":2}',
                                     '{"x":1.0}','{"x":NaN}','{"x":Infinity}'])
def test_trash_fails_cli(raw,tmp_path):
    controls = tmp_path/"controls.json"
    controls.write_text(raw,encoding="ascii")
    result = cli(controls,codec.HERE/"receipt.json")
    assert result.returncode!=0
    assert "verification failed" in result.stderr


def proof_closure():
    pending = ["Geometry.SourceReadSelectionAxiomAudit"]
    found = set()
    while pending:
        relative = "Lean/"+pending.pop().replace(".","/")+".lean"
        path = codec.ROOT/relative
        if relative in found or not path.is_file():
            continue
        found.add(relative)
        pending.extend(module for line in re.findall(r"^import (.+)$",path.read_text(encoding="utf-8"),re.M)
                       for module in line.split())
    return found


def test_complete_trust_pins_and_declaration_audit():
    assert proof_closure() <= set(codec.PINS)
    audit = (codec.ROOT/"Lean/Geometry/SourceReadSelectionAxiomAudit.lean").read_text(encoding="utf-8")
    declarations = set()
    for module in ("SourceReadSelection","SourceConstrainedSelection","SourceSelectionControls","SourceConstrainedRead"):
        proof = (codec.ROOT/f"Lean/Geometry/{module}.lean").read_text(encoding="utf-8")
        declarations.update(f"OPH.{module}.{name}" for name in re.findall(r"^theorem (\w+)",proof,re.M))
    audited = set(re.findall(r"^audit_reusable_bus_axioms (OPH\.[\w.]+)$",audit,re.M))
    assert declarations==audited
    assert len(declarations)==48


@pytest.mark.parametrize("relative",codec.PINS)
def test_real_source_tampering_rejected(relative,tmp_path):
    for source in (*codec.PINS,"code/source_read_selection/controls.json","code/source_read_selection/receipt.json"):
        destination = tmp_path/source
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(codec.ROOT/source,destination)
    package = tmp_path/"code/source_read_selection"
    result = cli(package/"controls.json",package/"receipt.json",package/"verify.py")
    assert result.returncode==0,result.stderr
    with (tmp_path/relative).open("ab") as file:
        file.write(b"\n")
    result = cli(package/"controls.json",package/"receipt.json",package/"verify.py")
    assert result.returncode!=0
    assert "source pins" in result.stderr


def test_positive_transport_and_omitted_history_control(artifact):
    row = artifact["toy"]["rows"][3]
    # The sole three-step successful word is 0,1,2. Conditioning on it would
    # discard 26/27 of this exact schedule law, including eight cut-avoiders.
    assert row["words"]==27 and row["avoid_words"]==8
    assert row["failure_counts"]==[26,26]
    assert artifact["toy"]["rows"][0]["failure_counts"]==[1,1]
    assert row["scalar_means"]==162
    for case in artifact["captured_controls"][0]["cases"]:
        assert case["decoded"]==case["payload"]
    for case in artifact["captured_controls"][1]["cases"]:
        assert case["decoded"]==0


def test_uniform_marginals_do_not_supply_free_word_law():
    # Choose one edge uniformly and repeat it. Every position has uniform
    # edge marginals, but the far receiver never sees either source input.
    for edge in ((0,1),(1,2),(2,3)):
        for payload in (-1,1):
            local,_ = build.execute({0:F(payload),1:F(0),2:F(0),3:F(0)},[edge]*6,[3])
            assert local==["0"]
    assert build.toy(6)["failure_counts"][0] < 3**6


def test_probability_endpoints_and_positive_bound(artifact):
    rows = artifact["probabilities"]
    assert rows[0]["lower_units"]==rows[0]["upper_units"]==2**32
    assert rows[0]["balanced_error_lower"]=="1/2"
    assert F(rows[-1]["balanced_error_lower"])>F(29,100)
    assert verify.probability_row(3,0,16)["balanced_error_lower"]=="1/2"
    assert verify.probability_row(3,3,16)["lower_units"]==0


def test_error_bound_is_not_a_bound_for_each_input():
    # A decoder always guessing -1 has zero error for -1 and error one for +1.
    # Only their sum/independent balanced average supports the stated bound.
    loss = [int(-1 != payload) for payload in (-1,1)]
    assert loss==[0,1]
    assert F(sum(loss),2)==F(1,2)


def require_ci_audit(source):
    workflow = yaml.safe_load(source)
    step = next(s for s in workflow["jobs"]["build"]["steps"] if s.get("id")=="changed_modules")
    roots = step["run"].split("targets=(",1)[1].split(")",1)[0]
    assert "Geometry.SourceReadSelectionAxiomAudit" in shlex.split(roots,comments=True)


def test_ci_audit_is_unconditional_and_controls_cover_dependencies():
    source = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    require_ci_audit(source)
    workflow = yaml.safe_load((codec.ROOT/".github/workflows/source-read-selection.yml").read_text())
    assert workflow["jobs"]["controls"]["strategy"]["matrix"]["os"]==["ubuntu-latest"]
    steps = workflow["jobs"]["controls"]["steps"]
    assert any(step.get("run")=="python -m pytest -q code/source_read_selection" for step in steps)
    triggers = workflow.get("on",workflow.get(True))
    for event in ("push","pull_request"):
        paths = triggers[event]["paths"]
        for relative in (*codec.PINS,"Lean/Geometry.lean",".github/workflows/lean-ci.yml"):
            assert any(fnmatchcase(relative,pattern) for pattern in paths),relative


@pytest.mark.parametrize("replacement",['','# "Geometry.SourceReadSelectionAxiomAudit"',
                                         '"Geometry.SourceReadSelection"'])
def test_missing_or_commented_trust_gate_fails(replacement):
    source = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    with pytest.raises(AssertionError):
        require_ci_audit(source.replace('"Geometry.SourceReadSelectionAxiomAudit"',replacement))


@pytest.mark.parametrize("field",["controls_sha256","status","toy_histories","toy_scalar_means",
                                    "captured_scalar_means","source_cut","probabilities"])
def test_forged_receipt_rejected(field,tmp_path):
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt[field] = None
    path = tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(receipt))
    result = cli(codec.HERE/"controls.json",path)
    assert result.returncode!=0 and "receipt mismatch" in result.stderr
