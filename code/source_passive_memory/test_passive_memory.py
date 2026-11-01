"""Retained-evidence replay and deliberate attempts to smuggle in feedback."""
from fractions import Fraction
import importlib
import importlib.util
from pathlib import Path
import re
import sys

import pytest
import yaml

# Unique package namespace: repository suites already contain unrelated modules
# named verify/build/codec, and Python itself has a module called code.
HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "oph_passive_memory_tests", HERE/"__init__.py", submodule_search_locations=[str(HERE)])
PACKAGE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PACKAGE
SPEC.loader.exec_module(PACKAGE)
codec = importlib.import_module(SPEC.name+".codec")
build = importlib.import_module(SPEC.name+".build")
verify = importlib.import_module(SPEC.name+".verify")


@pytest.fixture
def retained():
    return codec.load(HERE/"controls.json")


def test_retained_controls_independently_replay(retained):
    actual = verify.verify(retained)
    assert codec.canonical(actual) == (HERE/"receipt.json").read_bytes()
    assert actual["executed_means"] == 65
    assert actual["fibre"]["equal_image_witnesses"] == 6


def test_builder_reproduces_retained_bytes():
    assert codec.canonical(build.build()) == (HERE/"controls.json").read_bytes()


def test_verifier_does_not_call_builder(monkeypatch, retained):
    def broken(*args, **kwargs):
        raise AssertionError("verifier called the generating state machine")
    for name in ("build","trace","graph_census","fibre_certificate"):
        monkeypatch.setattr(build,name,broken)
    assert verify.verify(retained)["executed_means"] == 65


@pytest.mark.parametrize("field,value", [
    ("index", True), ("index",1), ("op","reset"), ("op","copy"),
    ("edge",[0,3]), ("edge",[0,999]), ("inputs",["0","0"]),
    ("writers",["init:0","event:0"]), ("output","0"), ("loss","0"),
])
def test_event_mutations_rejected(retained,field,value):
    retained["executions"][2]["events"][0][field] = value
    with pytest.raises(ValueError,match="read/write"):
        verify.verify(retained)


@pytest.mark.parametrize("metric,value", [
    ("registers",1), ("initializations",1), ("means",3), ("reads",4),
    ("writes_including_initialization",8), ("touches",[4,0,0,0,0]),
    ("total_load","0"), ("initial_quadratic","0"),
    ("final_quadratic","0"), ("loss","0"),
])
def test_resource_undercounts_rejected(retained,metric,value):
    retained["executions"][2]["metrics"][metric] = value
    with pytest.raises(ValueError,match="resource metrics"):
        verify.verify(retained)


@pytest.mark.parametrize("mutation",["elide_mean","erase_ancilla","hide_ancilla",
                                    "stale_final_writer","widen_menu","claim_trajectory"])
def test_hidden_preparation_or_operation_rejected(retained,mutation):
    case = retained["executions"][2]
    if mutation == "elide_mean":
        case["events"].pop()
    elif mutation == "erase_ancilla":
        case["final"][1] = "0"
    elif mutation == "hide_ancilla":
        case["final"].pop()
    elif mutation == "stale_final_writer":
        case["final_writers"][0] = "init:0"
    elif mutation == "widen_menu":
        case["allowed_edges"].append([1,2])
    else:
        retained["fibre"]["kind"] = "finite_pair_mean_execution"
    with pytest.raises(ValueError):
        verify.verify(retained)


def test_consistent_wrong_experiment_rejected(retained):
    # This computes a valid native mean trace with the same source residual,
    # operation counts, total load and Q. It assigns the ancilla loads to
    # different ports, so it must not pass as the retained, specified word.
    word = [(0,2),(0,1),(0,3),(0,4)]
    wrong = build.trace("cold_star_4",[Fraction(3,2),0,0,0,0],word,word)
    original = retained["executions"][2]
    for key in ("means","registers","total_load","initial_quadratic","final_quadratic","loss"):
        assert wrong["metrics"][key] == original["metrics"][key]
    assert wrong["final"][0] == original["final"][0]
    retained["executions"][2] = wrong
    with pytest.raises(ValueError,match="read/write"):
        verify.verify(retained)


def test_initial_preparation_not_free(retained):
    word = [(2,4),(3,5),(2,3)]
    retained["executions"][6] = build.trace("prepared_balanced_copy_and_cleanup",
        [5,3,4,4,4,4],word,word,baseline=4)
    with pytest.raises(ValueError,match="initial state"):
        verify.verify(retained)


def test_captured_port_permutation_rejected(retained):
    retained["executions"][-1]["global_ports"][1:3] = [5,1]
    with pytest.raises(ValueError,match="port identity"):
        verify.verify(retained)


@pytest.mark.parametrize("field,value",[("edges",46049),("component_sizes",[7680,7680])])
def test_captured_graph_census_is_recomputed(retained,field,value):
    retained["support"][field] = value
    with pytest.raises(ValueError,match="support census"):
        verify.verify(retained)


@pytest.mark.parametrize("field,value",[
    ("right",["1","2","3","4"]), ("left",["1","2","3","0"]),
    ("common_image",["1","2","7","0"]), ("amount","5"), ("edge",[3,0]),
])
def test_fibre_witness_mutations_rejected(retained,field,value):
    retained["fibre"]["witnesses"][0][field] = value
    with pytest.raises(ValueError):
        verify.verify(retained)


def test_source_pin_drift_rejected(retained):
    retained["source_sha256"][codec.SUPPORT] = "0"*64
    with pytest.raises(ValueError,match="source pins"):
        verify.verify(retained)


@pytest.mark.parametrize("source", ['{"a":1,"a":2}', '{"a":1.0}', '{"a":NaN}', '{"a":Infinity}'])
def test_ambiguous_json_rejected(tmp_path,source):
    path = tmp_path/"bad.json"
    path.write_text(source,encoding="utf-8")
    with pytest.raises(ValueError):
        codec.load(path)


@pytest.mark.parametrize("value",[0, True, "0.5", "2/4", "1/0", "NaN"])
def test_rational_encoding_is_exact_and_unique(value):
    with pytest.raises(ValueError):
        codec.rational(value)


def test_small_words_keep_positive_registers_and_charge_all_loss():
    # Exhaust all three-edge words of length <= 4 on three registers, including
    # no-op means, and check each prefix with the independent replay engine.
    from itertools import product
    for length in range(5):
        for word in product([(0,1),(0,2),(1,2)],repeat=length):
            word = list(word)
            initial = list(map(Fraction,[2,0,1]))
            trace = build.trace("exhaustive",initial,word,word)
            verify.replay(trace,("exhaustive",initial,word,0,None))


def test_axiom_audit_covers_all_theorems():
    names = set()
    for module in codec.PROOFS:
        source = (codec.ROOT/f"Lean/Geometry/{module}.lean").read_text(encoding="utf-8")
        names.update(f"OPH.{module}.{name}" for name in re.findall(r"^theorem ([\w.]+)",source,re.M))
    audit = (codec.ROOT/"Lean/Geometry/SourcePassiveMemoryAxiomAudit.lean").read_text(encoding="utf-8")
    audited = re.findall(r"^audit_passive_axioms (OPH\.[\w.]+)$",audit,re.M)
    assert len(audited) == len(set(audited)) == 40
    assert set(audited) == names
    assert "import Geometry.SourcePassiveMemoryAxiomAudit" in (codec.ROOT/"Lean/Geometry.lean").read_text(encoding="utf-8")


def test_axiom_audit_runs_when_only_a_dependency_changes():
    workflow = yaml.load((codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8"),
                         Loader=yaml.BaseLoader)
    selector = next(step["run"] for step in workflow["jobs"]["build"]["steps"]
                    if step.get("id") == "changed_modules")
    # Only explicit default targets run when the audit module itself is absent
    # from git diff. Importing it into an unchanged umbrella is insufficient.
    defaults = re.search(r"(?ms)^targets=\(\n(.*?)^\)", selector)
    assert defaults is not None
    assert '"Geometry.SourcePassiveMemoryAxiomAudit"' in defaults.group(1)


def test_dedicated_ci_preserves_frozen_runner():
    projection = codec.load(codec.ROOT/"code/invariant_mining/outputs/source_projection.json")
    pin = next(p for p in projection["control_documents"] if p["path"] == "tools/run_mandatory_suite.py")
    assert "sha256:"+codec.digest(codec.ROOT/pin["path"]) == pin["sha256"]
    workflow = yaml.load((codec.ROOT/".github/workflows/source-passive-memory.yml").read_text(encoding="utf-8"),
                         Loader=yaml.BaseLoader)
    job = workflow["jobs"]["controls"]
    assert set(job["strategy"]["matrix"]["os"]) == {"ubuntu-latest"}
    assert any(s.get("run") == "python -m pytest -q code/source_passive_memory" for s in job["steps"])
    for trigger in ("push","pull_request"):
        assert "code/source_passive_memory/**" in workflow["on"][trigger]["paths"]
        assert "code/source_routing/support_w12_l3.json" in workflow["on"][trigger]["paths"]
        assert ".github/workflows/lean-ci.yml" in workflow["on"][trigger]["paths"]
