"""Independent replay, actual CLI attacks, input retention and noisy publication."""
from copy import deepcopy
from fractions import Fraction as F
from fnmatch import fnmatchcase
import inspect
from itertools import product
from pathlib import Path
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


@pytest.fixture(scope="module")
def edges():
    return verify.support_edges()


def cli(controls,receipt,script=None):
    return subprocess.run([sys.executable,str(script or codec.HERE/"verify.py"),
                           "--controls",str(controls),"--receipt",str(receipt)],
                          capture_output=True,text=True,timeout=120)


def test_full_producer_and_independent_cli(artifact):
    assert codec.canonical(build.build()) == (codec.HERE/"controls.json").read_bytes()
    result = cli(codec.HERE/"controls.json",codec.HERE/"receipt.json")
    assert result.returncode == 0,result.stderr
    assert "1280 histories and 78255 native means" in result.stdout


def change(path, value):
    def mutate(packet):
        target = packet
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return mutate


ATTACKS = [
    change(["schema"],"other"),change(["m1_derived"],True),change(["m1_derived"],0),
    change(["source_selected"],True),change(["source","error_clock"],"reset_each_episode"),
    change(["source","raw_input_amplitudes_immutable"],True),
    change(["source","accumulator_retired_at_start"],False),
    change(["source","one_time_unit_seed"],False),change(["source","request_alphabet"],[3]),
    change(["source","schedule"],"A3_selected"),change(["source","grid_Q"],2**512),
    change(["source","rails",0],[0,1]),change(["source","payloads"],[1]),
    change(["source","remote"],"arbitrary_repeated_long_distance_reads"),
    change(["source","logical_input_values_retained"],False),change(["pins"],{}),
    change(["exhaustive"],[]),change(["exhaustive",7,"histories"],5),
    change(["reuse","histories"],1),change(["examples"],[]),
    change(["examples",1,"sessions"],[[3,2]]),
    change(["examples",1,"payload"],2),
    change(["examples",1,"tape",0,3],7),  # off-schedule supported edge
    change(["examples",1,"tape",0,4],45), # unsupported jump
    change(["examples",1,"tape",4,5],0),  # pretend a fresh seed
    change(["examples",1,"tape",8,7],-1), # stale current writer
    change(["examples",1,"tape",8,9],0),  # injected raw reset
    change(["examples",1,"tape",8,2],"host_add"),
    change(["examples",1,"polls",6,3],0), # reset error clock
    change(["examples",1,"polls",6,4],0), # uncharged scale alignment
    change(["examples",1,"polls",6,8],999), # forged output
    change(["examples",1,"polls",6,5],100), # hidden magnitude oracle
    change(["examples",1,"scales",2],0),
    change(["examples",1,"odd_sums"],0),
    change(["examples",1,"final",0,2],-1),
    change(["exhaustive",0,"scalar_means"],0),
    change(["exhaustive",0,"preparation_writes"],0),
    change(["exhaustive",0,"scalar_samples"],0),
    change(["exhaustive",0,"scalar_mean_reads"],0),
    change(["exhaustive",0,"scalar_mean_writes"],0),
    change(["exhaustive",0,"all_histories_sha256"],"0"*64),
    change(["reuse","episodes"],5),change(["reuse","max_scale"],1),
    change(["reuse","scalar_means"],5*878-1),
]


@pytest.mark.parametrize("mutate",ATTACKS)
def test_resealed_semantic_forgeries_fail_real_cli(artifact,tmp_path,mutate):
    packet = deepcopy(artifact)
    mutate(packet)
    assert codec.digest(packet) != codec.digest(artifact),"ineffective mutation"
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt["controls_sha256"] = codec.digest(packet)
    controls, receipt_path = tmp_path/"controls.json", tmp_path/"receipt.json"
    controls.write_bytes(codec.canonical(packet))
    receipt_path.write_bytes(codec.canonical(receipt))
    result = cli(controls,receipt_path)
    assert result.returncode != 0
    assert "verification failed:" in result.stderr


@pytest.mark.parametrize("raw",[
    b"",b"null",b"[]",b"{}",b'{"x":1,"x":2}',b'{"x":NaN}',
    b'{"x":Infinity}',b'{"x":1.0}',b'{"x":1e0}',b'{"x":',b'false',
])
def test_trash_rejected_by_cli(tmp_path,raw):
    path = tmp_path/"trash.json"
    path.write_bytes(raw)
    result = cli(path,codec.HERE/"receipt.json")
    assert result.returncode != 0
    assert "verification failed:" in result.stderr


@pytest.mark.parametrize("field",[
    "histories_replayed","scalar_means_replayed","scalar_samples_replayed",
    "odd_sums_replayed","maximum_scale","preparation_writes_per_history",
    "dynamic_scalar_writes","instrument","error_model","cost_exclusions",
])
def test_forged_receipt_rejected(tmp_path,field):
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt[field] = None
    path = tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(receipt))
    result = cli(codec.HERE/"controls.json",path)
    assert result.returncode != 0
    assert "receipt mismatch" in result.stderr


@pytest.mark.parametrize("payload",range(-2,3))
@pytest.mark.parametrize("sessions",[
    [[]],[[3,2,3,3]],[[2,3],[],[3,3],[2,2]],build.STRESS,
])
def test_signed_preparation_updates_and_samples(sessions,payload,edges):
    ordinary = verify.replay(sessions,payload,edges)
    noisy = verify.replay(sessions,payload,edges,disturbed=True)
    assert [p[-1] for p in noisy["polls"]] == [p[-1] for p in ordinary["polls"]]
    assert any(a[5:7] != b[5:7] for a,b in zip(noisy["tape"],ordinary["tape"]))
    assert [e[:5] for e in noisy["tape"]] == [e[:5] for e in ordinary["tape"]]


def test_data_independent_schedule_and_stable_logical_inputs(edges):
    cases = [verify.replay(build.STRESS,p,edges) for p in range(-2,3)]
    for case in cases:
        assert [e[:5] for e in case["tape"]] == [e[:5] for e in cases[0]["tape"]]
        assert case["scales"] == cases[0]["scales"]
        assert len(case["tape"]) == 878
        assert case["scales"][2] > 1
        assert case["scales"][3] > 1
        assert all(p[-1] == 1 for p in case["polls"] if p[2] == 2)
        assert all(p[-1] == case["payload"] for p in case["polls"] if p[2] == 3)
        assert any(p[0] == 14 and p[3] > 800 for p in case["polls"])


@pytest.mark.parametrize("edge",[(0,1),(9,5),(0,5),(9,1),(7,0),(4,9),
                                (7,1),(4,5),(8,1),(11,5),(1,5),
                                (1,14),(5,40),(14,23),(40,38),(23,45),(38,39)])
def test_each_used_seam_is_required(edge,edges):
    missing = edges-{frozenset(edge)}
    with pytest.raises(ValueError,match="unsupported native mean"):
        verify.replay([[3,2,3],[2,3]],-2,missing)


def test_no_direct_accumulator_clear_edge(edges):
    assert frozenset((0,9)) not in edges
    case = verify.replay([[3],[3]],2,edges)
    assert all(frozenset(e[3:5]) != frozenset((0,9)) for e in case["tape"])
    assert len([e for e in case["tape"] if e[2] == "reset"]) == 8


@pytest.mark.parametrize("z,r,bound,want",[
    (F(0),F(1,4),2,0), (F(-2),F(1,4),2,-2), (F(2),F(1,4),2,2),
    (F(1,2),F(1,2),2,None), (F(1,2),F(1,4),2,None),
    (F(3),F(1,4),2,None), (F(9,4),F(1,4),2,2),
    (F(-9,4),F(1,4),2,-2), (F(1,2),F(1,2)-F(1,100),2,None),
    (F(0),F(1),2,None), (F(0),F(1),0,0),
])
def test_integer_interval_boundaries(monkeypatch,z,r,bound,want):
    # W=0, e=0 yields radius 1/(2G). Choose G to realize exact boundaries.
    q = 3+F(1,2*r)
    monkeypatch.setattr(build,"G",q-3)
    monkeypatch.setattr(verify,"GRID",q)
    plus,minus = (q-3)*z,-(q-3)*z
    assert build.publish(plus,minus,0,0,bound) == want
    assert verify.decode(plus,minus,0,0,bound) == want


@pytest.mark.parametrize("grid",[4,5,19])
@pytest.mark.parametrize("scale",[0,1,3])
@pytest.mark.parametrize("count",[0,1,7])
@pytest.mark.parametrize("bound",[-1,0,1,3])
def test_clipped_decoder_matches_explicit_feasible_set(monkeypatch,grid,scale,count,bound):
    # Enumerate candidates directly; do not use floor/ceiling or the producer.
    # Small grids include wide, overlapping intervals beyond the availability regime.
    monkeypatch.setattr(build,"G",grid-3)
    monkeypatch.setattr(verify,"GRID",grid)
    radius = F((4+5*count)*2**scale,8*(grid-3))
    for quarter in range(-24,25):
        obs = F(quarter,4)
        plus = obs*(grid-3)/2**scale
        minus = -plus
        feasible = [k for k in range(-max(bound,0),max(bound,0)+1)
                    if abs(k)<=bound and abs(obs-k)<=radius]
        expected = feasible[0] if len(feasible) == 1 else None
        assert build.publish(plus,minus,scale,count,bound) == expected
        assert verify.decode(plus,minus,scale,count,bound) == expected


@pytest.mark.parametrize("bad",[[],[[0]],[[1]],[[4]],[[True]],[[3.0]],[["3"]]])
def test_invalid_request_aliases_fail(bad,edges):
    with pytest.raises(ValueError):
        verify.replay(bad,1,edges)
    with pytest.raises(ValueError):
        build.execute(bad,1)


def test_verifier_independence_and_local_decoder():
    source = inspect.getsource(verify)
    assert not re.search(r"(?:from|import)\s+.*\bbuild\b",source)
    assert tuple(inspect.signature(verify.decode).parameters) == (
        "plus","minus","exponent","count","limit")
    assert tuple(inspect.signature(build.publish).parameters) == (
        "p","m","scale","means","bound")
    shared = (codec.HERE/"codec.py").read_text()
    assert "def execute" not in shared and "def replay" not in shared


@pytest.mark.parametrize("grid",[8,16,64,256])
def test_insufficient_fixed_precision_fails(monkeypatch,edges,grid):
    monkeypatch.setattr(verify,"GRID",grid)
    with pytest.raises(ValueError,match="finite precision"):
        verify.replay(build.STRESS,2,edges)


@pytest.mark.parametrize("relative",[
    codec.SUPPORT,"Lean/Geometry/SourceAccumulatorProgram.lean",
    "Lean/Geometry/SourceNativeRecords.lean","docs/AXIOM_REFERENCE.md",
    "code/source_read_acceptance/codec.py",
])
def test_modified_dependency_rejected_before_replay(tmp_path,relative):
    paths = set(codec.pins()) | {"code/source_native_accumulator/controls.json",
                               "code/source_native_accumulator/receipt.json"}
    for path in paths:
        target = tmp_path/path
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(codec.ROOT/path,target)
    with (tmp_path/relative).open("ab") as stream:
        stream.write(b"\n")
    package = tmp_path/"code/source_native_accumulator"
    result = cli(package/"controls.json",package/"receipt.json",package/"verify.py")
    assert result.returncode != 0 and "custody" in result.stderr


def require_ci_gate(source):
    workflow = yaml.safe_load(source)
    step = next(s for s in workflow["jobs"]["build"]["steps"] if s.get("id") == "changed_modules")
    assert "if" not in step
    targets = step["run"].split("targets=(",1)[1].split(")",1)[0]
    assert "Geometry.SourceAccumulatorAxiomAudit" in shlex.split(targets,comments=True)


def test_all_new_theorems_have_transitive_axiom_gate():
    declarations = set()
    for module in ("SourceNativeAccumulator","SourceAccumulatorProgram","SourceAccumulatorDecoder"):
        source = (codec.ROOT/f"Lean/Geometry/{module}.lean").read_text(encoding="utf-8")
        namespace = re.search(r"^namespace (\S+)",source,re.M).group(1)
        declarations.update(namespace+"."+name for name in re.findall(r"^theorem (\w+)",source,re.M))
        code = re.sub(r"/-.*?-/|--[^\n]*","",source,flags=re.S)
        assert not re.search(r"\b(sorry|admit|axiom|native_decide|unsafe)\b",code)
    audit = (codec.ROOT/"Lean/Geometry/SourceAccumulatorAxiomAudit.lean").read_text()
    assert set(re.findall(r"^audit_reusable_bus_axioms (\S+)",audit,re.M)) == declarations
    guard = (codec.ROOT/"Lean/Geometry/SourceReusableBusAxiomAudit.lean").read_text(encoding="utf-8")
    assert "Lean.collectAxioms name" in guard and "sorryAx" in guard and "Lean.ofReduceBool" in guard
    assert "SourceAccumulatorAxiomAudit" in (codec.ROOT/"Lean/Geometry.lean").read_text(encoding="utf-8")
    require_ci_gate((codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("replacement",['','# "Geometry.SourceAccumulatorAxiomAudit"',
                                         '"Geometry.SourceAccumulatorProgram"'])
def test_disabled_ci_trust_gate_rejected(replacement):
    source = (codec.ROOT/".github/workflows/lean-ci.yml").read_text(encoding="utf-8")
    with pytest.raises(AssertionError):
        require_ci_gate(source.replace('"Geometry.SourceAccumulatorAxiomAudit"',replacement))


def test_pinned_sources_have_portable_bytes():
    for path in codec.pins():
        assert b"\r\n" not in (codec.ROOT/path).read_bytes(),path


def test_workflow_covers_pinned_sources():
    path = codec.ROOT/".github/workflows/source-native-accumulator.yml"
    workflow = yaml.safe_load(path.read_text())
    triggers = workflow.get("on",workflow.get(True))
    for event in ("push","pull_request"):
        globs = triggers[event]["paths"]
        for source in codec.pins():
            assert any(fnmatchcase(source,g) for g in globs),source
    job = workflow["jobs"]["controls"]
    assert job["strategy"]["matrix"]["os"] == ["ubuntu-latest"]
    assert any("pytest -q code/source_native_accumulator" in s.get("run","") and "if" not in s
               for s in job["steps"])
