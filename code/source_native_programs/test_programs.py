"""Independent native replay, re-sealed CLI attacks and transitive proof gates."""
from copy import deepcopy
from fractions import Fraction as F
from fnmatch import fnmatchcase
from itertools import product
from pathlib import Path
import hashlib
import re
import shlex
import subprocess
import sys
import os

import pytest
import yaml

from . import build,check_routes,codec,routes,verify,family_build,family_verify,lean_control


@pytest.fixture(scope="module")
def packet():
    return codec.load_artifact(codec.HERE/"controls.json")


@pytest.fixture(scope="module")
def allowed():
    return check_routes.edges(3)


def cli(path, receipt):
    env = dict(os.environ,PYTHONPATH=str(codec.ROOT/"code"))
    return subprocess.run([sys.executable,"-m","source_native_programs.verify",
                           "--controls",str(path),"--receipt",str(receipt)],
                          cwd=codec.ROOT,env=env,capture_output=True,text=True,timeout=120)


def test_complete_producer_independent_replay_and_committed_receipt(packet):
    assert codec.canonical(build.controls()) == (codec.HERE/"controls.json").read_bytes()
    result = verify.verify(packet)
    assert result == codec.load_artifact(codec.HERE/"receipt.json")
    assert result["histories"] == 9 and result["scalar_means"] == 16183386


def changed(path,value):
    def alter(packet):
        target = packet
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return alter


ATTACKS = [
    changed(["schema"],"other"),changed(["pins"],{}),changed(["source_selected"],True),
    changed(["m1_derived"],True),changed(["m1_derived"],0),changed(["scope"],"production"),
    changed(["evaluation"],"all_repetitions_expanded"),changed(["program"],[]),
    changed(["program",1,0],[0,1]),changed(["cases"],[]),changed(["routes"],[]),
    changed(["routes",0,"pairs",0],[17,17]),changed(["routes",0,"pairs",1],[1,2]),
    changed(["routes",0,"pairs",0],[True,21]),changed(["routes",0,"pairs",1],[-1,20]),
    changed(["routes",0,"pairs",1],[100000,20]),changed(["routes",0,"direction"],"write"),
    changed(["routes",1,"pairs",0],[9,4]),changed(["routes",0,"carrier"],True),
    changed(["plan","grid_bits"],8),changed(["plan","cleanup_blocks"],1),
    changed(["plan","scalar_means"],1),changed(["plan","maximum_scale"],0),
    changed(["plan","residual_numerator"],0),changed(["plan","cleanup_stages"],1),
    changed(["plan","residual_denominator_exponent"],200),
    changed(["plan","segments",0,"blocks",0,"repeat"],0),
    changed(["plan","segments",0,"blocks",0,"edges"],[[6,2],[7,1]]),
    changed(["plan","segments",0,"blocks",1,"edges"],[[2,1]]),
    changed(["plan","segments",0,"blocks",1,"cleanup"],True),
    changed(["plan","segments",0,"polls",0,"version"],[1,0,-1]),
    changed(["plan","segments",0,"polls",0,"scale"],9),
    changed(["plan","segments",0,"polls",0,"cap"],999),
    changed(["plan","segments",0,"polls",0,"pair"],[6,7]),
    changed(["cases",0,"payload"],[2,2]),changed(["cases",0,"scalar_means"],0),
    changed(["cases",0,"evaluated_means"],0),changed(["cases",0,"stationary_means"],0),
    changed(["cases",0,"final_state_sha256"],"0"*64),
    changed(["cases",0,"history_sha256"],"0"*64),
    changed(["cases",0,"stored_values",0],123),
    changed(["cases",0,"stored_values"],[]),
]


@pytest.mark.parametrize("mutation",ATTACKS)
def test_resealed_forgeries_rejected_by_real_cli(packet,tmp_path,mutation):
    forged = deepcopy(packet)
    mutation(forged)
    assert codec.digest(forged) != codec.digest(packet)
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt["controls_sha256"] = codec.digest(forged)
    path,receipt_path = tmp_path/"controls.json",tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(forged))
    receipt_path.write_bytes(codec.canonical(receipt))
    result = cli(path,receipt_path)
    assert result.returncode != 0, result.stdout


@pytest.mark.parametrize("suffix",[[None],[None,{"ignored":"suffix"}],[False],[{}],[[]],[0]])
def test_resealed_extra_route_records_rejected_by_real_cli(packet,tmp_path,suffix):
    assert packet["pins"] == codec.pins(), "baseline custody must be valid before mutation"
    forged = deepcopy(packet)
    forged["routes"].extend(suffix)
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt["controls_sha256"] = codec.digest(forged)
    path,receipt_path = tmp_path/"controls.json",tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(forged))
    receipt_path.write_bytes(codec.canonical(receipt))
    result = cli(path,receipt_path)
    assert result.returncode != 0 and "extra route witness" in result.stderr


@pytest.mark.parametrize("mutation",["null","hidden_suffix","extra_route","missing","duplicate","swap"])
def test_complete_stream_inventory_rejects_extra_missing_and_reordered_rows(packet,tmp_path,mutation):
    rows = deepcopy(packet["routes"])
    path = tmp_path/"routes.jsonl"
    path.write_bytes("".join(codec.compact(row)+"\n" for row in rows).encode("ascii"))
    baseline = check_routes.check(check_routes.read_rows(path),2,3)
    assert baseline["routes"] == 8
    if mutation == "null":
        rows.append(None)
    elif mutation == "hidden_suffix":
        rows.extend([None,{"ignored":"suffix"}])
    elif mutation == "extra_route":
        rows.append(deepcopy(rows[-1]))
    elif mutation == "missing":
        rows.pop()
    elif mutation == "duplicate":
        rows[-1] = deepcopy(rows[-2])
    else:
        rows[-1],rows[-2] = rows[-2],rows[-1]
    path.write_bytes("".join(codec.compact(row)+"\n" for row in rows).encode("ascii"))
    with pytest.raises(ValueError):
        check_routes.check(check_routes.read_rows(path),2,3)


def test_family_cli_rejects_null_followed_by_hidden_route_suffix(tmp_path):
    family = codec.load_artifact(codec.HERE/"families.json")
    assert family["pins"] == codec.pins(), "baseline family custody must be valid"
    rows = list(routes.generate(27,3))
    assert check_routes.check(rows,27,3) == family["routes"][0]
    path = tmp_path/"routes-27.jsonl"
    path.write_bytes("".join(codec.compact(row)+"\n" for row in rows+[None,{}]).encode("ascii"))
    env = dict(os.environ,PYTHONPATH=str(codec.ROOT/"code"))
    result = subprocess.run([sys.executable,"-m","source_native_programs.family_verify",
                             "--workdir",str(tmp_path)],cwd=codec.ROOT,env=env,
                            capture_output=True,text=True,timeout=120)
    # It must reject the inventory itself, before asking for an execution plan
    # or later witnesses; an unrelated missing-file error would mask the bug.
    assert result.returncode != 0 and "extra route witness" in result.stderr


@pytest.mark.parametrize("raw",[b"",b"null",b"[]",b"{}",b"true",b'{"x":1,"x":2}',
                              b'{"x":1.0}',b'{"x":NaN}',b'{"x":Infinity}',b'{"x":1e0}'])
def test_trash_rejected_by_real_cli(tmp_path,raw):
    path = tmp_path/"trash.json"
    path.write_bytes(raw)
    assert cli(path,codec.HERE/"receipt.json").returncode != 0


@pytest.mark.parametrize("field",["histories","scalar_means","stationary_means","scalar_samples",
                                  "grid_bits","preparation_writes_per_history","large_family_scope"])
def test_resealed_receipt_fields_rejected(tmp_path,field):
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt[field] = None
    path = tmp_path/"receipt.json"
    path.write_bytes(codec.canonical(receipt))
    result = cli(codec.HERE/"controls.json",path)
    assert result.returncode != 0 and "receipt mismatch" in result.stderr


def test_verifier_does_not_call_producer(packet,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("producer called by verifier")
    monkeypatch.setattr(build,"compile_program",forbidden)
    monkeypatch.setattr(build,"execute",forbidden)
    monkeypatch.setattr(routes,"generate",forbidden)
    verify.verify(packet)


def test_stationary_evaluation_matches_full_expansion(packet,allowed):
    plan = packet["plan"]
    fast = verify.replay(plan,[-2,2],allowed)
    full = verify.replay(plan,[-2,2],allowed,accelerate=False)
    assert full["stationary_means"] == 0
    assert full["evaluated_means"] == full["scalar_means"]
    for key in ("checkpoints","block_states_sha256","final_state_sha256","scalar_means"):
        assert fast[key] == full[key]
    assert fast["stationary_means"] > 0


@pytest.mark.parametrize("mode",["preparation_error","disturbed"])
def test_signed_preparation_idle_port_disturbance_and_samples(packet,allowed,mode):
    result = verify.replay(packet["plan"],[-2,2],allowed,**{mode:True})
    assert result["published"] == [0,2]
    assert result["scalar_means"] == packet["plan"]["scalar_means"]


def test_banks_retain_versions_and_publish_computed_inputs(packet,allowed):
    result = verify.replay(packet["plan"],[-2,2],allowed)
    committed = {}
    for checkpoint in result["checkpoints"]:
        for sample in checkpoint["observations"]:
            key = tuple(sample["version"])
            if key in committed:
                assert committed[key] == sample["decoded"]
            committed[key] = sample["decoded"]
    assert [committed[t,i,-1] for t in (1,2,3) for i in (0,1)] == [1,-3,-1,1,0,2]


def test_actual_opposite_operand_polarity_and_reversed_bank_orientation(packet,allowed):
    rows = deepcopy(packet["routes"])
    adjacency = routes.graph(3)
    protected = {12*c+r for c in range(1,5) for r in (0,1)}
    found = 0
    for row in rows:
        c = row["carrier"]
        if row["direction"] == "read":
            row["pairs"] = routes.paired(adjacency,(12*c+5,12*c+9),(5,3),protected|{1,2,4,6,7,9})
        else:
            row["pairs"] = routes.paired(adjacency,(4,9),(12*c+1,12*c),
                                         (protected-{12*c,12*c+1})|{1,2,3,5,6,7})
        found += row["pairs"][-1][0] > row["pairs"][-1][1]
    assert found > 0
    _,mapping = check_routes.check(rows,2,3,True)
    plan = build.compile_program(build.LAYERS,rows)
    verify.certify_plan(plan,build.LAYERS,mapping)
    codec.equal(build.execute(plan,[-2,2]),verify.replay(plan,[-2,2],allowed),"opposite polarity replay")


@pytest.mark.parametrize("program",[[[[]]],[[[0]],[[0,0]],[[0]],[[0,0,0]]],
                                     [[[0],[0]],[[1],[0,1]],[[0,1,0],[]],[[1],[1]]]])
def test_variable_finite_programs(program,allowed):
    n = len(program[0])
    rows = list(routes.generate(n,3))
    _,mapping = check_routes.check(rows,n,3,True)
    for bound in (0,1,7):
        plan = build.compile_program(program,rows,bound)
        verify.certify_plan(plan,program,mapping)
        payload = [-bound]*n
        codec.equal(build.execute(plan,payload),verify.replay(plan,payload,allowed),"generic program")


def test_interval_comparator_matches_complete_candidate_set():
    for plus,minus,shift,cap,count in product(range(-3,4),range(-3,4),range(3),range(3),range(3)):
        result,center,radius = verify.local_decode(plus,minus,shift,cap,count,16,F(1,64))
        candidates = [i for i in range(-cap,cap+1) if abs(center-i) <= radius]
        assert result == (candidates[0] if len(candidates) == 1 else None)


@pytest.mark.parametrize("q",[3,13])
def test_complete_metric_is_independent_and_matches_existing_custody(q):
    assert family_build.menu(q) == family_verify.metric_menu(q)
    if q == 3:
        assert sum(map(len,family_verify.metric_menu(q))) == 311


def require_ci_gate(source):
    workflow = yaml.safe_load(source)
    job = workflow["jobs"]["build"]
    assert "if" not in job and not job.get("continue-on-error",False)
    step = next(s for s in job["steps"] if s.get("id") == "changed_modules")
    assert "if" not in step and not step.get("continue-on-error",False)
    targets = step["run"].split("targets=(",1)[1].split(")",1)[0]
    assert "Geometry.SourceNativeProgramsAxiomAudit" in shlex.split(targets,comments=True)


def test_complete_transitive_axiom_inventory_and_ci_gate():
    declarations = set()
    for module in ("SourceNativeShuttle","SourceNativeProgramError","SourceNativeCore",
                   "SourceNativeStoredProgram","SourceNativeProgramBudget",
                   *(path.stem for path in sorted((codec.ROOT/"Lean/Geometry").glob("SourceBank*.lean")))):
        source = (codec.ROOT/f"Lean/Geometry/{module}.lean").read_text(encoding="utf-8")
        code = re.sub(r"/-.*?-/|--[^\n]*","",source,flags=re.S)
        namespace = re.search(r"^namespace (\S+)",code,re.M).group(1)
        declarations.update(namespace+"."+name for name in re.findall(r"^(?:theorem|lemma) (\w+)",code,re.M))
        assert not re.search(r"\b(sorry|admit|axiom|native_decide|unsafe)\b",code)
    audit = (codec.ROOT/"Lean/Geometry/SourceNativeProgramsAxiomAudit.lean").read_text(encoding="utf-8")
    assert set(re.findall(r"^audit_reusable_bus_axioms (\S+)",audit,re.M)) == declarations
    assert "SourceNativeProgramsAxiomAudit" in (codec.ROOT/"Lean/Geometry.lean").read_text(encoding="utf-8")
    require_ci_gate((codec.ROOT/".github/workflows/lean.yml").read_text(encoding="utf-8"))


def test_kernel_control_certificate_is_complete_and_current(packet,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("producer called by kernel control exporter")
    monkeypatch.setattr(build,"compile_program",forbidden)
    monkeypatch.setattr(build,"execute",forbidden)
    monkeypatch.setattr(routes,"generate",forbidden)
    assert lean_control.TARGET.read_bytes() == lean_control.render(packet)
    assert len(packet["plan"]["segments"]) == 42


@pytest.mark.parametrize("mutation",[
    changed(["program"],[]),changed(["routes"],[]),changed(["plan","segments"],[]),
    changed(["plan","segments",0,"blocks",0,"repeat"],1),
    changed(["plan","segments",1,"polls",-1,"scale"],0),
    changed(["plan","segments",1,"polls",-1,"version"],[0,1,-1]),
    changed(["plan","segments",-1,"kind"],"read"),
    changed(["plan","scalar_means"],0),changed(["cases"],[]),
    changed(["cases",0,"stored_values"],[]),
    changed(["cases",0,"stored_values",-1],True),
    changed(["cases",0,"stored_values",-1],12345),
])
def test_kernel_control_cli_rejects_forged_or_stale_correspondence(packet,tmp_path,mutation):
    # The intact packet must reach and pass the real CLI before an attack.
    env = dict(os.environ,PYTHONPATH=str(codec.ROOT/"code"))
    path = tmp_path/"controls.json"
    command = [sys.executable,"-m","source_native_programs.lean_control",
               "--controls",str(path),"--check"]
    path.write_bytes(codec.canonical(packet))
    good = subprocess.run(command,cwd=codec.ROOT,env=env,capture_output=True,text=True,timeout=120)
    assert good.returncode == 0,good.stderr
    forged = deepcopy(packet)
    mutation(forged)
    path.write_bytes(codec.canonical(forged))
    bad = subprocess.run(command,cwd=codec.ROOT,env=env,capture_output=True,text=True,timeout=120)
    assert bad.returncode != 0,bad.stdout


@pytest.mark.parametrize("mutation,reason",[
    (changed(["m1_derived"],True),"unproved M1 promotion"),
    (changed(["pins"],{}),"source/proof custody"),
    (changed(["cases",0,"stored_values",0],12345),"native state, version, observation or work differs"),
    (changed(["cases",0,"history_sha256"],"0"*64),"native state, version, observation or work differs"),
])
def test_kernel_control_cli_authenticates_evidence_outside_projection(packet,tmp_path,mutation,reason):
    forged = deepcopy(packet)
    mutation(forged)
    # These fields do not alter any represented theorem. Projection freshness
    # alone therefore cannot authenticate the supplied native evidence.
    assert lean_control.render(forged) == lean_control.render(packet)
    env = dict(os.environ,PYTHONPATH=str(codec.ROOT/"code"))
    path = tmp_path/"controls.json"
    command = [sys.executable,"-m","source_native_programs.lean_control",
               "--controls",str(path),"--check"]
    path.write_bytes(codec.canonical(packet))
    good = subprocess.run(command,cwd=codec.ROOT,env=env,capture_output=True,text=True,timeout=120)
    assert good.returncode == 0,good.stderr
    path.write_bytes(codec.canonical(forged))
    bad = subprocess.run(command,cwd=codec.ROOT,env=env,capture_output=True,text=True,timeout=120)
    assert bad.returncode != 0 and reason in bad.stderr,bad.stdout+bad.stderr


@pytest.mark.parametrize("replacement",['','# "Geometry.SourceNativeProgramsAxiomAudit"',
                                         '"Geometry.SourceNativeStoredProgram"'])
def test_disabled_ci_axiom_gate_is_rejected(replacement):
    source = (codec.ROOT/".github/workflows/lean.yml").read_text(encoding="utf-8")
    with pytest.raises(AssertionError):
        require_ci_gate(source.replace('"Geometry.SourceNativeProgramsAxiomAudit"',replacement))


def require_native_ci_gate(workflow):
    triggers = workflow.get("on",workflow.get(True))
    for event in ("push","pull_request"):
        for path in codec.pins():
            assert any(fnmatchcase(path,g) for g in triggers[event]["paths"]),path
    for name,commands in {
        "controls":["python -m pytest -q code/source_native_programs"],
        "families":[
            "python -m source_native_programs.family_build --workdir temp/native-program-witnesses --output temp/native-program-families.json",
            "python -m source_native_programs.family_verify --workdir temp/native-program-witnesses --families temp/native-program-families.json",
        ],
    }.items():
        job = workflow["jobs"][name]
        assert "if" not in job and not job.get("continue-on-error",False)
        for command in commands:
            assert any(step.get("run","").strip() == command and "if" not in step and
                       not step.get("continue-on-error",False) for step in job["steps"])


@pytest.mark.parametrize("key",["scalar_means_upper","grid_bits_sufficient","maximum_scale_upper",
                                "complete_metric_reads","scope","input_bound"])
def test_resealed_family_bounds_rejected(packet,key):
    menu = family_verify.metric_menu(3)
    rows = list(routes.generate(27,3))
    summary = check_routes.check(rows,27,3)
    bound = family_build.resource_bound(3,menu,summary["max_paired_cells"])
    family_verify.check_bound(bound,3,menu,summary)
    bound[key] = 0
    with pytest.raises(ValueError,match="resource bound"):
        family_verify.check_bound(bound,3,menu,summary)


@pytest.mark.parametrize("menu",[
    [[] for _ in range(27)],
    [[0]*27 for _ in range(27)],
    [[(17*i+7*j)%27 for j in range(i%8)] for i in range(27)],
],ids=["empty","concentrated_duplicate_reads","uneven_fanout"])
def test_family_bound_covers_actual_compilation_and_full_cleanup_word(menu):
    witnesses = list(routes.generate(27,3))
    summary,mapping = check_routes.check(witnesses,27,3,True)
    plan = build.compile_program([menu,menu],witnesses,28)
    verify.certify_plan(plan,[menu,menu],mapping)
    bound = family_build.resource_bound(3,menu,summary["max_paired_cells"])
    family_verify.check_bound(bound,3,menu,summary)
    assert plan["maximum_scale"] <= bound["maximum_scale_upper"]
    assert plan["cleanup_blocks"] <= bound["cleanup_blocks_sufficient"]
    # Charge the actual emitted blocks at the bound's larger cleanup count.
    # Comparing only the producer's smaller word would miss an unpaid-cleanup
    # term in the analytic resource formula.
    k = bound["cleanup_blocks_sufficient"]
    steps = 0
    for stage in plan["segments"]:
        for block in stage["blocks"]:
            repeat = block["repeat"]
            if block["cleanup"]:
                assert repeat % plan["cleanup_blocks"] == 0
                repeat = repeat // plan["cleanup_blocks"] * k
            steps += len(block["edges"]) * repeat
    assert steps <= bound["scalar_means_upper"]
    grid = 1 << bound["grid_bits_sufficient"]
    residual = F(28*plan["cleanup_stages"]*plan["maximum_bus_cells"]**2,1 << k)
    _,_,radius = verify.local_decode(0,0,plan["maximum_scale"],0,steps,grid,residual)
    assert 2*radius < 1


@pytest.mark.parametrize("field,value",[("source_selected",True),("m1_derived",True),
                                       ("production_native_replay",True),("q3_native_replay",False),
                                       ("routes",[]),("bounds",[]),("pins",{}),
                                       ("placement","Morton_addresses_preserved")])
def test_family_scope_and_inventory_forgeries_fail_before_replay(tmp_path,field,value):
    family = codec.load_artifact(codec.HERE/"families.json")
    assert family["pins"] == codec.pins(), "baseline family custody must be valid before mutation"
    family[field] = value
    with pytest.raises(ValueError):
        family_verify.verify(family,tmp_path)


@pytest.mark.parametrize("mutation",["crlf","extra_space","duplicate_key","float_port"])
def test_route_witness_bytes_are_strict(packet,tmp_path,mutation):
    row = codec.compact(packet["routes"][0])+"\n"
    if mutation == "crlf":
        row = row.replace("\n","\r\n")
    elif mutation == "extra_space":
        row = " "+row
    elif mutation == "duplicate_key":
        row = row.replace('"carrier":1,','"carrier":1,"carrier":1,')
    else:
        row = row.replace('"carrier":1,','"carrier":1.0,')
    path = tmp_path/"routes.jsonl"
    path.write_bytes(row.encode("ascii"))
    with pytest.raises(ValueError):
        list(check_routes.read_rows(path))
