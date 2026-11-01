"""Full replay, observation ambiguity, stopping accounting and resealed attacks."""
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

from . import build,codec,verify


@pytest.fixture(scope="module")
def artifact():
    return codec.load_artifact(codec.HERE/"controls.json")


def cli(controls,receipt,script=None):
    return subprocess.run([sys.executable,str(script or codec.HERE/"verify.py"),
                           "--controls",str(controls),"--receipt",str(receipt)],
                          capture_output=True,text=True,timeout=600)


def test_full_reproduction_and_real_cli(artifact):
    codec.equal(build.build(),artifact,"producer reproduction")
    result = cli(codec.HERE/"controls.json",codec.HERE/"receipt.json")
    assert result.returncode==0,result.stderr
    assert "174762 histories and 2642144 native means" in result.stdout


def setitem(path,value):
    def change(packet):
        target = packet
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return change


MUTATIONS = [
    setitem(["schema"],"other"), setitem(["m1_derived"],True),
    setitem(["full_axiom_source_contract"],True), setitem(["m1_derived"],0),
    setitem(["source","source_selection"],"A3_derived"),
    setitem(["source","reference"],"conditioned_on_success"),
    setitem(["source","cover"],"receiver_only"),
    setitem(["source","feasible"],"successful_histories_only"),
    setitem(["source","rails",4,0],44),setitem(["source","grid_Q"],32768),
    setitem(["source","payloads"],[-1]),setitem(["source","macro_alphabet"],[[0,1]]),
    setitem(["exhaustive"],[]), setitem(["examples"],[]),setitem(["path_laws"],[]),
    setitem(["exhaustive",8,"aborts_per_payload"],0),
    setitem(["exhaustive",8,"unconditional_expected_means"],"8"),
    setitem(["exhaustive",8,"scalar_means_both_payloads"],0),
    setitem(["exhaustive",8,"words_per_payload"],0),
    setitem(["exhaustive",0,"all_stopped_histories_sha256"],"0"*64),
    setitem(["examples",0,"result"],1),setitem(["examples",0,"steps"],0),
    setitem(["examples",0,"polls",0,3],-1),
    setitem(["examples",0,"word"],[3,2,1,0]),
    setitem(["path_laws",8,"abort_probability"],"0"),
    setitem(["path_laws",8,"expected_scalar_means"],"1"),
    setitem(["path_laws",8,"sufficient_grid_Q"],"1"),
    lambda p: p.update(unchecked_success=True),
]


@pytest.mark.parametrize("mutation",MUTATIONS)
def test_resealed_mutations_fail_real_cli(artifact,mutation,tmp_path):
    altered = deepcopy(artifact)
    mutation(altered)
    receipt = codec.load_artifact(codec.HERE/"receipt.json")
    receipt["controls_sha256"] = codec.digest(altered)
    controls_path,receipt_path = tmp_path/"controls.json",tmp_path/"receipt.json"
    controls_path.write_bytes(codec.canonical(altered))
    receipt_path.write_bytes(codec.canonical(receipt))
    with pytest.raises((ValueError,KeyError,TypeError)):
        verify.verify(altered)
    result = cli(controls_path,receipt_path)
    assert result.returncode!=0
    assert "verification failed" in result.stderr


@pytest.mark.parametrize("field",range(7))
def test_consumed_native_event_forgery(artifact,field):
    altered = deepcopy(artifact)
    altered["examples"][0]["tape"][0][field] += 1
    with pytest.raises(ValueError,match="native tape"):
        verify.verify(altered)


@pytest.mark.parametrize("raw",["null\n","[]\n","{}\n",'{"a":1,"a":2}\n',
                                     '{"x":1.0}\n','{"x":NaN}\n','{"x":Infinity}\n',"true\n"])
def test_trash_rejected(raw,tmp_path):
    path = tmp_path/"controls.json"
    path.write_text(raw,encoding="ascii")
    result = cli(path,codec.HERE/"receipt.json")
    assert result.returncode!=0 and "verification failed" in result.stderr


def test_interval_policy_is_exact_singleton_rule():
    # Continuous uncertainty models, sampled at every endpoint and on both sides.
    # This independently checks candidate sets, not a copy of the decision branches.
    for t in range(9):
        E,A = build.budget(t),F(build.AMPLITUDE,build.Q)
        endpoints = [-A-E,-E,E,A+E]
        values = {F(0),*endpoints}
        for x in endpoints:
            values.update([x-F(1,100000),x+F(1,100000)])
        for z in values:
            candidates = []
            if -E<=z<=A+E:
                candidates.append(1)
            if -A-E<=z<=E:
                candidates.append(-1)
            expected = candidates[0] if len(candidates)==1 else 0
            assert build.publish(2+z,2-z,t)==expected
            assert verify.local_decision((2+z)*build.Q,(2-z)*build.Q,t)==expected
    assert list(inspect.signature(build.publish).parameters)==[
        "positive_sample","negative_sample","checkpoint"]


def test_policies_on_ambiguous_and_impossible_observations():
    possible = {"blank":{-1,1},"negative":{-1},"positive":{1},"impossible":set()}
    observations = list(possible)
    sound = []
    for outputs in product((-1,0,1),repeat=4):
        policy = dict(zip(observations,outputs))
        if all(result==0 or all(v==result for v in possible[o]) for o,result in policy.items()):
            sound.append(policy)
    assert all(policy["blank"]==0 for policy in sound)
    assert any(policy["impossible"]==1 for policy in sound)  # soundness alone is vacuous here
    canonical = {o:next(iter(values)) if len(values)==1 else 0 for o,values in possible.items()}
    assert canonical=={"blank":0,"negative":-1,"positive":1,"impossible":0}


def test_actual_rounding_and_no_payload_timing_leak(artifact):
    seen_odd = False
    for minus,plus in zip(artifact["examples"][::2],artifact["examples"][1::2]):
        assert minus["steps"]==plus["steps"]
        assert minus["result"]==-plus["result"]
        assert [row[:2] for row in minus["tape"]]==[row[:2] for row in plus["tape"]]
        for row in minus["tape"]+plus["tape"]:
            if (row[2]+row[3])%2:
                seen_odd = True
                assert abs(F(row[6])-F(row[2]+row[3],2))==F(1,2)
    assert seen_odd


def test_uniform_marginals_do_not_supply_independence():
    # Each time marginal is uniform under this four-word correlated law.
    # None of its words delivers: the IID arrival formula cannot be reused.
    for edge in range(4):
        for payload in (-1,1):
            assert verify.replay((edge,)*8,payload)["result"]==0
    assert verify.arrival_law(4,8)[1]<4**8


def test_suffixes_are_not_executed_but_their_mass_is_retained(artifact):
    # Successful first prefix 0,1,2,3 has 4^4 continuations at H=8.
    assert artifact["exhaustive"][8]["first_publication_counts"][4]==4**4
    case = verify.replay((0,1,2,3,3,3,3,3),1)
    assert case["steps"]==4 and len(case["tape"])==8
    assert len(case["polls"])==5 and case["result"]==1
    failure = verify.replay((0,)*8,1)
    assert failure["steps"]==8 and len(failure["tape"])==16 and failure["result"]==0


@pytest.mark.parametrize("d",range(1,9))
@pytest.mark.parametrize("h",range(21))
def test_general_counting_and_precision(d,h):
    codec.equal(build.counting(d,h),verify.path_law(d,h),"binomial/first-hit recurrence")


def test_finite_precision_cannot_be_extrapolated():
    # The fixed experimental error budget overtakes the sufficient signal bound.
    assert 2*build.budget(8)<F(build.AMPLITUDE,build.Q)/2**8
    assert 2*build.budget(16)>F(build.AMPLITUDE,build.Q)/2**16


@pytest.mark.parametrize("relative",[
    codec.SUPPORT,"Lean/Geometry/SourceReadAcceptance.lean",
    "Lean/ObserverPatchHolography/RepairWordSchedule.lean","docs/AXIOM_REFERENCE.md"])
def test_changed_source_fails_before_replay(relative,tmp_path):
    paths = set(codec.pins()) | {"code/source_read_acceptance/controls.json",
                               "code/source_read_acceptance/receipt.json"}
    for path in paths:
        destination = tmp_path/path
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(codec.ROOT/path,destination)
    with (tmp_path/relative).open("ab") as stream:
        stream.write(b"\n")
    package = tmp_path/"code/source_read_acceptance"
    result = cli(package/"controls.json",package/"receipt.json",package/"verify.py")
    assert result.returncode!=0 and "custody" in result.stderr


def require_ci_audit(text):
    workflow = yaml.safe_load(text)
    step = next(s for s in workflow["jobs"]["build"]["steps"] if s.get("id")=="changed_modules")
    targets = step["run"].split("targets=(",1)[1].split(")",1)[0]
    assert "Geometry.SourceReadAcceptanceAxiomAudit" in shlex.split(targets,comments=True)


def test_transitive_proof_audit_and_ci():
    audit = (codec.ROOT/"Lean/Geometry/SourceReadAcceptanceAxiomAudit.lean").read_text(encoding="utf-8")
    names = set()
    for module in ("SourceReadAcceptance","SourceReadAcceptanceSchedule"):
        proof = (codec.ROOT/f"Lean/Geometry/{module}.lean").read_text(encoding="utf-8")
        names.update("OPH."+module+"."+n for n in re.findall(r"^theorem (\w+)",proof,re.M))
    assert names==set(re.findall(r"^audit_reusable_bus_axioms ([\w.]+)$",audit,re.M))
    assert "Lean/Geometry/SourceReusableBusAxiomAudit.lean" in codec.proof_paths()
    require_ci_audit((codec.ROOT/".github/workflows/lean-ci.yml").read_text())
    workflow = yaml.safe_load((codec.ROOT/".github/workflows/source-read-acceptance.yml").read_text())
    assert workflow["jobs"]["controls"]["strategy"]["matrix"]["os"]==["ubuntu-latest"]
    assert any(s.get("run")=="python -m pytest -q code/source_read_acceptance"
               for s in workflow["jobs"]["controls"]["steps"])
    triggers = workflow.get("on",workflow.get(True))
    for event in ("push","pull_request"):
        for path in codec.pins():
            assert any(fnmatchcase(path,pattern) for pattern in triggers[event]["paths"]),path


@pytest.mark.parametrize("replacement",['','# "Geometry.SourceReadAcceptanceAxiomAudit"',
                                         '"Geometry.SourceReadAcceptance"'])
def test_disabled_audit_gate_is_rejected(replacement):
    source = (codec.ROOT/".github/workflows/lean-ci.yml").read_text()
    with pytest.raises(AssertionError):
        require_ci_audit(source.replace('"Geometry.SourceReadAcceptanceAxiomAudit"',replacement))
