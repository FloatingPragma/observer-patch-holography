"""Hostile semantic inputs, independent custody and complete interventions."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import pytest

from . import check_process, process, verify

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


@pytest.fixture
def events():
    return list(process.execute(2, 2))


def index(events, kind):
    return next(i for i, row in enumerate(events) if row[0] == kind)


@pytest.mark.parametrize("kind", ["prepare", "fork", "flight", "wait", "accumulate", "commit", "checkpoint"])
def test_every_event_kind_is_required(events, kind):
    events.pop(index(events, kind))
    with pytest.raises(ValueError):
        check_process.check(2, 2, None, events)


@pytest.mark.parametrize("kind", ["prepare", "fork", "flight", "wait", "accumulate", "commit", "checkpoint"])
def test_duplicate_events_rejected(events, kind):
    i = index(events, kind)
    events.insert(i, deepcopy(events[i]))
    with pytest.raises(ValueError):
        check_process.check(2, 2, None, events)


@pytest.mark.parametrize("kind", ["prepare", "fork", "flight", "wait", "accumulate", "commit", "checkpoint"])
def test_boolean_indices_are_not_integers(events, kind):
    events[index(events, kind)][1] = False
    with pytest.raises(ValueError, match="integer"):
        check_process.check(2, 2, None, events)


@pytest.mark.parametrize("mutation", ["address", "preparation", "fork_value", "flight_value",
    "flight_duration", "negative_wait", "wait_clock", "wrong_writer", "stale_version",
    "sum", "commit", "checkpoint", "hidden_metadata", "premature_commit", "arrival_order",
    "future_layer", "post_checkpoint", "truncated"])
def test_false_execution_rejected(events, mutation):
    if mutation == "address":
        events[0][2][0] = 9
    elif mutation == "preparation":
        events[0][3] += 1
    elif mutation in {"fork_value", "flight_value"}:
        events[index(events, mutation.split("_")[0])][-1] += 1
    elif mutation == "flight_duration":
        events[index(events, "flight")][4] = "0"
    elif mutation == "negative_wait":
        events[index(events, "wait")][5] = "1"
    elif mutation == "wait_clock":
        events[index(events, "wait")][4] = "1/3"
    elif mutation == "wrong_writer":
        events[index(events, "accumulate")][2] = 7
    elif mutation == "stale_version":
        row = next(r for r in events if r[0] == "accumulate" and r[1] == 1)
        row[4] = row[2]+1
    elif mutation == "sum":
        events[index(events, "accumulate")][-1] += 1
    elif mutation == "commit":
        events[index(events, "commit")][-1] += 1
    elif mutation == "checkpoint":
        events[index(events, "checkpoint")][-1][0] += 1
    elif mutation == "hidden_metadata":
        events.insert(index(events, "accumulate"), ["scheduler_remote_read", 0, 7, 0, 8])
    elif mutation == "premature_commit":
        row = events.pop(index(events, "commit"))
        events.insert(index(events, "accumulate"), row)
    elif mutation == "arrival_order":
        i = index(events, "flight")
        j = max(k for k, row in enumerate(events) if row[0] == "flight" and row[1] == 0)
        events[i], events[j] = events[j], events[i]
        # q=2 admits only one nonzero flight length. Cross the phase boundary
        # instead of pretending that equal-duration arrivals have a strict order.
        i = index(events, "wait")
        events.insert(index(events, "flight"), events.pop(i))
    elif mutation == "future_layer":
        events[index(events, "fork")][1] = 1
    elif mutation == "post_checkpoint":
        events.append(events[0])
    else:
        events.pop()
    with pytest.raises(ValueError):
        check_process.check(2, 2, None, events)


@pytest.mark.parametrize("value", [True, 0.25, "01/4", "2/8", "NaN", "Infinity", "-0", "1"*129])
def test_malformed_clock_rejected(events, value):
    events[index(events, "flight")][4] = value
    with pytest.raises((ValueError, TypeError)):
        check_process.check(2, 2, None, events)


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "bool", "extra", "scope", "layers"])
def test_capture_cannot_filter_or_relabel_attempts(mutation):
    packet = verify.strict_load(HERE/"capture.json")
    if mutation == "missing":
        packet["cases"].pop()
    elif mutation == "duplicate":
        packet["cases"].append(deepcopy(packet["cases"][0]))
    elif mutation == "bool":
        packet["cases"][1]["intervention"] = False
    elif mutation == "extra":
        packet["cases"][0]["passed"] = True
    elif mutation == "scope":
        packet["scope"] = "physical source verified"
    else:
        packet["cases"][0]["layers"] = 1
    with pytest.raises(ValueError):
        verify.verify_capture(packet)


def test_resealed_corrupt_producer_does_not_pass(events, monkeypatch):
    packet = verify.strict_load(HERE/"capture.json")
    events[index(events, "flight")][4] = "0"
    packet["cases"][0]["stream_sha256"] = hashlib.sha256("".join(
        json.dumps(row, separators=(",", ":"))+"\n" for row in events).encode()).hexdigest()
    original = process.execute
    monkeypatch.setattr(process, "execute", lambda q, layers, intervention, stage=-1:
                        iter(events) if (q, intervention) == (2, None) else original(q, layers, intervention, stage))
    with pytest.raises(ValueError, match="unit-flight duration"):
        verify.verify_capture(packet)


def test_independent_checker_imports_no_producer(tmp_path, events):
    path = tmp_path/"events.json"
    path.write_text(json.dumps(events), encoding="utf-8")
    code = '''
import importlib.abc, json, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name in {"source_record_gluing.process", "source_record_gluing.build"}:
            raise RuntimeError("producer import forbidden")
sys.meta_path.insert(0, Block())
from source_record_gluing.check_process import check
with open(sys.argv[1]) as source:
    result = check(2, 2, None, json.load(source))
assert result["events"] == 234 and result["reads_per_layer"] == 32
'''
    result = subprocess.run([sys.executable, "-c", code, str(path)], cwd=ROOT,
                            env=dict(os.environ, PYTHONPATH=str(ROOT/"code")), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr


def test_all_q2_interventions_have_exact_consumed_path_effects():
    baseline = check_process.check(2, 2, None, process.execute(2, 2))
    for i in range(8):
        changed = check_process.check(2, 2, i, process.execute(2, 2, i))
        assert changed["parents"] == baseline["parents"]
        for target in range(8):
            direct = sum(parent == (-1, i) for parent in baseline["parents"][0, target])
            two_hop = sum((-1, i) in baseline["parents"][parent] for parent in baseline["parents"][1, target])
            assert changed["layers"][0][target]-baseline["layers"][0][target] == direct
            assert changed["layers"][1][target]-baseline["layers"][1][target] == two_hop


def test_all_intermediate_versions_are_independently_probed():
    baseline = check_process.check(2, 2, None, process.execute(2, 2))
    for source in range(8):
        result = check_process.check(2, 2, source, process.execute(2, 2, source, 0), 0)
        for target in range(8):
            assert result["layers"][0][target]-baseline["layers"][0][target] == int(source == target)
            assert result["layers"][1][target]-baseline["layers"][1][target] == int(
                (0, source) in baseline["parents"][1, target])


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "value", "wrong_site", "late", "bool_stage"])
def test_commit_intervention_cannot_be_forged(mutation):
    events = list(process.execute(2, 2, 0, 0))
    i = index(events, "intervene")
    stage = 0
    if mutation == "missing":
        events.pop(i)
    elif mutation == "duplicate":
        events.insert(i, deepcopy(events[i]))
    elif mutation == "value":
        events[i][-1] += 1
    elif mutation == "wrong_site":
        events[i][2] = 1
    elif mutation == "late":
        events[i], events[i+1] = events[i+1], events[i]
    else:
        stage = False
    with pytest.raises(ValueError):
        check_process.check(2, 2, 0, events, stage)


def test_initial_input_agreement_does_not_certify_version_dependencies():
    # Both versions hold x+1. Replacing a read of a with b agrees on *all*
    # initial x, but loses the actual a->output dependency. A commit-local
    # intervention at a distinguishes them without altering b or x.
    good = lambda x, da, db: (x+1+da, x+1+db, x+1+da)
    bad = lambda x, da, db: (x+1+da, x+1+db, x+1+db)
    assert all(good(x, 0, 0) == bad(x, 0, 0) for x in range(-20, 21))
    assert good(3, 1, 0) != bad(3, 1, 0)


def test_unequal_flight_arrivals_must_respect_the_reference_clock():
    events = list(process.execute(3, 1))
    i = index(events, "flight")
    j = next(k for k, row in enumerate(events) if row[0] == "flight" and row[4] != events[i][4])
    events[i], events[j] = events[j], events[i]
    with pytest.raises(ValueError, match="arrival chronology"):
        check_process.check(3, 1, None, events)


def test_grid_clock_clearance_and_deleted_rotation_control():
    for q in range(2, 50):
        assert all((2*q+1)*m != 2*q*q for m in range(3*(q-1)**2+1))
    # At speed one, diagonal (3/5,4/5,0) is delivered by one rotated unit
    # flight in one time unit; any axis-only word needs at least 7/5.
    v = (Q(3, 5), Q(4, 5), Q(0))
    assert sum(x*x for x in v) == 1 and sum(abs(x) for x in v) > 1


def test_actual_cli_rejects_false_receipt(tmp_path):
    receipt = verify.strict_load(HERE/"receipt.json")
    receipt["total_reference_events"] -= 1
    path = tmp_path/"forged-receipt.json"
    path.write_text(json.dumps(receipt), encoding="utf-8")
    result = subprocess.run([sys.executable, "-m", "source_record_gluing.verify", "--receipt", str(path)],
                            cwd=ROOT, env=dict(os.environ, PYTHONPATH=str(ROOT/"code")),
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "source custody disagrees with replay" in result.stderr


def test_actual_cli_rejects_changed_scientific_source(tmp_path):
    clone = tmp_path/"source"
    for relative in list(verify.pins())+["code/source_record_gluing/receipt.json"]:
        target = clone/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/relative, target)
    changed = clone/"code/source_selection_model/RECORD_GLUING.md"
    changed.write_bytes(changed.read_bytes()+b"\nTransport law changed.\n")
    result = subprocess.run([sys.executable, "-m", "source_record_gluing.verify"], cwd=clone,
                            env=dict(os.environ, PYTHONPATH=str(clone/"code")), capture_output=True, text=True)
    assert result.returncode != 0
    assert "source custody disagrees with replay" in result.stderr


def test_every_lean_theorem_is_axiom_audited_and_ci_gated():
    source = (ROOT/"Lean/Geometry/SourceRecordGluing.lean").read_text(encoding="utf-8")
    audit = (ROOT/"Lean/Geometry/SourceRecordGluingAxiomAudit.lean").read_text(encoding="utf-8")
    names = re.findall(r"^theorem (\w+)", source, re.M)
    assert len(names) == 20
    assert {"OPH.SourceRecordGluing."+name for name in names} == set(
        re.findall(r"^audit_reusable_bus_axioms (\S+)$", audit, re.M))
    ci = (ROOT/".github/workflows/lean.yml").read_text(encoding="utf-8")
    assert re.search(r'^\s*"Geometry.SourceRecordGluingAxiomAudit"\s*$', ci, re.M)
