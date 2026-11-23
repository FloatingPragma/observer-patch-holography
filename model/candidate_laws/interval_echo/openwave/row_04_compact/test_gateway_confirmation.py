from copy import deepcopy
from base64 import b64decode
from hashlib import sha256
import json
import subprocess
import sys

import pytest

from . import gateway_confirmation as held
from .gateway_release import ARMS, verify


def test_owner_can_reproduce_the_exact_registration_bytes():
    task = held.HERE.parent/"submission"/"M10_4_TASK.txt"
    proposed = task.read_bytes().replace(b"\r\n", b"\n")
    assert sha256(proposed).hexdigest() == held.TASK_SHA256
    assert held.PREREG_SHA256.encode() in proposed
    assert sha256(held.PROTOCOL.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == held.PREREG_SHA256


def test_offline_first_parent_registration_selects_only_the_earliest_copy(tmp_path, monkeypatch):
    def git(*args):
        return subprocess.run(["git", "-C", str(tmp_path), *args], check=True,
                              capture_output=True, text=True).stdout.strip()

    git("init", "--quiet")
    git("config", "user.name", "Test")
    git("config", "user.email", "test@example.invalid")
    git("commit", "--allow-empty", "-m", "baseline")
    base = git("rev-parse", "HEAD")
    protocol = tmp_path/"local_prereg.json"
    protocol.write_bytes(json.dumps({"prior_platform_revision": base}).encode())
    digest = sha256(protocol.read_bytes()).hexdigest()
    task = tmp_path/held.TASK
    task.parent.mkdir(parents=True)
    active = ("registered "+digest+"\n").encode()
    monkeypatch.setattr(held, "PROTOCOL", protocol)
    monkeypatch.setattr(held, "PREREG_SHA256", digest)
    monkeypatch.setattr(held, "TASK_SHA256", sha256(active).hexdigest())

    def publish(data):
        task.write_bytes(data)
        git("add", "--", held.TASK)
        git("commit", "-qm", "task version")
        head = git("rev-parse", "HEAD")
        git("update-ref", "refs/remotes/origin/main", head)
        return head

    prior = publish(b"proposal\n")
    assert held.registration_anchor(tmp_path) is None
    with pytest.raises(ValueError, match="absent"):
        held.verify_anchor(tmp_path, prior)
    git("switch", "-c", "approved-task")
    side = publish(active)
    git("switch", "-c", "platform", prior)
    git("merge", "--no-ff", "approved-task", "-m", "register reviewed task")
    first = git("rev-parse", "HEAD")
    git("update-ref", "refs/remotes/origin/main", first)
    publish(b"revision\n")
    duplicate = publish(active)
    assert held.registration_anchor(tmp_path) == first
    held.verify_anchor(tmp_path, first)
    with pytest.raises(ValueError, match="first registered"):
        held.verify_anchor(tmp_path, side)
    with pytest.raises(ValueError, match="first registered"):
        held.verify_anchor(tmp_path, duplicate)
    protocol.write_bytes(b"{}")
    with pytest.raises(ValueError, match="gates"):
        held.verify_anchor(tmp_path, first)


def test_held_verifier_rejects_changed_control_failure_and_bin_permutation():
    # Synthetic validator fixture. This executes no registered acquisition.
    settings = dict(home=0, start=0, direction=1, exposure=15)
    def distribution(p):
        return [p, *([0.0]*11), 1-p, 0.0]
    bins = {arm: distribution(p) for arm, p in zip(ARMS, (.99, .6, .3, .2), strict=True)}
    peer = dict(anchor="synthetic", settings=settings, bins=bins)
    actual = dict(peer, floor_ppt=990_000_000_000)
    held.check(actual, peer)
    for arm, p in (("coherent", .93), ("dephased", .98), ("stochastic_swap", .85)):
        both = deepcopy(peer)
        both["bins"][arm] = distribution(p)
        with pytest.raises(ArithmeticError, match="scientific gates"):
            held.check(dict(both, floor_ppt=actual["floor_ppt"]), both)
    with pytest.raises(ArithmeticError, match="scientific gates"):
        held.check(dict(actual, floor_ppt=940_000_000_000), peer)
    for floor in (True, 990_000_000_000., float("nan"), -1, held.SCORE+1):
        with pytest.raises(ValueError, match="integer enclosure"):
            held.check(dict(actual, floor_ppt=floor), peer)
    invalid = deepcopy(peer)
    invalid["bins"]["coherent"][0] -= .01
    with pytest.raises(ValueError, match="independent terminal law"):
        held.check(actual, invalid)
    with pytest.raises(ArithmeticError, match="normalization"):
        held.check(dict(invalid, floor_ppt=actual["floor_ppt"]), peer)
    for value in (float("nan"), True, -.001):
        false = deepcopy(peer)
        false["bins"]["coherent"][0] = value
        with pytest.raises(ValueError):
            held.check(actual, false)
    false = deepcopy(peer)
    false["bins"]["coherent"].reverse()
    with pytest.raises(ArithmeticError, match="differs"):
        held.check(actual, false)
    for garbage in ({"a": [False]}, {"a": [float("nan")]}, {"a": []}, {"wrong": [0.5]}):
        with pytest.raises(ValueError):
            verify({"a": [0.5]}, garbage)


@pytest.mark.parametrize("anchor", ["", "a"*39, "A"*40, "0"*40+" ", False])
def test_invalid_anchors(anchor):
    with pytest.raises(ValueError):
        held.selection(anchor)


def test_held_gate_failure_retains_original_evidence(tmp_path, monkeypatch):
    # No flight is executed in this synthetic file-custody test.
    anchor = "0"*40
    settings = held.selection(anchor)
    bins = {arm: [.5, *([0.0]*11), .5, 0.0] for arm in ARMS}
    peer = {"anchor": anchor, "settings": settings, "bins": bins}
    reference = tmp_path/"independent.json"
    reference.write_text(json.dumps(peer), encoding="utf-8")
    receipt = tmp_path/"held.json"
    monkeypatch.setattr(sys, "argv", ["gateway_confirmation", anchor, "--openwave-checkout", str(tmp_path),
                                      "--compare", str(reference), "--record", str(receipt)])
    monkeypatch.setattr(held, "verify_anchor", lambda *args: None)
    monkeypatch.setattr(held, "acquisition", lambda _: dict(peer, floor_ppt=500_000_000_000))
    with pytest.raises(ArithmeticError, match="scientific gates"):
        held.main()
    original = receipt.read_bytes()
    result = json.loads(original)
    assert result["passed"] is False and result["acquisition"]["bins"] == bins
    assert result["status"] == "rejected"
    assert result["comparison"] == {"absolute": 0., "arm": "coherent", "outcome": 0}
    assert result["independent_sha256"] == sha256(reference.read_bytes()).hexdigest()
    assert b64decode(result["independent_base64"], validate=True) == reference.read_bytes()
    with pytest.raises(FileExistsError):
        held.main()
    assert receipt.read_bytes() == original


@pytest.mark.parametrize("data", [None, b"\xff", b'{"bins":{},"bins":{}}', b'{"p":NaN}', b"[]"])
def test_incomplete_independent_calculation_retains_evidence(tmp_path, monkeypatch, data):
    anchor = "0"*40
    reference, receipt = tmp_path/"peer", tmp_path/"receipt"
    if data is not None:
        reference.write_bytes(data)
    monkeypatch.setattr(sys, "argv", ["gate", anchor, "--openwave-checkout", str(tmp_path),
                                      "--compare", str(reference), "--record", str(receipt)])
    monkeypatch.setattr(held, "verify_anchor", lambda *args: None)
    bins = {arm: [1., *([0.]*13)] for arm in ARMS}
    monkeypatch.setattr(held, "acquisition", lambda _: {"anchor": anchor, "settings": {}, "bins": bins,
                                                       "floor_ppt": held.SCORE})
    with pytest.raises((ValueError, OSError)):
        held.main()
    result = json.loads(receipt.read_bytes())
    assert result["status"] == "incomplete" and result["passed"] is None
    assert result["acquisition"]["anchor"] == anchor
    if data is not None:
        assert b64decode(result["independent_base64"], validate=True) == data


def test_nonfinite_native_output_keeps_an_incomplete_receipt(tmp_path, monkeypatch):
    anchor, receipt = "0"*40, tmp_path/"receipt"
    monkeypatch.setattr(sys, "argv", ["gate", anchor, "--openwave-checkout", str(tmp_path),
                                      "--compare", str(tmp_path/"peer"), "--record", str(receipt)])
    monkeypatch.setattr(held, "verify_anchor", lambda *args: None)
    monkeypatch.setattr(held, "acquisition", lambda _: {"bins": [float("nan")]})
    with pytest.raises(ValueError):
        held.main()
    result = json.loads(receipt.read_bytes())
    assert (result["status"], result["stage"], result["passed"]) == ("incomplete", "acquisition", None)
    assert result["invalid_acquisition_text"] == '{"bins": [NaN]}'
