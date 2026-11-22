from copy import deepcopy
from hashlib import sha256
import json
import subprocess
import sys

import pytest

from . import gateway_confirmation as held
from .gateway_release import ARMS, verify


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
    with pytest.raises(ValueError, match="absent"):
        held.verify_anchor(tmp_path, prior)
    first = publish(active)
    publish(b"revision\n")
    duplicate = publish(active)
    held.verify_anchor(tmp_path, first)
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
    assert result["comparison"] == {"absolute": 0., "arm": "coherent", "outcome": 0}
    assert result["independent_sha256"] == sha256(reference.read_bytes()).hexdigest()
    with pytest.raises(FileExistsError):
        held.main()
    assert receipt.read_bytes() == original
