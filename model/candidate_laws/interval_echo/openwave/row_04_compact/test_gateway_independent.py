from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from . import gateway_independent as peer


def test_complete_published_pilot():
    # Four published appointments, all fourteen outcomes per arm.
    # No registration commit is selected or hashed here.
    public = json.loads(Path(__file__).with_name("gateway_results.json").read_text(encoding="utf-8"))
    expected = public["terminal_by_pulse"]["1.0"]
    assert public["home"] == 13 and public["exposures"] == [0, 1, 7, 37]
    for index, t in enumerate(public["exposures"]):
        obtained = peer.probabilities(peer.BOARD, home=13, start=0, direction=1, exposure=t)
        assert obtained.keys() == expected.keys()
        for arm, vector in obtained.items():
            assert len(vector) == 14
            # Public JSON stores ten decimals. Do not scale tolerance with bin size.
            np.testing.assert_allclose(vector, expected[arm][index], atol=5.1e-11, rtol=0)


@pytest.mark.parametrize("rejected", [np.linspace(.11, .96, 12), np.zeros(12), np.ones(12)])
@pytest.mark.parametrize("start,direction", [(2, -1), (9, 1)])
def test_full_sweep_weights_against_sequential_failure_histories(rejected, start, direction):
    # Synthetic independent route calculation, including immediate and total exhaustion.
    weights, exhausted = peer.source_weights(rejected, start, direction)
    expected = np.zeros(12)
    reach = 1.
    port = start
    for _ in range(48):
        expected[port] += reach
        reach *= rejected[port]
        port = (port+direction) % 12
    np.testing.assert_allclose(weights, expected, atol=1e-14, rtol=0)
    assert abs(exhausted-reach) < 1e-14
    assert abs(weights @ (1-rejected)+exhausted-1) < 1e-14


def test_independent_reader_rejects_different_board(tmp_path):
    board = tmp_path/"unrelated.npz"
    board.write_bytes(b"unverified")
    with pytest.raises(ValueError, match="scored board"):
        peer.probabilities(board, home=13, start=0, direction=1, exposure=7)


def test_registration_precedes_selection_and_receipt(tmp_path, monkeypatch):
    # Synthetic Git history, separate from the registered baseline and task.
    def git(*args):
        return subprocess.run(["git", "-C", str(tmp_path), *args], check=True,
                              capture_output=True, text=True).stdout.strip()

    git("init", "--quiet")
    git("config", "user.name", "Test")
    git("config", "user.email", "test@example.invalid")
    git("commit", "--quiet", "--allow-empty", "-m", "baseline")
    base = git("rev-parse", "HEAD")
    content = b"synthetic approved version\n"
    monkeypatch.setattr(peer, "BASE", base)
    monkeypatch.setattr(peer, "TASK_HASH", sha256(content).hexdigest())
    task = tmp_path/peer.TASK
    task.parent.mkdir(parents=True)

    def publish(version):
        task.write_bytes(version)
        git("add", "--", peer.TASK)
        git("commit", "-qm", "install task version")
        head = git("rev-parse", "HEAD")
        git("update-ref", "refs/remotes/origin/main", head)
        return head

    earlier = publish(b"proposal\n")
    output = tmp_path/"output.json"
    monkeypatch.setattr(sys, "argv", ["peer", earlier, "--openwave-checkout", str(tmp_path),
                                      "--output", str(output)])
    monkeypatch.setattr(peer, "settings_for", lambda _: pytest.fail("drew an unregistered setting"))
    with pytest.raises(ValueError, match="No exact registration"):
        peer.main()
    assert not output.exists()

    registered = publish(content)
    publish(b"different version\n")
    duplicate = publish(content)
    peer.first_registration(tmp_path, registered)
    with pytest.raises(ValueError, match="first exact registration"):
        peer.first_registration(tmp_path, duplicate)


def test_independent_failure_is_retained_without_overwrite(tmp_path, monkeypatch):
    # Exercise custody on an invalid board, before any source or flight runs.
    output = tmp_path/"independent.json"
    board = tmp_path/"wrong.npz"
    board.write_bytes(b"wrong board")
    settings = dict(home=13, start=0, direction=1, exposure=7)
    monkeypatch.setattr(peer, "first_registration", lambda *_: None)
    monkeypatch.setattr(peer, "settings_for", lambda _: settings)
    monkeypatch.setattr(sys, "argv", ["peer", "synthetic", "--openwave-checkout", str(tmp_path),
                                      "--board", str(board), "--output", str(output)])
    with pytest.raises(ValueError, match="Incorrect scored board"):
        peer.main()
    original = output.read_bytes()
    assert json.loads(original) == {"anchor": "synthetic", "settings": settings,
                                    "error": "Incorrect scored board"}
    with pytest.raises(FileExistsError):
        peer.main()
    assert output.read_bytes() == original
