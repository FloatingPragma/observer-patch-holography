import ast
import hashlib
import inspect
import json
from copy import deepcopy

import numpy as np
import pytest

from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import federation

from . import channel, decoder, experiment, headline
from .acquisition import counts
from .channel import HERE, SPEC, encode, readout, validate_law
from .clock_join import banks, hardware, probe


def evidence():
    return json.loads((HERE / "receipt.json").read_text())


def test_replay_and_independent_calculation():
    saved = evidence()
    actual = experiment.acquire()
    actual["clock_attempt"] = probe()
    experiment.compare(saved, actual)
    worst, rows = headline.check(SPEC, saved)
    assert worst < 1e-12 and len(rows) == 10


def test_replay_accepts_only_complete_source_versions(tmp_path, monkeypatch):
    old = hashlib.sha256(b"old").hexdigest()
    new = hashlib.sha256(b"new").hexdigest()
    monkeypatch.setattr(channel, "ROOT", tmp_path)
    monkeypatch.setattr(channel, "SPEC", {"sources": {"a": old, "b": old, "reader": old}})
    monkeypatch.setattr(channel, "REPLAY_UPDATES", {"a": new, "b": new})
    monkeypatch.setattr(channel, "REPLAY_70E58E24", {})
    monkeypatch.setattr(channel, "REPLAY_REGIONAL_CONSENSUS", {})
    monkeypatch.setattr(channel, "REPLAY_MODULAR_STANDARD", {})
    monkeypatch.setattr(channel, "REPLAY_GRAVITY_CAP", {})
    monkeypatch.setattr(channel, "REPLAY_NESTED_CAP", {})
    monkeypatch.setattr(channel, "REPLAY_JOINT_READ", {})
    monkeypatch.setattr(channel, "REPLAY_DIRECTIONAL_CHARGE", {})
    monkeypatch.setattr(channel, "REPLAY_DARK_OBSERVATIONS", {})
    for name in ("a", "b", "reader"):
        (tmp_path / name).write_bytes(b"old")
    channel.sources()
    (tmp_path / "a").write_bytes(b"new")
    with pytest.raises(ValueError, match="source: b"):
        channel.sources()
    (tmp_path / "b").write_bytes(b"new")
    channel.sources()
    latest = hashlib.sha256(b"latest").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_70E58E24", {"a": latest, "b": latest})
    (tmp_path / "a").write_bytes(b"latest")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "b").write_bytes(b"latest")
    channel.sources()
    regional = hashlib.sha256(b"regional").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_REGIONAL_CONSENSUS", {"a": regional, "b": regional})
    (tmp_path / "b").write_bytes(b"regional")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "a").write_bytes(b"regional")
    channel.sources()
    standard = hashlib.sha256(b"standard").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_MODULAR_STANDARD", {"a": standard, "b": standard})
    (tmp_path / "b").write_bytes(b"standard")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "a").write_bytes(b"standard")
    channel.sources()
    gravity = hashlib.sha256(b"gravity").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_GRAVITY_CAP", {"a": gravity, "b": gravity})
    (tmp_path / "a").write_bytes(b"gravity")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "b").write_bytes(b"gravity")
    channel.sources()
    nested = hashlib.sha256(b"nested").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_NESTED_CAP", {"a": nested, "b": nested})
    (tmp_path / "b").write_bytes(b"nested")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "a").write_bytes(b"nested")
    channel.sources()
    joint = hashlib.sha256(b"joint").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_JOINT_READ", {"b": joint})
    (tmp_path / "b").write_bytes(b"joint")
    channel.sources()
    directional = hashlib.sha256(b"directional").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_DIRECTIONAL_CHARGE", {"b": directional})
    (tmp_path / "b").write_bytes(b"directional")
    channel.sources()
    dark = hashlib.sha256(b"dark").hexdigest()
    monkeypatch.setattr(channel, "REPLAY_DARK_OBSERVATIONS", {"b": dark})
    (tmp_path / "b").write_bytes(b"dark")
    channel.sources()
    (tmp_path / "a").write_bytes(b"gravity")
    with pytest.raises(ValueError):
        channel.sources()
    (tmp_path / "a").write_bytes(b"nested")
    (tmp_path / "reader").write_bytes(b"changed")
    with pytest.raises(ValueError):
        channel.sources()


def test_receiver_inputs_and_isolation():
    tree = ast.parse(inspect.getsource(decoder))
    assert not any(isinstance(node, (ast.Import, ast.ImportFrom, ast.Attribute,
                                     ast.Global, ast.Nonlocal)) for node in ast.walk(tree))
    assert set(decoder.decode.__code__.co_names) <= {"type", "int", "bytes", "len", "max", "ValueError"}
    table = bytes([0, 1] + [1] * 11)
    trials = [dict(truth=b, port=p, seed=i, outcome=12)
              for i, (b, p) in enumerate([(0, 82), (1, 74)])]
    assert [decoder.decode(r["outcome"], table) for r in trials] == [1, 1]
    for extra in (True, -1, 13, 1., trials[0], (12, 0)):
        with pytest.raises(ValueError):
            decoder.decode(extra, table)
    for bad in (bytearray(table), table[:-1], bytes([2]) + table[1:]):
        with pytest.raises(ValueError):
            decoder.decode(12, bad)
    with pytest.raises(TypeError):
        decoder.decode(12, table, truth=1)


def test_quiet_removal_and_missing_complement():
    row = evidence()["stations"][0]["windows"][-1]
    p = np.array(row["probabilities"])
    with pytest.raises(ValueError, match="incomplete"):
        validate_law(p[:, :-1])
    for value in (-1e-15, np.nan):
        changed = p.copy()
        changed[0, 0] = value
        with pytest.raises(ValueError):
            validate_law(changed)
    bits = np.array([0, 0, 1, 1])
    outcomes = np.array([1, 12, 12, 0])
    keep = (bits != 1) | (outcomes != 12)
    with pytest.raises(ValueError, match="enrolled"):
        counts(bits[keep], outcomes[keep], 2)
    saved = evidence()
    saved["stations"][0]["windows"][-1]["native"]["counts"][1][-1] -= 1
    with pytest.raises(ValueError, match="enrolment"):
        headline.check(SPEC, saved)


def test_physical_controls_and_mutated_gate():
    saved = evidence()
    for site in saved["stations"]:
        rows = site["windows"]
        by_s = {r["exposure"]: r for r in rows}
        for row in rows:
            assert row["zero_kappa_max_delta"] < SPEC["numeric_tolerance"]
            for null in (row["native"]["permuted"], row["no_encoding"]):
                assert null["interval"][0] <= .5 <= null["interval"][1]
        slow = by_s[2]["half_exchange_probabilities"]
        np.testing.assert_allclose(slow, by_s[1]["probabilities"], atol=1e-12, rtol=0)
        assert by_s[2]["half_exchange_max_delta"] >= SPEC["minimum_law_shift"]
        assert site["planted_success"] == 1
    site = saved["stations"][0]
    row = site["windows"][-1]
    true = deepcopy(row["probabilities"])
    row["probabilities"] = row["half_exchange_probabilities"]
    with pytest.raises(AssertionError):
        headline.check(SPEC, saved)
    row["probabilities"] = [np.roll(p[:12], 1).tolist() + p[12:] for p in true]
    with pytest.raises(AssertionError):
        headline.check(SPEC, saved)


def test_stationary_clock_and_protection():
    result = probe()
    for run in result["isolated_clock_trials"]:
        count = 0
        for event in run["events"]:
            count += event["tick"]
            assert count == event["elapsed"]
        assert count == run["deadline"] and len(run["events"]) <= 8
        assert run["committed_answer"] is None
    station = SPEC["stations"][0]
    board = federation(0)
    cells = banks(board, station["ports"][0] // 12)
    timeout = hardware(board, station, cells, 9, 31)
    assert timeout["status"] == "timeout" and len(timeout["events"]) == 8
    assert result["joint_message_trials"] == 0
    np.testing.assert_allclose(result["off_block_by_bit"], [-.5, -.5], atol=1e-12, rtol=0)
    free, _ = readout(echo.generator(board), encode(board, station), SPEC["control_exposure"])
    assert np.max(abs(result["protected_one_load_law"] - free)) > SPEC["minimum_law_shift"]


def test_receipt_custody_and_failed_verdict(tmp_path, monkeypatch):
    saved = evidence()
    path = tmp_path / "receipt.json"
    path.write_bytes(b"preserve")
    monkeypatch.setattr(experiment, "acquire", lambda: deepcopy(saved))
    monkeypatch.setattr("sys.argv", ["experiment", "--record", str(path)])
    with pytest.raises(FileExistsError):
        experiment.main()
    assert path.read_bytes() == b"preserve"
    changed = deepcopy(saved)
    changed["stations"][0]["windows"][0]["native"]["quiet"] -= 1
    with pytest.raises(ValueError):
        experiment.compare(saved, changed)
    saved["stations"][0]["verdict"] = "BELOW_REGISTERED_SENSITIVITY"
    target = tmp_path / "negative.json"
    monkeypatch.setattr("sys.argv", ["experiment", "--record", str(target)])
    with pytest.raises(SystemExit) as stop:
        experiment.main()
    assert stop.value.code == 1 and json.loads(target.read_text())["stations"][0] == saved["stations"][0]


def test_first_hand_advance_contact_changes_free_port_records():
    from scipy.linalg import expm
    from scipy.sparse.linalg import expm_multiply
    from model.candidate_laws.interval_echo import clock

    board = federation(0)
    launch = encode(board, SPEC["stations"][0])
    state = expm_multiply(-2j * echo.generator(board), launch)
    a, b = clock.RING[-1], clock.RING[-2]
    assert b in board.neighbours[a]
    # First addressed pi gate of post(), applied in isolation to the signal.
    law = echo.native_law
    one = [law.SRC.index((1, 0)), law.SRC.index((0, 1))]
    gate = expm(-1j * np.pi * law.seam_term(echo.KAPPA))[np.ix_(one, one)]
    after = gate @ state[[a, b]]
    changes = np.max(abs(abs(after)**2 - abs(state[[a, b]])**2), axis=0)
    np.testing.assert_allclose(changes, [0.00711113, 0.00146173], atol=5e-9, rtol=0)
