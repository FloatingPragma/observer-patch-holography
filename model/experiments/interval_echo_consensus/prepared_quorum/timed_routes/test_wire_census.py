"""Independent wire-byte and serialized-appointment audits on the public deposits."""

import gzip
import hashlib
import json
from pathlib import Path


def deposit(name):
    with gzip.open(Path(__file__).with_name(name)) as stream:
        return json.load(stream)


def public_bytes(obj):
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def with_local_epoch(obj, epoch):
    if not isinstance(obj, list):
        return obj
    body = [with_local_epoch(x, epoch) for x in obj]
    if obj and isinstance(obj[0], str):
        if obj[0] in ("PC", "LC", "DC"):
            index = 1
        elif obj[0] in ("change", "propose", "prepare", "accept", "commit"):
            index = 2
        else:
            raise AssertionError("Unknown transmitted structure")
        assert type(body[index]) is int and body[index] == epoch
        body[index] = None
    return body


def test_current_frame_bytes_and_source_failures():
    bank = deposit("local_view.json.gz")
    nominal = {1: 14478, 10: 14574, 100: 14670, 1000: 14766}
    for board, cases in bank.items():
        assert set(map(int, cases)) == set(nominal)
        for key, case in cases.items():
            view = int(key)
            if board == "0" and view != 10:
                assert "insufficient_stock" in case and "events" not in case
                continue
            total, absolute = 0, 0
            for event in case["events"]:
                packet = case["packets"][event["packet"]]
                word = public_bytes(packet)
                assert hashlib.sha256(word).hexdigest() == event["packet"]
                frame = public_bytes(with_local_epoch(packet, view))
                assert len(frame) == event["compact_bytes"] <= 8192
                assert len(word) == event["absolute_bytes"] <= 8192
                total += len(frame)
                absolute += len(word)
            assert (len(case["events"]), total, absolute) == (24, 14766, nominal[view])


def test_exact_clock_service_from_each_signed_byte_and_route():
    levels = deposit("records.json.gz")["levels"]
    outcomes = []
    for board in levels.values():
        calibration = board["clock"]
        assert calibration["read_appointments"] == 14
        for case in board["cases"].values():
            cursor = 0
            service = 0
            for event in case["traffic"]:
                assert event["start"] == cursor
                mode = "parked" if event["parked"] else "native"
                route = case["columns"][mode][f"{event['from']}/{event['to']}"]
                bits = 8 * len(public_bytes(case["packets"][event["packet"]]))
                cursor += bits * route["slots"]
                assert event["stop"] == cursor
                assert (
                    event["clock"]["terminal"] == calibration["tick"][(cursor - 1) % 2]
                )
                if not event["parked"]:
                    service = cursor
            assert service == case["communication_slots_to_finalization"]
            outcomes.append(service)
    assert (len(outcomes), min(outcomes), max(outcomes)) == (12, 4864608, 7788160)
