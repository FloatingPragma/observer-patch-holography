import json

import pytest

from .apparatus import Checkpoint, primary_read
from .measure import (
    HELD,
    HERE,
    IDENTITY,
    PACKET,
    SLOTS,
    access,
    bits,
    bits_word,
    certify,
)


def test_complete_deposited_tables_and_erasure_codes():
    archive = json.loads((HERE / "records.json").read_text())
    packet = json.loads(PACKET.read_text())
    assert archive["authorized_sets"] == packet["publicness_policy"]
    for bank in archive["boards"]:
        assert len(set(bank["owners"])) == 12
        assert (
            certify(bank["complete_words"], packet)
            == bank["checked_compositions"]
            == 1600
        )
    tables = archive["boards"][0]["complete_words"]
    certificates = access(tables)
    assert certificates == archive["interfaces"]
    assert [row["maximum_code_size"] for row in certificates] == [2] * 10 + [24, 24]
    assert (
        sum(
            sum(row["fixed_identity_apertures_by_distinct_outputs"].values())
            for row in certificates
        )
        == 4095
    )
    tables[IDENTITY][0] ^= 1
    with pytest.raises(ValueError, match="material joint data"):
        certify(tables, packet)


def test_native_held_routes_complete_and_missing_are_independent():
    bank = json.loads((HERE / "records.json").read_text())["boards"][0]
    bench = Checkpoint(0)
    for i in (0, 7, 23):
        bench.prepare(SLOTS[i])
        baseline = bench.idle_snapshot()
        bench.apply(HELD)
        bench.audit_closure(baseline)
        assert bits_word(bench.reads()) == bank["complete_words"][HELD][i]
        full = bench.relay()
        omitted = bench.relay(omit=True)
        assert full == bank["controls"]["service"][i]
        assert omitted == bank["controls"]["first_handoff_omitted"][i]
        assert len({primary_read(row) for row in full}) == 1
        assert sum(primary_read(row) is None for row in omitted) == 7
        assert len({primary_read(row) for row in omitted}) == 2
        bench.prepare(SLOTS[i])
        bench.apply(HELD, parked=True)
        assert bits_word(bench.reads()) == bank["controls"]["parked_routes"][i]
    assert primary_read(bits(bank["complete_words"][IDENTITY][0])[:-1]) is None
