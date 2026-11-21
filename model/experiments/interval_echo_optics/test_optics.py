import json

import numpy as np
import pytest

from model.executable_model import federation
from .experiment import HERE, N, analyze, decision, flags, material, unpack


def test_native_report_requires_far_trigger():
    board = federation(1)
    far = min(b // 12 for a, b in board.external if a // 12 == 0)
    on = material(board, far, (3, 5), 9991)
    off = material(board, far, (3, 5), 9991, trigger_enabled=False)
    assert [i for i, e in enumerate(on["events"], 1) if e["report"]] == [4]
    assert not any(e["report"] for e in off["events"])
    assert off["material_words"][3] == 0
    assert on["material_words"][3] > 0
    assert [e["handoff"] for e in on["events"]] == [e["handoff"] for e in off["events"]]
    assert on["material_words"][4:] == off["material_words"][4:]


def test_preserved_count_packets_and_paired_capacity_response():
    receipt = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "packets.npz", allow_pickle=False) as src:
        for p in receipt["panels"]:
            k = p["id"] + "_"
            values = analyze(src[k+"keys"], src[k+"counts"], src[k+"offsets"])
            assert values == p["estimates"] and src[k+"counts"].sum() == N
            assert p["decision"] == decision(values["residual"]["value"], values["residual"]["paired_block_se"])
            if p["intervals"] == [8,8]:
                rows = unpack(src[k+"keys"])
                if p["width"] == 3:
                    assert np.array_equal(rows[:,5],rows[:,6]) and p["decision"] == "disagrees"
                else:
                    assert (rows[:,5] <= rows[:,6]).all() and np.any(rows[:,5] < rows[:,6])
                    assert p["decision"] == "agrees" and p["raw_decision"] == "disagrees"
                    assert values["filtered_minus_raw"]["value"] == pytest.approx(-0.007836193147769954)


def test_owned_events_and_independent_optical_inversion():
    receipt = json.loads((HERE / "results.json").read_text())
    for witness in receipt["material_witnesses"] + receipt["five_bit_material_witnesses"]:
        guide = witness["guide"]
        home = {0,*witness["private_home"], *witness["home_word_banks"]}
        assert 0 in home and guide[-1] not in home
        assert set(guide) & home == {0}
        o, i = witness["intervals"]
        assert [t for t,(m,_,_) in enumerate(flags(witness),1) if m] == [o,i+o]
        assert [t for t,(_,r,_) in enumerate(flags(witness),1) if r] == [o+1]
    with np.load(HERE/"packets.npz", allow_pickle=False) as archive:
        for p in receipt["panels"]:
            key = p["id"]
            means = np.average(unpack(archive[key+"_keys"]), weights=archive[key+"_counts"], axis=0)
            e, c, s, signal, reunion, traveler, _ = means
            distance = (c - e) * .5
            outward_time = signal - s - distance
            inward_time = reunion - s - outward_time
            expected = sum(np.sqrt(t*t-distance**2) for t in (outward_time, inward_time))
            residual = (traveler - expected)/(reunion - s)
            assert p["estimates"]["residual"]["value"] == pytest.approx(residual,abs=1e-14)
