import numpy as np
import pytest

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.native_tidal_response import measure
from . import gateway_delivery as delivery


def test_station_ancestry_and_one_load_controls():
    board = federation(0)
    blanks, _ = measure.local_calibration()
    _, parked, _, tidal_mask, _ = measure.build_station(board, blanks)
    assert [delivery.station(board, k)[0] for k in range(board.cells)] == parked
    delivery.exact_contact_checks()
    report = delivery.compare_nomination(board)
    assert (report["source"], report["relays"], report["receiver"]) == (156, (157, 166), 163)
    assert report["tested_masks"] == 32768 and report["accepted_masks"] == 1
    mask = set(report["neutral_writes"])
    assert len(mask) == 11 and len(mask.intersection(board.external)) == 2
    assert not tidal_mask.intersection(board.external)
    before = report["challenges_before_and_after"]["blank"]
    assert delivery.switches_from_records(156, before[::-1]) == ((157, 166), 163, mask)
    with pytest.raises(ArithmeticError, match="duplicate"):
        delivery.switches_from_records(156, before + before[:1])
    inconsistent = [(p, q, r, 0 if (p, q) == (156, 157) else hit) for p, q, r, hit in before]
    with pytest.raises(ArithmeticError, match="inconsistent reciprocal"):
        delivery.switches_from_records(156, inconsistent)
    assert report["addressed_blank_delivery"] == [(163,), (163,)]
    assert report["stop_on_receipt_with_direct_link"] == ["17/75", "32/75", "26/75", "0"]
    assert report["two_attempt_bins_source_relays_receiver_other"] == {
        "passive": ["32/75", "38/75", "1/15", "0"],
        "blank": ["14/75", "4/25", "1/15", "44/75"],
        "first_relay_unavailable": ["101/150", "22/75", "1/30", "0"],
        "second_relay_unavailable": ["101/150", "22/75", "1/30", "0"],
    }
    p = report["port_probabilities_at_half_exposure"]
    assert len(p) == board.n and np.isclose(sum(p), 1, rtol=0, atol=1e-11)
    np.testing.assert_allclose(p[163], 0.000040929832, atol=1e-11, rtol=0)


@pytest.mark.parametrize(
    ("level", "unique", "refused", "eleven"), [(0, 12, 8, 12), (1, 80, 0, 64), (2, 320, 0, 224)]
)
def test_complete_home_census(level, unique, refused, eleven):
    results = delivery.survey(federation(level))
    assert len(results["homes"]) == results["cells"] == unique + refused
    assert results["counts"].get("wrong relay count", 0) == refused
    assert results["counts"]["unique"] == unique
    sizes = [layout["writes"] for row in results["homes"] for layout in row["layouts"]]
    assert sizes.count(11) == eleven and sizes.count(10) == unique - eleven


def test_missing_write_crosstalk_and_invalid_echo():
    board = federation(0)
    source, relays, (sink,) = delivery.station(board, 13)
    _, cut = delivery.wiring(board, source, relays, sink)
    for edge in cut:
        assert delivery.delivery_sink(board, source, relays, cut - {edge}) is None
    with pytest.raises(ArithmeticError, match="nonvacant"):
        delivery.replay_installation(board, (156, 157), [("echo", (156, 163))])
    with pytest.raises(ArithmeticError, match="previously written"):
        delivery.replay_installation(board, (156, 157), [("echo", (152, 166))] * 2)
    with pytest.raises(ArithmeticError, match="previously written"):
        delivery.replay_installation(board, (156, 157), [("swap", (157, 156))])


def test_pair_control_with_direct_receiver(monkeypatch):
    from . import open_pair_check

    board = federation(0)
    report = delivery.compare_nomination(board)
    monkeypatch.setattr(open_pair_check, "TIMES", (1,))
    monkeypatch.setattr(open_pair_check, "STRENGTHS", (0,))
    (row,) = delivery.direct_link_pair_control(board, report)
    assert len(row["bins"]) == 31
    np.testing.assert_allclose(sum(row["bins"]), 1, rtol=0, atol=2e-10)
    np.testing.assert_allclose(1 - row["bins"][-1], 0.982571099, rtol=0, atol=1e-9)
