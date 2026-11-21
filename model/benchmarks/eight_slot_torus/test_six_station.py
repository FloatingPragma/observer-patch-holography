from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from . import six_station as radar


def test_exact_reflected_edges_and_rank_one_query():
    radar.check_source()
    start = radar.STATIONS[0]
    x, denom = radar.Dyadic.launch(start)
    hits, residual = radar.Dyadic.monitor((x, denom), start, radar.HORIZON-1, 0)
    # The private-axis slots stay at this address and the other six leave.
    assert hits[-1] == Fraction(1, 16)
    assert sum(hits) + radar.Dyadic.norm(residual) == 1
    assert sum(residual[0][start]) == 0
    assert len(radar.STATIONS) == len(set(radar.STATIONS)) == 6
    assert radar.MARK not in radar.STATIONS


@pytest.mark.parametrize("marked", [0, 1])
def test_retained_histories_and_source_map(marked):
    source, target = 2, 4
    ranges, tags = radar.exact_menus(source, target, marked)
    assert sum(ranges.flat) == 1
    assert np.all(tags.sum(axis=(-1, -2)) == 1)
    np.testing.assert_allclose(ranges.astype(float), radar.literal_history(source, target, marked),
                               rtol=0, atol=radar.SOURCE_TOLERANCE)
    for c in (radar.ACTIVE, radar.SILENT, radar.DISABLED):
        for bit in range(2):
            charge = radar.literal_history(source, target, marked, c, bit, True)
            np.testing.assert_allclose(tags[c, bit].astype(float), charge,
                                       rtol=0, atol=radar.SOURCE_TOLERANCE)
    assert np.all(tags[radar.SILENT, :, 1:, 1:] == 0)
    assert np.all(tags[radar.DISABLED, :, 1:, :] == 0)
    assert np.array_equal(tags[radar.ACTIVE, :, 0, :], tags[radar.SILENT, :, 0, :])
    assert np.array_equal(ranges[0], tags[radar.ACTIVE, 0, 0])
    assert not np.array_equal(tags[radar.DISABLED, 0, 0], tags[radar.SILENT, 0, 0])
    assert tags[radar.SILENT, 0, 0, 1:].sum() > 0
    assert ranges[0, 0] > 0 and ranges[1:, 0].sum() > 0 and ranges[1:, 1:].sum() > 0
    assert not np.any(tags[:, 1, 1:radar.SLOT+1])
    assert not np.any(tags[radar.ACTIVE, :, radar.SLOT+1:, 1:radar.SLOT+1])


def test_pulse_sham_and_controlled_latch():
    sham, sham_tags = radar.exact_menus(0, 2, 0)
    pulse, pulse_tags = radar.exact_menus(0, 2, 1)
    # The marked address and the queried port are distinct. A local sign
    # cannot change that port's simultaneous click, nor any earlier one.
    np.testing.assert_array_equal(sham.sum(axis=1)[1:radar.PULSE_TICK+1],
                                  pulse.sum(axis=1)[1:radar.PULSE_TICK+1])
    assert np.any(sham != pulse)
    # When both launches are late the marked services act on the vacuum in
    # both released branches. The no-B branch crosses the second pulse.
    np.testing.assert_array_equal(sham_tags[radar.ACTIVE, 1, 1:],
                                  pulse_tags[radar.ACTIVE, 1, 1:])


def test_exact_apportionment_and_support():
    law = np.array([Fraction(1, 3)]*3 + [Fraction(0)], object)
    weights = radar.integer_weights(law)
    assert weights[0] == weights[1]+1 == weights[2]+1
    assert weights.sum() == radar.GRID and weights[-1] == 0
    assert all(abs(Fraction(int(w), radar.GRID)-q) < Fraction(1, radar.GRID)
               for w, q in zip(weights, law))
    with pytest.raises(ValueError, match="incomplete"):
        radar.integer_weights(law[:2])


def test_source_support_and_count_conservation():
    shape = (1, 2, 2, 2, 2)
    laws, counts = {}, {}
    for name, shots in (("range", radar.RANGE_SHOTS), ("tag", radar.TAG_SHOTS)):
        weight = np.zeros(shape, np.int64)
        weight[:, 0, 1, 0, 1] = weight[:, 1, 0, 0, 1] = radar.GRID
        laws[name+"_time_weights"] = weight
        counts[name+"_time_counts"] = weight // radar.GRID * shots
    radar.validate_counts(laws, counts)
    mutated = counts["range_time_counts"]
    mutated[0, 0, 1, 0, 1] -= 1
    with pytest.raises(ValueError, match="loss"):
        radar.validate_counts(laws, counts)
    mutated[0, 0, 1, 0, 0] += 1
    with pytest.raises(ValueError, match="support"):
        radar.validate_counts(laws, counts)


def test_consistency_detects_redistribution_inside_source_support():
    laws, counts = {}, {}
    for name, shots in (("range", radar.RANGE_SHOTS), ("tag", radar.TAG_SHOTS)):
        weight = np.full((1, 1, 1, 1, 2), radar.GRID//2, np.int64)
        laws[name+"_time_weights"] = weight
        counts[name+"_time_counts"] = np.full(weight.shape, shots//2, np.int64)
    radar.validate_counts(laws, counts)
    counts["tag_time_counts"][0, 0, 0, 0] += np.array([11000, -11000])
    with pytest.raises(ValueError, match="simultaneous source-count bound"):
        radar.validate_counts(laws, counts)


def test_archive_checksum_and_duplicate_menus(monkeypatch):
    payload = BytesIO()
    with ZipFile(payload, "w") as archive:
        for kind in ("range", "tag"):
            for statistic in ("counts", "weights"):
                archive.writestr(f"{kind}_time_{statistic}.npy", b"")
        with pytest.warns(UserWarning, match="Duplicate"):
            archive.writestr("range_time_counts.npy", b"")
    expected = "0"*64
    monkeypatch.setattr(Path, "read_text", lambda self, encoding: expected)
    monkeypatch.setattr(Path, "read_bytes", lambda self: payload.getvalue())
    with pytest.raises(ValueError, match="checksum"):
        radar.verify_archive()
    expected = sha256(payload.getvalue()).hexdigest()
    with pytest.raises(ValueError, match="duplicate"):
        radar.verify_archive()
