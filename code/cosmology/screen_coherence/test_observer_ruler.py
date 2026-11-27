"""Check the published skies against forward pulse propagation and source clocks."""

from collections import deque
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("strict_history", reason="needs the oph-model-sim checkout on PYTHONPATH")

from observer_ruler import compare, sky_catalog
from read_ruler_clock import pulse_reading
from strict_history.wide_replay import WideReplay

HERE = Path(__file__).resolve().parent


def forward_array(links, source, emission, limit):
    queue = deque([(source, emission)])
    visited = {source}
    received = [None] * 12
    while queue:
        u, t = queue.popleft()
        if u < 12:
            received[u] = t - emission
        if t - emission == limit:
            continue
        for v, closure in links(u):
            if v not in visited and (not closure or t < closure):
                visited.add(v)
                queue.append((v, t + 1))
    return received


@pytest.mark.parametrize("level,pairs", [(13, 422), (202, 227)])
def test_read_only_sky_reconstruction(level, pairs):
    survey = json.loads((HERE / f"level{level}_ruler.json").read_bytes())
    gram = np.array(survey["local_direction_gram"])
    for name in ("strict", "transparent"):
        panel = survey[name]
        rebuilt = compare(sky_catalog(panel["records"]), gram)
        assert rebuilt == {k: v for k, v in panel.items() if k != "records"}
    strict = survey["strict"]
    assert strict["unique_certified_pairs"] == pairs
    for sky in strict["skies"]:
        n = sky.get("length_certified", 0)
        assert sky.get("same_bin", 0) == n
        assert sky.get("complete_arrays", 0) == n
    assert all(row[2] == 0 for row in strict["distinct_pairs_by_reading"])
    control = survey["transparent"]
    assert control["unique_certified_pairs"] == 671
    assert set(control["pulse_counts_by_extra_ticks"]) == {"0"}
    assert control["examples"]["different_bins"][
        "bin_centre_span_rad"
    ] == pytest.approx(2.03444, abs=1e-5)


def test_independent_level202_broadcasts():
    survey = json.loads((HERE / "level202_ruler.json").read_bytes())
    model = WideReplay(202, 110, 20261126)
    assert survey["strict"]["pulse_counts_by_extra_ticks"] == {
        "0": 7844,
        "1": 48,
        "2": 1,
    }
    # Reissue both crests separately at the source and read every home clock.
    for source, times, expected, shift in [
        ((3, 1), (2, 6), [7, 11], "0"),
        ((3, 0), (2, 6), [6, 11], "1/4"),
        ((2, 2), (2, 6), [8, 14], "1/2"),
        ((3, 0), (4, 8), [8, 13], "1/4"),
        ((3, 0), (3, 9), [7, 14], "1/6"),
        ((3, 0), (2, 10), [6, 15], "1/8"),
        ((3, 0), (82, 86), [87, 91], "0"),
    ]:
        result = pulse_reading(survey, source, times)
        assert (result["received"], result["redshift"]) == (expected, shift)
        i = survey["source_labels"].index([str(source[0]), source[1]])
        for e in times:
            received = forward_array(model.links, 12 * source[0] + source[1], e, 48)
            segment = next(
                row for row in reversed(survey["strict"]["records"][i]) if row[0] <= e
            )
            assert segment[1] == received
            # At e each listed partner has sent the reply on the second tick.
            actual = []
            for peer, close in model.links(12 * source[0] + source[1]):
                peer_label = [str(peer // 12), peer % 12]
                if peer_label in survey["source_labels"] and (
                    not close or e - 2 < close
                ):
                    back = dict(model.links(peer))[12 * source[0] + source[1]]
                    if not back or e - 1 < back:
                        actual.append(survey["source_labels"].index(peer_label))
            assert segment[2] == sorted(actual)
    assert all(
        row.get("identical_complete_array_patterns", 0)
        == row.get("length_certified", 0)
        for row in survey["strict"]["skies"]
    )


def test_ruler_reader_retains_misses_and_clock_resolution():
    survey = json.loads((HERE / "level202_ruler.json").read_bytes())
    missed = next(
        i
        for i, rows in enumerate(survey["strict"]["records"])
        if rows[-1][1][0] is None
    )
    answer = pulse_reading(survey, survey["source_labels"][missed], (82, 86))
    assert answer["received"] == [None, None] and answer["redshift"] is None
    with pytest.raises(ValueError):
        pulse_reading(survey, [3, 0], (4, 4))
    with pytest.raises(ValueError):
        pulse_reading(survey, [3, 0], (0, 96))
