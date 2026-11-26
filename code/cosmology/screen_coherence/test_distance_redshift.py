"""Independent time-of-flight and first-hit checks for observer zero."""

from collections import defaultdict
from fractions import Fraction
import json
from math import inf, isfinite, nextafter
import os
from pathlib import Path

import numpy as np
import pytest

from distance_redshift import FlashMeter, gateway_image, inward, timed_sources
from observer_zero_radar import ArchiveRadar

HERE = Path(__file__).resolve().parent


def forward_hits(links, station, launch, deadlines):
    """Enumerate transit probability forward, permanently absorbing arrivals at zero."""
    live = {station: 1.0}
    hits = 0.0
    readings = []
    for tick in range(launch, max(deadlines) + 1):
        hits += live.pop(0, 0.0)
        if tick in deadlines:
            readings.append(hits)
        if tick == max(deadlines):
            break
        after = defaultdict(float)
        for source, mass in live.items():
            count = 0
            for target, closure in links(source):
                if not closure or tick < closure:
                    after[target] += mass / 6
                    count += 1
            after[source] += mass * (6 - count) / 6
        live = after
    assert abs(hits + sum(live.values()) - 1) < 1e-12
    return readings


def test_closure_direction_and_bearings():
    edges = [
        {1: 2, 6: 0},
        {0: 2, 2: 0},
        {1: 0, 3: 0},
        {2: 0, 4: 0},
        {3: 0, 5: 0},
        {4: 0, 6: 0},
        {5: 0, 0: 0},
    ]
    r = timed_sources([0, 12, 13, 14, 15, 16, 17], edges, [14], 0)[0]
    assert (r["reflection_at"], r["echo_at"]) == (3, 7)
    assert inward(edges, 3, 0, 4) == 4
    assert inward(edges, 3, 0, 3) is None
    assert inward(edges, 1, -1, 4) is None
    gates = [{2: 2}, {2: 0}, {0: 2, 1: 0}, {}]
    assert gateway_image([1, 3, 12, 13], gates, 2, 1) == [1, 3]
    assert gateway_image([1, 3, 12, 13], gates, 2, 2) == [3]
    assert gateway_image([1, 3, 12, 13], gates, 3, 0) == []
    meter = FlashMeter(edges)
    for start in (0, 1, 2, 4):
        deadlines = [start, start + 4, start + 8, start + 32]
        probabilities = meter.probabilities([(d, 3, start, d) for d in deadlines])
        np.testing.assert_allclose(
            [probabilities[d] for d in deadlines],
            forward_hits(lambda u: edges[u].items(), 3, start, deadlines),
            atol=1e-14,
            rtol=1e-13,
        )


def test_survey_receipts_and_cohort():
    data = json.loads((HERE / "distance_redshift.json").read_bytes())
    bank = json.loads((HERE / "observer_zero_radar.json").read_bytes())
    reference = data["local_reference"]
    assert reference["radar_distance"] == 1 and reference["z_pulse"] == 0
    samples = {}
    for name, rows in data["panels"].items():
        sample = {}
        for i, row in enumerate(rows):
            source = tuple(row["source"])
            assert source == tuple(bank["reflectors"][i][:2])
            wait = (
                2 * bank["reflectors"][i][2]
                if name == "transparent"
                else bank["echo_waits"][bank["emission_ticks"].index(int(name))][i]
            )
            assert row["radar_distance"] == (None if wait is None else wait / 2)
            sent, received = row["beacon_sent"], row["beacon_received"]
            if sent is not None:
                assert sent[1] - sent[0] == data["pulse_period"]
            z = (
                None
                if received is None or None in received
                else ((received[1] - received[0]) / (sent[1] - sent[0]) - 1)
            )
            assert row["z_pulse"] == z
            if wait is not None and z is not None:
                sample[source] = z
                flash = row["flash"]
                assert (
                    0
                    <= flash["counts"][0]
                    <= flash["counts"][1]
                    <= data["launches_per_source"]
                )
                for gate, count in enumerate(flash["counts"]):
                    distance = flash["equivalent_distances"][gate]
                    if not count:
                        assert distance is None
                        continue
                    # Test the count ratio directly. Different libm pow implementations
                    # can round the archived square root to adjacent binary64 values.
                    assert isfinite(distance) and distance >= 0
                    ratio = Fraction(reference["flash"]["counts"][gate], count)
                    lower = Fraction(nextafter(distance, 0.0))
                    upper = Fraction(nextafter(distance, inf))
                    assert lower**2 <= ratio <= upper**2
        samples[name] = sample
    assert [len(samples[e]) for e in ("0", "4", "96", "transparent")] == [
        156,
        144,
        144,
        274,
    ]
    early_red = {k for k, z in samples["0"].items() if z > 0}
    assert len(early_red) == 14 and early_red <= samples["96"].keys()
    assert set(samples["96"].values()) == set(samples["transparent"].values()) == {0}


def test_forward_counts_on_archive():
    path = os.environ.get("OPH_LEVEL13_ARCHIVE")
    if not path:
        pytest.skip(
            "Set OPH_LEVEL13_ARCHIVE to verify first arrivals on the archived universe"
        )
    archive = ArchiveRadar(Path(path))
    data = json.loads((HERE / "distance_redshift.json").read_bytes())

    def transparent_links(u):
        return ((v, 0) for v, _ in archive.links(u))

    for name in ("0", "96", "transparent"):
        row = next(r for r in data["panels"][name] if r["source"] == [67108867, 2])
        owner, port = row["source"]
        links = transparent_links if name == "transparent" else archive.links
        expected = forward_hits(
            links,
            12 * owner + port,
            row["reflection_at"],
            [row["echo_at"] + g for g in data["gates"]],
        )
        np.testing.assert_allclose(
            row["flash"]["first_hit_probabilities"], expected, rtol=1e-12, atol=1e-15
        )
