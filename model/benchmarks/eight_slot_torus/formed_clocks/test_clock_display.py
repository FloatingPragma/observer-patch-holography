"""The portable replay uses complete archived local counts with physical denominators."""
import base64
import gzip
import json

import numpy as np

from . import clock_display
from .. import runs


def test_offline_replay_carries_original_launch_counts(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("display attempted to acquire another sample")
    monkeypatch.setattr(np.random, "default_rng", forbidden)
    monkeypatch.setattr(clock_display.loaded_clock, "counts_for", forbidden)
    loaded = clock_display.loaded_document()
    driven = clock_display.driven_document()
    records = loaded["observer"]
    assert records["port_order"] == ["+X", "+Y", "-X", "-Y"]
    assert len(records["counts"]) == 11
    for i, word in enumerate(records["counts"]):
        assert len(word) == 7 and sum(word) == 2**24
        for port in range(4):
            assert loaded["series"][port]["points"][i] == [records["station_counts"][i], 100 * word[port] / 2**24]
    particle = driven["observer"]
    assert len(particle["words"]) == 108
    for word in particle["words"].values():
        original = (word["unmarked"] + word["other_first_marker"] + sum(n for _, n in word["first_failures"])
                    + sum(word["ports"]) + word["missed_corner"])
        assert original == particle["trials_per_visit"] == 2**44
    received = particle["words"]["forward_4_sieved_3_booked"]
    home = particle["words"]["forward_4_sieved_3_home"]
    assert received["station_count"] == home["station_count"] == 141
    assert received["station"] != home["station"]
    assert sum(received["ports"]) > 0 and sum(home["ports"]) == 0
    paths = [runs.write(doc, tmp_path / name) for doc, name in
             ((loaded, clock_display.LOADED_RUN), (driven, clock_display.TRANSPORT_RUN))]
    for doc, path in zip((loaded, driven), paths):
        assert runs.read(path) == doc
