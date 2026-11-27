"""Independent checks of the complete timing survey and distant-patch bound."""

import json
from pathlib import Path

import pytest

pytest.importorskip("strict_history", reason="needs the oph-model-sim checkout on PYTHONPATH")

from observer_zero_radar import local_graph, outward
from distance_redshift import inward
from strict_history.wide_replay import WideReplay

HERE = Path(__file__).resolve().parent


def test_individual_clock_paths_and_existing_aperture():
    data = json.loads((HERE / "level202_echo_reach.json").read_text())
    old = json.loads((HERE / "level202_distance_redshift.json").read_text())
    labels = [12 * int(owner) + port for owner, port in data["source_labels"]]
    # Verify the recorded replies by two separate earliest-arrival searches.
    ports, _, edges = local_graph(WideReplay(202).links, 0, 120)
    index = {p: i for i, p in enumerate(ports)}
    for epoch in (0, 4, 16, 96):
        row = next(
            r for r in reversed(data["surveys"]) if r["first_send_tick"] <= epoch
        )
        logged = dict(zip(labels, row["echo_wait_by_label"], strict=True))
        out = outward(edges, epoch, 60)
        for p, wait in logged.items():
            i = index[p]
            back = None if out[i] < 0 else inward(edges, i, epoch + out[i], 60 - out[i])
            assert wait == (None if back is None else back - epoch)
        for r in old["panels"][str(epoch)]:
            u, p = r["source"]
            wait = logged.get(12 * int(u) + p)
            actual = None if r["echo_at"] is None else r["echo_at"] - epoch
            assert actual == (wait if wait is not None and wait <= 48 else None)
    assert data["stationary_from_send"] == 15
    assert data["surveys"][-1]["responding_remote_ports"] == 168
    bound = data["antipodal_fresh_challenge"]
    assert bound["target_vertex_heights"] == [2, 2, 2]
    assert int(bound["minimum_roundtrip_ticks"]) == 12 * 2**202 - 2
    assert bound["target_echoes_ever"] == 0
