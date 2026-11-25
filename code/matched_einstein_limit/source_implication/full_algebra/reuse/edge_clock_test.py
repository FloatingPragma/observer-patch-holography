import json
from pathlib import Path

import numpy as np

from . import edge_clock, edge_clock_history, full_clock


def test_charge_rows_and_small_middle_asymmetry():
    record = edge_clock.run()
    assert record == json.loads(Path(edge_clock.__file__).with_suffix(".json").read_text())
    assert record["ports"] == [0, 1, 38, 39]
    for edge in record["edges"]:
        for menu in ("forward", "reverse", "orientations_pooled"):
            assert edge["menus"][menu]["joint"]["rank"] == 20
            assert edge["menus"][menu]["common_dimension"] == 0
        for menu in ("both_tagged", "opposite_assignments"):
            assert edge["menus"][menu]["joint"]["rank"] == 31
            assert edge["menus"][menu]["common_dimension"] == 9
        assert edge["common_dimension_all_four_charge_assignments"] == 0
    sharp = record["middle_precision"]
    assert all(3e-25 < x < 6e-25 for x in sharp["raw_coordinate_l_infinity_lower"])


def test_joint_timestamps_against_independent_native_solver():
    record = edge_clock_history.run()
    assert record == json.loads(Path(edge_clock_history.__file__).with_suffix(".json").read_text())
    assert all(r["through_tick_9"]["rank"] == r["through_tick_8_with_relay"]["rank"] == 31
               for r in record["reported_routes"])
    assert all(r["rank"] == r["through_tick_11_with_relay"]["rank"] == 63
               for r in record["unknown_charge_joint_state"])
    station = full_clock.station
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    directions = station.directions()
    directions[-1] = np.stack((-np.eye(4), np.eye(4)))/4
    for menu in record["native_affine_clock"]:
        q = menu["charge"]
        i, j = menu["route"]
        effects = np.zeros((77, 2, 4, 4), complex)
        for parity in (0, 1):
            h = station.iso.ranging.generator(board, ports, q, 8, bool(parity))
            kernel = station.iso.ranging.kernels(h, ports, .5, 12)
            forward = station.iso.ranging.first_hits(kernel, i, j)
            reverse = station.iso.ranging.first_hits(kernel, j, i)
            for k, (a, t) in enumerate(record["event_labels"]):
                amplitude = reverse[t-a]@forward[a]
                effects[k, parity] = amplitude.conj().T@amplitude
            for k in range(11):
                effects[66+k, parity] = forward[k+1].conj().T@forward[k+1]
        denom = int(record["affine_denominator"])
        ar = np.asarray(menu["matrix"]+menu["relay_matrix"], dtype=float)/denom
        off = np.asarray(menu["offsets"]+menu["relay_offsets"], dtype=float)/denom
        np.testing.assert_allclose(ar, np.einsum("rbxy,dbyx->rd", effects, directions).real,
                                   rtol=0, atol=1e-13)
        np.testing.assert_allclose(off, effects.trace(axis1=2, axis2=3).sum(axis=1).real/8,
                                   rtol=0, atol=1e-13)
