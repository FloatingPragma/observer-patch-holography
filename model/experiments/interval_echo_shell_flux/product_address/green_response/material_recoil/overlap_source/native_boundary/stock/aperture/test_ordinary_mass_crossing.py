"""Independent checks of the load inventory and the protected-link counterflow."""

import json
import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    occupied_mass_exchange as exchange,
)


def test_regional_partition_is_the_native_board_energy():
    words = exchange.home.source.thermal.apparatus()[0]
    board, homes, seams, _, _, full = exchange.home.source.balls.energy.returned(words)
    raw = exchange.home.source.thermal.pair.stock.kernels()[2]
    mask = np.arange(board.n) % 3 == 0
    for index in range(0, len(words), max(len(words) // 11, 1)):
        word = words[index]
        residents = np.zeros(board.n, dtype=int)
        residents[homes] = word[:4]
        local, outside = exchange.regional_energy(
            raw, board, seams, word, residents, mask
        )
        assert local + outside == 2 * full[index]
    for left in range(3):
        for right in range(3):
            lo, hi = exchange.bond_lines(raw, left, right)
            assert (lo + hi) / 2 == raw.diagonal()[3 * left + right] / 4


def test_counterflow_keeps_members_and_transfers_local_energy():
    data = json.loads((exchange.HERE / "occupied_mass_exchange.json").read_bytes())
    assert data["enrolled_eleven_marker_backgrounds"] == 0
    for row in data["boundary_edges"]:
        for bg in row.get("backgrounds", []):
            s = bg["complete_shell"]
            i, j = row["inside_outside"]
            assert i in s["within_radius"] and j not in s["within_radius"]
            assert s["whole_board_native_energy_gain"] == 0
            assert s["inside_rest_gain"] > s["inside_native_energy_gain"] > 0
            assert s["protected_timing"]["tested_links"] == sum(
                s["protected_timing"]["measured_period_counts"]
            )
    lo, hi = data["regional_calorimeter"]["inside_change_interval"]
    assert lo < 120 < hi and data["enclosed_record_entropy_increase_interval"][0] > 0
