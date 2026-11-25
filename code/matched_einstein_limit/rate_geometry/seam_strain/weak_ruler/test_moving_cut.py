import json

import numpy as np
from scipy.linalg import expm

from . import moving_cut as moving


def test_time_derivative_of_moving_charge_in_shifted_hamiltonian():
    h = np.array([[1, 1j, 0], [-1j, 0, 2], [0, 2, -1]], complex)
    left, zone = np.diag([1, 0, 0]), np.diag([0, 1, 0])
    rho = np.outer([1, 0.5j, 1], [1, -0.5j, 1]) / 2.25
    c, v = 2.25, -0.7
    p = left + zone / 2

    def energy(op, g):
        return (op @ g + g @ op) / 2

    nflux = -1j * (h @ p - p @ h)
    eflux = -1j * (h @ h @ p - p @ h @ h) / 2
    middle = eflux + c * nflux - v * energy(zone, h) - v * c * zone

    def level(t):
        g = h + c * np.eye(3)
        u = expm(-1j * t * g)
        return np.trace(energy(p + t * v * zone, g) @ u @ rho @ u.conj().T).real

    dt = 1e-5
    rate = (level(-dt) - level(dt)) / (2 * dt)
    assert abs(np.trace(middle @ rho).real - rate) < 1e-8


def test_physical_four_arm_rank_and_blind_lottery():
    regenerated = moving.run()
    archived = json.loads((moving.face.hinge.HERE / "moving_cut.json").read_text())

    def overlap(a, b):
        assert max(a[0], b[0]) <= min(a[1], b[1])

    for record in (regenerated, archived):
        assert record["pump_age"] == 4 and record["shared_ports"] == 490
        for key in ("finite_preparations", "infinitesimal_rotations"):
            result = record[key]
            assert result["rank_four"]
            assert result["any_origin_and_sweep_max_absolute_residual_lower"] > 0
        recipe = record["finite_preparations"]
        gap = recipe["signed_energy_current_gap"]
        assert -7.14e-7 < gap[0] <= gap[1] < -7.13e-7
        assert recipe["any_origin_and_sweep_max_absolute_residual_lower"] > 6.50e-7
        assert (
            record["infinitesimal_rotations"][
                "any_origin_and_sweep_max_absolute_residual_lower"
            ]
            > 2.47e-6
        )
        assert recipe["exactly_matched"] == [
            "number_current",
            "zone_energy",
            "zone_number",
        ]
        for row in recipe["lotteries"].values():
            assert all(0 <= lo <= hi <= 1.00000000000001 for lo, hi in row.values())
            assert (
                sum(lo for lo, hi in row.values())
                <= 1
                <= sum(hi for lo, hi in row.values())
            )
    for arm, cols in regenerated["contrasts"].items():
        for k, band in cols.items():
            overlap(band, archived["contrasts"][arm][k])
    overlap(
        regenerated["finite_preparations"]["signed_energy_current_gap"],
        archived["finite_preparations"]["signed_energy_current_gap"],
    )
