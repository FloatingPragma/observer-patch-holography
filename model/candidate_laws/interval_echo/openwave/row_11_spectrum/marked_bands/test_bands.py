import json

import numpy as np
from scipy.linalg import expm

from . import material, measure, paired


def test_full_seam_and_orientation_source():
    law = material.LAW
    x, c = law.local_term_parts()
    h = x + 8 * c
    fill = law.SRC.index((1, 1))
    host = np.kron(h[fill::9, fill::9], np.eye(9))
    occupied = np.diag(np.tile([int(w != (1, 1)) for w in law.SRC], 11))
    assert np.max(abs(host @ h - h @ host)) < 2e-12
    for term in (33 * occupied - (h - host), 9 * occupied + (h - host)):
        assert min(np.linalg.eigvalsh(term)) > -1e-12
    # Using an odd contrast head is an independently achievable one-unit detuning.
    for parity, shift in [(1, 32), (-1, 33)]:
        v = sum(
            np.eye(99)[:, law.SRC.index(word)] * amp
            for word, amp in [((0, 2), 1 / np.sqrt(2)), ((2, 0), parity / np.sqrt(2))]
        )
        filled = law.echo() @ v
        assert np.linalg.norm(h @ filled) < 3e-13
        for word in [(0, 1), (2, 1)]:
            inserted = np.zeros((11, 9))
            inserted[:, law.SRC.index(word)] = filled.reshape(11, 9)[:, fill]
            np.testing.assert_allclose(
                h @ inserted.ravel(), shift * inserted.ravel(), atol=1e-13
            )
    assert material.formation_certificate()["absolute_energy_per_defect_ceiling"] == 198


def test_frequency_is_read_from_positions_and_spacing_is_physical():
    rows = json.loads((material.HERE / "results.json").read_text())["cases"]
    for row in rows:
        b = row["band"]
        q = row["rest_q"]

        def energy(theta):
            edge = np.exp(1j * theta)
            h = [
                [17, -0.5, -0.5 * edge.conjugate()],
                [-0.5, 49, -0.5],
                [-0.5 * edge, -0.5, 17],
            ]
            return np.linalg.eigvalsh(h)[b]

        target = energy(q)
        lo, hi = row["energy"]["energy_interval"]
        assert lo < target < hi
        gap = energy(q + measure.DELTA) - target
        lo, hi = row["recoil"]["gap_interval"]
        assert lo < gap < hi
        assert all(sum(cohort) == measure.SHOTS for cohort in row["recoil"]["by_cell"])
        assert (
            row["flight"][0]["displacement"] < -1.8
            and row["flight"][1]["displacement"] > 1.8
        )
    assert (
        rows[0]["energy"]["energy_interval"][1]
        < rows[1]["energy"]["energy_interval"][0]
    )
    assert rows[0]["recoil"]["gap_interval"][1] < rows[1]["recoil"]["gap_interval"][0]
    assert material.enclosure(-1) == material.enclosure(1)


def test_three_modes_from_separated_sources_on_real_twin_tracks():
    r = json.loads((material.HERE / "paired.json").read_text())
    assert r["remote_store"] == paired.store()
    assert len(set(p for m in r["side_contacts"] for p in m)) == 24
    for row in r["cases"]:
        h = paired.guide(row["sign"])[: paired.N, : paired.N]
        # Independently propagate by powering a short full-matrix exponential.
        propagator = np.linalg.matrix_power(expm(-1j * h), row["time"])
        local = np.linalg.eigh(h[np.ix_(paired.LEFT, paired.LEFT)])[1][:, row["band"]]
        arrived = (
            abs(local[::-1] @ propagator[np.ix_(paired.RIGHT, paired.LEFT)] @ local)
            ** 2
        )
        assert abs(arrived - row["channel"]["measured_mixed_arrival"]) < 5e-5
        assert (
            row["initial_receiver_probability"] == 0 and row["doublet_weight"] > 0.99975
        )
        assert row["remote_orbitals_and_miss"][row["band"]] / paired.SHOTS > 0.9995
        assert (
            len(row["remote_orbitals_and_miss"]) == 4
            and sum(row["remote_orbitals_and_miss"]) == paired.SHOTS
        )
        assert row["channel"]["unnormalized_partial_transpose_min"] < -0.498
        assert all(
            dark == 0
            for bright, dark, miss in row["spatial_orders"][
                "perpendicular_half_turn_orders"
            ]
        )
        assert all(
            bright == 0
            for bright, dark, miss in row["spatial_orders"][
                "buffer_swap_and_local_phase_orders"
            ]
        )
        assert abs(np.diff(row["two_lines"])[0] - row["beat_from_pilot_counts"]) < 6e-10
