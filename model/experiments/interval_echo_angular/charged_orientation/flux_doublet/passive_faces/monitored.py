"""A charged magnetic doublet from prepared direction heads and a supplied monitor."""

from functools import cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ...live_trimer.measure import LAW, blocks
from ..measure import CELL, faces, face_edges, induced_hamiltonian, links, spatial_order

HERE = Path(__file__).resolve().parent
SHOTS = 2**20


@cache
def pair(sign, angle):
    # Same head spinor for vacancy and excess. Addresses of the two block columns
    # exchange when the sign changes; endpoint reversal conjugates the connection.
    chi = np.array([1, np.exp(-1j * angle)]) / np.sqrt(2)
    u, v = map(LAW.SRC.index, ((1 + sign, 1), (1, 1 + sign)))
    records = [u, v, 9 + u, 9 + v]
    absolute = blocks()[(1, 1)][np.ix_(records, records)] / 8
    filled = LAW.SRC.index((1, 1))
    np.testing.assert_array_equal(
        blocks()[(1, 1)][filled::9, filled::9] / 8, 8 * np.eye(2)
    )
    h = absolute - 8 * np.eye(4)
    embed = np.kron(chi[:, None], np.eye(2))
    z = embed.conj().T @ h @ embed
    leakage = embed.conj().T @ (h @ h) @ embed - z @ z
    # The rational local polynomial certificates follow from the directed -4
    # dimer and its matched host. Their compressed entries contain unit phases.
    np.testing.assert_array_equal(h @ h, -8 * h)
    np.testing.assert_allclose(z @ z, -4 * z, atol=1e-13)
    np.testing.assert_allclose(leakage, -4 * z, atol=1e-13)
    return h, embed, z, leakage


def monitored_step(sign, chart, slices, time):
    dt = time / (2 * slices)
    step = np.eye(12, dtype=complex)
    for a, b in [*CELL.seams, *reversed(CELL.seams)]:
        h, embed, _, _ = pair(sign, 2 * np.pi * int(chart[a, b]) / 20)
        # A filled pair outside the two addresses has zero relative generator.
        checked = embed.conj().T @ expm(-1j * dt * h) @ embed
        step[[a, b]] = checked @ step[[a, b]]
    return np.linalg.matrix_power(step, slices)


def material_face(face, chart, sign):
    # The same stored even head is used for both species. Read the full 99-state
    # absolute crossing, and leave all heads unread after the ordinary one lap.
    native = blocks()[(1, 1)] / 8
    pulse = expm(-1j * np.pi / 8 * native)
    first = 1.0 + 0j
    survival = 1.0
    two = 1.0 + 0j
    checked = 1.0 + 0j
    for a, b in face_edges(face):
        u, v = sorted((a, b))
        angle = 2 * np.pi * chart[u, v] / 20
        chi = np.array([1, np.exp(-1j * angle)]) / np.sqrt(2)
        loads = (1 + sign, 1) if a == u else (1, 1 + sign)
        i, j = map(LAW.SRC.index, (loads, loads[::-1]))
        transfer = pulse[j::9, i::9]
        hop = transfer @ chi
        first *= np.vdot(chi, hop)
        survival *= np.vdot(hop, hop).real
        two *= np.vdot(chi, transfer @ hop)
        checked *= np.vdot(chi, hop) ** 2
    np.testing.assert_allclose(two, 0.0, atol=1e-13)
    return first, survival, checked


def material_face_records(chart, rng):
    data = []
    for sign in (-1, 1):
        for face in faces():
            first, survival, checked = material_face(face, chart, sign)
            phase = np.exp(-1j * sign * np.pi / 10)
            np.testing.assert_allclose(
                [first, survival, checked],
                [phase / 8, 1 / 8, phase**2 / 64],
                atol=3e-14,
            )
            banks = []
            for amplitude, p in (
                (first, survival),
                (0.0, 0.0),
                (checked, abs(checked) ** 2),
            ):
                shots = []
                for angle in (0, np.pi / 2):
                    x = (amplitude * np.exp(-1j * angle)).real
                    prob = np.maximum(
                        [(1 + p + 2 * x) / 4, (1 + p - 2 * x) / 4, (1 - p) / 2], 0.0
                    )
                    shots.append(rng.multinomial(SHOTS, prob / sum(prob)).tolist())
                banks.append(shots)
            xy = [(r[0] - r[1]) / SHOTS for r in banks[0]]
            inferred = round(20 * np.arctan2(xy[1], xy[0]) / (2 * np.pi))
            assert inferred == -sign
            data.append([sign, list(face), inferred, *banks])
    return data


def experiment():
    chart = links(1)
    rng = np.random.default_rng(314591)
    rows = []
    for sign in (-1, 1):
        z = np.zeros((12, 12), complex)
        leak = np.zeros_like(z)
        for a, b in CELL.seams:
            blockslice = np.ix_((a, b), (a, b))
            _, _, term, loss = pair(sign, 2 * np.pi * int(chart[a, b]) / 20)
            z[blockslice] += term
            leak[blockslice] += loss
        comparison = 4 * induced_hamiltonian(-sign * chart) - 180 * np.eye(12)
        np.testing.assert_allclose(z, comparison, atol=1e-12)
        np.testing.assert_allclose(leak, -4 * z, atol=1e-12)
        energies, cols = np.linalg.eigh(z)
        levels = np.unique(np.round(energies, 10), return_counts=True)
        assert np.array_equal(levels[1], [2, 4, 6])
        assert spatial_order(-sign * chart) == [10] * 15
        assert all(
            (sum((-sign * chart[a, b]) for a, b in face_edges(face)) % 20)
            == (-sign) % 20
            for face in faces()
        )
        ground = cols[:, :2]
        time = 0.5
        exact = expm(-1j * time * z)
        appointments = []
        for slices in (8, 32, 128, 512):
            actual = monitored_step(sign, chart, slices, time)
            error = float(np.linalg.norm(actual - exact, 2))
            state = actual[:, 0]
            stats = abs(state) ** 2
            probs = np.r_[stats, max(0.0, 1 - sum(stats))]
            counts = rng.multinomial(SHOTS, probs / sum(probs))
            out = actual @ ground
            retained = np.linalg.eigvalsh(out.conj().T @ out)
            departure = out - ground @ (ground.conj().T @ out)
            appointments.append(
                dict(
                    slices=slices,
                    per_edge_time=time / (2 * slices),
                    head_projections=60 * slices,
                    transfer_operator_error=error,
                    sharp_port_0_yield=float(sum(stats)),
                    sharp_port_0_counts_with_refusals=counts.tolist(),
                    lower_doublet_yield_range=retained.tolist(),
                    lower_doublet_band_leakage=float(np.linalg.norm(departure, 2) ** 2),
                )
            )
        assert all(
            a["transfer_operator_error"] > b["transfer_operator_error"]
            for a, b in zip(appointments, appointments[1:])
        )
        rows.append(
            dict(
                charge=sign,
                energy_levels=levels[0].tolist(),
                multiplicities=levels[1].tolist(),
                autonomous_out_of_prepared_heads_norm_squared=(-4 * levels[0]).tolist(),
                perpendicular_negative_order_pairs=15,
                appointments=appointments,
            )
        )
    return dict(
        phase_origin="fixed common thirty-head product at links(+1); the positive charge conjugates it",
        source="sharp port zero; lowest-doublet survival also evaluated on the two eigenvectors",
        service="after every addressed contact a supplied coherent rank-one direction-head monitor must accept",
        face_count_order=[
            "charge",
            "face",
            "inferred_flux_from_first_lap",
            "one_lap_unread_XY",
            "two_laps_unread_XY",
            "two_laps_monitored_XY",
        ],
        phase_bank_bins=["bright", "dark", "route_or_monitor_refusal"],
        material_face_counts=material_face_records(chart, rng),
        exposure=time,
        starts_per_port_bank=SHOTS,
        seed=314591,
        records=rows,
    )


if __name__ == "__main__":
    result = experiment()
    data = result.pop("material_face_counts")
    body = json.dumps(result, indent=2)[:-2]
    body += (
        ',\n  "material_face_counts": [\n    '
        + ",\n    ".join(json.dumps(row) for row in data)
        + "\n  ]\n}\n"
    )
    (HERE / "monitored.json").write_text(body)
