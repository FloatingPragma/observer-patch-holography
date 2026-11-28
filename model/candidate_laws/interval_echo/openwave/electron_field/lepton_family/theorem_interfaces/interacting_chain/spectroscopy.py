"""Momentum-resolved energy lines and a spatially separated local CAR receiver."""

from functools import cache

import numpy as np
from scipy.linalg import eigh

from . import ring


@cache
def spectrum(length, up, down, k):
    basis = ring.fourier_basis(length, up, down, k)
    h = ring.hamiltonian(length, up, down)
    block = (basis.conj().T @ h @ basis).toarray() / 16
    np.testing.assert_allclose(block, block.conj().T, atol=2e-14)
    energies, frame = eigh(block, driver="evd")
    loss = np.linalg.norm(block @ frame - frame * energies, ord="fro")
    orth = np.linalg.norm(frame.conj().T @ frame - np.eye(len(energies)), ord="fro")
    assert max(loss, orth) < 3e-11
    return energies, frame, basis


def floor(length):
    half = length // 2
    lowest = []
    for k in range(length):
        energies = spectrum(length, half, half, k)[0]
        lowest.extend((e, k) for e in energies[:2])
    lowest.sort()
    e0, k0 = lowest[0]
    energies, frame, basis = spectrum(length, half, half, k0)
    ground = basis @ frame[:, 0]
    assert np.linalg.norm(ring.spin_plus(length, half, half) @ ground) < 1e-10
    neel = sum(1 << (2 * x + (x % 2)) for x in range(length))
    index = ring.sector(length, half, half)[1]
    return e0, k0, ground, float(abs(ground[index[neel]]) ** 2), lowest[1][0] - e0


def visible_lines(length, sign, e0, k0, ground):
    half = length // 2
    # Both columns have Sz=+1/2. A positive-spin hole removes a down spin.
    up, down = (half + 1, half) if sign == 1 else (half, half - 1)
    spin, add = int(sign == -1), sign == 1
    created = ring.field(length, half, half, 0, spin, add) @ ground
    raising = ring.spin_plus(length, up, down)
    assert np.linalg.norm(raising @ created) < 1e-10
    lines = []
    for relative_k in range(length):
        k = (relative_k + k0) % length
        energies, frame, basis = spectrum(length, up, down, k)
        local = frame.conj().T @ (basis.conj().T @ created)
        start = 0
        for end in range(1, len(energies) + 1):
            if end < len(energies) and energies[end] - energies[start] < 1e-9:
                continue
            overlap = local[start:end]
            weight = float(np.vdot(overlap, overlap).real)
            if weight > 1e-12:
                eigenspace = frame[:, start:end]
                received = basis @ (eigenspace @ overlap)
                assert np.linalg.norm(raising @ received) < 1e-10
                separation = min(
                    (energies[start] - energies[start - 1]) if start else np.inf,
                    (energies[end] - energies[start])
                    if end < len(energies)
                    else np.inf,
                )
                lines.append(
                    {
                        "q_index": relative_k,
                        "sector_floor": float(energies[0]),
                        "energy": float(energies[start]),
                        "transfer": float(energies[start] - e0),
                        "reservoir_corrected_transfer": float(
                            energies[start] - e0 - sign / 2
                        ),
                        "dark_levels_below_in_spin_z_block": start,
                        "multiplicity_in_spin_z_block": end - start,
                        "nearest_other_energy_same_momentum": float(separation),
                        "local_field_weight": weight,
                    }
                )
                break
            start = end
        else:
            raise AssertionError("No accessible charged line at this momentum")
    return lines


def branch_bank(length, sign, e0, k0, ground, preparation):
    lines = visible_lines(length, sign, e0, k0, ground)
    weights = np.array([row["local_field_weight"] for row in lines])
    levels = np.array([row["transfer"] for row in lines])
    receipts = {}
    # The coherent pass/refusal energy test does not record the separate
    # momenta. Keep phase relations, including in the zero-wait control.
    for wait in (0, 16):
        amplitudes = np.fft.fft(weights * np.exp(-1j * wait * levels))
        # A local charge-conserving pump and drain each have factor 1/sqrt(2).
        receipts[str(wait)] = [float(preparation * abs(a) ** 2 / 4) for a in amplitudes]
    return {
        "sign": sign,
        "lines": lines,
        "ground_then_pump_then_branch": float(preparation * sum(weights) / 2),
        "per_original_spin_retrieval_by_separate_receiver": receipts,
    }


def lowest_window(length, sign, e0, k0, ground, preparation):
    half = length // 2
    up, down = (half + 1, half) if sign == 1 else (half, half - 1)
    added = ring.field(length, half, half, 0, int(sign < 0), sign > 0) @ ground
    catalogue = []
    for k in range(length):
        energies, frame, basis = spectrum(length, up, down, k)
        amplitudes = frame.conj().T @ (basis.conj().T @ added)
        catalogue.extend(
            (float(e), k, i, abs(amplitudes[i]) ** 2) for i, e in enumerate(energies)
        )
    catalogue.sort()
    groups = []
    for item in catalogue:
        if not groups or item[0] - groups[-1][0][0] > 1e-9:
            groups.append([])
        groups[-1].append(item)
    visible = [i for i, group in enumerate(groups) if sum(v[3] for v in group) > 1e-12][
        :2
    ]
    first, second = (groups[j] for j in visible)
    energies = np.array([first[0][0], second[0][0]])
    gap = energies[1] - energies[0]
    upper_neighbour = groups[visible[1] + 1][0][0]
    columns = []
    for group in (first, second):
        projected = np.zeros_like(added, dtype=complex)
        for k in sorted({row[1] for row in group}):
            _, frame, basis = spectrum(length, up, down, k)
            v = basis @ frame[:, [row[2] for row in group if row[1] == k]]
            projected += v @ (v.conj().T @ added)
        columns.append(projected)
    weights = [float(np.vdot(vector, vector).real) for vector in columns]
    accepted = sum(weights)
    words = np.array(ring.sector(length, up, down)[0])
    densities = np.array(
        [[(words >> (2 * x + s)) & 1 for s in range(2)] for x in range(length)]
    )
    sink = np.array(
        [
            ring.field(length, half, half, x, int(sign < 0), sign > 0) @ ground
            for x in range(length)
        ]
    )
    read = {}
    for label, phase in (
        ("parked", 1),
        ("wait_16", np.exp(-16j * gap)),
        ("half_beat", -1),
    ):
        wave = columns[0] + phase * columns[1]
        spatial = densities @ (abs(wave) ** 2 / accepted)
        single = ((densities[:, 0] - densities[:, 1]) ** 2) @ (
            abs(wave) ** 2 / accepted
        )
        returns = preparation * abs(sink.conj() @ wave) ** 2 / 4
        charge = spatial[:, 0] + spatial[:, 1] - 1
        spin = spatial[:, 0] - spatial[:, 1]
        np.testing.assert_allclose([sum(charge), sum(spin)], [sign, 1], atol=1e-10)
        read[label] = {
            "joint_spin_return": returns.tolist(),
            "local_charge_excess_given_pass": charge.tolist(),
            "local_twice_spin_z_given_pass": spin.tolist(),
            "local_single_occupation_given_pass": single.tolist(),
            "local_spin_retention_given_single": (spin / single).tolist(),
        }
    dephased = preparation * sum(abs(sink.conj() @ col) ** 2 for col in columns) / 4
    np.testing.assert_allclose(
        dephased,
        (
            np.array(read["parked"]["joint_spin_return"])
            + np.array(read["half_beat"]["joint_spin_return"])
        )
        / 2,
        atol=1e-12,
    )
    return {
        "corrected_levels": (energies - e0 - sign / 2).tolist(),
        "next_full_sector_energy": float(upper_neighbour - e0 - sign / 2),
        "nearest_excluded_above_band": float(upper_neighbour - energies[1]),
        "upper_gate_midpoint_in_corrected_units": float(
            (energies[1] + upper_neighbour) / 2 - e0 - sign / 2
        ),
        "half_beat_wait": float(np.pi / gap),
        "local_unscaled_weights": weights,
        "window_pass_per_original": float(preparation * accepted / 2),
        "levels_inside_gate": [
            {
                "corrected": float(g[0][0] - e0 - sign / 2),
                "spin_z_block_rank": len(g),
                "local_weight": float(sum(c[3] for c in g)),
            }
            for g in groups[: visible[1] + 1]
        ],
        "three_waits": read,
        "energy_resolved_then_forgotten_all_waits": dephased.tolist(),
    }


def experiment(length):
    e0, k0, g, preparation, gap = floor(length)
    banks = [branch_bank(length, sg, e0, k0, g, preparation) for sg in (1, -1)]
    # The staggered particle/hole map exchanges addition and removal and
    # shifts their relative momenta by pi. Their clock readings differ by U.
    for j, p in enumerate(banks[0]["lines"]):
        m = banks[1]["lines"][(j + length // 2) % length]
        np.testing.assert_allclose(
            [p["local_field_weight"], p["reservoir_corrected_transfer"]],
            [m["local_field_weight"], m["reservoir_corrected_transfer"]],
            atol=1e-10,
        )
    for wait in ("0", "16"):
        np.testing.assert_allclose(
            banks[0]["per_original_spin_retrieval_by_separate_receiver"][wait],
            banks[1]["per_original_spin_retrieval_by_separate_receiver"][wait],
            atol=1e-10,
        )
    for bank in banks:
        bank["two_energy_window"] = lowest_window(
            length, bank["sign"], e0, k0, g, preparation
        )
    for name in ("parked", "wait_16", "half_beat"):
        a, b = (bank["two_energy_window"]["three_waits"][name] for bank in banks)
        for key in a:
            factor = -1 if key == "local_charge_excess_given_pass" else 1
            np.testing.assert_allclose(a[key], factor * np.array(b[key]), atol=1e-10)
    return {
        "length": length,
        "half_filled_ground": float(e0),
        "ground_translation_index": k0,
        "half_filled_ground_gap": float(gap),
        "neel_ground_preparation_probability": preparation,
        "sign_banks": banks,
    }
