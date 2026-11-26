"""Read the first variation of interval entropy and independent bond energy on one ring."""

import argparse
import json
from functools import cache
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal, expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30.experiment import select_route
from model.executable_model.carrier import federation
from model.experiments.interval_echo_gauge.port_phases.apparatus import reference

HERE = Path(__file__).resolve().parent
SIZES = (4, 8, 16, 24)
CENTRES = (48, 52)


def native_contact_checks():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    binary = [law.SRC.index((a, b)) for a in (0, 1) for b in (0, 1)]
    inclusion = np.eye(law.DIM)[:, binary]
    hop = np.zeros((4, 4))
    hop[1, 2] = hop[2, 1] = 1
    angle = 0.37
    target = expm(1j * angle * hop)
    rows = {}
    for name, ex, k in (
        ("native", 1, candidate.KAPPA),
        ("exchange_removed", 0, candidate.KAPPA),
        ("contrast_halved", 1, candidate.KAPPA / 2),
    ):
        h = ex * exchange + k * contrast
        local = inclusion.T @ h @ inclusion
        leakage = np.linalg.norm(h @ inclusion - inclusion @ local)
        if leakage > 1e-12:
            raise ArithmeticError("Binary blank input leaves the declared subspace")
        port_phase = reference(k, np.pi / candidate.KAPPA)[np.ix_(binary, binary)]
        forward = expm(-1j * angle * local)
        backward = expm(-1j * (2 * np.pi - angle) * local)
        pulse = forward @ port_phase @ backward @ port_phase
        rows[name] = {
            "blank_leakage": leakage,
            "pulse_operator_error": np.linalg.norm(pulse - target, ord=2),
            "transfer_probability": abs(pulse[1, 2]) ** 2,
        }
    return rows


def ring_modes(n):
    half = n // 2
    # The physical exchange closing the ring crosses half-1 intervening fermions.
    shift = 0 if half % 2 else 0.5
    momenta = 2 * np.pi * (np.arange(-half, half) + shift) / n
    phase = np.exp(1j * np.outer(np.arange(n), momenta)) / np.sqrt(n)
    filled = np.cos(momenta) > 0
    if filled.sum() != half or min(abs(np.cos(momenta))) < 1e-12:
        raise ArithmeticError("Half-filled reference has an ambiguous Fermi point")
    return momenta, phase, filled


def tangents(momenta, phase, filled):
    below = np.flatnonzero(filled)[-2:][::-1]
    above = np.flatnonzero(momenta > np.pi / 2)[:2]
    h, p = below[0], above[1]
    node = (51 + 55) / 2 - 0.5
    fringe_phase = np.pi / 2 - (momenta[p] - momenta[h]) * node
    rotated = np.exp(1j * fringe_phase) * phase[:, p]
    # The difference of these two projectors is the derivative of an orbital rotation.
    minus = (phase[:, h] - rotated) / np.sqrt(2)
    plus = (phase[:, h] + rotated) / np.sqrt(2)
    return {
        "plane_wave": (phase[:, below[0]], phase[:, above[0]]),
        "coherent_pair": tuple(
            np.sum(phase[:, pair], axis=1) / np.sqrt(2) for pair in (below, above)
        ),
        "signed_rotation": (minus, plus),
    }


def compiled_modes(target):
    """Compile the source basis into adjacent two-mode rotations and port phases."""
    working = target.copy()
    n = len(target)
    rotations = []
    for col in range(n - 1):
        for row in range(n - 1, col, -1):
            a, b = working[row - 1, col], working[row, col]
            norm = np.hypot(abs(a), abs(b))
            if norm == 0:
                continue
            a, b = a / norm, b / norm
            x, y = working[row - 1, col:].copy(), working[row, col:].copy()
            working[row - 1, col:] = a.conjugate() * x + b.conjugate() * y
            working[row, col:] = -b * x + a * y
            rotations.append((row, a, b))
    phases = working.diagonal().copy()
    replay = np.diag(phases / abs(phases))
    for row, a, b in reversed(rotations):
        x, y = replay[row - 1].copy(), replay[row].copy()
        replay[row - 1] = a * x - b.conjugate() * y
        replay[row] = b * x + a.conjugate() * y
    if np.linalg.norm(replay - target) > 1e-10:
        raise ArithmeticError("Native rotation schedule gives the wrong orbitals")
    return replay, len(rotations)


def weights(n, length):
    wave = np.pi / n
    midpoints = np.arange(1, length)
    return np.sin(wave * midpoints) * np.sin(wave * (length - midpoints)) / (
        wave * np.sin(wave * length)
    )


def four_site_ratios(vacancies, occupations, eigenvalues):
    modular = np.log(vacancies) - np.log(occupations)
    positive = eigenvalues > 0
    return (modular[positive] / eigenvalues[positive]).tolist()


def interval_read(n, start, end, phase, filled, sources):
    stations = np.arange(start, end)
    length = end - start
    f = weights(n, length)
    actual, modes = eigh_tridiagonal(np.zeros(length), -np.pi * f)
    amplitudes = modes.T @ phase[stations]
    # Count both outcomes on the reference. Subtraction from one loses rare counts.
    occupations = np.sum(abs(amplitudes[:, filled]) ** 2, axis=1)
    vacancies = np.sum(abs(amplitudes[:, ~filled]) ** 2, axis=1)
    if min(*occupations, *vacancies) <= 0:
        raise ArithmeticError("The selected interval reference has a null mode")
    if max(abs(occupations + vacancies - 1)) > 1e-12:
        raise ArithmeticError("The two complementary probabilities fail normalization")
    modular = np.log(vacancies) - np.log(occupations)
    c = phase[stations][:, filled] @ phase[stations][:, filled].conj().T
    bonds = np.diag(-np.pi * f, 1) + np.diag(-np.pi * f, -1)
    flat = np.diag(-np.pi * np.full(length - 1, f.mean()), 1)
    flat += flat.T.copy()
    signals = {}
    for name, (hole, particle) in sources.items():
        hole_here, particle_here = hole[stations], particle[stations]
        delta_modes = abs(modes.T @ particle_here) ** 2 - abs(modes.T @ hole_here) ** 2
        delta_bonds = -np.real(
            particle_here[1:] * particle_here[:-1].conj()
            - hole_here[1:] * hole_here[:-1].conj()
        )
        entropy = float(modular @ delta_modes)
        signed = float(2 * np.pi * f @ delta_bonds)
        signals[name] = {
            "entropy": entropy,
            "weighted_energy": signed,
            "relative_difference": abs(entropy - signed) / abs(signed),
            "flat_weight_energy": float(2 * np.pi * f.mean() * delta_bonds.sum()),
        }
    row = {
        "stations": [int(stations[0]), int(stations[-1])],
        "commutator_norm": np.linalg.norm(c @ bonds - bonds @ c, ord=2),
        "uniform_control_commutator_norm": np.linalg.norm(c @ flat - flat @ c, ord=2),
        "minimum_mode_probability": float(min(*occupations, *vacancies)),
        "readings": signals,
    }
    if length == 4:
        row["positive_modular_to_energy_ratios"] = four_site_ratios(
            vacancies, occupations, actual
        )
    return row


def inverse_stress(n, start, end, interval, which, source):
    angle = np.pi / n
    length = end - start
    moments = [
        np.sin(angle * (length + shift)) * interval(start, end + shift)["readings"][source][which]
        for shift in (-1, 0, 1)
    ]
    return (moments[2] - 2 * np.cos(angle) * moments[1] + moments[0]) / (
        2 * n * np.sin(angle) * np.sin(angle * length)
    )


def count_bank(n, phase, filled, sources):
    shots = 2**30
    seed = 24649
    rng = np.random.default_rng(seed)
    part = slice(46, 50)
    baseline = phase[part][:, filled] @ phase[part][:, filled].conj().T
    hole, particle = (orbital[part] for orbital in sources["coherent_pair"])
    excited = baseline + np.outer(particle, particle.conj()) - np.outer(hole, hole.conj())
    f = weights(n, 4)
    _, modes = eigh_tridiagonal(np.zeros(4), -np.pi * f)

    def occupations(c):
        return np.real(np.diag(modes.T @ c @ modes))

    def interference(c):
        rows = []
        for i in range(3):
            block = c[i:i + 2, i:i + 2]
            both = np.linalg.det(block).real
            single = np.trace(block).real - 2 * both
            r = block[0, 1].real
            probs = [1 - single - both, single / 2 + r, single / 2 - r, both]
            rows.append(rng.multinomial(shots, probs).tolist())
        return rows

    calibration, reference_counts, replaced = (
        rng.binomial(shots, occupations(c)) for c in (baseline, baseline, excited)
    )
    bonds_reference, bonds_replaced = interference(baseline), interference(excited)
    # Twelve binary frequencies and 24 multinomial frequencies enter the common event.
    radius = np.sqrt(np.log(2 * 36 / 0.01) / (2 * shots))
    cal = calibration / shots
    delta = (replaced - reference_counts) / shots

    def log_odds(p):
        return np.log1p(-p) - np.log(p)

    lower, upper = log_odds(cal + radius), log_odds(cal - radius)
    products = np.stack([u * v for u in (lower, upper) for v in (delta - 2 * radius, delta + 2 * radius)])
    energy_delta = (np.array(bonds_replaced) - np.array(bonds_reference)) / shots
    energy = np.pi * f @ (energy_delta[:, 2] - energy_delta[:, 1])
    energy_halfwidth = 4 * np.pi * f.sum() * radius
    return {
        "seed": seed, "shots_per_setting": shots, "settings": 18,
        "reference_mode_calibration_occupied": calibration.tolist(),
        "reference_modes_occupied": reference_counts.tolist(), "replaced_modes_occupied": replaced.tolist(),
        "bond_order": ["empty", "symmetric", "antisymmetric", "double"],
        "reference_bonds": bonds_reference, "replaced_bonds": bonds_replaced,
        "simultaneous_confidence": 0.99, "frequency_radius": radius,
        "entropy_read": log_odds(cal) @ delta,
        "entropy_interval": [products.min(axis=0).sum(), products.max(axis=0).sum()],
        "weighted_energy_read": energy,
        "weighted_energy_interval": [energy - energy_halfwidth, energy + energy_halfwidth],
    }


def all_open(board, route):
    ring = set(route)
    packed = set(range(len(route) // 2))
    incident = sum(a in ring or b in ring for a, b in board.seams)
    walls = sum((a in packed) != (b in packed) for a, b in board.seams)
    deficit = sum(6 - int(board.degree[a]) for a in packed)
    receiver = next(
        q for q in range(board.n)
        if q not in ring and len(ring.intersection(board.neighbours[q])) == 1
    )
    source, = ring.intersection(board.neighbours[receiver])
    k = candidate.KAPPA
    return {
        "incident_seams": incident, "sea_energy_lower_bound": k * incident / 2,
        "same_filling_comparator_port_range": [min(packed), max(packed)],
        "comparator_domain_walls": walls, "comparator_degree_deficit": deficit,
        "same_filling_comparator_energy": (k + 0.5) * walls + k * deficit,
        "outside_receiver": receiver, "unique_ring_neighbour": source,
        "outside_receiver_quadratic_coefficient": 1 / 8,
    }


def run():
    board = federation(2)
    route, _ = select_route(board)
    n = len(route)
    if len(set(route)) != n or any(
        b not in board.neighbours[a] for a, b in zip(route, route[1:] + route[:1])
    ):
        raise ArithmeticError("The source includes a repeated station or a nonlocal bond")
    momenta, target, filled = ring_modes(n)
    phase, count = compiled_modes(target)
    sources = tangents(momenta, phase, filled)

    @cache
    def interval(start, end):
        return interval_read(n, start, end, phase, filled, sources)

    inversions = []
    for end in (51, 55):
        for name, (hole, particle) in sources.items():
            direct = -np.real(
                particle[end] * particle[end - 1].conj() - hole[end] * hole[end - 1].conj()
            )
            for length in SIZES:
                start = end - length
                energy = inverse_stress(n, start, end, interval, "weighted_energy", name)
                entropy = inverse_stress(n, start, end, interval, "entropy", name)
                inversions.append({
                    "left_endpoint": start,
                    "right_bond": end,
                    "source": name,
                    "bond_energy": direct,
                    "energy_inversion": energy,
                    "entropy_inversion": entropy,
                    "entropy_relative_difference": abs(entropy - direct) / abs(direct),
                })
    return {
        "contact_controls": native_contact_checks(),
        "route_stations": route,
        "ring_size": n,
        "all_open": all_open(board, route),
        "half_filling": int(filled.sum()),
        "adjacent_source_rotations": count,
        "source_orbital_maximum_error": np.max(abs(phase - target)),
        "fermi_velocity": 1,
        "intervals": [
            interval(centre - length // 2, centre + length // 2)
            for length in SIZES for centre in CENTRES
        ],
        "local_stress": inversions,
        "count_bank": count_bank(n, phase, filled, sources),
    }


def stable(value):
    if isinstance(value, dict):
        return {k: stable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [stable(x) for x in value]
    if isinstance(value, (float, np.floating)):
        if abs(value) < 1e-24:
            return 0.0
        return float(format(value, ".7g"))
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = stable(run())
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        if data["route_stations"] != saved["route_stations"]:
            raise ArithmeticError("Route reconstruction disagrees with stored stations")
        if data["adjacent_source_rotations"] != saved["adjacent_source_rotations"]:
            raise ArithmeticError("Compiled source schedule differs from the receipt")

        def compare(actual, old):
            if isinstance(actual, dict):
                if actual.keys() != old.keys():
                    raise ArithmeticError("Receipt fields differ")
                for key, val in actual.items():
                    compare(val, old[key])
            elif isinstance(actual, list):
                for val, prior in zip(actual, old, strict=True):
                    compare(val, prior)
            elif isinstance(actual, str):
                if actual != old:
                    raise ArithmeticError("Receipt source labels differ")
            else:
                np.testing.assert_allclose(actual, old, rtol=2e-5, atol=2e-9)

        compare(data["contact_controls"], saved["contact_controls"])
        compare(data["all_open"], saved["all_open"])
        compare(data["local_stress"], saved["local_stress"])
        compare(data["count_bank"], saved["count_bank"])
        for row, original in zip(data["intervals"], saved["intervals"], strict=True):
            if row["stations"] != original["stations"]:
                raise ArithmeticError("Aperture endpoints disagree with the receipt")
            compare(row, original)
    print("Verified the native ring, nested interval records and held energy responses.")


if __name__ == "__main__":
    main()
