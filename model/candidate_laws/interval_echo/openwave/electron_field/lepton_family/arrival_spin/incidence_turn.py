"""Repeat the winding order test with the first-order flight generator on each contact."""

import argparse
from itertools import product
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    winding_turn as winding,
)


def contact_gate(mass):
    weight = winding.circ.tree.guide.HOP
    incidence = np.sqrt(weight) * np.array([[-1, 1]])
    h = winding.circ.tree.guide.generator(incidence, np.array([[mass]]))
    energy = np.sqrt(mass * mass + 2 * weight)
    wait = np.pi / energy
    # Balance the two arms with idle evolution. Every idle vertex accrues exp(-i m t).
    balanced = np.exp(1j * mass * wait) * expm(-1j * wait * h)
    relative = -np.exp(1j * mass * wait)
    stay, move = (1 + relative) / 2, (1 - relative) / 2
    np.testing.assert_allclose(
        balanced[:2, :2], [[stay, move], [move, stay]], atol=1e-14
    )
    np.testing.assert_allclose(balanced[2, :2], 0, atol=1e-14)
    # Flipping the sign convention of the incidence grade leaves the two-site propagator intact.
    h_reversed = winding.circ.tree.guide.generator(-incidence, np.array([[mass]]))
    np.testing.assert_allclose(
        expm(-1j * wait * h_reversed)[:2, :2], expm(-1j * wait * h)[:2, :2], atol=1e-14
    )
    wrong_clock = expm(-1j * np.pi / (2 * weight) * h)[:, 0]
    return (
        stay,
        move,
        {
            "first_auxiliary_return_wait": float(wait),
            "balanced_contact_operator_distance_from_ideal_swap": float(
                abs(1 + relative)
            ),
            "incidence_return_probabilities_from_sharp_endpoint": (
                abs(balanced[:, 0]) ** 2
            ).tolist(),
            "incidence_probabilities_at_laplacian_swap_wait": (
                abs(wrong_clock) ** 2
            ).tolist(),
            "contact_outcomes": ["initial site", "opposite site", "auxiliary grade"],
        },
    )


def contact_sector(detector):
    # Retain misrouted amplitudes at every step. Each enabled contact acts on disjoint pairs.
    sharp = [
        detector.state(f[0], winds)
        for f in detector.faces
        for winds in product((0, 1), repeat=4)
    ]
    indices = {s: i for i, s in enumerate(sharp)}
    for state in sharp:
        for e in detector.bits:
            other = detector.pulse(state, e)
            if other not in indices:
                indices[other] = len(sharp)
                sharp.append(other)
    operators = {}
    for e in detector.bits:
        perm = np.array([indices[detector.pulse(state, e)] for state in sharp])
        np.testing.assert_array_equal(perm[perm], np.arange(len(perm)))
        operators[e] = perm
    return indices, operators


def propagate_order(source, sequence, ops, stay, move):
    wave = source.copy()
    for e in sequence:
        wave = stay * wave + move * wave[ops[e]]
    return wave


def order_record(arms, k):
    # The two source columns label the source reference. Kraus rows retain unmeasured exterior heads.
    outputs = [np.einsum("nw,ewr->ner", op, k) for op in arms]
    constructive = float(np.linalg.norm((outputs[0] + outputs[1]) / 2) ** 2 / 2)
    destructive = float(np.linalg.norm((outputs[0] - outputs[1]) / 2) ** 2 / 2)
    capture = float(np.linalg.norm(k) ** 2 / 2)
    np.testing.assert_allclose(constructive + destructive, capture, atol=1e-14)
    return {
        "capture_probability": capture,
        "constructive_per_enrollment": constructive,
        "destructive_per_enrollment": destructive,
        "destructive_given_capture": destructive / capture,
    }


def compare_axes(row, sign, spur, k, gates):
    faces = [tuple(f) for f in row["faces"]]
    detector = winding.FaceDetector(faces, sign, spur)
    index, ops = contact_sector(detector)
    a, b = map(tuple, row["axes_as_port_permutations"])
    start = faces[0][0]
    seqs = [
        detector.turn(start, x)[0] + detector.turn(x[start], y)[0]
        for x, y in ((a, b), (b, a))
    ]
    source = np.zeros((len(index), 2), complex)
    for i in (0, 1):
        source[index[detector.state(start, (i, 0, 0, 0))], i] = 1
    ideal = [propagate_order(source, seq, ops, 0, 1) for seq in seqs]
    np.testing.assert_array_equal(ideal[0], ideal[1])
    for arm, labels in zip(ideal, row["enrolled_arm_port_and_head_words"], strict=True):
        np.testing.assert_array_equal(
            arm.argmax(axis=0), [index[tuple(s)] for s in labels]
        )
    record = {}
    for label, (stay, move, diagnostics) in gates.items():
        arms = [propagate_order(source, seq, ops, stay, move) for seq in seqs]
        distances = [
            float(np.linalg.norm(op - ref, ord=2))
            for op, ref in zip(arms, ideal, strict=True)
        ]
        delta = diagnostics["balanced_contact_operator_distance_from_ideal_swap"]
        # Telescoping the unitary products gives a worst-case bound on every enrolled state.
        for distance, seq in zip(distances, seqs, strict=True):
            if distance > len(seq) * delta + 1e-13:
                raise ArithmeticError(
                    "Order error exceeds the full-sector operator bound"
                )
        result = order_record(arms, k)
        bound = min(1.0, (sum(map(len, seqs)) * delta / 2) ** 2)
        if result["destructive_given_capture"] > bound + 1e-14:
            raise ArithmeticError("Order fringe exceeds the calibrated all-input bound")
        record[label] = {
            **result,
            "operator_error_on_two_prepared_words_by_order": distances,
            "destructive_bound_for_every_state_with_commuting_ideal_orders": float(
                bound
            ),
            "supplied_anticommuting_spin_destructive_given_capture": result[
                "constructive_per_enrollment"
            ]
            / result["capture_probability"],
        }
    return {
        "axes_as_port_permutations": [a, b],
        "reachable_sharp_words": len(index),
        "pulse_counts_by_order": list(map(len, seqs)),
        "variants": record,
    }


def local_rest_records(flight, mass, region):
    """Compress all local matrix units to the rested guide copies, keeping exterior heads blind."""
    links = winding.along(flight.route)
    enabled = [link in links for link in flight.links]
    graph = nx.Graph()
    for row in flight.bonds[enabled]:
        graph.add_edge(*np.flatnonzero(row))
    copies = list(nx.connected_components(graph))
    h = winding.full_generator(flight, mass, enabled)
    mask = sum(bit for e, bit in flight.bits.items() if any(p in region for p in e))
    records = []
    for j, component in enumerate(copies):
        ports = [flight.ports[i] for i in component]
        if len(ports) != len(flight.route) or set(ports) != set(flight.route):
            raise ArithmeticError("A guide copy has different port support")
        rest = np.zeros(len(h))
        rest[list(component)] = (
            1  # common unnormalized vector; its squared norm is the route length
        )
        np.testing.assert_array_equal(h @ rest, mass * rest)
        for i in component:
            p, word = flight.states[i]
            if p in region:
                records.append((j, (p, word & mask), word & ~mask))
    inside = sorted({local for _, local, _ in records})
    exterior = sorted({ext for _, _, ext in records})
    tensor = np.zeros((len(inside), len(exterior), len(copies)), dtype=int)
    for j, local, ext in records:
        tensor[inside.index(local), exterior.index(ext), j] += 1
    for left, right in product(tensor, repeat=2):
        compression = left.T @ right
        np.testing.assert_array_equal(
            compression, compression[0, 0] * np.eye(len(copies), dtype=int)
        )
    return {
        "region": region,
        "rest_copies": len(copies),
        "local_port_and_incident_head_records": len(inside),
        "exterior_words_at_capture": len(exterior),
        "scalar_matrix_units_checked": len(inside) ** 2,
        "capture_probability_per_copy": len(region) / len(flight.route),
    }


def run():
    guide = winding.circ.tree.guide
    board, route = guide.longest_unpinned_route()
    effect, masses = guide.register_effect()
    bare = guide.generator(guide.path_incidence(board, route), effect)
    levels = np.linalg.eigvalsh(bare)
    multiplicities = [
        int(np.count_nonzero(abs(levels - mass) < 1e-10)) for mass in masses
    ]
    gates = {
        "light_character": contact_gate(float(masses[0])),
        "zero_mass_control": contact_gate(0.0),
    }
    # Reconstruct the input channel from the source and guide, independent of the archived fringe.
    archived = json.loads((HERE / "winding_turn.json").read_text())
    rows = []
    for sector in archived["sectors"]:
        spur = tuple(sector["receiver_spur"])
        sign = sector["sign"]
        source = winding.Flight(board, route, sign, True, spur)
        guide_read = winding.staged(source, float(masses[0]), board)
        read = next(r for r in guide_read["reads"] if r["time"] == 8)
        k, _ = winding.filters_and_flux(read)
        local_spur = tuple(p - route[-1] // 12 * 12 for p in spur)
        axes = [
            compare_axes(row, sign, local_spur, k, gates)
            for row in sector["eligible_turns"]
            if row["eligible"]
        ]
        rows.append(
            {
                "sign": sign,
                "axes": axes,
                "interior_rest_measurements": [
                    local_rest_records(source, float(masses[0]), region)
                    for region in (route[5:6], route[5:8], route[1:-1])
                ],
            }
        )
    return {
        "rest_energy_multiplicities_on_one_undressed_guide": multiplicities,
        "family_rest_energies": masses.tolist(),
        "single_contact": {key: values[2] for key, values in gates.items()},
        "flight_exposure": 8,
        "sectors": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
