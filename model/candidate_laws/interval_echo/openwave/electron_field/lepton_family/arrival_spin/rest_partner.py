"""Search for local stationary partners at the signed family energies."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    winding_turn as wind,
)

guide = wind.circ.tree.guide


def geometries(board, route, incident=None):
    """All eligible two-face traps adjacent to the ten interior guide ports."""
    off_guide = set(range(board.n)) - set(route)
    if incident is None:
        incident = face_index(board)
    for port in route[1:-1]:
        for hub in sorted(set(board.neighbours[port]) & off_guide):
            if board.degree[hub] != 5:
                continue
            for a, b in combinations(sorted(incident.get(hub, [])), 2):
                if set(a).intersection(b) == {hub} and set(a).union(b) <= off_guide:
                    yield (port, hub, a, b)


def face_index(board):
    index = {}
    for face in board.faces():
        for p in face:
            index.setdefault(p, []).append(tuple(sorted(face)))
    return index


def federation_census(masses):
    counts = []
    classes = {}
    for level in range(3):
        b = guide.federation(level)
        heads = face_index(b)
        graph = nx.Graph(b.seams).subgraph([p for p in range(b.n) if b.degree[p] == 6])
        routes = 0
        placements = 0
        occupied = set()
        for component in nx.connected_components(graph):
            ends = sorted(p for p in component if graph.degree[p] == 1)
            if len(ends) != 2:
                raise ArithmeticError("Charged degree-six component is not a path")
            route = tuple(nx.shortest_path(graph, *ends))
            routes += 1
            for fixture in geometries(b, route, heads):
                p, h, a, c = fixture
                leaf_degrees = tuple(
                    sorted(
                        tuple(sorted(int(b.degree[i]) for i in f if i != h))
                        for f in (a, c)
                    )
                )
                at = route.index(p)
                # Reversing the guide, swapping the two faces, or swapping either pair of tips
                # relabels the full kinetic graph. The vertex cost depends only on physical degree.
                key = (len(route), min(at, len(route) - 1 - at), leaf_degrees)
                measured = spectrum(Patch(b, route, fixture), masses)
                bucket = classes.setdefault(
                    key,
                    {
                        "guide_length": key[0],
                        "least_tip_distance": key[1],
                        "face_tip_degrees": key[2],
                        "placements": 0,
                        "representative_cells": b.cells,
                        "representative_route": route,
                        "representative_fixture": fixture,
                        **measured,
                    },
                )
                for name in (
                    "vertices",
                    "incidence_rows",
                    "proper_cell_rotations_preserving_hub_and_faces",
                    "left_kernel_dimension",
                ):
                    if bucket[name] != measured[name]:
                        raise ArithmeticError("Mismatch within a weighted geometry")
                np.testing.assert_allclose(
                    measured["kinetic_floor"],
                    bucket["kinetic_floor"],
                    atol=1e-12,
                    rtol=0,
                )
                bucket["placements"] += 1
                placements += 1
                occupied.add(p)
        counts.append(
            {
                "cells": b.cells,
                "maximal_degree_six_guides": routes,
                "interior_guide_ports_with_a_placement": len(occupied),
                "placements": placements,
            }
        )
    return {
        "boards": counts,
        "spectra_by_weighted_geometry": [classes[k] for k in sorted(classes)],
    }


class Patch:
    def __init__(self, board, route, fixture):
        port, hub, a, b = fixture
        self.hub, self.port = hub, port
        self.board = board
        origin = 12 * (hub // 12)
        shapes = {tuple(sorted(v - origin for v in face)) for face in (a, b)}
        self.frame_stabilizer = sum(
            turn[hub - origin] == hub - origin
            and {tuple(sorted(turn[v] for v in f)) for f in shapes} == shapes
            for turn in wind.circ.tree.rotations()
        )
        faces = (a, b)
        bridge = wind.edge(port, hub)
        edges = sorted(
            set(wind.along(route) + [bridge] + wind.face_edges(a) + wind.face_edges(b))
        )
        bits = {e: 1 << i for i, e in enumerate(edges)}
        face_masks = []
        for face in faces:
            mask = sum(bits[e] for e in wind.face_edges(face))
            cw = sum(
                bits[wind.edge(x, y)]
                for x, y in zip(face, (*face[1:], face[0]))
                if y < x
            )
            face_masks.append((mask, cw))
        middle = route.index(port)
        parents = [
            (hub, port),
            *zip(route[middle:], route[middle + 1 :]),
            *zip(route[middle:0:-1], route[middle - 1 :: -1]),
        ]
        # On a vacancy the written tree points away from the parked load.
        base = sum(bits[wind.edge(x, y)] for x, y in parents if y < x)
        x, y = face_masks
        states = [
            (
                hub,
                base
                | (x[1] if not u else x[0] ^ x[1])
                | (y[1] if not v else y[0] ^ y[1]),
            )
            for u, v in product((0, 1), repeat=2)
        ]
        index = {s: i for i, s in enumerate(states)}
        transitions = set()
        for i, (p, w) in enumerate(states):
            moves = []
            for e, bit in bits.items():
                x, y = e
                if p in e and ((p == x) != bool(w & bit)):
                    moves.append(
                        (
                            (y if p == x else x, w ^ bit),
                            "bridge" if e == bridge else "matter",
                        )
                    )
            for k, (mask, cw) in enumerate(face_masks):
                if w & mask in (cw, mask ^ cw):
                    moves.append(((p, w ^ mask), f"ring{k}"))
            for other, kind in moves:
                if other not in index:
                    index[other] = len(states)
                    states.append(other)
                j = index[other]
                transitions.add((kind, min(i, j), max(i, j)))
        self.n = len(states)
        rows, tags, anchors = [], [], []
        physical = sorted({x for e in edges for x in e})

        def inward(word):
            return np.array(
                [
                    sum((p == (x if word & bit else y)) for (x, y), bit in bits.items())
                    for p in physical
                ]
            )

        reference = inward(states[0][1]) + (np.array(physical) == hub)
        for p, w in states:
            np.testing.assert_array_equal(
                inward(w) - reference, -(np.array(physical) == p).astype(int)
            )
        for kind, i, j in sorted(transitions):
            p, q = states[i][0], states[j][0]
            if (
                kind == "bridge"
                and p != hub
                or p in route
                and q in route
                and route.index(p) > route.index(q)
            ):
                i, j = j, i
            row = np.zeros(self.n)
            rate = 0.5 if kind.startswith("ring") else guide.HOP
            row[i], row[j] = -np.sqrt(rate), np.sqrt(rate)
            rows.append(row)
            tags.append(kind)
            anchors.append(i)
        for i, (p, w) in enumerate(states):
            cost = guide.DEGREE_FIVE_COST * (6 - int(board.degree[p]))
            if cost:
                row = np.zeros(self.n)
                row[i] = np.sqrt(cost)
                rows.append(row)
                tags.append("cost")
                anchors.append(i)
        self.b = np.array(rows)
        self.tags = np.array(tags)
        self.ports = np.array([states[i][0] for i in [*range(self.n), *anchors]])
        self.region = np.array([fixture[1], *sorted((set(a) | set(b)) - {hub})])
        # Compare the independent boundary head with the load in the complete five-port region.
        sharp_inside = np.isin(self.ports, self.region)
        bw = np.array(
            [
                int(
                    (bridge[0] if states[i][1] & bits[bridge] else bridge[1])
                    in self.region
                )
                for i in [*range(self.n), *anchors]
            ]
        )
        np.testing.assert_array_equal(bw - bw[0] - 1, -sharp_inside.astype(int))
        # The reference charge was at the hub; its initial bridge points toward the guide.
        self.bridge_flux = bw
        reachable = {0}
        for _ in range(len(states)):
            size = len(reachable)
            for _, i, j in transitions:
                if i in reachable or j in reachable:
                    reachable.update((i, j))
            if len(reachable) == size:
                break
        if len(reachable) != len(states):
            raise ArithmeticError("Seed words split into disconnected matter sectors")
        if not (len(states) == len(index) and any(t == "bridge" for t in tags)):
            raise ArithmeticError("Incomplete connected charge sector")

    def local_space(self):
        chi = np.zeros(self.n)
        chi[:4] = [0.5, -0.5, -0.5, 0.5]
        columns = np.zeros((len(self.b) + self.n, 4))
        columns[: self.n, 0] = chi
        weights = []
        for j, kind in enumerate(("ring0", "ring1", "cost"), 1):
            raw = (self.b @ chi) * (self.tags == kind)
            norm = np.linalg.norm(raw)
            columns[self.n :, j] = raw / norm
            weights.append(norm)
            np.testing.assert_allclose(self.b.T @ (raw / norm), norm * chi, atol=2e-14)
        np.testing.assert_allclose(weights, [1, 1, np.sqrt(8)], atol=3e-14)
        np.testing.assert_allclose(columns.T @ columns, np.eye(4), atol=2e-14)
        if np.any(columns[self.ports != self.hub]):
            raise ArithmeticError(
                "A local column carries the charge outside the resting port"
            )
        normal = np.array(weights) / np.linalg.norm(weights)
        e1 = np.array([1.0, 0, 0])
        bright = e1 - normal * (e1 @ normal)
        bright /= np.linalg.norm(bright)
        dark = np.cross(normal, bright)
        pair = np.vstack((np.zeros(2), np.column_stack((bright, dark))))
        np.testing.assert_allclose(
            self.b.T @ (columns[self.n :, 1:] @ pair[1:]), 0, atol=2e-14
        )
        return columns, pair


def reduced(mass, rings):
    return guide.generator(
        np.array([[rings[0]], [rings[1]], [np.sqrt(8)]]), np.array([[mass]])
    )


def spectrum(patch, masses):
    columns, pair = patch.local_space()
    states = columns @ pair
    residuals = []
    for m in masses:
        hv = np.vstack(
            (
                m * states[: patch.n] + patch.b.T @ states[patch.n :],
                patch.b @ states[: patch.n] - m * states[patch.n :],
            )
        )
        residuals.append(float(np.max(abs(hv + m * states))))
    if max(residuals) > 2e-13:
        raise ArithmeticError(
            "Stationary signed family pair fails the full Hamiltonian"
        )
    kinetic = eigh(patch.b.T @ patch.b, eigvals_only=True)
    if kinetic[0] <= 0:
        raise ArithmeticError(
            "Connected sector visiting a penalized port must have positive kinetic floor"
        )
    return {
        "vertices": patch.n,
        "incidence_rows": len(patch.b),
        "left_kernel_dimension": len(patch.b) - patch.n,
        "proper_cell_rotations_preserving_hub_and_faces": patch.frame_stabilizer,
        "lowest_positive_energies": np.sqrt(masses * masses + kinetic[0]).tolist(),
        "kinetic_floor": float(kinetic[0]),
        "maximum_eigen_residual_by_character": residuals,
    }


def hold(mass, rings, time):
    """Exact four-column evolution with a fixed binary setting of the two face controls."""
    b = np.array([*rings, np.sqrt(8)])
    e = np.sqrt(mass**2 + b @ b)
    normal = b / np.linalg.norm(b)
    c, s = np.cos(e * time), np.sin(e * time) / e
    out = np.zeros((4, 4), complex)
    out[0, 0] = c - 1j * mass * s
    out[0, 1:] = out[1:, 0] = -1j * s * b
    idle = np.exp(1j * mass * time)
    out[1:, 1:] = idle * np.eye(3) + (c + 1j * mass * s - idle) * np.outer(
        normal, normal
    )
    return out


def control_loop(mass, duration, ticks):
    # On each side of the square only one face varies. The actual pulses are binary.
    corners = ((1, 1), (0, 1), (0, 0), (1, 0), (1, 1))
    out = np.eye(4, dtype=complex)
    for start, end in zip(corners, corners[1:]):
        for tick in range(ticks):
            share = (1 - np.cos(np.pi * (tick + 0.5) / ticks)) / 2
            step = duration / ticks
            first = hold(mass, start, step * (1 - share) / 2)
            out = first @ hold(mass, end, step * share) @ first @ out
    np.testing.assert_allclose(out.conj().T @ out, np.eye(4), atol=1e-11)
    return out


def guide_control(patch, route, masses):
    # Restrict to the charge on the original guide with the four parked face words.
    vertices = np.flatnonzero(np.isin(patch.ports[: patch.n], route))
    if len(vertices) != 4 * len(route):
        raise ArithmeticError(
            "The guide must contain four spectator words at each position"
        )
    path = np.count_nonzero(patch.b[:, vertices], axis=1) == 2
    results = []
    for rate in (0, 1, 2):
        enabled = path & (
            (patch.tags == "matter")
            | (patch.tags == "ring0") & (rate >= 1)
            | (patch.tags == "ring1") & (rate >= 2)
        )
        b = patch.b[np.ix_(enabled, vertices)]
        levels = np.linalg.eigvalsh(b.T @ b)
        one_path = guide.path_incidence(patch.board, route)
        mu = np.linalg.eigvalsh(one_path.T @ one_path)
        field_shifts = [
            int(rate >= 1) * u + int(rate >= 2) * v
            for u, v in product((0, 1), repeat=2)
        ]
        np.testing.assert_allclose(
            levels, np.sort(np.concatenate([mu + s for s in field_shifts])), atol=4e-14
        )
        results.append(
            {
                "faces_enabled": rate,
                "positive_mass_level_multiplicity": int(
                    np.count_nonzero(abs(levels) < 1e-11)
                ),
                "uniform_guide_face_mode_shifts": field_shifts,
            }
        )
    shifted = np.sqrt(masses**2 + 1)
    return {
        "bridge_open": False,
        "guide_vertices": len(vertices),
        "controls": results,
        "uniform_positive_two_face_doublet_energies": shifted.tolist(),
        "uniform_positive_two_face_doublet_ratios": (shifted / shifted[0]).tolist(),
    }


def prepare_with_flux(patch, mass, columns, pair, bridge_u, cycle, decoder, rng):
    """A local source starts in the vertex grade. Keep every failed downstream charge trial."""
    local_projection = (patch.ports == patch.hub)[:, None] * (bridge_u @ columns)
    no_click = columns.T @ local_projection
    np.testing.assert_allclose(local_projection, columns @ no_click, atol=3e-13)
    source = decoder @ np.eye(4)[:, 0]
    quarter = hold(mass, (1, 1), np.pi / (2 * np.sqrt(10 + mass**2)))
    accepted = source
    failures = []
    for _ in range(32):
        prev = float(np.vdot(accepted, accepted).real)
        accepted = no_click @ quarter @ accepted
        failures.append(max(0.0, prev - float(np.vdot(accepted, accepted).real)))
    success = float(np.vdot(accepted, accepted).real)
    fidelity = float(abs(np.vdot(pair[:, 0], accepted)) ** 2 / success)
    outcomes = []
    for repeats in (0, 7, 14):
        wave = (
            bridge_u
            @ columns
            @ decoder
            @ np.linalg.matrix_power(cycle, repeats)
            @ accepted
        )
        click = float(np.linalg.norm(wave[patch.ports == patch.port]) ** 2)
        trials = [*failures, click, max(0.0, success - click)]
        outcomes.append(
            {
                "cycles": repeats,
                "decoder_click_per_enrollment": click,
                "counts_failure_rounds_then_click_then_no_click": wind.sample(
                    trials, 2**20, rng
                ),
            }
        )
    return {
        "success_from_vertex_source": success,
        "bright_fidelity_given_success": min(1.0, fidelity),
        "failure_probabilities_by_round": failures,
        "shots_per_setting": 2**20,
        "postpreparation_reads": outcomes,
    }


def relative_energy(patch, route, mass, columns, pair, rng):
    """Two quadratures against the positive guide level on the same bridged-off patch."""
    n = patch.n
    b = patch.b * (patch.tags != "bridge")[:, None]
    ref = np.zeros(n + len(b))
    inds = np.flatnonzero(np.isin(patch.ports[:n], route))
    ref[inds] = 1 / np.sqrt(len(inds))
    target = columns @ pair[:, 0]
    h = guide.generator(b, np.array([[mass]]))
    np.testing.assert_allclose(h @ ref, mass * ref, atol=2e-14)
    np.testing.assert_allclose(h @ target, -mass * target, atol=2e-14)
    levels, eig = eigh(h)
    state = (ref + target) / np.sqrt(2)
    # Fix the analyzer convention: phase is target amplitude divided by reference amplitude.
    coeff = eig.T @ state
    amps = np.array([ref, target]) @ eig
    readings = []
    previous = 0.0
    for power in range(15):
        t = 2**power
        a, d = amps @ (np.exp(-1j * t * levels) * coeff)
        fringe = 2 * a.conjugate() * d
        pairs = [
            wind.sample([(1 + v) / 2, (1 - v) / 2], 2**22, rng)
            for v in (fringe.real, fringe.imag)
        ]
        z = complex(*[(plus - minus) / (plus + minus) for plus, minus in pairs])
        phase = np.angle(z)
        # Each step doubles the wait. Resolve the wrap from the preceding observed phase.
        unwrapped = phase + 2 * np.pi * round((2 * previous - phase) / (2 * np.pi))
        previous = unwrapped
        readings.append(
            {"time": t, "quadrature_counts": pairs, "unwrapped_phase": float(unwrapped)}
        )
    return {
        "shots_per_quadrature": 2**22,
        "estimated_target_minus_guide_frequency": float(-previous / t),
        "readings": readings,
    }


def observer(patch, masses, route):
    columns, pair = patch.local_space()
    rows = []
    rng = np.random.default_rng(191183)
    for mass in masses:
        for setting in product((0, 1), repeat=2):
            enabled = (
                (patch.tags == "cost")
                | ((patch.tags == "ring0") & bool(setting[0]))
                | ((patch.tags == "ring1") & bool(setting[1]))
            )
            h = guide.generator(patch.b * enabled[:, None], np.array([[mass]]))
            np.testing.assert_allclose(
                h @ columns, columns @ reduced(mass, setting), atol=2e-14
            )
            np.testing.assert_allclose(
                hold(mass, setting, 0.375),
                expm(-0.375j * reduced(mass, setting)),
                atol=2e-14,
            )
        cycle = control_loop(mass, 24, 256)
        decoder = hold(mass, (0, 1), np.pi / (2 * np.sqrt(9 + mass**2)))
        enabled = patch.tags == "bridge"
        read_h = guide.generator(patch.b * enabled[:, None], np.array([[mass]]))
        # Reconstruct the load click on the opposite physical port; suppressing the bridge is the sham.
        bridge_u = expm(-1j * np.pi / np.sqrt(8 + mass**2) * read_h)
        free = guide.generator(patch.b, np.array([[mass]]))
        eigen, vec = eigh(free)
        initial = columns @ pair
        free_reads = []
        for t in (0, 4, 8, 12, 32):
            wave = vec @ (np.exp(-1j * t * eigen)[:, None] * (vec.T @ initial))
            np.testing.assert_allclose(
                wave, np.exp(1j * mass * t) * initial, atol=3e-13
            )
            free_reads.append(
                {
                    "time": t,
                    "charge_outside_hub_by_input": np.sum(
                        abs(wave[patch.ports != patch.hub]) ** 2, axis=0
                    ).tolist(),
                }
            )
        reversed_mass = guide.generator(patch.b, np.array([[-mass]]))
        np.testing.assert_allclose(reversed_mass @ initial, mass * initial, atol=2e-14)
        settings = []
        for repeats in (0, 7, 14):
            state = np.linalg.matrix_power(cycle, repeats) @ pair
            compressed = pair.T @ state
            gram = compressed.conj().T @ compressed
            max_leak = max(0.0, float(1 - np.linalg.eigvalsh(gram)[0]))
            wrong_port = (
                np.linalg.norm((columns @ state)[patch.ports != patch.hub]) ** 2
            )
            if wrong_port > 1e-25:
                raise ArithmeticError("Rotation moved the parked charge")
            reads = []
            for readout in (decoder, np.eye(4)):
                decoded = columns @ readout @ state
                wave = bridge_u @ decoded
                select = patch.ports == patch.port
                flux_leaves = patch.bridge_flux == 1
                np.testing.assert_array_equal(
                    flux_leaves, ~np.isin(patch.ports, patch.region)
                )
                # Only the bridge is available during release, so the two detectors must agree on the launched state.
                clicks = np.sum(abs(wave[select]) ** 2, axis=0)
                np.testing.assert_allclose(
                    clicks, np.sum(abs(wave[flux_leaves]) ** 2, axis=0), atol=2e-13
                )
                reads.append(
                    {
                        "far_port_by_input": clicks.tolist(),
                        "no_click": (1 - clicks).tolist(),
                        "cut_bridge_by_input": np.sum(
                            abs(decoded[select]) ** 2, axis=0
                        ).tolist(),
                    }
                )
            settings.append(
                {
                    "cycles": repeats,
                    "dark_input_to_bright_state": float(abs(compressed[0, 1]) ** 2),
                    "bright_input_to_dark_state": float(abs(compressed[1, 0]) ** 2),
                    "worst_rest_space_leakage": max_leak,
                    "charge_outside_hub_during_controls": float(wrong_port),
                    "weight_in_positive_energy_space_after_rotation": np.sum(
                        abs(vec[:, eigen > 0].T @ (columns @ state)) ** 2, axis=0
                    ).tolist(),
                    "decoded": reads[0],
                    "no_decode": reads[1],
                }
            )
        # The parked degree-six bridge recipient has the same face bits, but no local cost row.
        chi = columns[: patch.n, 0]
        bridge_rows = patch.b[patch.tags == "bridge"]
        mate = -(bridge_rows.T @ bridge_rows @ chi) / guide.HOP
        mate[patch.ports[: patch.n] != patch.port] = 0
        np.testing.assert_allclose(mate @ mate, 1.0, atol=1e-14)
        strengths = [
            float(np.linalg.norm((patch.b @ mate) * (patch.tags == kind)))
            for kind in ("ring0", "ring1", "cost")
        ]
        np.testing.assert_allclose(strengths, [1, 1, 0], atol=2e-14)
        rows.append(
            {
                "signed_stationary_energy": -float(mass),
                "loop": "four sides, 24 native time units per side, 256 symmetric binary ticks per side",
                "displaced_local_face_and_cost_strengths": strengths,
                "settings": settings,
                "free_charge_reads": free_reads,
                "reversed_mass_control_energy": float(mass),
                "heralded_preparation": prepare_with_flux(
                    patch, mass, columns, pair, bridge_u, cycle, decoder, rng
                ),
                "relative_energy_read": relative_energy(
                    patch, route, mass, columns, pair, rng
                ),
            }
        )
    return rows


def run():
    board, route = guide.longest_unpinned_route()
    masses = guide.register_effect()[1]
    fixtures = sorted(geometries(board, route))
    if not fixtures:
        raise ArithmeticError("No admissible placement")
    census = []
    for fixture in fixtures:
        patch = Patch(board, route, fixture)
        census.append(
            {
                "guide_port": fixture[0],
                "hub": fixture[1],
                "faces": fixture[2:],
                **spectrum(patch, masses),
            }
        )
    primary = Patch(board, route, fixtures[0])
    return {
        "route": route,
        "family_effects": masses.tolist(),
        "magnitude_ratios": (masses / masses[0]).tolist(),
        "rest_pair_eigenvalues": (-masses).tolist(),
        "placements": census,
        "all_maximal_guides": federation_census(masses),
        "guide_control": guide_control(primary, route, masses),
        "local_rotation_and_charge_read": observer(primary, masses, route),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
