"""Turn a flux-captured charge and a directed face with addressed material contacts."""

import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random
import sys

import networkx as nx
import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    circulation as circ,
    tomography,
)


def edge(a, b):
    return tuple(sorted((a, b)))


def along(vertices):
    return [edge(a, b) for a, b in zip(vertices, vertices[1:])]


def face_edges(face):
    return along((*face, face[0]))


class Flight:
    """Incidence lift of the reachable two-loop written graph in one mass character."""

    def __init__(self, board, route, sign, spur, receiver_spur=None):
        fixture = circ.tree.WrittenTree(board, route, sign)
        model = circ.Circulation(fixture, board, spur)
        enabled = along(route)
        for _, _, face in model.loops:
            enabled.extend(face_edges(face))
        if spur:
            enabled.append(model.spur)
        states = model.states[:4].copy()
        if receiver_spur is not None:
            tip, leaf = receiver_spur
            new = edge(tip, leaf)
            if new in enabled or leaf in {x for e in enabled for x in e}:
                raise ValueError("Receiver spur must end at an unused port")
            enabled.append(new)
            bit = model.bits.setdefault(new, 1 << len(model.bits))
            states = [
                (p, (w | bit) if (tip < leaf) == (sign == 1) else (w & ~bit))
                for p, w in states
            ]
        index = {s: i for i, s in enumerate(states)}
        transitions = set()
        for i, state in enumerate(states):
            for kind, neighbors in (
                ("hop", (model.contact(state, e) for e in enabled)),
                ("ring", model.flips(state)),
            ):
                for other in neighbors:
                    if other == state:
                        continue
                    if other not in index:
                        index[other] = len(states)
                        states.append(other)
                    j = index[other]
                    transitions.add((kind, min(i, j), max(i, j)))
        addresses = sorted({p for e in enabled for p in e})
        baseline = model.inward(states[0][1], addresses)
        for p, word in states:
            np.testing.assert_array_equal(
                model.inward(word, addresses) - baseline,
                sign
                * (
                    (np.array(addresses) == p).astype(int)
                    - (np.array(addresses) == route[0]).astype(int)
                ),
            )
        self.states, self.bits = states, model.bits
        self.n = len(states)
        bonds, anchors, kinds, cuts, links = [], [], [], [], []
        for kind, i, j in sorted(transitions):
            # Use the original source-to-detector incidence sign and left anchor on the guide.
            link = edge(states[i][0], states[j][0])
            if link in along(route) and route.index(states[i][0]) > route.index(
                states[j][0]
            ):
                i, j = j, i
            anchor = i
            row = np.zeros(self.n)
            row[i], row[j] = -1, 1
            bonds.append(row * np.sqrt(circ.tree.guide.HOP if kind == "hop" else 0.5))
            anchors.append(anchor)
            kinds.append(kind)
            links.append(
                link
                if kind == "hop"
                else next(
                    face
                    for mask, _, face in model.loops
                    if states[i][1] ^ states[j][1] == mask
                )
            )
            cuts.append(link == fixture.entrance)
        for i, (p, _) in enumerate(states):
            cost = circ.tree.guide.DEGREE_FIVE_COST * (6 - int(board.degree[p]))
            if cost:
                row = np.zeros(self.n)
                row[i] = np.sqrt(cost)
                bonds.append(row)
                anchors.append(i)
                kinds.append("cost")
                cuts.append(False)
                links.append((p,))
        self.ports = np.array([states[i][0] for i in [*range(self.n), *anchors]])
        self.links = links
        self.route = route
        self.bonds = np.array(bonds)
        self.kinds = np.array(kinds)
        self.cuts = np.array(cuts)
        self.flux = np.array(
            [
                fixture.boundary_head(states[i][1]) - fixture.initial_inward == sign
                for i in [*range(self.n), *anchors]
            ]
        )
        np.testing.assert_array_equal(
            self.flux,
            [states[i][0] in fixture.receiver for i in [*range(self.n), *anchors]],
        )
        mask, cw, face = model.loops[1]
        self.face = tuple(p - fixture.base for p in face)
        self.pairs = []
        local_mask = sum(
            bit
            for (a, b), bit in model.bits.items()
            if a in fixture.receiver and b in fixture.receiver
        )
        backgrounds = sorted({w & ~mask for p, w in states if p == route[-1]})
        if len({w & local_mask for w in backgrounds}) != 1:
            raise ArithmeticError("Spectator heads depend on the admitted face state")
        for background in backgrounds:
            orientations = [cw, mask ^ cw]
            self.pairs.append([index[route[-1], background | o] for o in orientations])
        # The two kept words in each pair differ only on the receiving triangle.
        self.source = np.zeros((self.n, 2))
        self.source[:2, 0] = self.source[2:4, 1] = 2**-0.5

    def read(self, mass, cut=False, stop_ring=False):
        b = self.bonds.copy()
        if cut:
            b[self.cuts] = 0
        if stop_ring:
            b[self.kinds == "ring"] = 0
        costs, modes = np.linalg.eigh(b.T @ b)
        energies = np.sqrt(mass**2 + np.maximum(costs, 0))
        sources = modes.T @ self.source
        positive_weight = np.sum(
            abs(sources) ** 2 * ((1 + mass / energies) / 2)[:, None], axis=0
        )
        reads = []
        for time in circ.tree.TIMES:
            osc = np.sin(time * energies) / energies
            sites = modes @ (
                (np.cos(time * energies) - 1j * mass * osc)[:, None] * sources
            )
            aux = -1j * b @ (modes @ (osc[:, None] * sources))
            wave = np.vstack((sites, aux))
            np.testing.assert_allclose(wave.conj().T @ wave, np.eye(2), atol=3e-13)
            reads.append(dict(time=time, **terminal_read(self, wave)))
        return {
            "kinetic_floor": float(costs[0]),
            "lowest_positive_vertex_energy": float(energies[0]),
            "positive_energy_fraction_per_source": positive_weight.tolist(),
            "reads": reads,
        }


def routing_tree(faces, spur):
    graph = nx.Graph(circ.tree.cell().seams)
    graph.remove_edges_from(e for face in faces for e in face_edges(face))
    rng = random.Random(7451)
    for attempt in range(4096):
        if attempt:
            for a, b in sorted(graph.edges):
                graph[a][b]["weight"] = rng.random()
        tree = nx.minimum_spanning_tree(graph)
        if not nx.is_tree(tree):
            raise ArithmeticError("Disconnected face complement")
        if (
            spur
            and tree.has_edge(*spur)
            and nx.shortest_path(tree, faces[0][0], spur[1])[-2] != spur[0]
        ):
            continue
        forks = {}
        for i, j in permutations(range(4), 2):
            for start, tip in permutations(faces[i], 2):
                for end in faces[j]:
                    path = nx.shortest_path(tree, start, end)
                    if tip in path:
                        forks[i, j] = start, tip, end, path[path.index(tip) :]
                        break
                if (i, j) in forks:
                    break
        if len(forks) == 12:
            return tree, forks
    raise ArithmeticError("No spanning tree implements all directed winding gates")


class FaceDetector:
    def __init__(self, faces, sign, spur):
        self.faces = faces
        self.sign = sign
        self.tree, self.forks = routing_tree(faces, spur)
        edges = sorted(
            {edge(a, b) for a, b in self.tree.edges}
            | {e for f in faces for e in face_edges(f)}
        )
        self.bits = {e: 1 << i for i, e in enumerate(edges)}
        self.face_words = []
        for face in faces:
            forward = sum(
                self.bits[edge(a, b)]
                for a, b in zip(face, (*face[1:], face[0]))
                if b < a
            )
            mask = sum(self.bits[e] for e in face_edges(face))
            self.face_words.append((forward, mask ^ forward))

    def state(self, p, winds):
        word = sum(self.face_words[i][w] for i, w in enumerate(winds))
        for parent, child in nx.bfs_edges(self.tree, p):
            if (parent < child) == (self.sign == 1):
                word |= self.bits[edge(parent, child)]
        return p, word

    def pulse(self, state, contact):
        p, word = state
        a, b = contact
        bit = self.bits[contact]
        if p in contact and ((p == a) == bool(word & bit)) == (self.sign == 1):
            return (b if p == a else a), word ^ bit
        return state

    def apply(self, state, schedule):
        for contact in schedule:
            state = self.pulse(state, contact)
        return state

    def face_not(self, face_index, anchor):
        a, b = [x for x in self.faces[face_index] if x != anchor]
        first, second, third = edge(anchor, a), edge(a, b), edge(b, anchor)
        return [first, second, third, second, first]

    def controlled_flip(self, parking, control, target):
        anchor, tip, dest, path = self.forks[control, target]
        reach = along(nx.shortest_path(self.tree, parking, anchor))
        fork_edge = edge(anchor, tip)
        bit = self.bits[fork_edge]
        one_moves = (
            (anchor == fork_edge[0]) == bool(self.face_words[control][1] & bit)
        ) == (self.sign == 1)
        invert_control = [] if one_moves else self.face_not(control, anchor)
        branch = [fork_edge, *along(path)]
        enter = reach + invert_control + branch
        return enter + self.face_not(target, dest) + enter[::-1]

    def turn(self, p, rot):
        moved_faces = [self.faces.index(tuple(rot[x] for x in f)) for f in self.faces]
        pulses = []
        for i, j in enumerate(moved_faces):
            if i < j:
                for c, t in ((i, j), (j, i), (i, j)):
                    pulses.extend(self.controlled_flip(p, c, t))
        pulses.extend(along(nx.shortest_path(self.tree, p, rot[p])))
        return pulses, moved_faces


def rotations_for_face(face):
    all_turns = sorted(circ.tree.rotations())
    half = [
        r
        for r in all_turns
        if r != tuple(range(12)) and all(r[r[i]] == i for i in range(12))
    ]
    for a, b in combinations(half, 2):
        if all(a[b[i]] == b[a[i]] for i in range(12)):
            group = [tuple(range(12)), a, b, tuple(a[b[i]] for i in range(12))]
            yield a, b, [tuple(turn[x] for x in face) for turn in group]


def calibrate(face, sign, spur):
    rows = []
    for a, b, faces in rotations_for_face(face):
        row = {"axes_as_port_permutations": [a, b], "faces": faces}
        rows.append(row)
        if len({p for f in faces for p in f}) != 12:
            row["eligible"] = False
            continue
        detector = FaceDetector(faces, sign, spur)
        ports = [f[0] for f in faces]
        schedules = {}
        for p, rot in product(ports, (a, b)):
            schedule, moved = detector.turn(p, rot)
            schedules[p, rot] = schedule
            for winds in product((0, 1), repeat=4):
                rotated = [winds[moved.index(i)] for i in range(4)]
                actual = detector.apply(detector.state(p, winds), schedule)
                if actual != detector.state(rot[p], rotated):
                    raise ArithmeticError(
                        "Pulse fails to rotate the charge or an oriented face"
                    )
        finals = []
        arm_labels = [[], []]
        for winds in product((0, 1), repeat=4):
            state = detector.state(face[0], winds)
            arms = [
                detector.apply(state, schedules[face[0], x] + schedules[x[face[0]], y])
                for x, y in ((a, b), (b, a))
            ]
            if arms[0] != arms[1]:
                raise ArithmeticError("Rotation orders disagree on the final head word")
            finals.append(arms[0])
            if not any(winds[1:]):
                for labels, state in zip(arm_labels, arms):
                    labels.append(state)
        if len(set(finals)) != 16:
            raise ArithmeticError("Facet rotation erases an input word")
        row.update(
            eligible=True,
            sharp_turn_calibrations=len(schedules) * 16,
            shortest_turn_contacts=min(map(len, schedules.values())),
            longest_turn_contacts=max(map(len, schedules.values())),
            routing_tree=sorted(edge(x, y) for x, y in detector.tree.edges),
            identical_order_words=len(finals),
            enrolled_arm_port_and_head_words=arm_labels,
        )
    return rows


def full_generator(flight, mass, open_rows):
    b = flight.bonds * np.asarray(open_rows)[:, None]
    return np.block([[mass * np.eye(flight.n), b.T], [b, -mass * np.eye(len(b))]])


def propagate(h, waves, time):
    levels, vec = np.linalg.eigh(h)
    return vec @ (np.exp(-1j * time * levels)[:, None] * (vec.T @ waves))


def terminal_read(flight, wave):
    filters = wave[flight.pairs]
    vec = filters.reshape(len(filters), 4)
    choi = vec.T @ vec.conj() / 2
    partial = choi.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    capture = float(np.trace(choi).real)
    flux_effect = wave[flight.flux].conj().T @ wave[flight.flux]
    flux = float(np.trace(flux_effect).real / 2)
    return {
        "no_signed_flux_or_failed_launch": 1 - flux,
        "other_flux_arrival": flux - capture,
        "flux_effect_real": flux_effect.real.tolist(),
        "flux_effect_imag": flux_effect.imag.tolist(),
        "face_captured": capture,
        "reference_partial_transpose_min": float(np.linalg.eigvalsh(partial)[0]),
        "kraus_real": filters.real.tolist(),
        "kraus_imag": filters.imag.tolist(),
    }


def interior_read(flight, wave, ports):
    """Alternative vertex-grade measurement during bare flight. Keep failures unnormalized."""
    mask = sum(bit for e, bit in flight.bits.items() if all(p in ports for p in e))
    accepted = [(i, p, w) for i, (p, w) in enumerate(flight.states) if p in ports]
    outside = sorted({w & ~mask for _, _, w in accepted})
    local = sorted({(p, w & mask) for _, p, w in accepted})
    tensor = np.zeros((len(outside), len(local), 2), complex)
    for i, p, w in accepted:
        tensor[outside.index(w & ~mask), local.index((p, w & mask))] = wave[i]
    matrix = tensor.reshape(len(outside), -1)
    choi = matrix.T @ matrix.conj() / 2
    m = len(local)
    pt = choi.reshape(m, 2, m, 2).transpose(0, 3, 2, 1).reshape(2 * m, 2 * m)
    return {
        "ports": ports,
        "sharp_vertex_records": m,
        "capture": float(np.trace(choi).real),
        "reference_partial_transpose_min": float(np.linalg.eigvalsh(pt)[0]),
    }


def staged(flight, mass, board, cut=False):
    route, n = flight.route, flight.n
    guide_links = along(route)
    guide_rows = np.array([link in guide_links for link in flight.links])
    src = set(route[:3])
    dst = set(route[-3:])
    for link in flight.links:
        if len(link) == 3:
            if route[0] in link:
                src.update(link)
            if route[-1] in link:
                dst.update(link)
    # Include the dead-end spurs attached to these triangular tips, keeping the middle guide closed.
    for link in flight.links:
        if len(link) == 2 and link not in guide_links:
            if src.intersection(link):
                src.update(link)
            if dst.intersection(link):
                dst.update(link)
    encoder = full_generator(flight, mass, [set(link) <= src for link in flight.links])
    decoder = full_generator(flight, mass, [set(link) <= dst for link in flight.links])
    guide = full_generator(flight, mass, guide_rows)
    graph = nx.Graph()
    for row in flight.bonds[guide_rows]:
        i, j = np.flatnonzero(row)
        graph.add_edge(i, j)
    bare_b = circ.tree.guide.path_incidence(board, route)
    reference_h = circ.tree.guide.generator(bare_b, np.array([[mass]]))
    reference = np.linalg.eigvalsh(reference_h)
    project = np.zeros_like(guide)
    dims, errors = [], []
    for component in nx.connected_components(graph):
        verts = sorted(component, key=lambda i: route.index(flight.ports[i]))
        if {flight.ports[i] for i in verts} != set(route) or len(verts) != len(route):
            raise ArithmeticError("Incomplete bare guide component")
        active = np.flatnonzero(
            guide_rows & np.any(flight.bonds[:, verts] != 0, axis=1)
        )
        active = sorted(active, key=lambda idx: guide_links.index(flight.links[idx]))
        np.testing.assert_array_equal(flight.ports[n + np.array(active)], route[:-1])
        ids = np.r_[verts, n + np.array(active)]
        np.testing.assert_array_equal(guide[np.ix_(ids, ids)], reference_h)
        levels, states = np.linalg.eigh(guide[np.ix_(ids, ids)])
        np.testing.assert_allclose(levels, reference, atol=2e-13)
        errors.append(float(np.max(abs(levels - reference))))
        pos = states[:, levels > 0]
        # No occupied site or anchored auxiliary in the complete decoder window at launch.
        forbidden = np.isin(flight.ports[ids], list(dst))
        _, sv, vh = np.linalg.svd(pos[forbidden], full_matrices=True)
        rank = np.count_nonzero(sv > 1e-11)
        selected = pos @ vh[rank:].T
        dims.append(selected.shape[1])
        np.testing.assert_allclose(selected[forbidden], 0, atol=2e-14)
        project[np.ix_(ids, ids)] = selected @ selected.T
    initial = np.vstack((flight.source, np.zeros((len(guide) - n, 2))))
    encoded = propagate(encoder, initial, 4)
    packet = project @ encoded
    np.testing.assert_allclose(packet[np.isin(flight.ports, list(dst))], 0, atol=2e-14)
    if cut:
        cut_link = edge(route[-4], route[-3])
        guide = full_generator(
            flight,
            mass,
            [op and link != cut_link for op, link in zip(guide_rows, flight.links)],
        )
    reads = []
    for t in circ.tree.TIMES:
        in_flight = propagate(guide, packet, t)
        wave = propagate(decoder, in_flight, 4)
        np.testing.assert_allclose(
            wave.conj().T @ wave, packet.conj().T @ packet, atol=5e-13
        )
        reads.append(
            dict(
                time=t,
                **terminal_read(flight, wave),
                interior_apertures=[
                    interior_read(flight, in_flight, region)
                    for region in (route[5:6], route[5:8])
                ],
            )
        )
    return {
        "encoder_wait": 4,
        "decoder_wait": 4,
        "bare_positive_rest_energy": float(min(reference[reference > 0])),
        "guide_spectrum_max_errors": errors,
        "launch_dimensions_per_guide_copy": dims,
        "accepted_launch_by_source": np.sum(abs(packet) ** 2, axis=0).tolist(),
        "reads": reads,
    }


def filters_and_flux(read):
    k = np.array(read["kraus_real"]) + 1j * np.array(read["kraus_imag"])
    flux = np.array(read["flux_effect_real"]) + 1j * np.array(read["flux_effect_imag"])
    return k, flux


def sample(weights, shots, rng):
    np.testing.assert_allclose(np.sum(weights), 1, atol=1e-12)
    if min(weights) < -1e-12:
        raise ArithmeticError("Negative detector weight")
    prob = np.maximum(np.round(weights, 12), 0)
    return rng.multinomial(shots, prob / prob.sum()).tolist()


def coherence_counts(read, rng, panels):
    shots = 2**38
    k, flux = filters_and_flux(read)
    records, matrices = [], []
    settings = (
        np.eye(2),
        np.array([[1, 1], [1, -1]]) / np.sqrt(2),
        np.array([[1, -1j], [1, 1j]]) / np.sqrt(2),
    )
    for source in tomography.INPUTS:
        amplitudes = k @ source
        success = np.linalg.norm(amplitudes) ** 2
        click = np.vdot(source, flux @ source).real
        settings_counts = []
        for basis in settings:
            bins = np.sum(abs(amplitudes @ basis.T) ** 2, axis=0)
            settings_counts.append(
                sample([*bins, click - success, 1 - click], shots, rng)
            )
        records.append(settings_counts)
        diag, x, y = np.array(settings_counts)[:, :2] / shots
        off = ((x[0] - x[1]) - 1j * (y[0] - y[1])) / 2
        matrices.append(np.array([[diag[0], off], [off.conjugate(), diag[1]]]))
    zero, one, xp, xm, yp, ym = matrices
    off = (xp - xm + 1j * (yp - ym)) / 2
    pt = np.block([[zero, off.conj().T], [off, one]]) / 2
    minimum = float(np.linalg.eigvalsh(pt)[0])
    # Union bound for 4 outcomes x 3 settings x 6 input states in each counted panel.
    eps = np.sqrt(np.log(2 * panels * 4 * 3 * 6 / 0.01) / (2 * shots))
    bound = float(4 * np.sqrt(2) * (eps + 5e-12))
    return {
        "shots_per_setting": shots,
        "counts": records,
        "minimum_partial_transpose_eigenvalue": minimum,
        "simultaneous_99_percent_error": bound,
        "certifies_entanglement": bool(minimum + bound < 0),
    }


def order_counts(rows, read, rng):
    k, _ = filters_and_flux(read)
    xyz, edges = circ.tree.icosahedron_from_coordinates()
    xyz = xyz[circ.tree.isomorphism(circ.tree.cell().seams, edges)]
    counts = []
    for row in rows:
        if not row["eligible"]:
            counts.append(None)
            continue
        labels = row["enrolled_arm_port_and_head_words"]
        outputs = sorted({tuple(label) for arm in labels for label in arm})
        arms = []
        for arm in labels:
            embed = np.eye(len(outputs))[
                :, [outputs.index(tuple(label)) for label in arm]
            ]
            arms.append(np.einsum("ij,ejk->eik", embed, k))
        lifts = []
        for turn in row["axes_as_port_permutations"]:
            rot = np.linalg.lstsq(xyz, xyz[list(turn)], rcond=None)[0]
            axis = np.linalg.eigh((rot + rot.T) / 2)[1][:, -1]
            lifts.append(-1j * np.einsum("i,ijk->jk", axis, circ.tree.PAULI))
        a, b = lifts
        np.testing.assert_allclose(a @ b, -b @ a, atol=1e-13)
        spin_up = np.array([1, 0])
        panels = {}
        for label, left, right in (
            ("turn_orders", *arms),
            ("equal_order", arms[0], arms[0]),
            (
                "supplied_spin",
                arms[0][..., None] * (b @ a @ spin_up),
                arms[1][..., None] * (a @ b @ spin_up),
            ),
        ):
            # One half for the mixture of the two inputs, even for the unnormalized herald.
            plus = np.linalg.norm((left + right) / 2) ** 2 / 2
            minus = np.linalg.norm((left - right) / 2) ** 2 / 2
            weights = [
                float(plus),
                float(minus),
                read["other_flux_arrival"],
                read["no_signed_flux_or_failed_launch"],
            ]
            panels[label] = {
                "probabilities": weights,
                "counts": sample(weights, 2**20, rng),
            }
        counts.append(panels)
    return counts


def run():
    board, route = circ.tree.guide.longest_unpinned_route()
    effect, masses = circ.tree.guide.register_effect()
    _, chars = np.linalg.eigh(effect)
    mass = float(masses[0])
    np.testing.assert_allclose(effect @ chars[:, 0], mass * chars[:, 0], atol=1e-14)
    swap = np.array([[0, 1], [1, 0]])
    np.testing.assert_allclose(
        expm(-1j * np.pi * (np.eye(2) - swap) / 2), swap, atol=1e-14
    )
    rng = np.random.default_rng(191182)
    sectors = []
    for sign in (-1, 1):
        base = Flight(board, route, sign, True)
        origin = route[-1] // 12 * 12
        face = [p + origin for p in base.face]
        spurs = sorted(
            (tip, leaf)
            for tip in face[1:]
            for leaf in board.neighbours[tip]
            if leaf not in face and leaf not in route
        )
        spur = spurs[0]
        primary = Flight(board, route, sign, True, spur)
        rows = calibrate(base.face, sign, tuple(p - origin for p in spur))
        joined = primary.read(mass)
        record = staged(primary, mass, board)
        sample_times = (4, 8, 12)
        selected = [r for r in record["reads"] if r["time"] in sample_times]
        variants = []
        for choice in spurs:
            flight = Flight(board, route, sign, True, choice)
            variants.append({"receiver_spur": choice, "continuous": flight.read(mass)})
        sectors.append(
            {
                "sign": sign,
                "source_spur": True,
                "receiver_spur": spur,
                "number_of_vertices_and_incidence_rows": list(primary.bonds.shape)[
                    ::-1
                ],
                "eligible_turns": rows,
                "continuous_spur_census": variants,
                "continuous_symmetric_receiver": base.read(mass),
                "continuous_symmetric_source": Flight(
                    board, route, sign, False, spur
                ).read(mass),
                "continuous_stopped_rings": primary.read(mass, stop_ring=True),
                "continuous_cut": primary.read(mass, cut=True),
                "staged_positive_guide": record,
                "staged_cut_after_preparation": staged(primary, mass, board, cut=True),
                "positive_branch_face_tomography": [
                    {
                        "time": s["time"],
                        **coherence_counts(s, rng, 2 * len(sample_times)),
                    }
                    for s in selected
                ],
                "rotation_shots_per_setting": 2**20,
                "positive_branch_rotation_readings": [
                    {"time": s["time"], "per_axis_pair": order_counts(rows, s, rng)}
                    for s in record["reads"]
                ],
                "continuous_rotation_at_eight": order_counts(
                    rows, joined["reads"][2], rng
                ),
            }
        )
    return {
        "route": route,
        "family_mass_effects": masses.tolist(),
        "seed": 191182,
        "rotation_outcomes": [
            "constructive",
            "destructive",
            "other flux arrival",
            "miss or failed launch",
        ],
        "tomography_source_order": ["0", "1", "X+", "X-", "Y+", "Y-"],
        "tomography_bases": ["Z", "X", "Y"],
        "tomography_outcomes": [
            "plus",
            "minus",
            "other flux arrival",
            "miss or failed launch",
        ],
        "sectors": sectors,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
