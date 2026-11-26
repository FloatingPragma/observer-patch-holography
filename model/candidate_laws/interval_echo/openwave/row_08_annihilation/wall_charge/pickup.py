"""Native signed common-mode and gradient readouts on held material addresses."""

from fractions import Fraction
from itertools import product

import networkx as nx
import numpy as np

from model.executable_model import federation
from model.experiments.material_inertia.pair_spectrum.native_source_probe import (
    circuit,
    physical,
)
from model.experiments.material_inertia.pair_spectrum.native_source_probe.regional_census import (
    meter,
)
from ..wall_flight import extended_geometry

POINTER = meter.pointer
N_CAL = 2**20
QWORDS = tuple(product(range(3), repeat=2))


def read_word(ratio, gradient=False, reverse=False, second=False):
    basis = POINTER.neutral_basis(3, 4)
    sequence = list(basis)
    for target in (0, 1):
        pulse = POINTER.diagonal_word(target, 4, ratio)
        if reverse ^ (gradient and target == 0):
            pulse = POINTER.inverse(pulse)
        sequence += pulse
    if second:
        sequence.append(POINTER.pulse("reference", 4, 4, Fraction(1, 64)))
    return sequence + POINTER.inverse(basis)


def calibration():
    initial = np.zeros((243, 9), complex)
    for j, (a, b) in enumerate(QWORDS):
        initial[np.array([a, b, 1, 2, 0]) @ [81, 27, 9, 3, 1], j] = 1
    labels = 3 * physical.WORDS[:, 3] + physical.WORDS[:, 4]
    hist = np.zeros((2, 2, len(QWORDS), 17))
    rng = np.random.default_rng(831209)
    for gradient, reverse in product((0, 1), repeat=2):
        first = physical.primitive(
            read_word(Fraction(1, 65), gradient, reverse), initial
        )
        branch = first.copy()
        branch[labels != meter.DARK] = 0
        reset = [POINTER.pulse("contact", 3, 4, 1)]
        second = physical.primitive(
            reset + read_word(Fraction(1, 130), gradient, reverse, True), branch
        )
        for j, (a, b) in enumerate(QWORDS):
            counts = np.r_[
                [
                    sum(abs(first[labels == w, j]) ** 2)
                    for w in range(9)
                    if w != meter.DARK
                ],
                [sum(abs(second[labels == w, j]) ** 2) for w in range(9)],
            ]
            q = (b - a if gradient else a + b - 2) * (-1 if reverse else 1)
            p = (1 + np.cos(np.pi * q / 2)) / 2
            y = (1 + np.sin(np.pi * q / 4)) / 2
            predicted = np.zeros(17)
            predicted[5] = p  # first BRIGHT (20), remove first DARK (02)
            predicted[8 + meter.BRIGHT] = (1 - p) * y
            predicted[8 + meter.DARK] = (1 - p) * (1 - y)
            if np.max(abs(counts - predicted)) > 4e-12:
                raise ArithmeticError(
                    "Native ternary instrument differs from signed pickup"
                )
            hist[gradient, reverse, j] = predicted
    reads = np.array(
        [rng.multinomial(N_CAL, np.maximum(0, p)) for p in hist.reshape(-1, 17)]
    )
    return dict(
        probabilities=hist.tolist(),
        counts=reads.reshape(hist.shape).tolist(),
        calibration_shots=N_CAL,
    )


def full_dimer_calibration():
    """One coherent four-port charge reading on two complete flanking dimers."""
    code = np.eye(27)
    labels = circuit.CODE
    basis = circuit.apply(POINTER.neutral_basis(1, 2), code)[
        np.ix_(labels[1], labels[1])
    ]
    rotations = {}
    for ratio in (Fraction(1, 65), Fraction(1, 130)):
        native = circuit.apply(POINTER.diagonal_word(0, 2, ratio), code)
        rotations[ratio] = []
        for idx in labels:
            outside = np.delete(native[:, idx], idx, axis=0)
            if np.max(abs(outside)) > 2e-12:
                raise ArithmeticError(
                    "Pointer phase leaves the material or controller code"
                )
            rotations[ratio].append(native[np.ix_(idx, idx)])
    inputs = [(a, 2 - a, b, 2 - b) for a, b in product((0, 2), repeat=2)]
    # Genuine doubly charged standards use the same four sample ports.
    inputs += [(0, 0, 2, 0), (2, 2, 2, 0)]
    rng = np.random.default_rng(831210)
    out = []
    for loads in inputs:
        arms = []
        for ratio in rotations:
            state = basis[:, 0].copy()
            for value in loads:
                state = rotations[ratio][value] @ state
            if ratio == Fraction(1, 130):
                state = np.diag([1, -1j]) @ state
            arms.append(abs(basis.conj().T @ state) ** 2)
        first, second = arms
        joint = np.array([first[0], first[1] * second[0], first[1] * second[1]])
        charge = sum(loads) - 4
        expected = np.eye(3)[{0: 0, 2: 1, -2: 2}[charge]]
        np.testing.assert_allclose(joint, expected, rtol=0, atol=4e-12)
        counts = rng.multinomial(N_CAL, expected)
        out.append(dict(loads=list(loads), net_charge=charge, counts=counts.tolist()))
    return out


def routing():
    ring, _, wings = extended_geometry()
    board = federation(0)
    footprint = set(ring + wings[0] + wings[1])
    records = []
    for ports in (ring, wings[0][1:], wings[1][1:]):
        owner = ports[0] // 12
        local = set(range(12 * owner, 12 * owner + 12))
        graph = nx.Graph()
        graph.add_nodes_from(sorted(local))
        graph.add_edges_from(
            (a, b) for a, b in board.seams if a in local and b in local
        )
        for cross, park in sorted(
            (c, p)
            for c in local - footprint
            if board.degree[c] == 5
            for p in board.neighbours[c]
            if p in local - footprint
        ):
            for head in sorted(
                c
                for c in board.neighbours[cross]
                if c in local and c != park and board.degree[c] == 5
            ):
                subgraph = graph.subgraph(local - {park, cross})
                if all(nx.has_path(subgraph, p, head) for p in ports):
                    paths = [nx.shortest_path(subgraph, p, head) for p in ports]
                    for path in paths:
                        physical.validate_route(board, path, park, cross)
                    records.append(
                        dict(
                            owner=owner, park=park, cross=cross, head=head, paths=paths
                        )
                    )
                    break
            else:
                continue
            break
        else:
            raise ValueError("No native receiving head")
    return records
