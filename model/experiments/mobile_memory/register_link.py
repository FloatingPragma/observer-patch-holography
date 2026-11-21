"""Compile an addressed cage link into the scalar register algebra."""

from collections import deque
from itertools import combinations
from math import atan2, cos, hypot, pi, sin

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.scalar_channel import compiler
from model.experiments.cage import construction as cage, instruments, registers
from model.candidate_laws.scalar_channel.carrier import require
from .cage_transfer import _ring


TOLERANCE = 2e-11


def _chart(level, seam):
    full = cage.federation(level)
    a, b = seam
    require(a // 12 != b // 12 and tuple(sorted(seam)) in full["index"],
            "the register chart needs an external carrier seam")
    ports = [12 * cell + port for cell in (a // 12, b // 12) for port in range(12)]
    label = {port: index for index, port in enumerate(ports)}
    patch_edges = [
        (label[i], label[j]) for i, j in full["edges"]
        if i in label and j in label and (i // 12 == j // 12 or {i, j} == {a, b})
    ]
    patch = cage.incidences(patch_edges, 24, sparse=True)
    edges, wedges = patch["edges"], patch["wedges"]
    nsource = len(edges) + len(wedges)
    fuels = list(range(nsource, nsource + 24))
    references = list(range(nsource + 24, nsource + 48))
    owners = {}
    cells = [set(), set()]
    for mode, (i, j) in enumerate(edges):
        owners[(mode, fuels[i])] = i
        owners[(mode, fuels[j])] = j
        if i // 12 == j // 12:
            cells[i // 12].add(mode)
    for index, (i, j, k) in enumerate(wedges, start=len(edges)):
        for pair in ((i, j), (j, k)):
            owners[(patch["index"][tuple(sorted(pair))], index)] = j
        owners[(index, fuels[j])] = j
        if i // 12 == k // 12 == j // 12:
            cells[j // 12].add(index)
    for port, (fuel, reference) in enumerate(zip(fuels, references, strict=True)):
        owners[(fuel, reference)] = port
        cells[port // 12].update((fuel, reference))
    graph = compiler.encode_graph(nsource + 48, list(owners), owners)

    full_wedges = {wedge: len(full["edges"]) + i for i, wedge in enumerate(full["wedges"])}
    profiles = []
    for port in seam:
        ring, _ = _ring(full, port // 12)
        source, _, error = cage.modes(full, cage.ring_profile(full, ring))
        pair, block, patch_error = cage.modes(
            patch, cage.ring_profile(patch, tuple(label[v] for v in ring)),
        )
        image = np.zeros_like(source)
        for index in np.flatnonzero(np.any(pair, axis=1)):
            if index < len(edges):
                i, j = edges[index]
                mapped = full["index"][tuple(sorted((ports[i], ports[j])))]
            else:
                i, j, k = wedges[index - len(edges)]
                mapped = full_wedges[ports[i], ports[j], ports[k]]
            image[mapped] = pair[index]
        require(max(error, patch_error, np.max(abs(image - source))) < TOLERANCE,
                "the patch pair differs from the actual federation pair")
        frame = np.zeros((graph["n"], 2))
        frame[:nsource] = pair @ instruments.rotation(block[0, 0], block[0, 1]).T
        profiles.append(frame)
    centers = [label[a], label[b]]
    bridge = patch["index"][tuple(sorted(centers))]
    # The external seam mode is a spectator. Include one more spectator
    # in the literal check so both parities of the bridge are represented.
    anchor = next(mode for mode, ends in enumerate(edges) if centers[0] in ends and mode != bridge)
    core = [references[centers[0]], fuels[centers[0]], bridge,
            fuels[centers[1]], references[centers[1]], anchor]
    return graph, cells, profiles, core, ports


def _compress(graph, frame, available, target, column):
    """Rotate leaves into their parents along links of the addressed cell."""
    parents = {target: target}
    pending = deque([target])
    while pending:
        parent = pending.popleft()
        for edge in graph["inc"][parent]:
            i, j = graph["links"][edge]
            child = j if i == parent else i
            if child in available and child not in parents:
                parents[child] = parent
                pending.append(child)
    require(set(parents) == available, "the addressed cell route disconnected")
    gates = []
    for child in reversed(list(parents)[1:]):
        parent = parents[child]
        u, v = frame[[parent, child], column]
        if abs(v) <= 1e-15:
            continue
        angle = atan2(v, u)
        c, s = cos(angle), sin(angle)
        frame[[parent, child]] = np.array([[c, s], [-s, c]]) @ frame[[parent, child]]
        i, j = sorted((parent, child))
        edge = graph["ix"][i, j]
        orientation = 1 if parent == i else -1
        words = compiler.hopping_words(graph, edge, 1j * orientation)
        require(not any(compiler.anti(p, q) for (_, p), (_, q) in combinations(words, 2)),
                "the two Pauli rotations of a native hop do not commute")
        gates.append({
            "modes": [parent, child],
            "angle": angle,
            "words": [[coefficient, list(word)] for coefficient, word in words],
        })
        require(abs(frame[child, column]) <= TOLERANCE and
                abs(frame[parent, column] - hypot(u, v)) <= TOLERANCE,
                "a native rotation failed to eliminate the routed coordinate")
    return gates


def _lowering_words(graph, lower, upper):
    """Integer Pauli polynomial for four times f_lower† f_upper."""
    edge = graph["ix"][tuple(sorted((lower, upper)))]
    result = {}
    for value, word in compiler.hopping_words(graph, edge, 1):
        registers.add(result, (int(2 * value), 0), word)
    sign = 1 if lower < upper else -1
    for value, word in compiler.hopping_words(graph, edge, 1j * sign):
        registers.add(result, (0, -int(2 * value)), word)
    return result


def _adjoint(polynomial):
    adjoint = {}
    for (x, z), (real, imag) in polynomial.items():
        registers.add(adjoint, (real, -imag), (2 * (x & z).bit_count(), x, z))
    return adjoint


def _quartic_batches(graph, core):
    """Separate the commuting onsite and exchange batches of 16 times the link."""
    left, right = _lowering_words(graph, *core[:2]), _lowering_words(graph, core[4], core[3])
    left_dagger, right_dagger = _adjoint(left), _adjoint(right)
    onsite = registers.combination(
        (1, registers.product(left_dagger, left)),
        (1, registers.product(right_dagger, right)),
    )
    exchange = registers.combination(
        (-1, registers.product(left_dagger, right)),
        (-1, registers.product(right_dagger, left)),
    )
    for batch in (onsite, exchange):
        words = [(0, x, z) for x, z in batch]
        require(registers.hermitian(batch) and not any(
            compiler.anti(p, q) for p, q in combinations(words, 2)
        ), "the cage-link batch is not a commuting Hermitian polynomial")
    require(not registers.combination(
        (1, registers.product(onsite, exchange)), (-1, registers.product(exchange, onsite))
    ), "the onsite and exchange sums do not commute")
    return onsite, exchange


def _dense_pauli(word, qubits):
    phase, x, z = word
    matrix = np.zeros((1 << qubits, 1 << qubits), complex)
    for column in range(1 << qubits):
        matrix[column ^ x, column] = 1j ** phase * (-1) ** (column & z).bit_count()
    return matrix


def _literal_core_check(wait):
    """Compare the five-mode bridge plus its parity anchor against literal CAR."""
    core = compiler.encode_graph(6, [(0, 1), (1, 2), (1, 5), (2, 3), (3, 4)])
    batches = _quartic_batches(core, range(6))
    dimension = 1 << core["M"]
    code_basis = np.eye(dimension, dtype=complex)
    fock_basis = np.eye(64, dtype=complex)
    lifts = [_dense_pauli(word, core["M"]) for word in compiler.logical_pauli_generators(core)]
    destroy = instruments.fock_operators(6)
    gamma = [a + a.conj().T for a in destroy]
    number_signs = [np.eye(64) - 2 * a.conj().T @ a for a in destroy]
    flips = []
    for i in range(5):
        flip = -1j * gamma[i] @ gamma[5]
        for z in number_signs[i + 1:5]:
            flip = flip @ z
        flips.append(flip)
    columns = []
    for state in range(dimension):
        code, fock = code_basis[:, 0], fock_basis[:, 0]
        for k in compiler.bits(state):
            code, fock = lifts[k] @ code, flips[k] @ fock
        columns.append((code, fock))
    v, w = (np.column_stack(part) for part in zip(*columns, strict=True))
    d = destroy[0].conj().T @ destroy[1] - destroy[4].conj().T @ destroy[3]
    literal = w.conj().T @ (d.conj().T @ d) @ w
    encoded = np.zeros((dimension, dimension), complex)
    pulse = np.eye(dimension, dtype=complex)
    for batch in batches:
        for (x, z), coefficient in sorted(batch.items()):
            term = complex(*coefficient) / 16 * _dense_pauli((0, x, z), core["M"])
            encoded += term
            pulse = expm(-1j * wait * term) @ pulse
    residual = float(np.max(abs(v.conj().T @ encoded @ v - literal)))
    occupation_residual = max(
        float(np.max(abs(v.conj().T @ _dense_pauli(word, core["M"]) @ v
                         - w.conj().T @ sign @ w)))
        for word, sign in zip(core["B"], number_signs, strict=True)
    )
    require(max(residual, occupation_residual) < TOLERANCE,
            "the encoded core differs from literal fermions")

    exact = expm(-1j * wait * literal)
    unitary_residual = float(np.linalg.norm(v.conj().T @ pulse @ v - exact, ord=2))
    require(unitary_residual < TOLERANCE, "the commuting batches differ from the literal pulse")
    return {
        "literal_even_basis_columns": dimension,
        "literal_generator_residual": residual,
        "literal_occupation_residual": occupation_residual,
        "literal_unitary_residual": unitary_residual,
    }


def compile_link(level=0, seam=(10, 20), wait=pi / 2):
    """Return the two-cell register rotations and exact controlled pulse.

    The chart retains the supplied two-cell seam, fuel and reference modes.
    The rotation word, its inverse and the central pulse are addressed gates.
    Source terms present during their execution require a separate schedule.
    """
    require(np.isfinite(wait) and wait >= 0, "the pulse needs a finite nonnegative wait")
    graph, cells, frames, core, ports = _chart(level, seam)
    routes = []
    residuals = []
    for cell, frame, pair in zip(cells, frames, ((core[0], core[1]), (core[4], core[3])), strict=True):
        lower, upper = pair
        forward = _compress(graph, frame, cell, lower, 0)
        forward.extend(_compress(graph, frame, cell - {lower}, upper, 1))
        expected = np.zeros_like(frame)
        expected[lower, 0] = expected[upper, 1] = 1
        residuals.append(float(np.max(abs(frame - expected))))
        routes.append(forward)
    require(max(residuals) < TOLERANCE, "the native routes failed to compress the actual cage modes")
    batches = _quartic_batches(graph, core)
    loops = [word for _, _, word in compiler.fundamental_cycles(graph)[0]]
    quartic_words = [(0, x, z) for batch in batches for x, z in batch]
    route_words = [tuple(word) for route in routes for gate in route for _, word in gate["words"]]
    require(not any(compiler.anti(p, q) for p in quartic_words + route_words for q in loops),
            "the routed link leaves its BKSF cycle code")
    central_support = sorted({
        ports[owner] for word in quartic_words for owner in compiler.support_stars(graph, word)
    })
    require(central_support == sorted(seam), "the central link uses another observer's registers")
    return {
        "oriented_seam": list(seam),
        "patch_modes": graph["n"],
        "patch_registers": graph["M"],
        "cycle_checks": len(loops),
        "route_residual": max(residuals),
        "routing_gate_counts": [len(route) for route in routes],
        "central_support_ports": central_support,
        "routes": routes,
        "register_owner_ports": [ports[owner] for owner in graph["owner"]],
        "core_modes": core,
        "central_word_counts": [len(batch) for batch in batches],
        "quartic_batches_times_16": {
            name: [[x, z, *coefficient] for (x, z), coefficient in sorted(batch.items())]
            for name, batch in zip(("onsite", "exchange"), batches, strict=True)
        },
        "central_wait": float(wait),
        **_literal_core_check(wait),
    }
