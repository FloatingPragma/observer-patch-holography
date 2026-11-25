"""Compile a two-head comparison from isolated native contacts and blank rail pulses."""

from math import pi

import numpy as np
from scipy.linalg import expm

from .. import measure as screen

IDENTITY, X, Y, Z = screen.station.phase.AXES
PHI = np.array([1, 0, 0, 1]) / np.sqrt(2)
PROJECTOR = np.outer(PHI, PHI)


def contact(odd):
    # Existing block restriction on the actual 99-state marked excess.
    blocks = screen.station.iso.ranging.parity_blocks(1, 8, -1 if odd else 1)
    h = blocks[1:, 1:].transpose(0, 2, 1, 3).reshape(4, 4)
    psi = np.array([0, 1, 1, 0])
    np.testing.assert_array_equal(
        2 * h, (8 + odd) * (2 * np.eye(4) - np.outer(psi, psi))
    )
    return h


def blank():
    phase = screen.station.phase
    emission = phase.formation.emission
    cols = [emission.LAW.SRC.index((2, 1)), emission.LAW.SRC.index((1, 2))]
    h = emission.four_h(8)[np.ix_(cols, cols)] / 4
    np.testing.assert_array_equal(h, 8 * IDENTITY + (IDENTITY - X) / 2)
    swap = expm(-1j * pi * h)
    np.testing.assert_allclose(swap, X, atol=2e-13)
    half = expm(-1j * pi * h / 2)
    raw_axis = half.conj().T @ Z @ half

    def y_read(sign):
        # One signed port outcome at the end of a half blank exchange.
        orientation = np.trace(raw_axis @ Y).real / 2
        return (IDENTITY + sign * np.sign(orientation) * raw_axis) / 2

    phase_time = pi / 48
    reference = np.diag([np.exp(-24j * phase_time), 1])
    z = reference @ reference
    hadamard = reference @ expm(-1.5j * pi * h) @ reference
    ideal = (X + Z) / np.sqrt(2)
    np.testing.assert_allclose(abs(np.trace(hadamard.conj().T @ ideal)), 2, atol=2e-13)
    return z, hadamard, y_read


def extractor(h, odd):
    u = expm(-1j * pi * h / (8 + odd)).reshape(2, 2, 2, 2)
    # Blank input rail 0; the stored reference head ends in 0 with certainty.
    np.testing.assert_allclose(u[:, 1, 0, :], 0, atol=2e-13)
    return u[:, 0, 0, :]


def parity_read(odd):
    """Joint native parity Kraus operators on the unknown coin after erasure."""
    z, h, _ = blank()
    gate = expm(-1j * pi * contact(odd))
    flip = np.kron(z, IDENTITY)
    cycle = flip @ gate @ flip @ gate
    # First |+> rail and X analysis resolve parity. A second |0> rail
    # traverses exactly the same contact sequence, cancelling the coin Z.
    read = (np.kron(h.conj().T, IDENTITY) @ cycle @ np.kron(h, IDENTITY)).reshape(
        2, 2, 2, 2
    )[:, :, 0, :]
    erase = cycle.reshape(2, 2, 2, 2)[:, :, 0, :]
    np.testing.assert_allclose(erase[1], 0, atol=2e-13)
    result = erase[0] @ read
    np.testing.assert_allclose(result[1 - odd], 0, atol=2e-13)
    np.testing.assert_allclose(abs(np.trace(result[odd])), 2, atol=2e-13)
    return result


def bell_effect(reference_odd, target_odd, *, pulse=True):
    z, hadamard, y_read = blank()
    reference = extractor(contact(reference_odd), reference_odd)
    target_h = contact(target_odd)
    target = extractor(target_h, target_odd)
    half = expm(-0.5j * pi * target_h / (8 + target_odd)) if pulse else np.eye(4)
    flip = np.kron(z, IDENTITY)
    # Native z cancels the known extraction phase; Hadamard acts on the rail.
    start = np.kron(hadamard @ z @ reference, IDENTITY)
    decoder = flip @ half @ flip @ half @ start
    # Read and remove the output rail first. A fresh excess transfers the
    # target head to the blank read dimer. That extraction reverses rail Y.
    plus_y = y_read(1)
    target_plus_y = target.conj().T @ y_read(-1) @ target
    success = np.kron(plus_y, target_plus_y)
    return decoder.conj().T @ success @ decoder


def geometry():
    board = screen.camera.front.optics.pairing.formation.federation(1)
    route = screen.camera.front.optics.pairing.formation.region(board)
    processing = screen.station.phase.read_dimers(board, route)
    reserved = set(route) | set(processing[0]) | set(processing[1])
    references = []
    for x, y in processing:
        choices = [
            (u, v)
            for u in board.neighbours[x]
            for v in board.neighbours[y]
            if u not in reserved
            and v not in reserved
            and v in board.neighbours[u]
            and board.degree[u] == board.degree[v] == 5
            and u // 12 == v // 12 == x // 12
        ]
        references.append(min(choices))
    target = np.asarray(route).reshape(2, 2)
    triples = [
        (list(map(int, ref)), list(map(int, proc)), list(map(int, tar)))
        for ref, proc, tar in zip(references, processing, target, strict=True)
    ]
    assert len({s for mount in triples for pair in mount for s in pair}) == 12
    for mount in triples:
        for pair in mount:
            assert pair[1] in board.neighbours[pair[0]]
        for a, b in zip(mount, mount[1:]):
            assert all(y in board.neighbours[x] for x, y in zip(a, b, strict=True))
    return triples


def shuttle(mount):
    """Move a logical rail from its original to its auxiliary mark in canonical port order."""
    amp = np.eye(2, dtype=complex)
    original, processing, auxiliary = reversed(mount)
    for ports, next_ports in ((original, processing), (processing, auxiliary)):
        addresses = sorted(ports) + sorted(next_ports)
        h = np.zeros((4, 4), complex)
        for a, b in zip(ports, next_ports, strict=True):
            i, j = addresses.index(a), addresses.index(b)
            h[np.ix_([i, j], [i, j])] = 8 * IDENTITY + (IDENTITY - X) / 2
        pulse = expm(-1j * pi * h)
        np.testing.assert_allclose(pulse[:2, :2], 0, atol=2e-13)
        transfer = pulse[2:, :2]
        # The physical wires to the auxiliary mark cross in sorted-port
        # convention. A single blank processing-dimer swap removes this X.
        if np.argmax(abs(transfer[:, 0])) == 1:
            extra = expm(-1j * pi * (8 * IDENTITY + (IDENTITY - X) / 2))
            transfer = transfer @ extra if ports == processing else extra @ transfer
        amp = transfer @ amp
    np.testing.assert_allclose(amp, IDENTITY, atol=3e-13, rtol=0.0)
    return amp


def paired_effect(a, b):
    # Reference-left, reference-right, flight-left, flight-right in trace products.
    return (
        np.kron(a, b)
        .reshape((2,) * 8)
        .transpose(0, 2, 1, 3, 4, 6, 5, 7)
        .reshape(16, 16)
    )


def report():
    tag_errors = [float(np.linalg.norm(parity_read(odd)[1 - odd], 2)) for odd in (0, 1)]
    errors = [
        [float(np.linalg.norm(bell_effect(r, t) - PROJECTOR, 2)) for t in (0, 1)]
        for r in (0, 1)
    ]
    assert np.max(errors) < 3e-13
    no = bell_effect(0, 0, pulse=False)
    return dict(
        mounts=geometry(),
        corrected_shuttle_errors=[
            float(np.linalg.norm(shuttle(m) - IDENTITY, 2)) for m in geometry()
        ],
        parity_wrong_port_operator_norm=tag_errors,
        parity_tag_and_erase_dwell="four pi exposures with two auxiliary rails",
        local_bell_error_by_reversal=errors,
        parked_reference_gap=24,
        blank_swap_dwell="pi",
        rail_quarter_phase_dwell="pi/48",
        extraction_dwell_by_reversal=["pi/8", "pi/9"],
        entangling_dwell_by_reversal=["pi/16", "pi/18"],
        omitted_pulses_phi_acceptance=float(np.trace(no @ PROJECTOR).real),
    )
