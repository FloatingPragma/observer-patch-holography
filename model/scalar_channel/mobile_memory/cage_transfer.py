"""A record-controlled exchange of whole scalar cages across an external seam."""

from functools import lru_cache
from math import isfinite, pi, sin
from numbers import Real

import numpy as np
from scipy.linalg import block_diag, expm

from ..carrier import require
from ..cage import construction as cage
from ..cage.instruments import fock_operators, lift, pointers, rotation


TOLERANCE = 2e-11
PORTS_PER_CELL = 12


@lru_cache(maxsize=3)
def _carrier(level):
    return cage.federation(level)


def _ring(carrier, cell):
    ports = range(PORTS_PER_CELL * cell, PORTS_PER_CELL * (cell + 1))
    eligible = [
        ring
        for ring in cage.rings(carrier, ports)
        if all(port in ports and carrier["degree"][port] == 5 for port in ring)
    ]
    require(bool(eligible), "cell has no exact degree-five cage")
    ring = eligible[0]
    profile = cage.ring_profile(carrier, ring)
    incidence = carrier["C"]
    wedge_support = incidence.T @ profile
    require(
        not np.any(carrier["Q"] @ profile)
        and np.array_equal(incidence @ wedge_support, 10 * profile),
        "the ring lost its exact scalar pair",
    )
    require(
        all(all(port in ports for port in carrier["wedges"][i])
            for i in np.flatnonzero(wedge_support)),
        "the scalar pair reaches outside its cell",
    )
    _, _, block, residual = cage.pair(carrier, profile)
    require(residual < TOLERANCE, "the carrier and scalar pair do not intertwine")
    return ring, block


def _record(state, projector):
    probability = float(np.linalg.norm(projector @ state) ** 2)
    require(
        abs(float(np.linalg.norm(state)) - 1.0) < TOLERANCE
        and 0.0 <= probability <= 1.0 + TOLERANCE,
        "the detector has an invalid record law",
    )
    return pointers(min(1.0, probability)).tolist()


def transfer(level=0, seam=None, memory=True, coupling=1.0, wait=None):
    """Return all fuel-pointer records for a sharp cage sent toward an adjacent cell.

    The supplied link is one positive quartic term on the two ring patches.
    A classical memory bit enables the term after the source-on sender pulse.
    The term remains active through the wait and the finite receiver pulse.
    The reference fills the lower mode of both cage pairs. Spectator modes in
    the negative sea factor because the two pairs reduce the scalar law.
    """
    require(type(level) is int and level in (0, 1, 2), "unknown federation level")
    require(type(memory) is bool, "current memory must be a bit")
    require(
        isinstance(coupling, Real) and isfinite(coupling) and coupling >= 0,
        "coupling must be nonnegative and finite",
    )
    if wait is None:
        wait = pi / (2 * coupling) if coupling else pi / 2
    require(isinstance(wait, Real) and isfinite(wait) and wait >= 0, "invalid wait")
    carrier = _carrier(level)
    if seam is None:
        seam = next(
            (a, b)
            for a, b in carrier["edges"]
            if a // PORTS_PER_CELL != b // PORTS_PER_CELL
        )
    require(
        isinstance(seam, (list, tuple)) and len(seam) == 2
        and all(type(port) is int for port in seam),
        "an oriented seam must contain two port indices",
    )
    sender, receiver = seam
    sender_cell, receiver_cell = sender // PORTS_PER_CELL, receiver // PORTS_PER_CELL
    require(
        sender_cell != receiver_cell and tuple(sorted(seam)) in carrier["index"],
        "the selected seam does not connect two cells",
    )
    sender_ring, sender_block = _ring(carrier, sender_cell)
    receiver_ring, receiver_block = _ring(carrier, receiver_cell)
    require(
        np.allclose(sender_block, receiver_block, atol=TOLERANCE, rtol=0),
        "the two cage pairs have different free generators",
    )

    # The first four modes are the two exact scalar pairs. The last mode is
    # the receiver's blank fuel. Fermion signs are kept on the full 32 states.
    ops = fock_operators(5)
    free = lift(block_diag(sender_block, receiver_block, 0.0), ops)
    lower, upper = rotation(sender_block[0, 0], sender_block[0, 1])
    lower_ops = [lower[0] * ops[i] + lower[1] * ops[i + 1] for i in (0, 2)]
    upper_ops = [upper[0] * ops[i] + upper[1] * ops[i + 1] for i in (0, 2)]
    destroy = [lo.conj().T @ hi for lo, hi in zip(lower_ops, upper_ops)]
    difference = destroy[0] - destroy[1]
    link = difference.conj().T @ difference

    vacuum = np.zeros(len(free), dtype=complex)
    vacuum[0] = 1.0
    reference = lower_ops[1].conj().T @ lower_ops[0].conj().T @ vacuum
    frame = np.column_stack([exciton.conj().T @ reference for exciton in destroy])
    exact_block = np.array([[1.0, -1.0], [-1.0, 1.0]])
    require(
        np.linalg.norm(frame.conj().T @ frame - np.eye(2)) < TOLERANCE
        and np.linalg.norm(link @ reference) < TOLERANCE
        and np.linalg.norm(link @ frame - frame @ exact_block) < TOLERANCE
        and np.linalg.norm(free @ link - link @ free) < TOLERANCE,
        "the quartic link changes the scalar reference or leaves the cage sector",
    )

    # The declared sharp cage preparation is a finite potential on the
    # sender's seam profile, with the free scalar source running throughout.
    sender_number = ops[0].conj().T @ ops[0]
    preparation = expm(
        -1j * cage.T_SHARP * (free - (2 * cage.ENERGY ** 2 / cage.REST) * sender_number)
    )
    prepared = preparation @ reference
    require(
        abs(abs(np.vdot(frame[:, 0], prepared)) - 1.0) < TOLERANCE,
        "the source-on sender pulse did not prepare the sharp cage",
    )
    active = float(coupling) if memory else 0.0
    source = free + active * link
    evolve = expm(-1j * wait * source)
    received = evolve @ prepared
    reference_waited = evolve @ reference
    receiver_exciton = destroy[1].conj().T @ destroy[1]
    exciton_record = _record(received, receiver_exciton)
    reference_exciton_record = _record(reference_waited, receiver_exciton)
    require(
        np.linalg.norm(received - frame @ (frame.conj().T @ received)) < TOLERANCE
        and abs(exciton_record[3] - sin(active * wait) ** 2) < TOLERANCE
        and reference_exciton_record[3] < TOLERANCE,
        "the full Fock evolution differs from the bound-cage record law",
    )

    receiver_hop = ops[2].conj().T @ ops[4] + ops[4].conj().T @ ops[2]
    capture = expm(-1j * pi / (2 * cage.STRENGTH) * (source + cage.STRENGTH * receiver_hop))
    fuel_number = ops[4].conj().T @ ops[4]
    return {
        "level": level,
        "seam": list(seam),
        "sender_ring": list(sender_ring),
        "receiver_ring": list(receiver_ring),
        "memory": memory,
        "coupling": float(coupling),
        "wait": float(wait),
        "exciton_record": exciton_record,
        "reference_exciton_record": reference_exciton_record,
        "receiver_record": _record(capture @ received, fuel_number),
        "reference_record": _record(capture @ reference_waited, fuel_number),
    }
