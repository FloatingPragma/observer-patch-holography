"""Time a chosen-bit reply across observer zero's settled material boundary."""

import argparse
from collections import deque
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from observer_zero_radar import ArchiveRadar, MANIFEST_SHA
from model.evaluation.patch_spatial_readout import infer, transitions
from model.executable_model import cell
from model.executable_model.symmetry import rotations
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_28_acceleration.cut_signal import (
    deposited_collar,
    receive,
)

LAW = candidate.native_law
PI = math.pi


def decode_twice(collar, bit, direction, half, contacts):
    """Keep each pointer word and its unnormalised conditional seam amplitudes."""
    branches = [(0, collar)]
    for interaction in contacts:
        following = []
        for word, memory in branches:
            wave = receive(memory, bit, direction, interaction, half)
            for pointer in (0, 1):
                state = wave[:, bit, pointer, 1 - pointer].copy()
                following.append((2 * word + pointer, state))
        branches = following
    return branches


@lru_cache(maxsize=None)
def pulse_operators(kappa):
    exchange, contrast = LAW.local_term_parts()
    half = expm(-0.5j * PI * LAW.seam_term(candidate.KAPPA))
    contacts = [
        expm(-1j * PI * (exchange + kappa * contrast) / den) for den in (32, 33)
    ]
    return half, contacts


@lru_cache(maxsize=None)
def channel(sign, direction, coin=None, kappa=candidate.KAPPA):
    half, contacts = pulse_operators(kappa)
    opposite = "right" if direction == "left" else "left"
    output = []
    for challenge in (0, 1):
        first = decode_twice(
            deposited_collar(sign, coin), challenge, direction, half, contacts
        )
        for choice in (0, 1):
            joint = np.zeros((2, 2))
            for word, memory in first:
                decoded = int(word != 0)
                reply = decoded ^ choice
                for returned, state in decode_twice(
                    memory, reply, opposite, half, contacts
                ):
                    joint[decoded, int(returned != 0)] += np.vdot(state, state).real
            np.testing.assert_allclose(joint.sum(), 1, atol=5e-13, rtol=0)
            output.append(
                dict(challenge=challenge, remote_choice=choice, joint=joint.tolist())
            )
    return output


def path_tree(source):
    parent = {0: None}
    queue = deque([0])
    while queue:
        here = queue.popleft()
        for other, tick in sorted(source.links(here)):
            if not tick and other not in parent:
                parent[other] = here
                queue.append(other)
    return parent


def clock_read(reads, opposite, shots):
    if not reads or reads[0] != 0 or any(a == b for a, b in zip(reads, reads[1:])):
        raise ValueError("Invalid flip register")
    return len(reads) - 1 + 2 * math.asin(math.sqrt(opposite / shots)) / PI


def measure_clock(intervals, rng, shots=2**18):
    """Independent initially occupied dimer, inspected at complete flips and at return."""
    h = LAW.seam_term(candidate.KAPPA)
    a, b = (LAW.SRC.index(s) for s in ((1, 0), (0, 1)))
    full = expm(-1j * PI * h)
    # The full-window read is deterministic at each calibrated flip.
    reads = [0]
    for _ in range(math.floor(intervals)):
        incoming = (a, b)[reads[-1]]
        p_one = abs(full[b, incoming]) ** 2
        if min(abs(p_one), abs(p_one - 1)) > 1e-12:
            raise ValueError("Reference flip is ambiguous")
        reads.append(int(p_one > 0.5))
    fraction = intervals - (len(reads) - 1)
    terminal = expm(-1j * PI * fraction * h)
    probability = abs(terminal[(a, b)[1 - reads[-1]], (a, b)[reads[-1]]]) ** 2
    opposite = int(rng.binomial(shots, probability))
    return dict(
        complete_flip_reads=reads,
        shots=shots,
        opposite=opposite,
        measured_flip_intervals=clock_read(reads, opposite, shots),
    )


def local_gram():
    h = candidate.generator(cell()).toarray()
    _, _, projector, rank, gram = infer(transitions(h, 1 / 8))
    if rank != 3:
        raise ValueError("Incorrect spatial calibration rank")
    np.testing.assert_allclose(gram, 4 * projector, atol=1e-12)
    for time in (2 / 8, 3 / 8):
        _, _, _, held_rank, held = infer(transitions(h, time))
        if held_rank != rank:
            raise ValueError("Held rank does not transfer")
        np.testing.assert_allclose(held, gram, atol=1e-12)
    return gram


def load_word(row, epoch):
    initial = [(row[0] >> port) & 1 for port in range(12)]
    written = [int(0 < tick <= epoch) for tick in row[1:13]]
    return initial + written


def displacement(word, gram):
    delta = (1 - 2 * np.array(word[:12])) * np.array(word[12:])
    return dict(
        signed_changes=delta.tolist(),
        gram_length=math.sqrt(max(0, float(delta @ gram @ delta))),
    )


@lru_cache(maxsize=None)
def pointer_operators(bit, direction):
    # Linear maps of the four pointer words on the same eleven-state collar.
    half, contacts = pulse_operators(candidate.KAPPA)
    columns = [decode_twice(mem, bit, direction, half, contacts) for mem in np.eye(11)]
    return [np.column_stack([output[k][1] for output in columns]) for k in range(4)]


def send_word(word, sign, direction):
    mem = deposited_collar(sign, None)
    rho = np.outer(mem, mem.conj())
    decoded = []
    for bit in word:
        ops = pointer_operators(bit, direction)
        branches = [op @ rho @ op.conj().T for op in ops]
        weights = [branches[0].trace().real, sum(m.trace().real for m in branches[1:])]
        output = int(weights[1] > weights[0])
        decoded.append(output)
        # No renormalisation: the final trace includes all earlier decoding errors
        # as excluded branches. The original collar is prepared only once.
        rho = sum(branches[1:]) if output else branches[0]
    return decoded, float(rho.trace().real)


def midpoint_accelerations(bank):
    result = []
    stamps = [bank["origin_emission_ticks"].index(s) for s in (0, 4, 8)]
    for i, target in enumerate(bank["reflectors"]):
        waits = [bank["echo_waits"][index][i] for index in stamps]
        if any(w is None for w in waits):
            raise ValueError("Expected every terminal reflector to answer")
        r = [Fraction(w, 2) for w in waits]
        t = [s + v for s, v in zip((0, 4, 8), r)]
        velocities = [(r[j + 1] - r[j]) / (t[j + 1] - t[j]) for j in (0, 1)]
        acceleration = 2 * (velocities[1] - velocities[0]) / (t[2] - t[0])
        result.append(
            dict(
                target=target,
                return_waits=waits,
                middle_radar_time=float(t[1]),
                acceleration=str(acceleration),
                sign=(acceleration > 0) - (acceleration < 0),
            )
        )
    return result


def frame_choices(gram, near, far):
    """Relative proper frame maps if the one facing signal bears opposite port directions."""
    opposite = int(np.argmin(gram[near]))
    candidates = [r for r in rotations() if r[far] == opposite]
    second = next(p for p in range(12) if abs(gram[p, far]) < 0.99)
    return dict(
        maps_given_one_direction=len(candidates),
        second_far_port=second,
        possible_second_near_directions=sorted({r[second] for r in candidates}),
    )


def measure(root):
    source = ArchiveRadar(root)
    parent = path_tree(source)
    gram = local_gram()
    cut = sorted(
        (p, q, t) for p in parent for q, t in source.links(p) if q not in parent
    )
    blank = LAW.seam_term(candidate.KAPPA)
    half = expm(-1j * PI * blank / 2)
    np.testing.assert_allclose(expm(-1.5j * PI * blank), half.conj().T, atol=1e-13)
    swap = expm(-1j * PI * blank)
    for bit in (0, 1):
        np.testing.assert_allclose(
            abs(swap[LAW.SRC.index((0, bit)), LAW.SRC.index((bit, 0))]) ** 2,
            1,
            atol=1e-13,
        )
    rng = np.random.default_rng(131004)
    records, kernels = [], {}
    for inside, outside, tick in cut:
        if (inside, tick) not in source.links(outside) or tick <= 0:
            raise ValueError("The boundary reports disagree")
        left, right = sorted((inside, outside))

        def initial(port):
            owner, p = divmod(port, 12)
            return 2 * ((source.row(owner)[0] >> p) & 1)

        sign = initial(left) - initial(right)
        if abs(sign) != 2:
            raise ValueError("Expected a sharp signed deposit")
        # A spare blank internal contact at each endpoint can host the local rail.
        rails = []
        for p in (inside, outside):
            rails.append(
                next(
                    q
                    for q, t in sorted(source.links(p))
                    if not t and q // 12 == p // 12
                )
            )
        direction = "left" if inside == left else "right"
        key = f"{sign:+d}/{direction}"
        if key not in kernels:
            kernels[key] = dict(
                coherent=channel(sign, direction),
                read_coin_zero=channel(sign, direction, coin=0),
                read_coin_one=channel(sign, direction, coin=1),
                removed_contrast=channel(sign, direction, kappa=0),
            )
        route = [inside]
        while parent[route[-1]] is not None:
            route.append(parent[route[-1]])
        route.reverse()
        # Blank transfers and programmed pointer gates run serially. A pointer
        # uses a half swap, a short interaction and a forward-time inverse.
        hops = len(route) - 1
        duration = 2 * hops + 4 * (0.5 + 1.5) + 2 * (1 / 32 + 1 / 33)
        clock_pair = next(
            (p, q)
            for p in range(12)
            if p not in route and p not in rails
            for q, t in sorted(source.links(p))
            if not t and p < q < 12 and q not in route and q not in rails
        )
        remote_word = load_word(
            source.row(outside // 12), source.manifest["settlement_round"]
        )
        received, success = send_word(
            remote_word, sign, "right" if direction == "left" else "left"
        )
        if received != remote_word or abs(success - 1) > 1e-11:
            raise ValueError("Remote load record failed to transfer")
        records.append(
            dict(
                inside=divmod(inside, 12),
                outside=divmod(outside, 12),
                write_round=tick,
                kernel=key,
                origin_clock_ports=clock_pair,
                facing_port_frame_test=frame_choices(gram, inside % 12, outside % 12),
                received_load_word=received,
                probability_entire_word_correct=success,
                remote_local_displacement=displacement(received, gram),
                rail_ports=[divmod(p, 12) for p in rails],
                outward_blank_route=[divmod(p, 12) for p in route],
                immediate=measure_clock(duration, rng),
                delayed=measure_clock(duration + 0.5, rng),
            )
        )
    bank = json.loads(
        Path(__file__).with_name("observer_zero_capacity.json").read_text()
    )
    if bank["manifest_sha256"] != MANIFEST_SHA or bank["origin"] != [0, 0]:
        raise ValueError("Radar receipt belongs to another history or station")
    own = source.row(0)
    return dict(
        manifest_sha256=MANIFEST_SHA,
        origin=[0, 0],
        local_position_readout=dict(
            rank=3,
            units="Unit diagonal of the slow-band Gram.",
            records=[
                dict(archive_round=t, **displacement(load_word(own, t), gram))
                for t in sorted(set(own[1:13]))
            ],
        ),
        midpoint_radar_accelerations=midpoint_accelerations(bank),
        settled_material_ports=len(parent),
        clock_unit="One complete load transfer on an isolated blank reference dimer (pi inverse exchange-gap units).",
        apparatus="Matched settled preparations with archived signs, native unread coherent coin, earlier mobile loads parked, clean addressed tags, local private rails and pointer records. Selected gates are serviced sequentially. Other terms are held. Reads, parking, switching and classical processing carry zero assigned duration. The delayed arm holds for half a reference flip at the remote end.",
        signal="Two sequential pointer uses on each crossing. Remote returns its decoded challenge XOR an independent local choice. All four uses share the conditional written memory. Joint axes are remote challenge decode and origin reply decode.",
        source="https://dovybzrnzxfo1.cloudfront.net/manifest.json",
        clock_seed=131004,
        records=records,
        kernels=kernels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    print(json.dumps(measure(parser.parse_args().archive), indent=1, allow_nan=False))
