"""Local, nondestructive publication of one measured outcome to a tree of material readers."""

import json
from argparse import ArgumentParser
from collections import Counter, deque
from itertools import product
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.clock import (
    DIGITS,
    EXTRA_DIGITS,
    Registers,
    _dimer_pulse,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
WORDS = ("not received", "negative", "positive", "recorded refusal")
RETAINED = (DIGITS[0], EXTRA_DIGITS[0])
OUTBOX = (DIGITS[2], EXTRA_DIGITS[1])
SHOTS = 65536


def rails(cell, layout):
    return [tuple(12 * cell + x for x in triple) for triple in layout]


def tree(board):
    neighbors = [set() for _ in range(board.cells)]
    for a, b in board.external:
        u, v = a // 12, b // 12
        neighbors[u].add(v)
        neighbors[v].add(u)
    queue = deque([0])
    visited = {0}
    edges = []
    while queue:
        parent = queue.popleft()
        for child in sorted(neighbors[parent] - visited):
            queue.append(child)
            visited.add(child)
            edges.append((parent, child))
    if len(visited) != board.cells:
        raise ValueError("the reader graph is disconnected")
    return edges


class Draws:
    """Enumerate the two sides of each verified fair read threshold."""

    def __init__(self, bits):
        self.draws = iter(bits)

    def random(self):
        return 0.25 if next(self.draws) else 0.75


def enroll(board):
    device = Registers(board, 21)
    for cell in range(board.cells):
        device.initialize_digits(rails(cell, RETAINED + OUTBOX))
        device.loads[12 * cell + 2] = 2
    return device


def pointer(device, cell):
    accepted = device.fair(cell, {cell})
    positive = device.fair(cell, {cell})
    return (1 + positive) if accepted else 3


def publish(board, edges, bits, controlled):
    meter = enroll(board)
    if not all(abs(q - 0.5) < 1e-12 for q in meter.rates.values()):
        raise ValueError("the source is outside its fair calibration")
    meter.rng = Draws(bits)
    source = pointer(meter, 0)
    meter.put(rails(0, RETAINED), source, {0})
    inventory = Counter(meter.loads)
    leaf = edges[-1][1]
    for parent, child in edges:
        value = meter.get(rails(parent, RETAINED), {parent})
        if value == 0 or (controlled == "cut" and (parent, child) == edges[0]):
            continue
        meter.put(rails(parent, OUTBOX), value, {parent})
        patch = {parent, child}
        for left, right in zip(
            rails(parent, OUTBOX), rails(child, RETAINED), strict=True
        ):
            for i in (0, 1):
                for u, v in meter.route(left[i], right[i], patch):
                    meter.swap(u, v, patch)
    if controlled == "independent":
        # The comparison leaf destroys its received word, then runs a fresh local
        # instrument with the same calibrated input law and error/refusal menu.
        fresh = pointer(meter, leaf)
        meter.put(rails(leaf, RETAINED), fresh, {leaf})
    reads = [meter.get(rails(i, RETAINED), {i}) for i in range(board.cells)]
    # The source's second fair read is preserved on its own rails throughout.
    reference = int(meter.flag(3, 4, {0}))
    if inventory != Counter(meter.loads):
        raise AssertionError("transport changed the ternary inventory")
    return (source, reference, tuple(reads), int(meter.gates["intercell_pi_contacts"]))


def panel(board, edges, mode, seed):
    nbits = 4 if mode == "independent" else 2
    alternatives = list(product((0, 1), repeat=nbits))
    counts = np.random.default_rng(seed).multinomial(
        SHOTS, [1 / len(alternatives)] * len(alternatives)
    )
    exact = np.zeros((4, 4), dtype=int)
    observed = np.zeros((4, 4), dtype=int)
    side = np.zeros((2, 4, 4), dtype=int)
    complete = disagreements = 0
    missing_sets = set()
    traffic = []
    for bits, n in zip(alternatives, counts, strict=True):
        src, ref, outcomes, contacts = publish(board, edges, bits, mode)
        terminal = outcomes[edges[-1][1]]
        exact[src, terminal] += 1
        observed[src, terminal] += n
        side[ref, src, terminal] += n
        missing_sets.add(tuple(i for i, word in enumerate(outcomes) if word == 0))
        complete += int(n) * all(word != 0 for word in outcomes)
        disagreements += int(n) * any(word not in (0, src) for word in outcomes)
        traffic.append(contacts)
    assert len(missing_sets) == 1
    return dict(
        equally_likely_histories=len(alternatives),
        exact_source_last=exact.tolist(),
        source_last=observed.tolist(),
        private_reference_source_last=side.tolist(),
        readers_without_a_word=list(missing_sets.pop()),
        delivered_everywhere=complete,
        starts_with_conflicting_readings=disagreements,
        intercell_pi_contacts=sorted(set(traffic)),
    )


def private_phase_control():
    # A half contact prepares a superposition in the proved two-rail subspace.
    # After the dwell, the same contact run for the complementary time unprepares it.
    launch = _dimer_pulse(np.pi / 2)
    close = _dimer_pulse(3 * np.pi / 2)
    detectors = {}
    rng = np.random.default_rng(201935)
    for source in range(2):
        branch = launch[:, source]
        coherent = abs((close @ branch)[source]) ** 2
        # Read the occupied rail once and retain both orthogonal records.
        # The returned counts sum the two disjoint measured histories.
        resolved = sum(abs(close[source, j] * branch[j]) ** 2 for j in range(2))
        for setting, prob in [("unread", coherent), ("locally_read", resolved)]:
            prob = float(np.clip(prob, 0, 1))
            restored = int(rng.binomial(SHOTS, prob))
            detectors[f"{source}/{setting}"] = dict(
                restore_probability=round(prob, 12),
                restored=restored,
                changed=SHOTS - restored,
            )
    return detectors


def run():
    boards = []
    for level in (0, 1, 2):
        board = federation(level)
        edges = tree(board)
        # Evaluate each material branch. Sample terminal banks from the exact
        # Born masses once per history, never independently at remote readers.
        panels = {
            mode: panel(board, edges, mode, 71234 + 3 * level + j)
            for j, mode in enumerate(("broadcast", "cut", "independent"))
        }
        boards.append(
            dict(
                level=level,
                owners=board.cells,
                edges=edges,
                last_reader=edges[-1][1],
                cut_edge=edges[0],
                panels=panels,
            )
        )
    return dict(
        words=WORDS,
        starts_per_panel=SHOTS,
        boards=boards,
        private_phase=private_phase_control(),
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
