"""Acquire joint record tables, verify against the frozen packet, and certify every aperture."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from .apparatus import PORTS, Checkpoint, primary_read

HERE = Path(__file__).resolve().parent
PACKET = (
    HERE.parents[4]
    / "code/capacity_readback/runtime/source_derived_public_checkpoint_packet.json"
)
SLOTS = [f"{p}/{orientation}" for p in PORTS for orientation in ("write", "check")]
HELD = "r3_s1_a1_f1"
IDENTITY = "r0_s0_a0_f0"


def bits_word(bits):
    return sum(value << i for i, value in enumerate(bits))


def bits(mask):
    return tuple((mask >> i) & 1 for i in range(12))


def acquire(level, settings):
    machine = Checkpoint(level)
    tables, controls = {}, {}
    for setting in settings:
        outcomes = []
        for label in SLOTS:
            machine.prepare(label)
            if primary_read(machine.reads()) != label:
                raise ValueError("the source failed its sharp read")
            baseline = machine.idle_snapshot()
            machine.apply(setting)
            machine.audit_closure(baseline)
            outcomes.append(bits_word(machine.reads()))
        tables[setting] = outcomes
    if level == 0:
        service, interrupted, parked = [], [], []
        for label in SLOTS:
            machine.prepare(label)
            machine.apply(HELD)
            service.append(machine.relay())
            interrupted.append(machine.relay(omit=True))
            machine.prepare(label)
            machine.apply(HELD, parked=True)
            parked.append(bits_word(machine.reads()))
        controls = dict(
            service=service, first_handoff_omitted=interrupted, parked_routes=parked
        )
    return dict(
        owners=machine.owners,
        board_cells=machine.board.cells,
        complete_words=tables,
        controls=controls,
    )


def certify(tables, packet):
    """Use independent deposited joint atoms and their composition table for the predictions."""
    aliases = packet["public_section_aliases"]
    slots = {section: slot for slot, section in aliases.items()}
    actual = {key: [primary_read(bits(v)) for v in row] for key, row in tables.items()}
    expected = {}
    for channel in packet["global_checkpoint_kernels"]:
        if [channel["authorized_observers"]] != [list(PORTS)]:
            raise ValueError("changed authorization in a deposited channel")
        answers = []
        for label in SLOTS:
            row = channel["rows"][aliases[label]]
            if len(row) != 1 or next(iter(row.values())) != 1:
                raise ValueError("expected one certain joint section")
            answers.append(slots[next(iter(row))])
        expected[channel["continuation_id"]] = answers
    if actual != expected:
        raise ValueError("material joint data differ from the fixed checkpoint")
    compositions = packet["support_relation_composition"]["table"]
    for left, by_right in compositions.items():
        for right, product in by_right.items():
            computed = [actual[left][SLOTS.index(label)] for label in actual[right]]
            if computed != actual[product]:
                raise ValueError("reported continuation product failed")
    return sum(len(row) for row in compositions.values())


def access(tables):
    """For a pair, the identity of every agreeing site determines all confusable apertures."""
    # max_agreement[i,j] is over reported settings (the same setting in the two runs).
    max_agreement = {}
    for i, j in combinations(range(24), 2):
        max_agreement[i, j] = max(
            12 - (row[i] ^ row[j]).bit_count() for row in tables.values()
        )
    # Evaluate every named subset directly as an independent check at identity.
    identity = tables[IDENTITY]
    by_size = {k: Counter() for k in range(1, 13)}
    for mask in range(1, 1 << 12):
        by_size[mask.bit_count()][len({word & mask for word in identity})] += 1
    receipts = []
    for size in range(1, 13):
        conflicts = {p for p, overlap in max_agreement.items() if overlap >= size}
        # A lower bound needs pairwise distinguishable codewords. An upper
        # bound needs a partition into that many pairwise-confusable cliques.
        if size <= 10:
            code, cliques = [0, 1], [list(range(0, 24, 2)), list(range(1, 24, 2))]
        else:
            code, cliques = list(range(24)), [[i] for i in range(24)]
        if (
            any(pair in conflicts for pair in combinations(code, 2))
            or any(
                pair not in conflicts for c in cliques for pair in combinations(c, 2)
            )
            or sorted(i for c in cliques for i in c) != list(range(24))
            or len(code) != len(cliques)
        ):
            raise ValueError("access certificate disagrees with the measured records")
        receipts.append(
            dict(
                reader_count=size,
                maximum_code_size=len(code),
                code=code,
                clique_cover=cliques,
                confusable_input_pairs=len(conflicts),
                fixed_identity_apertures_by_distinct_outputs={
                    str(k): v for k, v in sorted(by_size[size].items())
                },
            )
        )
    return receipts


def read_reference():
    """Use the exact packet at the protocol's named source revision ab22d7439."""
    data = PACKET.read_bytes()
    if (
        sha256(data).hexdigest()
        != "3604147f715327abdc96ad56160fa5ff5384af767d0a57c49beae5b92f3ed96a"
    ):
        raise ValueError("the frozen collective checkpoint packet changed")
    return json.loads(data)


def run():
    packet = read_reference()
    policy = packet["publicness_policy"]
    if (
        policy != [list(PORTS)]
        or packet["publicness_policy_id"] != "universal-twelve-port-publicness/v1"
    ):
        raise ValueError("the pinned collective family has changed")
    settings = sorted(
        chan["continuation_id"] for chan in packet["global_checkpoint_kernels"]
    )
    boards = [acquire(level, settings) for level in range(3)]
    for bank in boards:
        bank["checked_compositions"] = certify(bank["complete_words"], packet)
    ref = boards[0]["complete_words"]
    if any(bank["complete_words"] != ref for bank in boards[1:]):
        raise ValueError("the receiving boards disagree")
    controls = boards[0]["controls"]
    expected = [primary_read(bits(word)) for word in ref[HELD]]
    for i, word in enumerate(controls["service"]):
        for n, local in enumerate(word):
            atom = packet["local_checkpoint_packets"][HELD][PORTS[n]][SLOTS[i]]
            if atom != {f"{PORTS[n]}::record::{primary_read(local)}": 1.0}:
                raise ValueError("the local full-record relay failed")
    refused = Counter(
        primary_read(local) is None
        for word in controls["first_handoff_omitted"]
        for local in word
    )
    if (
        not refused[True]
        or not refused[False]
        or any(
            primary_read(local) not in (None, expected[i])
            for i, word in enumerate(controls["first_handoff_omitted"])
            for local in word
        )
    ):
        raise ValueError("missing-word control did not isolate delivery")
    wrong = sum(
        primary_read(bits(word)) != label
        for word, label in zip(controls["parked_routes"], expected, strict=True)
    )
    if wrong == 0:
        raise ValueError("parked continuation was indistinguishable")
    return dict(
        policy_id=packet["publicness_policy_id"],
        authorized_sets=policy,
        slot_order=SLOTS,
        source_packet=packet["packet_sha256"],
        boards=boards,
        interfaces=access(ref),
        omitted_round_refusals=refused[True],
        omitted_round_complete=refused[False],
        held_parked_disagreements=wrong,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    serialized = json.dumps(run(), indent=1) + "\n"
    if parser.parse_args().write:
        (HERE / "records.json").write_text(serialized)
    else:
        print(serialized, end="")
