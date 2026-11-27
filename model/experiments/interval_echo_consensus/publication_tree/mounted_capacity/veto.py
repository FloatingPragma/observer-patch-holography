"""Measure whether the one-slot veto is determined by each published frame."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from .apparatus import AUDITORS, MARKS, PORTS, WORDS, Board
from .measure import HERE, MENUS, anchor, pack, read_planes

FULL = (1 << WORDS) - 1


def decode_digit(pair):
    a, b = pair
    if a ^ b != FULL or a & b:
        raise ValueError("Invalid read")
    return b


def operate(board, initial, g, mode="native"):
    loads = initial.copy()
    board.swaps(loads, board.program(g, False))
    outputs = board.device.rows[18][:2]
    for a, b in outputs:
        loads[a], loads[b] = 0, FULL
    survey = board.read(
        loads,
        17,
        "parked" if mode == "parked" else "native",
        FULL,
        board.slots[PORTS : 2 * PORTS],
    )
    source_read = None
    if all((a ^ b) == FULL and not (a & b) for a, b in survey):
        marks = [decode_digit(pair) for pair in survey]
        if marks.count(FULL) != 1 or marks.count(0) != PORTS - 1:
            raise ValueError("The sharp host is not uniquely appointed")
        i = MARKS[0] if mode == "fixed" else marks.index(FULL)
        source_read = board.read(loads, 17, "native", FULL, [board.slots[i]])[0]
        occupied = decode_digit(source_read)
        # On each launch a positive pointer condition requests one local pi exchange.
        for (a, b), selected in zip(outputs, (occupied, FULL ^ occupied), strict=True):
            changes = (loads[a] ^ loads[b]) & selected
            loads[a] ^= changes
            loads[b] ^= changes
    station = [board.read(loads, s, "native", FULL, outputs) for s in AUDITORS]
    if station[0] != station[1]:
        raise ValueError("Veto readers disagree")
    return {
        "survey": pack(survey),
        "payload_read": pack([source_read])[0] if source_read else None,
        "output": [pack(r) for r in station],
    }


def bits(hex_value):
    packed = int(hex_value, 16).to_bytes(WORDS // 8, "little")
    return np.unpackbits(np.frombuffer(packed, dtype=np.uint8), bitorder="little")


def public_decisions(payloads, anchors, decisions):
    keys = payloads + WORDS * np.array(anchors, dtype=np.int64)[:, None]
    unique, inv = np.unique(keys, return_inverse=True)
    masks = np.zeros(len(unique), dtype=np.uint8)
    np.bitwise_or.at(masks, inv.ravel(), decisions.ravel())
    at = masks[inv.reshape(keys.shape)]
    if not np.all(at == at[0]):
        raise ValueError("Decision ambiguity varies along a measured orbit")
    _, examples = np.unique(keys.min(axis=0), return_index=True)
    counts = np.bincount(at[0, examples], minlength=4)
    if counts[0]:
        raise ValueError("Public orbit without a decision")
    return dict(zip(("veto", "allow", "ambiguous"), map(int, counts[1:]), strict=True))


def report():
    physical = HERE / "records.json.gz"
    previous = json.loads(gzip.decompress(physical.read_bytes()))
    board = Board()
    initial = board.load(range(WORDS))
    rows = []
    payloads = []
    marks = []
    decisions = []
    for receipt in previous["receipts"]:
        g = tuple(receipt["mount"])
        measured = operate(board, initial, g)
        payload, h, t = read_planes(
            [[int(a, 16), int(b, 16)] for a, b in receipt["raw"][0]]
        )
        flags = [bits(pair[1]) for pair in measured["output"][0]]
        if not np.all(flags[0] + flags[1] == 1):
            raise ValueError("Incomplete publication")
        rows.append(measured)
        payloads.append(payload)
        marks.append((int(h[0]), int(t[0])))
        decisions.append(flags[0] + 2 * flags[1])
    data = np.stack(payloads)
    reported = np.stack(decisions)
    census = {}
    for m, want in zip(
        MENUS,
        [(1, 1, 94), (416, 416, 0), (528, 528, 1024), (2048, 2048, 0)],
        strict=True,
    ):
        result = public_decisions(data, [anchor(h, t, m) for h, t in marks], reported)
        if tuple(result.values()) != want:
            raise ValueError("Measured veto differs from prediction")
        census[m] = result
    k = next(i for i, (h, t) in enumerate(marks) if h != 1 << MARKS[0])
    g = tuple(previous["receipts"][k]["mount"])
    fixed = operate(board, initial, g, "fixed")
    lost = operate(board, initial, g, "parked")
    changed = (
        int(fixed["output"][0][0][1], 16) ^ int(rows[k]["output"][0][0][1], 16)
    ).bit_count()
    if changed != WORDS // 2 or any(int(pair[1], 16) for pair in lost["output"][0]):
        raise ValueError("The source-change or refusal control failed")
    return {
        "census": census,
        "joint_vetoes": int((reported == 1).sum()),
        "rows": rows,
        "controls": {"mount": g, "fixed": fixed, "changed": changed, "parked": lost},
        "sources": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (
                physical,
                Path(__file__),
                HERE / "VETO_PROTOCOL.md",
                HERE / "apparatus.py",
            )
        },
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    d = report()
    p = HERE / "veto.json.gz"
    raw = json.dumps(d, separators=(",", ":")).encode()
    if args.write:
        p.write_bytes(gzip.compress(raw, mtime=0))
    if raw != gzip.decompress(p.read_bytes()):
        raise ValueError("Veto report differs")
    print(d["census"])
