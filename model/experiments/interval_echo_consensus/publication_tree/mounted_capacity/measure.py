"""Acquire all sharp port words and compute zero-error alphabets from joint readouts."""

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from model.executable_model import rotations
from .apparatus import Board, MARKS, PORTS, WORDS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
MENUS = ("unmarked", "host", "edge", "ordered")


def read_planes(raw, width=WORDS):
    full = (1 << width) - 1
    if len(raw) != 36 or any(a ^ b != full or a & b for a, b in raw):
        raise ValueError("Unresolved pair")
    packed = np.array(
        [
            np.frombuffer(b.to_bytes((width + 7) // 8, "little"), dtype=np.uint8)
            for _, b in raw
        ]
    )
    planes = np.unpackbits(packed, axis=1, bitorder="little")[:, :width].reshape(
        3, PORTS, width
    )
    return ((1 << np.arange(PORTS)) @ planes).astype(np.int64)


def anchor(host, tip, menu):
    if menu == "unmarked":
        return 0
    if menu == "host":
        return host
    if menu == "edge":
        return host | tip
    if menu == "ordered":
        return host + WORDS * tip
    raise ValueError(menu)


def measured_alphabet(payloads, anchors):
    keys = payloads + WORDS * np.asarray(anchors, dtype=np.int64)[:, None]
    reps = keys.min(axis=0)
    # A raw word shared by two input messages must identify the same canonical receipt.
    readings, inverse = np.unique(keys, return_inverse=True)
    lookup = np.empty(len(readings), dtype=np.int64)
    lookup[inverse.ravel()] = np.tile(reps, len(keys))
    if not np.all(lookup[inverse.reshape(keys.shape)] == reps):
        raise ValueError("Overlapping supports have distinct candidate representatives")
    # Equal minima imply intersection. Check equality of the entire possible-readout sets.
    by_rep = {}
    for word, rep in enumerate(reps):
        possible = tuple(np.unique(keys[:, word]))
        if by_rep.setdefault(int(rep), possible) != possible:
            raise ValueError("Intersecting but unequal readout supports")
    weights = [
        len(set(reps[[i for i in range(WORDS) if i.bit_count() == k]]))
        for k in range(13)
    ]
    if sum(weights) != len(by_rep):
        raise ValueError("Classes mix input weights")
    return {
        "code_size": len(by_rep),
        "by_weight": weights,
        "support_sizes": dict(sorted(Counter(len(s) for s in by_rep.values()).items())),
    }


def fixed_words(group):
    by_weight = np.zeros(13, dtype=int)
    for g in group:
        fixed = np.zeros(13, dtype=int)
        fixed[0] = 1
        seen = set()
        for i in range(PORTS):
            if i in seen:
                continue
            j = i
            length = 0
            while j not in seen:
                seen.add(j)
                length += 1
                j = g[j]
            fixed[length:] += fixed[:-length].copy()
        by_weight += fixed
    if np.any(by_weight % len(group)):
        raise ValueError("Non-integral orbit count")
    return (by_weight // len(group)).tolist()


def pack(raw):
    return [[format(a, "x"), format(b, "x")] for a, b in raw]


def report():
    group = rotations()
    board = Board()
    initial = board.load(range(WORDS))
    data = []
    marks = []
    receipts = []
    for g in group:
        raw, length = board.acquire(initial, g, WORDS)
        payload, host, tip = read_planes(raw[0])
        if not (np.all(host == host[0]) and np.all(tip == tip[0])):
            raise ValueError("Marker depends on data")
        if host[0] != 1 << g[MARKS[0]] or tip[0] != 1 << g[MARKS[1]]:
            raise ValueError("Deployment mark mismatch")
        data.append(payload)
        marks.append((int(host[0]), int(tip[0])))
        receipts.append({"mount": g, "pulses": length, "raw": [pack(r) for r in raw]})
    data = np.stack(data)
    results = {}
    for menu, target in zip(MENUS, (96, 832, 2080, 4096), strict=True):
        anchors = [anchor(h, t, menu) for h, t in marks]
        reference = anchor(1 << MARKS[0], 1 << MARKS[1], menu)
        stabilizer = [g for g, a in zip(group, anchors, strict=True) if a == reference]
        measured = measured_alphabet(data, anchors)
        counted = fixed_words(stabilizer)
        if measured["by_weight"] != counted or measured["code_size"] != target:
            raise ValueError("Material alphabet and fixed-word census disagree")
        results[menu] = measured | {"stabilizer_order": len(stabilizer)}
    # Independent mark remounts multiply every data support by exactly the same 60 marks.
    if len(set(marks)) != 60:
        raise ValueError("Incomplete independent marker menu")
    separated = measured_alphabet(data, [0] * len(group))
    # A second identical tape embeds each original receipt injectively on the diagonal.
    # Each station actually produced its own full tape. Equality at acquisition,
    # followed by the diagonal embedding theorem, preserves every joint support.
    for r in receipts:
        if r["raw"][0] != r["raw"][1]:
            raise ValueError("Paired stations did not realize a diagonal copy")
    g = next(g for g in group if g[MARKS[0]] != MARKS[0])
    partial, _ = board.acquire(initial, g, WORDS, parked=True)
    blank, _ = board.acquire(initial, g, WORDS, mode="parked")
    if partial[0] == [
        [int(a, 16), int(b, 16)] for a, b in receipts[group.index(g)]["raw"][0]
    ] or any(a or b for a, b in blank[0]):
        raise ValueError("A parked control reproduced its ideal acquisition")
    files = [
        HERE / "PROTOCOL.md",
        HERE / "apparatus.py",
        Path(__file__),
        ROOT / "model/executable_model/symmetry.py",
        ROOT / "model/executable_model/carrier.py",
        ROOT / "model/executable_model/carriers/FEDERATION_LEVEL0.npz",
        ROOT / "model/experiments/interval_echo_consensus/boolean_nodes/material.py",
        ROOT / "model/candidate_laws/interval_echo/candidate.py",
        ROOT / "model/candidate_laws/contrast_memory_spherical_v1/law.py",
    ]
    return {
        "schema": "oph.mounted_capacity.v1",
        "encoding": "hexadecimal; bit k is the sharp click at launch k; pair order is first rail then second",
        "marks": MARKS,
        "auditors": [18, 19],
        "sizes": results,
        "independent_mark_capacity": separated["code_size"],
        "independent_joint_mounts": len(group) ** 2,
        "duplicated_flag_capacity": results["ordered"]["code_size"],
        "sources": {
            str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in files
        },
        "controls": {
            "mount": g,
            "parked_crossings": [pack(r) for r in partial],
            "parked_echo": [pack(r) for r in blank],
        },
        "receipts": receipts,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = report()
    encoded = json.dumps(result, separators=(",", ":")).encode()
    path = HERE / "records.json.gz"
    if args.write:
        path.write_bytes(gzip.compress(encoded, mtime=0))
    if encoded != gzip.decompress(path.read_bytes()):
        raise ValueError("Replay differs from the deposited receipt")
    print({k: v["code_size"] for k, v in result["sizes"].items()})
