"""Acquire an unread physical selector and compare the resulting checkpoint kernels."""

import argparse
from functools import cache
import gzip
import json
from math import log
from pathlib import Path
import sys

import numpy as np

from .apparatus import (
    READ_CODES,
    calibrate,
    decode,
    exchange_word,
    federation,
    pulse,
    rail_swaps,
    raw_word,
)
from .checkpoint import correct
from ..seam_word_meter.finite_program.machine import digit_instruments

HERE = Path(__file__).parent
SHOTS = 16384
SEED = 1262605
ARMS = (("one", 3), ("internal", 7), ("boundary", 8), ("parked", 0))
BINS = 262144


def published(stem):
    return json.loads(gzip.decompress((HERE / (stem + ".json.gz")).read_bytes()))


def layout_selectors(cfg):
    board = federation(0)
    excluded = set(cfg["owners"] + cfg["auditors"] + [cfg["boundary"][1] // 12])
    selected = []
    used = set()
    for a, b in sorted(tuple(sorted(e)) for e in board.seams):
        if a // 12 != b // 12 or a // 12 in excluded or {a, b} & used:
            continue
        if board.degree[a] == board.degree[b] == 5:
            selected.append((a, b))
            used.update((a, b))
        if len(selected) == 8:
            break
    if len(selected) != 8:
        raise ValueError("Missing independently mounted selector dimers")
    return selected


def audited_labels(cfg, dimers):
    owners = sorted({a // 12 for a, b in dimers})
    stations = calibrate({**cfg, "owners": owners})
    ports = {p for edge in dimers for p in edge}
    for station in stations:
        station["paths"] = [r for r in station["paths"] if r["port"] in ports]
    read, _ = digit_instruments()
    routes = [
        {r["port"]: exchange_word(r["path"]) for r in s["paths"]} for s in stations
    ]
    charts = []
    for physical in range(256):
        loads = [0] * federation(0).n
        for station in stations:
            lo, hi = station["references"]
            loads[lo], loads[hi] = 0, 2
        for i, (a, b) in enumerate(dimers):
            bit = (physical >> i) & 1
            loads[a], loads[b] = 1 - bit, bit
        original = loads.copy()
        tapes = []
        for head, routing in zip(stations, routes):
            flags = []
            for pair in dimers:
                for rail in pair:
                    path = routing[rail]
                    for a, b in path:
                        loads[a], loads[b] = loads[b], loads[a]
                    value = loads[head["head"]]
                    zero, double = (loads[r] for r in head["references"])
                    flags.append(2 * read[value, zero] + read[value, double])
                    for a, b in reversed(path):
                        loads[a], loads[b] = loads[b], loads[a]
            tapes.append(flags)
        if loads != original or tapes[0] != tapes[1]:
            raise ArithmeticError("Reading selector did not restore material")
        trits = [READ_CODES.index(k) for k in tapes[0]]
        if any(trits[i] + trits[i + 1] != 1 for i in range(0, 16, 2)):
            raise ArithmeticError("Unresolved selection")
        observed = sum(v << j for j, v in enumerate(trits[1::2]))
        charts.append(dict(raw=tapes, decoded=observed))
    return stations, charts


def native_plate():
    from model.candidate_laws.interval_echo import candidate

    law = candidate.native_law
    a, b = [law.SRC.index(v) for v in ((1, 0), (0, 1))]
    rail_swaps()
    transfers = {}
    for name, t in (("parked", 0), ("half", np.pi / 2), ("full", np.pi)):
        u = pulse(t)
        columns = abs(u[:, [a, b]]) ** 2
        np.testing.assert_allclose(columns[[a, b]].sum(axis=0), 1, atol=1e-12)
        transfers[name] = np.round(columns[[a, b]], 14).tolist()
    return transfers


def apply_control(cfg, raw_states, stations, labels, arm, width):
    # Apply the selected full exchanges directly to both physical endpoints
    # of every sharp source column. A full native exchange acts on *values*.
    n = len(cfg["owners"])
    ports = [12 * o + j for o in cfg["owners"] for j in range(12)] + [
        cfg["boundary"][1]
    ]
    edges = [
        (ports.index(a), ports.index(b)) for a, b in [*cfg["inside"], cfg["boundary"]]
    ]
    decode_meter = np.array([READ_CODES.index(v) for v in range(3)])
    terminals = []
    for pattern in range(1 << width):
        code = labels[pattern]["decoded"]
        if arm == "one":
            enabled = [] if code == 0 else [code - 1]
        else:
            enabled = [j for j in range(width) if code & (1 << j)]
        loads = raw_states.copy()
        for j in enabled:
            a, b = edges[j]
            loads[:, [a, b]] = loads[:, [b, a]]
        raw = [np.take(station["sharp"], loads) for station in stations]
        if not np.array_equal(*raw):
            raise ArithmeticError("Two terminal heads disagree")
        trits = decode_meter[raw[0]]
        word = (trits[:, :-1].reshape(256, n, 12).sum(axis=2) % 2) @ (1 << np.arange(n))
        terminals.append(word)
    return np.asarray(terminals).T


@cache
def evaluator():
    directory = str(HERE.parents[3] / "code" / "capacity_readback")
    sys.path.insert(0, directory)
    try:
        from correctable_public_record_capacity import compound_confusability_graph
    finally:
        sys.path.remove(directory)
    return compound_confusability_graph


def graphs(words, destination, arm):
    # Matrix rows use the sharp material input index. Public graph vertices
    # and counted outputs are the measured parity words only.
    labels = [str(v) for v in words]
    rows = {
        x: {str(int(y)): 1 / len(outs) for y in outs}
        for x, outs in zip(labels, destination)
    }
    if any(len(r) != destination.shape[1] for r in rows.values()):
        raise ArithmeticError("Physical label fiber collision")
    graph = evaluator()(
        labels, [{"name": "unreported physical selector", "rows": rows}]
    )
    indexed = arm == "one"
    if indexed:
        known = [
            {
                "name": f"known {j}",
                "rows": {x: {str(y): 1} for x, y in zip(labels, col)},
            }
            for j, col in enumerate(destination.T)
        ]
        if any(evaluator()(labels, known).values()):
            raise ArithmeticError("Known permutations lost a record")
    # In each arm the measured decoder controls give disjoint output sets and
    # each source uses the full registered number of possible receipts.
    if arm == "one":
        chosen = [i for i in range(256) if correct(i % 128) == i % 128]
    elif arm == "internal":
        chosen = [0, 128]
    elif arm == "boundary":
        chosen = [0]
    else:
        chosen = list(range(256))
    if any(
        str(words[j]) in graph[str(words[i])]
        for p, i in enumerate(chosen)
        for j in chosen[p + 1 :]
    ):
        raise ArithmeticError("Claimed code has a terminal collision")
    fixed = list(range(128))
    fixed_code = [i for i in chosen if i in fixed]
    fixed_outputs = set(destination[fixed].ravel().tolist())
    if (
        len(fixed_code) * destination.shape[1] != len(fixed_outputs)
        or len(chosen) * destination.shape[1] != 256
    ):
        raise ArithmeticError("Code misses the physical output packing bound")
    if any(len(set(col)) != 256 for col in destination.T):
        raise ArithmeticError("Revealing the material selector lost an input")
    return {
        "edges": sum(map(len, graph.values())) // 2,
        "code_words": [words[i] for i in chosen],
        "fixed_class_code_words": [words[i] for i in fixed_code],
        "fixed_class_distinct_receipts": len(fixed_outputs),
        "receipts_per_source": destination.shape[1],
        "known_indexed_one_swap_sizes": [128, 256] if indexed else None,
        "visible_log_sizes": [128, 256],
    }


def acquire():
    checkpoint = published("checkpoint")
    cfg = checkpoint["config"]
    dimers = layout_selectors(cfg)
    stations, labels = audited_labels(cfg, dimers)
    plate = native_plate()
    original = np.array([raw_word(cfg, i % 128, i // 128) for i in range(256)])
    words = [decode(raw, 8)[0] for raw in checkpoint["raw"]]
    if len(set(words)) != 256 or any(r["decoded"] != j for j, r in enumerate(labels)):
        raise ArithmeticError("Source or selector codebook incomplete")
    rng = np.random.default_rng(SEED)
    panels = {}
    checked = 0
    for name, width in ARMS:
        destination = apply_control(
            cfg, original, checkpoint["stations"], labels, name, width
        )
        # The isolated pulse and its sharp restored read act as a tensor
        # product over the disjoint selector dimers. Unused dimers stay put.
        size = 1 << width
        prior = np.array(
            [
                np.prod([plate["half"][(x >> j) & 1][0] for j in range(width)])
                for x in range(size)
            ]
        )
        assert np.array_equal(prior, np.full(size, 1 / size))
        samples = rng.multinomial(SHOTS, prior, size=256)
        terminal_counts = np.zeros((256, 256), int)
        np.add.at(terminal_counts, (np.arange(256)[:, None], destination), samples)
        expected = np.zeros_like(terminal_counts, dtype=float)
        np.add.at(expected, (np.arange(256)[:, None], destination), prior)
        populated = expected > 0
        if not np.array_equal(populated, terminal_counts > 0):
            raise ArithmeticError(
                "A permitted output absent or forbidden output observed"
            )
        ell = log(2 * BINS / 0.01)
        allowance = np.sqrt(2 * SHOTS * expected * (1 - expected) * ell) + 2 * ell / 3
        ratio = float(np.max(abs(terminal_counts - SHOTS * expected) / allowance))
        if ratio > 1:
            raise ArithmeticError(
                "A terminal column fails the frozen simultaneous bound"
            )
        checked += terminal_counts.size
        panels[name] = dict(
            width=width,
            terminal_by_measured_pattern=destination.tolist(),
            counts_by_measured_pattern=samples.tolist(),
            terminal_counts=terminal_counts.tolist(),
            max_allowance_fraction=round(ratio, 8),
            minimum_populated_count=int(terminal_counts[populated].min()),
            certificate=graphs(words, destination, name),
        )
    assert checked == BINS
    return dict(
        seed=SEED,
        shots_per_input_arm=SHOTS,
        checked_terminal_bins=checked,
        dimers=dimers,
        station_routes=stations,
        selector_truth_table=labels,
        plate=plate,
        public_inputs=words,
        provenance="Seeded multinomial draws from the native isolated selector product law. Both restored heads read the same sharp selector outcome. The terminal swap circuit is conditioned on those readings; a hidden receiver sees only cell parities.",
        arms=panels,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    write = cli.parse_args().write
    data = acquire()
    encoded = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    path = HERE / "selector.json.gz"
    if write:
        path.write_bytes(gzip.compress(encoded, compresslevel=6, mtime=0))
    elif gzip.decompress(path.read_bytes()) != encoded:
        raise ValueError("Material selector receiver banks differ")
    print(
        [
            (
                k,
                len(v["certificate"]["fixed_class_code_words"]),
                len(v["certificate"]["code_words"]),
                v["max_allowance_fraction"],
            )
            for k, v in data["arms"].items()
        ]
    )
