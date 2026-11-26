"""Retain physical contact calendars, all detector outcomes and each receiver's tapes."""

import argparse
from collections import Counter
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path

from . import compiler, hardware
from .compiler import support_inputs
from model.executable_model.carrier import federation
from .hardware import Lane, LocalState, domains

HERE = Path(__file__).resolve().parent
DYNAMIC = ("project", "shutter", "intervene")


def play(board, lanes, inputs, mode="native", serial=False):
    loads = [(i * 7 + 2) % 3 for i in range(board.n)]
    local = [LocalState(mode == "parked") for _ in lanes]
    for lane, (w, _) in zip(lanes, inputs, strict=True):
        for reg, (a, b) in enumerate(lane.pairs):
            value = (w >> reg) & 1 if reg < 3 else int(mode == "dirty" and reg == 8)
            loads[a], loads[b] = 2 * value, 2 * (1 - value)
        for z in lane.refs:
            loads[z] = 0
    before = loads.copy()
    lengths = [len(lane.word) for lane in lanes]
    trace = [dict() for _ in lanes]
    clocks = sum(lengths) if serial else max(lengths)
    shifts = [sum(lengths[:i]) if serial else 0 for i in range(len(lanes))]
    for t in range(clocks):
        for i, lane in enumerate(lanes):
            j = t - shifts[i]
            if 0 <= j < lengths[i]:
                lamp = local[i].step(loads, lane.word[j], inputs[i][1])
                if lamp is not None:
                    trace[i][j] = lamp
    checked = []
    writable = {p for lane in lanes for pair in lane.pairs[3:9] for p in pair}
    if any(
        a != b
        for i, (a, b) in enumerate(zip(before, loads, strict=True))
        if i not in writable
    ):
        raise ArithmeticError("An untouched source, reference or transit cache changed")
    for i, lane in enumerate(lanes):
        if local[i].saved is not None or local[i].pending is not None:
            raise ArithmeticError("Leaked collar")
        outcomes = trace[i]
        reads = {}
        for label, at, reg, lo, hi in lane.reads:
            if outcomes[lo] + outcomes[hi] != 1:
                raise ArithmeticError("Unbalanced receiving rails")
            bank = reads.setdefault(label, [{}, {}])
            bank[at][reg] = outcomes[lo]
        if any(a != b for a, b in reads.values()):
            raise ArithmeticError("Independent stations disagree")
        counts = Counter(op[0] for op in lane.word)
        openings = sum(
            outcomes[j]
            for j, op in enumerate(lane.word)
            if op[0] in ("shutter", "intervene")
        )
        costs = dict(
            slots=lengths[i],
            native_exchanges=counts["swap"] + counts["flip"] + openings,
            echo_or_inverse=sum(
                counts[k] for k in ("write", "unwrite", "probe_write", "probe_undo")
            ),
            projections=counts["project"],
            conditional_slots=counts["shutter"] + counts["intervene"],
            opened=openings,
            waits=counts["wait"],
        )
        checked.append(
            dict(
                flags="".join(str(x) for x in outcomes.values()),
                readbacks={k: list(a.values()) for k, (a, _) in reads.items()},
                costs=costs,
            )
        )
    return dict(inputs=inputs, calendar=clocks, branches=checked)


def report():
    records = []
    for level in (0, 1):
        board = federation(level)
        groups = domains(board)
        lanes = tuple(Lane(board, g).program() for g in groups)
        if not set(lanes[0].graph).isdisjoint(lanes[1].graph):
            raise ValueError("Parallel schedules would contend for physical hardware")
        bank = [
            play(board, lanes, inputs) for inputs in product(support_inputs(), repeat=2)
        ]
        for launch in bank:
            for (w, edit), result in zip(
                launch["inputs"], launch["branches"], strict=True
            ):
                bit = int(w.bit_count() >= 2)
                after = bit ^ edit
                if result["readbacks"]["end"] != [
                    *((w >> i) & 1 for i in range(3)),
                    after,
                    bit,
                    after,
                    1,
                    after ^ ((w >> 2) & 1),
                    0,
                ]:
                    raise ArithmeticError(
                        "An intervened truth-table circuit missed its native prediction"
                    )
        serial = play(board, lanes, ((2, 1), (7, 0)), serial=True)
        altered = tuple(Lane(board, g).program(stale=True) for g in groups)
        controls = {
            name: [
                play(
                    board,
                    altered if name == "stale" else lanes,
                    (condition, condition),
                    mode=name,
                )
                for condition in support_inputs()
            ]
            for name in ("stale", "parked", "dirty")
        }
        records.append(
            dict(
                level=level,
                lanes=[lane.address_map() for lane in lanes],
                bank=bank,
                serial=serial,
                stale_lanes=[lane.address_map() for lane in altered],
                controls=controls,
            )
        )
    from ..seam_word_meter.finite_program import machine

    sources = [
        HERE / "PROTOCOL.md",
        Path(__file__),
        Path(compiler.__file__),
        Path(hardware.__file__),
        Path(machine.__file__),
        Path(machine.LAW.__file__),
    ]
    return dict(
        provenance="Complete deterministic simulated native dual-rail observations with sharp projected collars.",
        trace_order="Left-to-right flags for every project, shutter and intervene operation in each physical word.",
        sources={
            str(f.resolve().relative_to(HERE.parents[3])): hashlib.sha256(
                f.read_bytes()
            ).hexdigest()
            for f in sources
        },
        records=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if options.output:
        with options.output.open("xb") as handle:
            handle.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != data:
        raise ValueError("Native schedules or observations differ from archive")
