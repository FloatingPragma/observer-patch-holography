"""Retire the first live material label using three balanced control dimers."""

import argparse
from functools import partial
import gzip
from itertools import product
import json
from pathlib import Path

from ...apparatus import layout, COMMIT, FAR
from . import dispatcher, measure

SCRATCH = dispatcher.PAIRS[9:]
HERE = Path(__file__).parent


def at_both(bench, pair, active):
    left, right = (bench.read_bit(pair, at, active) for at in (COMMIT, FAR))
    if left != right:
        raise ArithmeticError("Two control-latch readouts disagree")
    return left


def native_mismatches(bench):
    def station(at):
        sensed = []
        for s in (0, 1):
            ports = layout(s)
            ranks = [
                bench.full ^ bench.rail(p, at, bench.full) for p in ports["selector"]
            ]
            for i in range(3):
                mismatches = 0
                for rank, live in enumerate(ranks):
                    if live:
                        mismatches |= bench.compare(
                            (ports["boundary"][0], ports["banks"][rank][i][0]), at, live
                        )
                sensed.append(mismatches)
        return sensed

    first, second = station(COMMIT), station(FAR)
    assert first == second
    return first


def scan(bench, flags, parked, park_parity=False):
    parity, least, greatest = SCRATCH
    full = bench.full
    at = len(bench.reads)
    assert [at_both(bench, p, full) for p in SCRATCH] == [0] * 3
    clicks = native_mismatches(bench)
    if clicks != flags:
        raise ArithmeticError(
            "Direct boundary/word detector differs from decoded mismatch"
        )
    for click in clicks:
        if not park_parity:
            bench.swap([parity], click)
    even_exit, odd_exit = (bench.rail(rail, COMMIT, full) for rail in parity)
    if even_exit ^ odd_exit != full or at_both(bench, parity, full) != odd_exit:
        raise ArithmeticError("Unresolved physical direction exits")
    bench.swap([least], even_exit)
    bench.swap([greatest], odd_exit)
    groups = [0] * 8
    for seeker, order in ((least, range(8)), (greatest, reversed(range(8)))):
        for j in order:
            live = at_both(bench, seeker, full)
            fire = at_both(bench, dispatcher.SHELVES[j], live)
            if not parked:
                bench.swap([dispatcher.SHELVES[j]], fire)
            bench.swap([seeker], fire)
            groups[j] |= fire
        bench.swap([seeker], at_both(bench, seeker, full))
    bench.swap([parity], odd_exit)
    restored = [at_both(bench, p, full) for p in SCRATCH]
    if any(restored):
        raise ArithmeticError("Selection left an armed latch")
    return groups, {
        "scanner": {
            "span": [at, len(bench.reads)],
            "direction": format(odd_exit, "x"),
            "flags": dispatcher.encode(clicks),
            "restored": dispatcher.encode(restored),
        }
    }


def report():
    source = dict(
        scratch=SCRATCH,
        tree=dispatcher.run(list(product(range(3), repeat=8)), 8, scanner=scan),
        masks=[dispatcher.mask_table(q, scanner=scan) for q in (0, 8)],
        parked=dispatcher.run([(0,)], 8, park=True, scanner=scan),
        parity_off=dispatcher.run([(0,)], 8, scanner=partial(scan, park_parity=True)),
    )
    reference = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    for seq, old in [
        (source["tree"]["steps"], reference["tree"]["steps"]),
        (source["parked"]["steps"], reference["parked"]["steps"]),
        *[
            ([a["proposal"]], [b["proposal"]])
            for a, b in zip(source["masks"], reference["mask_probes"])
        ],
    ]:
        for now, prior in zip(seq, old, strict=True):
            assert {
                k: v for k, v in now.items() if k not in ("scanner", "read_span")
            } == {k: v for k, v in prior.items() if k != "read_span"}
    off = source["parity_off"]["steps"][0]
    assert measure.word(off["groups"], 0) == 1 and measure.word(off["before"], 0) == 255
    assert measure.word(off["flags"], 0) == 8 and off["scanner"]["direction"] == "0"
    source["summary"] = measure.decode(source["tree"])
    assert source["summary"] == reference["summary"]
    return source


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    data = report()
    raw = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    path = HERE / "scan_records.json.gz"
    if write:
        path.write_bytes(gzip.compress(raw, compresslevel=6, mtime=0))
    elif gzip.decompress(path.read_bytes()) != raw:
        raise ValueError("Regenerated records differ from the archived receiver bank")
    print(data["summary"]["distinct_orders"], data["summary"]["prefix_means"][-1])
