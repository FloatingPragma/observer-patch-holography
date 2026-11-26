"""Measure all root interventions through one persistent 135-appointment program."""

import argparse
from functools import cache
import gzip
import hashlib
import json
from pathlib import Path

from model.executable_model import carrier
from .....seam_word_meter.finite_program.layout import mount
from .....seam_word_meter.finite_program import machine
from .. import apparatus

HERE = Path(__file__).resolve().parent
ROOTS = 7
MASKS = 1 << ROOTS
SHELVES = ROOTS + MASKS
WIDTH = 9
BASIS = (0,) + tuple(1 << i for i in range(ROOTS))


def reading(raw):
    if len(raw) != WIDTH or any(n not in "69" for n in raw):
        return None
    return sum((n == "9") << k for k, n in enumerate(raw))


def measure(bench, word, mode):
    bench.start([0] * SHELVES)
    events = []
    for root in range(ROOTS):
        tag = root + 1
        published = apparatus.store(bench, root, (tag << 1) | (word >> root & 1))
        events.append([root, [], published])
    for mask in range(MASKS):
        reg = ROOTS + mask
        actual = 33 if mode == "altered" and mask == 65 else mask
        stamp_fault = mode == "stamp" and mask == 127
        if stamp_fault:
            bench.swaps(bench.loads, [bench.words[0][1]])
        samples = []
        parity = 0
        valid = True
        for i in range(ROOTS):
            if actual & (1 << i):
                a, first = apparatus.read(bench, i, 0)
                parked = mode == "shutter" and mask == 127 and i == 0
                b, second = apparatus.read(bench, i, 1, parked=parked)
                samples.append([i, first, second])
                if a is None or a != b or a >> 1 != i + 1:
                    valid = False
                else:
                    parity ^= a & 1
        published = (
            apparatus.store(bench, reg, ((reg + 1) << 1) | parity) if valid else None
        )
        events.append([reg, samples, published])
        if stamp_fault:
            bench.swaps(bench.loads, [bench.words[0][1]])
    final = [[apparatus.read(bench, a, at)[1] for a in range(SHELVES)] for at in (0, 1)]
    if final[0] != final[1] or bench.tape:
        raise ValueError("Final read disagreement")
    return dict(word=word, mode=mode, events=events, final=final)


def decode_log(log, physical):
    # Verify each read against its parent's actual deposited word and targeted register.
    snapshots = {}
    past = []
    bits = []
    if len(log["events"]) != SHELVES or len(log["final"]) != 2:
        raise ValueError("Incomplete run")
    for eid, (reg, reads, store) in enumerate(log["events"], 1):
        if reg != eid - 1 or reg in snapshots or store is None:
            raise ValueError("Unsealed or reused target")
        source_mask = 0
        for address, first, second in reads:
            observed = reading(first)
            if (
                address not in snapshots
                or address >= ROOTS
                or source_mask & (1 << address)
            ):
                raise ValueError("Unexpected input origin")
            if first != second or observed is None:
                raise ValueError("Source readers disagree")
            parent_reg, parent_word = snapshots[address]
            if (
                parent_reg != address
                or parent_word != observed
                or observed >> 1 != parent_reg + 1
            ):
                raise ValueError(
                    "Source is not the authentic writer and deposited value"
                )
            source_mask |= 1 << address
        raw = store["after"]
        value = reading(raw)
        if value is None or value >> 1 != eid or store["peers"] != [raw, raw]:
            raise ValueError("Independent material publication failed")
        if reading(store["before"]) != 0 or any(
            bank[reg] != raw for bank in log["final"]
        ):
            raise ValueError("Mutable or missing final record")
        toggled = {
            p
            for k, pair in enumerate(physical["words"][reg])
            if value >> k & 1
            for p in pair
        }
        if set(store["changed_ports"]) != toggled:
            raise ValueError("Frame is not the sealed store's physical footprint")
        snapshots[reg] = (reg, value)
        past.append(source_mask)
        bits.append(value & 1)
    if any(past[:ROOTS]) or log["word"] != sum(
        v << i for i, v in enumerate(bits[:ROOTS])
    ):
        raise ValueError("Root preparation or independence failed")
    return bits[ROOTS:], past[ROOTS:]


def qualify(logs, physical, full):
    table, ancestry = {}, {}
    for log in logs:
        word = log["word"]
        if word in table:
            raise ValueError("Duplicate input")
        table[word], ancestry[word] = decode_log(log, physical)
    if set(table) != set(range(MASKS) if full else BASIS):
        raise ValueError("Incomplete source interventions")
    pivotal = [
        sum(1 << i for i in range(ROOTS) if table[1 << i][m] != table[0][m])
        for m in range(MASKS)
    ]
    if any(ancestry[w] != pivotal for w in table):
        raise ValueError("The physical read graph disagrees with pivotality")
    if any(table[0]):
        raise ValueError("Blank source has a bright parity")
    flips = 0
    for w, row in table.items():
        for m, support in enumerate(pivotal):
            if row[m] != (w & support).bit_count() % 2:
                raise ValueError("Parity and material measurement differ")
            if full:
                for i in range(ROOTS):
                    differs = table[w ^ (1 << i)][m] != row[m]
                    if differs != bool(support & (1 << i)):
                        raise ValueError("Full conditional intervention differs")
                    flips += int(differs)
    return dict(
        inputs=len(table),
        roots=ROOTS,
        distinct_root_pasts=len(set(pivotal)),
        expected_parents=sum(s.bit_count() for s in pivotal),
        pivotal=pivotal,
        verified_positive_flips=flips if full else sum(s.bit_count() for s in pivotal),
    )


def report():
    bench = mount.Mount(carrier.federation(2)).bench(
        registers=SHELVES, axes=1, width=WIDTH
    )
    bench.route = cache(bench.route)
    physical = apparatus.layout(bench)
    order = [*BASIS, *(i for i in range(MASKS) if i not in BASIS)]
    native = [measure(bench, word, "native") for word in order]
    altered = [measure(bench, word, "altered") for word in BASIS]
    controls = [measure(bench, 65, mode) for mode in ("shutter", "stamp")]
    analysis = dict(
        native=qualify(native, physical, True),
        altered=qualify(altered, physical, False),
    )
    accepted = analysis["native"]
    if accepted["pivotal"] != list(range(MASKS)) or accepted["expected_parents"] != 448:
        raise ArithmeticError("Frozen native all-cuts prediction failed")
    if analysis["altered"]["pivotal"] != [33 if m == 65 else m for m in range(MASKS)]:
        raise ArithmeticError("Changed-root control failed")
    baseline = next(launch for launch in native if launch["word"] == 65)
    for log in controls:
        if (
            log["events"][:-1] != baseline["events"][:-1]
            or log["events"][-1][2] is not None
            or log["final"] != [baseline["final"][0][:-1] + ["6" * WIDTH]] * 2
        ):
            raise ArithmeticError("Held control damaged earlier shelves or published")
    pins = [
        HERE / "PROTOCOL.md",
        Path(__file__),
        Path(apparatus.__file__),
        Path(machine.__file__),
        Path(mount.__file__),
        Path(machine.LAW.__file__),
    ]
    return dict(
        provenance="Complete deterministic simulated native sharp material readouts and checked frame changes.",
        sources={
            str(f.resolve().relative_to(HERE.parents[6])): hashlib.sha256(
                f.read_bytes()
            ).hexdigest()
            for f in pins
        },
        mount=physical,
        native=native,
        altered=altered,
        controls=controls,
        analysis=analysis,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if options.output:
        with options.output.open("xb") as out:
            out.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != data:
        raise ValueError("Seven-root archive replay failed")
