"""One-link radar and fixed-home reunion of the registered Interval Echo hand."""
from __future__ import annotations

from collections import deque
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
SCHEDULES = ((3, 5), (5, 3), (7, 7), (8, 8))
FIELDS = ("emit", "calibration_receipt", "start", "far_report", "reunion", "traveller", "raw_primary_control")
BLOCKS = 64
N = 2**21


def paths(board, far, width):
    neighbors = {i: sorted({p // 12 for q in range(12 * i, 12 * i + 12)
                            for p in board.neighbours[q]} - {i}) for i in range(board.cells)}

    def chains(start, length, exclude):
        def extend(path):
            if len(path) == length:
                yield path
            else:
                for c in neighbors[path[-1]]:
                    if c not in path and c not in exclude:
                        yield from extend(path + [c])
        return extend([start])

    for home_path in chains(0, 5, {far}):
        private = home_path[1:]
        for tail in chains(0, 9 + width, {*private, far}):
            rail = tail[::-1] + [far]
            reservoirs = sorted(set(neighbors[far]) - {*private, *rail})
            if not reservoirs:
                continue
            excluded = {*private, *rail, reservoirs[0]}
            bank = []
            queue = deque(private)
            while queue and len(bank) < 18:
                current = queue.popleft()
                for c in neighbors[current]:
                    if c not in excluded:
                        excluded.add(c)
                        bank.append(c)
                        queue.append(c)
                        if len(bank) == 18:
                            break
            if len(bank) == 18:
                return rail, private, bank
    raise ValueError("A finite guide and home storage could not be placed")


def material(board, far, pattern, seed, *, width=3, trigger_enabled=True):
    shared = 0
    rail, private, bank = paths(board, far, width)
    reg = clock.Registers(board, seed)
    home = {shared, *private, *bank}
    near_owner, far_owner = {shared}, {far}
    home_front = private[0]
    count_home = [d for c in private[1:] for d in clock._digits(c)]
    book_digits = [d for c in bank for d in clock._digits(c)]
    book = [book_digits[i:i + 9] for i in range(0, len(book_digits), 9)]
    # The outbound port is 0 and the return port is 6. The far response double
    # waits in a separate cell until an arrival is read on the far marker.
    reg.prepare_front(shared)
    reg.prepare_front(home_front)
    reg.initialize_digits(count_home + book_digits)
    for c in rail[:-2]:
        reg.initialize_digits(clock._digits(c))
    def movement(c):
        return clock._digits(c, clock.EXTRA_DIGITS[:1])[0]

    reg.initialize_digits([movement(rail[0])])

    def put_double(src, dst, support):
        if reg.loads[src] != 2 or reg.loads[dst] != 0:
            raise ValueError("unprepared optical source or destination")
        for a, b in reg.route(src, dst):
            reg.swap(a, b, support)

    def detects(cell, rail_port):
        if reg.loads[12 * cell + 11] != 0:
            raise ValueError("optical vacancy is occupied")
        return reg.flag(12 * cell + rail_port, 12 * cell + 11, {cell})

    def detects_remote_response():
        word = reg.route(remote_in, remote_out + 9)
        for a, b in word:
            reg.swap(a, b, far_owner)
        value = reg.flag(remote_out + 9, remote_out + 4, far_owner)
        for a, b in reversed(word):
            reg.swap(a, b, far_owner)
        return value

    def read_home():
        return reg.reading(home_front, count_home, home)

    def store(slot):
        reg.put(book[slot], read_home(), home)
        if reg.get(book[slot], home) != read_home():
            raise ArithmeticError("home register lost the optical reading")

    def home_service():
        positive = reg.fair(home_front, home)
        if positive:
            reg.post(home_front, count_home, home)
        return positive

    both = near_owner | far_owner
    outbound = 12 * shared
    remote_out = 12 * far
    remote_in = remote_out + 6
    near_in = outbound + 6
    reg.loads[outbound] = 2
    store(0)
    calibration = []
    calibration.append(home_service())
    put_double(outbound, remote_out, both)
    if not detects(far, 0):
        raise ArithmeticError("far mirror missing")
    put_double(remote_out, remote_in, far_owner)
    calibration.append(home_service())
    put_double(remote_in, near_in, both)
    if not reg.flag(near_in, outbound + 10, near_owner):
        raise ArithmeticError("home mirror receipt missing")
    store(1)
    # Reset the calibration double and supply a separate held-quiet remote arm.
    reg.loads[near_in] = 0
    # Use the optical response double from a private remote reservoir only on the
    # far marker read. It is protected in a distinct cell alongside the far station.
    parked = next(c for c in sorted({p // 12 for i in range(12) for p in board.neighbours[remote_out + i]})
                  if c not in {*rail, *home})
    reservoir = 12 * parked
    reg.loads[reservoir] = 2
    store(2)
    position = 0
    events = []
    returned = False
    for duration in pattern:
        for opportunity in range(1, duration + 1):
            primary_home = home_service()
            owner = set(rail[position:position + len(rail) - 1])
            prior = rail[position + len(rail) - 2]
            if not reg.flag(12 * prior + 8, 12 * prior + 9, owner):
                raise ValueError("departure marker absent")
            # The armed response is delivered at the *following* optical service.
            # In these frozen schedules that precedes the return handoff.
            response = bool(detects_remote_response()) if position == 1 else False
            if response:
                put_double(remote_in, near_in, both)
                # During reception port 6 of the original shared cell is vacant
                # beside the traveller’s live bank. The read compares it to port 7.
            if opportunity == duration:
                reg.shift(rail, right=(position == 0))
            did_move = not reg.flag(12 * prior + 8, 12 * prior + 9, owner)
            if did_move:
                position = 1 - position
            cells = rail[position:position + len(rail) - 1]
            owner = set(cells)
            front = cells[-1]
            arrival = reg.flag(12 * front + 8, 12 * front + 9, owner)
            if not arrival:
                raise ValueError("no material arrival")
            moved, = reg.signed_record(movement(cells[0]), int(did_move), owner, (owner,))
            if moved and front == far and trigger_enabled and reg.flag(remote_out + 8, remote_out + 9, far_owner):
                put_double(reservoir, remote_in, {parked, far})
            if response:
                if not reg.flag(near_in, outbound + 7, near_owner):
                    raise ArithmeticError("far report missing at home")
                store(3)
                # Return the optical double immediately to its remote reservoir,
                # leaving the sharp memory rail in its initial state.
                put_double(near_in, remote_in, both)
                put_double(remote_in, reservoir, {parked, far})
            count, fields = clock._fields(cells[:-1], width)
            primary = reg.fair(front, owner)
            offers = []

            def fair():
                offers.append(reg.fair(front, owner))
                return offers[-1]

            tick = clock._material_update(reg, owner, fields, moved, primary, width, fair)
            bits = iter(offers)
            check = clock.closing_tick(opportunity, primary, bits.__next__, age_bits=width) if moved else primary
            if tick != check or next(bits, None) is not None:
                raise ValueError("registered finite decision mismatch")
            if tick:
                reg.post(front, count, owner)
            reunion = bool(moved and front == shared and arrival)
            if reunion:
                store(4)
                returned = True
            events.append({"cell": front, "handoff": bool(moved), "report": response, "reunion": reunion,
                           "home_primary": primary_home, "traveller_primary": primary, "tick": int(tick),
                           "extra": offers, "elapsed": reg.reading(front, count, owner)})
    if not returned or sum(e["report"] for e in events) != int(trigger_enabled):
        raise ValueError("incomplete travelling packet")
    seam = (12 * shared + 8, 12 * shared + 9, 12 * shared + 6)
    public = []
    for name, data, writer in (("traveller", count, owner), ("home", count_home, home)):
        agreed = []
        for i, d in enumerate(data):
            a, b = reg.signed_record(seam, reg.read(d, writer), writer, (home, owner))
            if a != b:
                raise ValueError("owners disagree")
            agreed.append(a)
            if name == 'traveller':
                reg.assign(book[5][i], a, home)
        public.append(sum(bit << i for i, bit in enumerate(agreed)))
    words = [reg.get(b, home) for b in book]
    if public[1] != words[4] or public[0] != events[-1]["elapsed"] or public[0] != words[5]:
        raise ValueError("material reunion differs from signed register read")
    return {"carrier_seams": seam_digest(board), "cells": board.cells, "seed": seed, "guide": rail,
            "private_home": private, "home_word_banks": bank, "home_patch_cells": len(home),
            "far_response_reservoir": parked, "physical_link": next([a,b] for a,b in board.external
                if {a//12,b//12} == {shared,far}), "width": width, "intervals": list(pattern), "calibration_primaries": calibration,
            "events": events, "material_words": words, "operator_counts": dict(reg.gates)}


def estimator(means):
    emit, cal, start, echo, reunion, traveller, raw = means
    d = (cal - emit) / 2
    legs = np.array([echo - start - d, reunion - echo + d])
    if not (d > 0 and min(legs) > d):
        raise ValueError("Reciprocal optical estimate lies outside the timed domain")
    prediction = np.sqrt(legs**2 - d**2).sum() / legs.sum()
    observed, control = np.array([traveller, raw]) / legs.sum()
    return np.array([d, *(d / legs), prediction, observed, observed - prediction,
                     control, control - prediction, observed - control])


LABELS = ("distance", "out_speed", "in_speed", "optical_forecast", "ratio", "residual",
          "raw_ratio", "raw_residual", "filtered_minus_raw")


def decision(point, se):
    if abs(point) + 4 * se <= 0.004:
        return "agrees"
    if abs(point) - 4 * se > 0.004:
        return "disagrees"
    return "unresolved"


def sample_block(events, width, seed):
    from model.candidate_laws.interval_echo.cohorts import suppression

    rng = np.random.default_rng(seed)
    n = N // BLOCKS
    words = np.zeros((n, 7), dtype=np.uint8)
    home = rng.binomial(2, .5, size=n).astype(np.uint8)
    words[:, 1] = words[:, 2] = home
    traveler = np.zeros(n, dtype=np.uint8)
    raw = np.zeros(n, dtype=np.uint8)
    age = 0
    for moved, reported, reunion in events:
        home += rng.integers(2, size=n, dtype=np.uint8)
        primary = rng.integers(2, size=n, dtype=np.uint8)
        raw += primary
        age += 1
        if moved:
            tick = primary * (rng.random(n) >= suppression(age, width))
            age = 0
        else:
            tick = primary
        traveler += tick
        if reported:
            words[:, 3] = home
        if reunion:
            words[:, 4] = home
            words[:, 5] = traveler
            words[:, 6] = raw
    # Each nine-bit home register is below 32 for these journeys. Store the exact
    # joint multiplicity of its observed seven-word packets, including zero.
    shifts = np.arange(7, dtype=np.uint64) * 5
    packed = np.sum(words.astype(np.uint64) << shifts, axis=1)
    codes, counts = np.unique(packed, return_counts=True)
    return codes, counts.astype(np.uint32)


def unpack(codes):
    return (codes[:, None] >> (np.arange(7, dtype=np.uint64) * 5)) & 31


def analyze(keys, counts_by_block, offsets):
    totals = []
    for a, b in zip(offsets[:-1], offsets[1:]):
        weights = counts_by_block[a:b]
        if weights.sum() != N // BLOCKS:
            raise ValueError("Block enrollment changed")
        totals.append(weights @ unpack(keys[a:b]))
    totals = np.array(totals, dtype=float)
    whole = totals.sum(axis=0)
    point = estimator(whole / N)
    leave = np.array([estimator((whole - x) / (N - N // BLOCKS)) for x in totals])
    errors = np.sqrt((BLOCKS - 1) * np.mean((leave - leave.mean(axis=0))**2, axis=0))
    return {key: {"value": float(x), "paired_block_se": float(error)} for key,x,error in zip(LABELS,point,errors)}


def cohort(trace, width, seed):
    keys, counts, offsets = [], [], [0]
    for block in range(BLOCKS):
        k, c = sample_block(trace, width, np.random.SeedSequence([seed, block]))
        keys.append(k)
        counts.append(c)
        offsets.append(offsets[-1] + len(k))
    records = {"keys": np.concatenate(keys), "counts": np.concatenate(counts), "offsets": np.array(offsets)}
    measures = analyze(records["keys"], records["counts"], records["offsets"])
    return {"width": width, "base_seed": seed, "enrolled": N, "blocks": BLOCKS, "estimates": measures,
            "decision": decision(measures["residual"]["value"], measures["residual"]["paired_block_se"]),
            "raw_decision": decision(measures["raw_residual"]["value"], measures["raw_residual"]["paired_block_se"])}, records


def flags(row):
    return tuple(tuple(e[k] for k in ("handoff", "report", "reunion")) for e in row["events"])


def measure():
    witnesses, traces = [], {}
    for level in (1, 2):
        board = federation(level)
        targets = sorted(y // 12 for x, y in board.external if x // 12 == 0)
        if len(targets) != 3:
            raise ValueError("Expected three native home directions")
        for far in targets:
            for i, schedule in enumerate(SCHEDULES):
                seed = 9201900 + len(witnesses)
                witness = material(board, far, schedule, seed)
                witnesses.append(witness)
                trace = flags(witness)
                if traces.setdefault(i, trace) != trace:
                    raise ValueError("A different timeline requires a separate statistical cohort")
    broad = []
    for w in witnesses:
        if tuple(w["intervals"]) == (8, 8):
            far = w["guide"][-1]
            board = federation(1 if w["cells"] == 80 else 2)
            row = material(board, far, (8, 8), w["seed"] + 4019, width=5)
            if flags(row) != flags(w):
                raise ValueError("Capacity comparison changes the physical event timeline")
            broad.append(row)
    arrays, panels = {}, []
    for name, pattern, width, trace, seed in [
            *[(f"q{a}_{b}", (a,b), 3, traces[i], 551200 + i) for i,(a,b) in enumerate(SCHEDULES)],
            ("capacity8", (8,8), 5, traces[3], 996173)]:
        panel, data = cohort(trace, width, seed)
        for key, value in data.items():
            arrays[name + "_" + key] = value
        panels.append({"id": name, "intervals": list(pattern), "installations": 6, **panel})
    return {"candidate": "interval_echo", "fields": list(FIELDS), "protocol_sha256": sha256((HERE/"PROTOCOL.md").read_bytes()).hexdigest(),
            "capacity_protocol_sha256": sha256((HERE/"CAPACITY_PROTOCOL.md").read_bytes()).hexdigest(),
            "material_witnesses": witnesses, "five_bit_material_witnesses": broad, "panels": panels}, arrays


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write data to a fresh destination.")
    args = parser.parse_args()
    receipt, arrays = measure()
    text = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        paths = [args.output / "results.json", args.output / "packets.npz"]
        if any(path.exists() for path in paths):
            raise FileExistsError("Retained data already exist")
        np.savez_compressed(paths[1], **arrays)
        paths[0].write_text(text)
    else:
        if json.loads((HERE / "results.json").read_text()) != receipt:
            raise ValueError("Stored measurements differ")
        with np.load(HERE / "packets.npz", allow_pickle=False) as saved:
            if set(saved.files) != set(arrays):
                raise ValueError("Count fields differ")
            for key, data in arrays.items():
                np.testing.assert_array_equal(data, saved[key])
    for p in receipt["panels"]:
        print(p["id"],p["decision"],p["estimates"]["residual"],p["raw_decision"])
