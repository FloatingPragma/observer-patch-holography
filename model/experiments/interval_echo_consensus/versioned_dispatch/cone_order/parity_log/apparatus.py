"""Run all 37 appointments without resetting the roots or earlier outputs."""

from functools import cache

from model.executable_model import carrier
from ....seam_word_meter.finite_program.layout.mount import Mount
from ....seam_word_meter.finite_program.machine import decode_word

ROOTS = 5
MASKS = 1 << ROOTS
SHELVES = ROOTS + MASKS
WIDTH = 7


def tape_word(raw):
    """One hex digit retains both two-reference pointer codes on a balanced dimer."""
    return "".join(format(4 * a + b, "x") for a, b in raw)


def enroll():
    board = Mount(carrier.federation(1)).bench(registers=SHELVES, axes=1, width=WIDTH)
    board.route = cache(board.route)
    return board


def layout(board):
    return {
        "width": board.width,
        "cells": board.graph.cells,
        "capacity": board.graph.capacity,
        "heads": board.heads,
        "references": board.references,
        "words": list(board.words),
    }


def read(board, address, station, parked=False):
    before = board.loads.copy()
    board.sensor = "parked" if parked else "native"
    value = board.read(address, station)
    raw = board.tape.pop()["raw"]
    board.sensor = "native"
    if board.loads != before:
        raise ArithmeticError("Reading the shelf changed material stock")
    return value, tape_word(raw)


def store(board, address, value):
    """Seal a fresh target. Check the physical frame across the entire appointment."""
    start = board.loads.copy()
    station = address % 2
    board.put(address, value, station)
    command = board.tape.pop()
    if decode_word(command["before"]) != 0:
        raise ArithmeticError("The event target was already published")
    received = [read(board, address, at) for at in (0, 1)]
    if any(v != value for v, _ in received):
        raise ArithmeticError("Peer failed to receive the sealed value and writer")
    footprint = {p for pair in board.words[address] for p in pair}
    changed = [
        i for i, (a, b) in enumerate(zip(start, board.loads, strict=True)) if a != b
    ]
    if not set(changed) <= footprint:
        raise ArithmeticError("Sealed appointment changed another material port")
    return {
        "at": station,
        "before": tape_word(command["before"]),
        "after": tape_word(command["after"]),
        "peers": [raw for _, raw in received],
        "changed_ports": changed,
    }


def run(board, word, mode):
    board.start([0] * SHELVES)
    events = []
    for address in range(ROOTS):
        bit = (word >> address) & 1
        publication = store(board, address, ((address + 1) << 1) | bit)
        events.append(
            dict(id=address + 1, target=address, reads=[], publication=publication)
        )

    for mask in range(MASKS):
        target = ROOTS + mask
        identifier = target + 1
        actual = 3 if mode == "altered" and mask == 5 else mask
        wrong_stamp = mode == "stamp" and mask == 31
        if wrong_stamp:
            board.swaps(board.loads, [board.words[0][1]])

        reads = []
        parity = 0
        valid = True
        for source in range(ROOTS):
            if not (actual & (1 << source)):
                continue
            left, first = read(board, source, 0)
            shutter = mode == "shutter" and mask == 31 and source == 0
            right, second = read(board, source, 1, parked=shutter)
            reads.append([source, first, second])
            if left is None or left != right or left >> 1 != source + 1:
                valid = False
            else:
                parity ^= left & 1

        publication = (
            store(board, target, (identifier << 1) | parity) if valid else None
        )
        events.append(
            dict(id=identifier, target=target, reads=reads, publication=publication)
        )
        if wrong_stamp:
            board.swaps(board.loads, [board.words[0][1]])

    final = [
        [read(board, address, at)[1] for address in range(SHELVES)] for at in (0, 1)
    ]
    if final[0] != final[1] or board.tape:
        raise ArithmeticError("Final material census or tape disagrees")
    return dict(input=word, mode=mode, events=events, final=final)
