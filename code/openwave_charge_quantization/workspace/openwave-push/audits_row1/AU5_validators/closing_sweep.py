#!/usr/bin/env python3
"""The closing contrast sweep: an appended instrument on the native source.

After a source history the observer holds the known preparation G and a retained tape of current records.
Its reconstruction of the loads is G - B F (currents only).  The sweep appends one native event on every seam
of the carrier, in a stated order with stated coins, and retains the full record of each event:

    (seam, coin, contrast, current)      contrast d = n_a - n_b read by the event, current j = j(d, coin)

Each received contrast is compared with the contrast of the reconstruction on the same seam; the reconstruction
is then advanced by the received current.  No load of the source is read by the comparison.

Callable interface (frozen):

    closing_sweep(G, retained_tape, source_state, order, coins) -> dict

        G              twelve integers, the known preparation
        retained_tape  the retained source records: (seam, current), (seam, coin, current) or
                       {'seam': .., 'current': ..}; seam is an index 0..29 or a pair (a, b)
        source_state   twelve integers: the loads of the source when the sweep starts (the events act on it)
        order          the thirty seams of the sweep, each exactly once (indices or pairs)
        coins          thirty coins, one per sweep event

    returns {'records', 'comparison', 'source_before', 'source_after', 'reconstruction_before',
             'reconstruction_after', 'reconstruction_from', 'mismatch_positions', 'first_mismatch',
             'non_native_records', 'swept_graph_connected_and_spanning', 'sweep_moving_events',
             'source_changed_by_sweep', 'verdict', 'latch', 'audit_only'}

Verdicts:
    'consistent'                      no received contrast differs from the reconstruction: under the hypotheses
                                      of PROOF.md the reconstruction equals the source at latch L0 and at L1
    'reconstruction_error_detected'   at least one received contrast differs
    'sweep_record_not_native'         a received record violates current = j(contrast, coin)
    'inconclusive_sweep_not_spanning' only with allow_partial=True: the swept seams do not connect all ports

Latches.  L0 is the end of the retained tape, before the first sweep event; a regional charge decoded from the
retained counters refers to L0.  L1 is after the thirtieth sweep event.  The sweep changes the source: the
loads at L1 are the loads at L0 advanced by the thirty sweep currents, and counters have to be advanced by the
same thirty currents before they describe L1.

Standard library only; plain integers.  Nothing is read from disk.
"""
from __future__ import annotations

__all__ = ["closing_sweep", "sweep_step", "native_event", "reconstruct", "tape_per_seam", "seeded_order",
           "order_family", "connected_and_spanning", "SEAMS", "PORTS"]

SEAMS = ((0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 6),
         (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9), (4, 11), (5, 9), (5, 11), (6, 7),
         (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11))
SEAM_INDEX = {s: i for i, s in enumerate(SEAMS)}
PORTS = 12


def current_of(contrast, coin):
    """The native current as a function of the contrast and the coin."""
    return (contrast - (2 * coin - 1) * (abs(contrast) % 2)) // 2


def native_event(loads, seam, coin):
    """One native event on `seam`.  Returns (loads after, contrast read, current)."""
    a, b = seam
    d = loads[a] - loads[b]
    j = current_of(d, coin)
    out = list(loads)
    out[a] -= j
    out[b] += j
    return tuple(out), d, j


def _seam(x):
    if type(x) is int:
        if not 0 <= x < len(SEAMS):
            raise ValueError("seam index outside 0..29")
        return SEAMS[x]
    s = tuple(x)
    if s not in SEAM_INDEX:
        raise ValueError(f"{x!r} is not a seam of the carrier")
    return s


def _seam_and_current(record):
    if isinstance(record, dict):
        return _seam(record["seam"]), record["current"]
    record = tuple(record)
    if len(record) == 2:
        return _seam(record[0]), record[1]
    if len(record) == 3:
        return _seam(record[0]), record[2]
    raise ValueError("a retained record is (seam, current), (seam, coin, current) or a dictionary")


def tape_per_seam(retained_tape):
    """F: the sum of the retained currents on each of the thirty seams."""
    tape = [0] * len(SEAMS)
    for record in retained_tape:
        seam, j = _seam_and_current(record)
        if type(j) is not int:
            raise ValueError("a retained current is not an integer")
        tape[SEAM_INDEX[seam]] += j
    return tape


def reconstruct(G, retained_tape):
    """G - B F from the retained currents only."""
    loads = list(G)
    for (a, b), f in zip(SEAMS, tape_per_seam(retained_tape)):
        loads[a] -= f
        loads[b] += f
    return tuple(loads)


def connected_and_spanning(seams):
    """True when the graph of the given seams connects all twelve ports."""
    parent = list(range(PORTS))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for a, b in seams:
        parent[find(a)] = find(b)
    return len({find(v) for v in range(PORTS)}) == 1


def seeded_order(seed):
    """A permutation of the thirty seam indices from a linear congruential generator (Fisher-Yates), fixed by
    `seed` and independent of the interpreter's random module."""
    x = (seed * 2654435761 + 1013904223) % (1 << 32)
    order = list(range(len(SEAMS)))
    for i in range(len(order) - 1, 0, -1):
        x = (1664525 * x + 1013904223) % (1 << 32)
        k = (x >> 8) % (i + 1)
        order[i], order[k] = order[k], order[i]
    return order


def order_family():
    """The fixed family of twelve sweep orders: identity, reverse, and the seeded orders 1..10."""
    fam = {"identity": list(range(len(SEAMS))), "reverse": list(range(len(SEAMS) - 1, -1, -1))}
    for seed in range(1, 11):
        fam[f"seed{seed}"] = seeded_order(seed)
    return fam


def sweep_step(source, reconstruction, seam, coin, position=0, channel=None):
    """One sweep event.  The native event acts on the source and emits its record; the record passes through
    `channel` (identity when None); the received contrast is compared with the reconstruction's contrast; the
    reconstruction is advanced by the received current.

    returns (received record, source after, reconstruction after, mismatch, native_ok,
             contrast of the reconstruction on this seam before the event)"""
    after, d, j = native_event(source, seam, coin)
    record = {"position": position, "seam": [seam[0], seam[1]], "coin": coin, "contrast": d, "current": j}
    if channel is not None:
        record = channel(dict(record), tuple(reconstruction))
    a, b = seam
    predicted = reconstruction[a] - reconstruction[b]
    mismatch = record["contrast"] != predicted
    native_ok = (type(record["contrast"]) is int and record["coin"] in (0, 1)
                 and record["current"] == current_of(record["contrast"], record["coin"]))
    rec = list(reconstruction)
    if type(record["current"]) is int:
        rec[a] -= record["current"]
        rec[b] += record["current"]
    return record, after, tuple(rec), mismatch, native_ok, predicted


def closing_sweep(G, retained_tape, source_state, order, coins, *, channel=None, allow_partial=False):
    """Perform the appended sweep on the source and compare its contrast records with the reconstruction.
    See the module text for the interface, the verdicts and the two latches."""
    G = tuple(G)
    source_before = tuple(source_state)
    if len(G) != PORTS or len(source_before) != PORTS or any(type(v) is not int for v in G + source_before):
        raise ValueError("G and the source state are twelve integers")
    seams = [_seam(x) for x in order]
    coins = list(coins)
    if len(coins) != len(seams) or any(type(c) is not int or c not in (0, 1) for c in coins):
        raise ValueError("one native coin for every sweep event")
    full = sorted(seams) == sorted(SEAMS)
    if not full and not allow_partial:
        raise ValueError("the sweep order must contain every one of the thirty seams exactly once")
    spanning = connected_and_spanning(seams)

    reconstruction_before = reconstruct(G, retained_tape)
    source = source_before
    reconstruction = reconstruction_before
    records = []
    comparison = []
    mismatches = []
    non_native = []
    moving = 0
    for k, (seam, coin) in enumerate(zip(seams, coins)):
        record, source, reconstruction, mismatch, native_ok, predicted = sweep_step(
            source, reconstruction, seam, coin, k, channel)
        records.append(record)
        comparison.append({"position": k, "seam": [seam[0], seam[1]], "received_contrast": record["contrast"],
                           "reconstructed_contrast": predicted, "match": not mismatch})
        if mismatch:
            mismatches.append(k)
        if not native_ok:
            non_native.append(k)
        if record["current"]:
            moving += 1
    if non_native:
        verdict = "sweep_record_not_native"
    elif mismatches:
        verdict = "reconstruction_error_detected"
    elif not spanning:
        verdict = "inconclusive_sweep_not_spanning"
    else:
        verdict = "consistent"
    return {
        "schema": "oph.push.closing_contrast_sweep.v1",
        "records": records,
        "comparison": comparison,
        "source_before": list(source_before),
        "source_after": list(source),
        "reconstruction_before": list(reconstruction_before),
        "reconstruction_after": list(reconstruction),
        "reconstruction_from": {"G": list(G), "retained_records": len(list(retained_tape)),
                                "retained_tape_per_seam": tape_per_seam(retained_tape)},
        "mismatch_positions": mismatches,
        "first_mismatch": mismatches[0] if mismatches else None,
        "non_native_records": non_native,
        "swept_seams": len(seams),
        "swept_graph_connected_and_spanning": spanning,
        "sweep_moving_events": moving,
        "source_changed_by_sweep": tuple(source) != source_before,
        "verdict": verdict,
        "latch": {
            "L0": "end of the retained tape, before the first sweep event: the reconstruction compared against "
                  "(reconstruction_before) and any charge decoded from the retained counters refer to L0",
            "L1": "after the last sweep event: the loads are those of L0 advanced by the sweep currents "
                  "(reconstruction_after when the verdict is 'consistent'); counters must be advanced by the "
                  "sweep currents before a readout at L1",
            "charge_read_at": "L0",
        },
        "audit_only": {"reconstruction_equals_source_at_L0": reconstruction_before == source_before,
                       "reconstruction_equals_source_at_L1": tuple(reconstruction) == tuple(source),
                       "totals_equal": sum(G) == sum(source_before)},
    }
