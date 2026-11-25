"""Count material crossings, vertex balances, return times and proper marker changes."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json

import numpy as np

if __package__:
    from .carrier_routes import (
        AXES,
        CONFIG,
        HERE,
        POINTS,
        SIZE,
        Transport,
        coords,
        grid_edges,
    )
else:
    from carrier_routes import (
        AXES,
        CONFIG,
        HERE,
        POINTS,
        SIZE,
        Transport,
        coords,
        grid_edges,
    )
from background_clock.measure import census, scale, signal_readings
from background_clock.material import stream

D = SIZE - 1


def ray_steps(axis, i, j):
    slopes = [0, 0, 0]
    slopes[axis] = 24
    for k, rank in zip((k for k in range(3) if k != axis), (i, j), strict=True):
        slopes[k] = 2 * rank + 1
    crossings = sorted(
        (Fraction(24 * m, speed), k)
        for k, speed in enumerate(slopes)
        for m in range(1, D + 1)
        if 24 * m <= D * speed
    )
    return [AXES[k] for _, k in crossings]


def layer(site):
    return 0 if site == -1 else 1 + max(coords(site))


def shell_edges(edges, n):
    return [(a, b) for a, b in edges if layer(a) <= n < layer(b)]


def probe_incidence(net, off, src, edges, cut=None):
    """Read both endpoints of every offered link with fresh single-load challenges."""
    port = CONFIG["packet_port"]
    observed = []
    for a, b in edges:
        sender, receiver = src if a == -1 else off + a, off + b
        assert not net.read(sender, port) and not net.read(receiver, port)
        net.word[12 * sender + port] = 2
        net.hop(sender, receiver, port, cut=(a, b) == cut)
        at_sender, at_receiver = net.read(sender, port), net.read(receiver, port)
        assert at_sender + at_receiver == 1
        if at_receiver:
            observed.append((a, b))
        net.word[12 * (receiver if at_receiver else sender) + port] = 0
    return observed


def current(net, off, src, edges, mode="cycle"):
    port = CONFIG["packet_port"]
    hits = Counter()
    paired = Counter()
    edge_index = {e: k for k, e in enumerate(edges)}
    sinks = Counter()
    starts = 0
    for axis, i, j in product(range(3), range(12), range(12)):
        if mode == "single_face" and axis:
            continue
        starts += 1
        assert net.read(src, port) == 0
        net.word[12 * src + port] = (
            2  # one freshly prepared load at the injection owner
        )
        owner, site = src, -1
        crossings = {}
        for step in [None, *ray_steps(axis, i, j)]:
            nxt = 0 if step is None else site + step
            receiver = off + nxt
            assert net.read(owner, port) and not net.read(receiver, port)
            net.hop(owner, receiver, port)
            assert not net.read(owner, port) and net.read(receiver, port)
            if layer(site) < layer(nxt):
                crossings[layer(site)] = edge_index[site, nxt]
            hits[site, nxt] += 1
            site, owner = nxt, receiver
            if mode == "absorb" and site == AXES[0]:
                break
        if 2 in crossings and 4 in crossings:
            paired[crossings[2], crossings[4]] += 1
        sinks[site] += net.read(owner, port)
        net.word[12 * owner + port] = (
            0  # terminal read and reset, including the control drain
        )
    assert not hits.keys() - set(edges)
    balance = Counter()
    for (a, b), count in hits.items():
        balance[a] += count
        balance[b] -= count
    assert balance[-1] == starts and sum(sinks.values()) == starts
    assert all(balance[v] == -sinks[v] for v in POINTS)
    shells = []
    for n in range(1, D + 1):
        cut = shell_edges(edges, n)
        hist = Counter(hits[e] for e in cut)
        shells.append(
            dict(
                layer=n,
                boundary=len(cut),
                hits_by_multiplicity=sorted(hist.items()),
                total=sum(hits[e] for e in cut),
            )
        )
    return dict(
        injections=starts,
        shells=shells,
        paired_shells_2_4=[[a, b, h] for (a, b), h in sorted(paired.items())],
        edge_hits=[hits[e] for e in edges],
        vertex_outflow=[balance[v] for v in [-1, *POINTS]],
        terminal_reads=sorted(sinks.items()),
    )


def park_or_move_reflectors(net, off, before, after):
    port = CONFIG["range_port"]
    for axis, step in enumerate(AXES):
        if before is None:
            net.word[12 * (off + step) + port] = 2
            start = 1
        else:
            start = before[axis]
        stop = after[axis]
        direction = 1 if stop >= start else -1
        for k in range(start, stop, direction):
            a, b = off + k * step, off + (k + direction) * step
            assert net.read(a, port) and not net.read(b, port)
            net.hop(a, b, port)
        assert net.read(off + stop * step, port)


def range_read(net, off, axis, cut_at=None):
    port = CONFIG["packet_port"]
    step = AXES[axis]
    begin = net.service_clock
    assert not net.read(off, port)
    net.word[12 * off + port] = 2
    visited = [off]
    result = None
    for k in range(1, SIZE):
        a, b = off + (k - 1) * step, off + k * step
        net.hop(a, b, port, cut=k == cut_at, metered=True)
        if not net.read(b, port):
            assert net.read(a, port)
            break
        visited.append(b)
        if net.read(b, CONFIG["range_port"]):
            result = k
            break
    for a, b in zip(visited[:0:-1], visited[-2::-1], strict=True):
        net.hop(a, b, port, metered=True)
    assert net.read(off, port)
    net.word[12 * off + port] = 0
    elapsed = net.service_clock - begin
    if result is not None:
        assert elapsed == 2 * result * net.slot
    return dict(
        reflection=result is not None,
        elapsed_slots=elapsed,
        width=elapsed // (2 * net.slot) if result is not None else None,
    )


def marker_stroke(net, excluded=None):
    p = CONFIG["marker_port"]
    before = [net.read(owner, p) for owner in range(net.board.cells)]
    for owner in range(net.board.cells):
        if owner == excluded:
            continue
        a = 12 * owner + p
        net.word[a], net.word[a + 1] = net.word[a + 1], net.word[a]
    after = [net.read(owner, p) for owner in range(net.board.cells)]
    return np.not_equal(before, after).astype(int)


def publications(delta, off, widths):
    return [
        int(delta[off + s])
        for s in POINTS
        if all(0 < x <= hi for x, hi in zip(coords(s), widths, strict=True))
    ]


def history(net):
    horizon = 36
    offsets = CONFIG["box_offsets"]
    qs = [scale(2, t) for t in range(1, horizon + 1)]
    reference = stream(net.levels, net.board.cells, horizon, 0)
    pulses = np.zeros_like(reference)
    for owner in range(net.board.cells):
        net.word[12 * owner : 12 * owner + 2] = [0, 2] if owner == 0 else [2, 0]
    pulses[0] = np.array(net.word[1::12]) // 2
    for owner in range(net.board.cells):
        net.word[12 * owner + 8] = 2
    elapsed = [np.zeros(net.board.cells, dtype=int)]
    boxes = [dict(ranges=[], publications=[]) for _ in offsets]
    widths = [None for _ in offsets]
    for t, q in enumerate(qs, 1):
        # The parent-copy truth table uses the crest pair and clean ports 10/11.
        # Its transaction restores all other payloads. Routes finish between sweeps.
        for parents, children, tables in net.levels:
            bits = np.array(net.word[1::12]) // 2
            incoming = tables[
                np.arange(len(children)), 1, bits[parents], bits[children]
            ]
            for owner, bright in zip(children, incoming, strict=True):
                net.word[12 * owner : 12 * owner + 2] = [0, 2] if bright else [2, 0]
        net.word[:2] = [0, 2] if t in (2, 4, 6) else [2, 0]
        pulses[t] = np.array(net.word[1::12]) // 2
        if t == 1 or q != qs[t - 2]:
            for k, off in enumerate(offsets):
                park_or_move_reflectors(net, off, widths[k], [q] * 3)
                clocks_before = net.word[8::12]
                scans = [range_read(net, off, axis) for axis in range(3)]
                widths[k] = [s["width"] for s in scans]
                boxes[k]["ranges"].append(dict(layer=t, axes=scans))
                assert clocks_before == net.word[8::12]
        diff = sum(marker_stroke(net) for _ in range(q))
        elapsed.append(elapsed[-1] + diff)
        for box, off, w in zip(boxes, offsets, widths, strict=True):
            per_site = publications(diff, off, w)
            box["publications"].append(
                dict(layer=t, sites=len(per_site), changes=sum(per_site))
            )
        assert net.word[1::12] == (2 * pulses[t]).tolist()
    np.testing.assert_array_equal(pulses, reference)
    elapsed = np.stack(elapsed)
    independent, angular = census(2, qs, net.flip)
    np.testing.assert_array_equal(elapsed, independent)
    receivers = [int(np.flatnonzero(net.radii == r)[0]) for r in (8, 16)]
    signals = signal_readings(pulses, receivers, elapsed)
    controls = []
    for off in offsets:
        cut = range_read(net, off, 2, cut_at=4)
        muted = sum(marker_stroke(net, excluded=off + sum(AXES)) for _ in range(4))
        park_or_move_reflectors(net, off, [4] * 3, [2, 3, 4])
        scans = [range_read(net, off, axis) for axis in range(3)]
        delta = sum(marker_stroke(net) for _ in range(4))
        controls.append(
            dict(
                cut_range=cut,
                parked_changes=sum(publications(muted, off, [4] * 3)),
                rectangular_axes=scans,
                rectangular_changes=sum(
                    publications(delta, off, [s["width"] for s in scans])
                ),
            )
        )
    return dict(
        boxes=boxes,
        pulses=signals,
        receiver_radii=[8, 16],
        angular_control_per_face=[row[0] for row in angular],
        controls=controls,
    )


def run():
    if __package__:
        from .thermal_span import measure as thermal_span
    else:
        from thermal_span import measure as thermal_span
    net = Transport()
    spans = thermal_span(net)
    timed = history(net)
    protected = net.word.copy()
    boxes = []
    candidates = [(-1, 0), *grid_edges()]
    cut = shell_edges(candidates, 3)[0]
    cut_counts = []
    for off, src in zip(CONFIG["box_offsets"], CONFIG["injection_owners"], strict=True):
        edges = probe_incidence(net, off, src, candidates)
        assert edges == candidates
        cut_counts.append(
            len(shell_edges(probe_incidence(net, off, src, candidates, cut), 3))
        )
        boxes.append(
            {
                mode: current(net, off, src, edges, mode)
                for mode in ("cycle", "single_face", "absorb")
            }
        )
        assert net.word == protected
    return dict(
        carrier_owners=net.board.cells,
        slot_seam_exposures=net.slot,
        range_clock_owner=net.clock_owner,
        thermal_ranging=spans,
        edge_order=edges,
        boxes=boxes,
        history=timed,
        removed_edge_control=dict(edge=cut, shell=3, measured=cut_counts),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    serialized = json.dumps(run(), separators=(",", ":")) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(serialized)
    elif path.read_text() != serialized:
        raise SystemExit("archived records differ")
    print(f"Verified {len(serialized)} bytes of radial transport records")


if __name__ == "__main__":
    main()
