"""Register two earlier heat contrasts to distances in a separate ruler calibration."""

from fractions import Fraction
from hashlib import sha256
import json

import numpy as np

if __package__:
    from .carrier_routes import AXES, CONFIG, HERE, coords
else:
    from carrier_routes import AXES, CONFIG, HERE, coords


def input_challenges(net, probes, sweeps):
    # One sharp intervention in each column, all subsequent root writes parked.
    state = np.eye(net.board.cells, dtype=np.int8)
    seen = state.astype(bool)
    for _ in range(sweeps):
        for parents, children, tables in net.levels:
            state[children] = tables[
                np.arange(len(children))[:, None], 1, state[parents], state[children]
            ]
        state[0] = 0
        seen |= state.astype(bool)
    return [np.flatnonzero(seen[v]).tolist() for v in probes]


def distance(net, start, target):
    packet = CONFIG["packet_port"]
    beacon = CONFIG["range_port"]
    assert not net.read(start, packet) and not net.read(target, beacon)
    net.word[12 * start + packet] = net.word[12 * target + beacon] = 2
    begin = net.service_clock
    visited = [start]
    for x, y, stride in zip(coords(start), coords(target), AXES, strict=True):
        step = stride if y > x else -stride
        for _ in range(abs(y - x)):
            a = visited[-1]
            assert not net.read(a, beacon)
            net.hop(a, a + step, packet, metered=True)
            assert net.read(a + step, packet)
            visited.append(a + step)
    assert net.read(visited[-1], beacon)
    for a, b in zip(visited[:0:-1], visited[-2::-1], strict=True):
        net.hop(a, b, packet, metered=True)
    assert net.read(start, packet)
    net.word[12 * start + packet] = net.word[12 * target + beacon] = 0
    strokes = net.service_clock - begin
    assert strokes % (2 * net.slot) == 0
    return dict(
        start=start,
        target=target,
        return_strokes=strokes,
        distance=strokes // (2 * net.slot),
    )


def measure(net):
    original = (HERE.parent / "background_clock" / "records.json").read_bytes()
    heat = json.loads(original)["boards"][-1]["thermal"]
    probes = heat["cells"]
    assert len(set(probes)) == 4
    inputs = input_challenges(net, probes, heat["deadline"] - heat["reset_layer"])
    assert inputs == heat["pulse_inputs_per_probe"]
    assert all(0 <= v < 125 for group in inputs for v in group)
    incoming = [
        [distance(net, v, target) for v in origins]
        for origins, target in zip(inputs, probes, strict=True)
    ]
    cross = [distance(net, a, b) for a in probes[:2] for b in probes[2:]]
    radii = [max(item["distance"] for item in bank) for bank in incoming]
    shared = sorted(set().union(*inputs[:2]) & set().union(*inputs[2:]))
    words = heat["layer"]["energy_words"]
    means = [Fraction(0)] * 3
    for code, count in enumerate(words):
        a = (code & 1) - ((code >> 1) & 1)
        b = ((code >> 2) & 1) - ((code >> 3) & 1)
        for i, value in enumerate((a, b, a * b)):
            means[i] += Fraction(count * value, sum(words))
    observed = means[2] - means[0] * means[1]
    assert abs(float(observed) - heat["layer"]["observed"]) < 1e-14
    return dict(
        energy_archive_sha256=sha256(original).hexdigest(),
        probes=probes,
        reset_marker=heat["reset_marker"],
        read_marker=heat["marker"],
        input_challenges=net.board.cells,
        shared_inputs=shared,
        input_ranges=incoming,
        endpoint_radii=radii,
        cross_ranges=cross,
        energy_words=words,
        observed_covariance=str(observed),
        predicted_covariance=heat["layer"]["covariance"],
        simultaneous_99pct_halfwidth=heat["layer"]["simultaneous_99pct_halfwidth"],
    )
