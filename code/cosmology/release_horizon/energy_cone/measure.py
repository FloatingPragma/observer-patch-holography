"""Replay the material energy cone, joint thermal words and the held-state work response."""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from math import log, prod, sqrt

import numpy as np

if __package__:
    from .instruments import Bench, HERE, POINTS, PROBES
else:
    from instruments import Bench, HERE, POINTS, PROBES
from background_clock.material import stream
from background_clock.measure import scale, thermal
from radial_transport.carrier_routes import AXES, coords
from horizon import trigger_histories, upper
from radial_transport.thermal_span import distance
from relay import apply, transpose

CONFIG = json.loads((HERE / "protocol.json").read_text())
N = 2**50


def moments(words):
    total = sum(words)
    out = []
    for which in range(3):
        val = 0
        for k, n in enumerate(words):
            a, b = (k & 1) - ((k >> 1) & 1), ((k >> 2) & 1) - ((k >> 3) & 1)
            val += n * (a, b, a * b)[which]
        out.append(F(val, total))
    return out


def covariance(words):
    a, b, ab = moments(words)
    return ab - a * b


def product_words(p):
    return [
        prod(q if k & (1 << i) else 1 - q for i, q in enumerate(p)) for k in range(16)
    ]


def heat(net, reset, emitters, seed):
    weights = [F(0)] * 16
    source_bins = []
    strata = []
    origin_control = []
    for launch, weight in trigger_histories(reset):
        trace = stream(net.levels, net.board.cells, reset, launch)
        offers = []
        for owner in emitters:
            bright = np.flatnonzero(trace[:, owner])
            offers.append(int(bright[0]) if len(bright) else reset)
        conditional = [upper(n) for n in offers]
        law = product_words(conditional)
        weights = [a + weight * q for a, q in zip(weights, law, strict=True)]
        if launch in (1, reset + 1):
            arrival = [np.flatnonzero(trace[:, owner]).tolist() for owner in emitters]
            origin_control.append(
                dict(
                    root_launch=launch if launch <= reset else None,
                    first_crest=[a[0] if a else None for a in arrival],
                    contrast_means=list(map(str, moments(law)[:2])),
                )
            )
        source_bins.append((launch, weight, trace))
        strata.append([launch, str(weight), offers])
    # The independent counting path prepares all 256 probe/reservoir words
    # and draws each allowed native exchange, then reads the four probes.
    counts = thermal(source_bins, emitters, None, N, seed, "layer")
    sampled = counts["energy_words"]
    assert sum(weights) == 1 and sum(sampled) == N
    assert F(counts["covariance"]) == covariance(weights)
    left = [sum(weights[k] for k in range(16) if k % 4 == i) for i in range(4)]
    right = [sum(weights[k] for k in range(16) if k // 4 == i) for i in range(4)]
    split = [left[k % 4] * right[k // 4] for k in range(16)]
    assert covariance(split) == 0
    observed = covariance(sampled)
    exact = covariance(weights)
    # Six bounded moments in [-1,1], both windows together; union-bound Hoeffding.
    eps = sqrt(2 * log(12 / 0.01) / N)
    a, b, _ = moments(sampled)
    bound = eps * (1 + abs(float(a)) + abs(float(b))) + eps**2
    assert abs(float(observed - exact)) < bound
    return dict(
        trigger_strata=strata,
        root_interventions=origin_control,
        exact_energy_words=list(map(str, weights)),
        sampled_energy_words=sampled,
        exact_covariance=str(exact),
        observed_covariance=float(observed),
        simultaneous_99_halfwidth=bound,
        independent_trigger_words=list(map(str, split)),
    )


def trials(net, axes, live):
    # A local unitary transfers the held thermal dimer off its old ports at
    # isolation. It moves the unobserved controller to the vacated address.
    handoff = transpose(2, 4) + transpose(3, 5)
    labels = list(range(12))
    apply(labels, handoff)
    assert labels == [0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 10, 11]
    initial = np.array([i % 2 for i in POINTS], dtype=int)
    for i, bit in enumerate(initial):
        local = net.word[12 * i : 12 * i + 12]
        local[4 + bit] = 2
        apply(local, handoff)
        net.word[12 * i : 12 * i + 12] = local
    identity = np.eye(len(POINTS), dtype=np.int8)
    bank = identity.copy()
    for owner in POINTS:
        net.word[12 * owner + 8] = 2
    counts = np.zeros(len(POINTS), dtype=int)
    for t in range(1, 26 - len(axes)):
        for _ in range(scale(2, t)):
            counts += net.marker_stroke()
    marker = [counts[list(PROBES)].tolist()]
    publications = []
    interior = [s for s in POINTS if min(coords(s)) > 0]
    clocks = []
    for axis in axes:
        ticks, delta = net.frame(axis, bank, live)
        clocks.append(ticks)
        counts += delta
        publications.append(int(sum(delta[interior])))
        marker.append(counts[list(PROBES)].tolist())
        # Fresh energy interventions at this frame boundary are separate
        # columns from the unmodified states propagated from earlier boundaries.
        bank = np.hstack((bank, identity))
    by_boundary = bank.reshape(len(POINTS), -1, len(POINTS))
    cumulative = np.bitwise_or.reduce(by_boundary, axis=1)
    transfer = np.argmax(by_boundary[:, 0, :], axis=1)
    assert np.all(by_boundary.sum(axis=0) == 1) and np.all(by_boundary.sum(axis=2) == 1)
    for i, origin in enumerate(transfer):
        assert net.read(i, 2) == 1 - initial[origin]
        net.word[12 * i + 2 : 12 * i + 4] = [0, 0]
        net.word[12 * i + 8 : 12 * i + 10] = [0, 0]
    return (
        transfer[list(PROBES)].tolist(),
        [np.flatnonzero(cumulative[i]).tolist() for i in PROBES],
        clocks,
        marker,
        publications,
    )


def cone(net):
    cross = [distance(net, a, b) for a in PROBES[:2] for b in PROBES[2:]]
    incomplete, leak = net.energy_exchange(0, 25, cut=True)
    result = []
    for reset, axes in [(24, [1]), (23, [0, 1])]:
        orig, inputs, clocks, marker, publications = trials(net, axes, True)
        stopped, quiet, park_time, park_marks, park_counts = trials(net, axes, False)
        assert park_marks == marker and publications == park_counts
        assert park_time == clocks and all(t == 100 * net.slot for t in clocks)
        assert stopped == list(PROBES) and quiet == [[i] for i in PROBES]
        ranges = [
            [distance(net, start, target) for start in group]
            for group, target in zip(inputs, PROBES, strict=True)
        ]
        radii = [max(r["distance"] for r in group) for group in ranges]
        margins = [
            r["distance"] - radii[i] - radii[2 + j]
            for r, (i, j) in zip(cross, product(range(2), repeat=2), strict=True)
        ]
        result.append(
            dict(
                reset_layer=reset,
                read_layer=25,
                clock_strokes=clocks,
                marker_readings=marker,
                publications=publications,
                read_origins=orig,
                parked_origins=stopped,
                impulse_ranges=ranges,
                endpoint_radii=radii,
                separation_minus_radii=margins,
                thermal=heat(net, reset, orig, CONFIG["seed"] + 24 - reset),
            )
        )
    return dict(
        probes=list(PROBES),
        cross_ranges=cross,
        box_axis_returns=[distance(net, 0, 4 * d) for d in AXES],
        windows=result,
        single_rail_leakage=[round(float(v), 12) for v in leak],
        single_rail_decoded_weights=np.round(abs(incomplete) ** 2, 12).tolist(),
    )


def work(net):
    t0 = 2 * net.slot
    rows = []
    for t in [t0 + d for d in (-16, -8, -4, -2, 0, 2, 4, 8, 16)] + [2 * t0]:
        strokes = net.volume_returns(t)
        lengths = [F(v, 2 * t0) for v in strokes]
        vol = prod(lengths)
        row = dict(return_strokes=strokes, relative_volume=str(vol), cycle_ticks=t)
        for name, exposure in (
            ("fixed_action", np.pi / 3),
            ("fixed_duty", np.pi * t / (3 * t0)),
        ):
            prob, gap, occupations = net.spectroscopy(exposure, t)
            row[name] = dict(
                quadratures=prob,
                gap_in_reference_units=round(3 * t0 * gap, 12),
                energy_outcomes=np.round(occupations, 12).tolist(),
            )
        rows.append(row)
    secants = []
    for i in range(4):
        left, right = rows[i], rows[8 - i]
        dv = float(F(right["relative_volume"]) - F(left["relative_volume"]))
        slope = [
            (
                left[name]["gap_in_reference_units"]
                - right[name]["gap_in_reference_units"]
            )
            / dv
            for name in ("fixed_action", "fixed_duty")
        ]
        secants.append(
            dict(half_offset=t0 - left["cycle_ticks"], upper_pressures=slope)
        )
    compression = {}
    for name in ("fixed_action", "fixed_duty"):
        change = (
            rows[4][name]["gap_in_reference_units"]
            - rows[-1][name]["gap_in_reference_units"]
        )
        compression[name] = [round(pop * change, 12) for pop in (0, 0.2, 1)]
    return dict(
        reference_ticks=t0,
        rows=rows,
        secants=secants,
        compression_from_eight_to_one=compression,
    )


def run():
    net = Bench()
    return dict(
        rail_hop_ticks=net.slot,
        clock_owner=net.clock_owner,
        energy_cone=cone(net),
        work=work(net),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), separators=(",", ":")) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(data)
    elif path.read_text() != data:
        raise SystemExit("Energy-cone archive mismatch")
    print(f"Verified {len(data)} bytes of energy and work records")
