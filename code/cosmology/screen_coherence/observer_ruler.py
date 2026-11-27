"""A same-reception-time ruler survey with bilateral source echoes and twelve clocks."""

import argparse
from collections import Counter, defaultdict
import json
from math import acos, comb
from pathlib import Path

import numpy as np

from distance_redshift import PERIOD, TIMEOUT
from observer_zero_information import local_gram
from observer_zero_radar import ArchiveRadar, MANIFEST_SHA, local_graph


LAST = 96
FIRST = 2  # The two-tick source certificate can first travel on this crest.
SEED = 20261126


def arrival_words(edges, stations, home, last, limit):
    """Earliest tagged arrivals at each local clock, for every emission tick."""
    miss = limit + 1
    n = len(edges)
    # Each slot below contains at most one destination for each transmitting site.
    slots = [[] for _ in range(max(map(len, edges)))]
    for u, row in enumerate(edges):
        for slot, (v, close) in zip(slots, row.items()):
            slot.append((u, v, close))
    slots = [np.array(s, dtype=int).T for s in slots]
    empty = np.full((n, len(home)), miss, dtype=np.int16)
    empty[home, np.arange(len(home))] = 0
    following = empty
    words = np.empty((last + 1, len(stations), len(home)), dtype=np.int16)
    for tick in range(last + limit, -1, -1):
        current = empty.copy()
        for u, v, close in slots:
            live = (close == 0) | (tick < close)
            a, b = u[live], v[live]
            current[a] = np.minimum(current[a], 1 + following[b])
        if tick <= last:
            words[tick] = current[stations]
        following = current
    return words


def record_service(edges, stations, home, last=LAST, limit=TIMEOUT):
    words = arrival_words(edges, stations, home, last, limit)
    label = {u: i for i, u in enumerate(stations)}
    runs = []
    for i, station in enumerate(stations):
        # At e-2 challenge each named peer; the e-1 immediate return is due at e.
        # Any longer itinerary misses this deadline. Peers identify themselves.
        offers = [
            (label[v], close, edges[v][station])
            for v, close in edges[station].items()
            if v in label
        ]
        prior = None
        changes = []
        for e in range(FIRST, last + 1):
            waits = [int(d) if d <= limit else None for d in words[e, i]]
            cert = None
            if waits[0] is not None:
                cert = sorted(
                    j
                    for j, there, back in offers
                    if (not there or e - 2 < there) and (not back or e - 1 < back)
                )
            if (waits, cert) != prior:
                changes.append([e, waits, cert])
                prior = (waits, cert)
        runs.append(changes)
    return runs


def sky_catalog(records, last=LAST, period=PERIOD):
    """Decode the registered time series using only labels, stamps and returned lists."""
    skies = [dict() for _ in range(last + 1)]
    for label, changes in enumerate(records):
        crests = {}
        for n, (start, waits, cert) in enumerate(changes):
            stop = changes[n + 1][0] if n + 1 < len(changes) else last + 1
            if waits[0] is None:
                continue
            for e in range(start, min(stop, last + 1 - waits[0])):
                t = e + waits[0]
                crests[e] = (t, cert)
                if e - period not in crests:
                    continue
                previous, first_cert = crests[e - period]
                earliest = min(x for x in waits if x is not None)
                bearings = tuple(j for j, x in enumerate(waits) if x == earliest)
                pattern = tuple(None if x is None else x - earliest for x in waits)
                if label in skies[t]:
                    raise ValueError("Two source crests share a receiver stamp")
                skies[t][label] = (
                    e,
                    t - previous - period,
                    bearings,
                    pattern,
                    set(cert).intersection(first_cert),
                )
    return skies


def compare(skies, gram, period=PERIOD):
    counts = []
    distinct = defaultdict(set)
    examples = {}
    supported = set()
    pulse_counts = Counter()
    for t, sky in enumerate(skies):
        stats = Counter()
        bin_centres = set()
        observed = Counter()
        for a, (e, shift, bearings, pattern, neighbours) in sky.items():
            pulse_counts[shift] += 1
            if len(bearings) == 1:
                bin_centres.add(bearings[0])
            for b in sorted(neighbours):
                if b <= a or b not in sky:
                    continue
                f, other_shift, other_bearings, other_pattern, reply = sky[b]
                if a not in reply:
                    continue
                pair = (a, b)
                supported.add(pair)
                if abs(e - f) > 1:
                    raise ValueError(
                        "Common reception violates the bilateral flight bound"
                    )
                if len(bearings) != 1 or len(other_bearings) != 1:
                    category = "tied"
                    angle = None
                elif bearings == other_bearings:
                    category = "same_bin"
                    angle = 0.0
                else:
                    category = "different_bins"
                    cosine = float(gram[bearings[0], other_bearings[0]])
                    angle = round(acos(max(-1.0, min(1.0, cosine))), 8)
                delta = abs(shift - other_shift)
                stats[category] += 1
                stats["length_certified"] += 1
                if delta:
                    stats["different_redshifts"] += 1
                if None not in pattern and None not in other_pattern:
                    stats["complete_arrays"] += 1
                    if pattern == other_pattern:
                        stats["identical_complete_array_patterns"] += 1
                stats["nonzero_span_in_both"] += int(bool(delta and angle))
                key = (min(shift, other_shift), max(shift, other_shift), angle)
                observed[key] += 1
                distinct[key].add(pair)
                for name, present in (
                    ("different_redshifts", delta),
                    ("different_bins", angle),
                    ("nonzero_span_in_both", delta and angle),
                ):
                    if present and name not in examples:
                        examples[name] = dict(
                            receiver_tick=t,
                            endpoints=[a, b],
                            later_emissions=[e, f],
                            earlier_arrivals=[
                                t - period - shift,
                                t - period - other_shift,
                            ],
                            z=[shift / period, other_shift / period],
                            first_apertures=[list(bearings), list(other_bearings)],
                            bin_centre_span_rad=angle,
                        )
        counts.append(
            dict(
                receiver_tick=t,
                sources_with_period=len(sky),
                occupied_single_bins=sorted(bin_centres),
                **stats,
                span_counts=[
                    [a, b, ang, n]
                    for (a, b, ang), n in sorted(
                        observed.items(), key=lambda kv: str(kv[0])
                    )
                ],
            )
        )
    return dict(
        pulse_counts_by_extra_ticks={
            str(k): v for k, v in sorted(pulse_counts.items())
        },
        unique_certified_pairs=len(supported),
        examples=examples,
        # A pair can occur in several bins over its history.
        distinct_pairs_by_reading=[
            [a, b, ang, len(pairs)]
            for (a, b, ang), pairs in sorted(
                distinct.items(), key=lambda kv: str(kv[0])
            )
        ],
        skies=counts,
    )


def measure(links, provenance, gram):
    # Every path reaching local port j in <= limit steps stays within limit+3
    # of the reference port in the all-transmitting graph. Verify the extra 3.
    ports, depth, edges = local_graph(links, 0, 2 * (TIMEOUT + 3))
    index = {p: i for i, p in enumerate(ports)}
    home = [index[j] for j in range(12)]
    if max(depth[i] for i in home) > 3:
        raise ValueError("Home aperture extends beyond the dependency ball")
    labels = sorted(p for p, d in zip(ports, depth, strict=True) if p // 12 and d <= 12)
    stations = [index[p] for p in labels]
    result = dict(
        provenance=provenance,
        receiver="0",
        reference_port=0,
        last_reception=LAST,
        last_emission=LAST,
        period=PERIOD,
        maximum_crest_flight=TIMEOUT,
        first_emission=FIRST,
        enrollment="All remote tags inside the predeclared twelve-step transparent aperture.",
        source_labels=[[str(p // 12), p % 12] for p in labels],
        possible_pairs=comb(len(labels), 2),
        local_direction_gram=gram.round(10).tolist(),
        record_fields=[
            "first_emission_using_word",
            "arrival_delays_at_home_ports_0_to_11",
            "indices_of_two_tick_source_echoes_reported_to_port_0",
        ],
        span_count_fields=[
            "lower_z_times_period",
            "upper_z_times_period",
            "separation_of_bin_centres_in_radians",
            "count",
        ],
    )
    for name, service in (
        ("strict", edges),
        ("transparent", [{v: 0 for v in row} for row in edges]),
    ):
        records = record_service(service, stations, home)
        result[name] = dict(records=records, **compare(sky_catalog(records), gram))
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--level13-archive",
        type=Path,
        help="Use the manifest-verified level-13 archive; default is the level-202 replay.",
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if args.level13_archive:
        source = ArchiveRadar(args.level13_archive)
        origin = dict(level=13, manifest_sha256=MANIFEST_SHA)
    else:
        from strict_history.wide_replay import SCHEME, WideReplay

        source = WideReplay(202, 110, SEED)
        origin = dict(level=202, seed=SEED, key_scheme=SCHEME, precision=110)
    result = measure(source.links, origin, local_gram())
    args.output.write_text(
        json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
    )
