"""Paired one-way pulse periods, two-way ranges and gated standard flashes at observer zero."""

import argparse
from collections import defaultdict, deque
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from observer_zero_radar import ArchiveRadar, MANIFEST_SHA, local_graph, outward

HERE = Path(__file__).resolve().parent
PERIOD = 4
TIMEOUT = 48
GATES = (16, 32)
SHOTS = 2**20
SEED = 130413


def inward(edges, source, launch, limit):
    """First receipt of a labelled crest from this source, on observer zero's clock."""
    if launch < 0:
        return None
    seen = {source}
    queue = deque([(source, launch)])
    while queue:
        u, t = queue.popleft()
        if u == 0:
            return t
        if t - launch >= limit:
            continue
        for v, closure in edges[u].items():
            if v not in seen and (not closure or t < closure):
                seen.add(v)
                queue.append((v, t + 1))
    return None


def timed_sources(ports, edges, sources, epoch):
    index = {p: i for i, p in enumerate(ports)}
    outward_waits = outward(edges, epoch, TIMEOUT)
    rows = []
    for p in sources:
        station = index[p]
        delay = outward_waits[station]
        row = dict(
            source=list(divmod(p, 12)),
            query_at=epoch,
            reflection_at=None,
            echo_at=None,
            beacon_sent=None,
            beacon_received=None,
            radar_distance=None,
            z_pulse=None,
            flash=None,
        )
        if delay >= 0:
            hit = epoch + delay
            launches = [hit - PERIOD // 2, hit + PERIOD // 2]
            arrivals = [inward(edges, station, t, TIMEOUT) for t in launches]
            echo = inward(edges, station, hit, TIMEOUT - delay)
            row.update(
                reflection_at=hit,
                echo_at=echo,
                beacon_sent=launches,
                beacon_received=arrivals,
            )
            if echo is not None:
                row["radar_distance"] = (echo - epoch) / 2
            if all(t is not None for t in arrivals):
                row["z_pulse"] = (arrivals[1] - arrivals[0]) / PERIOD - 1
        rows.append(row)
    return rows


def gateway_image(ports, edges, station, launch):
    """First observer-zero aperture to receive a label, retaining simultaneous bearings."""
    queue = deque([(station, launch)])
    seen = {station}
    first = None
    directions = []
    while queue:
        u, t = queue.popleft()
        if first is not None and t > first:
            break
        if ports[u] // 12 == 0:
            first = t
            directions.append(ports[u] % 12)
            continue
        if t >= launch + TIMEOUT:
            continue
        for v, stop in edges[u].items():
            if v not in seen and (not stop or t < stop):
                seen.add(v)
                queue.append((v, t + 1))
    return sorted(directions)


def ruler_images(ports, edges, rows, gram):
    # All twelve labelled ports of the emitting cell form the same unit-radius
    # slow-band calibration figure. Its physical transverse projection is uncalibrated.
    index = {p: i for i, p in enumerate(ports)}
    records = {}
    for row in rows:
        time = row["reflection_at"]
        if time is None:
            continue
        owner = row["source"][0]
        key = f"{owner}/{time}"
        row["ruler_exposure"] = key
        if key in records:
            continue
        labels = [
            gateway_image(ports, edges, index[12 * owner + p], time) for p in range(12)
        ]
        bearings = sorted({d for hit in labels for d in hit})
        if any(not hit for hit in labels):
            verdict = "incomplete"
        elif any(len(hit) > 1 for hit in labels):
            verdict = "tied_bearing"
        elif len(bearings) == 1:
            verdict = "unresolved"
        else:
            verdict = "multiple_bins"
        span = None
        if verdict in ("multiple_bins", "unresolved"):
            span = max(
                float(np.arccos(np.clip(gram[a, b], -1, 1))) if a != b else 0.0
                for a in bearings
                for b in bearings
            )
        records[key] = dict(
            first_bearings_by_label=labels, readout=verdict, angular_span_radians=span
        )
    return records


class FlashMeter:
    """Exact first-detection probabilities for independent single-tag trials."""

    def __init__(self, edges):
        self.n = len(edges)
        self.u = np.array(
            [u for u, row in enumerate(edges) for _ in row], dtype=np.int32
        )
        self.v = np.array([v for row in edges for v in row], dtype=np.int32)
        self.stop = np.array([t for row in edges for t in row.values()], dtype=np.int32)
        self.last = int(max(self.stop, default=0))
        if max(map(len, edges)) > 6:
            raise ValueError(
                "The six-slot source routing requires at most six contacts"
            )
        self.matrix = lru_cache(maxsize=None)(self._matrix)

    def _matrix(self, t):
        allowed = (self.stop == 0) | (t < self.stop)
        u, v = self.u[allowed], self.v[allowed]
        # One of six equiprobable offers. A closed or unused slot leaves the tag
        # in place. Boundary truncation is immaterial inside the flight-time bound.
        diagonal = 1 - np.bincount(u, minlength=self.n) / 6
        return sparse.csr_matrix(
            (np.full(u.size, 1 / 6), (u, v)), shape=(self.n, self.n)
        ) + sparse.diags(diagonal)

    def probabilities(self, requests):
        by_deadline = defaultdict(list)
        for key, station, launch, deadline in requests:
            by_deadline[deadline].append((key, station, launch))
        result = {}
        for end, items in by_deadline.items():
            p = np.zeros(self.n)
            p[0] = 1
            by_launch = defaultdict(list)
            for key, station, launch in items:
                by_launch[launch].append((key, station))
            for t in range(end, min(by_launch) - 1, -1):
                if t < end:
                    p = self.matrix(min(t, self.last)) @ p
                    p[0] = 1
                for key, station in by_launch.get(t, ()):
                    result[key] = float(p[station])
        return result


def photometry(rows, id_index, meter, rng):
    requests = []
    for i, row in enumerate(rows):
        if row["echo_at"] is None:
            continue
        owner, port = row["source"]
        for gate in GATES:
            requests.append(
                (
                    (i, gate),
                    id_index[12 * owner + port],
                    row["reflection_at"],
                    row["echo_at"] + gate,
                )
            )
    hit = meter.probabilities(requests)
    for i, row in enumerate(rows):
        if (i, GATES[0]) not in hit:
            continue
        p1, p2 = [hit[i, g] for g in GATES]
        if not 0 <= p1 <= p2 <= 1:
            raise ValueError("First-hit probabilities are inconsistent")
        a, b, _ = rng.multinomial(SHOTS, [p1, p2 - p1, 1 - p2])
        row["flash"] = dict(
            first_hit_probabilities=[p1, p2], counts=[int(a), int(a + b)]
        )


def measure(root):
    from model.candidate_laws.interval_echo.openwave.row_28_acceleration.experiment import (
        calibrate,
    )
    from observer_zero_information import local_gram

    calibrated_tags = calibrate()
    gram = local_gram()
    catalog = json.loads((HERE / "observer_zero_radar.json").read_text())
    archive = ArchiveRadar(root)
    if (
        catalog["manifest_sha256"] != MANIFEST_SHA
        or catalog["origin_port"] != 0
        or catalog["origin_observer"] != 0
    ):
        raise ValueError("Unexpected source catalogue")
    sources = [12 * o + p for o, p, _ in catalog["reflectors"]]
    ports, _, edges = local_graph(archive.links, 0, 2 * (TIMEOUT + 3))
    blank = [{v: 0 for v in row} for row in edges]
    panels = {str(e): timed_sources(ports, edges, sources, e) for e in (0, 4, 96)}
    for e, rows in panels.items():
        expected = catalog["echo_waits"][catalog["emission_ticks"].index(int(e))]
        actual = [None if r["echo_at"] is None else r["echo_at"] - int(e) for r in rows]
        if actual != expected:
            raise ValueError(
                "Echo times disagree with the independent deposited receiver"
            )
    panels["transparent"] = timed_sources(ports, blank, sources, 0)
    if any(
        r["z_pulse"] != 0 or r["radar_distance"] != entry[2]
        for r, entry in zip(panels["transparent"], catalog["reflectors"], strict=True)
    ):
        raise ValueError("Transparent timing controls failed")
    # A first hit observed by the deadline cannot visit a port outside the
    # intact ball of radius deadline minus launch about the receiving port.
    radius = max(
        r["echo_at"] - r["reflection_at"]
        for rows in panels.values()
        for r in rows
        if r["echo_at"] is not None
    ) + max(GATES)
    fp, _, fe = local_graph(archive.links, 0, 2 * radius)
    lookup = {p: i for i, p in enumerate(fp)}
    native = FlashMeter(fe)
    transparent = FlashMeter([{v: 0 for v in row} for row in fe])
    rng = np.random.default_rng(SEED)
    ref = timed_sources(ports, blank, [1], PERIOD)[0]
    if ref["radar_distance"] != 1 or ref["z_pulse"] != 0:
        raise ValueError("The local reference must be one radar unit from the receiver")
    photometry([ref], lookup, transparent, rng)

    def calibrate_flash(row):
        if row["flash"] is not None:
            row["flash"]["equivalent_distances"] = [
                (a / b) ** 0.5 if b else None
                for a, b in zip(
                    ref["flash"]["counts"], row["flash"]["counts"], strict=True
                )
            ]

    for name, rows in panels.items():
        photometry(rows, lookup, transparent if name == "transparent" else native, rng)
        for row in rows:
            calibrate_flash(row)
    images = {
        name: ruler_images(ports, blank if name == "transparent" else edges, rows, gram)
        for name, rows in panels.items()
    }
    return dict(
        manifest_sha256=MANIFEST_SHA,
        calibration_tag_transfers=calibrated_tags,
        ruler="Twelve labeled flashes emitted on the source cell at the reflection time in separate matched trials. Register the earliest gateway(s) into observer zero for each label. Unit slow-band directions are the angular bins.",
        ruler_exposures=images,
        archive="https://dovybzrnzxfo1.cloudfront.net/manifest.json",
        receiver=[0, 0],
        clock="Equal local reference flip rates; one scheduled blank transfer per tick. Only clock differences enter the timing ratios.",
        enrollment="All remote port labels in the deposited initial 12-relay aperture at receiver [0,0]. All queries, timeouts and missing points retained.",
        timing="The echo gives half the round trip. Four preprogrammed period-four source clocks on matched trials cover the four launch phases. Select the labelled pair bracketing the remote reflection by two ticks on each side. Pulse redshift is receiver period / emitted period - 1. A survey groups a common receiver query time; the pulse receipt times are source dependent.",
        flash_protocol="In separate trials at reflection, launch one standard tag from the same source port. Each tick choose uniformly among six slots. A transmitting contact swaps the tag; an unused or written contact retains it. Count first hits at observer zero. Counts integrate from the first possible return through each receiver-clock gate. Source output, aperture and timing are fixed. Gates share the same shots.",
        distance_calibration="Relative gated count distance is sqrt(local reference counts / source counts), with the transparent reference [0,1] at one radar unit. This is a fluence-equivalent distance. A physical bolometric luminosity or transverse length scale requires additional calibration.",
        seed=SEED,
        launches_per_source=SHOTS,
        gates=GATES,
        pulse_period=PERIOD,
        timeout=TIMEOUT,
        flash_ball_radius=radius,
        local_reference=ref,
        panels=panels,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("archive", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.write_text(
        json.dumps(measure(args.archive), indent=1, allow_nan=False) + "\n"
    )
