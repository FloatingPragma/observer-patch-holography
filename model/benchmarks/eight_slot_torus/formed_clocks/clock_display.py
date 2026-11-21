"""Portable count replays for the loaded and driven matched-pair instruments."""
import hashlib
import io
import json
import time

import numpy as np

from . import clock_transport, loaded_clock
from .. import runs
from ..law import LAWS


PORTS = ("+X", "+Y", "-X", "-Y")
LOADED_RUN = "10-loaded-clock.run.json.gz"
TRANSPORT_RUN = "11-driven-clock.run.json.gz"


def loaded_document():
    """All seven original-start outcomes at one enrolled stationary station."""
    started = time.monotonic()
    score = loaded_clock.experiment()
    protocol = json.loads(loaded_clock.BOOK.read_text(encoding="utf-8"))
    trials = protocol["trials_per_setting"]
    settings = protocol["fixtures"][0]
    focus = protocol["focus_fractions"].index([.75, .5])
    member = loaded_clock.READS.index("exact")
    payload = loaded_clock.ARCHIVE.read_bytes()
    if hashlib.sha256(payload).hexdigest() != score["archive_sha256"]:
        raise ValueError("loaded clock archive changed during display preparation")
    with np.load(io.BytesIO(payload), allow_pickle=False) as archived:
        keys = archived["keys"]
        selected = np.all(keys[:, :4] == (0, member, focus, 0), axis=1)
        ages = keys[selected, 4].tolist()
        counts = archived["counts"][selected]
    if ages != protocol["terminal_ages"] or not np.all(counts.sum(axis=1) == trials):
        raise ValueError("incomplete enrolled root counter record")
    ticks = [age + 1 for age in ages]
    series = [{"name": f"Port {label}", "points": [[tick, 100 * int(word[i]) / trials]
               for tick, word in zip(ticks, counts)]} for i, label in enumerate(PORTS)]
    labels = [*PORTS, "Born, missed corner", "Other first marker", "Unmarked"]
    table = [["Local tick", "Detector outcome", "Count"]]
    for tick, word in zip(ticks, counts):
        table.extend([tick, name, f"{value:,}"] for name, value in zip(labels, word))
    doc = runs.scan_document(
        "Four ports read a stationary pair clock", "Clocks",
        f"One fixed station counts the relative rail phase of a pair formed in a collision. "
        f"Each point uses {trials:,} fresh ordinary launches. The first marker admits three in "
        "sixteen at the selected mouth. The 28-service local preparation puts one quarter of "
        "that pair at the observed corner. The four curves show counts as a percentage of all "
        "launches. The table includes both first-marker refusals and every missed corner. "
        "Ticks are services counted at the station.",
        "Station counter (services)", "Percent of original launches", series, [],
        LAWS["torus"], time.monotonic() - started, table)
    doc["observer"] = {"carrier": settings["carrier"], "station": settings["root"],
                       "focus": protocol["focus_fractions"][focus], "read": "exact",
                       "trials_per_visit": trials, "counts": counts.tolist(),
                       "outcomes": labels, "port_order": list(PORTS), "station_counts": ticks, "archive_sha256": score["archive_sha256"],
                       "gates_missed": score["failures"]}
    return doc


def driven_document():
    """Reunion hand differences from the local raw counts for both path orderings."""
    started = time.monotonic()
    summary, counts = clock_transport.published_replay()
    gates = summary["gates"]
    root = clock_transport.PROTOCOL.parent
    protocol = json.loads(clock_transport.PROTOCOL.read_text(encoding="utf-8"))
    trials = protocol["trials"]
    shape = (len(protocol["routes"]["forward"]), *protocol["carrier"], 8)
    records = {}
    for name, word in summary["words"].items():
        c = counts[name]
        failed = c[2:-5]
        if failed.size != np.prod(shape) or int(c.sum()) != trials:
            raise ValueError(f"incomplete driven-clock record: {name}")
        # The ordering is (first hop, x, y, z, private, slot). All omitted bins are zero.
        nz = np.flatnonzero(failed)
        records[name] = {"unmarked": int(c[0]), "other_first_marker": int(c[1]),
                         "first_failures": [[int(i), int(failed[i])] for i in nz],
                         "ports": c[-5:-1].tolist(), "missed_corner": int(c[-1]),
                         "station": word["receiver"], "station_count": word["station_count"]}
    tours = summary["counts"]["reunions"]
    series, table = [], [["Setting", "Reading", "Value"]]
    differences = {}
    def hand(counts):
        return np.arctan2(counts[1] - counts[3], counts[0] - counts[2])

    for name, route in protocol["routes"].items():
        points = []
        for spacing in protocol["service_spacings"]:
            moving = records[f"{name}_{spacing}_sieved_{len(route)}_home"]["ports"]
            parked = records[f"{name}_{spacing}_parked_{len(route)}_home"]["ports"]
            difference = float(np.angle(np.exp(1j * (hand(moving) - hand(parked)))))
            comparison = next(row for row in tours if row["route"] == name and row["spacing"] == spacing)
            if abs(difference - comparison["measured_difference"]) > 2e-12:
                raise ValueError("reported reunion hand differs from its four physical counters")
            points.append([spacing, difference * 1000])
            radius = comparison["halfwidth"] * 1000
            differences.setdefault(spacing, []).append((difference * 1000, radius))
            table.append([f"{name}, {spacing}", "Difference (mrad)", f"{difference*1000:.5f} ± {radius:.5f}"])
        series.append({"name": name.capitalize() + " circuit", "points": points})
    bounds = [{"name": label, "points": [[s, select(d + sign * r for d, r in readings)]
               for s, readings in differences.items()]} for label, select, sign in
              (("Lowest interval bound", min, -1), ("Highest interval bound", max, 1))]
    hops = len(protocol["routes"]["forward"])
    expected = 1000 * hops * (summary["transported_increment"] - summary["parked_increment"])
    bounds.append({"name": "Finite-law forecast", "points": [[s, expected] for s in differences]})
    # Give one complete displaced reading and its home null control in the visible table.
    for receiver in ("booked", "home"):
        word = records[f"forward_4_sieved_3_{receiver}"]
        station = "Midpoint" if receiver == "booked" else "Home"
        for label, count in zip(PORTS, word["ports"]):
            table.append([station, label, f"{count:,}"])
        for label, count in (("Missed corner", word["missed_corner"]),
                             ("Absorbed in transit", sum(n for _, n in word["first_failures"])),
                             ("Other first marker", word["other_first_marker"]), ("Unmarked", word["unmarked"])):
            table.append([station, label, f"{count:,}"])
    doc = runs.scan_document(
        "A driven pair clock returns to its starting station", "Clocks",
        f"A six-link rectangular route and its reverse carry the loaded clock back to the home "
        f"station. The graph subtracts the independently parked hand from the transported hand "
        f"at the same elapsed service count. The dashed lines show the finite forecast and the count "
        f"interval envelope. Each station word uses {trials:,} fresh ordinary launches. At spacing four, "
        "the midpoint receiver is three links from home. The table gives its counts and a separate "
        "home reading at local tick 141. Local steering plates, absorbing flags and station timing "
        "are supplied. The logs are compared after the run.",
        "Address services per directed link", "Reunion hand difference (mrad)", series, bounds,
        LAWS["torus"], time.monotonic() - started, table)
    doc["observer"] = {"carrier": protocol["carrier"], "home": protocol["root"],
                       "focus": protocol["focus"], "trials_per_visit": trials, "port_order": list(PORTS),
                       "first_failure_shape": list(shape), "first_failure_axes": ["hop", "x", "y", "z", "private", "slot"],
                       "words": records, "gates_missed": [k for k, v in gates.items() if not v]}
    return doc
