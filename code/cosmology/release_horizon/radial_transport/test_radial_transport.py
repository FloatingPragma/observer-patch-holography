"""Independent checks from deposited endpoint reads and directional marginals."""

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

R = json.loads(Path(__file__).with_name("records.json").read_text())
EDGES = R["edge_order"]


def xyz(v):
    return v // 25, (v // 5) % 5, v % 5


def radius(v):
    return 0 if v == -1 else 1 + max(xyz(v))


def crossing(n):
    return [idx for idx, (a, b) in enumerate(EDGES) if radius(a) <= n < radius(b)]


def parent_address(edge):
    return tuple(tuple(component // 2 for component in xyz(v)) for v in edge)


def test_complete_native_routes_and_independent_shell_cuts():
    assert (
        R["carrier_owners"] == 320 and len(EDGES) == len(set(map(tuple, EDGES))) == 301
    )
    assert EDGES[0] == [-1, 0]
    assert {tuple(e) for e in EDGES[1:]} == {
        (a, b)
        for a in range(125)
        for b in range(125)
        if all(x <= y for x, y in zip(xyz(a), xyz(b)))
        and sum(xyz(b)) - sum(xyz(a)) == 1
    }
    assert [len(crossing(n)) for n in range(1, 5)] == [3, 12, 27, 48]
    assert R["removed_edge_control"]["measured"] == [26, 26]


def test_every_edge_and_every_station_conserve_measured_pulses():
    for box in R["boxes"]:
        primary = box["cycle"]
        for n in range(1, 5):
            counts = [primary["edge_hits"][idx] for idx in crossing(n)]
            assert set(counts) == {432 // (3 * n * n)} and sum(counts) == 432
        for run in box.values():
            balance = Counter()
            for (a, b), f in zip(EDGES, run["edge_hits"], strict=True):
                balance[a] += f
                balance[b] -= f
            assert [balance[v] for v in range(-1, 125)] == run["vertex_outflow"]
            assert (
                balance[-1]
                == run["injections"]
                == sum(dict(run["terminal_reads"]).values())
            )
        assert all(
            f == 0
            for v, f in enumerate(primary["vertex_outflow"][1:])
            if radius(v) <= 4
        )
        assert all(radius(v) == 5 for v, _ in primary["terminal_reads"])
        assert set(box["single_face"]["edge_hits"][idx] for idx in crossing(1)) == {
            0,
            144,
        }
        assert box["absorb"]["vertex_outflow"][1 + 25] == -144
    assert R["boxes"][0] == R["boxes"][1]


def test_joint_crossings_against_complete_independent_cycle_product():
    for box in R["boxes"]:
        run = box["cycle"]
        pairs = run["paired_shells_2_4"]
        assert sum(h for _, _, h in pairs) == 432 and len(pairs) == 48
        assert all(
            tuple(map(xyz, EDGES[a])) == parent_address(EDGES[b]) for a, b, _ in pairs
        )
        independent_matches = sum(
            run["edge_hits"][a] * run["edge_hits"][b]
            for a in crossing(2)
            for b in crossing(4)
            if tuple(map(xyz, EDGES[a])) == parent_address(EDGES[b])
        )
        assert Fraction(independent_matches, 432**2) == Fraction(1, 12)


def test_census_from_return_records_and_independent_local_periods():
    h = R["history"]
    for box in h["boxes"]:
        measurements = {item["layer"]: item["axes"] for item in box["ranges"]}
        for pub in box["publications"]:
            axes = measurements[max(t for t in measurements if t <= pub["layer"])]
            widths = [
                r["elapsed_slots"] // (2 * R["slot_seam_exposures"]) for r in axes
            ]
            q = 1 if pub["layer"] <= 8 else 2 if pub["layer"] <= 16 else 4
            assert widths == [q] * 3 and all(r["reflection"] for r in axes)
            assert pub["sites"] == widths[0] * widths[1] * widths[2]
            assert pub["changes"] == pub["sites"] * q
            assert h["angular_control_per_face"][pub["layer"] - 1] == q**3
    signals = h["pulses"]
    for r, record in zip(
        (2, 4),
        sorted(signals["receivers"].values(), key=lambda x: x["arrival_layers"]),
        strict=True,
    ):
        received = record["local_clocks"]
        source = signals["emitter_clock"]
        assert [
            (b - a) // (d - c)
            for a, b, c, d in zip(received, received[1:], source, source[1:])
        ] == [r] * 3
        assert record["missing"] == 0
    for control in h["controls"]:
        assert control["cut_range"]["width"] is None
        assert not control["cut_range"]["reflection"]
        assert control["parked_changes"] == 252
        assert [r["width"] for r in control["rectangular_axes"]] == [2, 3, 4]
        assert control["rectangular_changes"] == 96


def test_thermal_contrasts_against_ranged_postreset_inputs():
    from hashlib import sha256

    t = R["thermal_ranging"]
    archive = Path(__file__).resolve().parents[1] / "background_clock" / "records.json"
    assert sha256(archive.read_bytes()).hexdigest() == t["energy_archive_sha256"]
    original = json.loads(archive.read_bytes())["boards"][-1]["thermal"]
    assert t["energy_words"] == original["layer"]["energy_words"]
    assert t["probes"] == original["cells"] and t["input_challenges"] == 320
    assert t["read_marker"] - t["reset_marker"] == 8
    assert [[a["start"] for a in g] for g in t["input_ranges"]] == original[
        "pulse_inputs_per_probe"
    ]
    assert not t["shared_inputs"]
    for record in [*t["cross_ranges"], *(r for g in t["input_ranges"] for r in g)]:
        d = sum(abs(a - b) for a, b in zip(xyz(record["start"]), xyz(record["target"])))
        assert record["return_strokes"] == 2 * d * R["slot_seam_exposures"]
    assert t["endpoint_radii"] == [3, 2, 3, 4]
    assert [r["distance"] for r in t["cross_ranges"]] == [3] * 4
    assert all(
        d <= max(a, b)
        for d, (a, b) in zip(
            [3] * 4,
            ((a, b) for a in t["endpoint_radii"][:2] for b in t["endpoint_radii"][2:]),
            strict=True,
        )
    )
    assert float(Fraction(t["observed_covariance"])) > t["simultaneous_99pct_halfwidth"]
