"""Driven translation of a loaded matched torus loop with common local plates.

The filter leaves all four successful corners coherent. Every occupied
rejected slot is absorbed and writes its site and first failed hop. The
prebooked plates do not use a global success bit during the journey.
"""
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from . import clock_loader as loader
from .. import law

PROTOCOL = Path(__file__).with_name("clock_transport_protocol.json")
PORTS = np.array([1, 1j, -1, -1j])
MEMBER = law.LAWS["torus"]
PUBLISHED_COUNTS = "b6e103a5c7a3e4534d84bbf363c9565ab387593ed828e3385ac550a789bf69aa"
PUBLISHED_LIKELIHOODS = "47bc95a34d77bab055f3f9769bff56a955eb9d7c7742a31436653be7187a9db3"


def common_ray(slot):
    """The optimal rail-independent input for the given routed output slot."""
    dest = np.eye(law.SLOTS, dtype=complex)[slot]
    w0, w1 = (law.contact_coin(dest, -phi) for phi in MEMBER.rail_contact_phases)
    overlap = np.vdot(w0, w1)
    ray = w0 + np.exp(-1j * np.angle(overlap)) * w1
    return ray / np.linalg.norm(ray)


def plate(start, end):
    """A common local unitary carrying an input ray to its phase-fixed image."""
    change = start - end
    if np.linalg.norm(change) < 1e-14:
        return np.eye(law.SLOTS)
    return np.eye(law.SLOTS) - np.outer(change, change.conj()) / (1 - np.vdot(end, start))


def _sites(carrier, root):
    loop = loader.circulation(carrier, root)
    return [(site, 2 * loop[site]) for site in loader.corners(carrier, root)]


def _step(states, age):
    for r in (0, 1):
        states[r] = law.step_maturing(states[r], MEMBER, age, r)


def translate(states, root, slot, sieve, age):
    """One physical mature service. Return both columns, root, all flag bins and slot power."""
    if age < MEMBER.maturation[2] or age >= 8192:
        raise ValueError("the controlled travel service must be mature and before the age wrap")
    carrier = states.shape[1:-1]
    incoming = common_ray(slot)
    outgoing = np.eye(law.SLOTS, dtype=complex)[int(law.PARTNER[slot])]
    for site, local in _sites(carrier, root):
        at = (slice(None), *site)
        states[at] = states[at] @ plate(local, incoming).T
    _step(states, age)
    root = loader.neighbour(root, slot, carrier)
    kept = np.zeros_like(states)
    sites = _sites(carrier, root)
    for site, _ in sites:
        at = (slice(None), *site, int(law.PARTNER[slot]))
        kept[at] = states[at]
    slot_power = float(np.vdot(kept, kept).real)
    if sieve:
        failed = np.sum(np.abs(states - kept) ** 2, axis=0)
        states = kept
    else:
        failed = np.zeros(states.shape[1:], dtype=float)
    for site, local in sites:
        at = (slice(None), *site)
        states[at] = states[at] @ plate(outgoing, local).T
    return states, root, failed, slot_power


def _read(states, site):
    rails = states[(slice(None), *site)]
    return np.sum(np.abs(rails[0] + PORTS[:, None] * rails[1]) ** 2, axis=-1) / 4


def _angle(ports):
    return float(np.angle((ports[0] - ports[2]) + 1j * (ports[1] - ports[3])))


def _wrap(value):
    return float(np.angle(np.exp(1j * value)))


def _outside(states, root):
    residual = states.copy()
    for site, ray in _sites(states.shape[1:-1], root):
        at = (slice(None), *site)
        residual[at] -= (states[at] @ ray.conj())[:, None] * ray
    return float(np.vdot(residual, residual).real)


def forecast():
    delta = np.diff(MEMBER.rail_contact_phases)[0] / 2
    d = np.exp(1j * delta) / 8 + 7 * np.exp(-1j * delta) / 8
    return float(delta), float(np.angle(d)), float((1 + abs(d)) / 2)


def _path(start, root, route, spacing, mode, protocol):
    states = start.copy()
    flag_shape = (len(route), *states.shape[1:])
    flags = np.zeros(flag_shape, dtype=float)
    actual, booked = root, root
    age = protocol["launch_age"]
    yield 0, age, actual, booked, states, flags, None, None
    for j, slot in enumerate(route):
        booked = loader.neighbour(booked, slot, states.shape[1:-1])
        slot_power = None
        if mode == "parked":
            _step(states, age)
        else:
            states, actual, flags[j], slot_power = translate(states, actual, slot, mode == "sieved", age)
        age += 1
        immediate = {"age": age, "booked_site": list(booked),
                     "booked_ports": _read(states, booked).tolist(),
                     "home_ports": _read(states, root).tolist(),
                     "surviving": float(np.vdot(states, states).real),
                     "first_failure": float(flags[j].sum())}
        for _ in range(spacing - 1):
            _step(states, age)
            age += 1
        yield j + 1, age, actual, booked, states, flags, slot_power, immediate


def replay(protocol=None):
    """Reproduce the local words and first failures for each independent terminal exposure."""
    if protocol is None:
        protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    carrier, root = tuple(protocol["carrier"]), tuple(protocol["root"])
    launch = loader.ordinary_launch(carrier, root, *protocol["focus"])
    # First the original, rail-conditioned natal loader. Travel uses only common plates.
    start, _ = loader.prepare_clock(carrier, root, "exact")
    start *= np.sqrt(launch["natal"] / 2)
    for age in range(loader.LOAD_AGE, protocol["launch_age"]):
        _step(start, age)
    delta, alpha, survival = forecast()
    initial = delta * (1 - (protocol["launch_age"] - loader.LOAD_AGE))
    arrays, words, first, arrivals = {}, {}, {}, {}
    offsets = protocol["appointments"]
    for name, route in protocol["routes"].items():
        for spacing in protocol["service_spacings"]:
            for mode in ("sieved", "unsieved", "parked"):
                prefix = f"{name}_{spacing}_{mode}"
                path = _path(start, root, route, spacing, mode, protocol)
                for j, age, actual, booked, states, flags, slot_power, immediate in path:
                    if immediate is not None:
                        if mode != "unsieved":
                            active = j if mode == "sieved" else 0
                            elapsed = immediate["age"] - protocol["launch_age"]
                            hand = initial - (elapsed - active) * delta + active * alpha
                            mass = launch["natal"] * survival ** active
                            worst = 0.
                            for site, label in ((booked, "booked_ports"), (root, "home_ports")):
                                incidence = mass / 4 if site in loader.corners(carrier, actual) else 0.
                                target = incidence * (1 + np.real(PORTS * np.exp(-1j * hand))) / 4
                                worst = max(worst, float(np.max(np.abs(immediate[label] - target))))
                            immediate["maximum_port_error"] = worst
                        arrivals[f"{prefix}_{j}"] = immediate
                    if j == 1 and mode != "parked":
                        first[prefix] = {"selected_slot_power": slot_power,
                                         "complement": launch["natal"] - slot_power}
                    if j not in offsets:
                        continue
                    remaining = float(np.vdot(states, states).real)
                    outside = _outside(states, actual)
                    for reader, site in (("booked", booked), ("home", root)):
                        key = f"{prefix}_{j}_{reader}"
                        ports = _read(states, site)
                        # Evaluate the miss from the unmonitored sites to avoid subtracting
                        # two nearly equal probabilities at an exactly localized source.
                        mask = np.ones(carrier, dtype=bool)
                        mask[site] = False
                        miss = float(np.sum(np.abs(states[:, mask, :]) ** 2))
                        # First two bins are the natal refusals. The final bin is terminal miss.
                        # All middle flag bins are indexed by (first hop, x,y,z,private,slot).
                        probabilities = np.r_[launch["unmarked"], launch["other"], flags.ravel(), ports, miss]
                        hand = _angle(ports) if ports.sum() > 1e-20 else None
                        arrays[key] = probabilities
                        separation = [abs(x - y) for x, y in zip(site, root)]
                        route_length = sum(min(d, period - d) for d, period in zip(separation, carrier))
                        word = {"age": age, "station_count": age + 1,
                                "earliest_causal_home_count": age + 1 + route_length,
                                "receiver": list(site), "occupied_root": list(actual),
                                "incidence": float(ports.sum()), "hand": hand,
                                "outside_local_loop_modes": outside, "surviving": remaining,
                                "flag_power_by_hop": np.sum(flags.reshape(len(route), -1), axis=1).tolist(),
                                "normalization_error": abs(float(probabilities.sum()) - 1)}
                        if mode != "unsieved":
                            forecast_hand = initial - j * spacing * delta
                            mass = launch["natal"]
                            if mode == "sieved":
                                forecast_hand += j * (alpha + delta)
                                mass *= survival ** j
                            incidence = mass / 4 if site in loader.corners(carrier, actual) else 0.
                            predicted = incidence * (1 + np.real(PORTS * np.exp(-1j * forecast_hand))) / 4
                            word["maximum_port_error"] = float(np.max(np.abs(ports - predicted)))
                        words[key] = word
    summary = {"natal": launch, "per_hop_survival": survival, "transported_increment": alpha,
               "parked_increment": -delta, "words": words, "arrivals": arrivals, "first_hop": first}
    return summary, arrays


def count_readout(summary, arrays, protocol, stored_counts=None):
    """Fresh multinomial panels. Use the four local counters for each reported angular interval."""
    rng = np.random.default_rng(protocol["seed"]) if stored_counts is None else None
    shots = protocol["trials"]
    total_bins = sum(len(p) for p in arrays.values())
    radius = np.sqrt(np.log(2 * total_bins / protocol["simultaneous_error"]) / (2 * shots))
    if stored_counts is not None and set(stored_counts) != set(arrays):
        raise ValueError("terminal count menu differs from the enrolled panel")
    counts = {}
    worst = 0.
    for key, p in arrays.items():
        if not np.isfinite(p).all() or p.min() < 0 or abs(p.sum() - 1) > 2e-11:
            raise ArithmeticError(f"incomplete original-start detector word: {key}")
        if stored_counts is None:
            c = rng.multinomial(shots, p)
        else:
            c = stored_counts[key]
        if (c.dtype.kind != "i" or c.shape != p.shape or c.min() < 0
                or c.max() > shots or c.sum() != shots):
            raise ValueError(f"invalid terminal counts: {key}")
        if np.any((p == 0) & (c != 0)):
            raise ValueError(f"count at a zero-probability outcome: {key}")
        counts[key] = c
        frequencies = c / shots
        worst = max(worst, float(np.max(np.abs(frequencies - p))))
        ports = frequencies[-5:-1]
        z = (ports[0] - ports[2]) + 1j * (ports[1] - ports[3])
        bound = 2 * np.sqrt(2) * radius
        angular = float(np.pi if abs(z) <= bound else np.arcsin(bound / abs(z)))
        summary["words"][key]["count_hand"] = None if z == 0 else float(np.angle(z))
        summary["words"][key]["count_halfwidth"] = angular
    comparisons = []
    for name, route in protocol["routes"].items():
        for spacing in protocol["service_spacings"]:
            key = f"{name}_{spacing}"
            a = summary["words"][f"{key}_sieved_{len(route)}_home"]
            b = summary["words"][f"{key}_parked_{len(route)}_home"]
            change = _wrap(a["hand"] - b["hand"])
            forecast_change = len(route) * (summary["transported_increment"] - summary["parked_increment"])
            x, y = a["count_hand"], b["count_hand"]
            read = None if x is None or y is None else _wrap(x - y)
            halfwidth = a["count_halfwidth"] + b["count_halfwidth"]
            comparisons.append({"route": name, "spacing": spacing, "exact_advance_difference": change,
                                "forecast_error": abs(_wrap(change - forecast_change)),
                                "measured_difference": read, "halfwidth": halfwidth,
                                "positive": bool(read is not None and halfwidth < read < np.pi - halfwidth)})
    summary["counts"] = {"exposures_per_word": shots, "family_bins": total_bins,
                         "simultaneous_radius": float(radius), "max_discrepancy": worst,
                         "all_bins_within_gate": bool(worst <= radius), "reunions": comparisons}
    return counts


def register(stored_counts=None):
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    summary, probabilities = replay(protocol)
    counts = count_readout(summary, probabilities, protocol, stored_counts)
    survival = summary["per_hop_survival"]
    natal = summary["natal"]["natal"]
    flag_errors = []
    for key, word in summary["words"].items():
        _, _, mode, appointment, _ = key.split("_")
        if mode == "sieved":
            for i, power in enumerate(word["flag_power_by_hop"]):
                expected = natal * survival ** i * (1 - survival) if i < int(appointment) else 0
                flag_errors.append(abs(power - expected))
    first = summary["first_hop"]
    first_errors = []
    for name in protocol["routes"]:
        for spacing in protocol["service_spacings"]:
            for mode in ("sieved", "unsieved"):
                row = first[f"{name}_{spacing}_{mode}"]
                first_errors.extend((abs(row["selected_slot_power"] - natal * survival),
                                     abs(row["complement"] - natal * (1 - survival))))
    max_word = lambda field: max(word.get(field, 0.) for word in summary["words"].values())
    reunions = summary["counts"]["reunions"]
    summary["gates"] = {
        "complete_words": max_word("normalization_error") <= 2e-11,
        "first_failures": max(flag_errors) <= 3e-12,
        "screened_and_parked_ports": max_word("maximum_port_error") <= 3e-12,
        "immediate_ports": max(row.get("maximum_port_error", 0.)
                               for row in summary["arrivals"].values()) <= 3e-12,
        "first_unsieved_complement": max(first_errors) <= 3e-12,
        "relative_signed_phase": max(row["forecast_error"] for row in reunions) <= 2e-10,
        "disjoint_midpoints": max(word["incidence"] for key, word in summary["words"].items()
                                  if "_3_home" in key and "_sieved_" in key or
                                  "_3_booked" in key and "_parked_" in key) <= 2e-12,
        "slow_reunion_counts": all(row["positive"] for row in reunions if row["spacing"] >= 4),
        "complete_count_gate": summary["counts"]["all_bins_within_gate"],
    }
    return summary, probabilities, counts


def _published_npz(path, sha256):
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != sha256:
        raise ValueError(f"published archive SHA-256 mismatch: {path.name}")
    return np.load(io.BytesIO(payload), allow_pickle=False)


def published_replay(seed_replay=False):
    """Re-evolve the wave and score archived counts. The optional draw checks the seed."""
    root = PROTOCOL.parent
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    saved = json.loads(gzip.decompress((root / "clock_transport_records.json.gz").read_bytes()))
    with _published_npz(root / "clock_transport_counts.npz", PUBLISHED_COUNTS) as stored:
        counts = dict(stored)
    summary, probabilities, _ = register(counts)
    with _published_npz(root / "clock_transport_likelihoods.npz", PUBLISHED_LIKELIHOODS) as stored:
        if set(stored.files) != set(probabilities):
            raise ValueError("changed terminal likelihood menu")
        rng = np.random.default_rng(protocol["seed"]) if seed_replay else None
        for name, generated in probabilities.items():
            booked = stored[name]
            if (booked.shape != generated.shape or booked.dtype != np.float64 or
                    not np.allclose(generated, booked, rtol=1e-11, atol=2e-13)):
                raise ValueError(f"changed terminal likelihoods: {name}")
            if np.any((booked == 0) & (counts[name] != 0)):
                raise ValueError(f"count outside the original sampling support: {name}")
            if seed_replay and not np.array_equal(rng.multinomial(protocol["trials"], booked), counts[name]):
                raise ValueError(f"changed original draw: {name}")

    def compare(value, expected):
        if type(value) is not type(expected):
            raise ValueError("changed transport record type")
        if isinstance(expected, dict):
            if value.keys() != expected.keys():
                raise ValueError("changed transport record fields")
            for key in expected:
                compare(value[key], expected[key])
        elif isinstance(expected, list):
            if len(value) != len(expected):
                raise ValueError("changed transport record length")
            for item, other in zip(value, expected, strict=True):
                compare(item, other)
        elif isinstance(expected, float):
            if not np.isclose(value, expected, atol=2e-12, rtol=1e-11):
                raise ValueError("changed transport reading")
        elif value != expected:
            raise ValueError("changed transport decision")

    # Use the public serialized representation for numeric scalar types.
    compare(json.loads(json.dumps(summary, allow_nan=False)), saved)
    return saved, counts


def verify_records(seed_replay=False):
    summary, _ = published_replay(seed_replay)
    return summary["gates"]


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify all retained local counters and reports")
    parser.add_argument("--seed-replay", action="store_true", help="also redraw counts from the published source")
    args = parser.parse_args(argv)
    result, _ = published_replay(args.seed_replay)
    print(json.dumps(result["gates"] if args.check else result, indent=2, allow_nan=False))
    return int(args.check and not all(result["gates"].values()))


if __name__ == "__main__":
    raise SystemExit(main())
