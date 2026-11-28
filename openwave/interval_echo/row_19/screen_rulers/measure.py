"""Compare two charge-specific clocks, receiver triangles and independently measured images."""

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path

from mpmath import iv
import numpy as np
from scipy.linalg import expm

from matched_einstein_limit.source_implication.material_clocks import measure as clocks
from ..selected_collar import measure as source, observer_records as camera

HERE = Path(__file__).parent
RECEIPT = HERE / "results.json"
SHOTS = 2**74
TIMES = [tick / 2 for tick in range(2, 13)] + [4.375]
ROUTES = tuple(permutations(range(3), 2))
STOCKS = ("reference", "ordinary_minus", "ordinary_plus", "held_minus", "held_plus")
PROBES = ("vacancy", "excess")
CORNERS = [0, 2, 4]


def preparations():
    ref, prior, *_ = camera.BEAM.aperture.source_and_changes()
    basis = source.joint_basis()
    held = (basis[0] + basis[4]) / np.sqrt(2)
    np.testing.assert_allclose(source.marginals(held), 0, atol=1e-14)
    stocks = {"reference": ref}
    steps = {"ordinary": 1 / 2048, "held": 1 / 128}
    for key, direction in (("ordinary", prior["ordinary"]), ("held", held)):
        for _ in range(21):
            endpoints = [ref + sign * steps[key] * direction for sign in (-1, 1)]
            if min(np.min(np.linalg.eigvalsh(x)) for x in endpoints) > 0:
                try:
                    recipes = [source.station.preparation(x) for x in endpoints]
                except ArithmeticError:
                    recipes = []
                if (
                    recipes
                    and min(min(r["min_kept"], r["min_dephased"]) for r in recipes) > 0
                ):
                    break
            if key == "ordinary":
                raise ValueError("Ordinary preparation differs from frozen source")
            steps[key] /= 2
        else:
            raise ValueError("No physical held endpoints")
        stocks[key + "_minus"], stocks[key + "_plus"] = endpoints
    recipes = {k: source.station.preparation(x) for k, x in stocks.items()}
    if max(r["max_error"] for r in recipes.values()) > 1e-11:
        raise ValueError("Native preparation error")
    return stocks, steps, recipes


def clock_probabilities(probe, state):
    if probe == "vacancy":
        h = clocks.blank_clock(-1, 8)
        return [abs(expm(-1j * t * h)[1, 0]) ** 2 for t in TIMES]
    parked = clocks.marginal_blocks(state, "right")
    return clocks.amplitudes(parked, clocks.contacts(0, 8), np.eye(1), TIMES, True)[
        :, 1, 0
    ]


def ranges(board, route, receivers, q, strength):
    ranging = source.station.iso.ranging
    # The extra gateway supplies the fourth kernel block. Only the three
    # receivers are interrogated in either leg of the monitored route.
    ports = [*receivers, route[2]]
    if len(set(ports)) != 4:
        raise ValueError("Ranging ports coincide")
    blocks = []
    for odd in (False, True):
        h = ranging.generator(board, route, q, strength, odd)
        u = ranging.kernels(h, ports, 0.5, 11)
        hits = {(i, j): ranging.first_hits(u, i, j) for i, j in ROUTES}
        edges = []
        for i, j in ROUTES:
            bins = []
            for k in range(2, 12):
                arms = np.array(
                    [hits[j, i][k - s] @ hits[i, j][s] for s in range(1, k)]
                )
                bins.append(np.einsum("kai,kaj->ij", arms.conj(), arms))
            edges.append(bins)
        blocks.append(edges)
    return np.stack(blocks)


def complete(rng, probs):
    a = np.asarray(probs)
    if a.ndim != 1 or min(a) < 0 or sum(a) > 1:
        raise ValueError("Invalid categorical probabilities")
    counts = camera.precision_draw(rng, SHOTS, a[:, None])[0]
    return counts, [*map(float, a), float(1 - a.sum())]


def collect(board, route, rings, stocks, strength, seed):
    rng = np.random.default_rng(seed)
    counts, population = {}, {}
    for q, probe in zip((-1, 1), PROBES, strict=True):
        receivers = [rings[1][j] for j in CORNERS]
        effects = ranges(board, route, receivers, q, strength)
        screen = camera.BEAM.effects(board, route, rings, q, TIMES[-1], strength)
        counts[probe], population[probe] = {}, {}
        for name, state in stocks.items():
            ps = list(map(float, clock_probabilities(probe, state)))
            count_panel = {
                "clock": [camera.precision_draw(rng, SHOTS, p) for p in ps],
                "range": [],
                "image": [],
            }
            prob_panel = {"clock": ps, "range": [], "image": []}
            signals = camera.BEAM.probabilities(screen, state)[np.ix_(CORNERS, CORNERS)]
            histogram = np.einsum("betij,bji->et", effects, state).real
            for key, rows in (("range", histogram), ("image", signals.T)):
                for row in rows:
                    bank, law = complete(rng, row)
                    count_panel[key].append(bank)
                    prob_panel[key].append(law)
            counts[probe][name] = count_panel
            population[probe][name] = prob_panel
    return {
        "counts": counts,
        "population": population,
        "strength": strength,
        "seed": seed,
    }


def probability(value, settings, counted, shots=SHOTS):
    if counted:
        if not isinstance(value, int) or not 0 <= value <= shots:
            raise ValueError("Invalid count")
        p = iv.mpf(value) / shots
        # Inverted empirical binomial tail, simultaneously over all named bins.
        tail = iv.ln(iv.mpf(2 * settings) / iv.mpf(".003"))
        r = (iv.sqrt(2 * iv.mpf(value) * tail) + 2 * tail) / shots
        bounds = [max(0, (p - r).a), min(1, (p + r).b)]
        return iv.mpf(bounds)
    if not 0 <= value <= 1:
        raise ValueError("Invalid model probability")
    return iv.mpf(value)


def hand(v, t):
    ascending = 2 * iv.atan2(iv.sqrt(v), iv.sqrt(1 - v))
    return ascending if t < np.pi else 2 * iv.pi - ascending


class UnresolvedMeasurement(ValueError):
    """A complete measurement whose finite records cannot define its geometric reading."""


def panel_read(panel, counted, ruler_shots=SHOTS):
    number = len(STOCKS) * len(PROBES)
    if len(panel["clock"]) != len(TIMES):
        raise ValueError("Missing clock appointment")
    observed = [
        hand(probability(v, number * len(TIMES), counted, ruler_shots), t)
        for v, t in zip(panel["clock"], TIMES, strict=True)
    ]
    if any((b - a).a <= 0 for a, b in zip(observed[:10], observed[1:11], strict=True)):
        raise UnresolvedMeasurement("Clock branches unresolved")

    def read_rows(tag, nrow, ncol):
        raw = panel[tag]
        if len(raw) != nrow:
            raise ValueError("Missing " + tag + " route or sender")
        for row in raw:
            shots = ruler_shots if tag == "range" else SHOTS
            target = shots if counted else 1
            if len(row) != ncol or (
                sum(row) != target if counted else abs(sum(row) - target) > 1e-13
            ):
                raise ValueError("Incomplete " + tag + " outcome bank")
            yield [probability(v, number * nrow * ncol, counted, shots) for v in row]

    directed = {}
    cap = observed[10] / 2
    for route, row in zip(ROUTES, read_rows("range", 6, 11), strict=True):
        directed[route] = cap + sum(row[k] * (observed[k] / 2 - cap) for k in range(10))
    a, b, c = [
        iv.sqrt(directed[i, j] * directed[j, i]) for i, j in ((0, 1), (0, 2), (1, 2))
    ]
    radicand = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)
    if radicand.a <= 0:
        raise UnresolvedMeasurement("Triangle inequality unresolved")
    area = iv.sqrt(radicand) / 4
    columns, through = [], []
    for row in read_rows("image", 3, 4):
        received = sum(row[:3])
        if received.a <= 0:
            raise UnresolvedMeasurement("No resolved receiver hits")
        through.append(received)
        columns.append([x / received for x in row[:3]])
    # A 3-by-3 matrix with unit column sums acts affinely on barycentric
    # coordinates. Det equals its signed area gain on the unit-sum plane.
    a, b, c = columns
    det = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    return {
        "radar": area,
        "image": area * det,
        "image_gain": det,
        "clock": observed[-1],
        "throughput": sum(through) / 3,
    }


def pack(x):
    return str(x) if x is not None else None


def decode(bank, steps, counted, *, ruler_shots=SHOTS):
    if ruler_shots not in (SHOTS, 2**80):
        raise ValueError("Unregistered ruling budget")
    reads = {}
    for probe in PROBES:
        if set(bank[probe]) != set(STOCKS):
            raise ValueError("Missing preparation")
        reads[probe] = {}
        for stock in STOCKS:
            try:
                reads[probe][stock] = panel_read(
                    bank[probe][stock], counted, ruler_shots
                )
            except UnresolvedMeasurement as exc:
                return {
                    "unresolved_stock": stock,
                    "unresolved_probe": probe,
                    "tests": {"reason": str(exc)},
                }
    tests = {}
    for area in ("radar", "image"):
        slopes, ratios, fractions = {}, {}, {}
        for probe, rows in reads.items():
            slopes[probe] = {
                k: (rows[k + "_plus"][area] - rows[k + "_minus"][area]) / (2 * h)
                for k, h in steps.items()
            }
            ordinary, held = slopes[probe]["ordinary"], slopes[probe]["held"]
            ratios[probe] = None if 0 in ordinary else held / ordinary
            ref = rows["reference"][area]
            fractions[probe] = None if 0 in ref else held / ref

        def compare(values):
            if any(x is None for x in values.values()):
                return {"difference": None, "decision": "unresolved calibration"}
            delta = values["excess"] - values["vacancy"]
            return {
                "difference": pack(delta),
                "decision": "different" if 0 not in delta else "unresolved",
            }

        tests[area] = {
            "slopes": {
                p: {k: pack(v) for k, v in s.items()} for p, s in slopes.items()
            },
            "held_over_ordinary": {p: pack(x) for p, x in ratios.items()},
            "comparison": compare(ratios),
            "held_over_reference": {p: pack(x) for p, x in fractions.items()},
            "fractional_comparison": compare(fractions),
        }
    return {
        "readings": {
            p: {s: {k: pack(v) for k, v in d.items()} for s, d in r.items()}
            for p, r in reads.items()
        },
        "tests": tests,
    }


def run():
    stocks, steps, recipes = preparations()
    f = camera.BEAM.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = camera.BEAM.screens(board, route)
    banks = {}
    for k, seed in ((8, 190824), (0, 190825)):
        bank = collect(board, route, rings, stocks, k, seed)
        bank["counted"] = decode(bank["counts"], steps, True)
        bank["calculated"] = decode(bank["population"], steps, False)
        banks[str(k)] = bank
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "receiver_triangle": [rings[1][i] for i in CORNERS],
        "sender_triangle": [rings[0][i] for i in CORNERS],
        "half_steps": steps,
        "shots_per_setting": SHOTS,
        "recipes": recipes,
        "banks": banks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        if (
            data["protocol_sha256"]
            != sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
            or data["shots_per_setting"] != SHOTS
        ):
            raise ValueError("Protocol differs from archive")
        for b in data["banks"].values():
            if b["counted"] != decode(b["counts"], data["half_steps"], True):
                raise ValueError("Decoded intervals differ from archive")
    else:
        data = run()
        if args.write:
            with RECEIPT.open("x") as f:
                json.dump(data, f, indent=2, allow_nan=False)
                f.write("\n")
        elif data != json.loads(RECEIPT.read_text()):
            raise ValueError("Acquisition differs from archive")
    print(
        json.dumps(
            {
                k: {s: v[s]["tests"] for s in ("counted", "calculated")}
                for k, v in data["banks"].items()
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
