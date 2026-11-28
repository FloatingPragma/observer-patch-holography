"""A native five-vertex clock tests stochastic common elapsed time on the same board."""

from fractions import Fraction as F
from itertools import pairwise
from math import cos, factorial, inf, log, nextafter, pi, sin, sqrt

from .native_columns import check as check_columns

TIMES = (1, 4, 8, 16)
TOLERANCE = F(1, 10**28)


def edge(a, b):
    return (min(a, b), max(a, b))


def graphs(fixture, left, right, first, second, install=True):
    n = 240
    ring = fixture["ring_ports"]
    aperture = set(ring)
    cycle = {edge(a, b) for a, b in pairwise(ring + [ring[0]])}
    physical = {edge(*e) for e in fixture["exterior_edges"]}
    active = set(physical)
    for c in range(n // 12):
        word = first + 1 if c == left else second + 1 if c == right else 0
        blocked = {edge(*e) for e in fixture["perfect_matchings"][word]}
        for a, b in fixture["cell_edges"]:
            bond = edge(a + 12 * c, b + 12 * c)
            physical.add(bond)
            opened = (
                bond in cycle
                if install and (aperture & set(bond))
                else (a, b) not in blocked
            )
            if opened:
                active.add(bond)
    assert cycle <= physical
    adj = [[] for _ in range(n)]
    physical_degree = [0] * n
    for a, b in physical:
        physical_degree[a] += 1
        physical_degree[b] += 1
    for a, b in sorted(active):
        adj[a].append(b)
        adj[b].append(a)
    # Twice H, including the actual unit-load port and blank-seam energies.
    potential = [
        16 * (6 - d + len(neighbors)) for d, neighbors in zip(physical_degree, adj)
    ]
    if install:
        for i, p in enumerate(ring):
            assert (
                set(adj[p]) == {ring[i - 1], ring[(i + 1) % 5]} and potential[p] == 48
            )
    return adj, potential


def coefficients(radius):
    n, err = 0, F(1)
    while err > TOLERANCE or n < radius * max(TIMES) / 2:
        n += 1
        err *= F(radius * max(TIMES), 2 * n)
    degree = n - 1
    denom = 2**degree * factorial(degree)
    weights = {}
    for t in TIMES:
        weights[t] = [
            (-1) ** (j // 2) * t**j * (denom // (2**j * factorial(j))) for j in range(n)
        ]
    return degree, denom, weights


def projections(adj, diagonal, ports, degree, weights):
    launch, receiver = ports
    v = [0] * len(adj)
    v[launch] = 1
    received = [[], []]
    for n in range(degree + 1):
        received[0].append(v[launch])
        received[1].append(v[receiver])
        if n != degree:
            take = v.__getitem__
            v = [a * x - sum(map(take, ns)) for a, x, ns in zip(diagonal, v, adj)]
    out = {}
    for t, c in weights.items():
        probs = []
        for m in received:
            terms = [a * b for a, b in zip(m, c)]
            probs.append(sum(terms[::2]) ** 2 + sum(terms[1::2]) ** 2)
        out[t] = probs
    return out


def enclosure(value, error):
    return [nextafter(float(value - error), -inf), nextafter(float(value + error), inf)]


def ring_checks():
    # Restrict to the calibrated subspace. D = twice the clock Hamiltonian
    # with the independently checked scalar removed.
    def apply(x):
        return [2 * x[i] - x[i - 1] - x[(i + 1) % 5] for i in range(5)]

    for s in range(5):
        v = [int(j == s) for j in range(5)]
        x = apply(v)
        y = apply(x)
        assert apply(y) == [5 * b - 5 * a for a, b in zip(x, y)]
    gaps = [(5 - sqrt(5)) / 4, (5 + sqrt(5)) / 4]
    # Compute the actual detector probabilities from both evolved mode vectors.
    z0 = [1 / sqrt(5)] * 5
    reports = []
    for k, g in enumerate(gaps, start=1):
        wave = [sqrt(2 / 5) * cos(2 * pi * k * j / 5) for j in range(5)]
        assert max(abs(a - 2 * g * b) for a, b in zip(apply(wave), wave)) < 2e-15
        rows = {}
        for t in TIMES:
            phase = complex(cos(g * t), -sin(g * t))
            p = [abs(z0[0] + c * wave[0] * phase) ** 2 / 2 for c in (1, 1j, -1, -1j)]
            reconstructed = complex(p[0] - p[2], p[1] - p[3]) / (2 * z0[0] * wave[0])
            assert abs(reconstructed - phase.conjugate()) < 2e-15
            rows[str(t)] = [round(a, 12) for a in p]
        reports.append(rows)
    # Deliberate random one-period slip of clock 1, with probability one half.
    jump = 2 * pi / gaps[0]
    control = [(1 - cos(g * jump)) / 2 for g in gaps]
    assert abs(control[0]) < 1e-14 and control[1] > 0.1
    return reports, control


def run(fixture):
    checked = check_columns()
    reports, control = ring_checks()
    fresh, plus = (fixture["stock_counts"][key] for key in ("fresh", "plus"))
    delta = [[b - a for a, b in zip(row_a, row_b)] for row_a, row_b in zip(fresh, plus)]
    assert all(sum(row) == 0 for row in delta) and all(
        sum(row) == 0 for row in zip(*delta)
    )
    reference = fixture["reference_cost_ticket_table"]
    # Every fixed-reference price is -log(count/total). The additive log(total)
    # cancels in a signed contrast with total weight zero.
    charge = (
        -sum(d * log(n) for ds, ns in zip(delta, reference) for d, n in zip(ds, ns))
        / fixture["denominator"]
    )
    assert abs(charge - 130 * log(39 / 35) / 20736) < 1e-16
    sites = fixture["source_pairs"]
    first_pair = sites[0]
    adj, potential = graphs(fixture, *first_pair, 0, 0)
    neighbor_counts = list(map(len, adj))
    lo = min(potential)
    hi = max(u + 2 * d for u, d in zip(potential, neighbor_counts))
    center = (lo + hi) // 2
    radius = max(hi - center, center - lo)
    degree, denom, weights = coefficients(radius)
    diagonal = [u + d - center for u, d in zip(potential, neighbor_counts)]
    totals = {t: [0, 0] for t in TIMES}
    active = 0
    for left, right in sites:
        for a, row in enumerate(delta):
            for b, count in enumerate(row):
                adj, potential_now = graphs(fixture, left, right, a, b)
                assert potential_now == potential and neighbor_counts == list(
                    map(len, adj)
                )
                if not count:
                    continue  # Equal ticket counts cancel exactly in the contrast.
                active += 1
                raw = projections(adj, diagonal, fixture["ports"], degree, weights)
                for t, values in raw.items():
                    for k in range(2):
                        totals[t][k] += count * values[k]
    l1 = F(sum(abs(i) for row in delta for i in row), fixture["denominator"])
    denominator = denom**2 * len(sites) * fixture["denominator"]
    result = {}
    for t in TIMES:
        err = F(radius * t, 2) ** (degree + 1) / factorial(degree + 1)
        prob_err = l1 * err * (2 + err)
        ranges = [enclosure(F(z, denominator), prob_err) for z in totals[t]]
        assert all(low * high > 0 for low, high in ranges)
        result[str(t)] = dict(zip(("near_return", "transmission"), ranges))
    # Optional finite-precision cut for |tau-t|<=1 and the t=16 transmission.
    # The (5+sqrt5)/4 clock gives 1-cos(gu)>=2g²u²/pi² in this window.
    g = (5 + sqrt(5)) / 4
    spatial_lipschitz = sqrt(neighbor_counts[fixture["ports"][1]]) / 2
    c = spatial_lipschitz * pi / (g * sqrt(2))
    lower_shift = min(abs(v) for v in result["16"]["transmission"])
    phi = (1 + sqrt(5)) / 2
    # F_34 turns of the low beat leave only phi**(-34) residual turns on the high beat.
    recurrence = 2 * sin(pi * phi ** (-34)) ** 2
    return {
        "native_columns_checked": checked,
        "charge_change_nats": round(charge, 15),
        "unbounded_near_recurrence": {
            "low_full_turns": 5702887,
            "high_phase_deficit": recurrence,
        },
        "geometry_checks": len(sites) * 81,
        "nonzero_signed_configs": active,
        "time_units": "native exchange",
        "clock_ring": fixture["ring_ports"],
        "gaps_numerical": [(5 - sqrt(5)) / 4, g],
        "four_clock_port_probabilities": reports,
        "half_probability_full_turn_deficits": control,
        "center_twice_H": center,
        "radius_twice_H": radius,
        "taylor_degree": degree,
        "probability_change": result,
        "max_shift_per_sqrt_phase_deficit": c,
        "bounded_window_phase_deficit_needed": (lower_shift / c) ** 2,
    }
