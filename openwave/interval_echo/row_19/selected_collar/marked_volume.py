"""Two observers compare the same six marked cells while ordinary matter enters."""

from fractions import Fraction as F
from math import factorial, log, sqrt

from mpmath import iv
import networkx as nx
import numpy as np

from . import stationary_work as work

local = work.local
bounds = work.bounds
SEED = 28122928
TIME = F(1, 8)
ORDER = 36
N = work.N
HALF = sqrt(log(2 * 10000 / 0.01) / (2 * N)) + 5e-13


def interval(q):
    return iv.mpf(q.numerator) / q.denominator


def certified_flight(h, start):
    # A constant energy shift is an unobservable phase. Evaluate the finite
    # polynomial by exact integer matrix products and rational coefficients.
    n = len(h)
    shift = sorted(int(x) for x in h.diagonal())[n // 2]
    b = h - shift * np.eye(n, dtype=int)
    norm = max(sum(int(abs(x)) for x in row) for row in b)
    size = TIME * norm / 2
    assert size < ORDER + 2
    tail = size ** (ORDER + 1) / factorial(ORDER + 1) / (1 - size / (ORDER + 2))
    power = np.zeros(n, dtype=object)
    power[start] = 1
    real = np.zeros(n, dtype=object)
    imag = real.copy()
    for k in range(ORDER + 1):
        coefficient = (-1) ** ((k + 1) // 2) * (TIME / 2) ** k / factorial(k)
        if k % 2:
            imag += coefficient * power
        else:
            real += coefficient * power
        power = b @ power
    probability = [a * a + b * b for a, b in zip(real, imag, strict=True)]
    error = 2 * tail + tail * tail
    assert abs(sum(probability) - 1) <= n * error
    return np.array([float(x) for x in probability]), probability, error, shift, norm


def score(counts, values):
    center = sum(iv.mpf(k) * int(v) / N for k, v in zip(counts, values, strict=True))
    width = (max(values) - min(values)) * HALF
    return center + iv.mpf([-width, width])


def spectra():
    law = local.survey.routes.volume.candidate.native_law
    eigen, modes = np.linalg.eigh(law.seam_term(8))
    energies = np.rint(eigen).astype(int)
    assert np.max(abs(eigen - energies)) < 1e-12
    result = {}
    for load in (1, 2):
        sharp = [modes[law.SRC.index((load, other))] ** 2 for other in (0, 1)]
        labels = [
            int(e)
            for e in sorted(set(energies))
            if any(p[energies == e].sum() > 1e-13 for p in sharp)
        ]
        result[load] = (
            labels,
            [np.array([p[energies == e].sum() for e in labels]) for p in sharp],
        )
    return result


def clock_read(pos, neighbors, spectrum, rng):
    energy = []
    raw = {}
    for q, (labels, (empty, occupied)) in spectrum.items():
        total = iv.mpf(0)
        raw[q] = []
        for port in neighbors:
            chance = pos[port]
            k = work.coherent.draws((1 - chance) * empty + chance * occupied, rng)
            total += score(k, labels)
            raw[q].append(k)
        energy.append(total)
    return energy[1] - energy[0], raw


def reference_weights(near, far):
    # A trial starts equally often at each port. In the first protocol the
    # negative and positive changes of shell+center/2 cancel.
    cn, sn = near
    cf, sf = far
    xn = sn + cn / 2
    xf = sf + cf / 2
    assert xn.b < 0 < xf.a and cn.a > 0 and cf.a > 0 and sf.a > 0
    u = xf / (xf - xn)
    central = (u * cn + (1 - u) * cf) / 2
    mass = central / 2
    v = 2 * mass / (sf + cf)
    assert u.a > 0 and u.b < 1 and v.a > 0 and v.b < 1
    return u, v, mass


def run():
    iv.dps = 55
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = work.steady.symmetric_instrument(board)
    _, _, ports, gates, partners, deg, edges, marks = setup
    g = gates[0]
    near = partners[0]
    graph = nx.Graph(edges)
    far = min(
        (
            i
            for i, d in nx.single_source_shortest_path_length(graph, g).items()
            if d == 2
        ),
        key=ports.__getitem__,
    )
    nb = sorted(graph[g])
    ball = [g, *nb]
    V = len(ball)
    n = len(ports)
    old = local.inputs("restored_balls")["surveys"][0]
    assert (
        sorted(ball) == old["ball_sites"][1] and old["volumes"][1] == V and near in nb
    )
    h = work.steady.native_energies(deg, edges, marks)[0, 0].h
    start = {}
    end = {}
    reads = {}
    exact = []
    measured = []
    for site in (near, far):
        initial = np.eye(n)[site]
        p, poly, err, shift, norm = certified_flight(h, site)
        start[site] = initial
        end[site] = p
        k = work.coherent.draws(p, rng)
        dc = interval(poly[g])
        ds = interval(sum(poly[i] for i in nb)) - int(site in nb)
        e = interval(err) * iv.mpf([-1, 1])
        exact.append((dc + e, ds + e * len(nb)))
        measured.append(
            (
                score(k, [int(i == g) for i in range(n)]),
                score(k, [int(i in nb) for i in range(n)]) - int(site in nb),
            )
        )
        reads[str(ports[site])] = {
            "position_counts": k,
            "certified_center_gain": bounds(dc + e),
            "certified_shell_change": bounds(ds + e * len(nb)),
        }
    u, v, mass = reference_weights(*exact)
    um, vm, _ = reference_weights(*measured)
    advance_near = float(um.mid)
    advance_far = float(vm.mid)
    common = (start[near] + start[far]) / 2

    def released(a, b):
        return (
            common + (a * (end[near] - start[near]) + b * (end[far] - start[far])) / 2
        )

    arms = {
        "baseline": common,
        "mixed_release": released(advance_near, 1 - advance_near),
        "far_release": released(0, advance_far),
    }
    spectrum = spectra()
    f0, cal0 = clock_read(start[far], nb, spectrum, rng)
    f1, cal1 = clock_read(start[near], nb, spectrum, rng)
    ratio_a = f1 / f0
    clocks_b, ratio_b = work.coherent.balls.proper_clock(rng)
    assert ratio_a.b < 1 and ratio_a.a > 0 and ratio_b.b < 1 and ratio_b.a > 0
    held = {}
    est = {}
    for key, p in arms.items():
        freq, k = clock_read(p, nb, spectrum, rng)
        hit = float(sum(p[i] for i in ball))
        local_count = work.coherent.draws([hit, 1 - hit], rng)
        census = work.coherent.draws(p, rng)
        est[key] = {
            "freq": freq,
            "mass_a": score(local_count, [1, 0]),
            "mass_b": score(census, [int(i in ball) for i in range(n)]),
            "shell_b": score(census, [int(i in nb) for i in range(n)]),
        }
        held[key] = {
            "frequency": bounds(freq),
            "local_clock_counts": k,
            "inside_outside": local_count,
            "port_counts": census,
        }

    def contrast(key, base):
        a, b = est[key], est[base]
        # The rate has exactly two sharp values. Interpolating its cube is
        # exact; cubing the ensemble-mean frequency would give a different read.
        da = V * (1 + ratio_a + ratio_a**2) * (a["freq"] - b["freq"]) / f0
        db = V * (ratio_b**3 - 1) * (a["shell_b"] - b["shell_b"])
        return {
            k: bounds(v)
            for k, v in {
                "clock_volume_change": da,
                "position_volume_change": db,
                "local_mass_change": a["mass_a"] - b["mass_a"],
                "census_mass_change": a["mass_b"] - b["mass_b"],
            }.items()
        }

    observations = {k: contrast(k, "baseline") for k in arms if k != "baseline"}
    return {
        "seed": SEED,
        "starts": N,
        "time": str(TIME),
        "taylor_order": ORDER,
        "generator_shift": shift,
        "shifted_norm": norm,
        "probability_error_bound": str(err),
        "ports": ports,
        "center": ports[g],
        "neighbors": [ports[i] for i in nb],
        "near": ports[near],
        "far": ports[far],
        "launches": reads,
        "certified_near_release_fraction": bounds(u),
        "certified_far_only_release_fraction": bounds(v),
        "measured_near_release_fraction": bounds(um),
        "measured_far_only_release_fraction": bounds(vm),
        "set_release_fractions": [advance_near, advance_far],
        "exact_match_mass_gain": bounds(mass),
        "predicted_matched_volume_changes": {
            "mixed_release": bounds(V * (1 - ratio_a**3) * mass),
            "far_release": bounds(V * (ratio_a**3 - 1) * v * exact[1][1] / 2),
        },
        "clock_energy_labels": {str(q): z[0] for q, z in spectrum.items()},
        "clock_calibration": [cal0, cal1],
        "frequencies": [bounds(f0), bounds(f1)],
        "position_observer_clocks": clocks_b,
        "held": held,
        "contrasts": observations,
        "between_releases": contrast("mixed_release", "far_release"),
        "inputs": work.steady.exp.save_inputs(
            "MARKED_VOLUME.md", "restored_balls.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    work.steady.exp.emit("marked_volume", result)
    print(
        result["exact_match_mass_gain"], result["contrasts"], result["between_releases"]
    )
