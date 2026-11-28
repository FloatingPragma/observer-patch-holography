"""Test the regional entropy first law using phase-sensitive ordinary-matter records."""

from fractions import Fraction as F
from itertools import combinations
from math import log, sqrt

from mpmath import iv
import numpy as np

from . import local_trapping as trap
from . import restored_balls as balls

steady = trap.steady
local = trap.local
SEED = 28122919
BATCH = 2**56
REPEATS = 4096
N = BATCH * REPEATS


def draws(weights, rng):
    p = np.clip(weights, 0, 1)
    assert abs(sum(p) - 1) < 1e-11
    batches = rng.multinomial(BATCH, p / p.sum(), size=REPEATS)
    return [sum(map(int, col)) for col in batches.T]


def stationary_layouts(setup):
    _, _, ports, _, _, deg, edges, marks = setup
    ops = steady.native_energies(deg, edges, marks)
    n = len(ports)
    initial = np.ones((n, n), dtype=object), n
    layouts = {}
    for flags, op in ops.items():
        x, den = ops[0, 0].pinch(op.pinch(initial))
        layouts[flags] = np.vectorize(lambda v: float(F(v, den)))(x)
    return layouts, ops[0, 0].h


def stationary_states(board, setup):
    _, cells, ports, _, sites, _, _, _ = setup
    n = len(ports)
    layouts, _ = stationary_layouts(setup)
    deposit = np.zeros((n, n))
    deposit[sites[0], sites[0]] = 1
    forward, reverse, _ = steady.parity(board, cells)
    for direction, law in (("forward", forward), ("reversed", reverse)):
        change = sum(float(w) * (layouts[k] - layouts[0, 0]) for k, w in law.items())
        yield (
            direction,
            {
                "fresh": layouts[0, 0] + (25 / 36) * change,
                "shared": layouts[0, 0] + (101 / 144) * change,
                "ordinary": deposit,
            },
        )


def read_region(region, states, edges, rng, wait_error, trace_limits=None):
    ids = sorted(region)
    idx = {port: i for i, port in enumerate(ids)}
    internal = [(idx[a], idx[b]) for a, b in edges if a in idx and b in idx]
    fresh = states["fresh"][np.ix_(ids, ids)]
    _, modes = np.linalg.eigh(fresh)
    instrument, circuit = trap.local_rotation(modes.T, internal)
    rotated = {
        name: instrument @ data[np.ix_(ids, ids)] @ instrument.T
        for name, data in states.items()
    }
    # All bounds jointly use fewer than 100,000 individual probabilities
    # or signed quadratures. The latter have range two before division by two.
    half = sqrt(log(2 * 100000 / 0.01) / (2 * N)) + wait_error + 5e-13
    counts = {}
    prob = {}
    for key, x in rotated.items():
        counts[key] = draws([*x.diagonal(), 1 - np.trace(x)], rng)
        prob[key] = [
            iv.mpf([max(0, k / N - half), min(1, k / N + half)]) for k in counts[key]
        ]
    pair_reads = []
    row_sums = np.zeros(len(ids))
    x = rotated["fresh"]
    for i, j in combinations(range(len(ids)), 2):
        mean = (x[i, i] + x[j, j]) / 2
        pair = []
        for cross in (x[i, j].real, x[i, j].imag):
            k = draws([mean + cross, mean - cross, 1 - 2 * mean], rng)
            pair.append(k)
        observed = [abs(k[0] - k[1]) / (2 * N) + half for k in pair]
        bound = sum(observed)  # |Re|+|Im| bounds the modulus.
        row_sums[[i, j]] += bound
        pair_reads.append({"ports": [ids[i], ids[j]], "quadratures": pair})
    # Vacuum and one-load coherences are absent by the regional algebra.
    # Gershgorin bounds ||rho-diag rho||, then the integral formula for
    # log bounds ||log rho-log diag rho|| by e/(min(diag rho)-e).
    e = max(row_sums) * (1 + 1e-12) + 1e-17
    m = min(float(a.a) for a in prob["fresh"]) - e - 1e-16
    assert m > 0
    log_error = e / m + 1e-14
    reference = prob["fresh"]
    classical = -sum(v * iv.log(v) for v in reference)
    # S(diag rho)-S(rho)=D(rho||diag rho), in [0, log_error].
    entropy = classical + iv.mpf([-log_error, 0])
    slopes = {}
    if trace_limits is None:
        trace_limits = {"shared": F(1, 72), "ordinary": F(2)}
    for name, max_trace in trace_limits.items():
        # A pure admission change has trace norm at most twice its change
        # in probability. An arbitrary ordinary replacement is bounded by two.
        d = -sum(
            (a - b) * iv.log(b) for a, b in zip(prob[name], reference, strict=True)
        )
        bound = float(max_trace) * log_error
        slopes[name] = d + iv.mpf([-bound, bound])
    results = {
        "sites": ids,
        "circuit": circuit,
        "mode_and_outside_counts": counts,
        "fresh_phase_counts": pair_reads,
        "coherence_norm_upper": e,
        "reference_eigenvalue_lower": m,
        "log_norm_error_upper": log_error,
        "entropy": balls.bounds(entropy),
        "entropy_tangents": {k: balls.bounds(v) for k, v in slopes.items()},
    }
    return results, entropy, slopes


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = steady.symmetric_instrument(board)
    edges = setup[6]
    saved = local.inputs("restored_balls")["surveys"]
    expected = local.inputs("stationary_polarization")["mirror_matched"]
    error = float(F(expected["combined_trace_error_upper"]))
    regions, comparisons = {}, []
    for direction, states in stationary_states(board, setup):
        cached = {}
        for row in saved:
            if row["direction"] != direction:
                continue
            entropies, derivatives = [], []
            for sites in row["ball_sites"]:
                key = tuple(sites)
                if key not in cached:
                    result, ref, slopes = read_region(key, states, edges, rng, error)
                    cached[key] = ref, slopes
                    regions[f"{direction}_{','.join(map(str, key))}"] = result
                ref, slopes = cached[key]
                entropies.append(ref)
                derivatives.append(slopes)
            changes = {}
            for name, arm in (
                ("shared", "record_source"),
                ("ordinary", "ordinary_standard"),
            ):
                area = iv.mpf(row[arm]["area_tangent"])
                coin = iv.mpf(row[arm]["public_probability_tangent"])
                # At p=1/2 the public bit has zero entropy derivative.
                slope = sum(z[name] for z in derivatives) / 2 + coin * (
                    entropies[1] - entropies[0]
                )
                changes[name] = (area, slope)
            a, ds = changes["shared"]
            oa, os = changes["ordinary"]
            comparisons.append(
                {
                    "direction": direction,
                    "center": row["center"],
                    "hop_radii": row["hop_radii"],
                    "record_entropy_tangent": balls.bounds(ds),
                    "ordinary_entropy_tangent": balls.bounds(os),
                    "record_gain": balls.bounds(-ds / a) if 0 not in a else None,
                    "ordinary_gain": balls.bounds(-os / oa) if 0 not in oa else None,
                    "common_gain_determinant": balls.bounds(a * os - oa * ds),
                    "required_additional_entropy_contrast": (
                        balls.bounds(a * os / oa - ds) if 0 not in oa else None
                    ),
                }
            )
    return {
        "starts_per_setting": N,
        "seed": SEED,
        "regions": regions,
        "comparisons": comparisons,
        "inputs": trap.exp.save_inputs(
            "COHERENT_BALLS.md", "restored_balls.json", "stationary_polarization.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    trap.exp.emit("coherent_balls", result)
    for row in result["comparisons"]:
        print(
            row["direction"],
            row["center"],
            row["hop_radii"],
            row["record_gain"],
            row["ordinary_gain"],
            row["common_gain_determinant"],
        )
