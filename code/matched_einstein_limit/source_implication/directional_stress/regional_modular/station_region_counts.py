"""Read the untagged four-station density from occupations and six pair interferometers."""

import argparse
import json
from itertools import combinations

import numpy as np
from mpmath import iv
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from . import measure, station_region

c = measure.c


def draw(rng, trials, frequencies):
    # Integer totals beyond int64 are collected over independent, bounded
    # multinomials; summing as Python integers preserves each exact count.
    batch = min(2**60, trials)
    return (
        rng.multinomial(batch, frequencies, size=trials // batch)
        .astype(object)
        .sum(axis=0)
    )


def records(blocks, trials, rng):
    population = np.r_[
        1 - sum(b.trace().real for b in blocks), *[b.diagonal().real for b in blocks]
    ]
    assert (
        len(population) == 9
        and 0 < min(population)
        and abs(sum(population) - 1) < 1e-12
    )
    pos = draw(rng, trials, population)
    est = [np.diag(v.astype(complex) / trials) for v in (pos[1:5], pos[5:9])]
    hist = {"loads": pos.tolist()}
    for i, j in combinations(range(4), 2):
        for phase, typ in ((1, "real"), (-1j, "imaginary")):
            # The positive eigenvector of Y is (|i>+i|j>)/sqrt(2);
            # its expectation uses -Im(rho_ij).
            response = []
            for mat in blocks:
                total = (mat[i, i] + mat[j, j]).real
                score = 2 * (mat[i, j] * phase.conjugate()).real
                response.extend(((total + score) / 2, (total - score) / 2))
            response.append(1 - sum(response))
            assert min(response) > 0
            counts = draw(rng, trials, response)
            hist[f"{i}-{j}-{typ}"] = counts.tolist()
            for z, m in enumerate(est):
                m[i, j] += phase * (counts[2 * z] - counts[2 * z + 1]) / (2 * trials)
                m[j, i] = m[i, j].conjugate()
    return hist, block_diag([[pos[0] / trials]], *est)


def run():
    iv.dps = 32
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    launch = measure.bridge.beam.indices(sum(rings, []))
    ends = measure.bridge.beam.indices(route)
    fixed = measure.source()[0]
    state = np.array(
        [np.array(m, complex) * p for m, p in zip(fixed, (0.8657, 0.1343))]
    )
    waves = []
    for q in (-1, 1):
        parity = []
        for odd in (False, True):
            h = -2j * measure.bridge.native.generator(board, route, q, odd)
            initial = np.zeros((h.shape[0], 40))
            initial[launch] = np.eye(40)
            parity.append(
                expm_multiply(h, initial, traceA=h.diagonal().sum())[ends].reshape(
                    4, 4, 10, 4
                )
            )
        waves.append(np.array(parity))
    arms = station_region.controls()
    origin = {"ports": route}
    for label, shots, seed, ends in (
        ("original_equal_work", 2**60, 4013839, [-0.01, 0.01]),
        ("fixed_all_number_weights", 2**72, 4013840, [0, 0.01]),
    ):
        ds = np.array([np.array(v, complex) for v in arms[label][0]])
        out = compare(waves, state, ds, shots, seed, ends, label)
        if label == "original_equal_work":
            origin.update(out)
        else:
            origin[label] = out
    return origin


def compare(waves, state, ds, trials, seed, endpoints, label):
    rng = np.random.default_rng(seed)
    reads = {}
    # 9 population and 4 signed-bin probabilities per 12 pair settings,
    # at both endpoints with family confidence >=99%.
    bins = 2 * (9 + 12 * 4)
    eps = iv.sqrt(iv.ln(c.ri(2 * bins) / c.ri(".01")) / (2 * c.ri(trials))) + c.ri(
        "1e-14"
    )
    # One outside and 8 diagonal probabilities have error eps. 24 ordered
    # off-diagonals (across both sectors) have complex modulus error sqrt(2) eps.
    radius = iv.sqrt(c.ri(57)) * eps
    vals = []
    for key, t in zip(
        ("reference" if endpoints[0] == 0 else "minus", "plus"), endpoints
    ):
        rho = state + t * ds
        assert np.linalg.eigvalsh(rho).min() > 0
        occupied = [
            np.einsum("brhxi,bij,bshxj->rs", w, rho, w.conj()) / 20 for w in waves
        ]
        hist, m = records(occupied, trials, rng)
        # The two occupied blocks and one outside scalar exhaust the record.
        # Enclose the numerical eigenvalues of each empirical block using the
        # same rational polar residual as the population certificate.
        intervals = [c.ri(float(m[0, 0].real))]
        for inds in (slice(1, 5), slice(5, 9)):
            approx = c.prior.sym_matrix(m[inds, inds])
            intervals.extend(station_region.logarithm(approx, c.ri(0))[2])
        measured = -sum(v * iv.ln(v) for v in intervals)
        lower = min(x.a for x in intervals) - radius
        assert lower > 0
        # For trace-one densities |dS|<= ||d rho||_1 ||log rho||_op.
        # At most nine eigenmodes, so the trace norm is bounded by 3 times HS.
        error = -3 * iv.ln(lower) * radius
        vals.append(measured + iv.mpf([-1, 1]) * error)
        reads[key] = dict(setting_counts=hist, entropy_nats_99pct=str(vals[-1]))
    contrast = vals[1] - vals[0]
    assert contrast.a > 0 if label == "original_equal_work" else contrast.b < 0
    return dict(
        setting_count=13,
        charge_tag_retained=False,
        seed=seed,
        trials_per_setting=trials,
        simultaneous_coverage=0.99,
        arms=reads,
        source_control=label,
        source_endpoints=endpoints,
        signed_two_arm_entropy_difference_nats_99pct=str(contrast),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "station_region_counts.json").write_text(text)
    else:
        print(text, end="")
