"""Six directed round-trip clock panels and their sharp reciprocal calibration."""

from itertools import permutations
import json
from math import log, sqrt

from mpmath import iv
import numpy as np
from scipy import sparse

from ..measure import leg_kernels, pairing, ranging
from .geometry import ALL, HERE, certify, interval_numbers, raw_bounds

SHOTS = 2**31
SEED = 603927
ORDER = tuple(permutations(ALL, 2))


def return_histograms(exchange=True):
    common = json.loads((HERE.parents[1] / "full_algebra" / "results.json").read_text())
    state = np.asarray(common["x_real"]) + 1j * np.asarray(common["x_imag"])
    board = pairing.formation.federation(1)
    source = pairing.formation.region(board)
    ports = tuple(source) + (5, 6)
    bins = np.zeros((len(ORDER), 11))
    for odd, block in enumerate(state):
        for charge in (-1, 1):
            gen = ranging.generator(board, source, charge, 8, bool(odd))
            if not exchange:
                # Park every transfer between different spatial ports; retain
                # all native four-by-four Hamiltonian blocks at each port.
                coo = gen.tocoo()
                keep = coo.row // 4 == coo.col // 4
                gen = sparse.csr_matrix(
                    (coo.data[keep], (coo.row[keep], coo.col[keep])), shape=gen.shape
                )
            waves = leg_kernels(gen, ports)
            first = {(i, j): ranging.first_hits(waves, i, j) for i, j in ORDER}
            for k, (i, j) in enumerate(ORDER):
                for n in range(2, 12):
                    for s in range(1, n):
                        a = first[j, i][n - s] @ first[i, j][s]
                        bins[k, n - 2] += np.trace(a @ block @ a.conj().T).real / 2
    bins[:, -1] = 1 - bins[:, :-1].sum(axis=1)
    if bins.min() < 0:
        raise ValueError("incomplete terminal clock distribution")
    return bins


def clock_floor(means):
    """Half the maximum directed log-ratio cycle mean; all simple cycles suffice."""
    ratios = {(i, j): iv.ln(means[i, j] / means[j, i]) for i, j in ORDER}
    candidates = []
    for n in range(2, 7):
        for path in permutations(ALL, n):
            if path[0] != min(path):
                continue
            value = sum(ratios[path[k - 1], path[k]] for k in range(n)) / (2 * n)
            candidates.append((path, value))
    strongest = max(candidates, key=lambda item: float(item[1].a))
    largest = max(candidates, key=lambda item: float(item[1].b))
    bounds = iv.mpf([strongest[1].a, largest[1].b])
    # With a0=0, every optimal clock ai differs from zij (j=0) by at most 2*delta.
    gain = [iv.mpf(0)] + [
        ratios[i, 0] + iv.mpf([-2, 2]) * bounds.b for i in range(1, 6)
    ]
    energy = {
        (i, j): means[i, j] * means[j, i] * iv.exp(-gain[i] - gain[j]) for i, j in ORDER
    }
    if max(square.b for square in energy.values()) >= 4 * min(
        square.a for square in energy.values()
    ):
        raise ValueError("a reciprocal fit may violate the direct smooth chamber bound")
    error = max(sum(abs(energy[i, j] - 9).b for j in ALL if j != i) for i in ALL)
    gram_lower = (9 - error) / 2
    if gram_lower.a <= 0:
        raise ValueError("a reciprocal fit may leave the six-point simplex domain")
    return dict(
        witness_indices=list(strongest[0]),
        log_error=interval_numbers(bounds),
        all_optimal_gram_lower=interval_numbers(gram_lower)[0],
    )


def collect():
    iv.dps = 35
    original = raw_bounds()
    data = certify(original)
    data["clock_fit"] = clock_floor(original)
    observed = return_histograms()
    parked = return_histograms(exchange=False)
    if np.any(parked[:, :-1]):
        raise ValueError("parked control returned before the cap")
    rng = np.random.default_rng(SEED)
    populations = np.stack([rng.multinomial(SHOTS, p) for p in observed])
    times = np.arange(2, 13) / 4
    expected = observed @ times
    readings = populations @ times / SHOTS
    radius = (3 - 0.5) * sqrt(log(2 * len(ORDER) / 0.01) / (2 * SHOTS))
    # 1.8e-6 covers the normalized state-center error at the actual stationary
    # source. The two propagation schemes are checked against its original
    # independent interval calculation for all thirty ranges.
    err = iv.mpf(str(radius)) + iv.mpf("0.00000181")
    measured = np.full((6, 6), iv.mpf(0), dtype=object)
    for m, (i, j) in enumerate(ORDER):
        if not (original[i, j].a <= expected[m] <= original[i, j].b):
            raise ValueError("sparse return replay misses the integer-kernel enclosure")
        if abs(expected[m] - readings[m]) > radius:
            raise ValueError("fresh six-station histogram misses its simultaneous band")
        measured[i, j] = iv.mpf(float(readings[m])) + iv.mpf([-1, 1]) * err
    sampled = certify(measured)
    data["sampled"] = dict(
        seed=SEED,
        starts_per_pair=SHOTS,
        clock_radius=radius,
        terminal_counts=populations.tolist(),
        clock_fit=clock_floor(measured),
        **sampled,
    )
    data["parked_return_ranges"] = (parked @ times).tolist()
    return data


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = collect()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    elif data != json.loads(path.read_text()):
        raise ValueError(
            "six-station replay differs from saved interval bounds or counts"
        )
    print("Replayed six-station geometry and 30 independently counted clock panels.")
