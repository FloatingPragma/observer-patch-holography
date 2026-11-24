"""Neutral-pair dressing of a period-three wall on ternary Interval Echo ports."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate, currents
from . import belt
from ..written_junction import polarization

LENGTH = 10
PRIME = 65521
SEED = 121313


def words(excess, holes):
    for zeros in combinations(range(LENGTH), holes):
        available = set(range(LENGTH)) - set(zeros)
        for twos in combinations(sorted(available), excess + holes):
            yield tuple(
                0 if i in zeros else 2 if i in twos else 1 for i in range(LENGTH)
            )


def addresses(word, excess):
    """Two disjoint exterior three-port clocks and the intervening four-port core."""
    start = 1 if excess == 3 else -1
    high = (0, 4, 7) if excess == 3 else (0, 1, 4, 7)
    reference = tuple(1 + int(i in high) for i in range(LENGTH))
    outside = [(start + i) % LENGTH for i in (-3, -2, -1, 4, 5, 6)]
    return [
        a
        for a in range(LENGTH)
        if all(word[(i + a) % LENGTH] == reference[i] for i in outside)
    ]


def necklace(sector):
    orbital = {}
    reps = []
    for word in sector:
        if word in orbital:
            continue
        rotated = word
        orbit = []
        while rotated not in orbital:
            orbital[rotated] = len(reps), len(orbit)
            orbit.append(rotated)
            rotated = rotated[-1:] + rotated[:-1]
        reps.append(orbit)
    return orbital, reps


def contacts(orbitals, reps):
    edges = belt.ring(LENGTH).seams
    moves = []
    for orbit in reps:
        word = orbit[0]
        hops = []
        cost = 2 * sum(x * x for x in word)
        for i, j in edges:
            a, b = word[i], word[j]
            cost += (a - b) ** 2
            if a != b:
                swapped = list(word)
                swapped[i], swapped[j] = b, a
                hops.append(orbitals[tuple(swapped)])
        moves.append((cost, hops))
    return moves


def modular_rank(h, seed):
    rng = np.random.default_rng(seed)
    v = rng.integers(1, PRIME, size=h.shape[0], dtype=np.int64)
    pivots = []
    for _ in range(len(v)):
        col = v.copy()
        for idx, pivot in pivots:
            col = (col - col[idx] * pivot) % PRIME
        entries = np.flatnonzero(col)
        if not len(entries):
            break
        idx = int(entries[0])
        pivots.append((idx, col * pow(int(col[idx]), -1, PRIME) % PRIME))
        v = (h @ v) % PRIME
    return len(pivots)


def reservoir_read(reps, probabilities, excess):
    """At a registered full-swap pulse, count both reservoir and outgoing wall."""
    result = {}
    delta = 1 if excess == 3 else -1
    for reservoir in range(3):
        records = {}
        for orbit, weight in zip(reps, probabilities, strict=True):
            for word in orbit:
                output = (reservoir,) + word[1:]
                change = reservoir - word[0]
                accepted = change == delta and bool(addresses(output, excess + delta))
                label = f"{change},{output.count(0)},{int(accepted)}"
                records[label] = records.get(label, 0.0) + weight / len(orbit)
        assert abs(sum(records.values()) - 1) < 1e-10
        result[str(reservoir)] = records
    return result


def scan(excess, holes, all_histograms=None):
    sector = list(words(excess, holes))
    orbital, reps = necklace(sector)
    moves = contacts(orbital, reps)
    decoded = [addresses(o[0], excess) for o in reps]
    # A sharp native wall with one or more neutral substitutions in its four-site core.
    source = next((o[0] for o, d in zip(reps, decoded) if d), None)
    source_index = orbital[source][0] if source is not None else None
    if holes == 0:
        high = (0, 4, 7) if excess == 3 else (0, 1, 4, 7)
        source = tuple(1 + int(i in high) for i in range(10))
        source_index = orbital[source][0]
    occupied = [len(o) for o in reps]
    prime_root = next(
        z
        for z in range(2, PRIME)
        if pow(z, 10, PRIME) == 1 and pow(z, 5, PRIME) != 1 and pow(z, 2, PRIME) != 1
    )
    records = {}
    stationary = {}
    stationary_probs = None
    if holes == 1:
        chart_start = 1 if excess == 3 else -1
        high = (0, 4, 7) if excess == 3 else (0, 1, 4, 7)
        dressed = [1 + int(i in high) for i in range(10)]
        neutral = [
            i % 10 for i in range(chart_start, chart_start + 4) if dressed[i % 10] == 1
        ]
        dressed[neutral[0]], dressed[neutral[1]] = 0, 2
        source = tuple(dressed)
        source_index = orbital[source][0]
        assert addresses(source, excess)
    k_ranks = []
    for k in range(6):
        subset = [j for j, length in enumerate(occupied) if k * length % 10 == 0]
        index = {j: i for i, j in enumerate(subset)}
        dim = len(subset)
        h_int = np.zeros((dim, dim), np.int64)
        hop = np.zeros((dim, dim), complex)
        cost = np.zeros(dim)
        for i, oid in enumerate(subset):
            c, swaps = moves[oid]
            h_int[i, i] = len(swaps) + 2 * int(candidate.KAPPA) * c
            hop[i, i] = len(swaps) / 2
            cost[i] = c
            for target, turn in swaps:
                if target not in index:
                    continue
                j = index[target]
                h_int[i, j] -= pow(prime_root, (k * turn) % 10, PRIME)
                hop[i, j] -= (
                    0.5
                    * np.sqrt(occupied[oid] / occupied[target])
                    * np.exp(2j * np.pi * k * turn / 10)
                )
        np.testing.assert_allclose(hop, hop.conj().T, atol=1e-12)
        block = sparse.csr_matrix(h_int % PRIME)
        if all_histograms is not None:
            all_histograms[k].append(block)
        rank = modular_rank(block, SEED + k)
        k_ranks.append(
            dict(
                momentum=k,
                dimension=dim,
                krylov_rank=rank,
                max_energy_multiplicity=dim - rank + 1,
            )
        )
        if source_index is None or source_index not in index:
            continue
        mask = np.array([bool(decoded[j]) for j in subset])
        for coupling in (candidate.KAPPA, 0):
            e, v = np.linalg.eigh(hop + coupling * np.diag(cost))
            # A sharp orbit word has weight 1/orbit_length per allowed momentum.
            weights = (abs(v[mask]) ** 2).sum(axis=0)
            # At exchange alone accidental degeneracies can mix these vectors.
            if coupling == candidate.KAPPA:
                peak = int(np.argmax(weights))
                if stationary.get("probability", -1) < weights[peak]:
                    stationary = dict(
                        momentum=k, energy=e[peak], probability=weights[peak]
                    )
                    stationary_probs = np.zeros(len(reps))
                    stationary_probs[subset] = abs(v[:, peak]) ** 2
            times = (0, 0.25, 1, 4, 16)
            phases = np.exp(-1j * (e - e[0])[:, None] * times)
            data = records.setdefault(str(coupling), np.zeros(len(times)))
            # Compute both momentum directions at the same positive exposure.
            # Their amplitudes can differ for an asymmetric sharp core.
            for modes in [v] if k in (0, 5) else [v, v.conj()]:
                coefficients = modes[index[source_index]].conj()
                output = modes[mask] @ (coefficients[:, None] * phases)
                data += np.sum(abs(output) ** 2, axis=0) / occupied[source_index]
    return dict(
        total=10 + excess,
        holes=holes,
        dimension=len(sector),
        admitted_words=sum(len(o) for o, d in zip(reps, decoded) if d),
        sharp_source=source,
        core_read_times=[0, 0.25, 1, 4, 16],
        joint_probabilities={k: v.tolist() for k, v in records.items()},
        native_max_stationary_joint_read=stationary,
        reservoir_delta_holes_accepted=(
            reservoir_read(reps, stationary_probs, excess)
            if stationary_probs is not None
            else {}
        ),
        momentum_certificate=k_ranks,
    )


def local_checks():
    # Histogram conservation on every seam, including all coherent memory labels.
    op0, op1 = law.local_term_parts()
    hole_count = np.tile([int(a == 0) + int(b == 0) for a, b in law.SRC], 11)
    diff = hole_count[:, None] - hole_count
    swap, _ = law._load_ops()
    target = np.vstack([swap, np.zeros((90, 9))])
    pulse_errors = {
        str(k): float(
            np.linalg.norm(expm(-1j * np.pi * (op0 + k * op1))[:, :9] - target, 2)
        )
        for k in (candidate.KAPPA, 0.0)
    }
    return dict(
        reservoir_swap_error=pulse_errors,
        exchange_commutator=float(np.max(abs(op0 * diff))),
        contrast_commutator=float(np.max(abs(op1 * diff))),
        addressed_echo_commutator=float(np.max(abs(law.echo() * diff))),
    )


def probe_energies(kappa):
    hx, hc = law.local_term_parts()
    # A stationary sample and parked reference arms select diagonal seam energies.
    diag = np.diag(hx + kappa * hc)[:9].reshape(3, 3)
    arms = diag[2] - diag[0]
    return arms - arms[1]


def read_phase(delta, rng):
    shots = 262144
    counts = [
        int(rng.binomial(shots, p))
        for p in ((1 + np.cos(delta)) / 2, (1 + np.sin(delta)) / 2)
    ]
    x, y = 2 * np.array(counts) / shots - 1
    return dict(
        bright=counts, shots_per_quadrature=shots, phase=float(np.arctan2(y, x))
    )


def long_cores():
    path, _, _ = polarization.layout()
    results = []
    rng = np.random.default_rng(SEED)
    gains = {k: probe_energies(k) for k in [candidate.KAPPA, 0.0]}
    time = np.pi / (6 * abs(gains[candidate.KAPPA][2]))
    for shift in (1, -2):
        x = polarization.X
        ref = 1 + (x % 3 == 0).astype(int)
        sample = ref.copy()
        for pos in np.flatnonzero((ref == 2) & (x >= 0) & (x < 120)):
            sample[pos] = 1
        for pos in np.flatnonzero((ref == 2) & (x >= 0) & (x < 120)):
            sample[pos + shift] = 2
        # Three overlapping windows give three counts for a calibrated unit.
        # Put decorations inside the flat part of this fixed aperture.
        aperture = abs(x) <= 18

        def score(v):
            return int(np.sum((np.roll(v, 1) + v + np.roll(v, -1) - 4)[aperture]))

        standard = ref.copy()
        standard[np.flatnonzero(x == 1)[0]] += 1
        calibration = score(standard) - score(ref)
        left = np.flatnonzero((x >= -23) & (x <= -21))
        right = np.flatnonzero((x >= 21) & (x <= 23))

        def clock(ids):
            return int(x[ids[sample[ids] == 2][0]] % 3)

        for holes in range(5):
            measurements = {}
            for coupling, gain in gains.items():
                field = gain[sample] - gain[ref]
                integral = np.sum(
                    (np.roll(field, 1) + field + np.roll(field, -1))[aperture]
                )
                # Opposite arm order makes an added unit a positive quarter turn.
                measurements[str(coupling)] = read_phase(-time * integral, rng)
            results.append(
                dict(
                    shift=shift,
                    phase_records=measurements,
                    holes=holes,
                    left_clock=clock(left),
                    right_clock=clock(right),
                    unit=calibration,
                    core_window_excess=score(sample) - score(ref),
                    counted_vacancies=int(np.sum(sample[aperture] == 0)),
                )
            )
            if holes < 4:
                a = next(
                    i
                    for i in np.argsort(abs(x))
                    if abs(x[i]) < 15 and sample[i] == sample[i + 1] == 1
                )
                sample[a], sample[a + 1] = 0, 2
                assert abs(x[a] - x[a + 1]) == 1 and path[a] != path[a + 1]
    return dict(
        cases=results,
        energies={str(k): g for k, g in gains.items()},
        exposure=time,
        seed=SEED,
        empty=read_phase(0, rng),
        independent_unit=read_phase(-3 * time * gains[candidate.KAPPA][2], rng),
    )


def run():
    readings = []
    combined = []
    for excess in (3, 4):
        histograms = [[] for _ in range(6)]
        readings.extend(scan(excess, holes, histograms) for holes in range(4))
        for k, blocks in enumerate(histograms):
            full = sparse.block_diag(blocks, format="csr")
            rank = modular_rank(full, SEED + k)
            combined.append(
                dict(
                    total=10 + excess,
                    momentum=k,
                    dimension=full.shape[0],
                    krylov_rank=rank,
                    max_energy_multiplicity=full.shape[0] - rank + 1,
                )
            )
    return currents._portable(
        dict(
            local_histogram=local_checks(),
            long_cores=long_cores(),
            prime=PRIME,
            seed=SEED,
            belts=readings,
            all_histograms=combined,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(result)
