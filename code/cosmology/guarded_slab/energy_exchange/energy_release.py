"""Sharp energy on a material release surface, with every enrollment retained."""

import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from energy_routes import HERE, ROOT, apparatus, dimer, public_capacity, rotate
from guarded_transfer import observations, source


def filter_labels(depth):
    _, decoder, leak = dimer.gates("native")
    observed = sum(dimer.outcomes(decoder, leak, i, depth) for i in (0, 1)) / 2
    expected = [F(1, 2 ** (j + 1)) for j in range(1, depth + 1)] + [
        F(1, 2),
        F(1, 2 ** (depth + 1)),
    ]
    np.testing.assert_allclose(observed, list(map(float, expected)), atol=2e-13, rtol=0)
    return expected


def certify_gates():
    k, _ = source()
    np.testing.assert_allclose(k, [[0.75, 0.25], [0.25, 0.75]], atol=2e-13, rtol=0)
    law = observations.native_law
    h0, h1 = law.local_term_parts()
    idx = [law.SRC.index((b, a)) for a, b in law.SRC]
    for h, goal in [
        (h0 + observations.KAPPA * h1, idx),
        (observations.KAPPA * h1, list(range(9))),
    ]:
        np.testing.assert_allclose(
            expm(-1j * np.pi * h)[:, :9], np.eye(law.DIM)[:, goal], atol=2e-13, rtol=0
        )
    h, decoder, _ = dimer.gates("native")
    levels, modes = np.linalg.eigh(h[:2, :2])
    np.testing.assert_allclose(levels, [32, 33], atol=2e-13, rtol=0)
    np.testing.assert_allclose(
        abs(decoder[:2, :2] @ modes) ** 2, np.eye(2), atol=2e-13, rtol=0
    )
    # Projectors of rail and decoded energy furnish Z and X, up to energy sign.
    z = np.diag([1.0, -1.0])
    x = decoder[:2, :2].conj().T @ z @ decoder[:2, :2]
    np.testing.assert_allclose(x, [[0, 1], [1, 0]], atol=2e-13, rtol=0)
    assert (
        np.linalg.matrix_rank(np.stack([v.ravel() for v in (np.eye(2), x, z, z @ x)]))
        == 4
    )
    # A point input distinguishes a swap from its parked exchange control.
    return dict(
        energy_levels=[32, 33], native_swap_columns=9, one_slot_algebra_dimension=4
    )


def stroke_records():
    """Propagate native complex amplitudes and recover the slice from the marker reads."""
    law = observations.native_law
    e, v = law.local_term_parts()
    words = list(product(range(2), repeat=3))
    pairs = list(product(range(2), repeat=2))
    rails = [(2, 0), (0, 2)]
    h, decoder, _ = dimer.gates("native")
    h = h[:2, :2]
    z = np.diag([1.0, -1.0])
    x = decoder[:2, :2].conj().T @ z @ decoder[:2, :2]
    flight = expm(-1j * np.pi / 2 * h)
    y = flight.conj().T @ z @ flight
    np.testing.assert_allclose(abs(y[0, 1].imag), 1, atol=1e-13)
    axes = (z, x, y)
    inputs = (
        np.array([1.0, 0.0]),
        decoder[:2, :2].conj().T[:, 1],
        flight.conj().T[:, 0],
    )
    reference = np.array(
        [[np.vdot(psi, axis @ psi).real for axis in axes] for psi in inputs]
    )
    state0 = np.kron(np.kron(inputs[0], inputs[1]), inputs[2])
    data = {}
    rng = np.random.default_rng(58125)
    for name, gen in [
        ("native", e + observations.KAPPA * v),
        ("exchange_off", observations.KAPPA * v),
    ]:
        pulse = expm(-1j * np.pi * gen)[:9, :9]
        # Two physical rail transpositions act on disjoint pairs of ports.
        swap = np.zeros((4, 4), complex)
        for j, src in enumerate(pairs):
            for i, dst in enumerate(pairs):
                swap[i, j] = np.prod(
                    [
                        pulse[
                            law.SRC.index(tuple(rails[b][r] for b in dst)),
                            law.SRC.index(tuple(rails[b][r] for b in src)),
                        ]
                        for r in range(2)
                    ]
                )

        def embed(u, v):
            return np.array(
                [
                    [
                        swap[pairs.index((a[u], a[v])), pairs.index((b[u], b[v]))]
                        if a[3 - u - v] == b[3 - u - v]
                        else 0
                        for b in words
                    ]
                    for a in words
                ]
            )

        stroke = embed(0, 2) @ embed(0, 1)
        idx = [law.SRC.index(pair) for pair in rails]
        meter = pulse[np.ix_(idx, idx)]
        position = 0
        elapsed = 0
        state = state0.copy()
        log = []
        for _ in range(6):
            state = stroke @ state
            prob = abs(meter[:, position]) ** 2
            read = int(np.argmax(prob))
            np.testing.assert_allclose(prob, np.eye(2)[read], atol=2e-13, rtol=0)
            elapsed += int(read != position)
            position = read
            means = []
            for slot in range(3):
                local = []
                for axis in axes:
                    factors = [np.eye(2) for _ in range(3)]
                    factors[slot] = axis
                    op = np.kron(np.kron(*factors[:2]), factors[2])
                    local.append(float(np.vdot(state, op @ state).real))
                means.append(local)
            np.testing.assert_allclose(
                means, rotate(list(reference), elapsed), atol=3e-12, rtol=0
            )
            # The prepared Pauli axes have exact plus probabilities 0, 1/2 or 1.
            # A tiny BLAS error near a deterministic endpoint must not advance
            # the random stream on only some machines.
            half_units = np.rint(1 + np.array(means))
            np.testing.assert_allclose(
                1 + np.array(means), half_units, atol=3e-12, rtol=0
            )
            counts = rng.binomial(4096, half_units / 2)
            log.append(
                dict(clock_rail=read, clock_count=elapsed, plus_counts=counts.tolist())
            )
        data[name] = log
    return dict(
        axes=["position", "energy", "quarter-period position"],
        shots_per_setting=4096,
        banks=data,
    )


def pair_after(cold, warm, steps, parked=False):
    odd = F(0) if parked else (1 - F(1, 2) ** steps) / 2
    weights = [
        (1 - cold) * (1 - warm),
        (1 - cold) * warm,
        cold * (1 - warm),
        cold * warm,
    ]
    return [
        weights[0],
        (1 - odd) * weights[1] + odd * weights[2],
        odd * weights[1] + (1 - odd) * weights[2],
        weights[3],
    ]


def cumulants(means):
    # Mixture of three conditionally independent binary readouts.
    avg = [sum(w * v[i] for w, v in means) for i in range(3)]
    cov = [
        sum(w * (v[i] - avg[i]) * (v[j] - avg[j]) for w, v in means)
        for i, j in combinations(range(3), 2)
    ]
    third = sum(
        w * (v[0] - avg[0]) * (v[1] - avg[1]) * (v[2] - avg[2]) for w, v in means
    )
    return dict(means=avg, cross_covariance=cov, connected_third=third)


def exact(cold, warm, cap, delays, parked=False):
    ages = [F(1, 4) * F(3, 4) ** (t - 1) for t in range(1, cap + 1)]
    hit = sum(ages)
    means = []
    words = list(product(range(4), repeat=3))
    law = [F(0)] * len(words)
    for t, w in enumerate(ages, 1):
        pairs = [pair_after(cold, warm, t + d, parked) for d in delays]
        means.append((w / hit, [p[2] + p[3] for p in pairs]))
        for i, (a, b, c) in enumerate(words):
            law[i] += w * pairs[0][a] * pairs[1][b] * pairs[2][c]
    assert sum(law) == hit
    moments = cumulants(means)
    # With independent hit ages, each seen output has its own independent age.
    independent = []
    for d in delays:
        independent.append(
            [
                sum(
                    w * pair_after(cold, warm, t + d, parked)[s]
                    for t, w in enumerate(ages, 1)
                )
                for s in range(4)
            ]
            + [1 - hit]
        )
    assert all(sum(p) == 1 for p in independent)
    return dict(
        hit=hit,
        timeout=1 - hit,
        words=law,
        independent_age_marginals=independent,
        delivered=moments,
    )


def unbounded(cold, warm, delays):
    def age_power(k):
        z = F(1, 2) ** k
        return z / (4 - 3 * z)

    v = age_power(2) - age_power(1) ** 2
    third = age_power(3) - 3 * age_power(1) * age_power(2) + 2 * age_power(1) ** 3
    a = (cold - warm) / 2
    b = (cold + warm) / 2
    return dict(
        means=[b + a * age_power(1) / 2**d for d in delays],
        cross_covariance=[a * a * v / 2 ** (i + j) for i, j in combinations(delays, 2)],
        connected_third=a**3 * third / 2 ** sum(delays),
    )


def draw_bank(rng, shots, depths, cap, delays):
    # First leaks carry their assay index and leak appointment. Other starts
    # proceed with a projective energy label restored by the inverse decoder.
    enrolled = np.ones(shots, dtype=bool)
    energy = np.zeros((shots, 6), dtype=np.uint8)
    failures = []
    for i, depth in enumerate(depths):
        alive = np.flatnonzero(enrolled)
        label = rng.choice(
            depth + 2, len(alive), p=list(map(float, filter_labels(depth)))
        )
        failures.append([int((label == j).sum()) for j in range(depth)])
        enrolled[alive[label < depth]] = False
        energy[alive, i] = label == depth + 1
    energy = energy[enrolled]
    # First root switch after arming, a geometric first-arrival sampler.
    age = rng.geometric(0.25, len(energy))
    timeout = int((age > cap).sum())
    energy, age = energy[age <= cap], age[age <= cap]
    before = energy.copy()
    odd = np.zeros((len(age), 3), dtype=np.uint8)
    for t in range(1, cap + max(delays) + 1):
        offer = (rng.random(odd.shape) < 0.25) & (t <= age[:, None] + np.array(delays))
        odd ^= offer.astype(np.uint8)
        cold, warm = energy[:, ::2].copy(), energy[:, 1::2].copy()
        energy[:, ::2] = np.where(offer, warm, cold)
        energy[:, 1::2] = np.where(offer, cold, warm)
    assert np.array_equal(
        energy[:, ::2] + energy[:, 1::2], before[:, ::2] + before[:, 1::2]
    )
    initial = before @ 2 ** np.arange(6)
    parity = odd @ 2 ** np.arange(3)
    hist = np.bincount(
        (age - 1) * 512 + initial * 8 + parity, minlength=cap * 512
    ).reshape(cap, 64, 8)
    assert int(hist.sum()) + timeout + sum(map(sum, failures)) == shots
    return hist.astype(np.uint32), dict(
        first_filter_leaks=failures,
        timeout=timeout,
        delivered=int(hist.sum()),
        launches=shots,
    )


def encode(obj):
    if isinstance(obj, F):
        return str(obj)
    raise TypeError(type(obj).__name__)


def run(out):
    protocol = json.loads((HERE / "protocol.json").read_text())
    cap = protocol["trigger"]["max_reads"]
    delays = protocol["trigger"]["capture_delay"]
    depths = [
        protocol["bodies"]["cold_filter_depth"],
        protocol["bodies"]["warm_filter_depth"],
    ] * 3
    gate = certify_gates()
    ready = F(1)
    marginal = []
    for depth in depths:
        weights = filter_labels(depth)
        ready *= sum(weights[-2:])
        marginal.append(weights[-1] / sum(weights[-2:]))
    cold, warm = marginal[:2]
    native = exact(cold, warm, cap, delays)
    controls = {
        name: exact(c, w, cap, delays, park)
        for name, c, w, park in [
            ("equal_temperature", warm, warm, False),
            ("parked", cold, warm, True),
            ("sharp", F(0), F(1), False),
        ]
    }
    rng = np.random.default_rng(protocol["sample"]["seed"])
    counts = {}
    boards = []
    for level in protocol["board_levels"]:
        board = apparatus(level)
        counts[f"L{level}"], board["sample"] = draw_bank(
            rng, protocol["sample"]["original_launches_per_board"], depths, cap, delays
        )
        boards.append(board)
    out.mkdir(exist_ok=True, parents=True)
    archive = out / "energy_counts.npz"
    np.savez_compressed(archive, **counts)
    result = dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "protocol.json").read_bytes()
        ).hexdigest(),
        input_hashes={
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                Path(__file__),
                HERE / "energy_routes.py",
                HERE.parent / "guarded_transfer.py",
                HERE.parent / "guarded_release.py",
                Path(dimer.__file__),
                Path(observations.native_law.__file__),
                Path(observations.__file__),
                ROOT / "model/executable_model/carrier.py",
            )
        },
        gates=gate,
        ring_tomography=stroke_records(),
        filter_acceptance=ready,
        upper_probabilities=[cold, warm],
        original_delivery=ready * native["hit"],
        cut_after_enrollment_miss=ready,
        native=native,
        controls=controls,
        unbounded=unbounded(cold, warm, delays),
        public_capacity=public_capacity(protocol),
        boards=boards,
        counts_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
    )
    (out / "energy.json").write_text(
        json.dumps(result, indent=2, default=encode) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=HERE)
    run(parser.parse_args().out)
