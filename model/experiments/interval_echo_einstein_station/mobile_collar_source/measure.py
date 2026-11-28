"""Native formation of two intercell record correlations by freely moving matter."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.photon_identity.native_emission import (
    measure as emission,
)
from model.executable_model import federation


HERE = Path(__file__).resolve().parent
TIMES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0)
PROBE_TIMES = (0.5, 1.0, 2.0, 4.0)
SHOTS = 2**42
SEED = 763011
# Hoeffding per cell for 100,000 readings (more than the enrolled counts).
EPS = np.sqrt(np.log(2 * 100000 / 0.01) / (2 * SHOTS))
COIN_IDENTITY = np.eye(2, dtype=complex)
PAULI = (
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]]),
    np.diag([1.0, -1.0]),
)
ANALYZERS = tuple(tuple((COIN_IDENTITY + k * s) / 2 for k in (1, -1)) for s in PAULI)
SIGNS = np.array([1.0, -1.0])
ANALYZER = np.array(
    [[np.kron(p, q) for p in a for q in b] for a in ANALYZERS for b in ANALYZERS]
).reshape(3, 3, 2, 2, 4, 4)


def region(board):
    a, b = min(e for e in board.external if e[0] // 12 == 0)
    left = min(j for j in board.neighbours[a] if j // 12 == a // 12)
    right = min(j for j in board.neighbours[b] if j // 12 == b // 12)
    return left, a, b, right


def local_parts(charge, strength):
    """Embed both sign-even coins at each possible oriented endpoint word."""
    native = emission.four_h(strength)
    words = ((1, 1), (1 + charge, 1), (1, 1 + charge))
    coins = np.column_stack(
        [emission.memory_vector([(c, -1, 1), (c, 1, 1)]) for c in (0, 1)]
    ).real
    # Integer embeddings allow exact invariance and orientation tests.
    basis = np.zeros((emission.LAW.DIM, 6), int)
    for i, word in enumerate(words):
        basis[emission.columns(word), 2 * i : 2 * i + 2] = np.rint(
            np.sqrt(2) * coins
        ).astype(int)
    scaled = basis.T @ native @ basis
    np.testing.assert_array_equal(2 * native @ basis, basis @ scaled)
    blocks = (scaled / 8).reshape(3, 2, 3, 2).transpose(0, 2, 1, 3)
    np.testing.assert_allclose(blocks[0, 0], strength * COIN_IDENTITY)
    return blocks.real


def hamiltonian(board, route, charge, strength, *, center=True):
    arms = [tuple(sorted(route[:2])), tuple(sorted(route[2:]))]
    blank = np.array([e for e in board.seams if e not in arms])
    a, b = blank.T
    degree = np.bincount(blank.ravel(), minlength=board.n)
    adjacency = sparse.csr_matrix(
        (np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(board.n, board.n)
    )
    energies = (0.5 + strength) * degree + strength * (2 * charge + charge**2) * (
        6 - board.degree
    )
    position = sparse.diags(energies) - 0.5 * adjacency
    total = sparse.kron(position, sparse.eye(4), format="csr")
    local = local_parts(charge, strength)
    for slot, (i, j) in enumerate(arms):
        for dst, src in ((0, 0), (1, 1), (0, 1), (1, 0)):
            block = local[1 + dst, 1 + src].copy()
            if dst == src:
                block -= local[0, 0]
            active = (
                np.kron(block, np.eye(2)) if slot == 0 else np.kron(np.eye(2), block)
            )
            ends = (i, j)
            dyad = sparse.csr_matrix(
                ([1], ([ends[dst]], [ends[src]])), shape=(board.n, board.n)
            )
            total += sparse.kron(dyad, active, format="csr")
    # Uniform-material terms are multiples of the identity. Remove a total scalar.
    if center:
        total -= sparse.eye(total.shape[0], format="csr") * float(np.mean(total.diagonal()))
    return total


def starting_coins(route, charge):
    return tuple(int((a > b) ^ (charge == 1)) for a, b in (route[:2], route[2:]))


def flow(h, start, times):
    ptr = 0.0
    states = []
    for t in times:
        if t > ptr:
            start = expm_multiply((-1j * (t - ptr)) * h, start, traceA=0.0)
        states.append(start.copy())
        ptr = t
    return states


def reduce(state):
    matrix = state.reshape(-1, 4)
    return matrix.T @ matrix.conj()


def marginals(rho):
    tensor = rho.reshape(2, 2, 2, 2)
    return np.trace(tensor, axis1=1, axis2=3), np.trace(tensor, axis1=0, axis2=2)


def negativity(rho):
    transpose = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    return float(max(0.0, -np.linalg.eigvalsh(transpose)[0]))


def probabilities(rho):
    return np.einsum("ijabkl,lk->ijab", ANALYZER, rho).real


def tomography(rho, rng):
    prob = probabilities(rho)
    counts = np.empty((3, 3, 2, 2), np.int64)
    for i in range(3):
        for j in range(3):
            p = prob[i, j].ravel().clip(0, 1)
            counts[i, j] = rng.multinomial(SHOTS, p / p.sum()).reshape(2, 2)
    frequencies = counts / SHOTS
    # All individual moments use their own z-partner pilot.
    left = np.einsum("iab,a->i", frequencies[:, 2], SIGNS)
    right = np.einsum("jab,b->j", frequencies[2], SIGNS)
    joint = np.einsum("ijab,a,b->ij", frequencies, SIGNS, SIGNS)
    fit = np.eye(4, dtype=complex)
    for i in range(3):
        fit += left[i] * np.kron(PAULI[i], COIN_IDENTITY) + right[i] * np.kron(COIN_IDENTITY, PAULI[i])
        for j in range(3):
            fit += joint[i, j] * np.kron(PAULI[i], PAULI[j])
    fitted_cov = joint - np.outer(left, right)
    return dict(
        counts=counts.tolist(),
        zz_covariance=float(fitted_cov[2, 2]),
        negative_partial_transpose=negativity(fit / 4),
    ), fitted_cov


def source(h, route, charge, rng):
    coin = starting_coins(route, charge)
    basis = 2 * coin[0] + coin[1]
    initial = np.zeros(h.shape[0], complex)
    initial[4 * route[0] + basis] = 1
    rows = []
    states = {}
    for t, state in zip(TIMES, flow(h, initial, TIMES), strict=True):
        rho = reduce(state)
        states[t] = rho.copy()
        first, second = marginals(rho)
        delta = rho - np.kron(first, second)
        energy_cov = float(np.trace(np.kron(PAULI[2], PAULI[2]) @ delta).real)
        read, cov = tomography(rho, rng)
        if t == 4:
            held = cov.copy()
        rows.append(
            dict(
                exposure=t,
                generated_zz=energy_cov,
                state_minus_product_trace_distance=float(
                    np.linalg.norm(delta, ord="nuc") / 2
                ),
                negativity=negativity(rho),
                **read,
            )
        )
    return rows, states, held


def probe_effects(h, route):
    initial = np.zeros((h.shape[0], 4), complex)
    r = route[3]
    initial[4 * r : 4 * r + 4] = np.eye(4)
    for t, waves in zip(PROBE_TIMES, flow(h, initial, PROBE_TIMES), strict=True):
        column = waves.reshape(-1, 4, 4)
        for name, arrivals in (
            ("home_port", column[[route[0]]]),
            ("home_owner", column[:12]),
        ):
            matrix = arrivals.reshape(-1, 4)
            yield t, name, matrix.conj().T @ matrix


def pilots(effect, rng):
    # Probability at each independently prepared product of two Pauli eigenstates.
    probs = probabilities(effect).clip(0, 1)
    raw = rng.binomial(SHOTS, probs)
    coefficients = np.einsum("ijab,a,b->ij", raw / SHOTS, SIGNS, SIGNS) / 4
    return raw, coefficients


def read_probe(effects, rho, cov, rng):
    product = np.kron(*marginals(rho))
    rows = []
    for t, kind, e in effects:
        raw, coefficients = pilots(e, rng)
        p = np.trace(e @ rho).real.clip(0, 1)
        pp = np.trace(e @ product).real.clip(0, 1)
        counts = rng.binomial(SHOTS, [p, pp])
        predicted = float(np.sum(coefficients * cov))
        # Each Pauli product or marginal is a difference of binary parity read frequencies;
        # bounding cells individually costs at most 4 eps, product covariance 12 eps.
        # The independent effect pilots cost eps per coefficient after division by four.
        forecast_bound = float(
            12 * EPS * np.abs(coefficients).sum()
            + EPS * np.abs(cov).sum()
            + 108 * EPS**2
        )
        rows.append(
            dict(
                exposure=t,
                detector=kind,
                pilot_counts=raw.tolist(),
                raw_yes_joint_and_product=counts.tolist(),
                simulated_joint_probability=float(p),
                simulated_product_probability=float(pp),
                exact_joint_minus_product=float(p - pp),
                forecast=predicted,
                forecast_radius_99=forecast_bound,
                measured=float((counts[0] - counts[1]) / SHOTS),
            )
        )
    return rows


def uniform_controls(rng, arms):
    native = emission.four_h(emission.K)
    emission.certification(native)
    labels = np.repeat(
        [-1, *(abs(d) for _, d in emission.LAW.LABELS)], emission.LAW.N_SRC
    )
    blank = (labels == -1).astype(int)
    blank_commutator = (blank[:, None] - blank[None, :]) * native
    magnitude_commutator = (labels[:, None] - labels[None, :]) * native
    rows = []
    for coins in arms:
        v = np.eye(4)[2 * coins[0] + coins[1]]
        # Both endpoint words are (1,1). Their exact sign-even blocks are scalar.
        h = local_parts(-1, emission.K)[0, 0]
        uniform_h = np.kron(h, COIN_IDENTITY) + np.kron(COIN_IDENTITY, h)
        np.testing.assert_array_equal(uniform_h, 2 * emission.K * np.eye(4))
        rho = np.outer(v, v)
        rows.append(
            dict(
                coins=coins,
                uniform_seam_energy=float(uniform_h[0, 0].real),
                appointments=[
                    dict(exposure=t, **tomography(rho, rng)[0]) for t in TIMES
                ],
            )
        )
    return dict(
        blank_commutator_nonzero=int(np.count_nonzero(blank_commutator)),
        magnitude_commutator_nonzero=int(np.count_nonzero(magnitude_commutator)),
        arms=rows,
    )


def modular_read(states, rng):
    ref = np.diag(states[4.0]).real.reshape(2, 2).clip(0, 1)
    calibration = rng.multinomial(SHOTS, ref.ravel()).reshape(2, 2)
    freq = calibration / SHOTS
    left = freq.sum(axis=1)
    right = freq.sum(axis=0)

    # Invert the two-sided binomial Bernstein bound relative to the observed count.
    # |s-mu| <= sqrt(2 mu x) + 2x/3, so mu is at most the positive root squared.
    def radius(counts):
        x = np.log(2 * 100000 / 0.01)
        root = (np.sqrt(2 * x) + np.sqrt(2 * x + 4 * (counts + 2 * x / 3))) / 2
        return (np.sqrt(2 * x) * root + 2 * x / 3) / SHOTS

    err = radius(calibration)
    err_left = radius(calibration.sum(axis=1))
    err_right = radius(calibration.sum(axis=0))
    result = dict(
        calibration_counts=calibration.tolist(),
        min_joint_probability=float(freq.min()),
        min_joint_lower_99=float((freq - err).min()),
        faithful_99=bool(np.all(freq > err)),
    )
    if not result["faithful_99"]:
        return result
    score = -np.log(freq / np.outer(left, right))
    # Every marginal is itself a registered binary read.
    score_error = (
        err / (freq - err)
        + err_left[:, None] / (left[:, None] - err_left[:, None])
        + err_right[None, :] / (right[None, :] - err_right[None, :])
    )
    result["connected_surprise_nats"] = score.tolist()
    rows = []
    for t in (2.0, 4.0, 8.0):
        p = np.diag(states[t]).real.reshape(2, 2).clip(0, 1)
        prod = np.outer(p.sum(axis=1), p.sum(axis=0))
        counts = np.array(
            [rng.multinomial(SHOTS, v.ravel()).reshape(2, 2) for v in (p, prod)]
        )
        diff = (counts[0] - counts[1]) / SHOTS
        diff_error = radius(counts[0]) + radius(counts[1])
        bound = np.sum(
            diff_error * np.abs(score) + (np.abs(diff) + diff_error) * score_error
        )
        rows.append(
            dict(
                exposure=t,
                joint_and_product_counts=counts.tolist(),
                measured_nats=float(np.sum(diff * score)),
                radius_99_nats=float(bound),
            )
        )
    result["held_reads"] = rows
    return result


def run():
    results = []
    for level in (0, 1):
        board = federation(level)
        route = region(board)
        ref = []
        for q in (-1, 1):
            rng = np.random.default_rng([SEED, level, q + 2])
            h = hamiltonian(board, route, q, emission.K)
            rows, states, cov = source(h, route, q, rng)
            ref.append((rows, states, cov))
            effects = list(probe_effects(h, route))
            if q == -1:
                vac = effects
            else:
                excess = effects
        arms = []
        for q, (rows, states, cov) in zip((-1, 1), ref, strict=True):
            rng = np.random.default_rng([SEED, level, 10 + q])
            rho = states[4]
            controls, _, _ = source(hamiltonian(board, route, q, 0), route, q, rng)
            arms.append(
                dict(
                    charge=q,
                    initial_coins=starting_coins(route, q),
                    source=rows,
                    modular=modular_read(
                        states, np.random.default_rng([SEED, level, 30 + q])
                    ),
                    changed_law_max_joint_distance=max(
                        r["state_minus_product_trace_distance"] for r in controls
                    ),
                    changed_law_zz_reads=[r["zz_covariance"] for r in controls],
                    probes={
                        name: read_probe(e, rho, cov, rng)
                        for name, e in (("vacancy", vac), ("excess", excess))
                    },
                )
            )
        sham = uniform_controls(
            np.random.default_rng([SEED, level, 40]),
            [starting_coins(route, q) for q in (-1, 1)],
        )
        results.append(dict(cells=board.cells, route=route, arms=arms, controls=sham))
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        modular_protocol_sha256=hashlib.sha256(
            (HERE / "MODULAR_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        native_contrast_strength=emission.K,
        native_law_sha256=hashlib.sha256(
            Path(emission.LAW.__file__).read_bytes()
        ).hexdigest(),
        shots=SHOTS,
        per_frequency_radius_99=float(EPS),
        covariance_radius_99=float(12 * EPS),
        negativity_radius_99=float(15 * EPS),
        probe_difference_radius_99=float(2 * EPS),
        experiments=results,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
