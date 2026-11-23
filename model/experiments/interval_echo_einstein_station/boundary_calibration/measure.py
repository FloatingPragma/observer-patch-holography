"""Active exterior readings of the equal-energy orbital bit under the full native pair flow."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply
from scipy.special import entr

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    pair_index,
)
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
    compile_wave,
)
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.measure import (
    angular_frame,
)
from model.executable_model import federation
from model.experiments.interval_echo_einstein_station.hidden_entropy.measure import (
    two_dark_modes,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**35
DEPTH = 25
STEP = 0.125
SEED = 514283
TIMES = (0.25, 0.5)
AXES = np.eye(3)
REFERENCE = np.array([0.0, 0.0, 0.5])
CALIBRATION = np.vstack([sign * axis for axis in AXES for sign in (-1, 1)])
HELD = REFERENCE + STEP * CALIBRATION
PHASE_STATES = np.array([[0, 0.5, 0], [0, 0.5 - STEP, 0], [0, 0.5 + STEP, 0]])


def inputs(board):
    cell, frame, _, _ = angular_frame()
    edges = sorted((a, b) for a, b in board.external if a < 12)
    ports, launches = np.array(edges).T
    dark = two_dark_modes(frame, ports)
    # Every eigen-preparation used in a classical mixture compiles on the same cell.
    for bloch in np.vstack([CALIBRATION, REFERENCE, HELD, PHASE_STATES]):
        x, y, z = bloch
        rho = np.array([[1 + z, x - 1j * y], [x + 1j * y, 1 - z]]) / 2
        probabilities, vectors = np.linalg.eigh(rho)
        prepared = np.column_stack([compile_wave(cell, dark @ v) for v in vectors.T])
        np.testing.assert_allclose(
            (prepared * probabilities) @ prepared.conj().T,
            dark @ rho @ dark.T,
            atol=1e-13,
        )
    free = np.zeros((board.n, 2))
    free[:12] = dark
    np.testing.assert_allclose(
        candidate.generator(board) @ free,
        (6 * candidate.KAPPA + 3) * free,
        rtol=0,
        atol=2e-13,
    )
    return dark, launches


def product_ball(board, launches):
    distance = shortest_path(
        board.adjacency(), unweighted=True, indices=np.r_[np.arange(12), launches]
    )
    left, right = distance[:12].min(0), distance[12:].min(0)
    sites = np.flatnonzero(np.minimum(left, right) < DEPTH)
    first, second = np.triu_indices(len(sites), 1)
    first, second = sites[first], sites[second]
    retained = (
        np.minimum(left[first] + right[second], right[first] + left[second]) < DEPTH
    )
    return first[retained], second[retained]


def operators(board, first, second):
    n = board.n
    keys = first * n + second
    adjacent = (board.table[first] == second[:, None]).any(axis=1)
    entries, columns = [], []
    for mover, partner in ((first, second), (second, first)):
        for destinations in board.table[mover].T:
            legal = (destinations >= 0) & (destinations != partner)
            old = np.flatnonzero(legal)
            dest, other = destinations[legal], partner[legal]
            new_key = np.minimum(dest, other) * n + np.maximum(dest, other)
            row = np.searchsorted(keys, new_key)
            matches = (row < len(keys)) & (
                keys[np.minimum(len(keys) - 1, row)] == new_key
            )
            entries.append(row[matches])
            columns.append(old[matches])
    rows, cols = np.concatenate(entries), np.concatenate(columns)
    hops = sparse.coo_matrix(
        (np.full(len(rows), -0.5), (rows, cols)), shape=(len(keys), len(keys))
    ).tocsr()
    diagonal = (board.degree[first] + board.degree[second] - 2 * adjacent) / 2
    kinetic = hops + sparse.diags(diagonal)
    # Every source/body omitted from the product ball is reached after >= DEPTH hops.
    return kinetic, sparse.diags(-2 * candidate.KAPPA * adjacent), keys


def portals(first, second, launch):
    inside_a = first // 12 == launch // 12
    inside_b = second // 12 == launch // 12
    alphabet = np.zeros(len(first), dtype=int)
    alphabet[inside_a & ~inside_b] = 1 + first[inside_a & ~inside_b] % 12
    alphabet[inside_b & ~inside_a] = 1 + second[inside_b & ~inside_a] % 12
    both = inside_a & inside_b
    alphabet[both] = 13 + pair_index(12, first[both] % 12, second[both] % 12)
    return alphabet


def effects(packet, bins):
    first, second = packet.T
    diagonal = [
        np.bincount(bins, weights=abs(v) ** 2, minlength=79) for v in (first, second)
    ]
    interference = first.conj() * second
    real = np.bincount(bins, weights=interference.real, minlength=79)
    imag = np.bincount(bins, weights=interference.imag, minlength=79)
    a, b = diagonal
    response = np.column_stack([(a + b) / 2, real, -imag, (a - b) / 2])
    np.testing.assert_allclose(response.sum(axis=0), [1, 0, 0, 0], rtol=0, atol=2e-13)
    return response


def seed_packet(keys, n, dark, launches):
    initial = np.zeros((len(keys), 6), complex)
    for j, launch in enumerate(launches):
        row = np.searchsorted(keys, n * np.arange(12) + launch)
        np.testing.assert_array_equal(keys[row], n * np.arange(12) + launch)
        initial[row, 2 * j : 2 * j + 2] = dark
    return initial


def evolve(hamiltonian, packet, bins):
    rows = []
    step = -0.25j * hamiltonian
    for _ in TIMES:
        packet = expm_multiply(step, packet, traceA=step.diagonal().sum())
        rows.append(
            [
                effects(packet[:, 2 * j : 2 * j + 2], alphabet)
                for j, alphabet in enumerate(bins)
            ]
        )
    return np.array(rows)


def sample(response, targets, rng):
    laws = np.einsum(
        "tgbi,hi->htgb", response, np.column_stack([np.ones(len(targets)), targets])
    )
    counts = np.zeros(laws.shape, dtype=np.int64)
    for h, t, gate in np.ndindex(laws.shape[:-1]):
        law = laws[h, t, gate]
        assert law.min() > -1e-14
        cleaned = np.maximum(0.0, law)
        counts[h, t, gate] = rng.multinomial(SHOTS, cleaned / cleaned.sum())
    return counts


def h_bloch(vector):
    radius = min(1.0, np.linalg.norm(vector))
    return float(entr((1 + radius) / 2) + entr((1 - radius) / 2))


def regional_read(dark, states, rng):
    cell = angular_frame()[0]
    first, second = dark.T
    readers = np.column_stack(
        [
            compile_wave(cell, wave)
            for wave in (
                (first + second) / np.sqrt(2),
                (first + 1j * second) / np.sqrt(2),
                first,
            )
        ]
    )
    overlaps = dark.T @ readers
    counts = []
    for x, y, z in states:
        density = np.array([[1 + z, x - 1j * y], [x + 1j * y, 1 - z]]) / 2
        yes = np.sum(overlaps.conj() * (density @ overlaps), axis=0).real
        tallies = rng.binomial(SHOTS, yes)
        counts.append(np.stack([tallies, SHOTS - tallies], axis=-1))
    counts = np.array(counts)
    measured = 2 * counts[:, :, 0] / SHOTS - 1
    entropy = np.array([h_bloch(u) for u in measured])
    radii = np.linalg.norm(measured, axis=1)
    derivative = -np.arctanh(radii)[:, None] * measured / radii[:, None]
    variance = np.sum(derivative**2 * (1 - measured**2), axis=-1) / SHOTS
    slopes = np.diff(entropy[1:].reshape(-1, 2), axis=1).ravel() / (2 * STEP)
    errors = np.sqrt(variance[1:].reshape(-1, 2).sum(axis=1)) / (2 * STEP)
    return dict(
        counts=counts.tolist(),
        observed_xyz=measured.tolist(),
        entropy_nats=entropy.tolist(),
        slopes=slopes.tolist(),
        standard_errors=errors.tolist(),
    )


def noiseless_derivatives(response):
    late = response[-1]
    slopes = late[:, :, 1:]
    reference = late[:, :, 0] + slopes @ REFERENCE
    logs = np.zeros_like(reference)
    np.log(reference, out=logs, where=reference > 0)
    entropy = -np.einsum("gb,gbi->i", logs, slopes) / late.shape[0]
    # Whiten the population-and-coherence response on the same expected-count support
    # as the empirical inverse, without a noisy axis calibration.
    mask = response[:, :, :, 0] * SHOTS >= 100
    white = response[:, :, :, 1:][mask] / np.sqrt(response[:, :, :, 0][mask, None])
    _, singular, right = np.linalg.svd(white, full_matrices=False)
    return dict(
        escape_xyz=slopes[:, 0].mean(axis=0).tolist(),
        terminal_entropy_xyz=entropy.tolist(),
        weakest_instrument_axis=right[-1].tolist(),
        weighted_singular_values=singular.tolist(),
        bulk_xyz=[0.0, 0.0, float(-np.arctanh(REFERENCE[2]))],
    )


def primary(counts, meter, regional):
    """Aggregate independent gateway outcomes at the final appointment, retaining quiet."""
    p = counts[:, -1] / SHOTS
    if meter == "escape":
        value = p[:, :, 0]
        variance = value * (1 - value) / SHOTS
    else:
        value = entr(p).sum(axis=-1)
        score = np.zeros_like(p)
        np.log(p, out=score, where=p > 0)
        variance = (np.sum(p * score**2, axis=-1) - value**2) / SHOTS
    mean = value.mean(axis=1)
    gradients = np.diff(mean.reshape(3, 2), axis=1).ravel() / (2 * STEP)
    errors = np.sqrt(variance.sum(axis=1).reshape(3, 2).sum(axis=1)) / (3 * 2 * STEP)
    bulk_slope = regional["slopes"][2]
    eta = -bulk_slope / gradients[2] if abs(gradients[2]) > 5 * errors[2] else None
    return dict(
        slopes_xyz=gradients.tolist(),
        standard_errors=errors.tolist(),
        bulk_population_central_difference=bulk_slope,
        calibrated_coefficient=eta,
        held_residuals_xy=None
        if eta is None
        else (eta * gradients[:2] + np.asarray(regional["slopes"][:2])).tolist(),
    )


def spectral_inverse(calibration, held):
    cal = calibration.astype(float) / SHOTS
    offset = cal.mean(axis=0)
    slopes = ((cal[1::2] - cal[::2]) / 2).reshape(3, -1).T
    included = offset.ravel() * SHOTS >= 100
    root = np.sqrt(offset.ravel()[included])
    normalized = slopes[included] / root[:, None]
    obs = (held / SHOTS - offset).reshape(len(held), -1)[:, included] / root
    parameters, _, _, singular = np.linalg.lstsq(normalized, obs.T, rcond=None)
    fit = parameters.T
    real_fit = np.linalg.lstsq(normalized[:, [0, 2]], obs.T, rcond=None)[0].T
    return dict(
        diagnostic_recovered_xz=real_fit.tolist(),
        diagnostic_maximum_xz_error=np.max(
            abs(real_fit - HELD[:, [0, 2]]), axis=0
        ).tolist(),
        retained_bins=int(included.sum()),
        weighted_singular_values=singular.tolist(),
        recovered_xyz=fit.tolist(),
        maximum_xyz_error=float(np.max(abs(fit - HELD))),
        maximum_entropy_error=max(
            abs(h_bloch(u) - h_bloch(v)) for u, v in zip(fit, HELD, strict=True)
        ),
    )


def channel_capacity(calibration):
    """The gateway is uniformly chosen, independent of the bit, and recorded as a setting."""
    branch = calibration[-2:, -1].astype(float) / SHOTS
    conditional_entropy = entr(branch).sum(axis=-1).mean(axis=-1)

    def information(prior):
        mixture = (1 - prior) * branch[0] + prior * branch[1]
        return float(
            entr(mixture).sum(axis=-1).mean()
            - np.dot([1 - prior, prior], conditional_entropy)
        )

    result = minimize_scalar(
        lambda prior: -information(prior), bounds=(0.0, 1.0), method="bounded"
    )
    return dict(
        maximizing_prior=float(result.x), channel_capacity_nats=-float(result.fun)
    )


def phase_challenge(response, rng, norm_bound):
    targets = PHASE_STATES[1:]
    counts = sample(response, targets, rng)
    observed_tv = np.abs(counts[1] - counts[0]).sum(axis=-1) / (2 * SHOTS)
    predicted_tv = STEP * np.abs(response[:, :, :, 2]).sum(axis=-1)
    # TV is at most pure-state norm error, and the same bound carries through mixtures.
    full_bound = predicted_tv + 2 * norm_bound
    return dict(
        held_counts=counts.tolist(),
        sampled_total_variations=observed_tv.tolist(),
        compressed_total_variations=predicted_tv.tolist(),
        maximum_full_tv_bound=float(full_bound.max()),
        complete_experiment_tv_bound=min(1.0, float(full_bound.sum() * SHOTS)),
    )


def measure(level):
    board = federation(level)
    dark, launches = inputs(board)
    first, second = product_ball(board, launches)
    kinetic, agreement, keys = operators(board, first, second)
    packet = seed_packet(keys, board.n, dark, launches)
    bins = [portals(first, second, j) for j in launches]
    regional = regional_read(
        dark, np.vstack([REFERENCE, HELD]), np.random.default_rng([SEED, level, 4])
    )
    phase_regional = regional_read(
        dark, PHASE_STATES, np.random.default_rng([SEED, level, 5])
    )
    summary = dict(
        regional=regional,
        phase_regional=phase_regional,
        cells=board.cells,
        port_pairs_full=board.n * (board.n - 1) // 2,
        port_pairs_used=len(keys),
        receiver_ports=launches.tolist(),
        compression_depth=DEPTH,
        state_norm_error_bound=2 * 3.0**DEPTH / math.factorial(DEPTH),
        studies={},
    )
    for name, gain in (("registered", 1.0), ("half_exchange", 0.5)):
        response = evolve(gain * kinetic + agreement, packet, bins)
        if level == 0 and name == "registered":
            full = PairSector(board)
            dense_keys = full.aa * board.n + full.bb
            original = evolve(
                full.hamiltonian(candidate.KAPPA),
                seed_packet(dense_keys, board.n, dark, launches),
                [portals(full.aa, full.bb, port) for port in launches],
            )
            discrepancy = float(np.abs(original - response).max())
            assert discrepancy < 2e-12
            summary["full_board_maximum_effect_discrepancy"] = discrepancy
        rng = np.random.default_rng([SEED, level, int(gain * 2)])
        calibration = sample(response, CALIBRATION, rng)
        reference = sample(response, REFERENCE[None], rng)
        held = sample(response, HELD, rng)
        summary["studies"][name] = dict(
            calibration_counts=calibration.tolist(),
            reference_counts=reference.tolist(),
            held_counts=held.tolist(),
            escape=primary(held, "escape", regional),
            terminal_entropy=primary(held, "entropy", regional),
            derivatives=noiseless_derivatives(response),
            phase_challenge=phase_challenge(
                response,
                np.random.default_rng([SEED, level, 8, int(gain * 2)]),
                summary["state_norm_error_bound"],
            ),
            inverse=spectral_inverse(calibration, held),
            channel=channel_capacity(calibration),
        )
    return summary


def run():
    return dict(
        seed=SEED,
        shots_per_ensemble=SHOTS,
        exposures=TIMES,
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        phase_challenge_sha256=hashlib.sha256(
            (HERE / "PHASE_CHALLENGE.md").read_bytes()
        ).hexdigest(),
        phase_challenge_bloch=PHASE_STATES.tolist(),
        central_bloch=REFERENCE.tolist(),
        signed_axis_calibrations=CALIBRATION.tolist(),
        held_bloch=HELD.tolist(),
        boards=[measure(level) for level in range(3)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.coarse_response import (
            check_replay,
        )

        check_replay(
            json.loads((HERE / "results.json").read_text()),
            json.loads(json.dumps(result)),
        )
    for board in result["boards"]:
        rows = board["studies"]["registered"]
        print(
            board["cells"],
            rows["escape"],
            rows["terminal_entropy"],
            rows["inverse"]["maximum_xyz_error"],
            rows["channel"],
        )
