"""Free flight of one ordered excess and vacancy in the pinned filled screen."""

from math import sqrt

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import free_seam
from .federation import propagate

TIMES = (0.0, 1.0, 3.0, 8.0, 16.0, 32.0)


def ordered_index(x, y, n):
    return x * (n - 1) + y - (y > x)


def generators(board):
    """All positions of distinguishable conserved excess (2) and vacancy (0)."""
    n = board.n
    x = np.repeat(np.arange(n), n - 1)
    offset = np.tile(np.arange(n - 1), n)
    y = offset + (offset >= x)
    ids = np.arange(n * (n - 1))
    rows, cols = [], []
    adj = np.zeros_like(x, dtype=bool)
    for moving_excess in (True, False):
        moving, parked = (x, y) if moving_excess else (y, x)
        for k in range(board.table.shape[1]):
            dest = board.table[moving, k]
            crossing = dest == parked
            if moving_excess:
                adj |= crossing
                rows.append(ordered_index(parked[crossing], moving[crossing], n))
                cols.append(ids[crossing])
            keep = (dest >= 0) & ~crossing
            u, v = (
                (dest[keep], parked[keep])
                if moving_excess
                else (parked[keep], dest[keep])
            )
            rows.append(ordered_index(u, v, n))
            cols.append(ids[keep])
    row, col = np.concatenate(rows), np.concatenate(cols)
    off = sparse.coo_matrix(
        (np.full(len(row), -0.5), (row, col)), shape=(len(ids), len(ids))
    )
    degrees = board.degree[x] + board.degree[y] - adj
    kinetic = (off + sparse.diags(degrees / 2)).tocsr()
    if (kinetic != kinetic.T).nnz or np.max(np.abs(kinetic.sum(axis=1))) != 0:
        raise ValueError("Incorrect ordered-defect exchange graph")
    # Diagonals relative to the uniform background minus the common 12*kappa.
    cost = 2 * KAPPA * (board.degree[y] - board.degree[x] + adj)
    h = (kinetic + sparse.diags(cost)).tocsr()
    seam = free_seam(KAPPA)
    for excess, vacancy in ((0, 1), (1, 0), board.external[0], (0, n - 1)):
        start = ordered_index(excess, vacancy, n)
        occupation = np.ones(n, dtype=int)
        occupation[excess], occupation[vacancy] = 2, 0
        direct = np.zeros(len(ids))
        direct[start] = KAPPA * ((6 - board.degree) @ (occupation**2 - 1)) - 12 * KAPPA
        for a, b in board.seams:
            left, right = occupation[a], occupation[b]
            i = law.SRC.index((left, right))
            direct[start] += seam[i, i]
            if left == right:
                continue
            swap = {a: b, b: a}
            dest = ordered_index(
                swap.get(excess, excess), swap.get(vacancy, vacancy), n
            )
            direct[dest] += seam[law.SRC.index((right, left)), i]
        np.testing.assert_allclose(
            h[:, start].toarray().ravel(), direct, atol=2e-13, rtol=0
        )
    return x, y, adj, kinetic, h, cost


def high_component(h, sharp, degree):
    """Legendre polynomial approximation of the upper spectral projector.

    With Z=H-30 the spectrum of Z² is in [2²,46²]. Mapping this interval
    to [-1,1] allows the generating series for Legendre polynomials to
    approximate (Z²)^(-1/2) with a geometric uniform error bound.
    """
    radius, gap = 46.0, 2.0
    ratio = (radius - gap) / (radius + gap)
    z = h - 30 * sparse.eye(h.shape[0], format="csr")

    def argument(v):
        return ((radius**2 + gap**2) * v - 2 * (z @ (z @ v))) / (radius**2 - gap**2)

    previous, current = sharp, argument(sharp)
    polynomial = sharp + ratio * current
    for j in range(1, degree):
        following = ((2 * j + 1) * argument(current) - j * previous) / (j + 1)
        polynomial += ratio ** (j + 1) * following
        previous, current = current, following
    remainder = radius / (2 * gap) * ratio ** (degree + 1)
    return sharp / 2 + z @ polynomial / (radius + gap), remainder


def address_coupling(board):
    """Row and column sums of the exchange across the high-address partition."""
    multiplicities = {}
    high = {
        (a, b) if board.degree[a] == 5 else (b, a)
        for a, b in board.seams
        if board.degree[a] != board.degree[b]
    }
    largest_out = 0
    for a, b in high:
        adjacent = [(u, b) for u in board.neighbours[a] if u != b]
        adjacent += [(a, v) for v in board.neighbours[b] if v != a]
        adjacent.append((b, a))
        crossing = [pair for pair in adjacent if pair not in high]
        largest_out = max(largest_out, len(crossing))
        for pair in crossing:
            multiplicities[pair] = multiplicities.get(pair, 0) + 1
    largest_in = max(multiplicities.values())
    coupling = sqrt(largest_out * largest_in) / 2
    # At every time, Q(H-37)psi has norm at most the initial standard
    # deviation sqrt(10)/2, and the upper edge of QHQ is at most 28.
    delta, sigma = 9.0, sqrt(10) / 2
    q = (delta * sigma + coupling * sqrt(delta**2 + coupling**2 - sigma**2)) / (
        delta**2 + coupling**2
    )
    return {
        "cells": board.cells,
        "upper_band_rank": len(high),
        "maximum_exits": largest_out,
        "maximum_entries": largest_in,
        "cross_block_norm_bound": coupling,
        "variance_binding_floor": 1 - q * q,
    }


def phase_read(h, source, x, y):
    """Reconstruct the complex launch amplitude from two local phase quadratures.

    A coherent local pulse puts the two-port state in an equal mixture of
    amplitudes for the undisturbed filled medium and the emitted ordered pair.
    Signed +/- outcome differences remove any incoherent returns to the
    filled local configuration while the defects reside elsewhere.
    """
    # The local native echo maps these written inputs to blank states with
    # either unchanged filled ports or excess and vacancy at them.
    echo = law.echo()
    blank = np.eye(law.DIM, law.N_SRC)
    controls = {(1, 1): 0, (2, 0): 2}
    for ports, contrast in controls.items():
        memory = np.zeros(1 + law.N_LAB)
        for coin in (0, 1):
            memory[1 + law.LABELS.index((coin, contrast))] = 1 / sqrt(2)
        written = np.kron(memory, np.eye(law.N_SRC)[:, law.SRC.index((1, 1))])
        np.testing.assert_allclose(
            echo @ written, blank[:, law.SRC.index(ports)], atol=2e-13, rtol=0
        )
    # This phase is coupled to the independently specified sum of native port
    # loads above the filled reference. A reference rail carries any compensating
    # unit for the single charged control arms.
    charge = np.array([a + b - 2 for a, b in law.SRC])
    angle = np.pi / 2
    response = {}
    for name, port_state in (
        ("excess_vacancy", (2, 0)),
        ("excess", (2, 1)),
        ("vacancy", (0, 1)),
    ):
        phase = np.exp(1j * angle * charge[law.SRC.index(port_state)])
        response[name] = [
            float(abs((1 + phase) / 2) ** 2),
            float(abs((1 - 1j * phase) / 2) ** 2),
        ]
    np.testing.assert_allclose(
        list(response.values()), [[1, 0.5], [0.5, 1], [0.5, 0]], rtol=0, atol=1e-14
    )
    index = int(np.argmax(source))
    full = h + 12 * KAPPA * sparse.eye(h.shape[0], format="csr")

    def read_at(t, shots, seed):
        state = propagate(full, source.astype(complex), t)
        amplitude = state[index]
        untouched = ~np.isin(x, [x[index], y[index]]) & ~np.isin(
            y, [x[index], y[index]]
        )
        # Partial trace on all ports outside the marked seam. Configurations
        # with neither defect here contribute to the filled-state diagonal.
        reduced = (
            np.array(
                [
                    [1 + np.sum(abs(state[untouched]) ** 2), amplitude.conjugate()],
                    [amplitude, abs(amplitude) ** 2],
                ]
            )
            / 2
        )
        quadratures = []
        for number, phase in enumerate((1, 1j)):
            local_modes = np.array([[1, 1], [phase, -phase]]) / np.sqrt(2)
            observed = np.diag(local_modes.conj().T @ reduced @ local_modes).real
            probabilities = np.r_[observed, 1 - observed.sum()]
            rng = np.random.default_rng([seed, number])
            counts = rng.multinomial(shots, probabilities)
            quadratures.append(
                {"probabilities": probabilities.tolist(), "counts": counts.tolist()}
            )
        signed = [r["probabilities"][0] - r["probabilities"][1] for r in quadratures]
        np.testing.assert_allclose(
            signed, [amplitude.real, amplitude.imag], atol=2e-13, rtol=0
        )
        measured = -signed[1] / t
        positive, negative, _ = quadratures[1]["counts"]
        sampled = (negative - positive) / shots / t
        se = np.sqrt((reduced.trace().real - signed[1] ** 2) / shots) / t
        norm_bound = 12 * KAPPA + 44
        error = norm_bound**3 * t * t / 6
        expected = float(full[index, index])
        if abs(measured - expected) > error:
            raise ValueError(
                "Filled-background energy differs from the local phase slope"
            )
        return {
            "exposure": t,
            "local_amplitude_real": float(amplitude.real),
            "local_amplitude_imag": float(amplitude.imag),
            "energy_relative_to_filled": expected,
            "finite_difference_energy": float(measured),
            "bias_bound": error,
            "phase_outcomes_plus_minus_quiet": quadratures,
            "sampled_energy": float(sampled),
            "sampling_standard_error": float(se),
        }

    first = read_at(1 / 4096, 2**26, 26112221)
    first["load_u1_xy_return_probabilities"] = response
    first["held_longer_exposure"] = read_at(1 / 2048, 2**26, 26112222)
    return first


def run():
    board = federation(0)
    x, y, adj, k, h, cost = generators(board)
    # Register the first gateway of cell 0 and the first regular neighbor.
    gateway, outside = min(e for e in board.external if e[0] // 12 == 0)
    regular = min(
        p for p in board.neighbours[gateway] if p // 12 == 0 and board.degree[p] == 5
    )
    initial = ordered_index(regular, gateway, board.n)
    source = np.zeros(len(x))
    source[initial] = 1
    high = cost == 4 * KAPPA
    if not high[initial] or np.max(cost[~high]) != 2 * KAPPA:
        raise ValueError("Source lies outside the registered split")
    mean = float(h[initial, initial])
    bound_k = float(max(np.asarray(abs(k).sum(axis=1)).flat))
    max_low = np.max(cost[~high]) + bound_k
    floor = (mean - max_low) / (float(np.max(cost)) - float(np.max(cost[~high])))
    projected, error = high_component(h, source, 260)
    lower, _ = high_component(h, source, 200)
    if np.linalg.norm(projected - lower) > 5e-7 or error > 2e-9:
        raise ValueError("Upper-band projection failed its precision check")
    weight = float(projected @ projected)
    lift = projected / sqrt(weight)
    # Only physical patch-relative occupations enter this menu. High addresses
    # are adjacent regular-excess and gateway-vacancy pairs.
    marked_port = y == gateway
    both_shifted = (y != gateway) & (x != regular)
    outside_cell = (x // 12 == outside // 12) & (y // 12 == outside // 12)
    in_x, in_y = x // 12 == 0, y // 12 == 0
    masks = {
        "any_adjacent": adj,
        "high_addresses": high,
        "reoriented_about_initial_vacancy": high & marked_port & (x != regular),
        "vacancy_changed_bound": high & ~marked_port,
        "both_addresses_changed_bound": high & both_shifted,
        "home_both_changed_bound": high & both_shifted & in_x & in_y,
        "receiver_cell_bound": high & outside_cell,
        "home_excess_without_vacancy": in_x & ~in_y,
        "home_vacancy_without_excess": in_y & ~in_x,
    }
    arms = {}
    for name, operator, prep in (
        ("native_sharp", h, source),
        ("upper_band", h, lift),
        ("exchange_only_sharp", k, source),
    ):
        psi = prep.astype(complex)
        prev = 0.0
        readings = []
        for appointment, t in enumerate(TIMES):
            psi = propagate(operator, psi, t - prev)
            prev = t
            weights = abs(psi) ** 2
            np.testing.assert_allclose(weights.sum(), 1, rtol=0, atol=7e-11)
            row = dict(
                exposure=t,
                **{key: float(weights[mask].sum()) for key, mask in masks.items()},
            )
            if name != "upper_band":
                local_bound = high & in_x & in_y
                bins = [
                    row["home_both_changed_bound"],
                    weights[local_bound & ~both_shifted].sum(),
                    weights[~local_bound].sum(),
                ]
                rng = np.random.default_rng(
                    [26112220, name == "exchange_only_sharp", appointment]
                )
                # The local number census may be quiet or leave the pair unordered.
                # Its three disjoint outcomes retain every sharp-start trial.
                draws = rng.random(131072)
                row["home_counts"] = np.bincount(
                    np.searchsorted(np.cumsum(bins)[:-1], draws), minlength=3
                ).tolist()
            readings.append(row)
        arms[name] = readings
    rank = [address_coupling(federation(level)) for level in (0, 1, 2)]
    return {
        "ports": board.n,
        "states": len(source),
        "source": [regular, gateway],
        "receiver_cell": int(outside // 12),
        "base_reference_offset": 12 * KAPPA,
        "shifted_spectral_ranges": [[-16, max_low], [32, 32 + bound_k]],
        "shifted_source_mean": mean,
        "upper_address_all_time_floor": float(floor),
        "upper_weight_of_sharp": weight,
        "upper_state_mean": float(lift @ (h @ lift)),
        "spectral_projector_degree": 260,
        "projector_truncation_norm_bound": error,
        "upper_band_ranks": rank,
        "local_energy_read": phase_read(h, source, x, y),
        "home_detector_bins": [
            "high_pair_both_changed",
            "other_home_high_pair",
            "quiet",
        ],
        "arms": arms,
    }
