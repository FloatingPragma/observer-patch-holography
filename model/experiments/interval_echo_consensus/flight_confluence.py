"""A permanent occupied-address fork with fixed priority and the same free flight."""

import argparse
from fractions import Fraction
from math import ceil, factorial, floor
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    doubled_hamiltonian,
)
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from .__main__ import local_ops
from .asynchronous import SharpService, priority
from .coherent_sweep import service_column
from .strict_pair import verify_columns


def winner_map(board):
    incident = [[] for _ in range(board.n)]
    for i, (a, b) in enumerate(board.seams):
        incident[a].append(i)
        incident[b].append(i)
    return np.array(
        [min(edges, key=lambda j: priority(board.seams[j])) for edges in incident]
    )


def native_double():
    exchange, contrast = native_law.local_term_parts()
    columns = np.eye(len(exchange))[
        :, [native_law.SRC.index(w) for w in ((2, 0), (0, 2))]
    ]
    swap = np.array([[0.0, 1.0], [1.0, 0.0]])
    ids = [native_law.SRC.index(w) for w in ((2, 0), (0, 2))]
    # Exact integer write columns. Conjugating these identities by the
    # native echo proves both generators preserve their entire range.
    written = native_law._write_pattern()[:, ids]
    reverse, absolute_squared = native_law._label_ops()
    np.testing.assert_array_equal(
        np.kron(reverse, np.eye(9, dtype=int)) @ written, written[:, ::-1]
    )
    np.testing.assert_array_equal(
        np.kron(absolute_squared, np.eye(9, dtype=int)) @ written, 4 * written
    )
    np.testing.assert_allclose(
        exchange @ columns, columns @ ((np.eye(2) - swap) / 2), rtol=0, atol=5e-14
    )
    np.testing.assert_allclose(contrast @ columns, 4 * columns, rtol=0, atol=5e-14)
    # An isolated pi pulse implements SWAP on these two blank inputs. The
    # squared-contrast coefficient contributes the phase exp(-i 32 pi)=1.
    pulse = expm(-1j * np.pi * (exchange + KAPPA * contrast))
    hit = local_ops()[2] @ columns
    np.testing.assert_allclose(pulse @ columns, columns @ swap, rtol=0, atol=8e-14)
    np.testing.assert_allclose(pulse @ hit, hit, rtol=0, atol=8e-14)
    return float(np.linalg.norm(pulse @ columns - columns @ swap, ord=2))


def paths(board, all_ports):
    """Integer walk powers on the intact native doubled-load sector."""
    winner = winner_map(board)
    twice = board.laplacian().astype(np.int64) - 6 * sparse.eye(
        board.n, dtype=np.int64, format="csr"
    )
    assert np.max(abs(twice).sum(axis=1)) == 6
    # A doubled load costs four contrast units at every incident
    # contact and at each missing degree-reference contact.
    endpoints = np.array(board.seams).ravel()
    contrast = 4 * np.bincount(endpoints, minlength=board.n) + 4 * (6 - board.degree)
    np.testing.assert_array_equal(contrast, np.full(board.n, 24))
    assert set(board.degree) == {5, 6}
    for strength in (0, KAPPA):
        physical = doubled_hamiltonian(board, strength)
        error = 2 * physical - (24 * strength + 6) * sparse.eye(board.n) - twice
        assert not error.count_nonzero()
    source = np.arange(board.n) if all_ports else np.array([0])
    mate = []
    for x in source:
        e = winner[x]
        a, b = board.seams[e]
        y = a + b - x
        mate.append(y if winner[y] == e else x)
    mate = np.array(mate)
    bad = sum(
        (np.count_nonzero(winner[np.array(board.neighbours[x])] != winner[x]) < 4)
        for x in range(board.n)
    )
    assert bad == 0
    v = np.zeros((board.n, len(source)), dtype=np.int64)
    cols = np.arange(len(source))
    v[source, cols] = 1
    powers = []
    degree = 22
    assert 6**degree < 2**63
    for j in range(degree + 1):
        powers.append((v[source, cols].copy(), v[mate, cols].copy()))
        if j < degree:
            v = twice @ v
    results = {}
    for time in (Fraction(1, 4), Fraction(1)):
        weights = [
            (-1) ** (k // 2) * (time / 2) ** k / factorial(k) for k in range(degree + 1)
        ]
        error = (3 * time) ** (degree + 1) / factorial(degree + 1)
        probability_error = (2 + error) * error
        ranges = []
        for j, x in enumerate(source):
            stay = Fraction(0)
            for side in range(1 + int(mate[j] != x)):
                for parity in (0, 1):
                    component = sum(
                        weights[k] * int(powers[k][side][j])
                        for k in range(parity, degree + 1, 2)
                    )
                    stay += component**2
            ranges.append((1 - stay - probability_error, 1 - stay + probability_error))

        def enclose(bounds):
            a, b = bounds
            return [floor(a * 10**12) / 10**12, ceil(b * 10**12) / 10**12]

        results[str(time)] = dict(
            launch_zero_other_address_probability=enclose(ranges[0]),
            all_tested_minimum=enclose(tuple(map(min, zip(*ranges)))),
            all_tested_maximum=enclose(tuple(map(max, zip(*ranges)))),
        )
    return dict(
        ports=int(board.n),
        source_columns=len(source),
        winning_addresses=len(np.unique(winner)),
        times=results,
    )


def guard_effects(service, winner):
    sites = np.arange(len(winner))
    ready = np.array([np.isin(sites, edge) for edge in service.edges])
    effects = []
    for i, blocked_by in enumerate(service.higher):
        permit = ready[i].copy()
        for j in blocked_by:
            permit &= ~ready[j]
        np.testing.assert_array_equal(permit, winner == i)
        effects.append(permit)
    effects = np.array(effects)
    np.testing.assert_array_equal(effects.sum(axis=0), np.ones(len(sites), dtype=int))
    return effects


def first_attempts(board, service, winner, h, rng, shots):
    # Rebuild all guarded hit effects using the readiness of competing
    # addresses. An independent fair edge is issued after free flight.
    guard = guard_effects(service, winner)
    trials = []
    zero = np.eye(1, board.n)[0]
    for index, (t, phase) in enumerate(((0.0, 0), (0.25, 1), (1.0, 1j))):
        wave = zero.copy().astype(complex)
        if index:
            wave[board.neighbours[0][0]] = phase
            wave /= np.linalg.norm(wave)
        state = expm_multiply(-1j * t * h, wave)
        conditional = guard @ abs(state) ** 2
        joint = conditional / len(guard)
        np.testing.assert_allclose(sum(joint), 1 / len(guard), rtol=0, atol=3e-16)
        hits = rng.multinomial(shots, np.r_[joint, 1 - joint.sum()])
        trials.append(dict(time=t, hit=int(hits[:-1].sum()), quiet=int(hits[-1])))
    return dict(
        atomic_guards=len(guard),
        identity_effect_diagonal=1,
        exact_one_offer_probability=str(Fraction(1, len(guard))),
        expected_offers=len(guard),
        simulated_first_offers=trials,
    )


def run():
    # Rebuild all 99-state exact guarded-write and binary flight columns,
    # including both signed-memory parity sectors on any occupied edge.
    verify_columns()
    pulse_error = native_double()
    board = federation(0)
    winner = winner_map(board)
    service = SharpService(range(board.n), board.seams)
    order = sorted(range(len(board.seams)), key=service.rank.__getitem__)
    for x in range(board.n):
        source = np.zeros(board.n, dtype=int)
        source[x] = 2
        record, _, _ = service_column(service, order, source, [0] * len(board.seams))
        assert [j for j, outcome in enumerate(record) if outcome == "write"] == [
            winner[x]
        ]
    x = 0
    y = next(q for q in board.neighbours[x] if winner[q] != winner[x])
    counts = []
    h = doubled_hamiltonian(board, KAPPA) - 12 * KAPPA * sparse.eye(board.n)
    initial = np.zeros(board.n)
    initial[x] = 1
    # Change only the exchange coefficient to zero during free flight.
    # The squared-contrast and vacancy-contact terms retain the address
    # symmetry, while the same strict repair and clock still run.
    degree_native = np.bincount(np.array(board.seams).ravel(), minlength=board.n)
    unchanged_contrast = KAPPA * 4 * (degree_native + 6 - board.degree)
    changed = expm_multiply(-1j * sparse.diags(unchanged_contrast), initial)
    changed_other = float(np.sum(abs(changed[winner != winner[x]]) ** 2))
    assert changed_other == 0.0
    rng = np.random.default_rng(715123)
    shots = 2**20
    for time in (0.25, 1.0):
        arrived = expm_multiply(-1j * time * h, initial)
        probs = np.bincount(
            winner, weights=abs(arrived) ** 2, minlength=len(board.seams)
        )
        np.testing.assert_allclose(probs.sum(), 1.0, atol=1e-13)
        draws = rng.multinomial(shots, probs / probs.sum())
        counts.append(
            dict(
                time=time,
                pre_sweep_occupied_each_seam=draws.tolist(),
                other_address_count=int(shots - draws[winner[x]]),
                other_address_probability=float(1 - probs[winner[x]]),
                post_sweep_other_address_count=0,
                post_sweep_home_count=shots,
            )
        )
    return dict(
        random_guard=first_attempts(board, service, winner, h, rng, shots),
        all_sphere_short_time_lower="529/16384",
        launch=x,
        installed_winner=list(board.seams[winner[x]]),
        isolated_pi_pulse=dict(
            edge=[x, int(y)],
            other_winner=list(board.seams[winner[y]]),
            output_address_total_variation=1,
        ),
        exchange_off_other_address_probability=changed_other,
        raw_contact_pi_pulse_error=pulse_error,
        zero_contrast_same_blank_double_walk=True,
        native_panels=[paths(federation(level), level == 0) for level in range(3)],
        shots_per_schedule=shots,
        complete_first_address_counts=counts,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(report)
    else:
        print(report, end="")
