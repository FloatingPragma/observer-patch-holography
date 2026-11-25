"""Two rest clocks, local charged orbital flight, and their joint readings on one source."""

from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import block_diag, expm

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging import (
    measure as ranging,
)
from ..full_algebra.native_precision import prior

HERE = Path(__file__).resolve().parent
EMISSION = ranging.formation.emission
X = np.array([[0, 1], [1, 0]])
IDENTITY = np.eye(2)
TIMES = np.pi / 3 * np.arange(5)
SHOTS = 2**24


def density_center():
    record = json.loads((HERE.parent / "full_algebra/results.json").read_text())
    radii = json.loads((HERE.parent / "full_algebra/certificate.json").read_text())
    center = np.array(record["x_real"]) + 1j * np.array(record["x_imag"])
    # Use the rationalized, normalized Hermitian center used by the equilibrium proof.
    sym = [prior.sym_matrix(block) for block in center]
    sym = [(block + block.H) / 2 for block in sym]
    norm = sum(s.trace() for s in sym)
    exact = [s / norm for s in sym]
    p = Fraction(exact[1].trace())
    coin_x = Fraction(
        sum(
            block[2 * a + b, 2 * (1 - a) + b].as_real_imag()[0]
            for block in exact
            for a in range(2)
            for b in range(2)
        )
    )
    radius = Fraction(str(radii["state_HS_error_upper"]))
    assert radius <= Fraction(2, 10**7)
    return np.array([np.array(s, dtype=complex) for s in exact]), p, coin_x, radius


def contacts(charge, k):
    """Leakage and exact action are checked directly on the 99 native columns."""
    native = EMISSION.four_h(k)
    law = EMISSION.LAW
    words = [(1, 1)] if charge == 0 else [(1 + charge, 1), (1, 1 + charge)]
    raw = np.zeros((law.DIM, 4 * len(words)), dtype=int)
    for idx, parity in enumerate((1, -1)):
        for pos, word in enumerate(words):
            for c in (0, 1):
                for sign in (-1, 1):
                    row = EMISSION.columns(word)[1 + law.LABELS.index((c, sign))]
                    raw[row, (idx * len(words) + pos) * 2 + c] = (
                        parity if sign == -1 else 1
                    )
    projected = raw.T @ native @ raw
    np.testing.assert_array_equal(2 * native @ raw, raw @ projected)
    np.testing.assert_array_equal(raw.T @ raw, 2 * np.eye(raw.shape[1], dtype=int))
    if charge == 0:
        bright = np.eye(2)
    else:
        bell = np.array([0, 1, 1, 0] if charge == 1 else [1, 0, 0, 1])
        bright = np.eye(4) - np.outer(bell, bell) / 2
    np.testing.assert_array_equal(
        projected, 8 * block_diag(k * bright, (k + 1) * bright)
    )
    return projected / 8


def blank_clock(charge, k):
    law = EMISSION.LAW
    columns = [law.SRC.index(word) for word in ((1 + charge, 1), (1, 1 + charge))]
    h = EMISSION.four_h(k)
    basis = np.eye(law.DIM, dtype=int)[:, columns]
    reduced = h[np.ix_(columns, columns)]
    np.testing.assert_array_equal(h @ basis, basis @ reduced)
    np.testing.assert_array_equal(reduced, 4 * k * IDENTITY + 2 * (IDENTITY - X))
    return reduced / 4


def marginal_blocks(state, side):
    if side == "left":
        return np.array([s.reshape(2, 2, 2, 2).trace(axis1=1, axis2=3) for s in state])
    right = state.sum(axis=0).reshape(2, 2, 2, 2).trace(axis1=0, axis2=2)
    return np.array([right, 0 * right])


def amplitudes(source, h, rail, times, start_clock):
    """Born table for reported endpoint and whether the two sign readings disagree."""
    # Basis parity x (rail x coin). Trace the distant head after the original preparation.
    rho = block_diag(*(np.kron(rail, b) for b in source))
    spatial = rail.shape[0]
    position = [np.kron(np.diag(np.eye(spatial)[i]), IDENTITY) for i in range(spatial)]
    signs = [(IDENTITY + s * X) / 2 for s in (-1, 1)]
    start = [np.kron(sign, np.eye(2 * spatial)) for sign in signs]
    result = np.zeros((len(times), 2, spatial))
    for n, t in enumerate(times):
        unitary = expm(-1j * t * h)
        if start_clock:
            for s, projection in enumerate(start):
                branch = unitary @ projection @ rho @ projection @ unitary.conj().T
                for e, sign in enumerate(signs):
                    for j, place in enumerate(position):
                        result[n, int(s != e), j] += np.trace(
                            np.kron(sign, place) @ branch
                        ).real
        else:
            final = unitary @ rho @ unitary.conj().T
            for j, place in enumerate(position):
                result[n, 0, j] = np.trace(np.kron(IDENTITY, place) @ final).real
    np.testing.assert_allclose(result.sum(axis=(1, 2)), 1, atol=2e-13, rtol=0)
    return result


def orbit(blocks, k, times):
    sharp = np.diag([1, 0])
    return (
        sum(
            amplitudes(blocks, contacts(q, k), sharp, times, False)[:, 0, 1]
            for q in (-1, 1)
        )
        / 2
    )


def recurrence(reads):
    y0, y1, y2, y3 = reads[:4]
    first = 2 * (y1 - y0)
    second = y2 - 2 * y1 + y0
    third = y3 - 2 * y2 + y1
    return first * (third + first) - 2 * second**2


def arbitrary_phase_defect(reads):
    """Eliminate the frequency and phase from consecutive first differences."""
    a, b, c, d = np.diff(reads)
    return b * (b + d) - c * (a + c)


def mixed_joint(blocks, times):
    sharp = np.diag([1, 0])
    return (
        sum(amplitudes(blocks, contacts(q, 8), sharp, times, True) for q in (-1, 1)) / 2
    )


def current_reads(blocks, orientation):
    rail = np.array([[1, -1j * orientation], [1j * orientation, 1]]) / 2
    endpoint = np.kron(IDENTITY, np.diag([0, 0, 1, 1]))
    source = block_diag(*(np.kron(rail, b) for b in blocks))
    # First sharp sign read, with the recorded sign retained for the clock comparison.
    swap = np.kron(X, np.eye(4))
    started = (source + swap @ source @ swap) / 2
    sign_effect = 0
    curves = 0
    side = (-1, 1)
    for q in side:
        h = contacts(q, 8)
        sign_effect += np.trace(1j * (h @ endpoint - endpoint @ h) @ started).real / 2
        curves += amplitudes(blocks, h, rail, TIMES, True)[:, 1, :].sum(axis=1) / 2
    # Fresh copies at short forward times give an operational endpoint-only current secant.
    step = 1e-5
    samples = 0
    for q in side:
        samples += (
            amplitudes(blocks, contacts(q, 8), rail, [0, step, 2 * step], True).sum(
                axis=1
            )[:, 1]
            / 2
        )
    derivative = (-3 * samples[0] + 4 * samples[1] - samples[2]) / (2 * step)
    return dict(
        orientation=orientation,
        current=float(sign_effect),
        endpoint_rate=float(derivative),
        all_flip=curves.tolist(),
    )


def clock_motion_spectra():
    """An additive moving clock requires a four-energy sumset in the native spectrum."""

    def multiplicities(h):
        twice = (2 * h).astype(int)
        np.testing.assert_array_equal(twice, 2 * h)
        return {
            str(level / 2): count
            for level, count in sp.Matrix(twice).eigenvals().items()
        }

    native = [multiplicities(contacts(q, 8)) for q in (-1, 1)]
    assert native[0] == native[1] == {"0": 2, "8": 3, "9": 3}
    blank = blank_clock(1, 8)
    independent = multiplicities(np.kron(blank, IDENTITY) + np.kron(IDENTITY, blank))
    assert independent == {"16": 1, "17": 2, "18": 1}
    return dict(native=native[0], two_independent_blank_dimers=independent)


def run():
    state, p, coin_x, eps = density_center()
    left = marginal_blocks(state, "left")
    right = marginal_blocks(state, "right")
    p_low, p_high = p - 2 * eps, p + 2 * eps
    assert Fraction(1, 5) < p_low < p_high < Fraction(1, 2)
    lower = Fraction(3, 8) * p_low * (1 - p_low)
    x_lower = coin_x - 3 * eps
    assert lower > Fraction("0.0605") and 17 * x_lower / 8 > Fraction("0.05637")
    rest = amplitudes(left, contacts(0, 8), np.ones((1, 1)), TIMES, True)[:, 1, 0]
    blank = [abs(expm(-1j * t * blank_clock(1, 8))[1, 0]) ** 2 for t in TIMES]
    np.testing.assert_allclose(rest, blank, rtol=0, atol=2e-14)
    rows = {}
    for k in (8, 0):
        left_row = orbit(left, k, TIMES)
        right_row = orbit(right, k, TIMES)
        expected = (
            (1 - float(p)) * (1 - np.cos(k * TIMES))
            + float(p) * (1 - np.cos((k + 1) * TIMES))
        ) / 4
        np.testing.assert_allclose(left_row, expected, atol=5e-14, rtol=0)
        rows[k] = dict(
            left=left_row.tolist(),
            right=right_row.tolist(),
            left_defect=recurrence(left_row),
            right_defect=recurrence(right_row),
            left_arbitrary_phase=arbitrary_phase_defect(left_row),
            right_arbitrary_phase=arbitrary_phase_defect(right_row),
        )
    rng = np.random.default_rng(197023)
    samples = [int(rng.binomial(SHOTS, v)) for v in rows[8]["left"]]
    scored = recurrence(np.array(samples) / SHOTS)
    error = 0.0006
    phase_scored = arbitrary_phase_defect(np.array(samples) / SHOTS)
    assert 10 * np.exp(-2 * SHOTS * error**2) < 0.01
    assert scored - 64 * error > 0 and phase_scored - 16 * error > 0
    joints = mixed_joint(left, TIMES)
    clocks = 0.75 * np.sin(TIMES / 2) ** 2
    np.testing.assert_allclose(joints[:, 1, :].sum(axis=1), clocks, atol=2e-14, rtol=0)
    np.testing.assert_allclose(joints[:, 1, 1], clocks / 6, atol=2e-14, rtol=0)
    travel = [current_reads(left, v) for v in (-1, 1)]
    for r in travel:
        np.testing.assert_allclose(
            r["current"], r["orientation"] * 17 * float(coin_x) / 8, atol=1e-13, rtol=0
        )
        np.testing.assert_allclose(r["endpoint_rate"], r["current"], atol=1e-8, rtol=0)
        np.testing.assert_allclose(r["all_flip"], clocks, atol=2e-14, rtol=0)
    return dict(
        appointments_in_pi_over_3=list(range(len(TIMES))),
        left_odd_center=float(p),
        left_coin_X_center=float(coin_x),
        enclosed_odd=[float(p_low), float(p_high)],
        exact_equilibrium_defect_lower="0.0605",
        exact_equilibrium_abs_current_lower="0.05637",
        exact_equilibrium_arbitrary_phase_defect_lower="0.03025",
        exact_level_multiplicities=clock_motion_spectra(),
        blank_dimer=blank,
        parked_head=rest.tolist(),
        orbital=rows,
        independent_simulated_endpoint_counts=dict(
            shots_per_endpoint=SHOTS,
            far=samples,
            near=[SHOTS - k for k in samples],
            recurrence_defect=float(scored),
            arbitrary_phase_defect=float(phase_scored),
            arbitrary_phase_99pct_lower=float(phase_scored - 16 * error),
            simultaneous_99pct_lower=float(scored - 64 * error),
        ),
        mobile_joint_same_flip_by_home_away=joints.tolist(),
        directed=travel,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
