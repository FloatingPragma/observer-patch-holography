"""Density fluctuation and driven transport on the unchanged blank binary Interval Echo boards."""

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import cell, federation

from ..measure import clock_seam, spectrum

HERE = Path(__file__).resolve().parent
BASE = 1 / 4
ORDINARY = 1 / 2
MIX = (1 / 8, 3 / 8)
APPOINTMENTS = (1 / 32, 1 / 8, 1 / 4)


def lowest_profile(board):
    """Project port zero onto the first positive Laplacian eigenspace."""
    lap = board.laplacian()
    if board.n == 12:
        values, modes = np.linalg.eigh(lap.toarray())
    else:
        values, modes = eigsh(
            lap,
            k=8,
            sigma=-1e-5,
            which="LM",
            v0=np.sin(0.53 * np.arange(board.n) + 0.1),
            tol=1e-13,
        )
    keep = abs(values - values[1]) < 1e-9
    chosen = modes[:, keep] @ modes[0, keep]
    chosen -= chosen.mean()
    chosen /= np.linalg.norm(chosen)
    rayleigh = float(chosen @ (lap @ chosen))
    return chosen, dict(
        eigenvalue=rayleigh,
        multiplicity=int(keep.sum()),
        next_eigenvalue=float(values[np.flatnonzero(keep)[-1] + 1]),
        residual=float(np.linalg.norm(lap @ chosen - rayleigh * chosen)),
    )


def registered_profile(board):
    payload = (HERE / "profiles.json").read_bytes()
    z = json.loads(payload)[str(board.n)]
    desired, spectral_read = lowest_profile(board)
    norm2 = sum(x * x for x in z)
    g = np.asarray(z, dtype=float) / np.sqrt(float(norm2))
    assert sum(z) == 0 and np.linalg.norm(desired - g) < 2e-6
    lap_z = [
        len(neighbours) * z[i] - sum(z[j] for j in neighbours)
        for i, neighbours in enumerate(board.neighbours)
    ]
    rayleigh = Fraction(sum(x * y for x, y in zip(z, lap_z)), norm2)
    resid_sq = Fraction(sum(y * y for y in lap_z), norm2) - rayleigh**2
    return (
        z,
        g,
        dict(
            sha256=hashlib.sha256(payload).hexdigest(),
            exact_rayleigh=str(rayleigh),
            exact_residual_squared=str(resid_sq),
            selector=spectral_read,
        ),
    )


def swap_amplitude(neighbours, i, j, g, kappa):
    """Four times the integer-profile swap entry, with i filled and j empty."""
    a, b = neighbours[i], neighbours[j]
    common = a & b
    coefficient = {}
    const = 0
    if j in a:
        gradient = g[i] - g[j]
        const += gradient * (len(a) - len(b))
        for k in (a ^ b) - {i, j}:
            coefficient[k] = 2 * gradient * (1 + 2 * kappa) * (1 if k in b else -1)
    for k in common:
        term = g[i] + g[j] - 2 * g[k]
        const += term
        coefficient[k] = -2 * term
    return const, coefficient


def coefficients(board, z, kappa):
    """Exact squared two-commutator for the registered profile at every filling."""
    norm2 = sum(v * v for v in z)
    lap_z = [
        len(ns) * z[i] - sum(z[j] for j in ns) for i, ns in enumerate(board.neighbours)
    ]
    second = Fraction(sum(x * y for x, y in zip(z, lap_z)), 2 * norm2)
    polynomial = [4 * sum(x * x for x in lap_z), 0, 0]
    neighbours = [set(v) for v in board.neighbours]
    pairs = set(board.seams)
    for ns in board.neighbours:
        pairs.update(combinations(ns, 2))
    for i, j in sorted(pairs):
        constant, coeff = swap_amplitude(neighbours, i, j, z, kappa)
        linear = sum(coeff.values())
        squares = sum(v * v for v in coeff.values())
        polynomial[0] += 2 * constant**2
        polynomial[1] += 2 * (2 * constant * linear + squares)
        polynomial[2] += 2 * (linear**2 - squares)
    _, coeff = swap_amplitude(neighbours, *clock_seam(), z, kappa)
    # A reset on either endpoint erases the signed coherent seam. On a
    # spectator, it differentiates that Bernoulli bit, picking its coefficient.
    local = Fraction(sum(z[k] * v for k, v in coeff.items()), 2 * norm2)
    return second, [Fraction(p, 16 * norm2) for p in polynomial], local


def mixture_read(fillings, second, polynomial):
    f = [Fraction(x) for x in fillings]
    chi = [x * (1 - x) for x in f]
    fourth = [polynomial[0] + x * (polynomial[1] + x * polynomial[2]) for x in f]
    return dict(
        susceptibility=sum(chi) / len(f),
        passive_frequency_squared=sum(c * v for c, v in zip(chi, fourth))
        / (sum(chi) * second),
        driven_frequency_squared=sum(fourth) / (len(f) * second),
    )


def same_source_tangents(second, polynomial):
    """Ordinary filling and fixed-filling correlation at one faithful equilibrium."""
    f, spread = Fraction(1, 4), Fraction(1, 8)
    x = (f, f - spread, f + spread)
    weights = (Fraction(1, 2), Fraction(1, 4), Fraction(1, 4))
    weight_slopes = (-1, Fraction(1, 2), Fraction(1, 2))
    chi = [y * (1 - y) for y in x]
    chi_prime = [1 - 2 * y for y in x]
    a, b, c = [v / second for v in polynomial]
    q = [a + y * (b + c * y) for y in x]
    q_prime = [b + 2 * c * y for y in x]

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v, strict=True))

    susceptibility = dot(weights, chi)
    passive = sum(w * z * v for w, z, v in zip(weights, chi, q)) / susceptibility
    ordinary_passive = (
        sum(
            w * (z * dv + dz * (v - passive))
            for w, z, dv, dz, v in zip(weights, chi, q_prime, chi_prime, q)
        )
        / susceptibility
    )
    held_passive = (
        sum(dw * z * (v - passive) for dw, z, v in zip(weight_slopes, chi, q))
        / susceptibility
    )
    ordinary_driven, held_driven = dot(weights, q_prime), dot(weight_slopes, q)
    expected = ordinary_driven * held_passive / ordinary_passive
    values = dict(
        ordinary_passive=ordinary_passive,
        ordinary_driven=ordinary_driven,
        fixed_marginal_passive=held_passive,
        fixed_marginal_driven=held_driven,
        calibrated_prediction=expected,
        observed_over_predicted=held_driven / expected,
    )
    return {name: str(value) for name, value in values.items()}


def finite_appointment(board, kappa, poly):
    """Outward rational signal bands with explicit finite preparations and dwell."""
    width = (1 + kappa) * len(board.seams) + kappa * (
        6 * board.n - 2 * len(board.seams)
    )
    moments = {}
    for name, menu in source_menus().items():
        w = [Fraction(weight) for weight, _ in menu]
        x = [Fraction(f) for _, f in menu]
        chi = [v * (1 - v) for v in x]
        m4 = [poly[0] + v * (poly[1] + v * poly[2]) for v in x]
        moments[name] = {
            "driven": sum(a * b for a, b in zip(w, m4)),
            "passive": sum(a * c * b for a, c, b in zip(w, chi, m4))
            / sum(a * c for a, c in zip(w, chi)),
        }
    measures = {}
    for arm in ("ordinary", "collar"):
        for meter in ("passive", "driven"):
            low, high = [
                moments[arm + "_" + side][meter] for side in ("lower", "upper")
            ]
            measures[arm + "_" + meter] = (high - low, width**2 * (high + low))
    # The order-six remainder is <= t^6 * W^2 (m4_+ + m4_-)/720.
    # Require it below half the absolute t^4/24 contribution in every arm.
    if any(delta == 0 for delta, _ in measures.values()):
        raise ValueError(
            "a vanishing fourth-order difference does not give a signed appointment"
        )
    t = Fraction(1)
    while any(2 * t * t * tail >= 30 * abs(delta) for delta, tail in measures.values()):
        t /= 2
    bounds = {}
    for key, (delta, tail) in measures.items():
        center = delta * t**4 / 24
        radius = tail * t**6 / 720
        bounds[key] = [str(center - radius), str(center + radius)]
    return dict(
        dwell=str(t), hamiltonian_width_upper=width, absolute_signal_differences=bounds
    )


def compare_board(board):
    z, _, calibration = registered_profile(board)
    output = {}
    for kappa in (int(KAPPA), 0):
        second, polynomial, local = coefficients(board, z, kappa)
        sources = {
            name: mixture_read(f, second, polynomial)
            for name, f in (
                ("reference", [BASE]),
                ("ordinary", [ORDINARY]),
                ("fixed_marginals", MIX),
            )
        }

        def change(arm, meter):
            return sources[arm][meter] - sources["reference"][meter]

        passive, driven = "passive_frequency_squared", "driven_frequency_squared"
        comparisons = dict(
            ordinary_passive_change=change("ordinary", passive),
            ordinary_driven_change=change("ordinary", driven),
            fixed_marginal_passive_change=change("fixed_marginals", passive),
            fixed_marginal_driven_change=change("fixed_marginals", driven),
            localized_seam_driven_second_time_derivative=local,
            localized_seam_passive_change=Fraction(0),
        )
        output[str(kappa)] = dict(
            mixed_reference_tangents=same_source_tangents(second, polynomial),
            finite_appointment=finite_appointment(board, kappa, polynomial),
            sources={
                name: {k: float(v) for k, v in row.items()}
                for name, row in sources.items()
            },
            normalized_second_moment=str(second),
            normalized_fourth_polynomial=list(map(str, polynomial)),
            comparisons={k: str(v) for k, v in comparisons.items()},
            readings={k: float(v) for k, v in comparisons.items()},
            baseline_fourth_over_second_squared=float(
                sources["reference"][passive] / second
            ),
        )
    return dict(ports=board.n, profile=calibration, by_contrast=output)


def source_menus():
    """(weight, filling) preparations for independently scheduled pulse and clock trials."""

    def mixed(f, e):
        return [(1 - e, f), (e / 2, f - 1 / 8), (e / 2, f + 1 / 8)]

    return dict(
        reference=[(1, BASE)],
        ordinary=[(1, ORDINARY)],
        fixed_marginals=[(0.5, MIX[0]), (0.5, MIX[1])],
        common=mixed(BASE, 0.5),
        ordinary_lower=mixed(BASE - 1 / 32, 0.5),
        ordinary_upper=mixed(BASE + 1 / 32, 0.5),
        collar_lower=mixed(BASE, 0.25),
        collar_upper=mixed(BASE, 0.75),
    )


def modular_read():
    """Sharp counts on the endpoints of a seam in the common stationary source."""
    from math import log, sqrt

    distribution = sum(
        w * np.kron([1 - f, f], [1 - f, f]) for w, f in source_menus()["common"]
    )
    odds = Fraction(float(distribution[1] * distribution[2])) / Fraction(
        float(distribution[0] * distribution[3])
    )
    assert odds == Fraction(529, 657)
    trials = 2**30
    counts = np.random.default_rng(81702).multinomial(trials, distribution)
    p = counts / trials
    err = sqrt(log(2 * 4 / 0.01) / (2 * trials))
    lower = log((p[1] - err) * (p[2] - err) / ((p[0] + err) * (p[3] + err))) / 64
    upper = log((p[1] + err) * (p[2] + err) / ((p[0] - err) * (p[3] - err))) / 64
    return dict(
        probabilities=distribution.tolist(),
        source_counts=counts.tolist(),
        trials=trials,
        seed=81702,
        collar_modular_per_mixing_change="log(529/657)/64",
        signed_modular_slope=log(float(odds)) / 64,
        independent_99pct_band=[lower, upper],
    )


def cell_records(kappa):
    """Passive spectral weights and separate active emitter/receiver traces on all 4096 words."""
    board = cell()
    a, b = clock_seam()
    _, g, _ = registered_profile(board)
    menus = source_menus()
    fills = sorted({f for mixture in menus.values() for _, f in mixture})
    passive = np.zeros((len(fills), len(APPOINTMENTS)))
    driven = np.zeros_like(passive)
    checks = np.zeros((len(fills), 2))
    local = np.zeros(len(APPOINTMENTS))
    for number, (energies, vec, masks, index) in enumerate(spectrum(kappa)):
        bits = (masks[:, None] >> np.arange(board.n)) & 1
        observable = bits @ g
        coefficients = (vec.T @ (observable[:, None] * vec)) ** 2
        gaps = energies[:, None] - energies[None, :]
        powers = [np.sum(coefficients * gaps**p) for p in (2, 4)]
        weights = np.array([f**number * (1 - f) ** (board.n - number) for f in fills])
        checks += weights[:, None] * powers
        pulses = np.zeros((len(masks), len(fills)))
        arm = np.empty(0)
        for i, f in enumerate(fills):
            # Independent full and empty rail release at a randomly enrolled
            # port. The signed postprocessing g_i extracts the chosen mode.
            full = bits * f ** (number - 1) * (1 - f) ** (board.n - number)
            empty = (1 - bits) * f**number * (1 - f) ** (board.n - number - 1)
            pulses[:, i] = (full - empty) @ g
        left = np.flatnonzero((bits[:, a] == 1) & (bits[:, b] == 0))
        right = [index[int(masks[i]) ^ (1 << a) ^ (1 << b)] for i in left]
        spectators = np.array([j for j in range(board.n) if j not in (a, b)])
        if left.size:
            rest = bits[left][:, spectators]
            k = number - 1
            full = rest * BASE ** (k - 1) * (1 - BASE) ** (board.n - 2 - k)
            empty = (1 - rest) * BASE**k * (1 - BASE) ** (board.n - 3 - k)
            arm = -((full - empty) @ g[spectators])
        for i, t in enumerate(APPOINTMENTS):
            passive[:, i] += weights * np.sum(coefficients * np.cos(t * gaps))
            u = (vec * np.exp(-1j * t * energies)) @ vec.T
            detector = observable @ (abs(u) ** 2)
            driven[:, i] += detector @ pulses
            if left.size:
                offdiag = np.sum(
                    observable[:, None] * u[:, left].conj() * u[:, right], axis=0
                )
                local[i] += 2 * offdiag.real @ arm
    chi = np.array([f * (1 - f) for f in fills])
    sources = {}
    for name, spec in menus.items():
        ids = [fills.index(f) for _, f in spec]
        w = np.array([weight for weight, _ in spec])
        susceptibility = float(w @ chi[ids])
        sources[name] = dict(
            susceptibility=susceptibility,
            passive_normalized=(w @ passive[ids] / susceptibility).tolist(),
            driven=(w @ driven[ids]).tolist(),
        )
    return dict(
        times=list(APPOINTMENTS),
        material_clock_clicks=(np.sin(np.array(APPOINTMENTS) / 2) ** 2).tolist(),
        product_fillings=fills,
        direct_spectral_moments=(checks / chi[:, None]).tolist(),
        sources=sources,
        localized_seam_driven_change=local.tolist(),
    )


def binary_records(time_records, g, seed):
    """Compress complete shots using enrolled random ports and one binary score per shot."""
    from math import log, sqrt

    trials = 2**50
    probability_grid = 2**38
    local_step = 1 / 8
    scale = float(abs(g).sum() ** 2)
    rng = np.random.default_rng(seed)
    sources = time_records["sources"]
    rows = {}
    for name, panel in sources.items():
        rows[name] = (
            panel["susceptibility"],
            panel["passive_normalized"],
            panel["driven"],
        )
    ref = sources["reference"]
    for label, amplitude in (("local_plus", local_step), ("local_minus", -local_step)):
        moved = np.array(ref["driven"]) + amplitude * np.array(
            time_records["localized_seam_driven_change"]
        )
        rows[label] = (BASE * (1 - BASE), ref["passive_normalized"], moved)
    # Include the three independent dimer appointments in the same confidence family.
    settings = len(APPOINTMENTS) * (2 * len(rows) + 1)
    p_error = sqrt(log(2 * settings / 0.01) / (2 * trials)) + 0.5 / probability_grid

    def acquire(p):
        probability = round(p * probability_grid) / probability_grid
        return int(rng.binomial(trials, probability))

    counts = {}
    for name, (chi, passive, driven) in rows.items():
        counts[name] = {}
        for label, signals, factor in (
            ("passive", passive, 4 * chi),
            ("driven", driven, 1),
        ):
            raw = [acquire((1 + factor * y / scale) / 2) for y in signals]
            counts[name][label] = dict(
                positive_scores=raw,
                reconstructed=[(2 * x / trials - 1) * scale / factor for x in raw],
                simultaneous_99pct_halfwidth=2 * p_error * scale / factor,
            )
    return dict(
        trials_per_arm=trials,
        seed=seed,
        probability_grid=probability_grid,
        score_scale=scale,
        local_signed_step=local_step,
        settings=settings,
        dimer_positive_counts=[
            acquire(p) for p in time_records["material_clock_clicks"]
        ],
        click_probability_halfwidth=p_error,
        observations=counts,
    )


def run():
    _, g, _ = registered_profile(cell())
    cells = {}
    for kappa in (int(KAPPA), 0):
        data = cell_records(kappa)
        data["count_records"] = binary_records(data, g, 81724 + kappa)
        cells[str(kappa)] = data
    return dict(
        source="Interval Echo, blank conserved binary sector, all ports retained",
        active_contrast=int(KAPPA),
        local_seam=list(clock_seam()),
        ordinary_filling=ORDINARY,
        reference_filling=BASE,
        mixed_fillings=list(MIX),
        mixture_connected_covariance=1 / 64,
        common_stationary_seam_read=modular_read(),
        boards=[compare_board(b) for b in [cell(), *(federation(i) for i in range(3))]],
        cell_time_records=cells,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")
