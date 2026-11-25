"""Directed weak first-click bins recover the fixed-marginal source correlations.

Two trace-one head states with the same one-head marginals differ by
sum_j c_j P_j/4 for the nine products P_j = X,Y,Z tensor X,Y,Z. The
four-tick pump is linear. A bin probability per launch is the trace of the
post-pump state against its first-relay effect. Subtracting the reference
therefore gives a linear calibration M c, valid even for finite preparations.

Only the calibration axes enter the numerical pivot selection from the 144
candidate bins. The frozen nine scalar probabilities are the minimum number
that can span a nine-dimensional tangent space. Integer no-click propagation
independently certifies invertibility, both for the appointment instrument
and for the continuous absorber. The difference between survival operators
at a bin's endpoints is its unconditional effect in either instrument.

The retained modular covector is k_j = Tr[-log(rho_ref) pump(P_j/4)].
Solving M.T w=k gives a weighted sum of probabilities for its response.
Uniform trace and operator bounds control the approximation to these weights.
The count intervals use separate Bernoulli preparations, every launch included.
Geometric lengths, proper volume and stationarity are separate observables.
"""

from functools import lru_cache
from math import expm1

import numpy as np
from scipy.linalg import qr

from .. import hidden, live_cadence, live_radar, live_source
from . import flight


DIVISOR = 16
AXES = tuple(i + j for i in "XYZ" for j in "XYZ")
NEAR = ((0, 1), (1, 2), (2, 3))


@lru_cache(maxsize=1)
def directions():
    inputs = np.stack([np.array(a, dtype=complex) for a in hidden.correlation_axes()])
    return inputs, np.array([live_source.propagate_heads(x, True) for x in inputs])


def effects(charge, edges, divisor=DIVISOR):
    """Unconditional first-relay bin effects from the independent renewal series."""
    kern = flight.free(charge)[:: max(flight.DIVISORS) // divisor]
    table = []
    for a, b in edges:
        f = live_cadence.first_hits(kern[:, b, a], kern[:, b, b], -expm1(-1 / divisor))
        table.append(live_radar.gram(f)[1:].reshape(12, divisor // 2, 8, 8).sum(axis=1))
    return np.stack(table)


def menu(divisor=DIVISOR):
    ids = []
    obs = []
    for charge in (-1, 1):
        edges = sorted({pair for a, b in NEAR for pair in ((a, b), (b, a))})
        vals = effects(charge, edges, divisor)
        for edge, block in zip(edges, vals, strict=True):
            for bin_id, e in enumerate(block):
                ids.append((charge, *edge, bin_id))
                obs.append(e)
    return ids, np.stack(obs)


def matrix(obs):
    return np.einsum("kij,zji->kz", obs, directions()[1]).real


def select(m):
    """Propose nine calibration bins by column-pivoted QR on the transposed menu."""
    _, _, pivot = qr(m.T, pivoting=True)
    return pivot[:9]


# Frozen numerical pivot order, addressed by charge, sender index, relay index
# and the zero-based half-time bin. Selection uses only calibration derivatives.
SELECTED = (
    (1, 2, 1, 4),
    (-1, 0, 1, 1),
    (-1, 3, 2, 1),
    (-1, 1, 0, 0),
    (-1, 2, 3, 0),
    (-1, 2, 1, 5),
    (-1, 1, 2, 8),
    (-1, 1, 0, 7),
    (-1, 2, 3, 10),
)
UNIT = 1 << 50
DEGREE = 22


def gram(w):
    from ..source import certificate as cert

    g = cert.exact_gram(w)
    return g[:8, :8] + g[8:, 8:], g[:8, 8:] - g[8:, :8]


def contraction():
    """Dyadic no-click attenuation certified against the real exponential."""
    from fractions import Fraction as F
    from math import factorial

    x = F(-1, 2 * DIVISOR)
    lo = sum(x**k / factorial(k) for k in range(16))
    hi = lo + x**16 / factorial(16)
    rounded = round((lo + hi) * UNIT / 2)
    assert lo < hi and max(abs(q - F(rounded, UNIT)) for q in (lo, hi)) < F(1, UNIT)
    return rounded


def col_error(h, continuous=False):
    from fractions import Fraction as F
    from math import factorial, isqrt

    L = F(int(abs(h).sum(axis=1).max()) + 20 * continuous, 40 * DIVISOR)
    assert L <= F(3, 2)
    # At each Taylor recurrence, real and imaginary integer limbs are each
    # rounded within 1/(2 UNIT) entrywise. Partial sums are exact integers.
    rounding = F(1 + isqrt(2 * h.shape[0]), 2 * UNIT)
    err = total = 0
    for j in range(1, DEGREE + 1):
        err = L * err / j + rounding
        total += err
    tail = L ** (DEGREE + 1) / factorial(DEGREE + 1)
    # Exact input has norm one; as long as cumulative column error <1, the
    # numerical input norm is <=2, its freely propagated center norm <=3.
    # Only 8 doubled-complex entries incur rounding in the absorbing update.
    result = flight.DEADLINE * DIVISOR * (2 * tail + total + F(5, UNIT))
    assert result < 1
    return result


def bounded_effects(settings, continuous=False):
    """Exact integer survival Gram differences for the independently propagated clock bins."""
    from itertools import groupby
    from fractions import Fraction as F

    data = {}
    damp = contraction()
    for (charge, a, b), items in groupby(sorted(settings), key=lambda row: row[:3]):
        bins = [it[3] for it in items]
        ticks = {j * (DIVISOR // 2) for t in bins for j in (t, t + 1)}
        h, rows = flight.hardware(charge)
        norm = int(abs(h).sum(axis=1).max()) + 20 * continuous
        w = np.zeros((h.shape[0], 16), dtype=np.int64)
        w[rows[a], :8] = np.eye(8, dtype=np.int64) * UNIT
        snap = {0: gram(w)}
        for t in range(1, max(ticks) + 1):
            term = w.copy()
            for j in range(1, DEGREE + 1):
                den = 40 * DIVISOR * j
                assert norm * int(abs(term).max()) + den < 2**62
                v = h @ term
                re, im = v[:, 8:], -v[:, :8]
                if continuous:
                    re[rows[b]] -= 20 * term[rows[b], :8]
                    im[rows[b]] -= 20 * term[rows[b], 8:]
                term[:, :8] = (re + den // 2) // den
                term[:, 8:] = (im + den // 2) // den
                assert int(abs(w).max()) + int(abs(term).max()) < 2**62
                w += term
            if not continuous:
                w[rows[b]] = (w[rows[b]].astype(object) * damp + UNIT // 2) // UNIT
            if t in ticks:
                snap[t] = gram(w)
        delta = 3 * col_error(
            h, continuous
        )  # sqrt(8) bounds the bank of eight column errors.
        effect_error = 2 * (2 * delta + delta**2)
        for k in bins:
            start, end = (snap[j * (DIVISOR // 2)] for j in (k, k + 1))
            data[(charge, a, b, k)] = (
                tuple(x - y for x, y in zip(start, end, strict=True)),
                F(effect_error),
            )
    return data


@lru_cache(maxsize=1)
def exact_pump():
    import sympy as sp

    _, _, colerr, grams = live_source.pump_center(True)
    outputs = [
        live_source.rational_read(grams, ax, (2, 1)) for ax in hidden.correlation_axes()
    ]
    return outputs, sp.Rational(4 * colerr + 4 * colerr**2)


def parts(x):
    """Real Hilbert-Schmidt coordinates for exact Hermitian pairings."""
    import sympy as sp

    r, i = x.as_real_imag()
    return sp.Matrix([*r, *i])


@lru_cache(maxsize=1)
def exact_columns():
    import sympy as sp

    outs, _ = exact_pump()
    return sp.Matrix.hstack(*(parts(x) for x in outs))


def certified_response(settings, fixed):
    import sympy as sp

    outs, pump_error = exact_pump()
    mids = []
    error = sp.Rational(0)
    for s in settings:
        (re, im), effect_error = fixed[s]
        E = sp.Matrix([[*re.ravel().tolist(), *im.ravel().tolist()]])
        mids.append(list(E * exact_columns() / UNIT**2))
        # Exact E and true channel contract. Each signed Pauli axis has trace norm one.
        error = max(
            error,
            sp.Rational(effect_error) + (1 + sp.Rational(effect_error)) * pump_error,
        )
    return sp.Matrix(mids), error, outs, pump_error


def rounded_matrix(M):
    import sympy as sp

    return M.applyfunc(lambda x: sp.Rational(round(x * 10**14), 10**14))


def inverse_bound(M, delta):
    """Rational lower bound on the true least singular value via inverse norms."""
    inv = M.inv()
    inf = max(sum(abs(v) for v in inv.row(i)) for i in range(inv.rows))
    one = max(sum(abs(v) for v in inv.col(i)) for i in range(inv.cols))
    # sqrt(one*inf) <= (one+inf)/2, and ||entry error||_2 <= 9 delta.
    lower = 2 / (one + inf) - 9 * delta
    if lower <= 0:
        raise ValueError("bin errors can erase a source direction")
    return lower


def setting(row):
    q, a, b, j = row
    return dict(
        charge=q,
        sender=flight.STATIONS[a],
        relay=flight.STATIONS[b],
        start=j / 2,
        end=(j + 1) / 2,
    )


@lru_cache(maxsize=1)
def centered_front_operators():
    """Integer operators for the six independently spanning reciprocal fronts."""
    import sympy as sp
    from . import observability, threshold

    source = live_radar.source
    board = source.formation.federation(1)
    contacts = threshold.hop(board, source.formation.region(board), -1)
    ops = []
    for d, a, b in observability.SELECTED:
        block = threshold.blocks(contacts, np.array([a]), np.array([b]), d)[0]
        op = sp.Matrix(threshold.effect(block, d)[0])
        centered = 8 * op - sp.trace(op) * sp.eye(8)
        ops.append(centered * sp.sign(next(v for v in centered if v)))
    return tuple(ops)


@lru_cache(maxsize=1)
def front_response_matrix():
    """Exact source-center rows of the six rationally normalized front effects."""
    import sympy as sp
    from . import observability

    ops = centered_front_operators()
    return (
        sp.Matrix.vstack(*(parts(op).T / observability.norm_inf(op) for op in ops))
        * exact_columns()
    )


def supplement(late_mids, entry_error, pump_error):
    """Select and certify three late bins that span the blind front directions."""
    from math import nextafter, inf
    from . import observability

    fronts = front_response_matrix()
    rounded = rounded_matrix(fronts)
    numerical = np.array(rounded, float)
    blind = np.linalg.svd(numerical, full_matrices=True)[2][6:].T
    unseen = np.array(late_mids, float) @ blind
    _, _, pivot = qr(unseen.T, pivoting=True)
    choose = pivot[:3].tolist()
    combined = rounded.col_join(late_mids[choose, :])
    error = max(entry_error, pump_error + max(abs(z) for z in fronts - rounded))
    lower = inverse_bound(combined, error)
    return dict(
        front_links=[list(x) for x in observability.SELECTED],
        front_rows_are_centered_and_infinity_norm_normalized=True,
        added_settings=[setting(SELECTED[i]) for i in choose],
        added_response_in_blind_subspace_singular_values=np.linalg.svd(
            unseen[choose], compute_uv=False
        ).tolist(),
        mixed_rank=9,
        mixed_least_singular_value_lower=nextafter(float(lower), -inf),
    )


def held(ops, calibration, weights, model_residual):
    from .. import live_survey
    from ..source import measure as source
    import json

    board = source.formation.federation(1)
    formed, product = source.form(board, source.formation.region(board), -1)
    gs = live_survey.rotations(formed)

    def turn(side, angle):
        return np.cos(angle / 2) * np.eye(4) - 2j * np.sin(angle / 2) * gs[side]

    both = turn("left", np.pi / 2) @ turn("right", -np.pi / 2)
    states = {
        "reference": formed,
        "decorrelated": product,
        "opposite_quarter_turns": both @ formed @ both.conj().T,
        "unmarked_lottery": (
            product + turn("right", np.pi) @ formed @ turn("right", np.pi).conj().T
        )
        / 2,
    }
    post = {k: live_source.propagate_heads(rho) for k, rho in states.items()}
    probabilities = {k: np.einsum("kij,ji->k", ops, s).real for k, s in post.items()}
    log_data = json.loads((flight.HERE.parent / "live_modular.json").read_text())
    log = np.array(log_data["reference_log_real"]) + 1j * np.array(
        log_data["reference_log_imag"]
    )
    rates = 1 << 54
    # Independent Bernoulli acquisitions for each held preparation and bin.
    # Each unselected first-click time and every timeout is the other outcome.
    rng = np.random.Generator(np.random.PCG64(19418))
    counts = {
        k: [
            sum(
                int(z) for z in rng.binomial(1 << 45, float(p), size=rates // (1 << 45))
            )
            for p in ps
        ]
        for k, ps in probabilities.items()
    }
    half = np.sqrt(np.log(2 * len(counts) * len(ops) / 0.01) / (2 * rates))
    total_half = 2 * half * np.linalg.norm(weights, 1) + model_residual
    results = {}

    def entropy(rho):
        values = np.linalg.eigvalsh(rho)
        return -np.dot(values, np.log(values))

    for k, state in states.items():
        if k == "reference":
            continue
        dr = state - formed
        axes = directions()[0]
        coeff = np.einsum("aij,ji->a", 4 * axes, dr).real
        p = probabilities[k] - probabilities["reference"]
        recovery = np.linalg.solve(calibration, p)
        trial = post[k]
        s_delta = float(np.trace(log @ (trial - post["reference"])).real)
        dcount = (np.array(counts[k], dtype=float) - counts["reference"]) / rates
        pred = float(weights @ dcount)
        rebuilt = post["reference"] + np.einsum("k,kij->ij", recovery, directions()[1])
        residual = max(np.linalg.norm(m) for m in source.formation.marginals(dr))
        results[k] = dict(
            input_minimum_eigenvalue=float(np.linalg.eigvalsh(state)[0]),
            maximum_marginal_change=float(residual),
            true_input_coordinates=coeff.tolist(),
            maximum_coordinate_recovery_error=float(max(abs(coeff - recovery))),
            direct_retained_modular_change=s_delta,
            reconstructed_retained_modular_change=float(weights @ p),
            actual_entropy_change=float(entropy(trial) - entropy(post["reference"])),
            reconstructed_entropy_change=float(
                entropy(rebuilt) - entropy(post["reference"])
            ),
            counted_modular_interval=[pred - total_half, pred + total_half],
        )
    return dict(
        sampling_scope="Synthetic independent Bernoulli counts at the nine nominal bin probabilities. Each first click outside the bin and each timeout is counted in the complementary outcome.",
        positive_preparations=results,
        probabilities={k: v.tolist() for k, v in probabilities.items()},
        counts=counts,
        shots_per_bin=rates,
        simultaneous_confidence=0.99,
        count_mean_half_width=float(half),
        modular_half_width=float(total_half),
    )


def front_blind_read(result):
    """Detect the independently registered blind direction and its native lotteries."""
    import json
    from fractions import Fraction as F
    from math import nextafter, inf

    witness = json.loads((flight.HERE / "observability.json").read_text())
    c = list(map(F, witness["blind_direction_center"]))
    # The witness seed is dyadic; allow one extra float ulp in calibration export.
    correction = F(witness["maximum_coordinate_correction_upper"]) + F(1, 10**15)
    entry = F(result["entry_error_upper"]) + F(1, 10**15)
    row_index = SELECTED.index((-1, 2, 1, 5))
    m = list(map(F, result["response_matrix"][row_index]))
    center = sum(x * y for x, y in zip(m, c, strict=True))
    radius = sum(map(abs, m)) * correction + entry * (sum(map(abs, c)) + 9 * correction)
    scale = witness["native_lotteries"]["common_scale"]

    def bounds(s):
        return [
            nextafter(float((center - radius) * s), -inf),
            nextafter(float((center + radius) * s), inf),
        ]

    return dict(
        setting=setting(SELECTED[row_index]),
        blind_input_probability_slope=bounds(1),
        exact_native_lottery_probability_difference=bounds(F(1, scale)),
        lottery_scale=scale,
    )


def fine_counts(result):
    """Resolve the small quarter-turn contrast with an independent count budget."""
    n = 1 << 64
    name = "opposite_quarter_turns"
    probabilities = result["held"]["probabilities"]
    rng = np.random.Generator(np.random.PCG64(918733))
    counts = {
        k: [
            sum(int(z) for z in rng.binomial(1 << 49, float(p), size=n // (1 << 49)))
            for p in probabilities[k]
        ]
        for k in ("reference", name)
    }
    # 18 bounded Bernoulli settings, separate reference/held preparations.
    r = float(np.sqrt(np.log(2 * 18 / 0.01) / (2 * n))) * (1 + 1e-12)
    dp = np.array(
        [
            float(x - y) / n
            for x, y in zip(counts[name], counts["reference"], strict=True)
        ]
    )
    m = np.array(result["response_matrix"])
    weights = np.array(result["log_weights"])
    u = np.linalg.solve(m, dp)
    # Any proposed source-coordinate vector u is valid here, whether obtained
    # from an approximate solve or otherwise. Include its actual linear residual.
    residual = np.linalg.norm(dp - m @ u) + 1e-16
    error = 6 * r + 9 * result["entry_error_upper"] * np.linalg.norm(u) + residual
    source_norm = min(
        9.0,
        np.linalg.norm(u, 1)
        + 3 * error / result["certified_least_singular_value_lower"]
        + 1e-12,
    )
    norm_bias = source_norm * result["modular_coefficient_residual_upper"]
    uncertainty = 2 * r * result["weight_l1"] + norm_bias + 1e-12
    pred = float(weights @ dp)
    return dict(
        held=name,
        shots_per_bin=n,
        simultaneous_confidence=0.99,
        counts=counts,
        mean_half_width=r,
        coefficient_l1_upper_from_counts=float(source_norm),
        modular_bias_upper_from_counts=float(norm_bias),
        counted_modular_interval=[pred - uncertainty, pred + uncertainty],
    )


def run():
    import json
    import sympy as sp
    from math import nextafter, inf
    from ..source import certificate as cert

    ids, ops = menu()
    data = matrix(ops)
    indices = [ids.index(s) for s in SELECTED]
    fixed = bounded_effects(SELECTED)
    exact, error, outs, pump_err = certified_response(SELECTED, fixed)
    nominal = rounded_matrix(exact)
    error += max(abs(z) for z in exact - nominal)
    lower = inverse_bound(nominal, error)
    m = np.array(nominal, dtype=float)
    log_record = json.loads((flight.HERE.parent / "live_modular.json").read_text())
    log = np.array(log_record["reference_log_real"]) + 1j * np.array(
        log_record["reference_log_imag"]
    )
    log_q = cert.prior.sym_matrix(log)
    k = (parts(log_q).T * exact_columns()).T
    k0 = rounded_matrix(k)
    kerr = (
        sp.Rational(float(log_record["reference_log_op_error"]))
        + 7 * pump_err
        + max(abs(z) for z in k - k0)
    )
    w = nominal.T.inv() * k0
    residual = kerr + sum(abs(z) for z in w) * error
    coeff = np.array(w, dtype=float).ravel()

    # Any two normalized states have HS distance <=sqrt(2); the chosen Pauli/4
    # basis gives ||delta coefficient||_1 <=3*2sqrt(2)<9 on fixed marginals.
    def upper(x):
        return nextafter(float(x), inf)

    def lower_float(x):
        return nextafter(float(x), -inf)

    continuous_exact, ce, _, _ = certified_response(
        SELECTED, bounded_effects(SELECTED, continuous=True)
    )
    cm = rounded_matrix(continuous_exact)
    ce += max(abs(z) for z in continuous_exact - cm)
    clower = inverse_bound(cm, ce)
    csv = np.linalg.svd(np.array(cm, float), compute_uv=False)
    sv = np.linalg.svd(m, compute_uv=False)
    condition = sv[0] / sv[-1]
    result = dict(
        divisor=DIVISOR,
        gamma=1,
        deadline=6,
        scope="Fixed two one-head marginals, 2:1 rail lottery, and native four-tick vacancy pump. The identified covector is the retained reference modular Hamiltonian.",
        selection="rank-revealing QR on a 144-bin calibration menu, nine distinct scalar rows",
        axes=list(AXES),
        settings=[setting(s) for s in SELECTED],
        menu_singular_values=np.linalg.svd(data, compute_uv=False).tolist(),
        selected_singular_values=sv.tolist(),
        selected_condition=float(condition),
        rank=9,
        continuous_monitor=dict(
            method="Loss of the full unclicked semigroup norm between the same physical times; no time quadrature",
            rank=9,
            response_matrix=np.array(cm, float).tolist(),
            singular_values=csv.tolist(),
            condition=float(csv[0] / csv[-1]),
            entry_error_upper=upper(ce),
            least_singular_value_lower=lower_float(clower),
        ),
        certified_least_singular_value_lower=lower_float(lower),
        entry_error_upper=upper(error),
        response_matrix=m.tolist(),
        independent_first_hit_versus_survival_max_error=float(
            np.max(abs(data[indices] - m))
        ),
        three_bin_supplement=supplement(nominal, error, pump_err),
        log_weights=coeff.tolist(),
        weight_l1=float(np.linalg.norm(coeff, 1)),
        weight_l2=float(np.linalg.norm(coeff)),
        modular_coefficient_residual_upper=upper(residual),
        uniform_fixed_marginal_modular_bias_upper=upper(9 * residual),
        nominal_modular_covector=np.array(k0, dtype=float).ravel().tolist(),
        held=held(ops[indices], m, coeff, upper(9 * residual)),
    )
    result["front_blind_read"] = front_blind_read(result)
    result["fine_count_check"] = fine_counts(result)
    return result


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    s = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        Path(__file__).with_suffix(".json").write_text(s)
    else:
        print(s, end="")
