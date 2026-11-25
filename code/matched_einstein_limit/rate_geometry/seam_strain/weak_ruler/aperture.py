"""Interior weak-capture census and outward-cut measurements on one pumped source."""

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from math import factorial, isqrt
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from .. import live_modular, live_radar, live_source, live_survey
from ..source import certificate as cert
from ..source import measure as source

HERE = Path(__file__).resolve().parent
DURATION = 6
RATE = 1
WIDTH = 8


def enroll(board, radius):
    region = set(source.formation.region(board)[1:3])
    for _ in range(radius):
        region |= {j for i in region for j in board.neighbours[i]}
    cut = sorted((i, j) for i in region for j in board.neighbours[i] if j not in region)
    rim = sorted({j for i, j in cut})
    return sorted(region), rim, cut


def calibrations():
    # Unit efficiency for an isolated occupied detector. Unit outward crossing
    # is a two-port native blank seam with its true off-diagonal -1/2.
    g = np.array([[0, 0.5j], [0.5j, -RATE / 2]], complex)
    w = expm(DURATION * g)[:, 0]
    return -np.expm1(-RATE * DURATION), 1 - abs(w[0]) ** 2


def launch(board, origins):
    rows = live_radar.station_rows(board, origins).ravel()
    start = np.zeros((8 * board.n, len(rows)), complex)
    start[rows, np.arange(len(rows))] = 1
    return start


def captures(board, h, targets, origins, keep):
    mask = np.zeros(h.shape[0])
    mask[live_radar.station_rows(board, targets).ravel()] = 1
    g = -1j * h / 40 - sparse.diags(mask * RATE / 2)
    initial = launch(board, origins)
    final = expm_multiply(DURATION * g, initial, traceA=DURATION * g.diagonal().sum())
    # All-particle no-click sum. Absolute survival uses every port and retained word.
    by_origin = final.reshape(len(final), len(origins), WIDTH).transpose(1, 0, 2)
    out = {}
    for name, rows in keep.items():
        selected = by_origin[:, rows, :]
        effects = np.eye(WIDTH) - selected.conj().transpose(0, 2, 1) @ selected
        out[name] = {p: (m + m.conj().T) / 2 for p, m in zip(origins, effects)}
    return out


GRID = 8
DEGREE = 26
UNIT = 2**50


def bounded_error(h):
    """Per-column Euclidean error of the real/imaginary Taylor solver."""
    row = int(abs(h).sum(axis=1).max()) + 20
    length = Fraction(row, 40 * GRID)
    assert length < 3 and max(np.diff(h.indptr)) <= 10
    # At most ten native products, one damping product, their additions,
    # and one division. gamma bounds each component by gamma |G| |x|.
    u = Fraction(1, 2**53)
    gamma = 16 * u / (1 - 16 * u)
    add = DEGREE * u / (1 - DEGREE * u)
    assert (1 + gamma) * length < 3
    # exp((1+gamma)L)<21, and exp((1+gamma)L)-exp(L)<=gamma L exp((1+gamma)L).
    rnd = 21 * (add + gamma * length)
    tail = length ** (DEGREE + 1) / factorial(DEGREE + 1)
    local = rnd + tail + Fraction(1, 10**290)
    total = DURATION * GRID * local
    assert total < Fraction(1, 10**5)
    # All exact steps contract: norms of computed inputs grow by at most 1+local.
    return total / (1 - total) + Fraction(1 + isqrt(2 * h.shape[0]), 2 * UNIT)


def bounded_capture(board, h, targets, origins, keep):
    target = live_radar.station_rows(board, targets).ravel()
    rows = live_radar.station_rows(board, origins).ravel()
    n = rows.size
    wave = np.zeros((h.shape[0], 2 * n))
    wave[rows, np.arange(n)] = 1
    real_h = h.astype(float)
    for _ in range(DURATION * GRID):
        term = wave.copy()
        for k in range(1, DEGREE + 1):
            v = real_h @ term
            real, imag = v[:, n:], v[:, :n]
            imag *= -1
            real[target] -= 20 * term[target, :n]
            imag[target] -= 20 * term[target, n:]
            den = 40 * GRID * k
            term[:, :n] = real / den
            term[:, n:] = imag / den
            wave += term
    assert np.max(abs(wave)) < 1.01
    ix = np.rint(wave * UNIT).astype(np.int64)
    out = {}
    for name, rows in keep.items():
        selected = ix[rows]
        out[name] = {}
        for j, p in enumerate(origins):
            indices = np.r_[j * 8 : j * 8 + 8, n + j * 8 : n + j * 8 + 8]
            g = cert.exact_gram(selected[:, indices])
            r = np.eye(8, dtype=object) * (UNIT**2) - g[:8, :8] - g[8:, 8:]
            i = g[8:, :8] - g[:8, 8:]
            out[name][p] = (
                np.asarray(r, dtype=float) + 1j * np.asarray(i, dtype=float)
            ) / float(UNIT**2)
    return out


def effect_error(h):
    # One orthogonal 8-column launch packet per emitter. Gram is integer
    # through subtraction of the no-click from the identity.
    d = 3 * bounded_error(h)
    return 2 * d + d * d + Fraction(16, 2**52)


def native_effects(charge, certified=False):
    board = source.formation.federation(1)
    h = source.system(board, source.formation.region(board), charge)[0]
    solver = bounded_capture if certified else captures
    pv, pa = calibrations()
    result = []
    for radius in (0, 1):
        inner, rim, cut = enroll(board, radius)
        census = sum(
            solve_groups(
                board, h, solver, inner, sorted(inner + rim), {"click": slice(None)}
            )["click"].values()
        )
        ops = {"census": census / pv}
        region = live_radar.station_rows(board, inner).ravel()
        exterior = np.ones(h.shape[0], dtype=bool)
        exterior[region] = False
        for tag, emitters, targets, keep in (
            ("outward", [p for p, q in cut], rim, region),
            ("inward", [q for p, q in cut], sorted({p for p, q in cut}), exterior),
        ):
            effects = solve_groups(
                board,
                h,
                solver,
                targets,
                sorted(set(emitters)),
                {"crossing": keep, "click": slice(None)},
            )
            ops[tag] = sum(effects["crossing"].values()) / pa
            ops[tag + "_click"] = sum(effects["click"].values()) / pa
        result.append((radius, inner, rim, cut, ops))
    return result


def solve_groups(board, h, solver, targets, emitters, keep):
    combined = {name: {} for name in keep}
    for start in range(0, len(emitters), 3):
        out = solver(board, h, targets, emitters[start : start + 3], keep)
        for name, table in out.items():
            combined[name].update(table)
    return combined


def expectation(op, state):
    return float(np.trace(op @ state).real)


def entropy(state):
    w = np.linalg.eigvalsh(state)
    assert w.min() > 0
    return float(-np.dot(w, np.log(w)))


def panels(effects):
    rho, deltas, cal = live_survey.retained(True)
    prepared = live_source.postpump()
    arms = {name: prepared[name] for name in ("ordinary", "decorrelated", "reference")}
    out = []
    for radius, inner, rim, cut, ops in effects:
        means = {
            arm: {k: expectation(m, state) for k, m in ops.items()}
            for arm, state in arms.items()
        }
        deriv = {
            term: {k: expectation(m, d) for k, m in ops.items()}
            for term, d in deltas.items()
        }
        low, middle, high = sorted(arms, key=lambda s: means[s]["census"])
        fraction = (means[middle]["census"] - means[low]["census"]) / (
            means[high]["census"] - means[low]["census"]
        )
        mix = (1 - fraction) * arms[low] + fraction * arms[high]
        difference = mix - arms[middle]
        val, vec = np.linalg.eigh(rho)
        k = -(vec * np.log(val)) @ vec.conj().T
        equal = dict(
            middle=middle,
            low=low,
            high=high,
            high_probability=fraction,
            lottery_minus_middle={
                tag: expectation(m, difference) for tag, m in ops.items()
            },
            retained_modular=expectation(k, difference),
            retained_entropy=entropy(mix) - entropy(arms[middle]),
        )
        stationarity = {}
        for side in ("left", "right"):
            dv = deriv[side]["census"]
            if abs(dv) < 1e-14:
                continue
            matched = {}
            for arm in ("ordinary", "decorrelated", "hidden"):
                angle = -deriv[arm]["census"] / dv
                matched[arm] = dict(
                    angle=angle,
                    **{
                        direction: deriv[arm][direction]
                        + angle * deriv[side][direction]
                        for direction in ("outward", "inward")
                    },
                    retained_modular=cal["retained"][arm]
                    + angle * cal["retained"][side],
                )
            predictions = {}
            for direction in ("outward", "inward"):
                eta = (
                    -matched["ordinary"]["retained_modular"]
                    / matched["ordinary"][direction]
                )
                predictions[direction] = dict(
                    ordinary_eta=eta,
                    held_residual={
                        arm: matched[arm]["retained_modular"]
                        + eta * matched[arm][direction]
                        for arm in ("decorrelated", "hidden")
                    },
                )
            stationarity[side] = dict(
                census_per_radian=dv, matched=matched, transfer=predictions
            )
        out.append(
            dict(
                radius=radius,
                interior=inner,
                rim=rim,
                oriented_cut=cut,
                means=means,
                slopes=deriv,
                finite_equal_census=equal,
                matched_census=stationarity,
            )
        )
    return out


def short_window(board, charge, inner, cut):
    h = source.system(board, source.formation.region(board), charge)[0]
    out = np.zeros((8, 8), dtype=np.int64)
    inc = out.copy()
    for p, q in cut:
        i, j = live_radar.station_rows(board, (p, q))
        block = h[j][:, i].toarray()
        out += block.T @ block
        inc += block @ block.T
    # Each registered blank crossing has rate coefficient 1/4; 40H is integer.
    return dict(
        census_count=len(inner),
        outward_coefficient_times_400=out.tolist(),
        inward_coefficient_times_400=inc.tolist(),
        universal_blank_cut=(
            np.array_equal(out, 400 * len(cut) * np.eye(8, dtype=int))
            and np.array_equal(inc, 400 * len(cut) * np.eye(8, dtype=int))
        ),
    )


def synthetic_counts(row, seed):
    rng = np.random.default_rng(seed)
    trials = 2**40
    pv, pa = calibrations()
    total = len(row["interior"]) + len(row["rim"])
    surface = {
        k: len({edge[i] for edge in row["oriented_cut"]})
        for i, k in enumerate(("outward", "inward"))
    }
    draws = {}
    means = row["means"].copy()
    match = row["finite_equal_census"]
    f = match["high_probability"]
    means["matched_lottery"] = {
        k: (1 - f) * means[match["low"]][k] + f * means[match["high"]][k]
        for k in means["reference"]
    }
    for arm, mean in means.items():
        draws[arm] = {
            "census_click": int(rng.binomial(trials, mean["census"] * pv / total))
        }
        for k in ("outward", "inward"):
            click = mean[k + "_click"] * pa / surface[k]
            terminal = (mean[k] - mean[k + "_click"]) * pa / surface[k]
            draws[arm][k] = rng.multinomial(
                trials, [click, terminal, 1 - click - terminal]
            ).tolist()
    # Fewer than 128 Bernoulli probabilities across all charge/radius/arm
    # settings and both event categories. Union-bound 99% sampling coverage.
    eps = np.sqrt(np.log(2 * 128 / 0.01) / (2 * trials))
    middle = draws[match["middle"]]
    lottery = draws["matched_lottery"]
    empirical = {}
    for key in ("census", "outward", "inward"):
        if key == "census":
            d = (lottery["census_click"] - middle["census_click"]) / trials
            scale = total / pv
        else:
            d = (sum(lottery[key][:2]) - sum(middle[key][:2])) / trials
            scale = surface[key] / pa
        empirical[key] = [float((d - 2 * eps) * scale), float((d + 2 * eps) * scale)]
    return dict(
        trials_per_arm_and_setting=trials,
        simultaneous_confidence=0.99,
        matched_lottery_high_weight=f,
        counted_lottery_minus_middle=empirical,
        seed=seed,
        channels=[
            "weak click",
            "terminal crossing without click",
            "source-side survivor",
        ],
        draws=draws,
    )


def rational_mean(op, state):
    return sum(
        Fraction(float(a.real)) * Fraction(float(b.real))
        - Fraction(float(a.imag)) * Fraction(float(b.imag))
        for a, b in zip(op.ravel(), state.T.ravel())
    )


def encode(v):
    return [
        float(np.nextafter(float(v.a), -np.inf)),
        float(np.nextafter(float(v.b), np.inf)),
    ]


@lru_cache(maxsize=1)
def material():
    # Both certificates check the floating source preparations of this run.
    return live_source.certify(), live_modular.certify()


def bands(row, ops, charge, source_data, modular):
    ctx = MPIntervalContext()
    ctx.dps = 40

    def ri(v):
        t = Fraction(v)
        return ctx.mpf(t.numerator) / t.denominator

    def extend(value, error):
        return value + ctx.mpf([-1, 1]) * error

    def pack(tree):
        if isinstance(tree, dict):
            return {k: pack(v) for k, v in tree.items()}
        return encode(tree)

    pv, pa = calibrations()
    vcal = 1 - ctx.exp(-DURATION)
    angle = ctx.sqrt(3) * DURATION / 4
    acal = (
        1
        - ctx.exp(-ctx.mpf(DURATION) / 2)
        * (ctx.cos(angle) + ctx.sin(angle) / ctx.sqrt(3)) ** 2
    )
    for numerical, physical in ((pv, vcal), (pa, acal)):
        assert abs(physical - ri(numerical)).b < ri(Fraction(1, 10**14))
    board = source.formation.federation(1)
    route = source.formation.region(board)
    op_error = effect_error(source.system(board, route, charge)[0])
    arm_errors = source_data["channels"]["coupled"]["nominal_state_trace_norm_errors"]
    tangent_errors = modular["nominal_tangent_trace_norm_errors"]
    ref, variations, _ = live_survey.retained(True)
    physical = live_source.postpump()
    matrices = {arm: physical[arm] for arm in ("reference", "ordinary", "decorrelated")}

    def raw(states, err, norm):
        out = {}
        for name, m in ops.items():
            census = name == "census"
            count = (
                (len(row["interior"]) + len(row["rim"]))
                if census
                else len(
                    {
                        e[0 if name.startswith("outward") else 1]
                        for e in row["oriented_cut"]
                    }
                )
            )
            cal = vcal if census else acal
            # <=22 positive contractions are added and calibrated in binary64.
            # Per summand, gamma_22 times eight matrix rows, one division and
            # the checked isolated calibration together contribute <2e-13.
            assert count <= 22 and cal.a > 0.9
            ee = count * (ri(op_error) + ri(Fraction(2, 10**13))) / cal
            out[name] = {}
            for arm, state in states.items():
                e = ri(err[arm])
                v = ri(rational_mean(m, state))
                out[name][arm] = extend(v, ee * (ri(norm[arm]) + e) + count * e / cal)
        return out

    means = raw(matrices, arm_errors, {k: 1 for k in matrices})
    norm = {
        "ordinary": 2,
        "decorrelated": 2,
        "hidden": 1,
        **{k: x[1] for k, x in modular["rotation_output_trace_norms"].items()},
    }
    slopes = raw(variations, tangent_errors, norm)
    kval = {k: ctx.mpf(v) for k, v in modular["modular_slopes"].items()}
    c = slopes["census"]

    # Work on differences relative to the common physical reference, cancelling
    # the large O(N) census intercept before a finite lottery comparison.
    def contrast(b, a, k):
        if b == a:
            return ctx.mpf(0)
        changes = slopes[k] if k != "K" else kval
        zero = ctx.mpf(0)
        return changes.get(b, zero) - changes.get(a, zero)

    e = row["finite_equal_census"]
    low, middle, high = (e[k] for k in ("low", "middle", "high"))
    frac = contrast(middle, low, "census") / contrast(high, low, "census")
    assert frac.a > 0 and frac.b < 1
    equal = {
        k: contrast(low, middle, k) + frac * contrast(high, low, k)
        for k in ("outward", "inward", "K")
    }
    head_centers = {}
    for key, state in matrices.items():
        m = cert.prior.sym_matrix(state)
        head_centers[key] = (m + m.H) / 2

    def entropy_interval(m, err):
        # Orthogonal-polar spectrum proposal with a checked Hermitian residue.
        _, w, _, radius = live_modular.spectral_center(m)
        er = ri(err) + 8 * ri(radius)
        assert ri(min(w.diagonal())) - er > ri(Fraction(1, 10**4))
        assert ri(max(w.diagonal())) + er < 2
        center = -sum(ri(v) * ctx.log(ri(v)) for v in w.diagonal())
        # On [1e-4,2], |(-x log x)'| <=9. The trace path/Hermitian
        # trace-norm bound controls finite entropy even when eigenvectors rotate.
        return extend(center, 9 * er)

    f0 = Fraction(float(e["high_probability"]))
    f_lo, f_hi = map(Fraction, encode(frac))
    uncertainty = max(abs(f0 - f_lo), abs(f0 - f_hi))
    mix = head_centers[low] * (1 - f0) + head_centers[high] * f0
    mix_error = max(Fraction(arm_errors[k]) for k in (low, high)) + 2 * uncertainty
    equal["retained_entropy"] = entropy_interval(mix, mix_error) - entropy_interval(
        head_centers[middle], arm_errors[middle]
    )
    transfers = {}
    for side in ("left", "right"):
        if c[side].a <= 0 <= c[side].b:
            continue
        matched = {}
        for k in ("ordinary", "decorrelated", "hidden"):
            angle = -c[k] / c[side]
            matched[k] = {
                "angle": angle,
                "K": kval[k] + angle * kval[side],
                **{
                    d: slopes[d][k] + angle * slopes[d][side]
                    for d in ("outward", "inward")
                },
            }
        transfer = {}
        for d in ("outward", "inward"):
            if matched["ordinary"][d].a <= 0 <= matched["ordinary"][d].b:
                continue
            eta = -matched["ordinary"]["K"] / matched["ordinary"][d]
            transfer[d] = {
                "eta": eta,
                **{
                    k: matched[k]["K"] + eta * matched[k][d]
                    for k in ("decorrelated", "hidden")
                },
            }
        transfers[side] = dict(matched=matched, transfer=transfer)
    return dict(
        per_origin_effect_op_error=float(np.nextafter(float(op_error), np.inf)),
        census_calibration=encode(vcal),
        transmission_calibration=encode(acal),
        means=pack(means),
        slopes=pack(slopes),
        high_lottery_weight=encode(frac),
        finite_equal_census=pack(equal),
        fixed_census=pack(transfers),
    )


def common_level_curve(slope_intervals, material_intervals):
    """Implicitly hold the three regional detectors fixed on native source controls."""
    outputs = ("census", "outward", "inward")
    adjustable = ("ordinary", "left", "right")

    def midwidth(bounds):
        lo, hi = map(sp.Rational, bounds)
        return ((lo + hi) / 2, (hi - lo) / 2)

    midpoint = sp.Matrix(
        [[midwidth(slope_intervals[k][j])[0] for j in adjustable] for k in outputs]
    )
    width = sp.Matrix(
        [[midwidth(slope_intervals[k][j])[1] for j in adjustable] for k in outputs]
    )
    b = sp.Matrix([midwidth(slope_intervals[k]["decorrelated"])[0] for k in outputs])
    berr = sp.Matrix([midwidth(slope_intervals[k]["decorrelated"])[1] for k in outputs])
    # Column scaling only changes coordinates on the adjustment parameters.
    scales = sp.diag(*[1 / max(abs(v) for v in midpoint[:, j]) for j in range(3)])
    inv = (midpoint * scales).inv()
    residual = inv.applyfunc(abs) * width * scales
    neumann = max(sum(residual.row(i)) for i in range(3))
    assert neumann < sp.Rational(1, 100)
    trial = -inv * b
    errors = inv.applyfunc(abs) * (berr + width * scales * trial.applyfunc(abs))
    budget = max(errors) / (1 - neumann)
    errors = scales * (errors + residual * sp.ones(3, 1) * budget)
    trial = scales * trial
    # With decorrelation probability s>=0, the ordinary probability must also
    # increase at the simplex vertex. Both native local angles can take either sign.
    assert trial[0] - errors[0] > 0
    center, err = midwidth(material_intervals["decorrelated"])
    for tag, t, e in zip(adjustable, trial, errors):
        h, w = midwidth(material_intervals[tag])
        center += h * t
        err += abs(h) * e + w * (abs(t) + e)
    assert abs(center) > err
    # For G with eigenvalues ±1/2, rho(+pi/2)-rho(-pi/2)=-2i[G,rho].
    # Splitting the kernel into two convex lotteries thus gives *finite*
    # sources with identical three detector populations, without taking angles
    # or source fractions to zero.
    assert center - err > 0
    positive = {
        "decorrelated": (sp.S.One, sp.S.Zero),
        "ordinary": (trial[0], errors[0]),
    }
    negative = {"reference": (1 + trial[0], errors[0])}
    for name, x, e in zip(adjustable[1:], trial[1:], errors[1:]):
        assert abs(x) > e
        sign = "+" if x > 0 else "-"
        positive[name + sign] = (abs(x) / 2, e / 2)
        negative[name + ("-" if sign == "+" else "+")] = (abs(x) / 2, e / 2)
    den = sum(x for x, e in positive.values())
    denerr = sum(e for x, e in positive.values())

    def ratio_band(x, e):
        lo, hi = (x - e) / (den + denerr), (x + e) / (den - denerr)
        return [cert.enclosed(lo, 0)[0], cert.enclosed(hi, 0)[1]]

    gap_lo = (center - err) / (den + denerr)
    gap_hi = (center + err) / (den - denerr)
    t = sp.Rational(1, 2**24)
    # Along each segment (1-t)rho_ref+tL, lambda_min >=(1-t)*.0014.
    # The second derivative of entropy lies in [-2/((1-t)*.0014),0].
    # The same one-sided interval holds for both chord remainders.
    entropy_radius = t * t / ((1 - t) * sp.Rational(14, 10000))
    s_lo = t * gap_lo - entropy_radius
    s_hi = t * gap_hi + entropy_radius
    assert s_lo > 0
    finite = dict(
        positive={k: ratio_band(x, e) for k, (x, e) in positive.items()},
        negative={k: ratio_band(x, e) for k, (x, e) in negative.items()},
        rotation_magnitude_radians="pi/2",
        modular_gap_nats=[cert.enclosed(gap_lo, 0)[0], cert.enclosed(gap_hi, 0)[1]],
        lottery_weight_in_common_reference=float(t),
        finite_entropy_gap_nats=[cert.enclosed(s_lo, 0)[0], cert.enclosed(s_hi, 0)[1]],
    )
    return dict(
        fixed=["census", "outward", "inward"],
        parameter="decorrelation probability at the formed reference",
        invertibility_error_upper=cert.enclosed(neumann, 0)[1],
        control_slopes={
            k: cert.enclosed(t, e) for k, t, e in zip(adjustable, trial, errors)
        },
        retained_entropy_slope_nats=cert.enclosed(center, err),
        detector_rank=3,
        detector_and_entropy_rank=4,
        finite_native_lotteries=finite,
    )


def aggregate(panel_sets):
    source_data, modular = material()
    out = {}
    for q, values in panel_sets.items():
        rows = panels(values)
        for i, (row, (_, _, _, _, ops)) in enumerate(zip(rows, values)):
            row["certificate"] = bands(row, ops, q, source_data, modular)
            row["three_equal_readings"] = common_level_curve(
                row["certificate"]["slopes"], modular["modular_slopes"]
            )
            row["short_window"] = short_window(
                source.formation.federation(1), q, row["interior"], row["oriented_cut"]
            )
            row["simulated_counts"] = synthetic_counts(row, 337409 + 2 * (q + 1) + i)
        out[str(q)] = rows
    vp, ap = calibrations()
    return dict(
        cells=80,
        forming_charge=-1,
        pump_charge=-1,
        pump_ticks=4,
        duration=DURATION,
        detection_rate=RATE,
        isolated_census_probability=vp,
        isolated_transmission_probability=ap,
        panels=out,
    )


def run():
    return aggregate({q: native_effects(q, certified=True) for q in (-1, 1)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "aperture.json").write_text(result)
    else:
        print(result, end="")
