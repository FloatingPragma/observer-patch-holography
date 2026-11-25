"""Full reduced load state of all four radar stations on the independent material controls."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv
from scipy import sparse

from ...action_response.four_ticks import parse_interval
from . import measure

c = measure.c
R = sp.Rational
SIGNS = iv.mpf([-1, 1])


def controls():
    archive = json.loads((measure.HERE / "energy_odd.json").read_text())
    arms = {
        "original_equal_work": json.loads((measure.HERE / "control.json").read_text())
    }
    arms.update({k: v for k, v in archive.items() if isinstance(v, dict)})
    arms["fixed_all_number_weights"] = json.loads(
        (measure.HERE / "station_control.json").read_text()
    )
    rows = {}
    for key, arm in arms.items():
        rows[key] = (
            arm["center_coordinates"],
            False,
            R(arm["direction_error_with_serialization_upper"]),
        )
    old = json.loads(
        (measure.HERE.parent / "connection/force_entropy.json").read_text()
    )
    assert parse_interval(old["true_kernel_correction_HS_upper"]).b < c.ri("5.676e-7")
    rows["fixed_ranges_and_forces_correlation"] = (
        old["direction_center_even_then_odd"],
        True,
        R("5.68e-7"),
    )
    out = {}
    for key, (row, restricted, correction) in rows.items():
        frame = [c.prior.sym_matrix(z) for z in c.iso.basis(restricted)]
        ds = [
            sum((R(v) * a for v, a in zip(part, frame)), sp.zeros(4))
            for part in np.array(row).reshape(2, -1)
        ]
        out[key] = ds, correction
    return out


def analyzers():
    out = []
    for i in range(4):
        d = sp.zeros(4)
        d[i, i] = 1
        out.append(d)
    for i, j in combinations(range(4), 2):
        x, y = sp.zeros(4), sp.zeros(4)
        x[i, j] = x[j, i] = 1
        y[i, j], y[j, i] = -sp.I, sp.I
        out.extend([x, y])
    return out


def measurements(board, route, rings, q, axes, arms):
    srcs, radius, parity = measure.source()
    p = R("0.1343")
    assert max(abs(parity.a - c.ri(p)), abs(parity.b - c.ri(p))) <= c.ri(
        ".000300000001"
    )
    sectors = []
    for odd in (False, True):
        h = measure.bridge.native.generator(board, route, q, odd)
        w = measure.four_tick_columns(
            h,
            measure.bridge.beam.indices(route),
            measure.bridge.beam.indices(sum(rings, [])),
            measure.bridge.native.prior.weights()[7],
        )
        entry = []
        for a in axes:
            imaginary = any(sp.im(v) for v in a)
            real = a / (sp.I if imaginary else 1)
            op = sparse.csr_matrix(
                np.kron(np.array(real, dtype=int), np.eye(4, dtype=int))
            )
            entry.append(measure.compress(w, op, -sp.I if imaginary else 1))
        sectors.append(entry)
    # Four port blocks, each with operator error <1e-19, so ||dU||<2.1e-19;
    # every analyzer has norm 1. The averaged compression error is <1e-18.
    prop = R("1e-18")
    values = {"reference": []}
    errors = {"reference": []}
    for key in (*arms, "low", "high"):
        values[key], errors[key] = [], []
    for pair in zip(*sectors):
        d = [v - v.trace() * sp.eye(4) / 4 for v in pair]
        ns = [iv.sqrt(c.ri(c.inner(v, v))) for v in d]
        a, b = [c.inner(v, s) for v, s in zip(pair, srcs)]
        for key, mix in (("low", R(".134")), ("high", R(".1346"))):
            values[key].append((1 - mix) * a + mix * b)
            errors[key].append(c.ri(prop) + c.ri(radius) * max(x.b for x in ns))
        values["reference"].append((1 - p) * a + p * b)
        errors["reference"].append(
            c.ri(prop)
            + c.ri(radius) * max(x.b for x in ns)
            + c.ri(".000300000001") * abs(c.ri(a - b))
        )
        for key, (ds, corr) in arms.items():
            values[key].append(sum(c.inner(v, s) for v, s in zip(d, ds)))
            # Trace norm <= sqrt(8) times the certified tangent norm <1.01.
            errors[key].append(
                c.ri(3 * prop) + c.ri(corr) * iv.sqrt(sum(v**2 for v in ns))
            )
    return {key: (val, errors[key]) for key, val in values.items()}


def matrices(readings, axes, indices):
    parts = []
    for a in axes:
        v = a.extract(indices, indices)
        parts.append(v / c.inner(a, a))
    out = {}
    for key, (val, err) in readings.items():
        m = sum((n * a for n, a in zip(val, parts)), sp.zeros(len(indices)))
        e = iv.sqrt(sum(c.ri(c.inner(v, v)) * z**2 for v, z in zip(parts, err)))
        assert m == m.H
        out[key] = m, e
    return out


def norm(m):
    return iv.sqrt(c.ri(c.inner(m, m)))


def logarithm(m, radius):
    dim = m.rows
    w, v = np.linalg.eigh(np.array(m, complex))
    vec = c.prior.sym_matrix(v)
    diag = sp.diag(*(R(float(y)) for y in w))
    logs = sp.diag(*(R(float(np.log(y))) for y in w))
    for y, log_y in zip(diag.diagonal(), logs.diagonal()):
        assert abs(iv.ln(c.ri(y)) - c.ri(log_y)).b < c.ri("1e-12")
    defect = c.prior.one_norm(vec.H * vec - sp.eye(dim))
    resid = c.prior.one_norm(m - vec * diag * vec.H)
    assert defect < R("1e-13") and resid < R("1e-13")
    assert max(diag.diagonal()) < 1 and min(logs.diagonal()) > -21
    # Replace vec by its unitary polar factor. Changing V diag V* costs
    # <=3*defect, and changing V log(diag) V* costs <=63*defect.
    delta = c.ri(resid + 3 * defect) + radius
    spec = [c.ri(y) + SIGNS * delta for y in diag.diagonal()]
    lo = spec[0].a
    assert lo > 0
    error = delta / lo + c.ri(63 * defect + R("2e-12"))
    k = (-vec * logs * vec.H).applyfunc(sp.expand)
    return k, error, spec


def state_description(m):
    ref, error = m["reference"]
    k, log_error, spec = logarithm(ref, error)
    return ref, error, k, log_error, spec


def study(blocks, active, key):
    weight = R(1, len(active))
    dim = blocks[active[0]][0]["reference"][0].rows
    root_dim = iv.sqrt(c.ri(dim))
    p0 = 1 - weight * sum(blocks[q][1][0].trace() for q in active)
    outside_error = c.ri(weight) * root_dim * sum(blocks[q][1][1] for q in active)
    outside = c.ri(p0) + SIGNS * outside_error
    score, err, curvature = 0, 0, 0
    dout, edout = 0, 0
    spectrum = [outside]
    t = c.ri(".01")
    for q in active:
        table, (_, _, k, log_error, spec) = blocks[q]
        dm, edm = table[key]
        dnorm = norm(dm)
        lerr = log_error + outside_error / outside.a
        s = iv.ln(c.ri(p0 / weight))
        score += c.ri(weight) * (c.ri(c.inner(k, dm)) + s * c.ri(dm.trace()))
        # True score operator is log(p_out/(wB)); its spectral norm is
        # bounded by its values at the two spectral endpoints.
        gain = max(
            abs(iv.ln(outside / (c.ri(weight) * x))).b for x in (spec[0], spec[-1])
        )
        err += c.ri(weight) * root_dim * (dnorm * lerr + edm * gain)
        eig = spec[0].a - t * (dnorm + edm)
        assert eig > 0
        curvature += c.ri(weight) * (dnorm + edm) ** 2 / eig
        spectrum.extend(c.ri(weight) * x for x in spec)
        dout -= weight * dm.trace()
        edout += c.ri(weight) * root_dim * edm
    oslope = c.ri(dout) + SIGNS * edout
    assert outside.a - t * abs(oslope).b > 0
    curvature += abs(oslope).b ** 2 / (outside.a - t * abs(oslope).b)
    slope = score + SIGNS * err
    finite = iv.mpf([(t * slope - t * t * curvature / 2).a, (t * slope).b])
    entropy = -sum(x * iv.ln(x) for x in spectrum)
    return dict(
        entropy=str(entropy),
        spectrum=list(map(str, spectrum)),
        slope=str(slope),
        finite_change_0_to_0p01=str(finite),
        outside_probability_slope=str(oslope),
        second_derivative_magnitude_upper=str(curvature.b),
    )


def density_entropy(spec):
    return -sum(v * iv.ln(v) for v in spec)


def tight_response(tables, cache, active, arm):
    w = R(1, len(active))
    step = R(".01")
    dim = tables[active[0]][arm][0].rows
    rd = iv.sqrt(c.ri(dim))
    p0s, p1s = [], []
    bounds = {}
    for q in active:
        tab = tables[q]
        dm, ed = tab[arm]
        dmnorm = norm(dm)
        ev = []
        for parity in ("low", "high"):
            m, eb = tab[parity]
            for s in (0, step):
                code = (q, parity, arm if s else "reference")
                if code not in cache:
                    cache[code] = logarithm(m + s * dm, c.ri(0))
                ev.append(cache[code][2])
        base_err = max(tab[p][1].b for p in ("low", "high"))
        lower = min(v[0].a for v in ev)
        upper = max(v[-1].b for v in ev)
        anchor = (tab["low"][0] + tab["high"][0] + step * dm) / 2
        inverse = c.prior.sym_matrix(np.linalg.inv(np.array(anchor, complex)))
        inverse = (inverse + inverse.H) / 2
        # Neumann bound on exact inverse from rational matrix residual.
        residual = c.prior.one_norm(sp.eye(dim) - anchor * inverse)
        assert residual < R("1e-10")
        inv_error = c.ri(c.prior.one_norm(inverse) * residual / (1 - residual))
        weighted_sq = c.ri(c.inner(inverse * dm * inverse, dm))
        weighted = iv.sqrt(
            weighted_sq
            + (2 * c.ri(c.prior.one_norm(inverse)) * inv_error + inv_error**2)
            * dmnorm**2
        )
        swing = (
            (norm(tab["high"][0] - tab["low"][0]) + c.ri(step) * dmnorm) / 2
            + base_err
            + c.ri(step) * ed
        )
        fraction = swing / lower
        assert fraction.b < 1
        bounds[q] = dict(
            dlog_dual_norm=weighted / (1 - fraction),
            n=dmnorm,
            d=ed,
            r=base_err,
            lo=lower - base_err - c.ri(step) * ed,
            hi=upper + base_err + c.ri(step) * ed,
            variation=norm(tab["high"][0] - tab["low"][0]),
        )
    for parity in ("low", "high"):
        p0 = 1 - w * sum(tables[q][parity][0].trace() for q in active)
        p1 = p0 - step * w * sum(tables[q][arm][0].trace() for q in active)
        p0s.append(p0)
        p1s.append(p1)
    oe = c.ri(w) * rd * sum(z["r"] + c.ri(step) * z["d"] for z in bounds.values())
    outside_min = min(map(c.ri, p0s + p1s)).a - oe
    outside_max = max(map(c.ri, p0s + p1s)).b + oe
    assert outside_min > 0 and all(z["lo"] > 0 for z in bounds.values())
    # A single scalar p interpolates the conditional parity centers. The
    # interpolation error for each nominal log score on x in [0,1] is at
    # most sup|score''(x)|/8, with V=B(high)-B(low). Twice differentiating
    # the log resolvent gives ||log(B+xV)''|| <= ||V||^2/lambda_min^2.
    # Bound the nominal-to-physical score separately, uniformly in x and s.
    outside_variation = abs(c.ri(p0s[1] - p0s[0]))
    error, interp = 0, 0
    for q, z in bounds.items():
        gain = max(
            abs(iv.ln(o / (c.ri(w) * level))).b
            for o in (outside_min, outside_max)
            for level in (z["lo"], z["hi"])
        )
        # In the eigenbasis of B, every entry of D log_B has multiplier
        # at most 1/sqrt(lambda_i lambda_j). Self-adjointness of this map
        # lets the state error pair with the actual log derivative of dm.
        constant_occupation = arm == "fixed_all_number_weights"
        if constant_occupation:
            gain = iv.ln(z["hi"] / z["lo"]) / 2
        inside_error = (z["r"] + c.ri(step) * z["d"]) * z["dlog_dual_norm"]
        outside_score = 0 if constant_occupation else rd * z["n"] * oe / outside_min
        error += c.ri(w) * (inside_error + outside_score + rd * z["d"] * gain)
        interp += (
            c.ri(w)
            * rd
            * z["n"]
            * ((z["variation"] / z["lo"]) ** 2 + (outside_variation / outside_min) ** 2)
            / 8
        )
    slopes = []
    finite = []
    for p, p0, p1 in zip(("low", "high"), p0s, p1s):
        s0 = [c.ri(p0)]
        s1 = [c.ri(p1)]
        score = 0
        for q in active:
            dm = tables[q][arm][0]
            k, ke, spec = cache[q, p, "reference"]
            score += c.ri(w) * (
                c.ri(c.inner(k, dm))
                + iv.ln(c.ri(p0 / w)) * c.ri(dm.trace())
                + SIGNS * rd * bounds[q]["n"] * ke
            )
            s0.extend(c.ri(w) * x for x in spec)
            s1.extend(c.ri(w) * x for x in cache[q, p, arm][2])
        slopes.append(score)
        finite.append(density_entropy(s1) - density_entropy(s0))
    radius = error + interp
    return dict(
        slope=str(
            iv.mpf([min(v.a for v in slopes), max(v.b for v in slopes)])
            + SIGNS * radius
        ),
        finite_change_0_to_0p01=str(
            iv.mpf([min(v.a for v in finite), max(v.b for v in finite)])
            + SIGNS * c.ri(step) * radius
        ),
        uniform_log_score_perturbation_bound=str(error.b),
        parity_interpolation_bound=str(interp.b),
    )


def run():
    iv.dps = 34
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    assert route == (0, 1, 38, 39)
    rings = measure.bridge.beam.screens(board, route)
    axes, arms = analyzers(), controls()
    raw = {q: measurements(board, route, rings, q, axes, arms) for q in (-1, 1)}
    panels = {}
    centers = {}
    midpoint = {}
    for region, sites in (("stations", (0, 1, 2, 3)), ("bridge", (1, 2))):
        tables = {q: matrices(raw[q], axes, sites) for q in raw}
        for tab in tables.values():
            dm, ed = tab["fixed_all_number_weights"]
            # The true tangent has zero trace separately in each occupied
            # sector. Orthogonally projecting the nominal direction is contractive.
            tab["fixed_all_number_weights"] = (
                dm - sp.eye(dm.rows) * dm.trace() / dm.rows,
                ed,
            )
        blocks = {q: (m, state_description(m)) for q, m in tables.items()}
        if region == "stations":
            # Each mode can be read on either actual local total-load sector,
            # with misses returning 0. The fair weight is already 1/2.
            best = (0, None)
            for q, t in tables.items():
                dm, ed = t["fixed_all_number_weights"]
                for j, op in enumerate(axes[4:]):
                    v = c.ri(c.inner(op, dm)) / 2
                    e = norm(op) * ed / 2
                    lower = (abs(v) - e).a
                    if lower > best[0]:
                        best = (lower, (q, j, str(v + SIGNS * e)))
            lo, (q, j, band) = best
            i, k = list(combinations(route, 2))[j // 2]
            midpoint = dict(
                charge=q,
                ports=[i, k],
                quadrature=("real", "imaginary")[j % 2],
                fair_mode_slope=band,
                minimum_midpoint_entropy_deficit_at_symmetric_0p01=str(
                    (c.ri(".01") ** 2 * lo**2 / 2).a
                ),
            )

        centers[region] = {}
        for q, (m, _) in blocks.items():
            v, err = m["reference"]
            centers[region][str(q)] = dict(
                real=np.array(v, complex).real.tolist(),
                imaginary=np.array(v, complex).imag.tolist(),
                HS_error_upper_with_display_rounding=float(
                    np.nextafter(float(err.b) + 1e-14, np.inf)
                ),
            )
        panels[region] = {
            name: {key: study(blocks, active, key) for key in arms}
            for name, active in (
                ("vacancy", (-1,)),
                ("excess", (1,)),
                ("fair", (-1, 1)),
            )
        }
        cache = {}
        for tag, active in (("vacancy", (-1,)), ("excess", (1,)), ("fair", (-1, 1))):
            for key in arms:
                panels[region][tag][key].update(
                    tight_response(tables, cache, active, key)
                )
    for name, panel in panels["stations"].items():
        for key, data in panel.items():
            bridge = panels["bridge"][name][key]
            data["extra_conditional_entropy"] = str(
                parse_interval(data["entropy"]) - parse_interval(bridge["entropy"])
            )
            for value in ("slope", "finite_change_0_to_0p01"):
                # These two controls hold the full dimer density exactly.
                base = (
                    0
                    if key
                    in (
                        "fixed_work_control",
                        "fixed_work_second_control",
                        "fixed_all_number_weights",
                    )
                    else parse_interval(bridge[value])
                )
                data["extra_conditional_" + value] = str(
                    parse_interval(data[value]) - base
                )
    return dict(
        route=route,
        tick=4,
        launches=rings,
        untagged_charge_weights=["1/2", "1/2"],
        source_arms=list(arms),
        occupied_reference_matrices=centers,
        regions=panels,
        fixed_number_midpoint_witness=midpoint,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (Path(__file__).resolve().parent / "station_region.json").write_text(report)
    else:
        print(report, end="")
