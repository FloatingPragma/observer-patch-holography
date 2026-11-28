"""Infer a source-strength floor over arbitrary radial transport distributions."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.constants import c, hbar, proton_mass, parsec
from scipy.optimize import linprog
from scipy import sparse

from . import local_well as local
from . import well_certificate as certificate
from . import kernel_certificate as exact

HERE = local.HERE
WIDTH = local.WINDOW
LOG_WIDTH = np.log((1 + WIDTH) / (1 - WIDTH))
# IAU nominal solar mass parameter divided by CODATA 2022 G, in kg.
SOLAR_MASS = 1.3271244e20 / 6.67430e-11


def bounds(rows):
    v, e, b = [
        np.asarray(rows[k]) for k in ("observed_kms", "error_kms", "baryons_kms")
    ]
    return np.maximum(0, v - 2 * e) ** 2 - b * b, (v + 2 * e) ** 2 - b * b


def corners(radius, edges):
    ends = np.r_[radius * (1 - WIDTH), radius * (1 + WIDTH)]
    return np.unique(
        np.r_[
            0, ends, (ends[:, None] + edges).ravel(), abs(ends[:, None] - edges).ravel()
        ]
    )


def density(r, offsets, annuli):
    """Spherically averaged density for each offset, including the zero limit."""
    edges, weights = annuli
    a, b = edges[:-1], edges[1:]
    slope = weights / (b * b - a * a)
    primitive = np.r_[0, np.cumsum(slope * (b - a))]
    h = np.asarray(offsets)
    result = np.zeros_like(h, dtype=float)
    positive = h > 0
    hp = h[positive]
    overlap = np.interp(r + hp, edges, primitive) - np.interp(
        abs(r - hp), edges, primitive
    )
    result[positive] = overlap / (4 * np.pi * r * hp)
    # Direct edge-relative coordinates avoid cancellation for arbitrarily
    # small positive offsets around a catalogue corner.
    small = (h > 0) & (h < 1e-5 * r)
    if np.any(small):
        ratio = (edges - r) / h[small, None]
        result[small] = (np.diff(np.clip(ratio, -1, 1), axis=1) @ slope) / (
            4 * np.pi * r
        )
    side = ((a < r) & (r < b)).astype(float) + 0.5 * ((a == r) | (b == r))
    result[~positive] = side @ slope / (2 * np.pi * r)
    result[h >= r + edges[-1]] = 0
    return result


def bank(rows, annuli, offsets, mass):
    return mass * np.array(
        [
            (
                density((1 - WIDTH) * r, offsets, annuli)
                - density((1 + WIDTH) * r, offsets, annuli)
            )
            / LOG_WIDTH
            for r in rows["radius_kpc"]
        ]
    )


def one_ring(matrix, lo, hi, signed):
    top = matrix.max(axis=1).clip(min=0)
    bottom = matrix.min(axis=1).clip(max=0)
    upper = top - bottom if signed else top
    lower = top - bottom if signed else -bottom
    with np.errstate(divide="ignore", invalid="ignore"):
        values = np.where(lo > 0, lo / upper, np.where(hi < 0, -hi / lower, 0))
    index = int(np.argmax(values))
    return index, float(values[index])


def dual_ray(direction, row_scale):
    ray = direction / row_scale
    peak = np.max(abs(ray))
    if not peak:
        return ["0"] * len(ray)
    return [str(Fraction(round(v * 10**10), 10**10)) for v in ray / peak]


def local_inverse(matrix, lo, hi, signed):
    """Column generation prices every corner after each restricted LP."""
    ring, floor = one_ring(matrix, lo, hi, signed)
    if not np.isfinite(floor):
        direction = np.zeros(len(lo))
        direction[ring] = 1 if lo[ring] > 0 else -1
        return dict(
            feasible=False,
            single_ring=ring,
            direct_sign_conflict=True,
            dual_ray=dual_ray(direction, np.ones(len(lo))),
        )
    unit = max(1, floor)
    row_scale = np.maximum(1, np.maximum(abs(lo), abs(hi)))
    k = matrix * (unit / row_scale[:, None])
    lower, u = lo / row_scale, hi / row_scale
    size = k.shape[1]
    selected = set(
        np.r_[
            np.linspace(0, size - 1, min(size, 100), dtype=int),
            np.argmax(k, axis=1),
            np.argmin(k, axis=1),
        ].tolist()
    )
    full = np.column_stack((k, -k)) if signed else k
    if signed:
        selected.update(i + size for i in list(selected))
    for phase in ("feasibility", "budget"):
        for attempt in range(300):
            ix = np.array(sorted(selected))
            draw = full[:, ix]
            a = np.vstack((-draw, draw))
            rhs = np.r_[-lower, u]
            if phase == "feasibility":
                a = np.column_stack((a, -np.ones(len(rhs))))
            else:
                active = np.array([ix < size] + ([ix >= size] if signed else []))
                a = np.vstack(
                    (
                        np.column_stack((a, np.zeros(len(rhs)))),
                        np.column_stack((active, -np.ones(len(active)))),
                    )
                )
                rhs = np.r_[rhs, np.zeros(len(active))]
            fit = linprog(
                np.r_[np.zeros(len(ix)), 1],
                A_ub=a,
                b_ub=rhs,
                bounds=(0, None),
                method="highs",
                options={
                    "primal_feasibility_tolerance": 1e-9,
                    "dual_feasibility_tolerance": 1e-9,
                },
            )
            if not fit.success:
                raise RuntimeError(f"{phase}: {fit.message}")
            if phase == "feasibility" and fit.fun < 1e-9:
                break
            dual = fit.ineqlin.marginals
            direction = dual[len(lower) : 2 * len(lower)] - dual[: len(lower)]
            profit = full.T @ direction
            if phase == "budget":
                profit[:size] += dual[2 * len(lower)]
                if signed:
                    profit[size:] += dual[-1]
            profit[ix] = 0
            active = np.flatnonzero(profit > 2e-8)
            if not len(active):
                if phase == "feasibility":
                    return dict(
                        feasible=False,
                        single_ring=ring,
                        normalized_violation=fit.fun,
                        dual_ray=dual_ray(direction, row_scale),
                    )
                break
            selected.update(active[np.argsort(profit[active])[-24:]].tolist())
        else:
            raise RuntimeError("Continuous offset pricing did not converge")
    coeff = np.zeros(full.shape[1])
    coeff[ix] = fit.x[:-1] * unit
    coeff = coeff[:size] - coeff[size:] if signed else coeff
    nonzero = np.flatnonzero(abs(coeff) > 0)
    v = matrix[:, nonzero] @ coeff[nonzero]
    tol = 2e-6 * np.maximum(1, np.maximum(abs(lo), abs(hi)))
    assert np.all(v >= lo - tol) and np.all(v <= hi + tol)
    total = max(sum(coeff.clip(min=0)), sum((-coeff).clip(min=0)))
    return dict(
        feasible=True,
        single_ring=ring,
        witness_amplitude=float(total),
        dual_ray=dual_ray(direction, row_scale),
        knot_indices=nonzero.tolist(),
        weights=coeff[nonzero].tolist(),
        worst_squared_speed_violation=float(max(0, max(lo - v), max(v - hi))),
        predicted_extra_speed2=v.tolist(),
    )


def exact_ring(rows, annuli, mass, index):
    """Rational endpoint extrema of 4π log(51/49) times the unit response."""
    Q = Fraction
    edges, weights = [[certificate.decimal(x) for x in group] for group in annuli]
    r = certificate.decimal(rows["radius_kpc"][index])
    low = certificate.slope_events(Q(49, 50) * r, edges, weights, [(1, Q(1))])
    high = certificate.slope_events(Q(51, 50) * r, edges, weights, [(1, Q(1))])
    knots = sorted(set(low) | set(high))
    value = Q(0)
    slope = Q(0)
    responses = [Q(0)]
    mass = certificate.decimal(mass)
    for a, b in zip(knots, knots[1:]):
        slope += low.get(a, 0) - high.get(a, 0)
        if a == 0:
            responses.append(mass * slope)
        value += slope * (b - a)
        responses.append(mass * value / b)
    assert value == 0
    v, e, bar = [
        certificate.decimal(rows[k][index])
        for k in ("observed_kms", "error_kms", "baryons_kms")
    ]
    lo = max(Q(0), v - 2 * e) ** 2 - bar**2
    hi = (v + 2 * e) ** 2 - bar**2
    floor = max(Q(0), lo, -hi) / (max(responses) - min(responses))
    return dict(
        index=index,
        radius_kpc=float(r),
        scaled_extrema=[str(min(responses)), str(max(responses))],
        required_speed2=[str(lo), str(hi)],
        rational_factor=str(floor),
        signed_floor=float(floor) * 4 * np.pi * LOG_WIDTH,
    )


def train(banks, inputs, held, signed):
    rows = []
    targets = []
    n = len(inputs) - 1
    for name in sorted(set(inputs) - {held}):
        lower, u = bounds(inputs[name])
        scale = np.maximum(1, (u - lower) / 2) * len(lower) * n
        obs = np.asarray(inputs[name]["observed_kms"])
        bar = np.asarray(inputs[name]["baryons_kms"])
        rows.append(banks[name] / scale[:, None])
        targets.extend((obs**2 - bar**2) / scale)
    matrix = np.concatenate(rows)
    target = np.array(targets)
    # Condition columns without changing the empirical mass norm.
    scale = np.maximum(1e-14, np.max(abs(matrix), axis=0))
    k = sparse.csr_matrix(matrix / scale)
    if signed:
        k = sparse.hstack((k, -k), format="csr")
    eye = sparse.eye(len(target))
    a = sparse.vstack(
        (sparse.hstack((k, -eye)), sparse.hstack((-k, -eye))), format="csr"
    )
    b = np.r_[target, -target]
    objective = np.r_[np.zeros(k.shape[1]), np.ones(len(target))]
    first = linprog(objective, A_ub=a, b_ub=b, method="highs")
    if not first.success:
        raise RuntimeError(first.message)
    a = sparse.vstack((a, sparse.csr_matrix(objective)))
    norm = 1 / scale
    if signed:
        norm = np.r_[norm, norm]
    # Duplicate signed columns have identical scales and separate nonnegative masses.
    second = linprog(
        np.r_[norm / np.median(norm), np.zeros(len(target))],
        A_ub=a,
        b_ub=np.r_[b, first.fun + 1e-6],
        method="highs",
    )
    if not second.success:
        raise RuntimeError(second.message)
    x = second.x[: k.shape[1]]
    coeff = (x[: len(scale)] - x[len(scale) :]) / scale if signed else x / scale
    speed2 = np.asarray(inputs[held]["baryons_kms"]) ** 2 + banks[held] @ coeff
    return dict(weights=coeff.tolist(), squared_speed=speed2.tolist(), loss=first.fun)


def load_inputs():
    prior = HERE / "source_profile.json"
    inputs = json.loads(prior.read_bytes())["resolved_maps"]["inputs"]
    charts = local.profiles.mapped_annuli(local.profiles.maps.enroll()[0])
    meta, _ = local.profiles.galaxies.catalogue()
    masses = {
        k: 0.5 * meta[k]["L36_gsun"]
        + local.profiles.galaxies.sparc.GAS_HELIUM * meta[k]["MHI_gsun"]
        for k in inputs
    }
    return inputs, charts, meta, masses


def units(amplitude):
    # amplitude in (km/s)^2 kpc^3 per 10^9 solar masses. kappa=8.
    alpha_si = amplitude * 1e6 * (1000 * parsec) ** 3 / (1e9 * SOLAR_MASS)
    volume_per_kg = alpha_si / ((7 / 6) * c**2)
    site = (volume_per_kg * proton_mass) ** (1 / 3)
    proton_length = np.sqrt(192) * hbar / (proton_mass * c)
    return dict(
        required_volume_per_feed_mass_m3_per_kg=volume_per_kg,
        min_site_length_one_unit_per_proton_m=site,
        proton_probe_kinetic_rest_site_m=proton_length,
        required_to_proton_calibrated_strength=(site / proton_length) ** 3,
        max_probe_mass_for_proton_feed_eV_c2=(np.sqrt(192) * hbar / (site * c))
        * c**2
        / 1.602176634e-19,
    )


def run():
    inputs, charts, meta, masses = load_inputs()
    result = {}
    fixed = {}
    grid = np.r_[0, local.profiles.RATIOS]
    for name, rows in inputs.items():
        offsets = corners(np.asarray(rows["radius_kpc"]), charts[name][0])
        matrix = bank(rows, charts[name], offsets, masses[name])
        lo, hi = bounds(rows)
        arms = {}
        for label, signed in (("positive", False), ("signed", True)):
            contraction = np.minimum(0.001, (hi - lo) / 100)
            fit = local_inverse(matrix, lo + contraction, hi - contraction, signed)
            if fit["feasible"]:
                fit["offset_kpc"] = offsets[fit.pop("knot_indices")].tolist()
                fit["witness"] = exact.witness(rows, charts[name], masses[name], fit)
                fit["certificate"] = exact.dual(
                    rows, charts[name], masses[name], fit["dual_ray"], signed
                )
                assert fit["certificate"]["lower"] <= fit["witness"]["upper"] * (
                    1 + 1e-12
                )
            arms[label] = fit
        ex = exact_ring(rows, charts[name], masses[name], arms["signed"]["single_ring"])
        np.testing.assert_allclose(
            ex["signed_floor"], one_ring(matrix, lo, hi, True)[1], rtol=1e-9
        )
        result[name] = dict(
            corner_count=len(offsets),
            fits=arms,
            exact=ex,
            ordinary=units(arms["signed"]["certificate"]["lower"]),
        )
        fixed[name] = bank(
            rows, charts[name], grid * meta[name]["disk_scale_kpc"], masses[name]
        )
    held = {}
    for label, signed in (("positive", False), ("signed", True)):
        trial = {k: train(fixed, inputs, k, signed) for k in inputs}
        deviations = []
        hits = invalid = 0
        for k, row in trial.items():
            v = np.array(row["squared_speed"])
            err = abs(np.sqrt(np.maximum(v, 0)) - inputs[k]["observed_kms"])
            hits += int(
                np.sum((v >= 0) & (err <= 2 * np.array(inputs[k]["error_kms"])))
            )
            invalid += int(np.sum(v < 0))
            deviations.extend(err)
        held[label] = dict(
            galaxies=trial,
            median_abs_kms=float(np.median(deviations)),
            within_two=hits,
            no_orbits=invalid,
        )
    files = [
        HERE / "KERNEL_BUDGET.md",
        HERE / "source_profile.json",
        local.profiles.maps.HERE / "leroy_table7.dat",
        local.profiles.maps.HERE / "leroy_sample.dat",
        local.profiles.galaxies.sparc.DATA / "table1.dat",
        local.lab.HERE / "undisturbed_medium.json",
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        strength_units="(km/s)^2 kpc^3 per 10^9 solar masses",
        grid_disk_scales=grid.tolist(),
        individual=result,
        held=held,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "kernel_budget.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    for k, v in result["individual"].items():
        print(
            k,
            v["fits"]["positive"]["feasible"],
            v["fits"]["signed"].get("witness_amplitude"),
            v["exact"]["signed_floor"],
        )
    print(
        {k: (v["within_two"], v["median_abs_kms"]) for k, v in result["held"].items()}
    )
