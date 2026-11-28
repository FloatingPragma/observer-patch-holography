"""Find and certify common equal-curvature variations of the regional head records."""

import argparse
import json
from hashlib import sha256
from os.path import relpath
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts import (
    held as material,
)
from matched_einstein_limit.source_implication.full_algebra.enclosure import (
    geometry_jet,
)

from . import finite_dose_curvature as c
from . import finite_dose_material as beam

lab, iv, reg = c.lab, c.iv, c.reg
bounds = lab.six.bounds


def pack(x):
    return np.vectorize(lambda z: iv.mpf(float(z)), otypes=[object])(x)


def norm(x):
    return iv.sqrt(sum(q**2 for q in np.ravel(x)))


def certify(matrix, price):
    mid = np.vectorize(lambda z: float(z.mid))(matrix)
    scale = 2.0 ** (-np.floor(np.log2(np.linalg.norm(mid, axis=1))))
    A = matrix * scale[:, None]
    center = mid * scale[:, None]
    q = np.array([float(z.mid) for z in price])
    U, w, V = np.linalg.svd(center, full_matrices=True)
    right = pack((V[: len(w)].T / w) @ U.T)
    projection = V[len(w) :].T @ (V[len(w) :] @ q)
    d = pack(projection / np.linalg.norm(projection))
    coeff = pack(q @ np.linalg.pinv(center))
    identity = np.eye(len(w), dtype=int)
    error = identity - A @ right
    defect = norm(error)
    out = {
        "constraint_count": len(w),
        "right_inverse_defect": bounds(defect),
        "proposal": [float(v.mid) for v in d],
    }
    # Certify a box for z = A d + (I-A B)z. The numerical solve only
    # proposes radii. Strict interval inclusion proves contraction in the
    # resulting weighted maximum norm and encloses the true solution.
    rhs = A @ d
    absolute = np.vectorize(lambda v: float(abs(v).b))
    try:
        radii = 1.05 * np.linalg.solve(
            np.eye(len(w)) - absolute(error), absolute(rhs) + 1e-15
        )
    except np.linalg.LinAlgError:
        return out
    if np.min(radii) <= 0:
        return out
    z = np.array([iv.mpf([-t, t]) for t in radii])
    inclusion = rhs + error @ z
    if any(abs(v).b >= r for v, r in zip(inclusion, radii, strict=True)):
        return out
    correction = norm(right @ z)
    out["correction_norm_upper"] = float(correction.b)
    # d-Bz is a true common kernel vector. Projection on d bounds its
    # norm from below even when the perpendicular correction is large.
    denominator = ((sum(d * d) - (d @ right) @ z) / norm(d)).a
    if denominator <= 0:
        return out
    residual = price - coeff @ A
    numerator = residual @ d - (residual @ right) @ z
    band = numerator / iv.mpf([denominator, (norm(d) + correction).b])
    out["entropy_per_unit_true_tangent"] = bounds(band)
    return out


def derivatives(bins, x, hs, materials):
    basis = c.dose.heads.measure.joint_basis()
    deficit = np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    calibrated = 3 + c.station.pair(deficit, x)
    # 4x4 effects are converted once from bounded fixed-point bins. Bound
    # coefficient rounding by a conservative multiple of the entrywise amplitude.
    propagation = float(lab.six.survey.history.timing_error()[1])
    mass = np.einsum("rtbij,t->r", abs(bins), abs(reg.SCORE[:-1] - 3))
    arith = 1e-15 + 2e-14 * mass
    # Any two unit-trace preparations annihilate a common scalar identity.
    scalar = np.trace(deficit, axis1=-2, axis2=-1).sum(axis=1) / 8
    centered = deficit - scalar[:, None, None, None] * np.eye(4)
    radius = hs * np.linalg.norm(centered.reshape(30, -1), axis=1) + propagation + arith
    R = np.array(
        [reg.band(v - e, v + e) for v, e in zip(calibrated, radius, strict=True)],
        dtype=object,
    )
    J = np.array([c.station.pair(deficit, b) for b in basis]).T
    tolerances = 2 * propagation + arith
    J = np.array(
        [
            [reg.band(v - e, v + e) for v in row]
            for row, e in zip(J, tolerances, strict=True)
        ]
    )
    ref = c.along(R, [iv.mpf(0)] * len(R))
    grad = np.array([ref[t]["normalized_log_clock_gradient"] for t in c.FILLINGS])
    shape = grad @ (J / R[:, None])
    ids = [i for i, r in enumerate(lab.six.survey.ROUTES) if max(r) < 4]
    area_jet, volume_jet = geometry_jet(R[ids])
    volume = volume_jet.g @ J[ids]
    # Full entropy and connected charge have the same price on the fixed-marginal slice.
    q = -c.station.pair(basis, c.station.log_blocks(x))
    floor = np.linalg.eigvalsh(x).min() - hs - 1e-12
    logerr = 2e-12 + hs / floor
    price = np.array([reg.band(v - logerr, v + logerr) for v in q])
    grids = {
        "six": shape,
        "six_volume": np.vstack([shape, volume]),
        "six_volume_material": np.vstack([shape, volume, materials]),
    }
    result = {k: certify(m, price) for k, m in grids.items()}
    if hs:
        # At the actual equilibrium, dS=-eta dA on this common volume
        # level. This independent price avoids estimating log(X*) at X0.
        coeff = json.loads(
            Path(c.dose.enclosure.__file__).with_name("certificate.json").read_bytes()
        )
        area_price = -(area_jet.g @ J[ids]) * iv.mpf(coeff["gain_interval"])
        for key in ("six_volume", "six_volume_material"):
            result[key + "_stationarity"] = certify(grids[key], area_price)
    return result


def precision_tape(state, rng):
    # Use the same individual-head analyzers as the finite-dose source, with
    # independent batches that keep every subtotal exactly representable.
    shots = 2**68
    tape = []
    for a in c.dose.heads.AXES[1:]:
        row = []
        for b in c.dose.heads.AXES[1:]:
            effects = [
                np.kron((np.eye(2) + i * a) / 2, (np.eye(2) + j * b) / 2)
                for i in c.dose.heads.SIGNS
                for j in c.dose.heads.SIGNS
            ]
            p = np.einsum("kij,bji->bk", effects, state).real.ravel()
            order = np.argsort(p)
            counts = np.zeros(8, dtype=object)
            for _ in range(shots // (2**50 * 1024)):
                draws = rng.multinomial(2**50, p[order], size=1024)
                counts += draws.sum(axis=0).astype(object)
            row.append([int(v) for v in counts[np.argsort(order)]])
        tape.append(row)
    return tape


def run():
    iv.dps = 36
    source, states, *_ = c.trials()
    x, p = states["reference"], states["independent"]
    prior = material.measure.survey.HERE / "translation.json"
    mat = json.loads(prior.read_bytes())
    M = np.asarray(
        [[iv.mpf(v) / mat["scale"] for v in row] for row in mat["rows"]], dtype=object
    )
    cert = Path(c.dose.enclosure.__file__).with_name("certificate.json")
    hs = json.loads(cert.read_bytes())["state_HS_error_upper"]
    ports, bins = lab.six.survey.return_bins()
    refs = {
        key: derivatives(bins, x, e, M)
        for key, e in (("prepared", 0), ("exact_equilibrium_box", hs + 1e-12))
    }
    # Use the strongest supported certificate on the actual prepared reference.
    feasible = [
        (k, r)
        for k, r in refs["prepared"].items()
        if r.get("entropy_per_unit_true_tangent", [0, 0])[0] > 0
    ]
    if not feasible:
        return {"certificates": refs, "observations": None}
    key, selected = max(feasible, key=lambda item: item[1]["constraint_count"])
    d = np.tensordot(selected["proposal"], c.dose.heads.measure.joint_basis(), axes=1)
    for k in range(1, 27):
        step = 2.0**-k
        ends = [x - step * d, x + step * d]
        if min(np.linalg.eigvalsh(y).min() for y in ends) <= 0:
            continue
        recipes = [c.station.preparation(y) for y in ends]
        if (
            min(min(r["min_kept"], r["min_dephased"]) for r in recipes) > 0
            and max(r["max_error"] for r in recipes) < 1e-11
        ):
            break
    else:
        raise ValueError("Both rail lotteries required")
    rng = np.random.default_rng(28112088)
    K = c.station.log_blocks(p) - c.station.log_blocks(x)
    full = c.station.pair(c.dose.heads.measure.joint_basis(), K)
    readings = []
    times = {}
    for name, y in zip(("minus", "plus"), ends, strict=True):
        tape = precision_tape(y, rng)
        # Use the registered signs within each of the nine appointments.
        weights = np.outer(c.dose.heads.SIGNS, c.dose.heads.SIGNS).ravel().astype(int)
        count = np.array(tape, dtype=object).reshape(3, 3, 2, 4)
        signed = np.asarray(count @ weights, dtype=float).transpose(2, 0, 1).ravel() / (
            2 * 2**68
        )
        val = full @ signed
        delta = np.sqrt(2 * np.log(2 * 36 / 0.01) / float(2**68)) / 2
        err = abs(full).sum() * delta + 2e-13
        readings.append((val, err))
        prob = c.station.pair(bins, y)
        times[name] = {
            "time": [lab.counts(np.r_[v, 1 - v.sum()], rng) for v in prob],
            "pauli": tape,
        }
    contrast, err = readings[1][0] - readings[0][0], sum(r[1] for r in readings)
    acted = c.station.pair(K, 2 * step * d)
    assert abs(contrast - acted) <= err
    times_raw = {k: v["time"] for k, v in times.items()}
    means, _ = reg.timing_slopes(times_raw, lab.SHOTS, (("plus", "minus", 1),))
    diff = np.array(means["plus"]) - np.array(means["minus"])
    geom = c.along(means["minus"], diff)
    anchor = [i for i, r in enumerate(lab.six.survey.ROUTES) if max(r) < 4]
    box = np.array(
        [r + reg.band(0, 1) * z for r, z in zip(means["minus"], diff, strict=True)]
    )
    area = geometry_jet(box[anchor])[0].g @ diff[anchor]
    flight = material.material_slopes(ends[1] - ends[0])
    # Count each of the ten calibrated interior transfers with a separate
    # acquisition. Keep even channels for which the sign is unresolved.
    effects, scales, occupancy = beam.flight()
    internal = []
    exposure = 2**70
    seed = 28112089
    detector = np.random.default_rng(seed)
    for edge in material.measure.CHANNELS:
        scale = scales[edge]
        hits, means = [], []
        for y in ends:
            mean = c.station.pair(effects[edge], y) / scale
            active = c.station.pair(occupancy[edge], y) / scale
            bins_read, score = beam.port_read(
                detector,
                exposure,
                active,
                mean,
                scale,
                np.log(4 * 2 * len(material.measure.CHANNELS) / 0.01),
            )
            hits.append(bins_read)
            means.append(score)
        contrast_interval = means[1] - means[0]
        population = flight[8 + len(internal)]
        assert contrast_interval.a <= population <= contrast_interval.b
        internal.append(
            {
                "edge": list(edge),
                "plus_minus_miss_by_arm": hits,
                "contrast": bounds(contrast_interval),
            }
        )
    return {
        "certificates": refs,
        "interior_link_measurements": {
            "exposure": exposure,
            "seed": seed,
            "links": internal,
        },
        "strongest_prepared_test": key,
        "ordinary_energy_change": 0,
        "signed_recipe_step": step,
        "rail_lotteries": recipes,
        "tomography_shots": 2**68,
        "clock_shots": lab.SHOTS,
        "source_contrast_from_counts": bounds(reg.band(contrast - err, contrast + err)),
        "source_contrast_on_nominal_lottery": acted,
        "measured_area_contrast": bounds(area),
        "measured_curvature_contrasts": {
            str(t): bounds(r["d_normalized"]) for t, r in geom.items()
        },
        "independently_propagated_material_contrasts": flight,
        "tapes": times,
        "seed": 28112088,
        "inputs": {
            relpath(f, lab.HERE): sha256(f.read_bytes()).hexdigest()
            for f in (
                prior,
                cert,
                beam.cuts.HERE / "results.json",
                lab.HERE / "CURVATURE_SOURCE_NULL.md",
                lab.HERE / "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dest = lab.HERE / "curvature_source_null.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), result)
    print(
        {
            k: {
                s: r.get("entropy_per_unit_true_tangent", r["right_inverse_defect"])
                for s, r in v.items()
            }
            for k, v in result["certificates"].items()
        }
    )
    print(result.get("source_contrast_from_counts"))
