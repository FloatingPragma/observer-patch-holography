"""Transfer the frozen finite correlation source to six complete spatial reconstructions."""

import argparse
import json
from hashlib import sha256
from itertools import combinations

import numpy as np
from scipy.optimize import brentq

from . import finite_pairing_dose as dose

lab = dose.lab
station = dose.station
reg = dose.regional
iv = dose.iv
EDGES = tuple(combinations(range(6), 2))
TETRA = tuple(combinations(range(6), 4))
FILLINGS = tuple(
    t for t in combinations(range(6), 3) if 0 in t and ((4 in t) != (5 in t))
)


def unpack(recipe):
    heads = [(np.eye(2) + s * a) / 2 for a in station.phase.AXES[1:] for s in (-1, 1)]
    atoms = np.array([np.kron(a, b) for a in heads for b in heads])
    kept, mixed = (np.asarray(recipe[k]) for k in ("kept", "dephased"))
    assert min(min(kept), min(mixed)) > 0
    return np.tensordot(np.stack((kept + mixed / 2, mixed / 2)), atoms, axes=1)


def trials():
    source = json.loads((lab.HERE / "finite_pairing_dose.json").read_bytes())
    ordinary = json.loads((lab.HERE / "ordinary_stationarity.json").read_bytes())
    states = {k: unpack(v) for k, v in source["source_recipes"].items()}
    x, p, w = (states[k] for k in ("reference", "independent", "compensator"))
    t, u = source["doses"][-1]["dose"], source["selected_fraction"]
    states["held"] = (1 - t - u) * x + t * p + u * w
    states["ordinary"] = np.stack([x.sum(0) / 2] * 2)
    change = states["ordinary"] - x
    exact = station.entropy(states["held"]) - station.entropy(x)
    max_prob = 1e-4
    h = brentq(
        lambda z: station.entropy(x + z * change) - station.entropy(x) - exact,
        0,
        max_prob,
        xtol=1e-15,
    )
    states["matched_ordinary"] = x + h * change
    slope = reg.band(*ordinary["entropy_slope_from_counts"])
    # S''(x)[D,D] lies between -||D||_HS^2/lambda_min and zero.
    floor = (
        min(
            np.linalg.eigvalsh(x).min(), np.linalg.eigvalsh(x + max_prob * change).min()
        )
        - 1e-12
    )
    lip = np.linalg.norm(change) ** 2 / floor + 1e-12
    rate = slope + reg.band(-max_prob * lip, 0)
    interval = reg.band(*source["held_entropy_gain"]) / rate
    assert 0 < interval.a <= h <= interval.b < max_prob
    return source, states, h, interval


def clock_map():
    """Coefficients of log edge length in the thirty directed log ranges."""
    routes = lab.six.survey.ROUTES
    eye = np.eye(len(routes))
    direct = dict(zip(routes, eye, strict=True))

    def contrast(i, j):
        return direct[i, j] - direct[j, i] if i != j else 0

    gains = [sum(contrast(i, j) for j in range(4)) / 4 for i in range(4)]
    gains.extend(sum(contrast(i, j) + gains[j] for j in range(4)) / 4 for i in (4, 5))
    return np.array(
        [(direct[i, j] + direct[j, i] - gains[i] - gains[j]) / 2 for i, j in EDGES]
    )


def hinge_geometry(square, clock_rates, factor):
    """Schläfli derivative of sum(length * deficit), keeping the full path in the interval."""
    length = {e: square[e].sqrt() for e in EDGES}
    volumes, angles = {}, {}
    for tet in TETRA:
        i, j, k, fourth = tet

        def dot(a, b, c):
            return (square[a, b] + square[a, c] - square[b, c]) * 0.5

        a, b, c = (square[i, j], square[i, k], square[i, fourth])
        g, h, f = dot(i, j, k), dot(i, j, fourth), dot(i, k, fourth)
        six_v = (a * b * c + 2 * g * h * f - a * f * f - b * h * h - c * g * g).sqrt()
        volumes[tet] = six_v * (iv.mpf(1) / 6)
        for edge in combinations(tet, 2):
            i, j = edge
            k, fourth = (v for v in tet if v not in edge)
            normal_dot = square[i, j] * dot(i, k, fourth) - dot(i, j, k) * dot(
                i, j, fourth
            )
            angles[tet, edge] = iv.atan2(length[edge].v * six_v.v, normal_dot.v)
    for tri in FILLINGS:
        cells = [q for q in TETRA if len(set(q) & set(tri)) == 2]
        vol = sum(volumes[q] for q in cells)
        deficit = {
            e: 2 * iv.pi - sum(angles[q, e] for q in cells if set(e) <= set(q))
            for e in EDGES
        }
        price = sum(v * length[e].v for e, v in deficit.items())
        weights = []
        for edge in EDGES:
            dw = iv.mpf(0)
            for q in cells:
                if set(edge) <= set(q):
                    other = tuple(v for v in q if v not in edge)
                    angle = angles[q, other]
                    dw += square[edge].v * length[other].v * iv.cot(angle) / 6
            weights.append(dw)
        # Combine coefficients of the SAME physical clocks before propagating
        # their counting error. Each tetra has dV/dlog L_e = L_e² L_op cot(theta_op)/6.
        norm = vol.v ** (-iv.mpf(1) / 3)
        terms = [length[e].v * deficit[e] for e in EDGES]

        def differential(gradient):
            coeff = clock_map().T @ np.asarray(gradient, dtype=object)
            return factor * sum(a * b for a, b in zip(coeff, clock_rates, strict=True))

        change = differential(terms)
        shape_grad = [
            norm * (c - price * v / (3 * vol.v))
            for c, v in zip(terms, weights, strict=True)
        ]
        shape = differential(shape_grad)
        yield (
            tri,
            {
                "volume": vol.v,
                "integrated": price,
                "normalized": price * norm,
                "d_volume": differential(weights),
                "d_integrated": change,
                "d_normalized": shape,
                "normalized_log_clock_gradient": clock_map().T
                @ np.asarray(shape_grad, dtype=object),
            },
        )


def along(reference, delta, factor=1, at_reference=False):
    # All directed means are affine under a source mixture. Enclose dF/dz
    # throughout z in [0,1], so its interval also encloses the finite difference.
    path = [
        r + reg.band(0, 1) * factor * dr for r, dr in zip(reference, delta, strict=True)
    ]
    if at_reference:
        path = reference
    square = lab.six.squared_range_jets(path, np.array(delta, dtype=object)[:, None])
    return dict(
        hinge_geometry(
            square, [dr / r for dr, r in zip(delta, path, strict=True)], factor
        )
    )


def run():
    iv.dps = 35
    source, states, h, hband = trials()
    ports, bins = lab.six.survey.return_bins()
    seed = 28112082
    rng = np.random.default_rng(seed)
    times = {}
    for key, x in states.items():
        prob = station.pair(bins, x)
        times[key] = [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in prob]
    means, _ = reg.timing_slopes(times, lab.SHOTS, (("held", "reference", 1),))
    x = np.array(means["reference"], dtype=object)
    ds = {k: np.array(a, dtype=object) - x for k, a in means.items()}
    t, u = source["doses"][-1]["dose"], source["selected_fraction"]

    stock = {
        "held_from_counts": ds["held"],
        "held_from_components": t * ds["independent"] + u * ds["compensator"],
        "ordinary_from_counts": ds["matched_ordinary"],
        "ordinary_from_components": ds["ordinary"],
    }
    shapes = {
        k: along(x, v, hband if k == "ordinary_from_components" else 1)
        for k, v in stock.items()
    }
    ordinary_receipt = json.loads(
        (lab.HERE / "ordinary_stationarity.json").read_bytes()
    )
    h_tangent = reg.band(*source["held_connected_modular_gain"]) / reg.band(
        *ordinary_receipt["entropy_slope_from_counts"]
    )
    joint_tangent = along(x, stock["held_from_components"], at_reference=True)
    ordinary_tangent = along(x, ds["ordinary"], h_tangent, at_reference=True)
    # Cross-check the integral enclosure against direct endpoint geometry from
    # population probabilities, independently of the simulated radar counts.
    values = {
        k: 3 + station.pair(bins, states[k]) @ (reg.SCORE[:-1] - 3)
        for k in ("reference", "held", "matched_ordinary")
    }
    direct = {
        k: along([iv.mpf(float(z)) for z in r], [iv.mpf(0)] * len(r))
        for k, r in values.items()
    }
    result = []
    for tri in FILLINGS:
        row = {"missing_triangle": list(tri)}
        for obs in ("volume", "integrated", "normalized"):
            ref = direct["reference"][tri][obs]
            row[obs] = {"reference": lab.six.bounds(ref)}
            for name, v in shapes.items():
                est = v[tri]["d_" + obs]
                arm = "held" if name.startswith("held") else "matched_ordinary"
                diff = direct[arm][tri][obs] - ref
                assert est.a <= diff.b and diff.a <= est.b, (tri, obs, name, est, diff)
                row[obs][name] = lab.six.bounds(est)
                row[obs][arm + "_population_change"] = lab.six.bounds(diff)
            contrast = (
                shapes["held_from_components"][tri]["d_" + obs]
                - shapes["ordinary_from_components"][tri]["d_" + obs]
            )
            row[obs]["held_minus_equal_entropy_ordinary"] = lab.six.bounds(contrast)
        a = joint_tangent[tri]["d_normalized"]
        b = ordinary_tangent[tri]["d_normalized"]
        row["equal_entropy_first_variation"] = {
            k: lab.six.bounds(v)
            for k, v in (("joint", a), ("ordinary", b), ("residual", a - b))
        }
        result.append(row)
    return {
        "seed": seed,
        "shots_per_arm_route": lab.SHOTS,
        "ports": list(ports),
        "ordinary_fraction": h,
        "ordinary_fraction_for_equal_first_variation": lab.six.bounds(h_tangent),
        "ordinary_fraction_from_head_counts": lab.six.bounds(hband),
        "ordinary_parked_energy_increment": h
        * float(np.trace((states["ordinary"] - states["reference"])[1]).real),
        "fillings": result,
        "times": times,
        "inputs": {
            f: sha256((lab.HERE / f).read_bytes()).hexdigest()
            for f in (
                "FINITE_DOSE_CURVATURE.md",
                "EQUAL_SOURCE_TANGENT.md",
                "finite_pairing_dose.json",
                "ordinary_stationarity.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dest = lab.HERE / "finite_dose_curvature.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), result)
    for row in result["fillings"]:
        v = row["normalized"]
        print(
            row["missing_triangle"],
            v["held_from_components"],
            v["ordinary_from_components"],
        )
