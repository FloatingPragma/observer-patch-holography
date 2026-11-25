"""Resolve the oriented translation flux at individual cuts of the material chart."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .. import translation as regional

HERE = Path(__file__).parent
native = regional.translation
survey = regional.measure
iv = survey.iv
ri = regional.cert.ri
# Ten coordinates chosen from the material matrices before using an area or
# volume derivative. The choices are frozen for exact and held-source replays.
CHANNELS = (
    (1, 0, 3),
    (1, 2, 5),
    (-1, 2, 5),
    (-1, 2, 3),
    (1, 2, 3),
    (-1, 0, 3),
    (1, 3, 1),
    (1, 0, 16),
    (-1, 0, 4),
    (-1, 0, 16),
)


def cut_experiment():
    atlas = native.Atlas(1)
    patch = native.chosen(atlas)
    cells = set(patch["cells"])
    center = 0  # the source of the registered partial translation
    radius = 1 + max(atlas.dist[center][u] for u in cells)
    weights = {
        u: sp.Rational(radius**2 - atlas.dist[center][u] ** 2, 2 * radius)
        for u in cells
    }
    board = atlas.board
    route = regional.bridge.native.iso.pairing.formation.region(board)
    launch = regional.bridge.beam.indices(
        sum(regional.bridge.beam.screens(board, route), [])
    )
    edges = native.flux_on_edges(board, patch)
    # For internal edges keep the lower-numbered endpoint; the reverse is
    # the exact negative by flux_on_edges. Exterior edges point outward.
    keys = sorted(k for k in edges if k[::-1] not in edges or k[0] < k[1])
    pauli = [regional.cert.prior.sym_matrix(b) for b in regional.cert.iso.basis(True)]
    coeff = regional.bridge.native.prior.weights()[7]
    data = []
    for charge in (-1, 1):
        blocks = {key: [] for key in keys}
        for odd in (False, True):
            h, operators = native.op_table(board, route, patch, charge, odd)
            assert (sum(edges.values()) - operators["boundary"][0]).count_nonzero() == 0
            indices = np.unique(np.concatenate([edges[k].nonzero()[0] for k in keys]))
            w = regional.four_tick_columns(h, indices, launch, coeff)
            for key in keys:
                mat = edges[key]
                ids = np.unique(mat.nonzero()[0])
                sites = np.unique(ids // 4)
                assert len(sites) <= 4 and abs(mat).sum(axis=1).max() <= 4
                where = np.searchsorted(indices, ids)
                block = regional.region.compress(w[:, where, :], mat[ids, :][:, ids], 4)
                # Each of at most four 4x4 port blocks has propagation error
                # <1e-19. The selected propagator is a contraction, as is M/4.
                # 2*(2*(2e-19)+(2e-19)^2) < 1e-18 per HS Pauli coordinate.
                blocks[key].extend(
                    survey.pack(
                        regional.cert.ri_band(
                            regional.cert.inner(block, b), sp.Rational("1e-18")
                        )
                    )
                    for b in pauli
                )
        for u, v in keys:
            data.append(
                {
                    "channel": [charge, u, v],
                    "aperture_weight": str(weights[u] - weights.get(v, 0)),
                    "row": blocks[u, v],
                }
            )
    return {"radius_in_graph_steps": radius, "cuts": data}


def interval_rows(packed):
    return np.array(
        [[iv.mpf([lo, hi]) / survey.SCALE for lo, hi in row] for row in packed],
        dtype=object,
    )


def mid(array):
    return np.vectorize(lambda x: float(x.mid), otypes=[float])(array)


def fixed(array):
    return np.vectorize(lambda x: ri(f"{x:.16g}"), otypes=[object])(array)


def norm(array):
    return iv.sqrt(sum(v**2 for v in array.flat))


def exact_entropy_center():
    source = json.loads((survey.HERE.parent / "results.json").read_text())
    state = np.array(source["x_real"]) + 1j * np.array(source["x_imag"])
    blocks = [regional.cert.prior.sym_matrix(b) for b in state]
    blocks = [(b + b.H) / 2 for b in blocks]
    mass = sum(b.trace() for b in blocks)
    pauli = [regional.cert.prior.sym_matrix(b) for b in regional.cert.iso.basis(True)]
    result = []
    for block in blocks:
        coords, error = regional.cert.prior.entropy_coordinates(
            block / mass, 0, sp.Rational(".02"), pauli
        )
        result.extend(regional.cert.ri_band(v, 2 * error) for v in coords)
    return np.array(result, dtype=object)


def aperture_test(matrix, area):
    magnitudes = np.linalg.norm(mid(matrix), axis=1)
    scales = np.array(
        [iv.mpf(2) ** int(-np.floor(np.log2(w))) for w in magnitudes], dtype=object
    )
    b = matrix * scales[:, None]
    right = np.linalg.pinv(mid(b))
    c = mid(area) @ right
    residual = mid(area) - c @ mid(b)
    d = fixed(-residual / np.linalg.norm(residual))
    right = fixed(right)
    defect = norm(np.eye(len(b), dtype=int) - b @ right)
    assert defect.b < 1
    correction = norm(right) * norm(b @ d) / (1 - defect)
    score = area - fixed(c) @ b
    error = iv.mpf([-1, 1]) * correction
    slope = (-sum(score * d) + norm(score) * error) / (norm(d) + error)
    ent = iv.mpf(["23.1042361", "23.1042363"]) * slope
    assert ent.a > ri(".00070")
    return {
        "independent_constraints": len(b),
        "kernel_dimension": 18 - len(b),
        "slope_nats": survey.pack(ent),
        "right_inverse_defect": survey.pack(defect),
    }


def certify(data):
    full = json.loads((survey.HERE / "translation.json").read_text())
    summaries = interval_rows(full["rows"])
    cuts = interval_rows([p["row"] for p in data["cuts"]])
    labels = [tuple(p["channel"]) for p in data["cuts"]]
    material = np.vstack([summaries, cuts[[labels.index(i) for i in CHANNELS]]])
    b = mid(material)
    # Row rescaling improves the floating-point proposal only. All checks,
    # inverse norms and perturbation statements below use the raw native units.
    scale = 1 / np.linalg.norm(b, axis=1)
    right = fixed(np.linalg.inv(b * scale[:, None]) * scale[None, :])
    defect = norm(np.eye(18, dtype=int) - material @ right)
    if defect.b >= 1:
        raise ValueError("material channel rank is uncertified")
    lower = (1 - defect) / norm(right)
    weakest = fixed(np.linalg.svd(b)[2][-1])
    upper = norm(material @ weakest) / norm(weakest)
    assert lower.a > ri("6e-13") and upper.b < ri("9e-13")
    largest = max(norm(r).a for r in material)
    condition = largest / upper.b
    # No area, volume, gain or ordinary transfer enters these read weights.
    k = exact_entropy_center()
    weights = mid(k) @ mid(right)
    decimal = [f"{v:.16g}" for v in weights]
    frozen = np.array([ri(v) for v in decimal], dtype=object)
    residual = norm(k - frozen @ material)
    # Resolvent identity for log on the entire state ball in Hilbert-Schmidt norm.
    bounds = json.loads((survey.HERE.parent / "certificate.json").read_text())
    radius = ri(str(bounds["state_HS_error_upper"]))
    modular_error = residual + radius / (ri(".02") - radius)
    assert modular_error.b < ri("2e-5")
    # Only now compare to geometry and test the prespecified graph aperture.
    geo = json.loads((survey.HERE / "results.json").read_text())
    area, volume = [survey.unpack(geo[key]).ravel() for key in ("area", "volume")]
    aperture = [
        sum(
            (
                row * ri(sp.Rational(p["aperture_weight"]))
                for p, row in zip(data["cuts"], cuts)
                if p["channel"][0] == charge
            ),
            np.full(18, iv.mpf(0), dtype=object),
        )
        for charge in (-1, 1)
    ]
    return dict(
        native_cuts=len(cuts),
        resolved_channels=CHANNELS,
        rank=18,
        inverse_defect=survey.pack(defect),
        min_singular_lower=survey.pack(lower)[0],
        min_singular_upper=survey.pack(upper)[1],
        condition_lower_integer=int(float(condition.a)) - 1,
        spectrum_weights=decimal,
        l1_noise_gain_upper=int(float(sum(abs(v) for v in frozen).b)) + 1,
        modular_error_upper=survey.pack(modular_error)[1],
        aperture=aperture_test(np.vstack([volume, summaries, aperture]), area),
    )


def write(data, result):
    # Keep one oriented channel per line so that the small integer bounds
    # can be reviewed together in their even/odd Pauli ordering.
    content = ",\n".join(
        "    " + json.dumps(row, separators=(",", ":")) for row in data["cuts"]
    )
    (HERE / "effects.json").write_text(
        '{\n  "scale": '
        + str(survey.SCALE)
        + ',\n  "radius_in_graph_steps": '
        + str(data["radius_in_graph_steps"])
        + ',\n  "cuts": [\n'
        + content
        + "\n  ]\n}\n"
    )
    (HERE / "results.json").write_text(
        json.dumps({"scale": survey.SCALE, **result}, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    from . import held

    iv.dps = 36
    data = cut_experiment()
    result = certify(data)
    result["held"] = held.run(result["spectrum_weights"])
    if args.write:
        write(data, result)
    print(json.dumps(result))
