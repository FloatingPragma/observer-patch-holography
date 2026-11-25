"""Complete the reconstructing material reader with the 13 measured local coordinates."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import measure as cuts

HERE = cuts.HERE
regional = cuts.regional
cert = regional.cert
iv = cuts.iv


def local_basis():
    pauli = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    return [sp.kronecker_product(p, sp.eye(2)) / 2 for p in pauli] + [
        sp.kronecker_product(sp.eye(2), p) / 2 for p in pauli
    ]


def native_local_effects():
    atlas = cuts.native.Atlas(1)
    board = atlas.board
    patch = cuts.native.chosen(atlas)
    route = regional.bridge.native.iso.pairing.formation.region(board)
    launch = regional.bridge.beam.indices(
        sum(regional.bridge.beam.screens(board, route), [])
    )
    edges = cuts.native.flux_on_edges(board, patch)
    local = local_basis()
    channels = {}
    for q in (-1, 1):
        for odd in (False, True):
            h, ops = cuts.native.op_table(board, route, patch, q, odd)
            chosen = {(q, k): ops[k] for k in regional.NAMES}
            chosen.update(
                {key: (edges[key[1:]], 4) for key in cuts.CHANNELS if key[0] == q}
            )
            idx = np.unique(
                np.concatenate([op.nonzero()[0] for op, _ in chosen.values()])
            )
            w = regional.four_tick_columns(
                h, idx, launch, regional.bridge.native.prior.weights()[7]
            )
            for key, (op, unit) in chosen.items():
                keep = np.unique(op.nonzero()[0])
                indices = np.searchsorted(idx, keep)
                assert abs(op).sum(axis=1).max() / abs(unit) < 200
                if len(key) == 3:
                    assert (
                        len(np.unique(keep // 4)) <= 4
                        and abs(op).sum(axis=1).max() <= 4
                    )
                block = regional.region.compress(
                    w[:, indices, :], op[keep, :][:, keep], unit
                )
                channels.setdefault(key, []).append(block)
    basis_order = [(q, name) for q in (-1, 1) for name in regional.NAMES] + list(
        cuts.CHANNELS
    )
    rows = []
    for key in basis_order:
        blocks = channels[key]
        error = sp.Rational("2e-15" if len(key) == 2 else "5e-19")
        projected = [
            cert.ri_band(cert.inner(h, b), 2 * error) for h in blocks for b in local
        ]
        contrast = cert.ri_band((blocks[0] - blocks[1]).trace(), 8 * error) / iv.sqrt(8)
        rows.append([cuts.survey.pack(v) for v in projected + [contrast]])
    return rows


def local_log_center():
    d = json.loads((cuts.survey.HERE.parent / "results.json").read_text())
    x = np.array(d["x_real"]) + 1j * np.array(d["x_imag"])
    blocks = [cert.prior.sym_matrix(b) for b in x]
    blocks = [(b + b.H) / 2 for b in blocks]
    mass = sum(b.trace() for b in blocks)
    local = local_basis()
    blocks_log = []
    traces = []
    for b in blocks:
        coords, err = cert.prior.entropy_coordinates(
            b / mass, 0, sp.Rational(".02"), local + [sp.eye(4)]
        )
        blocks_log.extend(cert.ri_band(v, 2 * err) for v in coords[:-1])
        traces.append(cert.ri_band(coords[-1], 4 * err))
    return np.array(blocks_log + [(traces[0] - traces[1]) / iv.sqrt(8)], dtype=object)


def certify(local_effects):
    previous = json.loads((HERE / "results.json").read_text())
    allcuts = json.loads((HERE / "effects.json").read_text())["cuts"]
    summary = json.loads((cuts.survey.HERE / "translation.json").read_text())["rows"]
    by_key = {tuple(r["channel"]): r["row"] for r in allcuts}
    b = cuts.interval_rows(summary + [by_key[k] for k in cuts.CHANNELS])
    c = cuts.interval_rows(local_effects)
    material = np.array(
        [cert.ri(w) for w in previous["spectrum_weights"]], dtype=object
    )
    a = local_log_center()
    # dJ = B dp + C dm; the marginal measurements read dm directly.
    marginal = [f"{v:.16g}" for v in cuts.mid(a - material @ c)]
    fitted = np.array([cert.ri(s) for s in marginal], dtype=object)
    residual = cuts.norm(
        np.r_[cuts.exact_entropy_center() - material @ b, a - fitted - material @ c]
    )
    bounds = json.loads((cuts.survey.HERE.parent / "certificate.json").read_text())
    radius = cert.ri(str(bounds["state_HS_error_upper"]))
    error = residual + radius / (cert.ri(".02") - radius)
    assert error.b < cert.ri("3e-5")
    gain = sum(abs(v) for v in np.r_[material, fitted])
    return {
        "trace_one_rank": 31,
        "survey_volume_tangent_dimension": 30,
        "material_weights": previous["spectrum_weights"],
        "marginal_and_parity_weights": marginal,
        "full_modular_error_upper": cuts.survey.pack(error)[1],
        "joint_modular_error_upper": previous["modular_error_upper"],
        "full_reader_l1_upper": int(float(gain.b)) + 1,
    }


def held_volume_arms(report):
    from scipy.optimize import brentq
    from ... import station
    from . import held

    x0 = json.loads((cuts.survey.HERE.parent / "results.json").read_text())
    state = np.array(x0["x_real"]) + 1j * np.array(x0["x_imag"])
    local = np.array(local_basis(), dtype=complex)
    parity = np.array([np.eye(4), -np.eye(4)]) / np.sqrt(8)
    joint = station.iso.basis(True)
    source = (0.3 * parity).astype(complex)
    source[0] += local[0] + 0.25 * joint[1]
    source[1] += 0.7 * local[5]
    adjust = np.zeros_like(source)
    adjust[0] = joint[-1]
    effects = station.radar_menus()[6]
    r = station.pair(effects, state)
    ds = station.pair(effects, source)
    dz = station.pair(effects, adjust)
    area, volume = station.geometry(r)
    grad = station.range_gradient(r)[1]
    direction = source - adjust * (grad @ ds) / (grad @ dz)
    direction /= np.linalg.norm(direction)
    dr = station.pair(effects, direction)
    step = 4e-5
    displacements = []
    areas = []
    volumes = []
    for sign in (-1, 1):
        base = r + sign * step * dr
        shift = brentq(
            lambda t: station.geometry(base + t * dz)[1] - volume,
            -1e-6,
            1e-6,
            xtol=1e-16,
        )
        delta = sign * step * direction + shift * adjust
        assert (
            np.linalg.norm(delta) < 1e-4
            and np.linalg.eigvalsh(state + delta).min() > 0.02
        )
        displacements.append(delta)
        area_arm, volume_arm = station.geometry(station.pair(effects, state + delta))
        areas.append(float(area_arm))
        volumes.append(float(volume_arm))
    minus, plus = displacements
    delta = (plus - minus) / (2 * step)
    dm = np.r_[
        np.einsum("aij,bji->ba", local, delta).real.ravel(),
        np.einsum("bij,bji->", parity, delta).real,
    ]
    assert abs(dm[0]) > 0.1 and abs(dm[11]) > 0.1 and abs(dm[-1]) > 0.05
    dmats = held.material_slopes(delta)
    mat = float(np.array(report["material_weights"], float) @ dmats)
    heads = float(np.array(report["marginal_and_parity_weights"], float) @ dm)
    spectral = float(
        (station.entropy(state + plus) - station.entropy(state + minus)) / (2 * step)
    )
    areal = -x0["ordinary_eta"] * (areas[1] - areas[0]) / (2 * step)
    assert max(abs(v - volume) for v in volumes) < 1e-13
    assert max(abs(mat + heads - spectral), abs(areal - spectral)) < 3e-5
    return dict(
        arm_step=step,
        volume_target=float(volume),
        arm_volumes=volumes,
        marginal_and_parity_slopes=dm.tolist(),
        material_slopes=dmats,
        weighted_material=mat,
        weighted_marginals=heads,
        combined=mat + heads,
        direct_entropy=spectral,
        ordinary_calibrated_area=areal,
    )


def write(rows, report):
    body = ",\n".join("    " + json.dumps(row, separators=(",", ":")) for row in rows)
    (HERE / "local_effects.json").write_text(
        '{\n  "scale": '
        + str(cuts.survey.SCALE)
        + ',\n  "rows": [\n'
        + body
        + "\n  ]\n}\n"
    )
    Path(__file__).with_suffix(".json").write_text(
        json.dumps({"scale": cuts.survey.SCALE, **report}, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    iv.dps = 36
    rows = native_local_effects()
    report = certify(rows)
    report["held"] = held_volume_arms(report)
    if args.write:
        write(rows, report)
    print(json.dumps(report))
