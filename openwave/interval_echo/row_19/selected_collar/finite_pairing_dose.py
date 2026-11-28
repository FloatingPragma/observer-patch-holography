"""Test a finite, energy-preserving change of joint records at matched survey volume."""

import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.full_algebra import enclosure, station
from matched_einstein_limit.source_implication.shared_radar.measure import preparations
from scipy.optimize import brentq

from . import regional_balance as regional
from . import returned_sources as source

lab = regional.lab
heads = lab.heads
iv = enclosure.iv
HERE = lab.HERE


def geometry(ranges):
    """Interval area and volume from twelve measured ranges, without derivatives."""
    logs = np.zeros((4, 4), dtype=object)
    for (a, b), r in zip(zip(*station.ORDER, strict=True), ranges, strict=True):
        logs[a, b] = iv.ln(r)
    gains = [sum(logs[a, b] - logs[b, a] for b in range(4)) / 4 for a in range(4)]
    square = np.zeros_like(logs)
    for a in range(4):
        for b in range(a + 1, 4):
            square[a, b] = square[b, a] = iv.exp(
                logs[a, b] + logs[b, a] - gains[a] - gains[b]
            )

    def dot(a, b, c):
        return (square[a, b] + square[a, c] - square[b, c]) / 2

    area = sum(
        iv.sqrt(square[i, j] * square[i, k] - dot(i, j, k) ** 2) / 2
        for i, j, k in combinations(range(4), 3)
    )
    d, e, f = square[0, 1:]
    g, h, k = dot(0, 1, 2), dot(0, 1, 3), dot(0, 2, 3)
    volume = iv.sqrt(d * e * f + 2 * g * h * k - d * k * k - e * h * h - f * g * g) / 6
    return area, volume


def source_read(tape):
    blocks, err = source.observed(tape, lab.SHOTS, False)
    ent = 0.0
    width = 0.0
    for x, e in zip(blocks, err, strict=True):
        levels = np.linalg.eigvalsh(x)
        floor = levels[0] - 2 * e
        assert floor > 0 and levels[-1] + 2 * e < 1
        ent -= levels @ np.log(levels)
        # Pauli coordinate bounds give HS error <=2e and trace norm <=4e.
        width += 4 * e * (1 + abs(np.log(floor)))
    return regional.band(ent - width - 1e-12, ent + width + 1e-12)


def rotations(X):
    margins = heads.measure.marginals(X)
    for block in range(2):
        for head, marg in enumerate(margins):
            v = np.array([np.trace(marg[block] @ a).real for a in heads.AXES[1:]])
            axis = (
                sum(z * a for z, a in zip(v, heads.AXES[1:], strict=True))
                / np.linalg.norm(v)
                if np.linalg.norm(v)
                else heads.AXES[-1]
            )
            for k in range(1, 9):
                for sign in (1, -1):
                    theta = sign * np.pi * k / 16
                    u = np.cos(theta) * np.eye(2) - 1j * np.sin(theta) * axis
                    gate = np.kron(u, np.eye(2)) if head == 0 else np.kron(np.eye(2), u)
                    out = X.copy()
                    out[block] = gate @ X[block] @ gate.conj().T
                    yield out


def correlated_write(X, effects, observed, reference_ranges):
    basis = heads.measure.joint_basis()
    entropy = -station.pair(basis, station.log_blocks(observed))
    v = station.range_gradient(reference_ranges)[1]
    volume = np.array([v @ station.pair(effects, d) for d in basis])
    weights = entropy / np.linalg.norm(entropy) + volume / np.linalg.norm(volume)
    d = np.tensordot(weights / np.linalg.norm(weights), basis, axes=1)
    for k in range(1, 30):
        step = 2.0**-k
        y = X + step * d
        if np.linalg.eigvalsh(y).min() <= 0 or station.entropy(y) <= station.entropy(X):
            continue
        if (
            station.geometry(station.pair(effects, y))[1]
            <= station.geometry(station.pair(effects, X))[1]
        ):
            continue
        p = station.preparation(y)
        if min(p["min_kept"], p["min_dephased"]) > 0 and p["max_error"] < 1e-11:
            return y, step, p
    raise ValueError("No realizable jointly increasing source")


def run():
    iv.dps = 35
    certificate = Path(enclosure.__file__).with_name("certificate.json")
    cert = json.loads(certificate.read_bytes())
    gain = iv.mpf(cert["gain_interval"])
    target = iv.mpf(cert["volume_interval"])
    # Also include the actual nominal reference's volume, bounding its offset
    # by the independently enclosed certificate residual.
    offset = cert["volume_offset_upper"]
    target += regional.band(-offset, offset)
    previous = json.loads((HERE / "returned_head_lab.json").read_bytes())["tapes"][
        "reference"
    ]
    corrected, _ = source.observed(previous["incoming"], lab.SHOTS, False)
    idx = [i for i, z in enumerate(lab.six.survey.ROUTES) if max(z) < 4]
    calibrated = np.asarray(previous["time"])[idx] @ regional.SCORE / lab.SHOTS
    X = preparations()["equilibrium"]
    P = heads.product_reference(X)
    _, bins = lab.six.survey.return_bins()
    bins = bins[idx]
    effects = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, regional.SCORE[:-1] - 3)
    rng = np.random.default_rng(28112081)

    def acquire(x):
        probs = station.pair(bins, x)
        return {
            "time": [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in probs],
            "incoming": lab.pairs(x, rng, False),
        }

    def ranges(tape):
        means, _ = regional.timing_slopes(
            {"x": tape["time"]}, lab.SHOTS, (("x", "x", 1),)
        )
        return np.array(means["x"], dtype=object)

    def R(y):
        return station.pair(effects, y)

    def V(r):
        return station.geometry(r)[1]

    x0, rp = R(X), R(P)
    rot = np.array([R(y) for y in rotations(X)])
    # Retain every attempted dose and the largest far-end volume in the menu.
    rotation_rows = []
    for k in range(1, 21):
        t = 2.0**-k
        near = V((1 - t) * x0 + t * rp) - V(x0)
        top = max(V((1 - t) * r + t * rp) - V(x0) for r in rot)
        rotation_rows.append(
            {"dose": t, "near_volume_offset": near, "maximum_rotated_offset": top}
        )
    assert all(
        r["near_volume_offset"] < 0 and r["maximum_rotated_offset"] < 0
        for r in rotation_rows
    )
    W, step, recipe = correlated_write(X, effects, corrected, calibrated)
    inputs = {"reference": X, "independent": P, "compensator": W}
    for state in inputs.values():
        np.testing.assert_allclose(
            heads.measure.marginals(state),
            heads.measure.marginals(X),
            atol=1e-14,
            rtol=0,
        )
    tapes = {name: acquire(x) for name, x in inputs.items()}
    readings = {name: ranges(t) for name, t in tapes.items()}
    entropy = {name: source_read(t["incoming"]) for name, t in tapes.items()}
    assert entropy["compensator"].a > entropy["reference"].b
    doses = []
    for k in range(1, 21):
        t = 2.0**-k

        def r(u):
            return (
                (1 - t - u) * readings["reference"]
                + t * readings["independent"]
                + u * readings["compensator"]
            )

        v0 = geometry(r(0))[1] - target
        v1 = geometry(r(1 - t))[1] - target
        doses.append(
            {
                "dose": t,
                "near_volume_offset": lab.six.bounds(v0),
                "far_volume_offset": lab.six.bounds(v1),
            }
        )
        if v0.b < 0 < v1.a:
            break
    else:
        raise ValueError("No volume match")

    # A single derivative enclosure covers all mixtures at the selected dose.
    def hull(x, y):
        return regional.band(
            min(float(x.a), float(y.a)) - 1e-14, max(float(x.b), float(y.b)) + 1e-14
        )

    gradients = [
        j.g
        for j in enclosure.geometry_jet(
            [hull(x, y) for x, y in zip(r(0), r(1 - t), strict=True)]
        )
    ]

    def tight(box, shape):
        mid = np.asarray([iv.mpf(float(x.mid)) for x in box], dtype=object)
        return geometry(mid)[shape] + sum(gradients[shape] * (box - mid))

    target_mid = float(target.mid)
    u = brentq(
        lambda z: V(np.array([float(q.mid) for q in r(z)])) - target_mid, 0, 1 - t
    )
    # Enclose at least one exact matching fraction. Uniqueness is unnecessary.
    for spread in (2.0**-k for k in range(16, -1, -1)):
        a, b = max(0, u - spread), min(1 - t, u + spread)
        vl, vh = tight(r(a), 1) - target, tight(r(b), 1) - target
        if vl.b < 0 < vh.a:
            break
    assert vl.b < 0 < vh.a
    held = (1 - t - u) * X + t * P + u * W
    tapes["held"] = acquire(held)
    held_entropy = source_read(tapes["held"]["incoming"])
    measured = geometry(ranges(tapes["held"]))
    area0 = geometry(readings["reference"])[0]
    # Along the segment from the counted reference to the certified optimizer,
    # trace norm <= sqrt(8) HS radius. Eigenvalues remain above this known floor.
    hs = cert["state_HS_error_upper"]
    floor = np.linalg.eigvalsh(X).min() - hs - 1e-12
    # The difference has trace zero. Center -log(X) on the identity and
    # enclose the log derivative along the segment by 1/minimum eigenvalue.
    logs = np.log(np.linalg.eigvalsh(X))
    reference_slack = (
        hs * np.linalg.norm(logs - logs.mean()) + hs * hs / (2 * floor) + 1e-12
    )
    gainP = entropy["independent"] - entropy["reference"]
    gainW = entropy["compensator"] - entropy["reference"]
    lower = t * float(gainP.a) + a * float(gainW.a) - reference_slack
    assert lower > 0
    predicted = lower / float(gain.b) - 1e-12
    # At the same volume as the actual preparation, the Lagrange volume
    # terms cancel. Its objective gap to the exact global maximum is
    # second order in hs because the full Lagrangian is stationary there.
    curvature = cert["global_survey"]["trace_norm_curvature_limit"]
    objective_slack = (1 / floor + 8 * curvature) * hs**2 / 2 + 1e-14
    matched_prediction = (
        t * float(gainP.a) + a * float(gainW.a) - objective_slack
    ) / float(gain.b) - 1e-12
    assert matched_prediction > 0
    # With equal one-head marginals the log(product) price pairs to zero.
    # Read the frozen-reference nonadditive modular charge for a finite step.
    Xread, ex = source.observed(tapes["reference"]["incoming"], lab.SHOTS, False)

    def modular_gain(name):
        y, ey = source.observed(tapes[name]["incoming"], lab.SHOTS, False)
        q = 0.0
        err = 0.0
        for B, Y, e, f in zip(Xread, y, ex, ey, strict=True):
            levels, U = np.linalg.eigh(B)
            K = (U * (-np.log(levels))) @ U.conj().T
            diff = Y - B
            coeff = abs(np.einsum("ijkl,lk->ij", heads.PRODUCTS, K)).sum() / 4
            err += coeff * (e + f) + (
                np.linalg.norm(diff, ord="nuc") + 4 * (e + f)
            ) * 2 * e / (levels[0] - 2 * e)
            q += np.trace(diff @ K).real
        return regional.band(q - err - 1e-12, q + err + 1e-12)

    margin_families = {}
    for name, tape in tapes.items():
        y, ey = source.observed(tape["incoming"], lab.SHOTS, False)
        differences = np.einsum("ijkl,blk->bij", heads.PRODUCTS, y - Xread).real
        error = (ex + ey)[:, None]
        assert np.all(abs(differences[:, :, 0]) <= error) and np.all(
            abs(differences[:, 0, :]) <= error
        )
        margin_families[name] = float(
            max(abs(differences[:, :, 0]).max(), abs(differences[:, 0, :]).max())
        )
    qheld = modular_gain("held")
    qroot = t * modular_gain("independent") + regional.band(a, b) * modular_gain(
        "compensator"
    )
    # Independently enclose A(X*)-A(X) using the range-effect spectral widths
    # and the same certified trace distance, including the native propagation error.
    diam = np.ptp(np.linalg.eigvalsh(effects), axis=(1, 2))
    ce = cert["range_effect_error_upper"]
    widths = (diam / 2 + 2 * ce + 1e-12) * np.sqrt(8) * hs + 2 * ce + 1e-12
    ref_box = [
        regional.band(z - w - 1e-12, z + w + 1e-12)
        for z, w in zip(x0, widths, strict=True)
    ]
    grad = enclosure.geometry_jet(ref_box)[0].g
    gmid = np.array([float(z.mid) for z in grad])
    gw = np.array([float(abs(z - c).b) for z, c in zip(grad, gmid, strict=True)])
    # The same trace-zero displacement perturbs all clocks together.
    centered = effects - (np.trace(effects, axis1=-2, axis2=-1).real.sum(1) / 8)[
        :, None, None, None
    ] * np.eye(4)
    changes = (
        np.linalg.norm(np.tensordot(gmid, centered, axes=1))
        + gw @ np.linalg.norm(centered.reshape(12, -1), axis=1)
    ) * hs
    changes += (abs(gmid) + gw).sum() * 3e-12 + 1e-12
    reference_area = geometry([iv.mpf(float(z)) for z in x0])[0]
    ref_exact = reference_area + regional.band(-changes, changes)
    gain_observed = held_entropy - entropy["reference"]
    assert (measured[0] - area0).b < 0 and gain_observed.a > 0

    # Keep the same mixture fraction for all twelve directed clocks.
    mid = (a + b) / 2
    d_range = readings["compensator"] - readings["reference"]
    root_area = tight(r(mid), 0) + regional.band(a - mid, b - mid) * sum(
        gradients[0] * d_range
    )
    return {
        "input_sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in [
                certificate,
                HERE / "returned_head_lab.json",
                HERE / "FINITE_PAIRING_DOSE.md",
                HERE / "PAIRING_VOLUME_COMPENSATOR.md",
            ]
        },
        "seed": 28112081,
        "shots_per_setting": lab.SHOTS,
        "rotations": rotation_rows,
        "write_step": step,
        "rail_preparation": recipe,
        "source_recipes": {name: station.preparation(x) for name, x in inputs.items()},
        "doses": doses,
        "root_fraction_interval": [a, b],
        "root_endpoint_volume_offsets": list(map(lab.six.bounds, (vl, vh))),
        "selected_fraction": u,
        "reference_entropy_allowance": reference_slack,
        "entropy_gain_product": lab.six.bounds(gainP),
        "entropy_gain_compensator": lab.six.bounds(gainW),
        "minimum_entropy_gain_at_exact_volume": lower,
        "minimum_area_contraction_at_exact_volume": predicted,
        "radar_area_change_at_exact_volume": lab.six.bounds(root_area - ref_exact),
        "held_area_change": lab.six.bounds(measured[0] - area0),
        "held_volume_offset": lab.six.bounds(measured[1] - target),
        "held_entropy_gain": lab.six.bounds(gain_observed),
        "head_parity_max_measured_marginal_difference": margin_families,
        "modular_gain_at_bracketed_volume": lab.six.bounds(qroot),
        "held_connected_modular_gain": lab.six.bounds(qheld),
        "held_relative_entropy_to_original": lab.six.bounds(qheld - gain_observed),
        "reference_lagrangian_gap_upper": objective_slack,
        "minimum_area_contraction_same_prepared_volume": matched_prediction,
        "radar_area_change_same_prepared_volume": lab.six.bounds(root_area - area0),
        "ordinary_parked_energy_change": 0,
        "tapes": tapes,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    dest = HERE / "finite_pairing_dose.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        stored = json.loads(dest.read_bytes())
        assert (
            result["tapes"] == stored["tapes"]
            and result["input_sha256"] == stored["input_sha256"]
        )
        from cosmology.rar_deep_regime.lensing_profile import check

        check(stored, result)
    print(
        {
            k: v
            for k, v in result.items()
            if k
            in (
                "root_fraction_interval",
                "minimum_area_contraction_at_exact_volume",
                "radar_area_change_at_exact_volume",
                "held_area_change",
                "held_volume_offset",
                "held_entropy_gain",
            )
        }
    )
