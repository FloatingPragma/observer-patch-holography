"""Infer finite source and geometric changes from terminal positions of free probes."""

import argparse
import json
from hashlib import sha256
from os.path import relpath
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.full_algebra.enclosure import (
    geometry_jet,
)

from . import finite_dose_curvature as c

lab, iv, reg = c.lab, c.iv, c.reg
HERE = lab.HERE
AXES = c.dose.heads.PRODUCTS.reshape(16, 4, 4) / 2
ORIGIN = np.stack([np.eye(4) / 8] * 2)
BASIS = np.zeros((31, 2, 4, 4), complex)
BASIS[:15, 0] = AXES[1:]
BASIS[15:30, 1] = AXES[1:]
BASIS[-1] = np.stack([AXES[0], -AXES[0]]) / np.sqrt(2)
# Each bin probability is computed from one certified four-by-four kernel.
# This allowance includes conversion, Hermitian products and Pauli projection.
EFFECT_ERROR = 5e-14
SEED = 28112101


def coefficients(ops):
    return np.einsum("...bij,kbji->...k", ops, BASIS).real


def calibrate():
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    ports = (*region, 5, 6)
    columns = np.array([4 * p + i for p in ports for i in range(4)])
    weights = native.prior.weights()
    banks = []
    for sign in (-1, 1):
        by_block = []
        for odd in (False, True):
            h = native.generator(board, region, sign, odd)
            exact = native.kernel_panels(h, columns, weights)
            unitary = (
                exact[0].astype(float) + 1j * exact[1].astype(float)
            ) / native.prior.FRAC
            # Source port, observation tick, destination, final and initial heads.
            u = unitary[:11].transpose(2, 0, 1, 3, 4)
            by_block.append(u.conj().swapaxes(-1, -2) @ u)
        banks.append(np.stack(by_block, axis=3))
    effects = np.asarray(banks).reshape(-1, len(ports), 2, 4, 4)
    return ports, effects


def emit(probabilities, shots, rng):
    records = []
    batch = 2**52
    assert shots >= batch and shots % batch == 0
    for p in probabilities:
        full = np.append(p, 1 - p.sum())
        assert np.min(full) > -EFFECT_ERROR
        order = np.argsort(full, kind="stable")
        law = np.maximum(full[order], 0)
        law[-1] = 1 - law[:-1].sum()
        counts = np.zeros(len(law), dtype=object)
        batches = shots // batch
        for offset in range(0, batches, 1024):
            draws = rng.multinomial(batch, law, size=min(batches - offset, 1024))
            counts += draws.sum(axis=0).astype(object)
        records.append(counts[np.argsort(order)].tolist())
    return records


class PositionReader:
    def __init__(self, effects):
        flat = effects.reshape(-1, 2, 4, 4)
        self.baseline = c.station.pair(flat, ORIGIN)
        self.matrix = coefficients(flat)
        noise = np.sqrt(np.maximum(self.baseline, 1e-9))
        weighted = self.matrix / noise[:, None]
        self.inverse = np.linalg.pinv(weighted) / noise[None, :]
        product = self.inverse @ self.matrix
        rounding = 3e-13 * (abs(self.inverse) @ abs(self.matrix)) + 3e-14
        self.defect = float(np.linalg.norm(abs(np.eye(31) - product) + rounding))
        assert self.defect < 1
        self.singular_lower = (1 - self.defect) / np.linalg.norm(self.inverse)
        # Conversion to the true response matrix costs at most sqrt(m*31)
        # times the outward per-entry allowance in spectral norm.
        self.singular_lower -= np.sqrt(self.matrix.size) * EFFECT_ERROR
        assert self.singular_lower > 0

    def decode(self, records, shots):
        counts = np.array(records, dtype=object)
        if (
            counts.shape != (len(self.baseline) // 6, 7)
            or np.min(counts) < 0
            or not np.all(counts.sum(1) == shots)
        ):
            raise ValueError(
                "Every enrolled launch requires all six destinations and the complement"
            )
        freq = np.array(counts[:, :-1], dtype=float).ravel() / shots
        # The union bound covers every frequency of all four sources and exposures.
        width = reg.src.probability_error(freq, shots) + EFFECT_ERROR
        vector = self.inverse @ (freq - self.baseline)
        estimate = ORIGIN + np.tensordot(vector, BASIS, axes=1)
        residual = freq - self.baseline - self.matrix @ vector
        exact_residual = (
            abs(self.inverse @ residual)
            + 3e-13 * (abs(self.inverse) @ (abs(residual) + 1e-2))
            + 2e-13
        )
        radius = (
            np.linalg.norm(abs(self.inverse) @ width + exact_residual) + 1e-13
        ) / (1 - self.defect)
        return {
            "vector": vector,
            "matrix": estimate,
            "radius": radius,
            "width": width,
            "residual": residual,
        }

    def error(self, obs, price):
        # Bound Tr(price*(true-estimated)) using the same probability event,
        # including when the price itself was chosen from the measured records.
        k = coefficients(price)
        w = k @ self.inverse
        w_round = 1e-14 * (abs(k) @ abs(self.inverse) + 1)
        defect = (
            abs(k - w @ self.matrix)
            + 3e-13 * (abs(w) @ abs(self.matrix))
            + w_round @ abs(self.matrix)
            + 2e-13
        )
        r = obs["residual"]
        return (
            abs(w @ r)
            + (abs(w) + w_round) @ obs["width"]
            + (w_round + 3e-13 * abs(w)) @ (abs(r) + 1e-2)
            + np.linalg.norm(defect) * obs["radius"]
            + 2e-12
        )

    def change(self, x, y, price):
        center = c.station.pair(price, y["matrix"] - x["matrix"])
        err = self.error(x, price) + self.error(y, price)
        return reg.band(center - err, center + err)


def entropy(reader, z):
    mat = z["matrix"]
    rad = z["radius"]
    floor = np.linalg.eigvalsh(mat).min() - rad - 2e-12
    if floor <= 0:
        return None
    slope = reader.error(z, -c.station.log_blocks(mat))
    center = c.station.entropy(mat)
    return reg.band(
        center - slope - rad**2 / (2 * floor) - 2e-12, center + slope + 2e-12
    )


def joint_charge(reader, x, y):
    mat = x["matrix"]
    rad = x["radius"]
    left, right = c.dose.heads.measure.marginals(mat)
    levels = [
        (np.linalg.eigvalsh(mat).min(), rad),
        (np.linalg.eigvalsh(left).min(), np.sqrt(2) * rad),
        (np.linalg.eigvalsh(right).min(), np.sqrt(2) * rad),
        (np.trace(mat, axis1=1, axis2=2).real.min(), 2 * rad),
    ]
    if any(f - e <= 2e-12 for f, e in levels):
        return None
    log_error = sum(e / (f - e - 2e-12) for f, e in levels) + 1e-11
    trace_dist = np.sqrt(8) * (
        np.linalg.norm(y["matrix"] - mat) + rad + y["radius"] + 1e-12
    )
    prod = c.dose.heads.product_reference(mat)
    K = c.station.log_blocks(prod) - c.station.log_blocks(mat)
    return reader.change(x, y, K) + reg.band(-1, 1) * trace_dist * log_error


def survey_changes(reader, x, y, effects):
    xr = c.station.pair(effects, x["matrix"])
    yr = c.station.pair(effects, y["matrix"])
    ex = np.array([reader.error(x, E) for E in effects]) + 1e-13
    ey = np.array([reader.error(y, E) for E in effects]) + 1e-13
    box = np.array(
        [
            reg.band(min(a - u, b - v), max(a + u, b + v))
            for a, b, u, v in zip(xr, yr, ex, ey, strict=True)
        ]
    )
    outputs = []
    for jet in geometry_jet(box):
        gradient = np.array([float(v.mid) for v in jet.g])
        price = np.tensordot(gradient, effects, axes=1)
        spread = np.array(
            [float(abs(v - w).b) for v, w in zip(jet.g, gradient, strict=True)]
        )
        error = spread @ (abs(yr - xr) + ex + ey) + 1e-12 * (1 + abs(gradient).sum())
        outputs.append(reader.change(x, y, price) + reg.band(-error, error))
    return outputs


def analyze(reader, tapes, shots, effects, gain, multiplier, gap):
    states = {k: reader.decode(t, shots) for k, t in tapes.items()}
    x = states["reference"]
    sx = entropy(reader, x)
    result = {}
    for name, y in states.items():
        if name == "reference":
            continue
        sy = entropy(reader, y)
        ds = sy - sx if sy is not None and sx is not None else None
        area, volume = survey_changes(reader, x, y, effects)
        q = joint_charge(reader, x, y)
        bound = (
            None if ds is None else float(((ds - multiplier * volume - gap) / gain).a)
        )
        result[name] = {
            "reference_HS_radius": x["radius"],
            "held_HS_radius": y["radius"],
            "reference_smallest_eigenvalue": float(
                np.linalg.eigvalsh(x["matrix"]).min()
            ),
            "entropy_change": None if ds is None else lab.six.bounds(ds),
            "connected_charge_change": None if q is None else lab.six.bounds(q),
            "parked_energy_change": lab.six.bounds(
                reader.change(x, y, np.stack([np.zeros((4, 4)), np.eye(4)]))
            ),
            "inferred_area_change": lab.six.bounds(area),
            "inferred_volume_change": lab.six.bounds(volume),
            "minimum_contraction_from_entropy": bound,
        }
    return result


def run():
    iv.dps = 36
    prior, states, *_ = c.trials()
    x, y = states["reference"], states["held"]
    reverse = 2 * x - y
    recipe = c.station.preparation(reverse)
    assert (
        min(recipe["min_kept"], recipe["min_dephased"]) > 0
        and recipe["max_error"] < 1e-11
    )
    inputs = {"reference": x, "held": y, "reverse": reverse, "sham": x}
    ports, projectors = calibrate()
    reader = PositionReader(projectors)
    _, bins = lab.six.survey.return_bins()
    ids = [i for i, z in enumerate(lab.six.survey.ROUTES) if max(z) < 4]
    ranges = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins[ids], reg.SCORE[:-1] - 3)
    source = HERE / "link_source_contraction.json"
    old = json.loads(source.read_bytes())
    cert_file = Path(c.dose.enclosure.__file__).with_name("certificate.json")
    cert = json.loads(cert_file.read_bytes())
    gain = iv.mpf(cert["gain_interval"])
    multiplier = iv.mpf(old["ordinary_volume_multiplier"])
    prob = {k: c.station.pair(projectors, v) for k, v in inputs.items()}
    results = {}
    rng = np.random.default_rng(SEED)
    for power in (52, 60, 66, 70):
        n = 2**power
        tapes = {k: emit(v, n, rng) for k, v in prob.items()}
        rows = analyze(
            reader,
            tapes,
            n,
            ranges,
            gain,
            multiplier,
            prior["reference_lagrangian_gap_upper"],
        )
        results[str(power)] = {"records": tapes, "decoded": rows}
        if (
            power >= 66
            and (rows["held"]["minimum_contraction_from_entropy"] or 0) > 0
            and (rows["held"]["connected_charge_change"] or [-1])[0] > 0
        ):
            break
    return {
        "seed": SEED,
        "ports": list(ports),
        "ticks": list(range(1, 12)),
        "probe_signs": [-1, 1],
        "coordinate_count": 31,
        "left_inverse_defect_upper": reader.defect,
        "material_response_sigma_lower": reader.singular_lower,
        "exposures": results,
        "reverse_rail_lottery": recipe,
        "independent_held_radar_area": old["independent_area_change"],
        "independent_held_head_entropy": prior["held_entropy_gain"],
        "independent_held_head_charge": prior["held_connected_modular_gain"],
        "inputs": {
            relpath(f, HERE): sha256(f.read_bytes()).hexdigest()
            for f in (
                source,
                HERE / "finite_pairing_dose.json",
                HERE / "POSITION_SOURCE_PROTOCOL.md",
                cert_file,
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = HERE / "position_source.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(result["material_response_sigma_lower"])
    for power, data in result["exposures"].items():
        print(power, data["decoded"]["held"])
