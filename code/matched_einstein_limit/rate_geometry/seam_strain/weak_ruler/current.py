"""Signed conserved pump energy and number crossing each physical port cut."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext

from . import aperture, regional, threshold
from .. import live_modular, live_radar, live_source
from ..source import certificate as cert, measure as source

HERE = Path(__file__).parent


def compressed(cross, inside, hi, lo, denominator):
    """i (V_i† M_ie V_e − h.c.) in the eight-dimensional input space."""
    cut = cross.tocoo()
    assert not set(inside).intersection(cut.col)

    def entries(ids):
        x = hi[ids].astype(object) * (1 << 36) + lo[ids]
        return np.split(x, 2, axis=1)

    re, im = entries(inside)
    er, ei = entries(cut.col)
    left = np.zeros_like(re)
    right = np.zeros_like(im)
    np.add.at(left, cut.row, cut.data[:, None] * er)
    np.add.at(right, cut.row, cut.data[:, None] * ei)
    d = sp.Matrix(re.T @ left + im.T @ right) + sp.I * sp.Matrix(
        re.T @ right - im.T @ left
    )
    return (sp.I * (d - d.H)).applyfunc(sp.expand) / (denominator * (1 << 164))


def current_effects(ports):
    board = source.formation.federation(1)
    h = source.system(board, source.formation.region(board), -1)[0]
    inside = live_radar.station_rows(board, ports).ravel()
    keep = np.ones(h.shape[0], dtype=np.int64)
    keep[inside] = 0
    first = h[inside]
    hi, lo, colerr, _ = live_source.pump_center()
    norm = int(abs(h).sum(axis=1).max()) / sp.Integer(40)
    matrices = {}
    # H_R=(P H+H P)/2 assigns half of each boundary bond. The outward
    # currents are -i[H,P] and -i[H^2,P]/2 in exchange units.
    for k, m, den, bound in (
        ("number", first, 40, 2 * norm),
        ("energy", first @ h, 3200, norm**2),
    ):
        cross = m.multiply(keep)
        cross.eliminate_zeros()
        matrices[k] = (
            compressed(cross, inside, hi, lo, den),
            bound,
            (6 * colerr + 9 * colerr**2) * bound,
        )
    return matrices


def slopes(ports):
    (formed, product), ferr, _ = live_source.formation_center()
    right, err = live_modular.rotation_centers(formed, ferr)["right"]
    inputs = [("decorrelated", product - formed, 3 * ferr, 2), ("right", right, err, 1)]
    inputs.extend(
        (a + b, x, 0, 1)
        for (a, b), x in zip(
            ((a, b) for a in "XYZ" for b in "XYZ"),
            live_source.qubit_hidden.correlation_axes(),
            strict=True,
        )
    )
    out = {}
    for k, (effect, limit, flight_error) in current_effects(ports).items():
        row = {}
        for name, x, prep_error, norm in inputs:
            sigma = sp.kronecker_product(sp.diag(2, 1) / 3, x)
            v = sp.re(sp.trace(effect * sigma))
            row[name] = cert.enclosed(
                v, limit * prep_error + (norm + prep_error) * flight_error
            )
        out[k] = row
    return out


def clock_transfer(data, bank):
    """The same nine time bins as the entropy read, calibrated against native energy flow."""
    target = [data["energy"][k] for k in bank["axes"]]
    m = np.array(bank["response_matrix"])
    coefficients = np.linalg.solve(m.T, np.mean(target, axis=1))
    w = list(map(Fraction, coefficients))
    mid = [[Fraction(x) for x in r] for r in m]
    scale = sum(abs(x) for x in w)
    residual = max(
        abs(sum(w[j] * mid[j][i] for j in range(9)) - Fraction(y))
        for i, r in enumerate(target)
        for y in r
    )
    residual += scale * Fraction(bank["entry_error_upper"])
    trials = bank["held"]
    reads = {}
    for k, arm in trials["positive_preparations"].items():
        coords = np.array(arm["true_input_coordinates"])
        diff = (
            np.array(trials["probabilities"][k]) - trials["probabilities"]["reference"]
        )
        reads[k] = dict(
            direct=float(np.mean(target, axis=1) @ coords),
            from_clock_probabilities=float(coefficients @ diff),
        )
    # Each fixed-marginal coefficient is its Pauli expectation contrast, hence
    # at most two in magnitude for any two physical preparations. This crude
    # uniform bound also applies without knowing which held state was prepared.
    return dict(
        weights=coefficients.tolist(),
        coefficient_error_upper=np.nextafter(float(residual), np.inf),
        uniform_flux_bias_upper=np.nextafter(18 * float(residual), np.inf),
        held=reads,
    )


def universal_current_gap(rows):
    """Lower bound for max_R |ΔJ_E(R) + c ΔJ_N(R)| for every constant c."""
    left, right = rows
    n1, n2 = (
        Fraction(r["matched_front_lottery_difference"]["number"][0]) for r in rows
    )
    z1 = Fraction(left["front_energy_origin_for_zero_response"][1])
    z2 = Fraction(right["front_energy_origin_for_zero_response"][0])
    assert min(n1, n2) > 0 and z2 > z1
    # For positive n_i, a bound e on both residuals would bound the
    # distance between their two zeros by e*(1/n1+1/n2).
    limit = (z2 - z1) / (1 / n1 + 1 / n2)
    assert limit > Fraction(34, 10**12)
    return dict(
        uniform_over_energy_origin_lower=3.4e-11,
        best_possible_common_prediction_error_lower=1.7e-11,
    )


def run():
    board = source.formation.federation(1)
    fronts = json.loads((HERE / "observability.json").read_text())
    ctx = MPIntervalContext()
    ctx.dps = 40
    dz = threshold.z_intervals()
    alpha = -ctx.mpf(dz["decorrelated"]) / ctx.mpf(dz["right"])
    regions = []
    for rad in (1, 2):
        ports = aperture.enroll(board, rad)[0]
        data = slopes(ports)
        local = {
            k: aperture.encode(ctx.mpf(d["decorrelated"]) + alpha * ctx.mpf(d["right"]))
            for k, d in data.items()
        }
        lotteries = {k: regional.lottery_response(d, fronts) for k, d in data.items()}
        row = dict(
            radius=rad,
            ports=ports,
            derivatives=data,
            fixed_one_hop_metric=local,
            matched_front_lottery_difference=lotteries,
        )
        for name, d in (("local", local), ("front", lotteries)):
            num, energy = ctx.mpf(d["number"]), ctx.mpf(d["energy"])
            if num.a * num.b > 0:
                row[name + "_energy_origin_for_zero_response"] = aperture.encode(
                    -energy / num
                )
        bank = json.loads((HERE / "late_read.json").read_text())
        row["independent_clock_flux"] = clock_transfer(data, bank)
        regions.append(row)
    return dict(
        energy_assignment="half of each bond at each endpoint; absolute pump generator 40H from source.system",
        sign="positive is outward flow at the four-tick pump appointment",
        regions=regions,
        front_matched_current_gap=universal_current_gap(regions),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "current.json").write_text(data)
    else:
        print(data, end="")
