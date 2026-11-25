"""Compare every registered regional translation summary to the stationary entropy reader."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from ...directional_stress import bridge
from ...directional_stress.traction.certificate import four_tick_columns
from ...directional_stress.translation import certificate as region
from ...directional_stress.translation import measure as translation
from .. import certificate as cert
from . import measure

NAMES = ("shift", "domain_cut", "image_cut", "marked_force")


def material_rows():
    atlas = translation.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    launch = bridge.beam.indices(sum(bridge.beam.screens(board, route), []))
    patch = translation.chosen(atlas)
    weights = bridge.native.prior.weights()[7]
    pauli = [cert.prior.sym_matrix(b) for b in cert.iso.basis(True)]
    rows = []
    for charge in (-1, 1):
        blocks = {name: [] for name in NAMES}
        for odd in (False, True):
            h, ops = translation.op_table(board, route, patch, charge, odd)
            idx = np.unique(
                np.concatenate([ops[name][0].nonzero()[0] for name in NAMES])
            )
            w = four_tick_columns(h, idx, launch, weights)
            for name in NAMES:
                mat, unit = ops[name]
                assert abs(mat).sum(axis=1).max() / abs(unit) < 200
                block = region.compress(w, mat[idx, :][:, idx], unit)
                blocks[name].extend(cert.inner(block, axis) for axis in pauli)
        # The regional propagation bound is 2e-15 in operator norm per head
        # block. Every normalized Pauli product has trace norm two.
        for name in NAMES:
            rows.append(
                [
                    measure.pack(cert.ri_band(v, sp.Rational("4e-15")))
                    for v in blocks[name]
                ]
            )
    return rows


def certify(rows):
    iv = measure.iv
    source = json.loads((measure.HERE / "results.json").read_text())
    area, volume = [measure.unpack(source[key]).ravel() for key in ("area", "volume")]
    material = [[iv.mpf([lo, hi]) / measure.SCALE for lo, hi in row] for row in rows]
    constraints = np.vstack([volume, material])
    # Use central floating-point values only to propose a right inverse,
    # tangent, and linear combination. All bounds use the true interval rows.
    midpoint = np.vectorize(lambda x: float(x.mid), otypes=[float])
    norms = np.linalg.norm(midpoint(constraints), axis=1)
    exponents = [int(-np.floor(np.log2(n))) for n in norms]
    scales = [iv.mpf(2) ** exp for exp in exponents]
    matrix = constraints * np.array(scales, dtype=object)[:, None]
    center = midpoint(matrix)
    a = midpoint(area)
    inv = np.linalg.pinv(center)
    price = a @ inv
    residual = a - price @ center
    d = -residual / np.linalg.norm(residual)
    fixed = np.vectorize(lambda x: cert.ri(f"{x:.15g}"), otypes=[object])
    k, coefficients, d = fixed(inv), fixed(price), fixed(d)

    def norm(v):
        return iv.sqrt(sum(x**2 for x in v.flat))

    defect = norm(np.eye(len(matrix), dtype=int) - matrix @ k)
    if defect.b >= 1:
        raise ValueError(
            f"translation constraints have uncertified rank, defect={defect}"
        )
    correction = norm(k) * norm(matrix @ d) / (1 - defect)
    # D=d-K(BK)^(-1)Bd obeys every constraint; subtracting coefficients
    # times B from the area gradient keeps its exact contraction on D.
    score = area - coefficients @ matrix
    slope = -sum(score * d)
    length = norm(d)
    unit = (slope + iv.mpf([-1, 1]) * norm(score) * correction) / (
        length + iv.mpf([-1, 1]) * correction
    )
    gain = iv.mpf(["23.1042361", "23.1042363"])
    entropy = gain * unit
    if entropy.a <= 0:
        raise ValueError(
            f"no signed separation from regional translation, interval={entropy}"
        )
    assert entropy.a > cert.ri(".001593") and entropy.b < cert.ri(".001595")
    return dict(
        constraint_rank=len(matrix),
        constraint_row_scale_powers_of_two=exponents,
        kernel_dimension=len(d) - len(matrix),
        tangent_center=[measure.pack(v) for v in d],
        inverse_defect=measure.pack(defect),
        tangent_correction=measure.pack(correction),
        signed_area_loss_per_unit_HS=measure.pack(unit),
        signed_entropy_gain_per_unit_HS=measure.pack(entropy),
    )


def write(rows, report):
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(
            {"scale": measure.SCALE, "rows": rows, "order": list(NAMES), **report},
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    measure.iv.dps = 32
    rows = material_rows()
    report = certify(rows)
    if args.write:
        write(rows, report)
    print(json.dumps(report))
