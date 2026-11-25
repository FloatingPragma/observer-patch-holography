"""Equal traction change on both species at fixed radar and fixed complete marked forces."""

import argparse
import json
from math import ceil

import numpy as np
import sympy as sp
from mpmath import iv

from ...full_algebra.reuse import invisible
from ..connection import force_entropy, force_gradients
from . import measure

c = measure.c


def solve(rows, err, target=None):
    raw = np.array(rows, dtype=float)
    scales = np.linalg.norm(raw, axis=1)
    normalized = raw / scales[:, None]
    k = invisible.dyadic(np.linalg.pinv(normalized) / scales[None, :])
    norm = invisible.frobenius
    defect = norm(sp.eye(rows.rows) - rows * k) + err * norm(k)
    assert defect.b < c.ri(".01"), str(defect)
    if target is None:
        target = sp.zeros(rows.rows, 1)
        target[15] = -2
        target[19] = -2
    assert 22 <= rows.rows <= 30 and rows.cols == 30
    # Rows: 12 geometry, 2 marked forces, then 4 dimer effects per species.
    # X=-2 and occupation fixed gives unit whole-seam strain work on both.
    base = k * target
    scale = ceil(float(norm(base).b))
    d = base / scale
    # Solve for true D by correcting its actual residual (the nominal inverse
    # is already accurate on this target, sharper than ||K|| times ||target||).
    adjustment = (
        norm(k) * (norm(target / scale - rows * d) + err * norm(d)) / (1 - defect)
    )
    size = norm(d) + adjustment
    assert size.b < c.ri("1.01")
    return d, scale, adjustment, defect, size


def fair_region(out, regional):
    from ...action_response.four_ticks import parse_interval

    native = c.iso.ranging
    words = [(1, 1)] + [word for q in (-1, 1) for word in ((1 + q, 1), (1, 1 + q))]
    columns = [set(native.formation.emission.columns(w)) for w in words]
    assert len(set.union(*columns)) == sum(map(len, columns))
    for q in (-1, 1):
        for sign in (-1, 1):
            # Native verifier builds the basis in the three stated load words
            # and checks invariance under the full 99-column marked action.
            native.parity_blocks(q, 8, sign)
    spectra = [
        [parse_interval(v) for v in regional[q]["eigenvalues"]] for q in ("-1", "1")
    ]
    eigenvalues = [(spectra[0][0] + spectra[1][0]) / 2] + [
        v / 2 for row in spectra for v in row[1:]
    ]
    assert all(v.a > 0 for v in eigenvalues)
    slope = sum(parse_interval(out[q]["entropy_derivative"]) for q in ("-1", "1")) / 2
    lower = (
        sum(
            parse_interval(out[q]["finite_entropy_increase_at_0p01_lower"])
            for q in ("-1", "1")
        )
        / 2
    )
    # Region has one shared outside word. Along this exact control both outside
    # probabilities are individually fixed, so the pooling correction is constant.
    return dict(
        load_words=[list(w) for w in words],
        reference_spectrum=list(map(str, eigenvalues)),
        signed_entropy_slope=str(slope),
        finite_entropy_change_at_0p01=str(iv.mpf([lower.a, (c.ri(".01") * slope).b])),
    )


def constraint_rows():
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    means, mean_error = force_entropy.range_gradients(frame)
    force = sp.Matrix(
        [
            [c.inner(h, b) for h in pair for b in frame]
            for pair in force_gradients.force_effects()
        ]
    )
    rows = means.col_join(force)
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    weights = measure.bridge.native.prior.weights()[7]
    for q in (-1, 1):
        effects, _ = measure.local_effects(
            board, route, rings, q, weights, include_flow=False
        )
        rows = rows.col_join(
            sp.Matrix(
                [
                    [c.inner(h, b) for h in pair for b in frame]
                    for pair in effects["reference"]
                ]
            )
        )
    # ||delta row|| <= sqrt(8) times the largest 4x4 operator error, by orthogonal
    # projection on the normalized 30 basis coordinates and two parity blocks.
    err = iv.sqrt(
        c.ri(
            8
            * (
                12 * mean_error**2
                + 2 * sp.Rational("2e-15") ** 2
                + 8 * sp.Rational("3e-15") ** 2
            )
        )
    )
    return rows, err


def response(rows, err):
    d, scale, adj, defect, size = solve(rows, err)
    data = json.loads((measure.HERE / "regional.json").read_text())
    out = {}
    from ...action_response.four_ticks import parse_interval

    work = c.ri(1) / scale
    for q in ("-1", "1"):
        eq = data["charges"][q]
        weights = parse_interval(eq["work_basis_coefficients"]["strain_work"])
        slope = work * weights
        # dM/dt=-X/scale on the occupied two-dimensional block, hence
        # ||d rho_R||_F^2=2/scale^2, exactly. Both blocks stay faithful.
        eig = parse_interval(eq["eigenvalues"][2]).a - c.ri(".01") * work
        assert eig > 0
        gain = c.ri(".01") * slope.a - c.ri(".01") ** 2 * work**2 / eig
        assert gain.a > 0
        out[q] = dict(
            entropy_derivative=str(slope),
            finite_entropy_increase_at_0p01_lower=str(gain.a),
        )
    return dict(
        constraint_rank=rows.rows,
        fixed_parity_dimension=30,
        true_source_correction_upper=str(adj.b),
        right_inverse_error_upper=str(defect.b),
        true_source_norm_upper=str(size.b),
        work_derivative_exact=f"1/{scale}",
        area_volume_ranges_and_marked_forces_constant=True,
        center_coordinates=[float(z) for z in d],
        direction_error_with_serialization_upper=float(
            np.nextafter(float(adj.b) + 2e-14, np.inf)
        ),
        native_lottery_slack_lower=str((c.ri(".20") - 16 * c.ri(".01") * size).a),
        charges=out,
        fair_untagged_load_region=fair_region(out, data["charges"]),
    )


def run():
    iv.dps = 32
    return response(*constraint_rows())


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    d = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "control.json").write_text(d)
    else:
        print(d, end="")
