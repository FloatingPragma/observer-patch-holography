"""Choose a material control at fixed four-station number weights and bridge meters."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv
from scipy import sparse

from . import control, energy_odd, measure

c = measure.c


def run():
    iv.dps = 34
    b = measure.bridge
    board = b.native.iso.pairing.formation.federation(1)
    route = b.native.iso.pairing.formation.region(board)
    rings = b.beam.screens(board, route)
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    rows, err = control.constraint_rows()
    energy, _ = energy_odd.energy_rows()
    rows = rows.col_join(energy)
    occ = []
    for q in (-1, 1):
        row = []
        for odd in (False, True):
            h = b.native.generator(board, route, q, odd)
            w = measure.four_tick_columns(
                h,
                b.beam.indices(route),
                b.beam.indices(sum(rings, [])),
                b.native.prior.weights()[7],
            )
            e = measure.compress(w, sparse.eye(16, dtype=int, format="csr"), 1)
            row.extend(c.inner(e, a) for a in frame)
        occ.append(row)
    rows = rows.col_join(sp.Matrix(occ))
    total_err = iv.sqrt(
        err**2 + c.ri(16 * (sp.Rational("1e-15") ** 2 + sp.Rational("1e-18") ** 2))
    )
    ix = sp.zeros(1, 30)
    ix[0] = 1
    target = sp.zeros(27, 1)
    target[-1] = 1
    d, scale, adjust, defect, norm = control.solve(rows.col_join(ix), total_err, target)
    return dict(
        ports=route,
        region_occupations_derivatives=[0, 0],
        charge_unresolved_number_sector_probabilities_constant=True,
        fixed_readout_count=26,
        augmented_rank=27,
        material_drive="even IX/2",
        drive_derivative_exact=f"1/{scale}",
        right_inverse_defect_upper=str(defect.b),
        source_direction_norm_upper=str(norm.b),
        source_segment=[-0.01, 0.01],
        lottery_slack_lower=str((c.ri(".2") - 16 * c.ri(".01") * norm).a),
        center_coordinates=[float(v) for v in d],
        direction_error_with_serialization_upper=float(
            np.nextafter(float(adjust.b) + 2e-14, np.inf)
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    txt = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "station_control.json").write_text(txt)
    else:
        print(txt, end="")
