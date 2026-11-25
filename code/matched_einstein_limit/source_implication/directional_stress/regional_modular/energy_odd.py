"""Odd sharp-screen focusing at exactly fixed independently metered bridge energy flow."""

import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response.four_ticks import parse_interval
from .. import bridge
from . import control, measure, optical_certificate, screen_ward

c = measure.c


def energy_rows():
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    launches = bridge.beam.indices(sum(rings, []))
    ids = bridge.beam.indices(list(route[1:3]) + sum(rings, []))
    weights = bridge.native.prior.weights()[7]
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    ref = measure.source()
    effects = []
    reference = []
    for charge in (-1, 1):
        pair = []
        for odd in (False, True):
            h = bridge.native.generator(board, route, charge, odd)
            cells = bridge.operators(board, route, charge, odd)
            bridge.audit(board, cells, route, h)
            op = bridge.current(board, cells, route[1:3])
            cut = op[ids, :][:, ids]
            assert cut.count_nonzero() == op.count_nonzero()
            assert abs(cut).sum(axis=1).max() <= 8 * 500 and len(ids) == 12 * 4
            w = measure.four_tick_columns(h, ids, launches, weights)
            # Twice-action numerators; physical inward J_E = i [b,diff]/8.
            pair.append(measure.compress(w, cut, -8 * sp.I))
        effects.append([c.inner(v, a) for v in pair for a in frame])
        reference.append(measure.mean(pair, ref, sp.Rational("1e-15")))
    # For each launch 12 port amplitudes of radius 1e-19 give ||delta U||<4e-19.
    # ||J_E||<=500, hence 500*(2+4e-19)*4e-19 < 1e-15 per compressed effect.
    return sp.Matrix(effects), reference


def screen_reading(arm):
    sample = screen_ward.run(arm)
    data = optical_certificate.run(arm, sample, compare_charges=False)
    rays = data["fair_untagged_images"]
    odd = {
        key: parse_interval(rays["forward"][key]) - parse_interval(rays["reverse"][key])
        for key in rays["forward"]
    }
    assert 0 not in odd["source_focusing_slope"]
    assert 0 not in odd["source_expansion_slope"]
    return dict(
        fair_untagged_images=rays,
        odd_forward_minus_reverse={k: str(v) for k, v in odd.items()},
        required_common_clock_acceleration=str(
            -odd["source_focusing_slope"] / odd["source_expansion_slope"]
        ),
        independent_arrivals={
            k: v["extrapolated_arrivals"]
            for k, v in sample["fair_untagged_images"].items()
        },
    )


def fixed_control(rows, err, component):
    # Fix every previous mean, including both strain works. The extra control
    # drives only the first head coordinate (independently selected even-parity Pauli); no optical
    # statistic enters either preparation. This also tests a common clock.
    ix = sp.zeros(1, 30)
    ix[component] = 1
    target = sp.zeros(25, 1)
    target[-1] = 1
    center, scale, adjustment, defect, norm = control.solve(
        rows.col_join(ix), err, target
    )
    shear = dict(
        constraint_rank=25,
        even_pauli=("IX", "IY")[component],
        even_pauli_derivative_exact=f"1/{scale}",
        work_derivative_exact="0",
        right_inverse_error_upper=str(defect.b),
        true_source_norm_upper=str(norm.b),
        native_lottery_slack_lower=str((c.ri(".20") - 16 * c.ri(".01") * norm).a),
        center_coordinates=[float(x) for x in center],
        direction_error_with_serialization_upper=float(
            np.nextafter(float(adjustment.b) + 2e-14, np.inf)
        ),
    )
    shear["screen"] = screen_reading(shear)
    return shear


def run():
    iv.dps = 32
    rows, err = control.constraint_rows()
    energy, reference = energy_rows()
    rows = rows.col_join(energy)
    err = iv.sqrt(err**2 + c.ri(8 * 2 * sp.Rational("1e-15") ** 2))
    work = control.response(rows, err)
    work["screen"] = screen_reading(work)

    shear = fixed_control(rows, err, 0)
    other = fixed_control(rows, err, 1)
    first, second = [
        parse_interval(x["screen"]["required_common_clock_acceleration"])
        for x in (work, shear)
    ]
    assert first.b < second.a or second.b < first.a
    third = parse_interval(other["screen"]["required_common_clock_acceleration"])
    assert third.b < second.a or second.b < third.a
    return dict(
        central_bridge=[1, 38],
        flight_ticks=4,
        fixed_energy_current_derivatives_by_charge=[0, 0],
        same_symmetric_charge_and_ten_launch_ensemble=True,
        reference_energy_current_by_charge=list(map(str, reference)),
        fair_reference_energy_current=str(sum(reference) / 2),
        equal_work_control=work,
        fixed_work_control=shear,
        fixed_work_second_control=other,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "energy_odd.json").write_text(out)
    else:
        print(out, end="")
