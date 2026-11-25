"""Held completion slopes on the existing signed four-station survey."""

import numpy as np

from ...source_implication import geometry_gain as survey


def held_bridge_menu():
    board = survey.pairing.formation.federation(1)
    route = survey.pairing.formation.region(board)
    joint, independent, _ = survey.pairing.prepared_states(board, route, 1)
    pair = (joint / np.trace(joint).real, independent / np.trace(independent).real)
    shapes = survey.run()["panels"]
    area_ratios = [
        row["pairing_area_derivative"] / row["ordinary_area_derivative"]
        for row in shapes
    ]
    assert all(row["ordinary_area_derivative"] > 0 for row in shapes)
    rows = []
    for charge in (-1, 1):
        generators = [
            survey.ranging.generator(board, route, charge, 8, odd)
            for odd in (False, True)
        ]
        for reverse in (False, True):
            ends = route[1:3][::-1] if reverse else route[1:3]
            even, odd = [
                survey.pairing.material_read(h, ends, {1: pair})[0] for h in generators
            ]
            for index, t in enumerate(survey.pairing.TIMES):

                def value(state, kind):
                    return state[kind]["appointments"][index]["complete"]

                dark = value(even, "formed") - value(even, "independent")
                ordinary = value(odd, "independent") - value(even, "independent")
                assert dark < 0 and ordinary < 0
                rows.append(
                    dict(
                        probe=charge,
                        ports=ends,
                        exposure_per_leg=t,
                        pairing_completion_slope=dark,
                        ordinary_completion_slope=ordinary,
                        required_area_ratio=dark / ordinary,
                        transfer_fraction=[a * ordinary / dark for a in area_ratios],
                    )
                )
    return dict(
        area_slopes=[
            {
                k: p[k]
                for k in (
                    "tick",
                    "horizon",
                    "pairing_area_derivative",
                    "ordinary_area_derivative",
                )
            }
            for p in shapes
        ],
        area_ratio_range=[min(area_ratios), max(area_ratios)],
        completion=rows,
        maximum_convex_transfer_fraction=max(area_ratios)
        / min(row["required_area_ratio"] for row in rows),
    )
