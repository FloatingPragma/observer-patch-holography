"""Minimum extra response dimension for the measured regional entropy and two fluxes."""

import argparse
import json
from pathlib import Path

import sympy as sp

from . import late_read
from ..source import certificate as cert

HERE = Path(__file__).parent


def run():
    fronts = late_read.front_response_matrix()
    _, pump_error = late_read.exact_pump()
    center = late_read.rounded_matrix(fronts)
    widths = (fronts - center).applyfunc(abs) + sp.ones(6, 9) * pump_error
    region = json.loads((HERE / "regional.json").read_text())["regions"][0]
    cuts = json.loads((HERE / "current.json").read_text())["regions"]
    for row in [
        region["regional_entropy_slopes"],
        *(x["derivatives"]["energy"] for x in cuts),
    ]:
        intervals = [[sp.Rational(x) for x in row[k]] for k in late_read.AXES]
        center = center.col_join(sp.Matrix([[(a + b) / 2 for a, b in intervals]]))
        widths = widths.col_join(sp.Matrix([[(b - a) / 2 for a, b in intervals]]))
    # Each true matrix is center+E with |E| bounded entrywise by widths.
    # This infinity-norm Neumann test keeps the heterogeneous row errors.
    defect = center.inv().applyfunc(abs) * widths
    bound = max(sum(defect.row(i)) for i in range(9))
    assert bound < sp.Rational(12, 1000)
    return dict(
        front_rank=6,
        rank_with_regional_entropy_and_both_energy_currents=9,
        regional_radius=region["radius"],
        energy_cut_radii=[x["radius"] for x in cuts],
        source_axes=list(late_read.AXES),
        independent_additional_response_directions_required=3,
        inverse_defect_upper=cert.enclosed(bound, 0)[1],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "source_rank.json").write_text(data)
    else:
        print(data, end="")
