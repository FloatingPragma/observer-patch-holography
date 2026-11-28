"""Compare the maximum simultaneous native sharp energy with independent mass budgets."""

import argparse
from fractions import Fraction as F
import json
from math import ceil
from cosmology.rar_deep_regime.compensation import measure as sky
from . import patch, source
from ..closed_region.measure import check


def layers(ratios, stock):
    return [ceil(x / float(stock)) for x in ratios]


def run():
    standard = source.run_case(8)[1]
    # The largest b/N is three. A single occupied owner attains it.
    before = 72 + 3 * standard["ordinary_incident"]
    ceiling = 3 * standard["sharp"]["gain"] / before
    assert ceiling == F(17, 260)
    compact = patch.run()[-1]["patches"]["one face"]["gain_relative_to_initial"]
    fields = json.loads((sky.HERE / "inputs.json").read_text())["mistele_all"]["outer"]
    case = sky.mass_budget(fields)["cases"]["two_systematic_double_baryons"]
    ratios = case["min_exterior_debit_over_allowed_galaxy_baryons"]
    cosmic = sky.OMEGA_C_H2 / sky.OMEGA_B_H2
    return dict(
        radius_kpc=fields["radius_kpc"],
        favorable_lensing_ratios=ratios,
        perfect_delivery_max_sharp_mass_fraction=str(ceiling),
        compact_64_max_sharp_mass_fraction=str(compact),
        minimum_maximal_layers=layers(ratios, ceiling),
        minimum_compact_64_layers=layers(ratios, compact),
        planck_base_lcdm_cold_to_baryon=cosmic,
        planck_minimum_maximal_layers=layers([cosmic], ceiling)[0],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = run()
    path = source.HERE / "budget.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), result)
