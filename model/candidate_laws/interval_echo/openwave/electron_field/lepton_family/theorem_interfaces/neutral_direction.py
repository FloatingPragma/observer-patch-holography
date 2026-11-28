"""Which independent neutral-current comparison resolves the endpoint's charge direction?"""

import json
from fractions import Fraction

import sympy as sp

from .charged_endpoints.cube import MASKS, Q6
from .charged_endpoints.measure import CHARGE_READS, ROOT, read_archive


def verify():
    # The first row is the actual endpoint representation on the eight-cell probe.
    # The other two rows are proposed preparations on the same calibrated weak frame.
    dictionary = json.loads((ROOT / "code/defect_particle_classes/results.json").read_bytes())
    table = {row["mask"]: row for row in dictionary["patterns"]}
    singlet = MASKS[Q6]
    if singlet != 24 or (table[singlet]["n3"], table[singlet]["n2"]) != (0, 2):
        raise ValueError("Different charged endpoint representation")
    bank = json.loads(read_archive(CHARGE_READS))["selected_member_bound"]["result"]
    charge_menu = bank["charge_menu"]["admitted_local_transfer_alphabet"]
    if charge_menu != [-6, 0, 6] or [m for m in MASKS if table[m]["six_Y"] == 6] != [singlet]:
        raise ValueError("The archived charge-six singlet is not unique on the declared matter carrier")
    successes = [x["charge_shift"] for x in bank["source_probe_charge_records"] if not x.get("failure")]
    if len(successes) != 48 or any(sorted(x) != [-6, 0, 0, 0, 0, 0, 0, 6] for x in successes):
        raise ValueError("Different native endpoint transfer alphabet")
    masks = [singlet, 15, 23]
    weights = sp.Matrix([[sp.Rational(table[m][key]) for key in ("T3", "Y")] for m in masks])
    gauge_directions = weights[:2, :]
    # Reconstruct the generator from two signed readings, keeping the third for a prediction.
    s, u = sp.symbols("singlet neutral_partner", real=True)
    a, b = gauge_directions.inv() @ sp.Matrix([s, u])
    assert (a, b) == (s+2*u, s)
    assert (weights[2, :] @ sp.Matrix([a, b]))[0] == -s-u
    assert weights[0, :].rank() == 1 and gauge_directions.rank() == 2
    hypercharge = weights @ sp.Matrix([0, 1])
    electromagnetic = weights @ sp.Matrix([1, 1])
    assert list(hypercharge) == [1, -Fraction(1, 2), -Fraction(1, 2)]
    assert list(electromagnetic) == [1, 0, -1]
    return {"mask_order": masks, "endpoint_probes_bound_to_original_transfer_alphabet": len(successes),
            "pure_hypercharge_response": [str(x) for x in hypercharge],
            "electromagnetic_direction_response": [str(x) for x in electromagnetic],
            "two_read_reconstruction": ["a = s + 2u", "b = s"],
            "held_charged_doublet_response": "-s-u",
            "same_physical_neutral_direction_probe_executed": False}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
