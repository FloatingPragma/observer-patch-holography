"""Verified signs for the zero-or-one mobile with equilibrated, fair dipoles."""

import argparse
from fractions import Fraction as F
import json
from math import factorial, inf, nextafter
from pathlib import Path


HERE = Path(__file__).parent
BETA = F(1, 8)
PIN_BOUND = 17
SUPPORT_BOUND = 12
DIAGONAL_FLOOR = F(35, 2)
HOP_ROW_BOUND = 3
LOCAL_LOWER = F(3, 160000)
LOCAL_UPPER = BETA * SUPPORT_BOUND * PIN_BOUND
MOBILE_PORTS = 236


def uniform_range(distance):
    ratio = HOP_ROW_BOUND / DIAGONAL_FLOOR
    return (
        SUPPORT_BOUND * PIN_BOUND**2 / DIAGONAL_FLOOR
        * ratio ** (2 * distance) / (1 - ratio) ** 2
    )


def fixed_temperature_range(distance):
    hops = 2 * distance
    return (
        SUPPORT_BOUND * PIN_BOUND**2 * BETA
        * (HOP_ROW_BOUND * BETA) ** hops / factorial(hops)
    )


def exp_upper(x, last):
    """Positive Taylor series, with a geometric bound from the first omitted term."""
    if x < 0 or x >= last + 2:
        raise ValueError("the tail ratio must lie in [0,1)")
    partial = sum(x**j / factorial(j) for j in range(last + 1))
    return partial + x**(last + 1) / factorial(last + 1) / (1 - x / (last + 2))


def outward(low, high):
    return [nextafter(float(low), -inf), nextafter(float(high), inf)]


def run():
    # A no-hop visit to any affected vertex lower-bounds each single increment.
    assert exp_upper(F(133, 16), 40) < 5000
    archived = json.loads((HERE.parent / "unoriented.json").read_text())
    geometries = json.loads((HERE.parent / "matched_energy.json").read_text())
    distances = {
        row["distance"]: row["first_mixed"]["neighbourhood_distance"]
        for row in geometries["panels"] if row["mobile_charge"] == 1
        and not row["reversed_receiver"]
    }
    rows = []
    for panel in archived["panels"]:
        distance = distances[panel["distance"]]
        # Interpret each outward binary64 endpoint at its exact stored value.
        low, high = map(F, panel["mean_heat_trace"])
        assert 0 < low <= high < BETA * fixed_temperature_range(distance)
        denominator_floor = LOCAL_LOWER**2 - MOBILE_PORTS * high
        assert denominator_floor > 0
        critical_low = low / LOCAL_UPPER**2
        critical_high = high / denominator_floor
        signs = []
        for power in (72, 35):
            activity = F(1, 10**power)
            lower = low - activity * LOCAL_UPPER**2
            upper = high * (1 + activity * MOBILE_PORTS) - activity * LOCAL_LOWER**2
            sign = "unresolved"
            if lower > 0:
                sign = "negative"
            elif upper < 0:
                sign = "positive"
            signs.append(dict(activity=f"1e-{power}", sign_numerator=outward(lower, upper),
                              connected_grand_potential=sign))
        rows.append(dict(cell_distance=panel["distance"], support_distance=distance,
                         measured_heat_interval=panel["mean_heat_trace"],
                         fixed_beta_potential_slope_bound=outward(0, fixed_temperature_range(distance))[1],
                         uniform_beta_potential_slope_bound=outward(0, uniform_range(distance))[1],
                         zero_crossing_activity=outward(critical_low, critical_high),
                         signs=signs))
    return dict(ensemble="zero or one mobile; fair charge and orientations before logarithm",
                beta="1/8", local_increment_bounds=[str(LOCAL_LOWER), str(LOCAL_UPPER)],
                active_contrast=8, panels=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    content = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(content)
    else:
        print(content, end="")
