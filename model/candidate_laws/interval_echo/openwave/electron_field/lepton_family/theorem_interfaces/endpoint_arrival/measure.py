"""Original-source admission and energy bounds for a charged spin at a fresh port."""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

from ..charged_endpoints.cube import CELLS, LINKS
from ..charged_endpoints.measure import BOUNDS, CONTRACT, SOURCE_REVISION, read_archive
from ..neutral_direction import verify as verify_endpoint_weights
from . import spin_read


REGIONAL = "model/particle_primary_regional_fermion/results.json"
CONJUGATION = "model/particle_primary_charge_conjugation/results.json"
INPUTS = {
    CONTRACT: "7617133bb473e37250f31491d5fadf3c4275d7b97a4e8342321d9e108988ed31",
    BOUNDS: "d2b43d798437761904974f1ca622f55ec22e4c1cb0769b2dd0ddcac8a1205a41",
    REGIONAL: "0f9d7bd1511f49252cc76617c1ba08b9e90cf3cce41802fab17439305f76a314",
    CONJUGATION: "98ce0fb71731466f462ffb25c93d3df811028a87999404229c71493403949f3e",
}


def inputs():
    data = {}
    for name, digest in INPUTS.items():
        raw = read_archive(name)
        if sha256(raw).hexdigest() != digest:
            raise ValueError(f"Different source for {name}")
        data[name] = json.loads(raw)
    member = data[CONTRACT]["fixed_primary_member"]
    source = data[BOUNDS]["selected_member_bound"]["result"]
    regional, conjugation = (
        data[name]["members"][0]["result"] for name in (REGIONAL, CONJUGATION)
    )
    for receipt in (source, regional, conjugation):
        assert all(receipt[k] == v for k, v in member["ids"].items())
        assert receipt["member_id"] == member["member_id"]
    assert conjugation["selected_polynomial_density_preserved"]
    assert conjugation["original_H_and_energy_zero_preserved"]
    assert regional["coherent_charge_meter_radix"] == 5
    return source, regional


def edge(x, y):
    return tuple(sorted((x, y)))


def route_counts():
    neighbors = {x: [y for y in CELLS if edge(x, y) in LINKS] for x in CELLS}
    faces = [
        quad
        for quad in combinations(CELLS, 4)
        if all(len(set(neighbors[x]) & set(quad)) == 2 for x in quad)
    ]
    assert len(faces) == 6
    paths = [(y, x, z) for x in CELLS for y, z in permutations(neighbors[x], 2)]
    uses = Counter()
    for i in range(3):
        assert set(Counter(path[i] for path in paths).values()) == {6}
    for y, x, z in paths:
        route = (edge(y, x), edge(x, z))
        uses.update(route)
        assert edge(y, z) not in LINKS
        assert sum(any(set(link) <= set(face) for link in route) for face in faces) == 3
        # Each other incident native Wilson column has squared norm 1/2.
        # Two active spin rows per site cancel that 1/2 in the squared sum.
        assert (
            sum(len(neighbors[p]) - sum(p in link for link in route) for p in (y, x, z))
            == 5
        )
        for sign in (-1, 1):
            flux, charge = Counter(), Counter()
            for initial, final in ((y, x), (x, z)):
                low, high = edge(initial, final)
                flux[low, high] += sign if final == low else -sign
                charge[initial] -= 6 * sign
                charge[final] += 6 * sign
                divergence = Counter()
                for (a, b), head in flux.items():
                    divergence[a] += 6 * head
                    divergence[b] -= 6 * head
                assert all(charge[u] == divergence[u] for u in CELLS)
            assert (charge[y], charge[x], charge[z]) == (-6 * sign, 0, 6 * sign)
    assert len(paths) == 48 and set(uses.values()) == {8} and set(uses) == set(LINKS)
    return len(paths)


def calculate():
    source, regional = inputs()
    weights = verify_endpoint_weights()
    assert (
        weights["pure_hypercharge_response"][0]
        == weights["electromagnetic_direction_response"][0]
        == "1"
    )
    paths = route_counts()
    spin = spin_read.check()
    allocation = source["allocation"]
    k = F(allocation["combined_bad_event_coefficient"])
    assert F(3) * F(allocation["retained_q6_onsite_fraction"]) == k / 16
    assert F(105, 2) * F(allocation["retained_electric_fraction"]) == k * F(5, 18)
    deficit = (
        F(source["selected_mean_energy_upper"])
        + 720
        + F(allocation["dimensionless_scalar_loss_upper"])
    ) / k
    assert deficit == F(source["combined_bad_event_mean_upper"])
    assert F(source["original_success_weight_lower"]) == 1 - deficit
    admission = 1 - 3 * deficit
    assert admission > F("0.38339964")

    center = 6 + 2 * F(63, 4)
    assert center == F(75, 2) and F(source["exact_reference_frequency"]) == 6 + F(63, 4)
    assert regional["hopping_control"]["selected_bond_squared_norm_upper"] == 40
    # Rational upper square roots suffice for an outward bound on all link fields.
    selected, other = F(63246, 10000), F(22361, 10000)
    assert selected**2 > 40 and other**2 > 5
    residual = 2 * selected + other + F(6, 7)
    assert residual**2 < 248
    radius = F(30)
    joint = admission - F(248) / radius**2 - F(1, 32)
    assert F(regional["complete_two_clock_event_loss_upper"]) <= F(1, 32)
    assert joint > F("0.07659409")
    resolutions = F(1, 8)
    band = (center - radius - resolutions, center + radius + resolutions)
    return {
        "member_id": source["member_id"],
        "source_revision": SOURCE_REVISION,
        "paths_per_sign": paths,
        "per_original_sign_setting": {
            "signed_spin_arrival_lower": float(admission),
            "signed_spin_and_transfer_band_lower": float(joint),
        },
        "clock_units": {
            "center": str(center),
            "reported_total_energy_transfer_band": list(map(str, band)),
            "centered_second_moment_upper": 248,
        },
        "transfer_at_holder_origin_receiver_in_sixths": [[-6, 0, 6], [6, 0, -6]],
        "correlations_on_each_nonzero_accepted_clock_bin": spin,
        "conjugate_instruments_have_equal_clock_laws_after_sign_reversal": True,
        "exact_lower_bounds": {
            "signed_spin": str(admission),
            "signed_spin_and_clock": str(joint),
        },
        "source_input": "original polynomial preparation with local addressed gates",
        "full_energy_histogram_enumerated": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    result = calculate()
    if options.output:
        with options.output.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        recorded = Path(__file__).with_name("results.json")
        assert json.loads(recorded.read_text(encoding="utf-8")) == result
        print(
            "Verified both signed spin arrivals and their common source-clock bounds."
        )


if __name__ == "__main__":
    main()
