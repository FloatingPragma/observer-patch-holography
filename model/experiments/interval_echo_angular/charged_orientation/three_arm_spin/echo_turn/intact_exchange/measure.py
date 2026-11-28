"""Compare routed charge and winding exchange on the three-arm material."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import pairwise, permutations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ... import apparatus, rotation_domain
from ...measure import BASE
from ...recoil_search.guides import detour

HERE = Path(__file__).resolve().parent


def blank_pulse():
    """The full native blank pulse is the unsigned load swap at this strength."""
    words = apparatus.LAW.SRC
    swap = np.eye(len(words), dtype=int)[:, [words.index((b, a)) for a, b in words]]
    cost = np.diag([(a - b) ** 2 for a, b in words])
    h8 = apparatus.blocks()["blank"]
    kappa = int(apparatus.candidate.KAPPA)
    np.testing.assert_array_equal(
        h8, 4 * (np.eye(len(words), dtype=int) - swap) + 8 * kappa * cost
    )
    np.testing.assert_array_equal(swap @ cost, cost @ swap)
    assert kappa % 2 == 0
    # The contrast factor is exp(-i pi kappa cost)=I. The two eigenspaces
    # of the swap have eigenvalues 0 and 1 under (I-S)/2.
    np.testing.assert_allclose(expm(-1j * np.pi * h8 / 8), swap, atol=2e-14, rtol=0)
    return dict(kappa=kappa, exposure_in_pi_units=1, closed_load_words=len(words))


def paired_routes(route, back, spare):
    front = list(pairwise(route))
    reverse = list(pairwise(back))
    assert len(front) == len(reverse)
    near, far = route[0], back[0]
    swap = [step for both in zip(front, reverse, strict=True) for step in both]
    hold = [(far, spare), *front, *reverse, (spare, far)]
    assert Counter(swap) + Counter([(far, spare), (spare, far)]) == Counter(hold)
    assert (near, far) == (back[-1], route[-1])
    return swap, hold


def routed_operator(steps, hubs, allowed, index):
    """Follow the two vacancy coordinates and an independent sign-transport rule."""
    positions = list(hubs)
    antisymmetric = 1
    for source, target in steps:
        assert frozenset((source, target)) in allowed
        assert source in positions and target not in positions
        moving = positions.index(source)
        other = positions[1 - moving]
        # A hop through the other vacancy in the chosen *ordering* picks
        # up the wedge sign. This ordering is computational, not a route.
        if (index[source] < index[other]) != (index[target] < index[other]):
            antisymmetric *= -1
        positions[moving] = target
    return positions, antisymmetric


def idle_compensator(ports, original):
    board = apparatus.winding.flight.federation(0)
    used = set(ports) | set(original)
    isolated = {
        p
        for p in range(board.n)
        if p not in used and not used.intersection(board.neighbours[p])
    }
    a, b = next((a, b) for a, b in board.seams if {a, b} <= isolated)
    shutters = [
        list(edge)
        for edge in board.seams
        if ({a, b} & set(edge)) and set(edge) != {a, b}
    ]
    return dict(excess_at=a, spare_at=b, shutters=shutters)


def observer_banks():
    """Verify analyzer sensitivity and retain source and receiver refusals exactly."""
    circular = BASE @ np.array([[1, 1], [1j, -1j]]) / np.sqrt(2)
    analyzer = np.kron(circular, circular)
    identity = np.eye(9)
    winding_swap = identity.reshape(3, 3, 3, 3).transpose(1, 0, 2, 3).reshape(9, 9)
    p_both_sources = Fraction(1, 9)
    p_all_apertures = Fraction(1, 81)
    source_refusal = 1 - p_both_sources
    detector_refusal = p_both_sources - p_all_apertures

    def phase_bank(active_sign, reference_sign):
        bright = Fraction((active_sign + reference_sign) ** 2, 4)
        dark = Fraction((active_sign - reference_sign) ** 2, 4)
        return [
            str(p)
            for p in (
                p_all_apertures * bright,
                p_all_apertures * dark,
                detector_refusal,
                source_refusal,
            )
        ]

    labels = ("++", "+-", "-+", "--")
    records = {}
    for label in ("+-", "-+"):
        column = labels.index(label)
        data = {}
        for query, op, expected in (
            ("native_blank_routing", identity, label),
            ("control_swaps_winding_with_centre", winding_swap, label[::-1]),
        ):
            probs = abs(analyzer.conj().T @ op @ analyzer[:, column]) ** 2
            np.testing.assert_allclose(
                probs, np.eye(4)[labels.index(expected)], atol=2e-15, rtol=0
            )
            data[query] = [
                str(p_all_apertures * int(flag == expected)) for flag in labels
            ]
            data[query] += [str(detector_refusal), str(source_refusal)]
        records[label] = data
    return dict(
        phase_bins=["bright", "dark", "receiver_refused", "source_refused"],
        native_phase=phase_bank(1, 1),
        antisymmetric_changed_rule=phase_bank(-1, 1),
        winding_bins=[*labels, "receiver_refused", "source_refused"],
        winding_records=records,
    )


def exchange_chart(length):
    f = apparatus.fixture()
    _, ports, legs, _, _, _, _ = detour(length - 19)
    hubs = (f["hub"], f["tip"])
    routes = [[hubs[0], *leg, hubs[1]] for leg in legs]
    allowed = {frozenset(edge) for path in routes for edge in pairwise(path)}
    assert all(len(path) == length + 1 for path in routes)
    assert not allowed.intersection(map(frozenset, f["marked"]))
    p = sorted(ports)
    orderings = (p, list(reversed(p)), p[::2] + p[1::2])
    orders = [{v: i for i, v in enumerate(word)} for word in orderings]
    physical_routes = []
    for a, b, c in permutations(range(3)):
        spare = routes[c][-2]
        swap, hold = paired_routes(routes[a], routes[b][::-1], spare)
        for index in orders:
            assert routed_operator(swap, hubs, allowed, index) == (list(hubs[::-1]), -1)
            assert routed_operator(hold, hubs, allowed, index) == (list(hubs), 1)
        physical_routes.append(
            dict(
                cycle_arms=[a, b],
                parking_port=spare,
                exchange_intervals=len(swap),
                parked_intervals=len(hold),
                native_amplitudes=[1, 1],
                antisymmetric_amplitudes=[-1, 1],
                head_operator="identity_on_64_words",
            )
        )
    return dict(
        length=length,
        arm_paths=routes,
        vacancy_hubs=list(hubs),
        first_compensator=f["donor"],
        second_compensator=idle_compensator(ports, f["donor"]),
        fixed_active_head_edges=[list(edge) for edge in f["marked"]],
        matched_pairs=physical_routes,
        computational_orderings=len(orders),
    )


def run():
    reading = rotation_domain.run()
    receiver = reading["receiver"]
    # The intact head mask and terminal cell are the same in all refinements.
    charts = [exchange_chart(length) for length in (19, 20, 21)]
    for chart in charts:
        assert (
            sorted(map(sorted, chart["fixed_active_head_edges"][3:]))
            == receiver["absolute_one_edges"]
        )
    forbidden = [
        row["rotation"] for row in receiver["proper_rotations"] if not row["same_live"]
    ]
    assert len(forbidden) == 59
    return dict(
        blank_pulse=blank_pulse(),
        intact_turns=dict(
            receiver=receiver["global_port"],
            forbidden_rotation_indices=forbidden,
            shadow_rotation_probability="0",
            native_echo_writes_magnitude_one=reading["native"][
                "blank_to_magnitude_one_echo_probability"
            ],
        ),
        charts=charts,
        banks=observer_banks(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output is not None:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    else:
        reference = json.loads((HERE / "results.json").read_text())
        if result != reference:
            raise ValueError("Replayed exact exchange receipt differs")
        print(
            "18 matched routes give a positive native phase. Both winding records remain at their launch observers."
        )


if __name__ == "__main__":
    main()
