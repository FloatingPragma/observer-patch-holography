"""Check which observers can read circular winding on the charged guide."""

import argparse
from itertools import pairwise
import json
from pathlib import Path

import numpy as np

from ... import apparatus
from ...measure import BASE, ONE
from ...recoil_search.guides import detour
from ...recoil_search.threshold import HERE as CLOCK_DIRECTORY, ci

HERE = Path(__file__).resolve().parent


def fixed_arm_signature(length):
    f = apparatus.fixture()
    _, ports, legs, h, actions, _, _ = detour(length - 19)
    circular = np.zeros((h.shape[0], 2), dtype=complex)
    circular[ONE] = BASE @ np.array([[1, 1], [1j, -1j]]) / np.sqrt(2)
    rows = []
    for arm, rail in enumerate(legs):
        i, j = (k for k in range(3) if k != arm)
        other_two_exchange = list(range(3))
        other_two_exchange[i], other_two_exchange[j] = j, i
        reflection = actions[tuple(other_two_exchange)]
        assert (h[reflection, :][:, reflection] - h).nnz == 0
        exchanged = circular[reflection]
        overlaps = circular.conj().T @ exchanged
        np.testing.assert_allclose(abs(overlaps), [[0, 1], [1, 0]], atol=2e-15, rtol=0)
        path = [f["hub"], *rail, f["tip"]]
        mark = [f["marked"][arm], f["marked"][arm + 3]]
        all_edges = [*pairwise(path), *mark]
        for a, b in all_edges:
            for v in (a, b):
                location = ports.index(v)
                # The other head labels can be relabelled together. All 64
                # components at either endpoint still have that same location.
                assert np.all(
                    reflection[location * 64 : (location + 1) * 64] // 64 == location
                )
        rows.append(dict(arm=arm, individually_fixed_edges=len(all_edges)))
    return rows


def run():
    archive = json.loads((CLOCK_DIRECTORY / "threshold.json").read_text())
    starts = archive["shots_per_witness_axis"]
    results = []
    for layout in archive["lengths"]:
        length = layout["length"]
        first = layout["lowest_and_next_harmonic"][0]
        same, different, refused = first["winding"][
            "counts_same_different_all_refusals"
        ][1]
        assert same + different + refused == starts
        independent_sites = fixed_arm_signature(length)
        accepted = same + different
        interval = ci(different, accepted).tolist()
        assert interval[0] > 0.9
        results.append(
            dict(
                physical_steps=length,
                clock_center=first["energy"],
                fixed_arm_signatures=independent_sites,
                exact_single_arm_probability_difference=0,
                joint_y_same_different_refused=[same, different, refused],
                conditional_different_fraction=different / accepted,
                conditional_different_interval=interval,
            )
        )
    return dict(
        readings="threshold.json Pauli Y at the first extended clock setting of each length",
        confidence_tail=archive["tail"],
        original_starts_per_setting=starts,
        charts=results,
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
        reference = json.loads((HERE / "orientation.json").read_text())
        if result != reference:
            raise ValueError("The orientation replay differs from its receipt")
        print(
            "Single-arm equality and multi-arm winding readback reproduce on all three physical lengths."
        )


if __name__ == "__main__":
    main()
