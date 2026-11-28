"""Connect the signed cut, angular apparatus and stationary current to the long route."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from ..field_charge_records.paths import enclosing_records, joint_support, local_terms
from ..fan_rotation.measure import physical_pose
from .fixture import FanPath, archived_fixture
from .spectral import verified_inertia

HERE = Path(__file__).resolve().parent


def run():
    cases = []
    with archived_fixture() as archive:
        for size in (3, 4):
            for sign in (-1, 1):
                path = FanPath(archive, size, sign, asymmetric=True)
                indices, load, heads, edges, counts = joint_support(
                    path, path.source(), sign, local_terms(path)
                )
                bits, cells = enclosing_records(path, sign, load, heads, edges)
                np.testing.assert_array_equal(indices, np.arange(len(path.basis)))
                np.testing.assert_array_equal(cells, path.owner)
                # On every receiving word the whole-cell enclosing-head change is the signed material unit.
                np.testing.assert_array_equal(bits[:, -1], path.owner == size - 1)
                ring = path.field.ring_laplacian(path.words_at[0])
                code = path.code
                np.testing.assert_array_equal(ring @ code, code)
                for bit in path.field.bits.values():
                    head = np.array([int(bool(w & bit)) for w in path.words_at[0]])
                    np.testing.assert_array_equal(
                        code.T @ (ring * (head[None, :] - head[:, None])) @ code, 0
                    )
                # Each material transfer crosses two distinct sharp ports. Restricted to a fixed hub its current is zero.
                board_case = {"left": [path.bridges[-1][0]], "right": path.cells[-1]}
                pose = physical_pose(board_case, path.field.edges, path.board)
                words = list(path.words_at[0])
                outer = np.array([words.index(w ^ path.field.masks[0]) for w in words])
                rz = np.diag([-1, 1])
                np.testing.assert_array_equal(code[outer], code @ rz)
                head_phase = np.array(
                    [
                        (-1)
                        ** (
                            bool(w & path.field.bits[1, 2])
                            + bool(w & path.field.bits[4, 5])
                        )
                        for w in words
                    ]
                )
                x = code.T @ (head_phase[:, None] * code)
                np.testing.assert_array_equal(head_phase[:, None] * code, code @ x)
                np.testing.assert_array_equal(x @ x, np.eye(2))
                np.testing.assert_array_equal(rz @ x, -x @ rz)
                cuts = [Fraction(15, 4), Fraction(9, 2), Fraction(145, 32)]
                inertia = [verified_inertia(path.h, cut)["count_below"] for cut in cuts]
                assert inertia == [1, 2, 2]
                assert np.isrealobj(path.h.data) and (path.h != path.h.T).nnz == 0
                cases.append(
                    {
                        "sign": sign,
                        "physical_hubs": [c[0] for c in path.cells],
                        "complete_reachable_records": len(indices),
                        "microscopic_directed_transfers": dict(counts),
                        "signed_enclosing_head_change_given_charge_in_last_cell": sign,
                        "proper_poses_preserving_receiver": pose[
                            "preserve_aperture_and_written_support"
                        ],
                        "proper_poses_moving_hub": pose["move_gateway_aperture"],
                        "other_hub_fixed_poses": len(
                            pose["preserve_gateway_but_change_written_support"]
                        ),
                        "counts_below_3_75_4_5_and_4_53125": inertia,
                        "local_receiver_bonds_tested": len(path.field.bits),
                        "captured_ring_current_compressions_exact": 0,
                        "native_first_outer_ring_action": rz.tolist(),
                        "supplied_two_head_action": x.astype(int).tolist(),
                        "inverse_order_fringe_sign": -1,
                        "repeated_native_ring_fringe_sign": 1,
                    }
                )
    return cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    cases = run()
    if args.record:
        with args.record.open("x", encoding="utf-8") as handle:
            json.dump(cases, handle, indent=2)
            handle.write("\n")
    else:
        assert cases == json.loads((HERE / "arrival_meters.json").read_text())
        print(
            "Reproduced the charge, current and angular checks on all four signed preparations."
        )
