"""Can the three-arm charge send a local field through unoccupied guide material?"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..apparatus import LAW, fixture, native_graph
from ..measure import BASE, ONE
from ..echo_turn import pulse
from .measure import conjugation

HERE = Path(__file__).resolve().parent
STARTS = 2**37


def run():
    rng = np.random.default_rng(581027)
    chart = fixture()
    ports = chart["ports"]
    home, far = chart["hub"] // 12, chart["tip"] // 12
    owners = sorted({p // 12 for p in ports} - {home, far})
    for a, b in chart["marked"]:
        assert a // 12 == b // 12 and a // 12 in (home, far)
    assert not ({p // 12 for p in chart["donor"]} & set(owners))
    for a, b in chart["blank"]:
        assert a in ports and b in ports
    # Ready/spent orientation is physically discriminated by a terminus probe.
    # The full 99-state transfer is checked in both physical endpoint orders.
    for order in (False, True):
        pulse.primitive(int(order))
    assert len(LAW.SRC) == 9 and conjugation() == 99
    rows = []
    for sign in (-1, 1):
        generator = native_graph(sign)
        source = np.zeros((generator.shape[0], 2), complex)
        source[ONE] = BASE
        for time in (16, 32):
            operation = -1j * time * generator
            state = expm_multiply(operation, source, traceA=operation.diagonal().sum())
            weights = (abs(state.reshape(len(ports), 8, 8, 2)) ** 2).sum(
                axis=(2, 3)
            ) / 2
            np.testing.assert_allclose(weights.sum(), 1, atol=5e-11)
            total_by_port = weights.sum(axis=1)
            bulk = []
            for owner in owners:
                local = np.array([p // 12 == owner for p in ports])
                occupancy = float(total_by_port[local].sum())
                # Tracing the complement kills vacuum/one-defect cross terms.
                # If empty, each material tensor factor is in its unperturbed state.
                # The +/- occupied sectors are orthogonal, so unconditioned D equals P(occupied).
                draws = rng.multinomial(
                    STARTS, [occupancy / 3, (1 - occupancy) / 3, 2 / 3]
                )
                bulk.append(
                    dict(
                        owner=int(owner),
                        charge_probability_given_source=occupancy,
                        counts_charge_empty_refusal=draws.tolist(),
                    )
                )
            # Read the first physical winding head after excluding ANY occupant of the receiving owner.
            outside = np.array([p // 12 != far for p in ports])
            mask = (np.arange(8) & 1).astype(bool)
            changed = float(weights[outside][:, mask].sum())
            unchanged = float(weights[outside][:, ~mask].sum())
            receiver_charge = float(weights[~outside].sum())
            probs = np.array([changed / 3, unchanged / 3, receiver_charge / 3, 2 / 3])
            draws = rng.multinomial(STARTS, probs)
            rows.append(
                dict(
                    sign=sign,
                    time=time,
                    intermediate_owners=bulk,
                    receiver_edge=chart["marked"][3],
                    empty_receiver_changed_per_start=changed / 3,
                    counts_changed_unchanged_occupied_refused=draws.tolist(),
                )
            )
    return dict(
        starts_per_bank=STARTS,
        source_aperture=1 / 3,
        bulk_empty_conditional_trace_distance_from_idle=0,
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        (HERE / "empty_patch.json").write_text(json.dumps(result, indent=2) + "\n")
    for row in result["rows"]:
        print(
            row["sign"],
            row["time"],
            row["empty_receiver_changed_per_start"],
            row["counts_changed_unchanged_occupied_refused"][0],
        )
