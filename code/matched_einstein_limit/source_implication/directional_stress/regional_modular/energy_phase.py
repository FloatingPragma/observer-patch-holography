"""Independent phase-kick counts for the two fixed-energy-flow preparation controls."""

import argparse
import json

import numpy as np

from .. import bridge, phase_read
from . import measure


def run(proof=None):
    if proof is None:
        proof = json.loads((measure.HERE / "energy_odd.json").read_text())
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    pooled = None
    for q in (-1, 1):
        weights, bank = phase_read.sums_for_charge(board, route, rings, q)
        if pooled is None:
            pooled = {key: np.array(value) / 2 for key, value in bank.items()}
        else:
            for key, value in bank.items():
                pooled[key] += np.array(value) / 2
    ref, _, _ = measure.source()
    rho = np.array(
        [float(z) * np.array(v, complex) for z, v in zip((0.8657, 0.1343), ref)]
    )
    shots = 2**60
    seed = 2640138
    rng = np.random.default_rng(seed)
    magnitude = (60 + 24) * 64
    controls = ("equal_work_control", "fixed_work_control", "fixed_work_second_control")
    family = len(controls) * 2 * 19
    width = (
        magnitude
        * np.sqrt(2 * np.log(2 * family / 0.01) / shots)
        * np.abs(weights).sum()
    )
    groups = {}
    for key in controls:
        d = np.array(proof[key]["center_coordinates"]).reshape(2, 15)
        direction = np.einsum("bk,kij->bij", d, measure.c.iso.basis(False))
        slopes = {
            k: float(weights @ np.einsum("btij,bji->t", ef, direction).real)
            for k, ef in pooled.items()
        }
        assert abs(slopes["active"]) < 1e-8
        counts = []
        currents = []
        for s in (-0.01, 0.01):
            state = rho + s * direction
            assert min(np.linalg.eigvalsh(state).ravel()) > 0
            # Pool the untagged charge outcomes before the binary sampling.
            means = np.einsum("btij,bji->t", pooled["active"], state).real
            probability = (1 + means / magnitude) / 2
            assert np.all((0 < probability) & (probability < 1))
            yes = rng.binomial(shots, probability)
            counts.append(yes.tolist())
            currents.append(float(weights @ (magnitude * (2 * yes / shots - 1))))
        diff = currents[1] - currents[0]
        assert abs(diff - 0.02 * slopes["active"]) < 2 * width
        groups[key] = dict(
            pooled_phase_slope_at_numerical_control_center=slopes,
            pooled_current_estimates=currents,
            plus_minus_contrast_99pct=[
                float(diff - 2 * width),
                float(diff + 2 * width),
            ],
            plus_counts_at_minus_then_plus=counts,
        )
    return dict(
        trials_per_setting=shots,
        simultaneous_coverage=0.99,
        seed=seed,
        source_charge_tag_retained=False,
        controls=groups,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "energy_phase.json").write_text(output)
    else:
        print(output, end="")
