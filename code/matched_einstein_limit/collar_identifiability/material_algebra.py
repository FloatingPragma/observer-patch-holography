"""Modular responses after restriction to complete material outcome registers."""

from fractions import Fraction
import json

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_pairing import measure as paired
from ..seam_overlap import protected


def restricted_slope(perturbation, reference):
    """Reference and perturbation are already probability vectors on the same outcome register."""
    if np.any(reference <= 0) or abs(perturbation.sum()) > 1e-11:
        raise ValueError("the finite log requires positive reference and zero tangent trace")
    return float(-perturbation @ np.log(reference))


def run():
    # Each integer scalar column acts as an overall phase for every endpoint input.
    hidden = protected.audit()["protected_logical_registers"]
    assert hidden["max_integer_scalar_action_residual"] == 0
    assert Fraction(hidden["energy_contrast"]) == Fraction(1, 4)

    board = paired.formation.federation(1)
    route = paired.formation.region(board)
    joint = paired.prepared_states(board, route, 1)[0]
    joint = joint / np.trace(joint).real
    independent = np.kron(*paired.formation.marginals(joint))
    keys = ("relay_miss", "return_miss", "complete")
    output = []
    for strength in (8, 0):
        # Separate source preparation from the read pulse, so the same two arms feed both laws.
        h = paired.formation.hamiltonian(board, route, 1, strength)
        row, = paired.material_read(h, route[1:3], {1: (joint, independent)})
        first = row["formed"]["appointments"]
        second = row["independent"]["appointments"]
        for a, b in zip(first, second, strict=True):
            p = np.array([a[k] for k in keys])
            q = np.array([b[k] for k in keys])
            np.testing.assert_allclose([p.sum(), q.sum()], 1., atol=1e-13)
            ordinary = restricted_slope(p-q, p)
            subtracted = ordinary - restricted_slope(p-q, q)
            # This independent expression is the negative symmetric relative entropy.
            np.testing.assert_allclose(subtracted, -(p-q) @ np.log(p/q), atol=1e-15)
            output.append(dict(read_contrast=strength, exposure=a["time"], formed=p.tolist(),
                               independent=q.tolist(), modular_slope=ordinary,
                               independent_reference_subtracted=subtracted))
    return dict(hidden_material_slope=0, hidden_hamiltonian_energy=hidden["energy_contrast"],
                material_outcomes=keys, visible_returns=output)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
