"""Compare an equal-rest-energy axial sheet survey with a fixed Weyl source."""

import json
from math import pi
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo_partition_repair.published_strings import geometry
from .boosts import hopping
from .measure import DIPOLE, make_board, read_ball, regional_covariances


def continuum_stress(cutoff):
    """Covariant stress from the occupied massless ray moment of the measured spectrum."""
    nodes, weights = np.polynomial.legendre.leggauss(5)
    directions = np.concatenate((np.eye(3), -np.eye(3)))
    stress = np.zeros((4, 4))
    for node, weight in zip(nodes, weights, strict=True):
        momentum = cutoff * (1 + node) / 2
        base = 0.25 * (1 - (momentum / cutoff) ** 2) ** 2
        radial_measure = 2 * (cutoff * weight / 2) * momentum**3 / (2 * pi) ** 3
        for direction in directions:
            occupation = base * (1 + momentum * np.dot(DIPOLE, direction) / cutoff)
            covector = np.r_[-1.0, direction]
            stress += radial_measure * (4 * pi / 6) * occupation * np.outer(covector, covector)
    return stress


def signal_comparison():
    side, length, cutoff, floor, time = 4, 4.0, 2.5, 0.002, 0.72
    board = make_board(side, length, cutoff, floor)
    sites = np.indices((side,) * 3).reshape(3, -1).T
    reference, delta = regional_covariances(board, sites)
    source = reference + delta
    sender = 2 * np.ravel_multi_index((2, 2, 2), (side,) * 3)
    receiver = 2 * np.ravel_multi_index((2, 2, 3), (side,) * 3)
    sent = np.zeros(2 * side**3)
    sent[sender] = 1
    columns = np.column_stack((sent, (sent - reference[:, sender]), reference[:, sender],
                               (sent - source[:, sender]), source[:, sender]))
    evolved = expm_multiply(-1j * time * hopping(side, board.spacing), columns, traceA=0)
    target, vacuum_more, vacuum_less, source_more, source_less = evolved.T

    def complex_pair(z):
        return [float(z.real), float(z.imag)]

    return {
        "side": side, "time": time, "source_energy_density": board.energy,
        "region_entropy_change": read_ball(board, 1.2)["entropy_change"],
        "selected_receiver": {
            "vacuum_greater": complex_pair(vacuum_more[receiver]),
            "vacuum_lesser": complex_pair(vacuum_less[receiver]),
            "source_greater": complex_pair(source_more[receiver]),
            "source_lesser": complex_pair(source_less[receiver]),
            "free_signal": complex_pair(target[receiver]),
        },
        "maximum_source_change_in_lesser_signal": float(np.max(np.abs(vacuum_less - source_less))),
        "maximum_vacuum_sum_error": float(np.max(np.abs(vacuum_more + vacuum_less - target))),
        "maximum_source_sum_error": float(np.max(np.abs(source_more + source_less - target))),
    }


def audit():
    cutoff, rapidity = 2.5, 0.32
    stress = continuum_stress(cutoff)
    curvature = geometry.curvature(float(stress[0, 0]), 0)
    metric_tensor = np.asarray(geometry.einstein_tensor(curvature)) / geometry.KAPPA
    null = [np.array([1., 0., 0., s]) for s in (1., -1.)]
    timelike = [np.array([np.cosh(rapidity), 0., 0., s * np.sinh(rapidity)]) for s in (1., -1.)]

    def readings(tensor):
        return {
            "rest": float(tensor[0, 0]), "trace": float(np.sum(np.diag(tensor) * [-1, 1, 1, 1])),
            "axial_null": [float(k @ tensor @ k) for k in null],
            "mean_axial_boost_over_rest": float(np.mean([u @ tensor @ u for u in timelike]) / tensor[0, 0]),
        }

    return {"cutoff": cutoff, "rapidity": rapidity, "weyl_action": readings(stress),
            "sheet_einstein_over_coupling_at_equal_rest_energy": readings(metric_tensor),
            "finite_signal_on_two_preparations": signal_comparison()}


if __name__ == "__main__":
    Path(__file__).with_name("axial_response.json").write_text(json.dumps(audit(), indent=2) + "\n")
