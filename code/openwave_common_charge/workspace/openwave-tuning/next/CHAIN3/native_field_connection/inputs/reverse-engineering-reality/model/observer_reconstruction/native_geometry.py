"""Recover the native response geometry from committed current records.

The executor uses the existing twelve-port primitive and supplied calibration
schedule. The estimator receives only preparation labels, admitted seam labels,
committed coins/currents and their exact probabilities. It never sees output
loads, source coordinates or an incidence matrix. This is exact ensemble
identification on a declared finite apparatus, not statistical reconstruction
from a finite sample or identification of physical space with the scalar grid.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

import numpy as np


PORTS = 12
RECORD_FIELDS = {"preparation_port", "seam", "coin", "current", "probability"}


def _port(value):
    if type(value) is not int or not 0 <= value < PORTS:
        raise ValueError("a native port is an integer from zero through eleven")
    return value


def _probability(value):
    if type(value) is not str:
        raise ValueError("record probabilities must be exact rational strings")
    try:
        result = F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("record probabilities must be exact rational strings") from error
    if not 0 < result <= 1:
        raise ValueError("positive normalized record probabilities required")
    return result


def execute_calibration(seam_weights=None, permutation=None):
    """Execute all branches, returning observations without working-state data.

    Each preparation contains two loads at its labelled port. The fair native
    coin is retained even where it does not affect the current. Optional exact
    weights change the supplied scheduler; a permutation relabels the entire
    observable interface, including preparation and oriented seams.
    """
    from ..common_source import slow_geometry
    from ..observer import native

    seams = tuple(tuple(seam) for seam in slow_geometry.source().seams())
    if permutation is None:
        permutation = tuple(range(PORTS))
    if (not isinstance(permutation, (tuple, list)) or len(permutation) != PORTS
            or any(type(port) is not int for port in permutation)
            or set(permutation) != set(range(PORTS))):
        raise ValueError("a permutation of all twelve port labels is required")
    if seam_weights is None:
        weights = {seam: F(1, len(seams)) for seam in seams}
    else:
        if (not isinstance(seam_weights, dict) or not seam_weights
                or any(seam not in seams for seam in seam_weights)
                or any(type(weight) not in (int, F) or weight < 0
                       for weight in seam_weights.values())):
            raise ValueError("exact nonnegative weights on admitted source seams required")
        weights = {seam: F(weight) for seam, weight in seam_weights.items() if weight}
        if sum(weights.values(), F(0)) != 1:
            raise ValueError("the seam scheduler must have total probability one")
    records = []
    for preparation in range(PORTS):
        loads = tuple(2 * int(port == preparation) for port in range(PORTS))
        for seam, weight in sorted(weights.items()):
            for coin in (0, 1):
                _, record = native.native(loads, seam, coin)
                records.append({
                    "preparation_port": permutation[preparation],
                    "seam": [permutation[seam[0]], permutation[seam[1]]],
                    "coin": record[0], "current": record[1],
                    "probability": str(weight / 2),
                })
    return records


def _qadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def _qmul(left, right):
    return (left[0]*right[0] + 5*left[1]*right[1],
            left[0]*right[1] + left[1]*right[0])


def _qscale(value, scalar):
    return value[0]*scalar, value[1]*scalar


def _qmatrix_product(left, right):
    return [[sum((left[i][k][0]*right[k][j][0]
                  + 5*left[i][k][1]*right[k][j][1] for k in range(PORTS)), F(0)),
             sum((left[i][k][0]*right[k][j][1]
                  + left[i][k][1]*right[k][j][0] for k in range(PORTS)), F(0))]
            for i in range(PORTS) for j in range(PORTS)]


def _qproduct(left, right):
    flat = _qmatrix_product(left, right)
    return [[tuple(flat[PORTS*i+j]) for j in range(PORTS)] for i in range(PORTS)]


def _shift(matrix, scalar):
    return [[_qadd(matrix[i][j], _qscale(scalar, -int(i == j)))
             for j in range(PORTS)] for i in range(PORTS)]


def _encode_q5(matrix):
    return [[[str(a), str(b)] for a, b in row] for row in matrix]


def _band_certificate(laplacian):
    """Verify a candidate spectral band on the recovered operator exactly.

    The candidate roots 0, 5-sqrt(5), 6, 5+sqrt(5) are checked against the
    recovered matrix; a different scheduler cannot inherit this certificate.
    No incidence table or source-frame coordinates enter this calculation.
    """
    zero = (F(0), F(0))
    matrix = [[(value, F(0)) for value in row] for row in laplacian]
    numerator = _qproduct(_qproduct(matrix, _shift(matrix, (F(6), F(0)))),
                          _shift(matrix, (F(5), F(1))))
    denominator = _qmul((F(5), F(-1)), _qmul((F(-1), F(-1)), (F(0), F(-2))))
    norm = denominator[0]**2 - 5*denominator[1]**2
    inverse = denominator[0]/norm, -denominator[1]/norm
    projector = [[_qmul(value, inverse) for value in row] for row in numerator]
    square = _qproduct(projector, projector)
    image = _qproduct(matrix, projector)
    eigen = [[_qmul((F(5), F(-1)), value) for value in row] for row in projector]
    # This polynomial excludes an unexamined spectrum hidden off the band.
    polynomial = _qproduct(_qproduct(matrix, _shift(matrix, (F(6), F(0)))),
                           _qproduct(_shift(matrix, (F(5), F(-1))),
                                     _shift(matrix, (F(5), F(1)))))
    checks = {
        "full_spectral_polynomial_zero_exact": all(value == zero for row in polynomial for value in row),
        "projector_idempotent_exact": square == projector,
        "projector_symmetric_exact": all(projector[i][j] == projector[j][i]
                                         for i in range(PORTS) for j in range(PORTS)),
        "slow_eigenvalue_identity_exact": image == eigen,
        "projector_trace_three_exact": _qsum(projector[i][i] for i in range(PORTS)) == (F(3), F(0)),
        "gram_unit_diagonal_exact": all(value == (F(1, 4), F(0))
                                        for value in (projector[i][i] for i in range(PORTS))),
        "constant_mode_removed_exact": all(_qsum(row) == zero for row in projector),
    }
    return checks, projector


def _qsum(values):
    total = F(0), F(0)
    for value in values:
        total = _qadd(total, value)
    return total


def recover(records):
    """Estimate T and its geometry solely from complete observation laws.

    The readout dictionary supplies delta x = current*(e_b-e_a), and the
    preparation supplies x=2e_j. Hence the conditional mean output divided
    by two is column j of T. L=60*(I-T) is a declared normalization, whose
    unweighted-graph and rank-three consequences are independently tested.
    """
    if not isinstance(records, (tuple, list)) or not records:
        raise ValueError("a nonempty complete record table is required")
    transition = [[F(int(i == j)) for j in range(PORTS)] for i in range(PORTS)]
    masses = [F(0)] * PORTS
    schedules = [{} for _ in range(PORTS)]
    seen = set()
    for row in records:
        if type(row) is not dict or set(row) != RECORD_FIELDS:
            raise ValueError("only the declared observer record fields are accepted")
        preparation = _port(row["preparation_port"])
        seam = row["seam"]
        if not isinstance(seam, (tuple, list)) or len(seam) != 2:
            raise ValueError("an oriented two-port seam is required")
        a, b = (_port(port) for port in seam)
        if a == b:
            raise ValueError("a seam needs two distinct ports")
        coin, current = row["coin"], row["current"]
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError("a committed native coin is zero or one")
        if type(current) is not int or current not in (-1, 0, 1):
            raise ValueError("a committed total-two current is -1, zero or one")
        probability = _probability(row["probability"])
        key = preparation, a, b, coin
        if key in seen:
            raise ValueError("duplicate preparation/seam/coin record")
        seen.add(key)
        reconstructed = [2*int(i == preparation) for i in range(PORTS)]
        reconstructed[a] -= current
        reconstructed[b] += current
        if min(reconstructed) < 0:
            raise ValueError("the record contradicts the prepared conserved load")
        masses[preparation] += probability
        schedules[preparation][a, b, coin] = probability
        transition[a][preparation] -= probability * current / 2
        transition[b][preparation] += probability * current / 2
    if any(mass != 1 for mass in masses):
        raise ValueError("all twelve preparations require complete normalized outcomes")
    if any(schedule != schedules[0] for schedule in schedules[1:]):
        raise ValueError("the same calibration scheduler is required for every preparation")
    for (a, b, coin), probability in schedules[0].items():
        if schedules[0].get((a, b, 1-coin)) != probability:
            raise ValueError("both fair-coin outcomes must be retained")
        if (b, a, coin) in schedules[0]:
            raise ValueError("one orientation per seam is required")
    laplacian = [[60*(F(int(i == j))-transition[i][j]) for j in range(PORTS)]
                 for i in range(PORTS)]
    symmetric = all(laplacian[i][j] == laplacian[j][i]
                    for i in range(PORTS) for j in range(PORTS))
    graph = (symmetric and all(laplacian[i][i] == 5 for i in range(PORTS))
             and all(laplacian[i][j] in (0, -1) for i in range(PORTS)
                     for j in range(PORTS) if i != j))
    checks, projector = _band_certificate(laplacian)
    recovered = symmetric and graph and all(checks.values())
    gram = [[_qscale(value, 4) for value in row] for row in projector] if recovered else None
    numeric = ([[float(a) + math.sqrt(5)*float(b) for a, b in row] for row in gram]
               if recovered else None)
    return {
        "record_count": len(records), "preparation_count": PORTS,
        "transition_matrix_exact": [[str(value) for value in row] for row in transition],
        "laplacian_normalization": "L=60*(I-T); T acts on column load vectors",
        "laplacian_exact": [[str(value) for value in row] for row in laplacian],
        "symmetric_response_exact": symmetric,
        "degree_five_unweighted_graph_recovered_exact": graph,
        "rank_three_unit_gram_recovered": recovered,
        "band_certificate": checks,
        "slow_projector_qsqrt5": _encode_q5(projector) if recovered else None,
        "gram_qsqrt5": _encode_q5(gram) if recovered else None,
        "gram_numeric": numeric,
        "squared_distances_numeric": [[numeric[i][i]+numeric[j][j]-2*numeric[i][j]
                                        for j in range(PORTS)] for i in range(PORTS)] if recovered else None,
        "exact_number_encoding": "each pair [a,b] denotes a+b*sqrt(5)",
        "geometry_identified_up_to": "orthogonal isometry; port labels transport under interface relabeling",
        "hidden_output_loads_used": False,
        "source_incidence_or_coordinates_used_by_estimator": False,
    }


def report():
    """Retain the exact recovery and changed-scheduler/interface controls."""
    records = execute_calibration()
    result = recover(records)
    if not result["rank_three_unit_gram_recovered"]:
        raise AssertionError("record-only native response recovery failed")
    permutation = (7, 2, 10, 0, 5, 11, 8, 4, 1, 9, 3, 6)
    relabelled = recover(execute_calibration(permutation=permutation))
    relabel_error = max(abs(result["gram_numeric"][i][j]
                           - relabelled["gram_numeric"][permutation[i]][permutation[j]])
                       for i in range(PORTS) for j in range(PORTS))
    seams = sorted({tuple(row["seam"]) for row in records})
    biased_weights = {seam: F(1, 31) for seam in seams}
    biased_weights[seams[0]] = F(2, 31)
    biased = recover(execute_calibration(biased_weights))
    fixed = recover(execute_calibration({seams[0]: F(1)}))
    erased = [{**row, "current": 0} for row in records]
    tampered = recover(erased)
    # Coordinates are only an independent numerical witness to the recovered
    # Gram. Their eigensolver-dependent orientation is not retained as data.
    gram = np.asarray(result["gram_numeric"])
    values, vectors = np.linalg.eigh(gram)
    coordinates = vectors[:, -3:] * np.sqrt(values[-3:])
    rotation = np.array(((0., 1., 0.), (0., 0., -1.), (1., 0., 0.)))
    isometry_error = float(np.max(np.abs((coordinates@rotation)@(coordinates@rotation).T-gram)))
    controls = {
        "permutation": list(permutation),
        "relabelled_gram_covariance_error": relabel_error,
        "orthogonal_isometry_gram_error": isometry_error,
        "biased_scheduler": biased,
        "fixed_seam_scheduler": fixed,
        "erased_current_records": tampered,
    }
    checks = {
        "exact_native_rank_three_gram_recovered": result["rank_three_unit_gram_recovered"],
        "relabeling_transports_the_same_gram": relabel_error == 0,
        "orthogonal_coordinates_preserve_gram": isometry_error < 1e-12,
        "biased_scheduler_changes_response": biased["transition_matrix_exact"] != result["transition_matrix_exact"],
        "biased_scheduler_does_not_inherit_uniform_certificate": not biased["rank_three_unit_gram_recovered"],
        "fixed_seam_does_not_inherit_uniform_certificate": not fixed["rank_three_unit_gram_recovered"],
        "erased_currents_do_not_recover_geometry": not tampered["rank_three_unit_gram_recovered"],
    }
    return {
        "schema": "oph.observer_reconstruction.native_geometry.v1",
        "records": records, "recovery": result, "controls": controls,
        "source_execution": {
            "primitive": "model.observer.native.native",
            "native_branch_calls_main": len(records),
            "native_branch_calls_including_controls": 2*len(records)+720+24,
            "preparations_main": PORTS, "seams_main": 30,
            "coins_per_seam": 2, "operations_per_realized_history": 1,
            "native_quantum_compiler_words_executed": 0,
            "exact_record_ensemble": True, "finite_shots_sampled": False,
            "initialization_and_record_service_costs_supplied": True,
        },
        "scope": {
            "preparation_and_uniform_scheduler_supplied": True,
            "response_geometry_recovered_from_records": True,
            "geometry_is_field_grid_geometry": False,
            "physical_length_or_clock_identified": False,
            "propagation_or_quantum_field_law_inferred": False,
            "all_native_inputs_identified_by_twelve_preparations": False,
            "meaning": "exact finite ensemble response on twelve stated calibrations; no unique unrestricted repair law",
        },
        "checks": checks, "all_pass": all(checks.values()),
    }
