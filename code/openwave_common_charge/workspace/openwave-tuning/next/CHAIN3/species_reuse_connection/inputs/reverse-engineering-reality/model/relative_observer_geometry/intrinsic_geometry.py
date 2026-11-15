"""Exact single-carrier spatial readouts and their scheduling boundary.

The estimator sees labelled preparations and currents. It removes positive
seam exposure without importing coordinates. Even preparations also permit
exact recovery from finite native shots without a probability-law oracle.
The result concerns finite structural geometry, not physical distances.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from math import lcm
from pathlib import Path

import sympy as sp

from ..observer import native
from ..observer_reconstruction import native_geometry as ng
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
PORTS = 12
SOURCE_PATHS = (
    'model/relative_observer_geometry/intrinsic_geometry.py',
    'model/observer/native.py',
    'model/observer_reconstruction/native_geometry.py',
    'model/common_source/slow_geometry.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Geometry/SourceDerivedSpacetimeCarrier.lean',
)


def recover_structure(records):
    """Remove seam exposure from complete, possibly biased observation laws.

    Each preparation is 2 e_p. Per-preparation and per-seam schedules may
    differ, and coin probabilities may be biased or deterministic. Every
    preparation must expose the same seam support with positive mass. The
    estimator checks currents against that disclosed preparation. It does
    not use the hidden final loads or a source adjacency table.
    """
    if type(records) not in (list, tuple) or not records:
        raise ValueError('complete positive observation laws required')
    exposures = [{} for _ in range(PORTS)]
    currents = [{} for _ in range(PORTS)]
    masses = [F(0) for _ in range(PORTS)]
    seen = set()
    for row in records:
        if type(row) is not dict or set(row) != ng.RECORD_FIELDS:
            raise ValueError('only declared observation record fields accepted')
        p = row['preparation_port']
        seam = row['seam']
        coin, current = row['coin'], row['current']
        if type(p) is not int or p not in range(PORTS):
            raise ValueError('preparation port must be a native port')
        if (type(seam) not in (list, tuple) or len(seam) != 2
                or any(type(i) is not int or i not in range(PORTS) for i in seam)
                or seam[0] >= seam[1]):
            raise ValueError('one increasing orientation per seam required')
        a, b = seam
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError('a native coin is zero or one')
        if type(current) is not int or current not in (-1, 0, 1):
            raise ValueError('a native total-two current is integral')
        probability = ng._probability(row['probability'])
        key = p, a, b, coin
        if key in seen:
            raise ValueError('duplicate calibration outcome')
        seen.add(key)
        # In a 2e_p preparation every endpoint difference is even. The
        # native current equals 1_{p=a}-1_{p=b}, regardless of the coin.
        if current != int(p == a)-int(p == b):
            raise ValueError('current contradicts even-load native preparation')
        seam = a, b
        exposures[p][seam] = exposures[p].get(seam, F(0))+probability
        currents[p][seam] = currents[p].get(seam, F(0))+probability*current
        masses[p] += probability
    if any(mass != 1 for mass in masses):
        raise ValueError('each labelled preparation requires normalized complete laws')
    seams = tuple(sorted(exposures[0]))
    if not seams or any(set(exposure) != set(seams) for exposure in exposures):
        raise ValueError('the same positive seam support is required for every preparation')
    laplacian = [[F(0) for _ in range(PORTS)] for _ in range(PORTS)]
    for p in range(PORTS):
        for a, b in seams:
            conditional_current = currents[p][a, b]/exposures[p][a, b]
            laplacian[a][p] += conditional_current
            laplacian[b][p] -= conditional_current
    checks, projector = ng._band_certificate(laplacian)
    symmetric = all(laplacian[i][j] == laplacian[j][i]
                    for i in range(PORTS) for j in range(PORTS))
    graph = symmetric and all(laplacian[i][i] == 5 for i in range(PORTS))
    graph = graph and all(laplacian[i][j] in (0, -1)
                         for i in range(PORTS) for j in range(PORTS) if i != j)
    certified = graph and all(checks.values())
    return {
        'record_count': len(records), 'seams': [list(seam) for seam in seams],
        'laplacian_exact': [[str(x) for x in row] for row in laplacian],
        'seam_exposures_exact': [
            {f'{a},{b}': str(exposure[a, b]) for a, b in seams}
            for exposure in exposures],
        'rank_three_structural_gram_certified': certified,
        'band_certificate': checks,
        'gram_qsqrt5': ng._encode_q5([
            [ng._qscale(x, 4) for x in row] for row in projector]) if certified else None,
        'uniform_scheduler_or_fair_coins_required_for_this_readout': False,
        'full_response_slow_band_identified_for_the_actual_scheduler': False,
        'hidden_final_loads_or_source_coordinates_used': False,
    }


def execute_biased_calibration(seams):
    """Execute native endpoints with positive preparation-dependent exposures.

    Coin zero alone is used on even preparation labels; odd labels retain
    coins with probabilities 1/3 and 2/3. These are alternative declared
    protocols, not a change to the model's production scheduler.
    """
    records = []
    for p in range(PORTS):
        numerators = [1+(index+3*p) % 7 for index in range(len(seams))]
        total = sum(numerators)
        coins = ((0, F(1)),) if p % 2 == 0 else ((0, F(1, 3)), (1, F(2, 3)))
        for seam, weight in zip(seams, numerators, strict=True):
            for coin, coin_probability in coins:
                _, (_, current) = native.native(tuple(2*int(i == p) for i in range(PORTS)), seam, coin)
                records.append({'preparation_port': p, 'seam': list(seam),
                    'coin': coin, 'current': current,
                    'probability': str(F(weight, total)*coin_probability)})
    return records


def execute_single_shot_calibration(seams, coin_value=None):
    """Execute one native branch per preparation/seam, with fixed varied coins."""
    if coin_value is not None and (type(coin_value) is not int or coin_value not in (0, 1)):
        raise ValueError('fixed coin must be zero or one')
    records = []
    for p in range(PORTS):
        loads = tuple(2*int(i == p) for i in range(PORTS))
        for index, seam in enumerate(seams):
            coin = (p+index) % 2 if coin_value is None else coin_value
            _, (_, current) = native.native(loads, seam, coin)
            records.append({'preparation_port': p, 'seam': list(seam),
                            'coin': coin, 'current': current})
    return records


def recover_from_shots(records):
    """Recover exactly in the declared native model without probability inputs.

    Empirical frequencies are only bookkeeping weights. Conditional currents
    have zero variance on the known even preparations, so any observed coin
    frequencies give the exact conditional mean. Native rule, preparation and
    faithful record access remain premises; this is not a lab-noise theorem.
    """
    if type(records) not in (list, tuple) or not records:
        raise ValueError('nonempty native calibration shots required')
    counts, currents = Counter(), {}
    total = Counter()
    for row in records:
        if type(row) is not dict or set(row) != ng.RECORD_FIELDS-{'probability'}:
            raise ValueError('shots contain observations, not probability laws')
        p, seam = row['preparation_port'], row['seam']
        coin, current = row['coin'], row['current']
        if (type(p) is not int or p not in range(PORTS)
                or type(seam) not in (list, tuple) or len(seam) != 2
                or any(type(i) is not int or i not in range(PORTS) for i in seam)
                or type(coin) is not int or coin not in (0, 1)
                or type(current) is not int):
            raise ValueError('native preparation, seam, coin and current required')
        key = p, *seam, coin
        if key in currents and currents[key] != current:
            raise ValueError('inconsistent currents for a deterministic native calibration')
        currents[key] = current
        counts[key] += 1
        total[p] += 1
    weighted = [{'preparation_port': p, 'seam': [a, b], 'coin': coin,
                 'current': currents[p, a, b, coin],
                 'probability': str(F(count, total[p]))}
                for (p, a, b, coin), count in sorted(counts.items())]
    result = recover_structure(weighted)
    result['shot_count'] = len(records)
    result['probability_law_oracle_used'] = False
    result['empirical_frequencies_are_only_bookkeeping_weights'] = True
    result['conditional_current_variance_zero_under_native_model'] = True
    result['laboratory_noise_or_preparation_error_bounded'] = False
    return result


def _qinverse(x):
    denominator = x[0]*x[0]-5*x[1]*x[1]
    if denominator == 0:
        raise ZeroDivisionError('zero quadratic-field element')
    return x[0]/denominator, -x[1]/denominator


def qrank(matrix):
    """Gaussian rank over Q(sqrt(5)), with no floating threshold."""
    rows = [[tuple(x) for x in row] for row in matrix]
    if not rows or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('nonempty rectangular matrix required')
    pivot = 0
    zero = F(0), F(0)
    for column in range(len(rows[0])):
        index = next((i for i in range(pivot, len(rows)) if rows[i][column] != zero), None)
        if index is None:
            continue
        rows[pivot], rows[index] = rows[index], rows[pivot]
        inv = _qinverse(rows[pivot][column])
        rows[pivot] = [ng._qmul(x, inv) for x in rows[pivot]]
        for i in range(pivot+1, len(rows)):
            factor = rows[i][column]
            rows[i] = [ng._qadd(x, ng._qscale(ng._qmul(factor, y), -1))
                       for x, y in zip(rows[i], rows[pivot], strict=True)]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def scheduling_control(laplacian, seam):
    """Exactly certify the slow-band split under one doubled seam exposure."""
    a, b = seam
    edge = [F(int(i == a)-int(i == b)) for i in range(PORTS)]
    weighted = [[laplacian[i][j]+edge[i]*edge[j]
                 for j in range(PORTS)] for i in range(PORTS)]
    lam = F(5), F(-1)
    shifted = [[ng._qadd((weighted[i][j], F(0)), ng._qscale(lam, -int(i == j)))
                for j in range(PORTS)] for i in range(PORTS)]
    multiplicity = PORTS-qrank(shifted)
    if multiplicity != 2:
        raise AssertionError('one doubled seam must leave exactly two slow eigenvectors')
    # L + bb^T >= L and bb^T annihilates constants. Hence there is no
    # nonconstant eigenvalue below lambda_min(L), while exact nullity two
    # above proves that this smallest positive eigenvalue is attained twice.
    return {
        'doubled_seam': list(seam),
        'scheduler': 'weight 2/31 on one seam, 1/31 on each other seam',
        'actual_mean_response': 'I-(L+b*b^T)/62; b=e_a-e_b',
        'smallest_positive_weighted_laplacian_eigenvalue_qsqrt5': ['5', '-1'],
        'smallest_positive_weighted_laplacian_multiplicity_exact': multiplicity,
        'slowest_nonconstant_actual_response_qsqrt5': ['57/62', '1/62'],
        'perturbation_positive_semidefinite_exact': True,
        'normalized_long_attempt_response_gram_rank': multiplicity,
        'structural_gram_rank_after_exposure_removal': 3,
        'time_parameter': 'uniformly indexed native attempts with the declared biased seam scheduler',
    }


def _inverse(matrix):
    size = len(matrix)
    rows = [list(row)+[F(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        if pivot is None:
            raise ValueError('singular matrix')
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [x/scale for x in rows[column]]
        for i in range(size):
            if i != column:
                factor = rows[i][column]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[column], strict=True)]
    return [row[size:] for row in rows]


def inverse_response_metric(laplacian):
    """Compute all-mode inverse-response geometry on the constant-free space.

    For graph L, G=(L+J/12)^-1-J/12 is the positive inverse on 1^perp.
    R_ab=(e_a-e_b)^T G (e_a-e_b) is effective resistance. Its square root
    is a Euclidean distance with centered Gram G of rank eleven. This is
    another mathematical readout of the same response, not a rod apparatus.
    """
    inverse = _inverse([[x+F(1, PORTS) for x in row] for row in laplacian])
    gram = [[x-F(1, PORTS) for x in row] for row in inverse]
    centered_identity = [[sum(laplacian[i][k]*gram[k][j] for k in range(PORTS))
                          for j in range(PORTS)] for i in range(PORTS)]
    expected = [[F(i == j)-F(1, PORTS) for j in range(PORTS)] for i in range(PORTS)]
    if centered_identity != expected or any(sum(row) for row in gram):
        raise AssertionError('inverse response must invert exactly off the constant mode')
    rank = qrank([[(x, F(0)) for x in row] for row in gram])
    resistance = [[gram[i][i]+gram[j][j]-2*gram[i][j]
                   for j in range(PORTS)] for i in range(PORTS)]
    if any(resistance[i][j] > resistance[i][k]+resistance[k][j]
           for i in range(PORTS) for j in range(PORTS) for k in range(PORTS)):
        raise AssertionError('effective resistance triangle inequality failed')
    shells = Counter(str(resistance[i][j]) for i in range(PORTS) for j in range(i+1, PORTS))
    if rank != 11 or dict(shells) != {'11/30': 30, '7/15': 30, '1/2': 6}:
        raise AssertionError('unexpected exact carrier inverse-response metric')
    return {
        'gram_exact': [[str(x) for x in row] for row in gram],
        'squared_euclidean_distance_exact': [[str(x) for x in row] for row in resistance],
        'squared_distance_shell_counts': dict(sorted(shells.items())),
        'centered_gram_rank_exact': rank,
        'L_times_G_equals_constant_free_identity_exact': True,
        'resistance_triangle_inequality_all_triples_exact': True,
        'source': 'inverse of the recovered unweighted response on the constant-free space',
        'physical_distance_or_dimension_selected': False,
    }


def wave_density_metric(spatial_tensor, measure_density=1, time_weight=1):
    """Recover a comparison Lorentz metric from a four-dimensional action.

    Inputs are exact rational coefficients at one point, not inferred native
    continuum coefficients. The continuum argument in README obtains a
    measurable positive spatial tensor from the already proved local form.
    A smooth tensor, radial geometry and actual wave apparatus are distinct
    obligations. The time weight names the chosen wave-clock coordinate.
    """
    def rational(value):
        if type(value) not in (int, F) and not isinstance(value, sp.Rational):
            raise ValueError('exact rational coefficients required')
        return sp.Rational(value)

    try:
        rows = [list(row) for row in spatial_tensor]
    except TypeError as exc:
        raise ValueError('a three-dimensional spatial tensor is required') from exc
    if len(rows) != 3 or any(len(row) != 3 for row in rows):
        raise ValueError('a three-dimensional spatial tensor is required')
    a = sp.Matrix([[rational(value) for value in row] for row in rows])
    mu, weight = rational(measure_density), rational(time_weight)
    if mu <= 0 or weight <= 0:
        raise ValueError('positive measure density and time weight required')
    if a != a.T or any(a[:k, :k].det() <= 0 for k in range(1, 4)):
        raise ValueError('a symmetric positive definite spatial tensor is required')
    coefficient = mu*sp.diag(-weight, a)
    volume = sp.sqrt(-coefficient.det())
    metric = sp.simplify(volume*coefficient.inv())
    inverse = sp.simplify(coefficient/volume)
    if sp.simplify(metric.det()-coefficient.det()) != 0:
        raise ArithmeticError('four-dimensional density determinant identity failed')
    if sp.simplify(volume*inverse-coefficient) != sp.zeros(4):
        raise ArithmeticError('metric did not reproduce the complete action density')
    return {'action_density': coefficient, 'metric': metric,
            'inverse_metric': inverse, 'metric_volume_density': volume}


def product_wave_metric_control():
    """Exact tensor/control checks, separate from finite carrier readbacks."""
    # This anisotropic rational tensor is a test coefficient, not a fitted
    # value of the subsequential native continuum tensor.
    a = sp.Matrix([[1, 0, 0], [0, 2, 1], [0, 1, 3]])
    data = wave_density_metric(a.tolist(), F(3, 2))
    doubled = wave_density_metric(a.tolist(), 3)
    rate = sp.Rational(3, 2)
    rescaled = wave_density_metric(a.tolist(), rate*sp.Rational(3, 2), 1/rate**2)
    jacobian = sp.diag(rate, 1, 1, 1)  # old coordinate t = rate * new coordinate
    if sp.simplify(rescaled['metric']-jacobian.T*data['metric']*jacobian) != sp.zeros(4):
        raise ArithmeticError('clock reparameterization must transform the same metric')
    if sp.simplify(doubled['metric']-2*data['metric']) != sp.zeros(4):
        raise ArithmeticError('action normalization must retain its conformal factor')
    encode = lambda matrix: [[str(sp.simplify(value)) for value in row]
                             for row in matrix.tolist()]
    return {
        'coefficient_fixture': 'anisotropic rational product tensor; not measured native coefficients',
        'spatial_tensor': encode(a),
        'measure_density': '3/2',
        'action_density': encode(data['action_density']),
        'metric': encode(data['metric']),
        'metric_volume_density': str(data['metric_volume_density']),
        'signature': [-1, 1, 1, 1],
        'complete_action_density_recovered_exactly': True,
        'clock_reparameterization_rate': str(rate),
        'clock_reparameterization_is_tensor_covariant': True,
        'doubling_measure_preserves_cone_but_doubles_metric': True,
        'continuum_input': 'regular strongly local H1 form and chosen product wave action',
        'coefficient_regularities_proved': ['measurable', 'uniformly elliptic in local charts'],
        'smooth_coefficients_or_curvature_proved': False,
        'native_radial_transport_or_lab_proper_time_proved': False,
        'observed_spacetime_or_Einstein_dynamics_proved': False,
    }


@lru_cache(maxsize=3)
def _native_density_rows(level):
    """Integer kernels from total-one endpoint records, including all holds.

    This is the unchanged primitive with a disclosed total-one preparation,
    not a tagged-particle assumption about the fixed-number density field.
    Empty seams give two zero-current records; occupied seams are executed.
    """
    from ..federated_observer import source
    graph = source.federation(level)
    denominator = 2*len(graph.seams)
    rows = [Counter({p: denominator-2*len(graph.neighbours[p])})
            for p in range(graph.ports)]
    for a, b in graph.seams:
        for origin, before in ((a, (1, 0)), (b, (0, 1))):
            for coin in (0, 1):
                after, (_, current) = native.native(before, (0, 1), coin)
                decoded = before[0]-current, before[1]+current
                if after != decoded or sum(after) != 1 or min(after) < 0:
                    raise ArithmeticError('native current failed total-one reconstruction')
                rows[origin][(a, b)[after.index(1)]] += 1
    for p, row in enumerate(rows):
        expected = {p: denominator-len(graph.neighbours[p]),
                    **{q: 1 for q in graph.neighbours[p]}}
        if dict(row) != expected or sum(row.values()) != denominator:
            raise ArithmeticError('native records disagree with the full-port density kernel')
    return graph, tuple(tuple(sorted(row.items())) for row in rows), denominator


def _native_chart_dictionary(level):
    """Exact barycentric probes with a continuous cutoff on source face zero."""
    from ..physical_space_time import uniform_green as uniform
    graph, _, _ = _native_density_rows(level)
    vertices = ((0, 0, 0, 0, 0, 0), (2, 0, 0, 0, 0, 0),
                (0, 2, 0, 0, 0, 0))
    faces = ((0, 1, 2),)
    for _ in range(level):
        vertices, faces = uniform.refine(vertices, faces)
    _, whole_faces = uniform.affine_prefix(level)
    actual = {tuple(sorted((a//12, b//12))) for a, b in graph.inter_seams}
    if actual != set(uniform.dual_edges(whole_faces)):
        raise ArithmeticError('chart ancestry does not match the actual native seams')
    probes, weights = [], []
    for face in faces:
        x, y = (sum(F(vertices[v][i], 2**(level+1)) for v in face)/3
                for i in (0, 1))
        distance = min(x, y, 1-x-y)
        cutoff = min(F(1), 6*distance)
        probes.append((cutoff*x, cutoff*y))
        weights.append(max(F(0), 12*distance-3))
    probes.extend([(F(0), F(0))]*(graph.cells-len(probes)))
    weights.extend([F(0)]*(graph.cells-len(weights)))
    if sum(weights) <= 0:
        raise ArithmeticError('the declared interior chart has no sampled weight')
    return tuple(probes), tuple(weights)


def native_chart_increment(level=1, attempts=2):
    """Read a finite averaged tensor from full-port conditional increments.

    The time is attempts/(1260*16**level), not seconds or a modular parameter.
    These exact finite observations are not claimed to approximate the limiting
    tensor at this cap. Their continuum identification requires the ordered
    limits in README, with bounded H1 probes and the full constant mode.
    """
    if type(level) is not int or level not in (1, 2):
        raise ValueError('a native chart level of one or two is required')
    if type(attempts) is not int or not 1 <= attempts <= 32:
        raise ValueError('one through thirty-two native attempts required')
    graph, rows, denominator = _native_density_rows(level)
    probes, weights = _native_chart_dictionary(level)
    scale = lcm(*(v.denominator for pair in probes for v in pair))
    first = tuple(int(pair[0]*scale) for pair in probes for _ in range(12))
    second = tuple(int(pair[1]*scale) for pair in probes for _ in range(12))
    features = [first, second, tuple(x*x for x in first),
                tuple(x*y for x, y in zip(first, second, strict=True)),
                tuple(y*y for y in second)]
    for _ in range(attempts):
        features = [tuple(sum(count*values[q] for q, count in row) for row in rows)
                    for values in features]
    power = denominator**attempts
    time = F(attempts, denominator*4**level)
    tensor = [[F(0), F(0)], [F(0), F(0)]]
    for p in range(graph.ports):
        weight = weights[p//12]
        if not weight:
            continue
        for i, j, slot in ((0, 0, 2), (0, 1, 3), (1, 1, 4)):
            x, y = (first, second)[i][p], (first, second)[j][p]
            numerator = (features[slot][p]-x*features[j][p]
                         -y*features[i][p]+x*y*power)
            tensor[i][j] += weight*F(numerator, scale**2*power)
    norm = 24*sum(weights)*time  # twelve starting ports per cell and 2*time
    tensor = [[x/norm for x in row] for row in tensor]
    tensor[1][0] = tensor[0][1]
    if tensor[0][0] <= 0 or tensor[0][0]*tensor[1][1] <= tensor[0][1]**2:
        raise ArithmeticError('the executed two-coordinate increment tensor is degenerate')
    encode = lambda matrix: [[str(v) for v in row] for row in matrix]
    return {'level': level, 'native_attempts': attempts,
            'density_parameter': str(time), 'cells': graph.cells,
            'ports': graph.ports, 'native_seams': len(graph.seams),
            'nonzero_chart_weights': sum(w > 0 for w in weights),
            'averaged_screen_tensor': encode(tensor),
            'kernel': 'full-port native records; no cell reset or mean-zero projection',
            'preparation': 'disclosed total one; chart-weighted starting cell and uniform port',
            'finite_value_is_identified_with_continuum_tensor': False}


def native_increment_control():
    """Native chart values, source-record oracles and false-limit controls."""
    graph, rows, denominator = _native_density_rows(1)
    probes, _ = _native_chart_dictionary(1)
    cell = 3
    # Independently compose endpoint records for each of the twelve starting
    # ports. Remote seam/coin records retain their exact hold multiplicities.
    direct = sp.zeros(2)
    time = sp.Rational(2, denominator*4)
    endpoints = Counter()
    for origin in range(12*cell, 12*(cell+1)):
        for middle, first_count in rows[origin]:
            for target, second_count in rows[middle]:
                multiplicity = first_count*second_count
                endpoints[target//12] += multiplicity
                delta = sp.Matrix([sp.Rational(probes[target//12][i]-probes[cell][i])
                                   for i in range(2)])
                direct += multiplicity*(delta*delta.T)/(24*denominator**2*time)
    early = [native_chart_increment(level, n) for level in (1, 2) for n in (1, 2)]
    if direct != sp.Matrix(early[1]['averaged_screen_tensor']):
        raise ArithmeticError('native endpoint records disagree with polarized heat readout')
    bare = sp.Matrix([[sp.Rational(1, 36), -sp.Rational(1, 72)],
                      [-sp.Rational(1, 72), sp.Rational(1, 36)]])
    for row in early:
        factor = 1-sp.Rational(row['native_attempts']-1, 2*row['native_seams'])
        if sp.Matrix(row['averaged_screen_tensor']) != factor*bare:
            raise ArithmeticError('the exact early native chart response changed')
    # Resetting to a uniform port after each attempt fabricates a two-cross-
    # seam path: the true matching needs an intervening internal move.
    neighbors = {c: set() for c in range(graph.cells)}
    for a, b in graph.inter_seams:
        neighbors[a//12].add(b//12); neighbors[b//12].add(a//12)
    target = min(q for p in neighbors[cell] for q in neighbors[p]
                 if q != cell and q not in neighbors[cell])
    reset_paths = len(neighbors[cell] & neighbors[target])
    reset_probability = F(reset_paths, (12*denominator)**2)
    if endpoints[target] != 0 or reset_probability <= 0:
        raise ArithmeticError('the finite cell-reset countercontrol lost its distinction')
    # The one-attempt-first value persists with refinement. Actual continuum
    # energy has the strict upper bound 15/19 of this bare tensor, from the
    # existing terminal-star conductances and their native recovery maps.
    from ..physical_space_time import local_continuum as local
    conductances = local.weighted_dual_network(graph)[2]
    if any(w > F(15, 19) for _, _, w in conductances):
        raise ArithmeticError('the native terminal conductance upper bound changed')
    fixed_time = [early[0], native_chart_increment(2, 16)]
    if fixed_time[0]['density_parameter'] != fixed_time[1]['density_parameter']:
        raise ArithmeticError('refinement changed the declared observation separation')
    measured = sp.Matrix(fixed_time[1]['averaged_screen_tensor'])
    product = sp.diag(1, measured)
    metric = wave_density_metric(product.tolist(), F(1, 10))
    encode = lambda matrix: [[str(sp.simplify(v)) for v in row] for row in matrix.tolist()]
    return {'early_native_chart_readouts': early,
        'fixed_density_parameter_readouts': fixed_time,
        'two_attempt_endpoint_record_oracle_exact': True,
        'two_attempt_native_record_tapes_per_starting_port': denominator**2,
        'zero_current_records_included': True,
        'false_cell_reset': {'initial_cell': cell, 'terminal_cell': target,
            'actual_probability': '0', 'reset_probability': str(reset_probability)},
        'false_Schur_kernel': local.non_markov_control(),
        'bare_chart_tensor': encode(bare),
        'continuum_chart_tensor_lower_relative_to_bare': '3/7',
        'continuum_chart_tensor_upper_relative_to_bare': '15/19',
        'one_attempt_first_limit_equals_bare_not_effective_tensor': True,
        'finite_averaged_product_tensor': encode(product),
        'product_measure_density_in_barycentric_chart': '1/10',
        'finite_averaged_comparison_metric': encode(metric['metric']),
        'ordered_limits': ['spatial subsequence at fixed positive density parameter',
                           'density parameter to zero against bounded H1 probes',
                           'chart averaging radius to zero at Lebesgue points'],
        'limiting_tensor_pointwise_error_from_finite_run_bounded': False,
        'radial_factor_and_active_wave_law_are_declared': True,
        'native_time_is_identified_with_wave_or_lab_time': False}


def constant_optical_separation(spatial_tensor, displacement, time_weight=1):
    """Exact constant-coefficient optical duality, not a rough-metric ruler.

    For C=mu*diag(-w,a), the optical tensor is w*a^-1.  A coefficient
    averaged from finite native records can test this identity, but its
    constant extension is not the continuum coefficient field.
    """
    data = wave_density_metric(spatial_tensor, time_weight=time_weight)
    try:
        values = list(displacement)
    except TypeError as exc:
        raise ValueError('three exact rational displacement components required') from exc
    if len(values) != 3 or any(type(v) not in (int, F) and
                               not isinstance(v, sp.Rational) for v in values):
        raise ValueError('three exact rational displacement components required')
    vector = sp.Matrix([sp.Rational(v) for v in values])
    optical = sp.simplify(data['metric'][1:, 1:]/(-data['metric'][0, 0]))
    squared = sp.simplify((vector.T*optical*vector)[0])
    dual = optical*vector/sp.sqrt(squared) if squared else sp.zeros(3, 1)
    return {'optical_tensor': optical, 'distance_squared': squared,
            'maximizing_covector': sp.simplify(dual)}


def radial_front_readout(wave_parameter):
    """Exact product-wave comparison on its restored angular constant mode.

    Declare radial interval (0,4), Dirichlet endpoints, initial interval
    (1,3/2) and detector (2,5/2). For 0<=t<=1 no reflected wave reaches
    this detector. Values are unnormalized L2 pairings; angular mass is one.
    """
    if (type(wave_parameter) not in (int, F) and
            not isinstance(wave_parameter, sp.Rational)):
        raise ValueError('an exact wave parameter in [0,1] is required')
    t = F(wave_parameter)
    if not 0 <= t <= 1:
        raise ValueError('an exact wave parameter in [0,1] is required')
    overlap = max(F(0), min(F(5, 2), F(3, 2)+t)-max(F(2), F(1)+t))
    cosine = overlap/2  # the left-moving copy has no overlap in this window
    sine = max(F(0), t-F(1, 2))**2/4
    return {'wave_parameter': str(t), 'cosine_pairing': str(cosine),
            'retarded_sine_pairing': str(sine)}


def native_optical_front_control():
    """Native coefficient duality and controls on what a wave front means."""
    native_tensor = native_chart_increment(1, 2)
    a = sp.diag(1, sp.Matrix(native_tensor['averaged_screen_tensor']))
    displacement = sp.Matrix([1, 2, -1])
    separation = constant_optical_separation(a.tolist(), displacement)
    squared, dual = separation['distance_squared'], separation['maximizing_covector']
    if sp.simplify((dual.T*a*dual)[0]-1) != 0:
        raise ArithmeticError('optical covector must saturate the energy constraint')
    if sp.simplify((dual.T*displacement)[0]-sp.sqrt(squared)) != 0:
        raise ArithmeticError('optical dual readout must equal its separation')
    metric = wave_density_metric(a.tolist(), F(1, 10))['metric']
    null_vector = sp.Matrix([sp.sqrt(squared), *displacement])
    if sp.simplify((null_vector.T*metric*null_vector)[0]) != 0:
        raise ArithmeticError('optical separation disagrees with the wave metric')
    change = sp.Matrix([[1, 0, 0], [0, 2, 1], [0, 1, 1]])
    transformed = constant_optical_separation(
        (change.inv()*a*change.inv().T).tolist(), change.inv()*displacement)
    rate = sp.Rational(3, 2)
    rescaled = constant_optical_separation(a.tolist(), displacement, 1/rate**2)
    if transformed['distance_squared'] != squared or rescaled['distance_squared'] != squared/rate**2:
        raise ArithmeticError('optical separation lost chart or clock covariance')

    # This full-port eigenmode is compactly supported on carrier zero and
    # vanishes at its external seam ports.  A remote detector is dark even
    # though a positive local preparation crosses the same seam.
    graph, rows, denominator = _native_density_rows(0)
    profile = (0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0)+(0,)*(graph.ports-12)
    transfer = tuple(sum(count*profile[q] for q, count in row) for row in rows)
    if any(value != (denominator-6)*profile[p] for p, value in enumerate(transfer)):
        raise ArithmeticError('the actual full-port dark mode is not an eigenmode')
    crossing = next((p, q) for p, q in graph.inter_seams if p < 12 <= q)
    origin, target = crossing
    positive = sum(F(count, 12*denominator) for p, count in rows[target] if p < 12)
    if profile[origin] != 0 or positive <= 0 or transfer[target] != 0:
        raise ArithmeticError('the dark-mode control did not retain a transmitting seam')
    altered = constant_optical_separation((4*a).tolist(), displacement)
    if altered['distance_squared'] != squared/4:
        raise ArithmeticError('null-line representative control lost its length distinction')
    return {
        'native_coefficient_readout': native_tensor,
        'declared_radial_coefficient': '1',
        'constant_coefficient_displacement': [str(v) for v in displacement],
        'constant_coefficient_optical_distance_squared': str(squared),
        'optical_dual_energy_saturated_exactly': True,
        'metric_null_vector_identity_exact': True,
        'spatial_chart_covariance_exact': True,
        'wave_coordinate_change_t_equals_rate_times_tau': str(rate),
        'rescaled_optical_distance_squared': str(rescaled['distance_squared']),
        'density_measure_does_not_change_optical_tensor': True,
        'dark_native_mode': {
            'ports': graph.ports, 'native_seams': len(graph.seams),
            'first_carrier_profile': list(profile[:12]), 'other_ports_zero': True,
            'full_port_laplacian_eigenvalue': '6',
            'full_port_eigenvector_identity_exact': True,
            'remote_detector_port': target, 'transmitting_seam': list(crossing),
            'signed_mode_response_at_remote_detector_for_every_wave_parameter': '0',
            'uniform_first_carrier_one_attempt_probability_at_detector': str(positive),
            'continuous_matrix_cosine_cross_seam_t_squared_coefficient': '1/2',
            'finite_matrix_wave_has_a_positive_hard_arrival_delay': False},
        'declared_interval_constant_mode_front': {
            'radial_interval': ['0', '4'], 'boundary': 'Dirichlet',
            'initial_interval': ['1', '3/2'], 'detector_interval': ['2', '5/2'],
            'essential_set_distance': '1/2', 'angular_mode': 'restored normalized constant',
            'readouts': [radial_front_readout(t) for t in (F(1, 4), F(1, 2), F(5, 8), F(3, 4))],
            'mean_zero_angular_projection_erases_this_control': True},
        'null_set_representative_control': {
            'alteration': 'multiply the constant tensor by four only on the displacement line',
            'energy_form_heat_wave_and_intrinsic_distance_unchanged': True,
            'naive_curve_length_squared_on_that_line': str(altered['distance_squared']),
            'unmodified_constant_optical_distance_squared': str(squared)},
        'continuum_front_scope': {
            'complete_local_cosine_and_retarded_sine_onset_equals_intrinsic_set_distance': True,
            'onset_is_an_infimum_not_a_claim_of_response_at_the_front': True,
            'finite_coefficient_average_is_the_limiting_metric': False,
            'one_fixed_detector_must_see_every_signal': False,
            'smooth_null_geodesics_or_fastest_native_transport_identified': False,
            'radial_design_active_wave_and_clock_contracts_required': True}}


def _radial_adaptive_history_control(graph, laplacian, models):
    """Complete comparison histories on the invariant constant-radial code.

    Native timing tapes are executed; the spectral wave and angular
    instruments are comparison operations, not additional source pulses.
    The reference columns are retained through every stopped branch.
    """
    import numpy as np

    size = graph.ports
    eigenvalues, axes = np.linalg.eigh(np.array(laplacian, dtype=float))
    if abs(eigenvalues[0]) > 1e-10 or eigenvalues[1] <= 0:
        raise ArithmeticError('the angular comparison lost its unique constant mode')
    eigenvalues[0] = 0  # Exact conserved mode, not sqrt of roundoff.
    seam = next(edge for edge in graph.seams if edge[0] % 2 != edge[1] % 2)
    a, b = seam
    initial = np.zeros((size, 2), complex)
    initial[a, 0], initial[b, 1] = 1/np.sqrt(2), 1j/np.sqrt(2)
    mask = np.array([p % 2 == 0 for p in range(size)])
    swap = np.arange(size); swap[a], swap[b] = b, a

    # Same genuine native two-attempt clock, including zero-current marks.
    seams = tuple(tuple(edge) for edge in recover_structure(ng.execute_calibration())['seams'])
    before = (1,)+(0,)*(PORTS-1)
    marks, tape_rows, zero_incidents = Counter(), [], 0
    for edge in seams:
        for coin in (0, 1):
            middle, (_, current) = native.native(before, edge, coin)
            for second in seams:
                for next_coin in (0, 1):
                    final, (_, next_current) = native.native(middle, second, next_coin)
                    flags = int(0 in edge), int(0 in second)
                    marks[flags] += 1
                    zero_incidents += int(flags[0] and current == 0)+int(flags[1] and next_current == 0)
                    tape_rows.append((edge, coin, current, second, next_coin, next_current, final))
    if marks != Counter({(0, 0): 2500, (0, 1): 500, (1, 0): 500, (1, 1): 100}):
        raise ArithmeticError('the native comparison timer changed its complete law')
    weights = [F(sum(count for flags, count in marks.items() if sum(flags) == n), len(tape_rows))
               for n in range(3)]

    def histories(start, basis, spectrum, radial_copies, feedback=True):
        branches, out = {(): start}, []
        for step in range(3):
            out.append(branches)
            if step == 2:
                break
            following = {}
            for bits, state in branches.items():
                time = F(1, 2) if bits == (1,) else F(1, 4)
                advanced = basis@(np.exp(-1j*float(time)*np.sqrt(spectrum))[:, None]
                                  *(basis.conj().T@state))
                for bit in (0, 1):
                    terminal = advanced*np.tile(mask if bit == 0 else ~mask, radial_copies)[:, None]
                    if step == 0 and bit == 1 and feedback:
                        terminal = terminal.reshape(radial_copies, size, 2)[:, swap, :].reshape(-1, 2)
                    following[bits+(bit,)] = terminal
            branches = following
        return out

    base = histories(initial, axes, eigenvalues, 1)
    no_feedback = histories(initial, axes, eigenvalues, 1, False)
    probability = lambda state: float(np.linalg.norm(state)**2)
    feedback_tv = sum(float(weights[n])*sum(abs(probability(state)-probability(no_feedback[n][bits]))
                       for bits, state in branches.items())/2 for n, branches in enumerate(base))
    lift = lambda state: np.tile(state, (3, 1))/np.sqrt(3)
    checks = []
    for model in models:
        radial = np.array(model['radial_laplacian'], dtype=float)
        radial_values, radial_axes = np.linalg.eigh(radial)
        if abs(radial_values[0]) > 1e-12 or radial_values[1] <= 0:
            raise ArithmeticError('the radial comparison lost its unique constant mode')
        radial_values[0] = 0
        states = histories(lift(initial), np.kron(radial_axes, axes),
                           (radial_values[:, None]+eigenvalues).reshape(-1), 3)
        distance, mass = 0., 0.
        for n, branches in enumerate(states):
            for bits, state in branches.items():
                expected = lift(base[n][bits])
                # Half-trace bound on unnormalized pure work/reference outputs.
                distance += float(weights[n])*.5*(np.linalg.norm(state)+np.linalg.norm(expected))*np.linalg.norm(state-expected)
                mass += float(weights[n])*probability(state)
        if distance > 1e-10 or abs(mass-1) > 1e-11:
            raise ArithmeticError('an adaptive angular history distinguishes the constant radial code')
        checks.append({'radial_model': model['name'],
            'complete_work_reference_record_distance_upper_numerical': float(distance),
            'complete_output_probability': float(mass)})
    if feedback_tv < 1e-4 or not zero_incidents:
        raise ArithmeticError('the adaptive comparison lacks genuine feedback or null incidents')
    return {'complete_native_tapes': len(tape_rows),
        'complete_native_tape_sha256': hashlib.sha256(canonical(tape_rows).encode()).hexdigest(),
        'complete_tape_and_field_outcomes': sum(count*2**sum(flags) for flags, count in marks.items()),
        'clock_pattern_counts': [{'incident_marks': list(flags), 'count': count} for flags, count in sorted(marks.items())],
        'zero_current_incidents_retained': zero_incidents, 'reference_dimension': 2,
        'feedback_native_seam': list(seam), 'first_model_duration': '1/4',
        'second_duration_by_first_bit': ['1/4', '1/2'],
        'complete_angular_record_law': [{'incident_ticks': n, 'field_bits': list(bits),
            'probability': float(weights[n])*probability(state)} for n, branches in enumerate(base)
                                      for bits, state in branches.items()],
        'radial_model_comparisons': checks, 'feedback_removal_record_tv': float(feedback_tv),
        'theorem_covers_arbitrary_angular_reference_inputs_in_constant_radial_code': True,
        'native_timer_pauses_during_comparison_services': True,
        'source_pulses_for_comparison_operations_executed': False,
        'full_product_Gibbs_or_nonconstant_radial_inputs_covered': False}


def radial_nonidentifiability_control():
    """Connected local extensions agree on the retained angular readout.

    The source Laplacian comes from all total-one native endpoint records.
    The radial graphs and their normalized average are comparison choices,
    not extra native seams or an ordinary boundary trace. Their discarded
    radial contrasts are explicitly measurable within each enlarged model.
    """
    graph, rows, denominator = _native_density_rows(0)
    size = graph.ports
    lap = sp.SparseMatrix(size, size, {
        (p, q): (denominator if p == q else 0)-count
        for p, row in enumerate(rows) for q, count in row})
    lift = sp.SparseMatrix(3*size, size,
                          {(r*size+p, p): 1 for r in range(3) for p in range(size)})
    if lift.T*lift != 3*sp.eye(size):
        raise ArithmeticError('the normalized radial average lost its isometry')
    contrast = sp.Matrix([1, -1, 0])
    lifted_contrast = sp.Matrix([value for value in contrast for _ in range(size)])
    models = []
    for name, edges, weight in (
        ('path_weight_one', ((0, 1), (1, 2)), 1),
        ('path_weight_four', ((0, 1), (1, 2)), 4),
        ('cycle_weight_one', ((0, 1), (1, 2), (0, 2)), 1)):
        radial = sp.zeros(3)
        for a, b in edges:
            radial[a, a] += weight; radial[b, b] += weight
            radial[a, b] -= weight; radial[b, a] -= weight
        if radial*sp.ones(3, 1) != sp.zeros(3, 1) or radial.rank() != 2:
            raise ArithmeticError('the radial comparison graph must be connected and conservative')
        entries = {(r*size+i, r*size+j): value for r in range(3)
                   for (i, j), value in lap.todok().items()}
        for a in range(3):
            for b in range(3):
                for p in range(size):
                    key = a*size+p, b*size+p
                    entries[key] = entries.get(key, 0)+radial[a, b]
        product = sp.SparseMatrix(3*size, 3*size, entries)
        residual = product*lift-lift*lap
        if residual != sp.zeros(3*size, size):
            raise ArithmeticError('the angular source readout distinguishes the proposed extension')
        energy = (contrast.T*radial*contrast)[0]/2
        full_energy = (lifted_contrast.T*product*lifted_contrast)[0]/(2*size)
        if energy != full_energy:
            raise ArithmeticError('the distinguishing radial readout changed under angular lifting')
        spectrum = sorted(value for value, multiplicity in radial.eigenvals().items()
                          for _ in range(multiplicity))
        probe_time = sp.Rational(1, 100)
        center = 1-energy*probe_time**2/2
        remainder = max(spectrum)**2*probe_time**4/24
        models.append({'name': name, 'radial_edges': [list(edge) for edge in edges],
            'radial_edge_weight': weight,
            'radial_laplacian': [[str(v) for v in row] for row in radial.tolist()],
            'radial_spectrum': [str(v) for v in spectrum],
            'cycle_rank': len(edges)-2,
            'exact_product_intertwiner_residual_nonzeros': residual.nnz(),
            'radial_contrast_wave_second_derivative_at_zero': str(-energy),
            'nonconstant_radial_probe_time': str(probe_time),
            'nonconstant_radial_cosine_interval': [str(center-remainder), str(center+remainder)]})
    if len({row['radial_contrast_wave_second_derivative_at_zero'] for row in models}) != 3:
        raise ArithmeticError('the additional radial measurement failed to distinguish the models')
    probe_intervals = sorted(tuple(map(sp.Rational, row['nonconstant_radial_cosine_interval'])) for row in models)
    if any(left[1] >= right[0] for left, right in zip(probe_intervals, probe_intervals[1:])):
        raise ArithmeticError('the finite radial probe intervals do not distinguish the extensions')
    return {'native_ports': size, 'native_seams': len(graph.seams),
        'source': 'unchanged complete total-one native endpoint-record kernel',
        'product_ports': 3*size, 'radial_readout': '(1,1,1)/sqrt(3)',
        'exact_identity': 'B_R J = J L; J*J=I; B_R=I_3 tensor L+R tensor I',
        'models': models, 'complete_adaptive_angular_histories': _radial_adaptive_history_control(graph, lap, models),
        'retained_functional_calculus': ['heat', 'cosine wave', 'positive shifted resolvent'],
        'native_record_and_incident_mark_law_is_independent_of_radial_choice': True,
        'distinguishing_measurement': 'constant angular mode times radial contrast (1,-1,0)/sqrt(2)',
        'same_angular_clock_scale_used_for_all_comparisons': True,
        'radial_average_is_a_boundary_Dirichlet_to_Neumann_readout': False,
        'additional_radial_contrast_is_available_in_current_native_record_menu': False,
        'radial_models_are_executed_native_repair_extensions': False,
        'existing_Dirichlet_product_or_general_OPH_is_refuted': False}


def execute_radial_records(seams):
    """Complete total-one native branches; final loads are not exported."""
    records = []
    for p in range(PORTS):
        for seam in seams:
            for coin in (0, 1):
                _, (_, current) = native.native(tuple(int(i == p) for i in range(PORTS)), seam, coin)
                records.append({'preparation_port': p, 'seam': list(seam),
                                'coin': coin, 'current': current})
    return records


def _lumped_kernel(kernel, blocks):
    """Exact row-kernel lumping; averages cannot replace equal block sums."""
    if (not blocks or any(not block for block in blocks)
            or sorted(p for block in blocks for p in block) != list(range(kernel.rows))):
        raise ValueError('a partition of the complete kernel required')
    quotient = sp.Matrix([[sum(kernel[block[0], q] for q in target)
                           for target in blocks] for block in blocks])
    for i, block in enumerate(blocks):
        if any(sum(kernel[p, q] for q in target) != quotient[i, j]
               for p in block for j, target in enumerate(blocks)):
            raise ValueError('shell records are not strongly lumpable')
    return quotient


def recover_radial_records(records, observer=0):
    """Recover shell histories from complete labelled native branch records.

    Equal branch weights are the declared uniform-seam/fair-coin protocol,
    not empirical frequencies. Known total-one preparations and authentic
    complete current access are required. No hidden endpoint or graph is an
    input. The resulting radius is graph distance on this same carrier.
    """
    if type(observer) is not int or observer not in range(PORTS):
        raise ValueError('a native observer port required')
    if type(records) not in (list, tuple) or not records:
        raise ValueError('complete native branch records required')
    outcomes, seams = {}, set()
    for row in records:
        if type(row) is not dict or set(row) != ng.RECORD_FIELDS-{'probability'}:
            raise ValueError('only preparation, seam, coin and current records accepted')
        p, seam, coin, current = (row[name] for name in ('preparation_port', 'seam', 'coin', 'current'))
        if (type(p) is not int or p not in range(PORTS)
                or type(seam) not in (list, tuple) or len(seam) != 2
                or any(type(i) is not int or i not in range(PORTS) for i in seam)
                or seam[0] >= seam[1] or type(coin) is not int or coin not in (0, 1)
                or type(current) is not int or current not in (-1, 0, 1)):
            raise ValueError('native preparation, ordered seam, coin and current required')
        a, b = seam
        difference = int(p == a)-int(p == b)
        if current != (difference-(2*coin-1)*abs(difference))//2:
            raise ValueError('current contradicts the declared total-one native preparation')
        key = p, a, b, coin
        if key in outcomes:
            raise ValueError('duplicate native branch record')
        decoded = [int(i == p)+current*(int(i == b)-int(i == a)) for i in range(PORTS)]
        outcomes[key] = decoded.index(1)
        seams.add((a, b))
    if set(outcomes) != {(p, a, b, coin) for p in range(PORTS) for a, b in seams for coin in (0, 1)}:
        raise ValueError('complete preparations and both coins on every seam required')
    denominator = 2*len(seams)
    remote, incident = sp.zeros(PORTS), sp.zeros(PORTS)
    for (p, a, b, _), q in outcomes.items():
        (incident if observer in (a, b) else remote)[p, q] += sp.Rational(1, denominator)
    kernel = remote+incident
    shells, unseen = [(observer,)], set(range(PORTS))-{observer}
    while unseen:
        layer = tuple(sorted(q for q in unseen if any(kernel[p, q] > 0 for p in shells[-1])))
        if not layer:
            raise ValueError('the recovered native carrier is not connected')
        shells.append(layer); unseen.difference_update(layer)
    q0, q1 = (_lumped_kernel(matrix, shells) for matrix in (remote, incident))
    radial = denominator*(sp.eye(len(shells))-q0-q1)
    stopped = (sp.eye(len(shells))-q0).inv()*q1
    encode = lambda matrix: [[str(x) for x in row] for row in matrix.tolist()]
    return {'observer': observer, 'record_count': len(records), 'events_per_preparation': denominator,
        'shells': [list(block) for block in shells], 'shell_sizes': list(map(len, shells)),
        'native_row_kernel': encode(kernel), 'radial_generator': encode(radial),
        'remote_row_kernel': encode(q0), 'incident_row_kernel': encode(q1),
        'first_incident_row_kernel': encode(stopped),
        'all_marked_shell_histories_close_exactly': True,
        'hidden_final_loads_or_supplied_radial_graph_used': False,
        'uniform_seam_fair_coin_law_and_total_one_preparations_declared': True,
        'finite_sample_or_laboratory_noise_bound': False}


def native_radial_record_control():
    """Source shell quotient, marked-clock obstruction and exact arc readout."""
    structure = recover_structure(ng.execute_calibration())
    records = execute_radial_records(tuple(tuple(edge) for edge in structure['seams']))
    recovered = recover_radial_records(records)
    shells = recovered['shells']
    radial = sp.Matrix(recovered['radial_generator'])
    kernel = sp.Matrix(recovered['native_row_kernel'])
    q0, q1 = (sp.Matrix(recovered[name]) for name in ('remote_row_kernel', 'incident_row_kernel'))
    stopped = sp.Matrix(recovered['first_incident_row_kernel'])
    lift = sp.Matrix(PORTS, len(shells), lambda p, r: int(p in shells[r]))
    lap = 60*(sp.eye(PORTS)-kernel)
    mass = sp.diag(*recovered['shell_sizes'])
    if (recovered['shell_sizes'] != [1, 5, 5, 1] or lap*lift != lift*radial
            or mass*radial != radial.T*mass or lift.T*lift != mass
            or q0*sp.ones(4, 1) != sp.ones(4, 1)*sp.Rational(5, 6)):
        raise ArithmeticError('native radial histories lost their exact weighted quotient')
    commutator = stopped*radial-radial*stopped
    if commutator == sp.zeros(4):
        raise ArithmeticError('the stopped clock falsely became a scalar spectral continuation')
    wrong_shells = [list(block) for block in shells]
    wrong_shells[1][0], wrong_shells[2][0] = wrong_shells[2][0], wrong_shells[1][0]
    try:
        _lumped_kernel(kernel, wrong_shells)
    except ValueError:
        pass
    else:
        raise ArithmeticError('an inequitable shell partition falsely closed its histories')
    # Native arc reversal retains zero-current outcomes. Amplitudes are
    # stored multiplied by sqrt(60), so all operations stay rational.
    keys = [(row['preparation_port'], *row['seam'], row['coin']) for row in records]
    index = {key: i for i, key in enumerate(keys)}
    reversal = []
    for row in records:
        p, (a, b), coin, current = (row[name] for name in ('preparation_port', 'seam', 'coin', 'current'))
        q = a+b-p if current else p
        reversal.append(index[q, a, b, 1-coin if current else coin])
    if any(reversal[reversal[i]] != i for i in range(len(records))):
        raise ArithmeticError('native radial arc reversal is not an involution')
    vector = [sp.Integer(p == 0) for p, *_ in keys]
    previous, cosine = sp.eye(PORTS), kernel
    wave_rows = []
    for steps in range(4):
        compressed = sp.Matrix([sum(vector[i] for i, key in enumerate(keys) if key[0] == p)/60
                                for p in range(PORTS)])
        expected = sp.eye(PORTS)[:, 0] if steps == 0 else cosine[:, 0]
        if compressed != expected or any(compressed[p] != compressed[block[0]] for block in shells for p in block):
            raise ArithmeticError('native arc wave failed its radial polynomial readout')
        wave_rows.append({'incident_triggers': steps,
            'shell_compressed_amplitudes': [str(compressed[block[0]]) for block in shells],
            'Ramsey_plus_probability': str((1+compressed[0])/2),
            'Ramsey_minus_probability': str((1-compressed[0])/2)})
        means = [sum(vector[i] for i, key in enumerate(keys) if key[0] == p)/60 for p in range(PORTS)]
        reflected = [2*means[key[0]]-value for key, value in zip(keys, vector, strict=True)]
        vector = [reflected[reversal[i]] for i in range(len(keys))]
        if steps:
            previous, cosine = cosine, 2*kernel*cosine-previous
    return {**recovered, 'exact_native_arc_count': len(records),
        'weighted_self_adjoint_radial_generator': True,
        'radial_generator_characteristic_polynomial': str(radial.charpoly().as_expr()),
        'stopped_clock_commutator_entry_0_1': str(commutator[0, 1]),
        'exchanged_shell_ports_control_rejected': True,
        'first_incident_timeout_after_N_attempts': '(5/6)^N',
        'active_wave_complete_Ramsey_outputs': wave_rows,
        'active_wave_protocol': 'one protected native arc walk per incident trigger; timing source pauses during service',
        'passive_stopped_density_equals_active_wave': False,
        'source_pulses_for_this_radial_input_reexecuted': False,
        'radial_shells_are_an_independent_product_coordinate': False,
        'product_radial_factor_or_physical_length_selected': False}


def run():
    uniform_records = ng.execute_calibration()
    uniform = ng.recover(uniform_records)
    structural = recover_structure(uniform_records)
    seams = tuple(tuple(seam) for seam in structural['seams'])
    biased_records = execute_biased_calibration(seams)
    biased = recover_structure(biased_records)
    if not biased['rank_three_structural_gram_certified']:
        raise AssertionError('biased calibration failed exact structural recovery')
    if biased['gram_qsqrt5'] != uniform['gram_qsqrt5']:
        raise AssertionError('seam exposure changed the structural Gram')
    if biased['laplacian_exact'] != uniform['laplacian_exact']:
        raise AssertionError('seam exposure changed the structural response')
    shot_records = execute_single_shot_calibration(seams)
    shots = recover_from_shots(shot_records)
    if shots['gram_qsqrt5'] != uniform['gram_qsqrt5']:
        raise AssertionError('finite native shots failed exact structural recovery')
    laplacian = [[F(x) for x in row] for row in structural['laplacian_exact']]
    missing = recover_structure(ng.execute_calibration({seams[0]: F(1)}))
    if missing['rank_three_structural_gram_certified']:
        raise AssertionError('one exposed seam cannot inherit the carrier certificate')
    doubled_weights = {seam: F(1+int(seam == seams[0]), 31) for seam in seams}
    doubled_records = ng.execute_calibration(doubled_weights)
    actual_biased = ng.recover(doubled_records)
    for i in range(PORTS):
        for j in range(PORTS):
            edge_i = int(i == seams[0][0])-int(i == seams[0][1])
            edge_j = int(j == seams[0][0])-int(j == seams[0][1])
            expected = F(i == j)-(laplacian[i][j]+edge_i*edge_j)/62
            if F(actual_biased['transition_matrix_exact'][i][j]) != expected:
                raise AssertionError('weighted source response does not match rank-one perturbation')
    native_increment = native_increment_control()
    # Lazy imports avoid a cycle through the carrier calibration above. The
    # live receipt binds every dependency of the executed ancestry readout.
    from ..physical_space_time import uniform_green, local_continuum
    source_paths = sorted(set((*SOURCE_PATHS, *uniform_green.SOURCE_PATHS,
                              *local_continuum.SOURCE_PATHS)))
    return {
        'schema': 'oph.relative_observer_geometry.intrinsic.v1',
        'sources': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in source_paths},
        'single_carrier_recovery': biased,
        'same_structural_gram_as_uniform_protocol_exact': True,
        'finite_shot_recovery': {
            'native_branch_calls': len(shot_records),
            'shots_per_preparation_and_seam': 1,
            'probability_law_oracle_used': shots['probability_law_oracle_used'],
            'record_fields': sorted(ng.RECORD_FIELDS-{'probability'}),
            'same_structural_gram_exact': True,
            'conditional_current_variance_zero_under_native_model': True,
            'laboratory_noise_or_preparation_error_bounded': False},
        'biased_protocol': {'preparation_count': PORTS, 'native_branch_calls': len(biased_records),
            'preparation_dependent_positive_seam_exposure': True,
            'even_preparations_have_deterministic_coin_zero': True,
            'odd_preparations_have_coin_probabilities': ['1/3', '2/3'],
            'complete_exact_observation_laws_not_finite_shot_estimates': True},
        'actual_scheduler_slow_band_control': scheduling_control(laplacian, seams[0]),
        'full_response_metric': inverse_response_metric(laplacian),
        'product_wave_metric': product_wave_metric_control(),
        'native_increment_metric': native_increment,
        'native_optical_front': native_optical_front_control(),
        'radial_selection_control': radial_nonidentifiability_control(),
        'native_radial_record_quotient': native_radial_record_control(),
        'incomplete_exposure_control': {'exposed_seams': len(missing['seams']),
            'rank_three_structural_gram_certified': missing['rank_three_structural_gram_certified']},
        'scope': {'federation_needed_for_finite_carrier_result': False,
            'uniform_schedule_needed_for_structural_gram_recovery': False,
            'coin_fairness_needed_for_even_preparation_readout': False,
            'arbitrary_positive_schedule_preserves_rank_three_native_slow_band': False,
            'actual_native_slow_band_requires_its_own_schedule_certificate': True,
            'complete_ensemble_probability_law_needed_for_even_preparation_structure': False,
            'known_preparation_faithful_records_and_all_seam_coverage_still_required': True,
            'observed_seam_labels_disclose_interacting_endpoint_incidence': True,
            'carrier_architecture_derived_from_self_reference': False,
            'slow_band_definition_selects_one_metric_not_all_observer_readouts': True,
            'rank_three_proves_three_dimensional_physical_bulk': False,
            'specific_federation_architecture_required_for_continuum_claim': False,
            'extension_or_limit_with_compatible_observable_laws_required_for_continuum_claim': True,
            'physical_length_time_or_lab_calibration_established': False},
        'all_pass': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Single-carrier intrinsic geometry receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
