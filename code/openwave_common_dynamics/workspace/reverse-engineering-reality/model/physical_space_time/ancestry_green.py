"""Certified Green/clock observations on the actual source ancestry tower.

Finite cross-level comparisons use one fixed twenty-face readout dictionary.
They do not extrapolate a continuum limit from a finite sequence.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as F
from functools import lru_cache
import math
from types import MappingProxyType

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from . import source_green as elimination

SOURCE_PATHS = tuple(sorted(set((*elimination.SOURCE_PATHS,
    'model/physical_space_time/ancestry_green.py'))))
DIMENSION = 20
DYADIC_BITS = 40
RETAINED_OPERATOR_ERROR = F(1, 1000000)
# compare_receipt permits 2e-10 per floating entry; for a 20x20 matrix its
# operator allowance is at most 20 times that. Include it in the certificate.
REPLAY_MATRIX_ALLOWANCE = F(1, 250000000)


def _ceil_sqrt(integer):
    value = math.isqrt(integer)
    return value+int(value*value != integer)


def problem(level):
    """Check the actual source rows, local diameter, ancestry and graph bound."""
    if type(level) is not int or level not in range(4):
        raise ValueError('an original source level zero through three is required')
    builder = elimination.source._source()
    raw = builder.build_federation(level, 'port_pair')
    parents = tuple(map(int, builder.build_geodesic_icosahedral_tower(level).cell_ancestor_ids(
        fine_level=level, coarse_level=0)))
    scale = 4**level
    if (raw.carriers != 20*scale or raw.ports != 240*scale
            or parents != tuple(cell//scale for cell in range(raw.carriers))):
        raise ArithmeticError('the fixed twenty-face ancestry dictionary changed')
    local, _ = elimination.passage.recovered_carrier()
    local_set = set(local)
    intra = [set() for _ in range(raw.carriers)]
    cell_neighbors = [set() for _ in range(raw.carriers)]
    endpoints, edges = [], []
    for left, right in zip(raw.seam_a, raw.seam_b, strict=True):
        a, b = sorted((int(left), int(right)))
        edges.append((a, b))
        if a//12 == b//12:
            intra[a//12].add((a % 12, b % 12))
        else:
            endpoints.extend((a, b))
            cell_neighbors[a//12].add(b//12); cell_neighbors[b//12].add(a//12)
    if (any(found != local_set for found in intra) or len(set(edges)) != len(edges)
            or len(endpoints) != len(set(endpoints))
            or any(len(row) != 3 for row in cell_neighbors)):
        raise ArithmeticError('native carrier or matching hypotheses changed')
    # Independently rebuild the full integer Laplacian from actual source seams.
    row, col, values = [], [], []
    for a, b in edges:
        row.extend((a, b, a, b)); col.extend((a, b, b, a)); values.extend((1, 1, -1, -1))
    laplacian = sparse.coo_matrix((values, (row, col)),
        shape=(raw.ports, raw.ports), dtype=np.int64).tocsc()
    if (laplacian-raw.laplacian).nnz:
        raise ArithmeticError('source seam and operator declarations differ')
    local_neighbors = [set() for _ in range(12)]
    for a, b in local:
        local_neighbors[a].add(b); local_neighbors[b].add(a)
    def distances(neighbors, start):
        values = {start: 0}; pending = deque((start,))
        while pending:
            old = pending.popleft()
            for new in neighbors[old]:
                if new not in values:
                    values[new] = values[old]+1; pending.append(new)
        return values
    local_distances = tuple(distances(local_neighbors, start) for start in range(12))
    if any(len(row) != 12 or max(row.values()) != 3 for row in local_distances):
        raise ArithmeticError('the internal source diameter is not three')
    rooted = distances(cell_neighbors, 0)
    if len(rooted) != raw.carriers:
        raise ArithmeticError('the cell graph is disconnected')
    eccentricity = max(rooted.values())
    diameter_bound = 8*eccentricity+3
    inverse_bound = F(raw.ports*diameter_bound, 4)
    # Integer scores encode mean-free normalized base-face observations.
    score = np.full((raw.ports, DIMENSION), -1, dtype=np.int64)
    score[np.arange(raw.ports), np.repeat(parents, 12)] += DIMENSION
    if np.any(score.sum(axis=0)) or np.any(score.sum(axis=1)):
        raise ArithmeticError('ancestry observations lost their constant-free normalization')
    incident = tuple(edge for edge in edges if edge[0]//12 == 0 or edge[1]//12 == 0)
    if len(incident) != 33:
        raise ArithmeticError('the actual root incident clock changed')
    return raw, laplacian, score, incident, eccentricity, inverse_bound


def certify_solution(laplacian, score, candidate, inverse_bound, scale, *, bits=DYADIC_BITS):
    """Exact residual with the independently proved inverse bound from problem()."""
    if (type(bits) is not int or not 16 <= bits <= 44 or type(scale) is not int or scale < 1
            or type(inverse_bound) is not F
            or type(inverse_bound.numerator) is not int
            or type(inverse_bound.denominator) is not int
            or inverse_bound.numerator <= 0 or inverse_bound.denominator <= 0
            or math.gcd(inverse_bound.numerator, inverse_bound.denominator) != 1
            or candidate.shape != score.shape or not np.isfinite(candidate).all()):
        raise ValueError('finite candidate and declared positive exact bounds required')
    if (not sparse.issparse(laplacian) or laplacian.dtype.kind not in 'iu'
            or not isinstance(score, np.ndarray) or score.dtype.kind not in 'iu'
            or score.shape != (240*scale, 20) or laplacian.shape != (240*scale, 240*scale)
            or not np.isin(score, (-1, 19)).all() or np.any(score.sum(axis=0))
            or np.any(score.sum(axis=1)) or np.any(laplacian.data < -1)
            or np.any(laplacian.data > 6)
            or np.max(abs(laplacian).sum(axis=1)) > 12):
        raise ValueError('integer native operator and normalized integer ancestry scores required')
    denominator = 2**bits
    if np.max(np.abs(candidate))*denominator >= 2**56:
        raise ValueError('candidate exceeds the exact integer arithmetic allowance')
    encoded = np.rint(candidate*denominator).astype(np.int64)
    # Preserve the exact zero sum over the twenty observation channels.
    encoded[:, -1] = -encoded[:, :-1].sum(axis=1)
    maximum = max(abs(int(encoded.min())), abs(int(encoded.max())))
    if 12*maximum+19*denominator >= 2**63:
        raise ValueError('integer residual multiplication could overflow')
    residual = laplacian@encoded-denominator*score
    squared = sum(int(value)**2 for value in residual.flat)
    residual_norm_bound = F(_ceil_sqrt(squared), denominator)
    score_norm_bound = _ceil_sqrt(4800*scale)
    error = score_norm_bound*inverse_bound*residual_norm_bound/(4800*scale**2)
    return encoded, error, squared


def _retained_rows(matrix, exact_error):
    result = [[float(value) for value in row] for row in matrix]
    rounding = max(sum(abs(F(result[i][j])-matrix[i][j]) for j in range(20)) for i in range(20))
    if exact_error+rounding+REPLAY_MATRIX_ALLOWANCE > RETAINED_OPERATOR_ERROR:
        raise ArithmeticError('exact residual, conversion and replay exceed the retained covariance enclosure')
    return result


@lru_cache(maxsize=4, typed=True)
def covariance(level):
    """Return exact rational approximants and certified operator error bounds."""
    raw, laplacian, score, incident, eccentricity, inverse_bound = problem(level)
    scale = 4**level
    candidate = np.zeros(score.shape)
    candidate[1:] = splu(laplacian[1:, 1:].astype(float)).solve(score[1:].astype(float))
    candidate -= candidate.mean(axis=0)
    encoded, error, residual_squared = certify_solution(
        laplacian, score, candidate, inverse_bound, scale)
    total = encoded.sum(axis=0, dtype=object)
    groups = encoded.reshape(20, 12*scale, 20).sum(axis=1, dtype=object)
    pairing = 20*groups-total
    denominator = 4800*scale**2*2**DYADIC_BITS
    green = tuple(tuple(F(int(pairing[i, j]+pairing[j, i]), 2*denominator)
                        for j in range(20)) for i in range(20))
    local_pairing = np.zeros((20, 20), dtype=object)
    for a, b in incident:
        local_pairing += np.outer((score[a]-score[b]).astype(object),
                                 (encoded[a]-encoded[b]).astype(object))
    projector = tuple(tuple(F(int(i == j))-F(1, 20) for j in range(20)) for i in range(20))
    attempt = tuple(tuple(green[i][j]-projector[i][j]/(4*raw.seams*scale)
                          for j in range(20)) for i in range(20))
    clock = tuple(tuple(green[i][j]+projector[i][j]/(132*scale)
        -F(int(local_pairing[i, j]+local_pairing[j, i]), 132*denominator)
        for j in range(20)) for i in range(20))
    if any(sum(row) for matrix in (green, attempt, clock) for row in matrix):
        raise ArithmeticError('exact candidate covariances lost the constant mode')
    norm_bound = max(sum(map(abs, row)) for row in green)+error
    # ||W^T L_i W|| <= 1/(2s): only three cross seams can change cell values.
    root = _ceil_sqrt((norm_bound/2).numerator*(norm_bound/2).denominator)
    sqrt_bound = F(root, (norm_bound/2).denominator)
    clock_bound = (1+2*sqrt_bound)/(132*scale)
    return MappingProxyType({'level': level, 'ports': raw.ports, 'carriers': raw.carriers,
        'ancestry_scale': scale, 'cell_root_eccentricity': eccentricity,
        'full_port_diameter_upper': 8*eccentricity+3,
        'native_inverse_norm_upper': inverse_bound,
        'dyadic_bits': DYADIC_BITS, 'exact_residual_squared_numerator': residual_squared,
        'green': green, 'attempt': attempt, 'incident': clock,
        'green_error': error, 'attempt_error': error, 'incident_error': F(13, 11)*error,
        'incident_clock_correction_bound': clock_bound,
        'attempt_clock_correction_bound': F(1, 4*raw.seams*scale)})


def _comparison(first, second, key):
    left, right = first[key], second[key]
    delta = tuple(tuple(right[i][j]-left[i][j] for j in range(20)) for i in range(20))
    error = 2*RETAINED_OPERATOR_ERROR
    lower = max(abs(x) for row in delta for x in row)
    upper = max(sum(map(abs, row)) for row in delta)
    numeric = np.asarray(delta, dtype=float)
    values, vectors = np.linalg.eigh(numeric)
    vector = np.rint(vectors[:, int(np.argmax(np.abs(values)))]*2**24).astype(np.int64)
    norm = sum(int(x)**2 for x in vector)
    rayleigh = sum(delta[i][j]*int(vector[i])*int(vector[j])
                   for i in range(20) for j in range(20))/norm
    lower = max(lower, abs(rayleigh))
    # Outward conversion and 1e-8 extra margin cover the existing 2e-10
    # scalar replay tolerance without weakening that shared verifier.
    margin = F(1, 100000000)
    lower_bound = max(F(0), lower-error-margin)
    return {'field': key, 'source_levels': [first['level'], second['level']],
        'operator_difference_lower': max(0., math.nextafter(float(lower_bound), -math.inf)),
        'operator_difference_upper': math.nextafter(float(upper+error+margin), math.inf),
        'operator_difference_diagnostic': float(np.linalg.norm(numeric, 2)),
        'finite_comparison_does_not_bound_the_unexecuted_refinement_tail': True}


def run():
    values = tuple(covariance(level) for level in range(4))
    rows = []
    for value in values:
        row = {key: value[key] for key in ('level', 'ports', 'carriers', 'ancestry_scale',
            'cell_root_eccentricity', 'full_port_diameter_upper', 'dyadic_bits')}
        if max(sum(map(abs, row)) for row in value['green'])+value['green_error'] > 36:
            raise ArithmeticError('the executed Green norm exceeds its declared finite bound')
        clock_bound = F(19, 264*value['ancestry_scale'])
        if value['incident_clock_correction_bound'] > clock_bound:
            raise ArithmeticError('the source clock correction exceeds its outward rational enclosure')
        row.update({'native_inverse_norm_upper': str(value['native_inverse_norm_upper']),
            'exact_dyadic_residual_checked': True,
            'retained_covariance_operator_error': str(RETAINED_OPERATOR_ERROR),
            'floating_receipt_replay_operator_allowance_included': str(REPLAY_MATRIX_ALLOWANCE),
            'certified_Green_operator_norm_upper': '36',
            'incident_clock_correction_bound': str(clock_bound),
            'attempt_clock_correction_bound': str(value['attempt_clock_correction_bound'])})
        row.update({key+'_covariance': _retained_rows(value[key], value[key+'_error'])
                    for key in ('green', 'attempt', 'incident')})
        rows.append(row)
    return {'contract': {
        'readout': 'twenty fixed base-face ancestry counts, centered and normalized by sqrt(12*4^level)',
        'spatial_operator_scale': 'declared 4^level combinatorial normalization',
        'stationary_total': 'ports/3', 'full_initial_and_seam_current_disclosure_required': True,
        'native_attempt_and_root_cell_incident_clocks_compared': True,
        'time_average_limit_precedes_spatial_comparison': True,
        'physical_radial_coordinate_or_continuum_metric_supplied': False},
        'certified_source_levels': rows,
        'consecutive_same_readout_comparisons': [_comparison(a, b, key)
            for a, b in zip(values, values[1:]) for key in ('green', 'attempt', 'incident')],
        'acceptance': {'same_actual_source_ancestry_Green_and_clock_observations_certified': True,
            'bounded_Green_sequence_implies_vanishing_incident_clock_correction': True,
            'all_level_Green_Cauchy_tail_or_limit_proved': False,
            'interacting_QFT_or_Lorentzian_spacetime_identified': False}}
