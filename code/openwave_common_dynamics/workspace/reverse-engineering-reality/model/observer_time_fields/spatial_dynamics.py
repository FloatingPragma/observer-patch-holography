"""Native temporal correlations and the actual spatial Green operator.

The parameter counts native attempts; geometric discount is processing of
their correlations. Neither the logarithm nor its unitary continuation is
identified with physical time. All retained arithmetic is exact.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from types import MappingProxyType

import sympy as sp

from model.observer import native
from model.observer import coherent_continuation as coherent
from model.physical_space_time import source_green
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*source_green.SOURCE_PATHS, *coherent.SOURCE_PATHS,
    'model/__init__.py', 'model/observer_time_fields/__init__.py',
    'model/observer_time_fields/spatial_dynamics.py',
    'model/observer/__init__.py', 'model/observer/scalar_response/__init__.py',
    'model/observer/scalar_response/experiment.py', 'model/physical_space_time/__init__.py',
    'model/common_source/__init__.py', 'model/common_source/experiment.py',
    'model/common_source/gauge_matter.py', 'model/common_source/native_join.py',
    'model/common_source/records.py', 'model/common_source/scalar_coupling.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/source_scalar/__init__.py', 'model/observer/source_scalar/bounds.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/source_scalar/cutoff.py',
    'model/observer/source_scalar/experiment.py', 'model/observer/source_scalar/field.py',
    'model/observer_reconstruction/__init__.py', 'model/relative_observer_geometry/__init__.py',
    'tools/strict_json.py'))))


def _level(value):
    if type(value) is not int or not 0 <= value <= 100:
        raise ValueError('an original subdivision level from zero through one hundred required')
    return value


def _positive(value):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError('a positive exact rational required')
    return F(value)


def _fraction(value):
    return F(int(value.p), int(value.q)) if isinstance(value, sp.Rational) else F(value)


def _rows(matrix):
    return tuple(tuple(str(x) for x in row) for row in matrix.tolist())


def _digest(matrix):
    return hashlib.sha256(json.dumps(_rows(matrix), separators=(',', ':')).encode()).hexdigest()


@lru_cache(maxsize=1)
def native_density():
    """All3960 endpoint/coin branches, independent of a density closure."""
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    seams, laplacian = source_green.passage.recovered_carrier()
    counts = sp.zeros(66)
    branches = 0
    for column, before in enumerate(states):
        for edge in seams:
            for coin in (0, 1):
                after, (_, current) = native.native(before, edge, coin)
                decoded = tuple(before[j]+current*(int(j == edge[1])-int(j == edge[0]))
                                for j in range(12))
                if after != decoded:
                    raise ArithmeticError('the current tape lost field reconstruction')
                counts[index[after], column] += 1
                branches += 1
    z = sp.Matrix(12, 66, lambda i, j: states[j][i]-sp.Rational(1, 6))
    p0 = sp.eye(12)-sp.ones(12)/12
    chi = sp.Rational(5, 33)
    if (counts != counts.T or sp.ones(1, 66)*counts != 60*sp.ones(1, 66)
            or z*counts != (60*sp.eye(12)-laplacian)*z or z*z.T/66 != chi*p0):
        raise ArithmeticError('native exclusion, density closure or canonical covariance changed')
    # Exact source gap; this conservative lambda=2 is sufficient for all routing.
    _, pivots = (laplacian-2*p0)[:11, :11].LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(11)):
        raise ArithmeticError('the internal native gap no longer exceeds two')
    return MappingProxyType({'counts': sp.ImmutableMatrix(counts), 'field': sp.ImmutableMatrix(z),
        'laplacian': laplacian, 'center': sp.ImmutableMatrix(p0), 'chi': chi,
        'branches': branches})


def covariance(attempts):
    """Exact two-time field matrix and independent66-state native oracle."""
    if type(attempts) is not int or not 0 <= attempts <= 16:
        raise ValueError('zero through sixteen native attempts required')
    data = native_density()
    reduced = (sp.eye(12)-data['laplacian']/60)**attempts
    result = data['chi']*reduced*data['center']
    oracle = data['field']*data['counts']**attempts*data['field'].T/(66*60**attempts)
    if result != oracle:
        raise ArithmeticError('native complete configuration correlation disagrees')
    return sp.ImmutableMatrix(result)


def native_discount(z=1):
    """Exact geometric-discount resolvent in the original66-state sector."""
    z = sp.Rational(_positive(z))
    data = native_density()
    # (1/60) sum (60/(60+z))^(n+1) M^n = (60I+zI-C)^-1.
    oracle = data['field']*((60+z)*sp.eye(66)-data['counts']).inv(method='DM')*data['field'].T/66
    expected = data['chi']*(data['laplacian']+z*sp.eye(12)).inv()*data['center']
    if oracle != expected:
        raise ArithmeticError('native discounted correlation differs from the repair resolvent')
    return sp.ImmutableMatrix(expected)


def bounds(level, z=1, band=1, terms=1):
    """All-level exact analytic bounds; no screen allocation or fitted scale."""
    level = _level(level)
    z, band = _positive(z), _positive(band)
    if type(terms) is not int or not 1 <= terms <= 10000:
        raise ValueError('one through ten thousand correlation terms required')
    scale = 4**level
    seams = 630*scale
    normalization = 2*seams*scale
    incident_normalization = 66*scale
    if band >= incident_normalization:
        raise ValueError('the scaled low band must lie below the incident normalization')
    discount = F(normalization, normalization+z)
    return {'level': level, 'spatial_scale': scale, 'seams': seams,
        'native_attempt_normalization': normalization,
        'root_incident_normalization': incident_normalization,
        'discount': discount,
        'discounted_correlation_tail': discount**terms/z,
        'dynamic_Schur_resolvent_error': F(1, 8*scale),
        'scaled_attempt_log_error_on_band': band**2/(2*(normalization-band)),
        'scaled_incident_ordered_log_eigenvalue_error_on_band': band**2/(2*(incident_normalization-band))}


@lru_cache(maxsize=1)
def _base_blocks():
    graph, laplacian, raw_r, integer_f, denominator = source_green.matrices(0)
    local = native_density()['laplacian']
    d = sp.zeros(graph.ports, len(graph.inter_seams))
    for j, (a, b) in enumerate(graph.inter_seams):
        d[a, j], d[b, j] = 1, -1
    v = sp.Matrix(graph.ports, graph.cells, lambda i, j: int(i//12 == j))
    r = sp.Matrix(raw_r.tolist())
    f = sp.Matrix(integer_f.tolist())/denominator
    full = sp.ImmutableSparseMatrix(laplacian.shape[0], laplacian.shape[1],
        {(int(i), int(j)): int(laplacian[i, j]) for i, j in zip(*laplacian.nonzero(), strict=True)})
    if (d.T*d != 2*sp.eye(len(graph.inter_seams)) or d.T*v != r
            or full != sp.diag(*([local]*graph.cells))+d*d.T):
        raise ArithmeticError('actual native source block or matching identity changed')
    schur = r.T*f.inv(method='DM')*r/12
    if schur*sp.ones(graph.cells, 1) != sp.zeros(graph.cells, 1):
        raise ArithmeticError('the native Schur constant mode changed')
    return tuple(map(sp.ImmutableMatrix, (local, d, v, r, full, schur)))


@lru_cache(maxsize=4)
def _base_resolvent(z):
    local, d, v, r, full, schur = _base_blocks()
    z = sp.Rational(z)
    local_inverse = (local+z*sp.eye(12)).inv(method='DM')
    inverse_d = sp.Matrix.vstack(*(local_inverse*d[12*j:12*j+12, :] for j in range(20)))
    middle = (sp.eye(30)+d.T*inverse_d).inv(method='DM')
    solution = v/z-inverse_d*middle*r/z
    # An exact240port Dirichlet residual independently validates the reduction.
    if (full+z*sp.eye(240))*solution != v:
        raise ArithmeticError('the full native discounted solve has a nonzero exact residual')
    compressed = v.T*solution/12
    effective = (schur+z*sp.eye(20)).inv(method='DM')
    difference = effective-compressed
    if difference != difference.T or difference*sp.ones(20, 1) != sp.zeros(20, 1):
        raise ArithmeticError('dynamic elimination lost symmetry or constant cancellation')
    # A row bound is a rigorous (possibly loose) numerical witness in addition
    # to the all-level operator proof; its comparison is made over rationals.
    row_bound = max(sum(abs(x) for x in difference.row(i)) for i in range(20))
    if row_bound >= sp.Rational(1, 8):
        raise ArithmeticError('the retained base resolvent exceeds its theorem allowance')
    return MappingProxyType({'compressed': sp.ImmutableMatrix(compressed),
        'effective': sp.ImmutableMatrix(effective), 'difference': sp.ImmutableMatrix(difference),
        'row_bound': _fraction(row_bound)})


def base_resolvent(z=1):
    return _base_resolvent(_positive(z))


def memory_control():
    """Two actual attempts expose the nonclosure of cell-average fields."""
    _, _, v, _, full, _ = _base_blocks()
    compressed_l = v.T*full*v/12
    defect = v.T*full*full*v/12-compressed_l**2
    # This is C C*, with a measurable positive variance on a cell contrast.
    contrast = sp.Matrix([1, -1]+[0]*18)
    witness = (contrast.T*defect*contrast)[0]/1260**2
    if witness <= 0 or defect*sp.ones(20, 1) != sp.zeros(20, 1):
        raise ArithmeticError('the source cell-memory obstruction disappeared')
    return {'two_attempt_compressed_transfer_defect': _fraction(witness),
        'contrast': 'unnormalized cell0 minus cell1',
        'cell_average_native_process_is_closed_Markov_chain': False,
        'defect_trace_before_attempt_normalization': _fraction(sp.trace(defect))}


def incident_spectrum():
    """Exact generalized eigenvalue inclusion on native density observables."""
    data = native_density()
    seams, laplacian = source_green.passage.recovered_carrier()
    incident = tuple(edge for edge in seams if 0 in edge)
    li = sp.zeros(12)
    for a, b in incident:
        li[a, a] += 1; li[b, b] += 1; li[a, b] -= 1; li[b, a] -= 1
    remote = laplacian-li
    d = 2*len(incident)*sp.eye(12)+remote
    b = 2*len(incident)*sp.eye(12)-li
    # Generalized min-max: L/(a+L) <=_ordered D^-1/2 L D^-1/2 <=_ordered L/a.
    if d-b != laplacian:
        raise ArithmeticError('incident density transfer no longer has the native numerator')
    return MappingProxyType({'remote': sp.ImmutableMatrix(remote), 'incident': sp.ImmutableMatrix(li),
        'D': sp.ImmutableMatrix(d), 'B': sp.ImmutableMatrix(b),
        'kernel': sp.ImmutableMatrix(b*d.inv()), 'incident_seams': len(incident)})


def _log_interval(x, degree=48):
    """Outward rational enclosure of -log(1-x), 0<x<1."""
    lo, hi = map(F, x)
    if not 0 < lo <= hi < 1:
        raise ValueError('a strict positive subunit interval required')
    lower = sum((lo**j/j for j in range(1, degree+1)), F())
    upper = sum((hi**j/j for j in range(1, degree+1)), F())
    upper += hi**(degree+1)/((degree+1)*(1-hi))
    return lower, upper


def _outward(interval, scale=10**12):
    lo, hi = interval
    return (F(lo.numerator*scale//lo.denominator, scale),
            F(-((-hi.numerator*scale)//hi.denominator), scale))


def dispersion():
    """Three native bands reject a single mass/speed relativistic fit.

    Squared frequencies are normalized by the first band, so an arbitrary
    common positive clock-unit change cannot remove the residual.
    """
    local = native_density()['laplacian']
    identity = sp.eye(12)
    if local*(local-6*identity)*(local**2-10*local+20*identity) != sp.zeros(12):
        raise ArithmeticError('the native three-band Laplacian polynomial changed')
    if tuple(local.charpoly().all_coeffs()) != tuple(sp.Poly(
            sp.Symbol('lambda')*(sp.Symbol('lambda')-6)**5*
            (sp.Symbol('lambda')**2-10*sp.Symbol('lambda')+20)**3,
            sp.Symbol('lambda')).all_coeffs()):
        raise ArithmeticError('the native band multiplicities changed')
    grid = 1 << 80
    root = math.isqrt(5*grid*grid)
    a, b = F(root, grid), F(root+1, grid)
    eigen = ((5-b, 5-a), (F(6), F(6)), (5+a, 5+b))
    log = tuple(_log_interval((lo/60, hi/60)) for lo, hi in eigen)
    square = tuple((lo*lo, hi*hi) for lo, hi in log)
    normalized = tuple((lo/square[0][1], hi/square[0][0]) for lo, hi in square)
    k = ((eigen[2][0]-eigen[0][1])/(eigen[1][1]-eigen[0][0]),
         (eigen[2][1]-eigen[0][0])/(eigen[1][0]-eigen[0][1]))
    # y3-1-k(y2-1): every exact affine squared-frequency law makes this zero.
    defect = (normalized[2][0]-1-k[1]*(normalized[1][1]-1),
              normalized[2][1]-1-k[0]*(normalized[1][0]-1))
    error = defect[0]/(2*k[1])
    if error <= F(1, 10):
        raise ArithmeticError('the dimensionless multi-band obstruction lost its margin')
    normalized = ((F(1), F(1)), *normalized[1:])
    return {'native_laplacian_eigenvalues': ('5-sqrt(5)', '6', '5+sqrt(5)'),
        'multiplicities': (3, 5, 3), 'attempt_transfer': '1-lambda/60',
        'continued_frequency': '-log(1-lambda/60)',
        'squared_frequency_over_first_band': tuple(_outward(x) for x in normalized),
        'affine_fit_defect_interval': _outward(defect),
        'minimum_maximum_normalized_squared_frequency_error': _outward((error, error))[0],
        'common_clock_unit_rescaling_can_make_all_three_relativistic': False,
        'excludes_other_source_fields_interactions_or_physical_dictionaries': False}


def wave_bounds(level, band=1, parameter=1):
    """Coherent-walk low-band error, at t=n/sqrt(M*4^level).

    This is a different admitted operation from committed native waiting or
    exp(-it H). The relative bound needs no irrational numerical constant.
    """
    level = _level(level)
    band, parameter = _positive(band), _positive(parameter)
    scale = 4**level
    normalization = 630*scale**2
    if band >= 4*normalization:
        raise ValueError('the wave band must lie below four times its normalization')
    relative = band/(6*(4*normalization-band))
    # ceil(sqrt(band)) bounds the amplitude allowance over exact rational inputs.
    root = math.isqrt(band.numerator//band.denominator)
    root += int(F(root*root) < band)
    return {'level': level, 'squared_coherent_clock_normalization': normalization,
        'band_upper': band, 'parameter_upper': parameter,
        'relative_wave_frequency_error': relative,
        'cosine_propagator_error_on_band': min(F(2), parameter*root*relative),
        'squared_frequency_absolute_error_on_band': band*(2*relative+relative**2)}


def coherent_field():
    """Actual native arc execution against a12field Chebyshev recurrence."""
    data = native_density()
    field = tuple(6*state[0]-1 for state in coherent.STATES)
    arc = coherent.embed(field)
    coefficient = sp.Matrix([6]+[0]*11)
    previous = None
    kernel = sp.eye(12)-data['laplacian']/60
    # G[x,p]=x_p-1/6; coefficient6 at port0 reproduces the exact input.
    loads = data['field'].T
    if tuple(loads*coefficient) != field:
        raise ArithmeticError('coherent source sector ordering changed')
    rows = []
    for step in range(9):
        actual = coherent.compress(arc)
        expected = tuple(map(_fraction, loads*coefficient))
        if actual != expected or coherent.norm_squared(arc) != 330:
            raise ArithmeticError('native coherent field recurrence or norm failed')
        pairing = sum((F(x)*y for x, y in zip(field, actual, strict=True)), F())
        rows.append({'coherent_steps': step, 'normalized_return_amplitude': pairing/330,
            'Ramsey_plus': (1+pairing/330)/2, 'Ramsey_minus': (1-pairing/330)/2})
        following = kernel*coefficient
        if previous is not None:
            following = 2*following-previous
        previous, coefficient = coefficient, following
        arc = coherent.advance(arc)
    # Opposite reflection has the same first coherent event but W_plus^2=I.
    prepared = coherent.embed(field)
    if coherent.advance(prepared, 1) != coherent.advance(prepared, -1):
        raise ArithmeticError('first native coherent-event equality failed')
    plus_twice = coherent.advance(coherent.advance(prepared, 1), 1)
    if plus_twice != prepared or rows[2]['Ramsey_plus'] == 1:
        raise ArithmeticError('the admitted continuation-phase control disappeared')
    # A fixed phase on an entire committed record branch changes no CP map,
    # while changing the coherent full-transfer spectrum.
    plain, twisted = coherent.kernel_counts(), coherent.kernel_counts(True)
    moment_difference = F(coherent._trace_cube(twisted)-coherent._trace_cube(plain), 60**3)
    if moment_difference != -F(1, 900):
        raise ArithmeticError('same-instrument cross-record phase control changed')
    return {'native_arc_dimension': coherent.DIMENSION,
        'input': '6*n0-1, divided by sqrt(330) for unit amplitude',
        'all_return_and_complementary_Ramsey_outcomes': rows,
        'same_first_event_opposite_reflection_two_step_Ramsey_plus': F(1),
        'same_committed_instrument_coherent_trace_cube_difference': moment_difference,
        'reflection_and_coherent_uncommitted_record_services_required': True,
        'coherent_clock_parameter': 'n/sqrt(M*4^level)',
        'field_dispersion': 'sqrt(M*4^level)*arccos(1-lambda/(2M))',
        'density_log_relation': 'theta=arccos(exp(-H_density))',
        'active_three_class_modular_probe_is_this_density_subspace': False,
        'native_source_pulse_implementation_of_growing_walk': False}


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, (tuple, list)):
        return [_json(x) for x in value]
    if type(value) is dict:
        return {key: _json(item) for key, item in value.items()}
    if value is None or type(value) in (str, bool, int):
        return value
    raise TypeError('unexpected receipt value')


def run():
    data = native_density()
    correlations = [{'attempts': n, 'matrix_sha256': _digest(covariance(n)),
                     'port0_self_correlation': _fraction(covariance(n)[0, 0])}
                    for n in (0, 1, 2, 4, 8)]
    discounted = [{'discount_rate': z, 'matrix_sha256': _digest(native_discount(z))}
                  for z in (F(1, 2), F(1), F(4))]
    resolvents = []
    for z in (F(1, 2), F(1), F(4)):
        value = base_resolvent(z)
        resolvents.append({'discount_rate': z,
            'full_native_residual_exactly_zero': True,
            'compressed_sha256': _digest(value['compressed']),
            'effective_sha256': _digest(value['effective']),
            'difference_row_norm_upper': _outward((value['row_bound'], value['row_bound']))[1],
            'theorem_operator_bound': F(1, 8)})
    incident_spectrum()
    return _json({'schema': 'oph.observer_time_fields.spatial_dynamics.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'native': {'ports': 12, 'sector_dimension': 66, 'branches': data['branches'],
                   'stationary_susceptibility': F(5, 33)},
        'two_time_correlations': correlations, 'native_geometric_discount': discounted,
        'actual240port_dynamic_elimination': resolvents,
        'all_level_bounds': [bounds(level, terms=16) for level in (0, 1, 2, 4, 10)],
        'cell_memory_control': memory_control(), 'multi_band_clock_control': dispersion(),
        'coherent_wave_field': coherent_field(),
        'coherent_wave_bounds': [wave_bounds(level, band=16, parameter=4) for level in (0, 1, 2, 4, 10)],
        'scope': {'uniform_IID_native_seams_and_fair_endpoint_repairs': True,
            'declared_canonical_sector_stationarity': True,
            'full_initial_occupations_and_seam_current_tape_for_collective_fields': True,
            'exact_native_two_time_covariance_not_a_Gaussian_history_claim': True,
            'discount_is_record_processing_not_an_additional_scheduler': True,
            'same_source_spatial_Green_and_temporal_resolvent': True,
            'arbitrary_admitted_cross_port_matching_covered_by_bound': True,
            'incident_log_comparison_is_ordered_spectrum_not_full_clock_law': True,
            'weak_static_limit_implies_unique_temporal_semigroup': False,
            'existing_three_mode_compiler_implements_this_full_field_generator': False,
            'modular_continuation_equals_passive_native_waiting': False,
            'physical_relativistic_QFT_or_CFT_identified': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        from model.observer.scalar_response.experiment import compare_receipt
        compare_receipt(strict_json.load(args.check), result)
        print('Native temporal correlations, spatial resolvents and clock-dispersion controls replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, sort_keys=True); stream.write('\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
