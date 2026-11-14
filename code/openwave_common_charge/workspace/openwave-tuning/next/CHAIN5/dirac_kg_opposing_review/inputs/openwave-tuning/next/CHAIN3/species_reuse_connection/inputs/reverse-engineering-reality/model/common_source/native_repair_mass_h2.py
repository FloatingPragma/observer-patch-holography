"""Exact two-hop correlation control on the persistent geometric object.

The existing seam-weight cost is evaluated on the actual two equivariant
matching codes. Its raw diagonal action leaves matching syndromes; even a
separately declared matching-support average leaks from their geometric M6.
Complete first and second moments retain both orthogonal complements. No
compressed eigenvalue is programmed as a source phase or called a mass.

The cost is an additional coherent probe, not a term selected by corrected
C/E repair. This calculation supplies an exact attachment control for the
specified native object, not a physical family coupling or a new history.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
from pathlib import Path

import sympy as s

from . import native_gauge_object as objects
from . import parameter_selection as selection

ROOT = Path(__file__).resolve().parents[2]


def _matrix(value):
    return [[str(x) for x in value.row(i)] for i in range(value.rows)]


@lru_cache(maxsize=1)
def geometry():
    """Rational unknown-input frame and its exact, non-Euclidean metric."""
    projector = objects.geometric_projectors()['six']
    columns = projector.to_DM().rref()[1]
    frame = s.ImmutableMatrix(projector[:, columns])
    metric = s.ImmutableMatrix(frame.T*frame)
    if frame.cols != 6 or frame*metric.inv()*frame.T != projector:
        raise ArithmeticError('the actual geometric input frame changed')
    return {'projector': projector, 'frame': frame, 'metric': metric,
            'metric_inverse': metric.inv()}


@lru_cache(maxsize=2, typed=True)
def local_cost(power=2):
    """Derive local diagonal moments in native pair order, before projection."""
    if type(power) is not int or power not in (1, 2):
        raise ValueError('the fixed one-hop or two-hop source cost is required')
    source = selection.diagonal_cost(power)
    states = selection.code_geometry()['states']
    by_pair = {tuple(i for i, bit in enumerate(state) if bit): source['values'][j]
               for j, state in enumerate(states)}
    pairs = objects.geometry()['pairs']
    values = s.ImmutableMatrix([by_pair[pair] for pair in pairs])
    diagonals, variances = [], []
    for choice in (0, 1):
        groups = [[by_pair[pair] for pair in matching]
                  for matching in objects.equivariant_matchings(choice)]
        means = [sum(group)/6 for group in groups]
        variance = [s.factor(sum((v-m)**2 for v in group)/6)
                    for group, m in zip(groups, means, strict=True)]
        if not all(v.is_positive for v in variance):
            raise ArithmeticError('raw cost unexpectedly preserves a matching color')
        diagonals.append(s.ImmutableMatrix(s.diag(*means)))
        variances.append(s.ImmutableMatrix(s.diag(*variance)))
    return {'values': values, 'means': tuple(diagonals), 'variances': tuple(variances)}


@lru_cache(maxsize=4, typed=True)
def moments(power=2, sign=1):
    """Full native moments for H_left + sign*H_right on six unknown inputs.

    The frame has Gram G. A=G^-1 F^*HF is its compressed generator. The
    positive form L=F^*H^2F-(F^*HF)G^-1(F^*HF) retains every output outside
    the six-input domain. It splits into matching syndromes and the other
    thirty color coordinates. No projection is made during evolution.
    """
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('a declared equal or opposite local cost is required')
    data, g = local_cost(power), geometry()
    d0, d1 = data['means']; v0, v1 = data['variances']
    eye = s.eye(6)
    h = s.kronecker_product(d0, eye)+sign*s.kronecker_product(eye, d1)
    syndrome = s.kronecker_product(v0, eye)+s.kronecker_product(eye, v1)
    r, gi, p = g['frame'], g['metric_inverse'], g['projector']
    first = r.T*h*r
    second = r.T*(h*h+syndrome)*r
    color = r.T*h*(s.eye(36)-p)*h*r
    raw_syndrome = r.T*syndrome*r
    leakage = second-first*gi*first
    if leakage != color+raw_syndrome:
        raise ArithmeticError('orthogonal native complements do not exhaust the moment defect')
    if leakage.det() <= 0 or color.rank() != (0 if power == 1 else 6):
        raise ArithmeticError('the exact raw or color leakage rank changed')
    return {'color_generator': s.ImmutableMatrix(h), 'first': s.ImmutableMatrix(first),
            'second': s.ImmutableMatrix(second), 'leakage': s.ImmutableMatrix(leakage),
            'syndrome_leakage': s.ImmutableMatrix(raw_syndrome),
            'color_leakage': s.ImmutableMatrix(color)}


def complete_native_factors(power=2, sign=1):
    """All 66-by-66 native outputs, including zero rows, in exact factors.

    The binary total-load-four space has 10626 rows. A diagonal local cost leaves
    the 4356-row load-(2,2) sector invariant; the other 6270 rows are exactly
    zero for these inputs. F=(W0 tensor W1)R/6 and diagonal H therefore
    determine every amplitude for arbitrary unknown reference columns.
    """
    moments(power, sign)  # Validate the original domain before building factors.
    values = local_cost(power)['values']
    return {'left_encoder_integer': objects.color_matrix(0),
            'right_encoder_integer': objects.color_matrix(1),
            'frame': geometry()['frame'], 'encoder_divisor': s.Integer(6),
            'diagonal': s.ImmutableMatrix([values[a]+sign*values[b]
                                          for a in range(66) for b in range(66)])}


def preserving_coefficient_rank():
    """Test every linear combination of the two support-averaged H2 costs."""
    d0, d1 = local_cost(2)['means']; g = geometry()
    operators = (s.kronecker_product(d0, s.eye(6)), s.kronecker_product(s.eye(6), d1))
    columns = [((s.eye(36)-g['projector'])*h*g['frame']).reshape(216, 1)
               for h in operators]
    return s.Matrix.hstack(*columns).rank()


def report():
    g = geometry(); rows = []
    for power in (1, 2):
        data = local_cost(power)
        local = [{'matching_choice': choice,
                  'compressed_diagonal_exact': list(map(str, data['means'][choice].diagonal())),
                  'raw_syndrome_variances_exact': list(map(str, data['variances'][choice].diagonal()))}
                 for choice in (0, 1)]
        cases = []
        for sign in (1, -1):
            m = moments(power, sign)
            trace = lambda key: str(s.factor(s.trace(g['metric_inverse']*m[key])))
            cases.append({'right_cost_sign': sign,
                          'first_native_moment_exact': _matrix(m['first']),
                          'second_native_moment_exact': _matrix(m['second']),
                          'complete_native_leakage_rank': m['leakage'].rank(),
                          'complete_native_leakage_squared_frobenius': trace('leakage'),
                          'matching_syndrome_leakage_squared_frobenius': trace('syndrome_leakage'),
                          'remaining_color_leakage_rank': m['color_leakage'].rank(),
                          'remaining_color_leakage_squared_frobenius': trace('color_leakage')})
        rows.append({'seam_weight_power': power, 'local_costs': local, 'opposite_matching_pair': cases})
    rank = preserving_coefficient_rank()
    if rank != 2:
        raise ArithmeticError('a nonzero support-averaged H2 combination unexpectedly preserves M6')
    return {'classification': 'exact full-complement correlation-cost attachment control',
            'domain': 'binary 24-port total-load-four configurations; diagonal local cost preserves binary load-(2,2)',
            'native_input_columns': 6, 'complete_binary_native_total_load_four_rows': 10626,
            'diagonal_invariant_binary_load_two_two_rows': 4356,
            'retained_zero_binary_rows_in_other_load_sectors': 6270,
            'raw_color_complement_dimension': 30,
            'input_frame_exact': _matrix(g['frame']), 'input_metric_exact': _matrix(g['metric']),
            'costs': rows, 'averaged_H2_coefficient_constraint_rank': rank,
            'nonzero_averaged_H2_combination_preserves_M6': False,
            'compression_evolution_second_derivative_defect': '-G^-1 L, L=F^*H^2F-(F^*HF)G^-1(F^*HF)',
            'scope': {'coherent_cost_is_an_additional_architecture_input': True,
                      'raw_complement_removed_or_postselected': False,
                      'geometric_current_triplets_identified_as_families': False,
                      'compiled_source_history_executed_here': False,
                      'Higgs_family_coupling_or_physical_particle_mass_derived': False,
                      'measured_mass_or_programmed_pole_input': False}}


def source_bindings():
    paths = set(selection.source_bindings()) | set(objects.source_bindings()) | {
        'model/common_source/native_repair_mass_h2.py',
        'model/common_source/test_native_repair_mass_h2.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
