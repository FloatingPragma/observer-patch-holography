"""Exact complete-output enclosures for the canonical smooth continuation.

Polynomial approximants carry all sixteen records through the unbounded
Hamiltonian and both Ramsey instruments. They are verification objects, not
substitutes for the source-executed finite program. See SMOOTH_HISTORY.md.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
from pathlib import Path

from . import infinite_history as ih, protocol


@dataclass(frozen=True)
class Resolution:
    cutoff: int = 16
    preparation: int = 12
    degree: int = 16
    steps_per_unit: int = 16
    bits: int = 112

    def __post_init__(self):
        ih._integer(self.cutoff, 2, 32, 'retained total occupation')
        ih._integer(self.preparation, 0, self.cutoff, 'preparation cutoff')
        ih._integer(self.degree, 1, 32, 'Taylor degree')
        ih._integer(self.steps_per_unit, 2, 64, 'steps per unit')
        ih._integer(self.bits, 32, 192, 'coefficient precision')
        if self.preparation % 2 or self.steps_per_unit % 2:
            raise ValueError('even preparation cutoff and steps per unit required')


def _norm(vector, normalization):
    return F(ih._norm_numerator(vector.coefficients), vector.denominator**2) / normalization


def _compress(coefficients, dimension):
    if dimension is None:
        return coefficients
    return {key: value for key, value in coefficients.items()
            if key[0] < dimension and key[1] < dimension}


def _coordinate(coefficients, dimension=None):
    """B=a0+a0*, with factorial norms and compression after this action."""
    out = {}
    for (n0, n1, gauge), value in coefficients.items():
        ih._scaled_add(out, (n0+1, n1, gauge), value, 1)
        if n0:
            ih._scaled_add(out, (n0-1, n1, gauge), value, n0)
    return _compress(out, dimension)


def _series(vector, time, degree, action, denominator=1):
    """Taylor polynomial and spectral next-term bound for exp(-i time A)."""
    term, total = dict(vector.coefficients), dict(vector.coefficients)
    divisor = vector.denominator
    for power in range(1, degree+2):
        term = {key: (time.numerator*imag, -time.numerator*real)
                for key, (real, imag) in action(term).items()}
        factor = denominator*time.denominator*power
        divisor *= factor
        if power == degree+1:
            return ih._Polynomial(total, divisor//factor), ih._Polynomial(term, divisor)
        total = {key: (real*factor, imag*factor) for key, (real, imag) in total.items()}
        for key, (real, imag) in term.items():
            ih._add(total, key, real, imag)
    raise AssertionError('unreachable')


def _evolve(vector, time, normalization, resolution, dimension):
    steps = max(1, math.ceil(time*resolution.steps_per_unit))
    if dimension is None:
        return ih._advance(vector, time, steps, resolution.degree,
                           resolution.cutoff, resolution.bits, normalization)
    errors = {'Taylor': F(0), 'discarded_support': F(0), 'coefficient_rounding': F(0)}
    maximum = max((a+b for a, b, _ in vector.coefficients), default=0)

    def action(v):
        nonlocal maximum
        output = _compress(ih._hamiltonian_numerator(v), dimension)
        maximum = max(maximum, max((a+b for a, b, _ in output), default=0))
        return output

    for _ in range(steps):
        polynomial, remainder = _series(vector, time/steps, resolution.degree, action, 128)
        errors['Taylor'] += ih._root_upper(_norm(remainder, normalization))
        # The finite generator has already been compressed at EVERY power.
        # Its rectangle is retained in full; only coefficient rounding follows.
        vector, discarded, rounding = ih._cut_and_round(
            polynomial, 2*(dimension-1), resolution.bits, normalization)
        errors['discarded_support'] += discarded
        errors['coefficient_rounding'] += rounding
    return vector, errors, maximum


def _combine(left, right, right_phase):
    """(left + right_phase*i*right)/2, without branch normalization."""
    denominator = math.lcm(left.denominator, right.denominator)
    out = {key: (a*(denominator//left.denominator), b*(denominator//left.denominator))
           for key, (a, b) in left.coefficients.items()}
    for key, (a, b) in right.coefficients.items():
        factor = denominator//right.denominator
        ih._add(out, key, -right_phase*b*factor, right_phase*a*factor)
    return ih._Polynomial(out, 2*denominator)


def _pointer(vector, normalization, resolution, dimension):
    # Q0/2=B/sqrt(8). Choosing a rational LOWER coefficient is charged by
    # |delta alpha| ||Bv||, not treated as an exact irrational source angle.
    alpha_upper = ih._root_upper(F(1, 8))
    alpha = alpha_upper-F(1, ih.ROOT_SCALE)
    delta = F(1, ih.ROOT_SCALE)
    maximum = max((a+b for a, b, _ in vector.coefficients), default=0)

    def action(v):
        nonlocal maximum
        output = _coordinate(v, dimension)
        maximum = max(maximum, max((a+b for a, b, _ in output), default=0))
        return output

    minus, rem_minus = _series(vector, alpha, resolution.degree, action)
    plus, rem_plus = _series(vector, -alpha, resolution.degree, action)
    taylor = ih._root_upper((_norm(rem_minus, normalization)+_norm(rem_plus, normalization))/2)
    coordinate_vector = ih._Polynomial(action(vector.coefficients), vector.denominator)
    angle = delta*ih._root_upper(_norm(coordinate_vector, normalization))
    outputs, cuts, rounds = [], [], []
    for sign in (-1, 1):
        state = _combine(minus, plus, sign)
        state, cut, rounding = ih._cut_and_round(
            state, resolution.cutoff if dimension is None else 2*(dimension-1),
            resolution.bits, normalization)
        outputs.append(state)
        cuts.append(cut)
        rounds.append(rounding)
    errors = {'Taylor': taylor, 'coordinate_coefficient': angle,
              'discarded_support': ih._root_upper(sum((x*x for x in cuts), F(0))),
              'coefficient_rounding': ih._root_upper(sum((x*x for x in rounds), F(0)))}
    return tuple(outputs), errors, maximum


def _kick(vector, current):
    if current == 0:
        return vector
    out = {}
    for (a, b, gauge), (real, imag) in vector.coefficients.items():
        ih._add(out, (a, b, gauge), 3*real, 3*imag)
        ih._add(out, (a, b, gauge ^ 1), 4*current*imag, -4*current*real)
    return ih._Polynomial(out, 5*vector.denominator)


@dataclass(frozen=True)
class Enclosure:
    program: protocol.Program
    dimension: int | None
    resolution: Resolution
    normalization: F
    branches: dict
    stages: list
    vector_error: F
    preparation_error: F
    maximum_support: int
    seal: str = ''

    @property
    def output_error(self):
        return self.preparation_error+self.vector_error+self.vector_error**2/2


def _seal(enclosure):
    payload = [enclosure.program.specification(), enclosure.dimension,
               enclosure.resolution.__dict__, str(enclosure.normalization),
               enclosure.stages, str(enclosure.vector_error), str(enclosure.preparation_error),
               enclosure.maximum_support,
               [[list(record), str(state.denominator),
                 [[*key, str(a), str(b)] for key, (a, b) in sorted(state.coefficients.items())]]
                for record, state in enclosure.branches.items()]]
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _validate(enclosure):
    """Check an in-process proof object's integrity, not external authentication.

    Published receipts are verified by complete exact recomputation, never by
    trusting this public digest or importing user-supplied error claims.
    """
    if type(enclosure) is not Enclosure:
        raise ValueError('an executed exact enclosure required')
    if not enclosure.seal or enclosure.seal != _seal(enclosure):
        raise ValueError('enclosure integrity digest does not match its contents')
    reads = sum(op.kind in protocol.READS for op in enclosure.program.operations)
    if tuple(enclosure.branches) != tuple(product((0, 1), repeat=reads)):
        raise ValueError('every ordered record, including zero branches, is required')
    return enclosure


def enclose(program, *, dimension=None, resolution=None):
    """Enclose the whole declared program; None means infinite occupation.

    The preparation is the specified blank gauge/free scalar vacuum with a
    native-current kick. Arbitrary supplied work states are outside this API.
    """
    if resolution is None:
        resolution = Resolution()
    if type(program) is not protocol.Program or type(resolution) is not Resolution:
        raise ValueError('typed canonical Program and certificate Resolution required')
    program = protocol.validate_program(program)
    resolution = replace(resolution)
    if program.coupling != F(1, 4):
        raise ValueError('the exact common Hamiltonian requires coupling 1/4')
    if dimension is not None and (type(dimension) is not int or dimension != program.dimension):
        raise ValueError('finite comparison must use the program occupation cutoff')
    if any(op.kind == 'sharp_read' for op in program.operations):
        raise ValueError('sharp pointers require their separate nonpolynomial output theorem')
    if (sum((op.time for op in program.operations), F(0)) > 4
            or sum(op.kind == 'smooth_read' for op in program.operations) > 2
            or sum(op.kind in protocol.READS for op in program.operations) > 4):
        raise ValueError('certificate execution budget: model time <=4, two smooth reads, four records')
    initial = ih._initial(resolution.preparation if dimension is None else dimension-2)
    normalization = F(ih._norm_numerator(initial.coefficients), initial.denominator**2)
    preparation_error = F(0)
    if dimension is None:
        sqrt5_lower = F(math.isqrt(5*ih.ROOT_SCALE**2), ih.ROOT_SCALE)
        preparation_error = ih._root_upper(1-F(4, 9)*sqrt5_lower*normalization)
    branches = {(): _kick(initial, program.native_current)}
    stages, total_error = [], F(0)
    maximum = max((a+b for a, b, _ in initial.coefficients), default=0)
    for index, operation in enumerate(program.operations):
        following, local_errors, components = {}, [], {}
        for record, vector in branches.items():
            errors, support = {}, 0
            if operation.kind == 'evolve':
                state, errors, support = _evolve(vector, operation.time, normalization, resolution, dimension)
                following[record] = state
            elif operation.kind == 'smooth_read':
                outputs, errors, support = _pointer(vector, normalization, resolution, dimension)
                following.update({(*record, bit): state for bit, state in enumerate(outputs)})
            elif operation.kind == 'gauge_read':
                following.update({(*record, bit): ih._first_record(vector, bit, False) for bit in (0, 1)})
            elif operation.kind == 'feedback':
                following[record] = _kick(vector, int(operation.enabled and record[operation.record]))
            else:
                raise ValueError('unsupported certificate operation')
            local_errors.append(sum(errors.values(), F(0)))
            for key, value in errors.items():
                components.setdefault(key, []).append(value)
            maximum = max(maximum, support)
        increment = ih._root_upper(sum((x*x for x in local_errors), F(0)))
        total_error += increment
        branches = following
        stages.append({'index': index, 'kind': operation.kind,
                       'branches': len(branches), 'vector_error_increment_upper': str(increment),
                       'cumulative_vector_error_upper': str(total_error),
                       'component_l2_upper': {key: str(ih._root_upper(sum((x*x for x in values), F(0))))
                                              for key, values in sorted(components.items())}})
    result = Enclosure(program, dimension, resolution, normalization, branches, stages,
                       total_error, preparation_error, maximum)
    return replace(result, seal=_seal(result))


def _inner_squared(left, right):
    real, imag = 0, 0
    for (a, b, gauge), (lr, li) in left.coefficients.items():
        rr, ri = right.coefficients.get((a, b, gauge), (0, 0))
        weight = math.factorial(a)*math.factorial(b)
        real += weight*(lr*rr+li*ri)
        imag += weight*(lr*ri-li*rr)
    return F(real*real+imag*imag, (left.denominator*right.denominator)**2)


def distance_interval(left, right):
    """Rigorous complete CQ distance of the two exact normalized histories."""
    _validate(left)
    _validate(right)
    target = lambda x: (x.program.operations, x.program.coupling, x.program.native_current)
    if target(left) != target(right):
        raise ValueError('enclosures of the same continuous target program required')
    if left.branches.keys() != right.branches.keys():
        raise ValueError('every retained record must occur on both sides')
    lower = upper = F(0)
    for key in left.branches:
        u, v = left.branches[key], right.branches[key]
        p, q = _norm(u, left.normalization), _norm(v, right.normalization)
        overlap = _inner_squared(u, v)/(left.normalization*right.normalization)
        squared = (p+q)**2-4*overlap
        if squared < 0:
            raise ArithmeticError('exact rank-one distance violated Cauchy-Schwarz')
        root = ih._root_upper(squared)
        upper += root/2
        lower += max(F(0), root-F(1, ih.ROOT_SCALE))/2
    error = left.output_error+right.output_error
    return max(F(0), lower-error), min(F(1), upper+error)


def summary(enclosure):
    _validate(enclosure)
    rows = []
    for record, state in enclosure.branches.items():
        coefficients = [[*key, str(a), str(b)] for key, (a, b) in sorted(state.coefficients.items())]
        digest = hashlib.sha256(json.dumps([str(state.denominator), coefficients], separators=(',', ':')).encode()).hexdigest()
        rows.append({'record': list(record), 'polynomial_sha256': digest,
                     'nonzero_coefficients': len(coefficients),
                     'mass': ih._interval(_norm(state, enclosure.normalization), _norm(state, enclosure.normalization))})
    return {'dimension': enclosure.dimension, 'resolution': dict(enclosure.resolution.__dict__),
            'preparation_trace_distance_upper': str(enclosure.preparation_error),
            'polynomial_vector_error_upper': str(enclosure.vector_error),
            'complete_output_half_trace_norm_error_upper': str(enclosure.output_error),
            'maximum_expanded_total_occupation': enclosure.maximum_support,
            'stages': deepcopy(enclosure.stages), 'branches': rows}


def probability_interval(enclosure, slot, bit):
    """An unconditional recorded-bit probability, with all other records summed."""
    _validate(enclosure)
    reads = sum(op.kind in protocol.READS for op in enclosure.program.operations)
    if type(slot) is not int or not 0 <= slot < reads or type(bit) is not int or bit not in (0, 1):
        raise ValueError('an existing record slot and binary value required')
    mass = sum((_norm(state, enclosure.normalization) for record, state in enclosure.branches.items()
                if record[slot] == bit), F(0))
    error = enclosure.preparation_error+2*enclosure.vector_error+enclosure.vector_error**2
    return max(F(0), mass-error), min(F(1), mass+error)


def _execution_check(finite, temporal):
    """Independent source/product execution diagnostics, not roundoff bounds."""
    import numpy as np
    from . import composition
    _validate(finite)
    source = composition._executed(finite.program)
    independent = composition._executed(finite.program, 'product')
    gate = F(1, 10**8)

    def embedded(vector):
        target = np.zeros(finite.program.work_dimension, complex)
        for (n0, n1, gauge), (real, imag) in vector.coefficients.items():
            target[4*(n0+finite.dimension*n1)+gauge] = (
                complex(real, imag)/vector.denominator
                *math.sqrt(float(F(math.factorial(n0)*math.factorial(n1))/finite.normalization)))
        return target

    # Preparation is a separate service: an erroneous shared preparation
    # must not be hidden inside the much larger Hamiltonian splitting budget.
    prepared = embedded(ih._initial(finite.dimension-2)).reshape(-1, 1)
    for history, backend in ((source, 'source'), (independent, 'product')):
        protocol.validate_checkpoint(finite.program, history, backend)
        if (history.program != finite.program or history.preparation_mode != 'declared_vacuum'
                or history.cursor != len(finite.program.operations)
                or history.initial_state.shape != prepared.shape
                or not np.isfinite(history.initial_state).all()
                or float(np.linalg.norm(history.initial_state-prepared)) > float(gate)):
            raise ArithmeticError('executed preparation differs from the exact declared vacuum')
    comparison = protocol.compare_histories(source, independent)
    if any(not math.isfinite(value) or value < 0 or value > float(gate)
           for value in comparison.values()):
        raise ArithmeticError('emitted source maps failed the independent product comparison')
    squared = []
    for branch in independent.branches:
        target = embedded(finite.branches[branch.record])
        squared.append(float(np.linalg.norm(branch.state[:, 0]-target)**2))
    residual = math.sqrt(math.fsum(squared))
    if not math.isfinite(residual) or residual > float(temporal+finite.vector_error+gate):
        raise ArithmeticError('finite executed history is detached from the exact polynomial target')
    return {'source_vs_independent_product_residual_gate': str(gate),
            'independent_preparation_check_pass': True,
            'source_product_operator_and_branch_checks_pass': True,
            'product_vs_exact_continuous_polynomial_check_pass': True,
            'diagnostic_is_outward_roundoff_certificate': False,
            'realized_path_resource_maxima': {
                key: max(dict(branch.resources)[key] for branch in source.branches)
                for key in protocol.RESOURCE_KEYS},
            'expanded_tree_resources': dict(source.expanded_tree_resources)}


SOURCE_PATHS = (
    'model/common_source/smooth_history.py', 'model/common_source/infinite_history.py',
    'model/common_source/protocol.py', 'model/common_source/continuation_bounds.py',
    'model/common_source/composition.py', 'model/common_source/smooth_protocol.py',
    'model/common_source/scalar_coupling.py', 'model/common_source/experiment.py',
    'model/common_source/gauge_matter.py', 'model/common_source/quantum_native.py',
    'model/common_source/encoded_memory.py', 'model/common_source/encoded_observer.py',
    'model/common_source/native_routing.py', 'model/observer/source_scalar/compiler.py',
    'model/observer/source_scalar/field.py', 'model/observer/source_scalar/bounds.py',
    'model/observer/source_scalar/cutoff.py', 'model/requirements.txt',
    'model/common_source/composition_results.json',
    'model/common_source/quantum_native_results.json',
    'model/common_source/encoded_observer_results.json',
)


def run():
    """Recompute proof objects and execute the same six finite source programs."""
    from . import continuation_bounds as bounds
    rows, lifted = [], {2: [], 4: []}
    fine, coarse = Resolution(), Resolution(12, 8, 12, 12, 96)
    for current in (-1, 0, 1):
        program = protocol.continuation_program(2, 512, native_current=current)
        infinite = enclose(program, resolution=fine)
        finite_rows = []
        for dimension in (2, 4):
            finite_program = protocol.continuation_program(dimension, 512, native_current=current)
            finite = enclose(finite_program, dimension=dimension, resolution=fine)
            interval = distance_interval(infinite, finite)
            budget = bounds.program_bound(finite_program)
            temporal = F(budget['ideal_temporal_instrument_error_upper'])
            product_interval = max(F(0), interval[0]-temporal), min(F(1), interval[1]+temporal)
            if not 0 < interval[0] <= interval[1] < F(3, 10) or product_interval[1] >= 1:
                raise ArithmeticError('complete prepared-history occupation comparison is not nonvacuous')
            lifted[dimension].append(product_interval[1])
            finite_rows.append({
                'dimension': dimension, 'program': finite_program.specification(),
                'continuous_to_infinite_complete_output_distance_interval': ih._interval(*interval),
                'ideal_temporal_error_upper': str(temporal),
                'ideal_product_to_infinite_complete_output_distance_interval': ih._interval(*product_interval),
                'finite_polynomial_certificate': summary(finite),
                'source_execution': _execution_check(finite, temporal),
                'physical_service_coverage': bounds.service_error_budget(
                    product_interval[1], required_stage_keys=bounds.required_service_keys(finite_program.operations)),
            })
        row = {'native_current': current, 'infinite_polynomial_certificate': summary(infinite),
               'infinite_clock_one_probability_interval': ih._interval(*probability_interval(infinite, 3, 1)),
               'finite_comparisons': finite_rows}
        if current == 0:
            refinement = enclose(program, resolution=coarse)
            refinement_interval = distance_interval(infinite, refinement)
            if refinement_interval[0] != 0 or not infinite.output_error < refinement.output_error:
                raise ArithmeticError('same infinite-history refinement certificates disagree')
            row['coarser_polynomial_certificate'] = summary(refinement)
            row['same_target_refinement_distance_interval'] = ih._interval(*refinement_interval)
        rows.append(row)
    root = Path(__file__).resolve().parents[2]
    return {'schema': 'oph.common_source.smooth_history_enclosure.v1',
            'source_sha256': {path: hashlib.sha256((root/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
            'currents': rows,
            'complete_native_composition': [
                {'dimension': d, 'maximum_current_ideal_product_output_error_upper':
                 ih._interval(F(0), max(values))[1],
                 'physical_output_error_upper': None}
                for d, values in lifted.items()],
            'scope': {
                'preparation': 'independent blank gauge and declared free scalar vacuum; all three native-current kicks',
                'retained_records_per_work_history': 16,
                'complete_quantum_outputs_compared': True,
                'native_reference_extension': 'full-carrier C/E and protected/passive logical modes with independently prepared work',
                'arbitrary_initial_native_work_correlations_included': False,
                'uniform_moment_certificate_replaced': False,
                'source_splitting_in_product_comparison_included': True,
                'floating_or_physical_service_errors_certified': False,
                'spatial_refinement_or_physical_clock_calibration_claimed': False,
                'new_source_law_selected': False},
            'all_pass': True}


def main():
    import argparse
    from tools.strict_json import load
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    args = parser.parse_args()
    expected = load(args.check) if args.check else None
    result = run()
    if args.check:
        ih.compare_receipt(expected, result)
    if args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(ih.canonical(result))
    print(ih.canonical({'all_pass': result['all_pass'],
                        'complete_native_composition': result['complete_native_composition']}))


if __name__ == '__main__':
    main()
