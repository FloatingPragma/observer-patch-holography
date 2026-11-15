"""Source-executed local number noise on the continuing common register.

The supplied dimensionless variance is not a calibrated laboratory rate or
duration. Unobserved environment alternatives are summed without selecting a
branch; existing observer records and every untouched reference remain.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction as F
import math
import sys

import numpy as np

from ..anharmonic_transfer import environment as source
from ..observer.source_scalar import compiler
from ..scalar_observations.experimental_readout import canonical_source as canonical
from ..scalar_observations.experimental_readout import executable_source as native
from . import protocol


TOLERANCE = 3e-12
MAX_CERTIFICATE_INTEGER_DIGITS = 4000


def _serialization_limit():
    # Observe a caller's stricter interpreter policy without mutating it.
    active = sys.get_int_max_str_digits()
    return min(MAX_CERTIFICATE_INTEGER_DIGITS, active or MAX_CERTIFICATE_INTEGER_DIGITS)


def _check_fraction_size(value):
    limit = _serialization_limit()
    for integer in (value.numerator, value.denominator):
        # 30103/100000 > log10(2). This inexpensive integer-only estimate
        # conservatively bounds decimal length without formatting the value.
        digits = max(1, (abs(integer).bit_length()*30103+99999)//100000)
        if digits > limit:
            raise ValueError(f'noise exact certificate exceeds the {limit}-digit integer serialization budget')


def _certificate_preflight(dimension, event):
    """Check even the uncapped remainder retained by the inherited producer."""
    _check_fraction_size(event.variance)
    represented = F(float(event.variance))
    order = event.order
    entry = F(math.factorial(order), math.factorial(2*order))*(
        represented*(dimension-1)**2)**order
    coefficient = F((dimension-1)**2, 4)
    quadrature = min(F(1), F(dimension, 2)*entry)
    conversion = min(F(1), abs(event.variance-represented)*coefficient)
    combined = min(F(1), quadrature+conversion)
    effect = F(0) if order == 1 else min(F(1), represented*coefficient)
    for value in (represented, entry, quadrature, conversion, combined, effect,
                  F(0) if order == 1 else min(F(1), event.variance*coefficient)):
        _check_fraction_size(value)
    return combined, effect


@dataclass(frozen=True)
class PhaseNoise:
    after: int
    site: int
    variance: F
    order: int = 16

    def __post_init__(self):
        if type(self.after) is not int or self.after < 0:
            raise ValueError('noise cursor must be a nonnegative integer')
        if type(self.site) is not int or self.site not in (0, 1):
            raise ValueError('noise site must be zero or one')
        if type(self.variance) not in (int, F) or self.variance < 0:
            raise ValueError('noise variance must be a nonnegative exact integer or Fraction')
        variance = F(self.variance.numerator, self.variance.denominator)
        _check_fraction_size(variance)
        try:
            finite = math.isfinite(float(variance))
        except OverflowError:
            finite = False
        if not finite:
            raise ValueError('noise variance must have a finite source representation')
        if type(self.order) is not int or not 1 <= self.order <= 64:
            raise ValueError('noise quadrature order must be an integer from 1 through 64')
        object.__setattr__(self, 'variance', variance)


@dataclass(frozen=True)
class Environment:
    events: tuple[PhaseNoise, ...] = ()

    def __post_init__(self):
        if type(self.events) is not tuple or any(type(event) is not PhaseNoise for event in self.events):
            raise ValueError('environment requires an immutable tuple of PhaseNoise events')
        events = tuple(PhaseNoise(event.after, event.site, event.variance, event.order)
                       for event in self.events)
        addresses = tuple((event.after, event.site) for event in events)
        if addresses != tuple(sorted(set(addresses))):
            raise ValueError('noise events must have strictly sorted unique cursor/site addresses')
        object.__setattr__(self, 'events', events)

    def specification(self):
        checked = Environment(self.events)
        return {'schema': 'oph.common_source.environment.v1',
                'events': [{'after': event.after, 'site': event.site,
                            'variance': str(event.variance), 'order': event.order}
                           for event in checked.events]}


def validated(program, environment):
    """Validate and detach a schedule, including frozen objects mutated by a caller."""
    checked = protocol.validate_program(program)
    if type(environment) is not Environment:
        raise ValueError('a canonical Environment is required')
    detached = Environment(environment.events)
    if any(event.after > len(checked.operations) for event in detached.events):
        raise ValueError('noise cursor exceeds the complete program')
    totals = [F(0), F(0)]
    for event in detached.events:
        radii = _certificate_preflight(checked.dimension, event)
        totals = [min(F(1), previous+radius) for previous, radius in zip(totals, radii, strict=True)]
        for value in totals:
            _check_fraction_size(value)
    return detached


def _event(program, event):
    if type(event) is not PhaseNoise:
        raise ValueError('a PhaseNoise event is required')
    return validated(program, Environment((event,))).events[0]


def _snapshot(program, event):
    """Copy and check executed cache contents before using their Schur matrix.

    Independent expected phases are guards only. The returned matrix always
    comes from the copied source contraction, never the analytic target.
    """
    d, variance = program.dimension, float(event.variance)
    channel = source.channel(d, 1., variance, event.order)
    try:
        if type(channel) is not source.PhaseDiffusion or channel.dimension != d:
            raise ValueError('source environment dimension mismatch')
        if (not isinstance(channel.weights, np.ndarray) or channel.weights.dtype.kind not in 'fiu'
                or not isinstance(channel.schur, np.ndarray) or channel.schur.dtype.kind not in 'fiuc'
                or type(channel.unitaries) is not tuple or len(channel.unitaries) != event.order
                or any(not isinstance(unitary, np.ndarray) or unitary.dtype.kind not in 'fiuc'
                       for unitary in channel.unitaries)):
            raise ValueError('malformed source environment arrays')
        weights = np.array(channel.weights, dtype=float, copy=True)
        schur = np.array(channel.schur, dtype=complex, copy=True)
        unitaries = tuple(np.array(unitary, dtype=complex, copy=True) for unitary in channel.unitaries)
        metadata = deepcopy(channel.metadata)
        if (weights.shape != (event.order,) or schur.shape != (d, d)
                or not np.isfinite(weights).all() or not np.isfinite(schur).all()
                or np.any(weights <= 0) or abs(math.fsum(weights)-1) > TOLERANCE):
            raise ValueError('source environment weights or Schur matrix are invalid')
        if (metadata['dimension'] != d or metadata['quadrature_order'] != event.order
                or metadata['rate_per_second'] != 1. or metadata['duration_seconds'] != variance
                or len(metadata['branches']) != event.order
                or metadata['all_phase_unitaries_source_executed'] is not True
                or metadata['independent_classical_environment_draw_per_application'] is not True
                or metadata['environment_outcome_discarded_without_selection'] is not True
                or metadata['analytic_dephasing_schur_matrix_inserted'] is not False):
            raise ValueError('source environment metadata differs from requested channel')
        nodes, ideal_weights = source.canonical_quadrature(event.order)
        expected_weights = np.array([float(weight) for weight in ideal_weights])
        expected_weights /= math.fsum(expected_weights)
        if not np.array_equal(weights, expected_weights):
            raise ValueError('source environment changed its declared quadrature weights')
        context = canonical._context()
        sigma = context.sqrt(canonical._parameter(context, variance, 'variance'))
        bits = d.bit_length()-1
        phases = []
        expected_resources = None
        for index, (unitary, branch) in enumerate(zip(unitaries, metadata['branches'], strict=True)):
            phase = sigma*context.mpf(nodes[index])
            word = [compiler.pulse('I', float(-phase*(d-1)/2), [0])]
            word.extend(compiler.pulse('Z', float(phase*2**bit/2), [bit]) for bit in range(bits))
            resources = {key: compiler.resources(word)[key] for key in native.COUNT_KEYS}
            if expected_resources is None:
                expected_resources = resources
            target = np.diag(np.exp(-1j*float(phase)*np.arange(d)))
            if (unitary.shape != (d, d) or not np.isfinite(unitary).all()
                    or np.linalg.norm(unitary-target, 2) > TOLERANCE
                    or np.linalg.norm(unitary-np.diag(np.diag(unitary)), 2) > TOLERANCE
                    or np.max(np.abs(np.abs(np.diag(unitary))-1)) > TOLERANCE):
                raise ValueError('source environment phase branch is invalid')
            if (type(branch['branch']) is not int or branch['branch'] != index
                    or branch['standard_normal_node'] != nodes[index]
                    or branch['ideal_weight_decimal'] != ideal_weights[index]
                    or branch['weight'] != float(weights[index]) or branch['phase'] != float(phase)
                    or branch['source']['word_sha256'] != native._word_hash(word)
                    or any(type(branch['source']['resources'][key]) is not int
                           or branch['source']['resources'][key] != count
                           for key, count in resources.items())):
                raise ValueError('source environment branch binding changed')
            phases.append(np.diag(unitary))
        diagonal = np.asarray(phases)
        reconstructed = (diagonal.T*weights) @ diagonal.conj()
        if np.max(np.abs(reconstructed-schur)) > TOLERANCE:
            raise ValueError('source environment Schur matrix differs from all executed branches')
        for key in ('resources_per_realization', 'resources_weighted_mean_per_application',
                    'resources_maximum_per_application', 'resources_for_all_branch_contractions'):
            multiplier = event.order if key == 'resources_for_all_branch_contractions' else 1
            if (metadata[key] != {name: multiplier*value for name, value in expected_resources.items()}
                    or any(type(value) is not int for value in metadata[key].values())):
                raise ValueError('source environment resource accounting changed')
        if metadata['ideal_quadrature_error'] != source.quadrature_bound(d, 1., variance, event.order):
            raise ValueError('source environment ideal certificate changed')
    except (KeyError, TypeError, AttributeError, OverflowError) as error:
        raise ValueError('malformed source environment cache') from error
    return schur, metadata


def _original(value):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('Boolean values are not operator amplitudes')
    if isinstance(value, (list, tuple)):
        for item in value:
            _original(item)


def apply(rho, program, event, backend='source'):
    """Apply a complete local channel to a work-major/reference-minor operator.

    Work index is gauge+4*(n0+d*n1). Finite matrix units and unnormalized
    density blocks are admitted; no positivity projection or normalization
    changes the supplied operator. A source cache is validated on every call.
    """
    event = _event(program, event)
    if type(backend) is not str or backend not in ('source', 'gaussian'):
        raise ValueError('environment backend must be source or gaussian')
    if np.ma.isMaskedArray(rho):
        raise ValueError('masked operators are not admitted')
    _original(rho)
    raw = np.asarray(rho)
    size = program.work_dimension
    if (raw.ndim != 2 or raw.shape[0] != raw.shape[1] or not raw.shape[0]
            or raw.shape[0] % size or raw.dtype.kind not in 'fiuc'
            or not np.isfinite(raw).all()):
        raise ValueError('finite square operator on work times reference required')
    with np.errstate(over='ignore', invalid='ignore'):
        operator = np.array(raw, dtype=complex, copy=True)
    if not np.isfinite(operator).all():
        raise ValueError('operator must remain finite in the execution representation')
    if backend == 'source':
        schur, _ = _snapshot(program, event)
    else:
        # Independent analytic comparison only, evaluated from the exact
        # declared variance rather than an inserted source prediction map.
        def entry(gap):
            exponent = event.variance*gap*gap/2
            try:
                return math.exp(-float(exponent))
            except OverflowError:
                return 0.
        schur = np.array([[entry(i-j) for j in range(program.dimension)]
                         for i in range(program.dimension)])
    reference = len(operator)//size
    occupation = (np.arange(size)//(4*program.dimension**event.site)) % program.dimension
    occupation = np.repeat(occupation, reference)
    return operator*schur[occupation[:, None], occupation[None, :]]


def _bounds(program, event):
    d, variance, order = program.dimension, event.variance, event.order
    represented = F(float(variance))
    coefficient = F((d-1)**2, 4)
    conversion = min(F(1), abs(variance-represented)*coefficient)
    quadrature = F(source.quadrature_bound(d, 1., float(variance), order)
                   ['exact_half_diamond_error_upper'])
    # M=1 has its sole node at zero, hence its exact finite channel is identity.
    identity = F(0) if order == 1 else min(F(1), represented*coefficient)
    return {
        'after': event.after, 'site': event.site, 'quadrature_order': order,
        'variance_exact': str(variance), 'represented_variance_exact': str(represented),
        'ideal_quadrature_half_diamond_upper': str(quadrature),
        'variance_conversion_half_diamond_upper': str(conversion),
        'ideal_finite_to_gaussian_half_diamond_upper': str(min(F(1), quadrature+conversion)),
        'ideal_exact_variance_finite_to_identity_half_diamond_upper': str(
            F(0) if order == 1 else min(F(1), variance*coefficient)),
        'ideal_finite_to_identity_half_diamond_upper': str(identity),
    }


def event_metadata(program, event):
    event = _event(program, event)
    _, metadata = _snapshot(program, event)
    return {
        **_bounds(program, event),
        'source': metadata,
        'source_site_wires': list(range(2+event.site*(program.dimension.bit_length()-1),
                                        2+(event.site+1)*(program.dimension.bit_length()-1))),
        'resources_per_realization': deepcopy(metadata['resources_per_realization']),
        'resources_for_all_branch_contractions': deepcopy(metadata['resources_for_all_branch_contractions']),
        'independent_environment_draws': 1,
        'variance_units': 'dimensionless declared phase variance',
        'source_api_parameterization': 'rate=1 and duration=float(variance); no physical seconds assigned',
        'maximum_exact_certificate_integer_digits': _serialization_limit(),
        'all_environment_alternatives_summed_without_selection': True,
        'floating_execution_error_included': False,
        'physical_noise_parameters_or_service_errors_certified': False,
    }


def schedule_budget(program, environment):
    """Exact finite-history radii and source costs per realized observer path.

    Tree-level repeated applications across different visible records are a
    separate runtime count. Host matrix reuse removes none of these services.
    """
    environment = validated(program, environment)
    events = [event_metadata(program, event) for event in environment.events]
    totals = {key: str(min(F(1), sum((F(row[key]) for row in events), F(0))))
              for key in ('ideal_finite_to_gaussian_half_diamond_upper',
                          'ideal_finite_to_identity_half_diamond_upper')}
    resources = {key: {name: sum(row[key][name] for row in events) for name in native.COUNT_KEYS}
                 for key in ('resources_per_realization', 'resources_for_all_branch_contractions')}
    return {
        'schema': 'oph.common_source.open_environment_budget.v1', 'events': events,
        **totals, **resources, 'independent_environment_draws': len(events),
        'arbitrary_inputs_and_reference_systems_included': True,
        'all_visible_records_and_quantum_outputs_retained': True,
        'noise_draws_are_committed_observer_records': False,
        'floating_execution_error_included': False,
        'physical_noise_parameters_or_service_errors_certified': False,
        'simultaneous_lindblad_evolution_claimed': False,
        'scope': 'independent phase channels at declared checkpoints of the same finite program',
    }
