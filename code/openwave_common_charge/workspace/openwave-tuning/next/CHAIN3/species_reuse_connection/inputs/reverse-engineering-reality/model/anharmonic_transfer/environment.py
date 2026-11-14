"""Finite source-executed environment for the supplied number-dephasing law.

An independent classical draw selects a positive Gauss--Hermite branch. Each
branch executes exp(-i*phase*N) with native source pulses. Discarding the draw
gives a finite CPTP channel, including on states entangled with the spin. The
Gaussian law and its rate remain supplied apparatus; neither is derived here.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from numbers import Integral

import numpy as np

from ..observer.source_scalar import compiler as c
from ..scalar_observations.experimental_readout import canonical_source as canonical
from ..scalar_observations.experimental_readout import executable_source as native


DEFAULT_QUADRATURE_ORDER = 16
MAXIMUM_QUADRATURE_ORDER = 64


def _order(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or not 1 <= value <= MAXIMUM_QUADRATURE_ORDER):
        raise ValueError('quadrature order must be an integer from 1 through 64')
    return int(value)


def _parameters(dimension, rate, duration, quadrature_order):
    dimension = native._dimension(dimension)
    rate = native._real(rate, 'number dephasing rate')
    duration = native._real(duration, 'number dephasing duration')
    if rate < 0 or duration < 0:
        raise ValueError('nonnegative rate and duration required')
    return dimension, rate, duration, _order(quadrature_order)


def quadrature_bound(dimension, rate, duration, quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Uniform ideal half-diamond bound, valid with any reference system.

    For the standard-normal M-node Gaussian rule, the characteristic-function
    error at frequency k*sqrt(rate*duration) is at most
    M! * (rate*duration*k*k)**M / (2*M)!. The trace norm of the blockwise error
    is at most d times its largest scalar error, even with a reference system.
    Exact rational arithmetic is used on the binary64 input parameters.
    Floating source execution and quadrature construction are not certified
    by this ideal mathematical bound.
    """
    dimension, rate, duration, order = _parameters(dimension, rate, duration, quadrature_order)
    variance_span = Fraction(rate) * Fraction(duration) * (dimension-1)**2
    entry = Fraction(math.factorial(order), math.factorial(2*order))*variance_span**order
    half_diamond = min(Fraction(1), Fraction(dimension, 2)*entry)
    upper = float(half_diamond)
    if 0 < upper < 1:
        upper = min(1., math.nextafter(upper, math.inf))
    elif half_diamond > 0 and upper == 0:
        upper = math.nextafter(0., math.inf)
    return {
        'dimension': dimension, 'quadrature_order': order,
        'rate_per_second': rate, 'duration_seconds': duration,
        'half_diamond_error_upper': upper,
        'exact_half_diamond_error_upper': str(half_diamond),
        'exact_maximum_schur_entry_error_upper': str(entry),
        'formula': 'min(1, d*M!*(gamma*t*(d-1)^2)^M/(2*(2*M)!))',
        'arbitrary_reference_system_included': True,
        'scope': 'ideal Gaussian quadrature and exact native phase gates',
        'floating_execution_error_included': False,
        'physical_noise_law_or_rate_certified': False,
    }


def canonical_quadrature(quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Normal-weight Gaussian nodes from the fixed-precision Jacobi matrix.

    The off-diagonal entries sqrt(j) are the standard-normal orthonormal
    Hermite recurrence. Squared first eigenvector components give positive
    probability weights. No experimental observations enter this construction.
    """
    return _canonical_quadrature(_order(quadrature_order))


@lru_cache(maxsize=16)
def _canonical_quadrature(order):
    context = canonical._context()
    jacobi = context.matrix(order)
    for index in range(1, order):
        jacobi[index-1, index] = jacobi[index, index-1] = context.sqrt(index)
    nodes, vectors = context.eigsy(jacobi)
    weights = [vectors[0, index]**2 for index in range(order)]
    if any(weight <= 0 for weight in weights):
        raise ArithmeticError('Gaussian environment needs positive branch weights')
    return (tuple(context.nstr(value, canonical.DECIMAL_DIGITS) for value in nodes),
            tuple(context.nstr(value, canonical.DECIMAL_DIGITS) for value in weights))


@dataclass(frozen=True)
class PhaseDiffusion:
    """An enumerated classical environment and its contracted finite channel."""

    dimension: int
    weights: np.ndarray
    unitaries: tuple
    schur: np.ndarray
    metadata: dict

    @property
    def half_diamond_error_bound(self):
        return self.metadata['ideal_quadrature_error']['half_diamond_error_upper']

    def apply(self, rho):
        """Trace over the independently refreshed environment, keeping spin.

        The Schur matrix is assembled only from the executed branch phases.
        Tiling it represents identity on every retained reference register;
        this contraction is exactly the sum of the finite branch channels.
        """
        rho = np.asarray(rho, dtype=complex)
        if (rho.ndim != 2 or rho.shape[0] != rho.shape[1] or not rho.shape[0]
                or rho.shape[0] % self.dimension or not np.isfinite(rho).all()):
            raise ValueError('finite square density operator on reference times oscillator required')
        reference = rho.shape[0] // self.dimension
        return rho*np.tile(self.schur, (reference, reference))


def channel(dimension, rate, duration, quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Compile every branch, then contract a positive finite phase mixture.

    Repeated use shares the host matrix, but each application needs a fresh
    independent classical environment draw and incurs its full native word.
    The frozen zero-rate case uses the same declared branch schedule.
    """
    return _channel(*_parameters(dimension, rate, duration, quadrature_order))


@lru_cache(maxsize=64)
def _channel(dimension, rate, duration, order):
    nodes, weights_decimal = canonical_quadrature(order)
    context = canonical._context()
    sigma = context.sqrt(canonical._parameter(context, rate, 'rate')
                         * canonical._parameter(context, duration, 'duration'))
    weights = np.array([float(value) for value in weights_decimal])
    raw_weight_sum = math.fsum(weights)
    weights /= raw_weight_sum
    # Classical weights are rounded supplied probabilities, not a quantum
    # postselection. Preserve all positive branches and disclose the drift.
    bits = dimension.bit_length()-1
    phases, unitaries, records = [], [], []
    projection_changes = []
    for index, node in enumerate(nodes):
        phase = sigma*context.mpf(node)
        word = [c.pulse('I', float(-phase*(dimension-1)/2), [0])]
        word.extend(c.pulse('Z', float(phase*(2**bit)/2), [bit]) for bit in range(bits))
        target = np.diag(np.exp(-1j*float(phase)*np.arange(dimension)))
        executed, metadata = native._contract(word, bits, target,
                                              name='environment_number_phase')
        # Native Z/I words are exactly diagonal. The contraction's existing
        # polar correction can leave machine-size off-diagonal entries.
        # Represent that exact invariant explicitly, and retain its change.
        diagonal = np.diag(executed).copy()
        if np.any(np.abs(diagonal) == 0):
            raise ArithmeticError('source phase contraction lost a diagonal entry')
        diagonal /= np.abs(diagonal)
        represented = np.diag(diagonal)
        change = float(np.linalg.norm(represented-executed, 2))
        if not math.isfinite(change) or change > 2e-12:
            raise ArithmeticError('number-phase diagonal projection exceeds roundoff guard')
        represented.flags.writeable = False
        phases.append(diagonal)
        unitaries.append(represented)
        projection_changes.append(change)
        records.append({
            'branch': index, 'standard_normal_node': node,
            'ideal_weight_decimal': weights_decimal[index],
            'weight': float(weights[index]), 'phase': float(phase),
            'source': metadata, 'diagonal_roundoff_projection_norm': change,
        })
    phase_array = np.asarray(phases)
    # This is contraction of executed Kraus branches, not insertion of the
    # analytic exp(-gamma*t*(n-m)^2/2) law. The latter appears only in tests.
    schur = (phase_array.T*weights) @ phase_array.conj()
    if not np.isfinite(schur).all() or np.max(np.abs(np.diag(schur)-1)) > 3e-14:
        raise ArithmeticError('finite environment failed trace-preservation guard')
    weights.flags.writeable = schur.flags.writeable = False
    resources = {key: records[0]['source']['resources'][key] for key in native.COUNT_KEYS}
    if any(any(record['source']['resources'][key] != resources[key] for key in native.COUNT_KEYS)
           for record in records):
        raise ArithmeticError('phase branch resource schedules must match')
    all_resources = {key: order*value for key, value in resources.items()}
    metadata = {
        'schema': 'finite_source_number_dephasing_v1',
        'dimension': dimension, 'quadrature_order': order,
        'rate_per_second': rate, 'duration_seconds': duration,
        'branches': records, 'resources_per_realization': resources,
        'resources_weighted_mean_per_application': deepcopy(resources),
        'resources_maximum_per_application': deepcopy(resources),
        'resources_for_all_branch_contractions': all_resources,
        'independent_classical_environment_draw_per_application': True,
        'classical_environment_distribution_supplied': True,
        'environment_outcome_discarded_without_selection': True,
        'all_phase_unitaries_source_executed': True,
        'analytic_dephasing_schur_matrix_inserted': False,
        'gaussian_noise_law_or_rate_derived_from_observer': False,
        'ideal_quadrature_error': quadrature_bound(dimension, rate, duration, order),
        'numerical_diagnostics': {
            'rounded_weight_sum_before_normalization': raw_weight_sum,
            'rounded_weight_sum': math.fsum(weights),
            'maximum_diagonal_roundoff_projection_norm': max(projection_changes),
            'trace_preservation_maximum_entry_residual': float(np.max(np.abs(np.diag(schur)-1))),
            'outward_roundoff_certificate': False,
        },
        'canonical_compilation': canonical.policy(),
        'clock': 'supplied simulated Schrodinger-frame duration in seconds',
        'physical_service_error_certified': False,
    }
    return PhaseDiffusion(dimension, weights, tuple(unitaries), schur, metadata)


def phase_diffusion(rho, dimension, rate, duration, quadrature_order=DEFAULT_QUADRATURE_ORDER):
    """Apply the cached source-executed finite environment channel."""
    return channel(dimension, rate, duration, quadrature_order).apply(rho)
