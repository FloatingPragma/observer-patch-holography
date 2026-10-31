"""Source-operated full-Fock modular histories with retained read/abort outputs.

The convention is sigma_s(A)=rho**(is) A rho**(-is).  For rho=exp(-beta H)/Z,
its Heisenberg Hamiltonian parameter is t=-beta*s.  Beta, H, initialized bits
and source services are declared.  Neither s nor t is a laboratory duration.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from ..field_thermal.thermal import ThermalField
from ..full_fock_field.field import FieldGeometry, source_word
from ..full_fock_field.local_simulation import split_bound
from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import execute_complete


NUMERICAL_ALLOWANCE = 2e-9


def preparation_resources(word, bits):
    """Portable analytic budget for the existing nonnegative-state/phase compiler.

    The actual angle cost is retained. Only phases and rotations on tiny
    numerical amplitudes may vary between equivalent source preparations.
    No angle, source pulse, output or error is removed by this certificate.
    """
    if bits not in (4, 8):
        raise ValueError('four/eight preparation bits required by the field domain')
    dimension = 1 << bits
    amplitudes = ((dimension-1)*(dimension//2) if bits <= 6 else dimension-1)
    variable_expected = amplitudes+dimension
    amplitude_pulses = (3*bits*(dimension//2)*(dimension-1) if bits <= 6 else
                        3*((bits-1)*dimension+1))
    diagonal_pulses = 6+(3*bits-5)*dimension
    fixed_expected = F((2*bits-4)*dimension+4)
    fixed_expected += (F((dimension-1)*(4*bits-3)*(1 << (bits-2))) if bits <= 6 else
                       sum((F(1, 2)+2*k)*(1 << k) for k in range(bits)))
    if len(word) != amplitude_pulses+diagonal_pulses:
        raise ArithmeticError('preparation pulse schedule changed')
    fixed = F()
    variable, variable_cost = 0, 0.
    for event in word:
        kind, angle = event['kind'], event['angle']
        if kind in ('Z', 'I'):
            variable += 1
            variable_cost += abs(angle)
        else:
            coefficient = {'H': F(1, 2), 'J': F(1), 'X': F(1, 4), 'Y': F(1, 4)}[kind]
            if abs(angle) != float(coefficient)*math.pi:
                raise ArithmeticError('preparation fixed source scaffold changed')
            fixed += coefficient
    ceil_root = lambda n: math.isqrt(n)+(math.isqrt(n)**2 != n)
    amplitude_bound = (F(dimension-1, 2) if bits <= 6 else
                       F(sum(ceil_root(1 << k) for k in range(bits)), 2))
    variable_bound = amplitude_bound+ceil_root(dimension)
    if (variable != variable_expected or fixed != fixed_expected or
            variable_cost > float(variable_bound*F(22, 7))):
        raise ArithmeticError('preparation variable-angle schedule exceeds compiler budget')
    measured = math.fsum(abs(event['angle']) for event in word)
    upper = (fixed+variable_bound)*F(22, 7)
    if not math.isfinite(measured) or measured > float(upper):
        raise ArithmeticError('preparation processor norm exceeds compiler budget')
    return {'scheme': 'fixed-scaffold-plus-Walsh-Parseval-v1',
            'preparation_work_bits': bits,
            'fixed_scaffold_pi_coefficient': str(fixed),
            'variable_angle_sum_pi_coefficient_upper': str(variable_bound),
            'processor_norm_upper_rational': str(upper),
            'measured_preparation_processor_norm_upper': measured}


def validate_report_replay(retained, current):
    """Compare reports while certifying, rather than equating, variable prep cost.

    Complete error diagnostics and every pulse/service/memory count still
    compare normally. Each run reports its actual angle cost; its preparation
    cost must fit the same analytic certificate and add to the total cost.
    """
    variable_paths = {
        ('prepared_reference_readout', 'resources', 'integrated_processor_norm_upper'),
        ('prepared_reference_readout', 'preparation_resource_certificate',
         'measured_preparation_processor_norm_upper')}

    def compare(a, b, path=()):
        if type(a) is not type(b):
            raise ValueError(f'replay type mismatch at {path}')
        if isinstance(a, float):
            if not all(math.isfinite(x) for x in (a, b)):
                raise ValueError(f'nonfinite replay value at {path}')
            if any(isinstance(key, str) and 'probabilit' in key for key in path):
                if any(x < 0 or x > 1+NUMERICAL_ALLOWANCE for x in (a, b)):
                    raise ValueError(f'invalid replay probability at {path}')
        if isinstance(a, dict):
            if set(a) != set(b):
                raise ValueError(f'replay keys changed at {path}')
            for key in a:
                compare(a[key], b[key], path+(key,))
        elif isinstance(a, list):
            if len(a) != len(b):
                raise ValueError(f'replay length changed at {path}')
            for index, (x, y) in enumerate(zip(a, b, strict=True)):
                compare(x, y, path+(index,))
        elif path in variable_paths:
            if not all(math.isfinite(x) and x >= 0 for x in (a, b)):
                raise ValueError('nonfinite or negative measured preparation resource')
        elif isinstance(a, float):
            if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11):
                raise ValueError(f'replay numerical mismatch at {path}')
        elif a != b:
            raise ValueError(f'replay value changed at {path}')

    compare(retained, current)
    for result in (retained, current):
        readout = result['prepared_reference_readout']
        certificate = readout['preparation_resource_certificate']
        preparation = certificate['measured_preparation_processor_norm_upper']
        maximum = float(F(certificate['processor_norm_upper_rational']))
        fixed = float(F(certificate['fixed_scaffold_pi_coefficient']))*math.pi
        variable = float(F(certificate['variable_angle_sum_pi_coefficient_upper']))*math.pi
        total = readout['resources']['integrated_processor_norm_upper']
        continuation = result['full_algebra_certificate']['resources']['integrated_processor_norm_upper']
        if (preparation > maximum or preparation < fixed-NUMERICAL_ALLOWANCE or
                preparation-fixed > variable+NUMERICAL_ALLOWANCE or
                not math.isclose(total, preparation+continuation+2*math.pi,
                                 rel_tol=1e-12, abs_tol=1e-10)):
            raise ValueError('preparation resource budget or complete accounting failed')


def _rational(value, name):
    if type(value) not in (int, F):
        raise ValueError(f'exact rational {name} required')
    return F(value)


def _steps(value):
    if type(value) is not int or not 1 <= value <= 128:
        raise ValueError('one through 128 source subdivisions required')
    return value


def _complete_error(run, target):
    """Dilation operator-norm upper estimate, including every abort Kraus map.

    For a complete input basis, this controls arbitrary input/reference
    vectors, not only the prepared state. Floating values are diagnostics.
    """
    residual = float(np.linalg.norm(run.success-target, 2))
    failure = float(np.linalg.norm(run.failure_gram, 2))
    error = math.sqrt(residual**2+failure)
    if not all(math.isfinite(x) for x in (residual, failure, error)):
        raise ArithmeticError('nonfinite complete modular output')
    return error


@dataclass(frozen=True)
class ModularFlow:
    """Same prepared Gibbs reference, full algebra and native field action."""
    geometry: FieldGeometry = FieldGeometry((2, 1, 1))
    beta: F = F(1)

    def __post_init__(self):
        beta = _rational(self.beta, 'inverse temperature')
        if not 0 < beta <= 32:
            raise ValueError('inverse temperature must be in (0,32]')
        # Includes the existing support/conditioning and bounded preparation
        # admission. No eigenvalue clipping or number-sector deletion.
        ThermalField(self.geometry, float(beta)).target()
        object.__setattr__(self, 'beta', beta)

    @property
    def thermal(self):
        return ThermalField(self.geometry, float(self.beta))

    def duration(self, parameter):
        duration = -self.beta*_rational(parameter, 'dimensionless modular parameter')
        if abs(duration) > 4:
            raise ValueError('bounded Hamiltonian model duration has magnitude at most four')
        return duration

    def word(self, parameter, steps=4):
        return list(source_word(self.geometry, float(self.duration(parameter)), _steps(steps)))

    def bound(self, parameter, steps=4):
        return split_bound(self.geometry.sides, self.duration(parameter), _steps(steps),
                           F(self.geometry.speed), F(self.geometry.interaction), refined=True)

    def unitary(self, parameter):
        # This ideal reference is used only for validation, never as a source
        # primitive. The executed word consists entirely of native pulses.
        return expm(-1j*float(self.duration(parameter))*self.thermal.target()['hamiltonian'])

    def certificate(self, parameter=F(1, 4), steps=4):
        """Execute the complete input basis once, bounding every reference input."""
        dimension = 1 << self.geometry.modes
        word = self.word(parameter, steps)
        run = execute_complete(word, self.geometry.modes, np.eye(dimension, dtype=complex))
        error = _complete_error(run, self.unitary(parameter))
        bound = self.bound(parameter, steps)
        if error > float(bound)+NUMERICAL_ALLOWANCE:
            raise ArithmeticError('full-algebra modular source word exceeds its error budget')
        return run, {'dimensionless_modular_parameter': str(F(parameter)),
                     'Hamiltonian_model_duration': str(self.duration(parameter)),
                     'ideal_split_operator_bound': str(bound),
                     'complete_dilation_operator_error': error,
                     'abort_Gram_operator_norm': float(np.linalg.norm(run.failure_gram, 2)),
                     'resources': compiler.resources(word),
                     'all_input_columns_and_arbitrary_references_covered': True}

    def history(self, segments, *, read_mode=0):
        """Execute one to four increments, each followed by a complete number read.

        Pointer bits form a bounded in-register record tape. No read is
        selected, overwritten or normalized. All first export failures halt
        in retained abort slots; the source sender resets at every export.
        """
        modes, dimension = self.geometry.modes, 1 << self.geometry.modes
        if type(segments) not in (tuple, list) or not 1 <= len(segments) <= 4:
            raise ValueError('one through four modular increment/subdivision pairs required')
        if type(read_mode) is not int or not 0 <= read_mode < modes:
            raise ValueError('existing field read mode required')
        if any(type(pair) not in (tuple, list) or len(pair) != 2 for pair in segments):
            raise ValueError('each history segment is a modular increment/subdivision pair')
        word, ideal, bound = [], [np.eye(dimension, dtype=complex)], F()
        projectors = [np.diag([int(((i >> read_mode) & 1) == b) for i in range(dimension)])
                      for b in (0, 1)]
        for index, (parameter, steps) in enumerate(segments):
            word += self.word(parameter, steps)+compiler.cnot(read_mode, modes+index)
            u = self.unitary(parameter)
            ideal = [p@u@k for p in projectors for k in ideal]
            bound += self.bound(parameter, steps)
        bits = modes+len(segments)
        initial = np.zeros((1 << bits, dimension), complex)
        initial[:dimension] = np.eye(dimension)
        run = execute_complete(word, bits, initial)
        target = np.vstack(ideal)
        error = _complete_error(run, target)
        if error > float(bound)+NUMERICAL_ALLOWANCE:
            raise ArithmeticError('complete modular record history exceeds its error budget')
        return run, {'retained_binary_reads': len(segments), 'read_mode': read_mode,
                     'ideal_complete_history_bound': str(bound),
                     'complete_dilation_operator_error': error,
                     'resources': compiler.resources(word),
                     'initialized_record_bits': len(segments),
                     'history_outcomes': len(ideal), 'postselection_applied': False}

    def prepared_readout(self, parameter=F(1, 4), steps=4, *, read_mode=0):
        """Prepare this reference and retain its purifier during modular readout."""
        modes, dimension = self.geometry.modes, 1 << self.geometry.modes
        if type(read_mode) is not int or not 0 <= read_mode < modes:
            raise ValueError('existing field read mode required')
        preparation = self.thermal.word()
        preparation_certificate = preparation_resources(preparation, 2*modes)
        evolution = self.word(parameter, steps)
        word = preparation+evolution+compiler.cnot(read_mode, 2*modes)
        initial = np.zeros((2*dimension**2, 1), complex); initial[0, 0] = 1
        run = execute_complete(word, 2*modes+1, initial)
        resources = compiler.resources(word)
        measured = math.fsum(abs(angle) for _, angle, _ in run.word)
        if (resources['integrated_processor_norm_upper'] != measured or
                resources['native_pulses'] != len(run.word)):
            raise ArithmeticError('resource accounting differs from the executed source tape')
        target = self.thermal.target()
        columns = target['purification'].reshape(dimension, dimension).T
        columns = self.unitary(parameter)@columns
        expected = np.zeros_like(initial)
        evolved = columns.T.reshape(-1)
        for address, value in enumerate(evolved):
            expected[address+(((address >> read_mode) & 1)*dimension**2), 0] = value
        error = _complete_error(run, expected)
        bound = self.bound(parameter, steps)
        if error > float(bound)+2*NUMERICAL_ALLOWANCE:
            raise ArithmeticError('prepared modular observation exceeds its complete error budget')
        branches = run.success[:, 0].reshape(2, dimension, dimension)
        densities = [branch.T@branch.conj() for branch in branches]
        probability = [float(np.trace(rho).real) for rho in densities]
        rho = target['reference']
        return run, {'ideal_split_operator_bound': str(bound),
                     'complete_preparation_and_readout_error': error,
                     'pointer_probabilities': probability,
                     'first_export_abort_probability': run.failure_probability,
                     'complete_probability': sum(probability)+run.failure_probability,
                     'nonselective_read_reference_disturbance': float(np.linalg.norm(sum(densities)-rho, 'nuc')),
                     'resources': resources,
                     'preparation_resource_certificate': preparation_certificate,
                     'initialized_work_bits': 2*modes+1,
                     'retained_purifier_bits': modes, 'retained_pointer_bits': 1,
                     'postselection_or_success_normalization_applied': False}


def report():
    flow = ModularFlow()
    _, algebra = flow.certificate()
    _, history = flow.history(((F(1, 8), 4), (F(-1, 16), 4)))
    _, prepared = flow.prepared_readout()
    return {'schema': 'oph.complete_field_modular_flow.v1',
            'scope': 'FINITE_SOURCE_OPERATED_FULL_FOCK_MODULAR_HISTORIES',
            'geometry': list(flow.geometry.sides), 'full_Fock_dimension': 16,
            'declared_beta': str(flow.beta),
            'conventions': {'modular': 'sigma_s(A)=rho^(is) A rho^(-is)',
                            'Heisenberg': 'tau_t(A)=exp(itH) A exp(-itH)',
                            'common_parameter_map': 't=-beta*s'},
            'full_algebra_certificate': algebra, 'complete_read_history': history,
            'prepared_reference_readout': prepared,
            'physical_clock_calibration_or_regional_modular_covariance_claimed': False,
            'canonical_physical_premises_discharged': 0,
            'floating_residuals_are_outward_certificates': False}


if __name__ == '__main__':
    import json
    print(json.dumps(report(), indent=2, allow_nan=False))
