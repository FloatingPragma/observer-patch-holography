"""Same-reference regional comparison without imposing false covariance.

The global and regional modular automorphisms are compared by explicit
relative unitaries, compiled from the existing source services. The default
interacting dimer also admits an exact rational nonresonance certificate.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as sp

from . import intervals as q
from ..field_thermal.thermal import car_restriction, restrict
from ..full_fock_field.field import FieldGeometry, FullFock
from ..modular_flow.flow import ModularFlow, NUMERICAL_ALLOWANCE
from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import execute_complete


def exact_spectral_attachment():
    """Verify the spectrum and both reduced spectral weights from actual H.

    The characteristic polynomial identifies all eigenvalues/multiplicities.
    Nine partial-trace moment identities determine the nine spectral weights
    by an invertible Vandermonde system. No numerical eigenvector or assumed
    degeneracy enters the rational/algebraic check.
    """
    raw = FullFock(4).hamiltonian(FieldGeometry((2, 1, 1))).toarray()
    scaled = raw*16
    if not (np.array_equal(scaled.real, np.rint(scaled.real))
            and np.array_equal(scaled.imag, np.rint(scaled.imag))):
        raise ArithmeticError('actual source Hamiltonian changed its exact denominator')
    h = sp.Matrix(16, 16, lambda i, j:
                  (int(scaled[i, j].real)+sp.I*int(scaled[i, j].imag))/16)
    if h != h.conjugate().T:
        raise ArithmeticError('actual source Hamiltonian is not exactly Hermitian')
    x = sp.Symbol('x')
    t, u, d = sp.Rational(1, 16), sp.Integer(1), sp.sqrt(17)/4
    expected = x**4*(x*x-t*t)**2*((x-u)**2-t*t)**2*(x-u)*(x-2*u)*(x*x-u*x-4*t*t)
    if sp.expand(h.charpoly(x).as_expr()-expected) != 0:
        raise ArithmeticError('actual source characteristic polynomial differs')
    levels = [-t, t, 0, u-t, u+t, u, 2*u, (u-d)/2, (u+d)/2]
    half = sp.Rational(1, 2)
    weights = [[1, half, half, 0], [1, half, half, 0],
               [1, 3*half, 3*half, 0], [0, half, half, 1],
               [0, half, half, 1], [half, 0, 0, half], [0, 0, 0, 1]]
    weights += [[(1-sign*u/d)/4, (1+sign*u/d)/4,
                 (1+sign*u/d)/4, (1-sign*u/d)/4] for sign in (1, -1)]
    power = sp.eye(16)
    for degree in range(9):
        predicted = sp.diag(*[sp.simplify(sum(e**degree*w[a] for e, w in zip(levels, weights, strict=True)))
                              for a in range(4)])
        for site in (0, 1):
            # The signed CAR restriction has identical signs in each
            # nonzero matrix element here: H^k preserves total number and
            # the exterior address is equal on both sides of the trace.
            address = (lambda a, b: a+4*b) if site == 0 else (lambda a, b: b+4*a)
            reduced = sp.Matrix(4, 4, lambda i, j: sum(
                power[address(i, b), address(j, b)] for b in range(4)))
            if reduced != predicted:
                raise ArithmeticError(f'exact site {site} spectral moment {degree} differs')
        power *= h
    return {'actual_source_denominator': 16,
            'exact_source_Hermiticity_verified': True,
            'complete_Fock_dimension': 16,
            'distinct_global_energy_levels': 9,
            'characteristic_polynomial': str(sp.factor(expected)),
            'exact_partial_trace_moment_degrees': list(range(9)),
            'regional_site_moment_checks': 18,
            'regional_spectral_weights_proved_by_Vandermonde': True,
            'eigenvectors_or_floating_degeneracy_groups_used': False}


def frequency_certificate():
    """All nonzero local modular frequencies miss every global frequency.

    Fixed default dimer, beta=1 and mu=0. Inputs to the proof are exact
    rational/algebraic source coefficients. Exponentials/logarithms have
    rational outward enclosures from explicit positive series remainders.
    """
    attachment = exact_spectral_attachment()
    d, t = q.sqrt_integer(17)/4, F(1, 16)
    energies = [q.as_interval(-t), q.as_interval(t), q.as_interval(0),
                q.as_interval(1-t), q.as_interval(1+t), q.as_interval(1),
                q.as_interval(2), (1-d)/2, (1+d)/2]
    weights = [[1, F(1, 2), F(1, 2), 0], [1, F(1, 2), F(1, 2), 0],
               [1, F(3, 2), F(3, 2), 0], [0, F(1, 2), F(1, 2), 1],
               [0, F(1, 2), F(1, 2), 1], [F(1, 2), 0, 0, F(1, 2)],
               [0, 0, 0, 1]]
    weights += [[(1-sign/d)/4, (1+sign/d)/4, (1+sign/d)/4, (1-sign/d)/4]
                for sign in (1, -1)]
    raw = [sum((q.exp(-e)*w[a] for e, w in zip(energies, weights, strict=True)),
               q.Interval.point(0)) for a in range(4)]
    partition = sum(raw, q.Interval.point(0))
    probabilities = [p/partition for p in raw]
    # p_0 > p_up = p_down > p_double, certified without tolerance grouping.
    if not (raw[0].lo > raw[1].hi and raw[1].lo > raw[3].hi and raw[1] == raw[2]):
        raise ArithmeticError('regional spectral ordering not certified')
    frequencies = [q.log(raw[a]/raw[b]) for a, b in ((0, 1), (1, 3), (0, 3))]
    global_frequencies = [a-b for a in energies for b in energies]
    minimum = min(q.distance_from_zero(local-global_) for local in frequencies
                  for global_ in global_frequencies)
    if minimum <= F(1, 10000):
        raise ArithmeticError('global/regional modular nonresonance margin failed')
    return {'scope': 'EXACT_DEFAULT_DIMER_BETA1_MU0_COMMON_PARAMETER',
            'spectral_attachment': attachment,
            'regional_probability_intervals': [p.json() for p in probabilities],
            'positive_regional_frequency_intervals': [f.json() for f in frequencies],
            'global_energy_intervals': [e.json() for e in energies],
            'frequency_pair_checks': len(frequencies)*len(global_frequencies),
            'minimum_frequency_separation_lower': str(minimum),
            'simple_strict_separation_lower': '1/10000',
            'regional_centralizer_complex_dimension': 6,
            'every_covariant_linear_map_has_centralizer_range': True,
            'nontrivial_regional_clock_transmission_by_an_intertwiner': False,
            'other_references_parameters_or_clock_rescalings_excluded': False}


class ClockComparison:
    """Two one-site clocks and the global clock on one full-Fock field.

    C_AB(s)=W_A(s) W_B(s)^* compares explicitly implemented clock actions.
    It is an active correction, not a claim of passive regional covariance.
    """
    def __init__(self, geometry=FieldGeometry((2, 1, 1))):
        if not isinstance(geometry, FieldGeometry) or geometry.volume != 2:
            raise ValueError('one two-site full-Fock field required')
        self.flow = ModularFlow(geometry)
        self.geometry = geometry
        target = self.flow.thermal.target()
        self.reference = target['reference']
        self.hamiltonian = target['hamiltonian']
        self.local_references, self.local_generators = [], []
        for site in (0, 1):
            rho = restrict(self.reference, car_restriction(4, (2*site, 2*site+1)))
            diagonal = np.diag(rho).real
            if (np.min(diagonal) <= 0 or
                    np.linalg.norm(rho-np.diag(diagonal), 2) > 1e-12):
                raise ArithmeticError('local diagonal faithful reference certificate failed')
            # Scalar sector-independent gauge only; all local number coherences remain.
            generator = -np.log(diagonal)+math.log(diagonal[0])
            self.local_references.append(rho)
            self.local_generators.append(generator)

    @staticmethod
    def _parameter(s):
        if type(s) not in (int, F) or abs(s) > 4:
            raise ValueError('exact modular parameter in [-4,4] required')
        return F(s)

    @staticmethod
    def _label(label):
        if type(label) is not str or label not in ('global', 'site0', 'site1'):
            raise ValueError('global, site0 or site1 clock required')
        return label

    def generator(self, label):
        label = self._label(label)
        if label == 'global':
            return self.hamiltonian
        site = int(label[-1])
        return np.diag([self.local_generators[site][(mask >> (2*site)) & 3]
                        for mask in range(16)])

    def unitary(self, label, s):
        return expm(1j*float(self._parameter(s))*self.generator(label))

    def word(self, label, s, steps=8):
        label, s = self._label(label), self._parameter(s)
        if type(steps) is not int or not 1 <= steps <= 128:
            raise ValueError('one through 128 source subdivisions required')
        if label == 'global':
            return self.flow.word(s, steps)
        site = int(label[-1])
        local = compiler.compile_diagonal(float(s)*self.local_generators[site])
        return [{**event, 'targets': tuple(target+2*site for target in event['targets'])}
                for event in local]

    def comparison(self, target, origin, s):
        return self.unitary(target, s)@self.unitary(origin, s).conj().T

    def comparison_word(self, target, origin, s, steps=8):
        s = self._parameter(s)
        return self.word(origin, -s, steps)+self.word(target, s, steps)

    def certificate(self, target='site0', origin='global', s=F(1, 2), steps=8):
        word = self.comparison_word(target, origin, s, steps)
        run = execute_complete(word, 4, np.eye(16))
        expected = self.comparison(target, origin, s)
        bound = sum((self.flow.bound(s, steps) if label == 'global' else F())
                    for label in (origin, target))
        failure_norm = float(np.linalg.norm(run.failure_gram, 2))
        error = math.sqrt(float(np.linalg.norm(run.success-expected, 2))**2+failure_norm)
        if error > float(bound)+NUMERICAL_ALLOWANCE:
            raise ArithmeticError('complete clock comparison exceeds source error budget')
        rho = self.reference
        disturbance = float(np.linalg.norm(expected@rho@expected.conj().T-rho, 'nuc')/2)
        return run, {'target': target, 'origin': origin, 'modular_parameter': str(s),
                     'ideal_split_operator_bound': str(bound),
                     'complete_dilation_operator_error': error,
                     'failure_Gram_operator_norm': failure_norm,
                     'complete_input_Gram_residual': float(np.linalg.norm(
                         run.success.conj().T@run.success+run.failure_gram-np.eye(16), 2)),
                     'common_reference_half_trace_disturbance': disturbance,
                     'resources': compiler.resources(word),
                     'input_basis_columns': 16, 'all_number_sectors_retained': True,
                     'first_export_failure_slots_retained': len(word),
                     'postselection_or_reference_reset_applied': False}


def report():
    service = ClockComparison()
    records = [service.certificate(a, b)[1] for a, b in
               (('site0', 'global'), ('site1', 'global'), ('site1', 'site0'))]
    a, b = F(1, 3), F(-1, 5)
    c = service.comparison('site0', 'global', a)
    u = service.unitary('global', a)
    cocycle = c@u@service.comparison('site0', 'global', b)@u.conj().T
    cocycle_residual = float(np.linalg.norm(cocycle-service.comparison('site0', 'global', a+b), 2))
    composition = service.comparison('site1', 'site0', a)@service.comparison('site0', 'global', a)
    root = Path(__file__).resolve().parents[2]
    sources = ('Lean/QFT/ModularClockComparison.lean',
               'model/modular_clock_comparison/__main__.py',
               'model/modular_clock_comparison/comparison.py',
               'model/modular_clock_comparison/intervals.py',
               'model/modular_clock_comparison/contract.json',
               'model/full_fock_field/field.py', 'model/field_thermal/thermal.py',
               'model/modular_flow/flow.py', 'model/reference_recovery/recovery.py',
               'model/observer/source_scalar/compiler.py')
    return {'schema': 'oph.modular_clock_comparison.v1',
            'scope': 'FINITE_SOURCE_CLOCK_COCYCLE_AND_EXACT_REGIONAL_NONRESONANCE',
            'source_sha256': {path: hashlib.sha256((root/path).read_bytes()).hexdigest()
                              for path in sources},
            'law': 'same full-Fock Weyl/contact dimer, speed=1/8, U=1, beta=1, mu=0',
            'frequency_certificate': frequency_certificate(),
            'source_clock_comparisons': records,
            'cocycle_identity_numerical_residual': cocycle_residual,
            'regional_composition_numerical_residual': float(np.linalg.norm(
                composition-service.comparison('site1', 'global', a), 2)),
            'finite_premises_closed': ['actual full and reduced faithful reference',
                'complete source implementation of relative modular corrections',
                'arbitrary-input/reference complete-channel comparison',
                'transitive comparison and time-cocycle identities',
                'exact common-parameter modular nonresonance classification'],
            'physical_time_or_energy_calibration_established': False,
            'passive_regional_covariance_or_Einstein_premises_established': False,
            'canonical_physical_premises_discharged': 0}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', type=Path)
    args = parser.parse_args()
    result = report()
    text = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if args.write:
        args.write.write_text(text)
    else:
        print(text, end='')
