"""Prepare a full-Fock Gibbs marginal using the existing proper-code source.

The finite field action, inverse temperature and initialized work bits are
declared. The complete purification stays present. Faithfulness and stationarity
refer to the field algebra, not the joint field/purifier algebra or every read.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from ..full_fock_field import field
from ..observer.source_scalar import compiler
from ..reference_observation.observation import ReferenceObservation, frozen
from ..interacting_field.lattice import real


NUMERICAL_LIMIT = 2e-9


def execute_with_exports(word, bits, initial):
    """Retain the valid output and every first failed proper-code export.

    Failure terminates this finite program and resets the source sender.
    Each event contains separate Kraus vectors for the complement basis
    values; its coarse abort density is their sum of outer products, never
    the outer product of their coherent sum. A vector contains only the
    untouched register/reference axes; the common failure flag is implicit.
    """
    if type(bits) is not int or not 1 <= bits <= 9:
        raise ValueError('one through nine bounded work bits required')
    if np.ma.isMaskedArray(initial):
        raise ValueError('masked state forbidden')
    state = np.asarray(initial)
    if (state.ndim not in (1, 2) or state.shape[0] != 1 << bits
            or state.size > 2**18 or state.dtype.kind not in 'iufc'
            or not np.isfinite(state).all() or not state.size):
        raise ValueError('finite bounded complete work/reference columns required')
    vector = state.ndim == 1
    state = np.array(state, complex, copy=True).reshape(1 << bits, -1)
    word = compiler._validate(word, bits)
    retained_entries = sum((6-(1 << len(e['targets'])))*(1 << (bits-len(e['targets'])))*state.shape[1]
                           for e in word)
    if retained_entries > 2**24:
        raise ValueError('complete failure outputs exceed the bounded host storage domain')
    failures = []
    for event in word:
        targets = event['targets']
        dimension = 1 << len(targets)
        full = compiler.native_pulse(event['kind'], event['angle'])
        if (np.ma.isMaskedArray(full) or np.asarray(full).shape != (6, 6)
                or np.asarray(full).dtype.kind not in 'iufc' or not np.isfinite(full).all()
                or np.linalg.norm(full.conj().T@full-np.eye(6), 2) > 2e-11):
            raise ArithmeticError('source primitive is not an admitted finite full unitary')
        axes = [bits-1-j for j in targets]
        axes += [axis for axis in range(bits) if axis not in axes]
        axes += [bits]
        arranged = state.reshape([2]*bits+[state.shape[1]]).transpose(axes)
        failed = (full[dimension:, :dimension]@arranged.reshape(dimension, -1))
        failed = failed.reshape(6-dimension, -1, state.shape[1])
        failures.append(tuple(frozen(row[:, 0] if vector else row) for row in failed))
        state = compiler._apply(state, full[:dimension, :dimension], targets, bits)
    return frozen(state[:, 0] if vector else state), tuple(failures)


def failure_probability(events):
    """Sum the coarse first-failure event weights over their Kraus branches."""
    return math.fsum(float(np.vdot(branch, branch).real) for event in events for branch in event)


def _matrix(value, dimension):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked matrices are not admitted')
    raw = np.asarray(value)
    if (raw.shape != (dimension, dimension) or raw.dtype.kind not in 'iufc'
            or not np.isfinite(raw).all()):
        raise ValueError('finite matrix on the complete declared carrier required')
    return np.asarray(raw, complex)


def entropy(density):
    """Entropy of an admitted faithful density; no spectral clipping."""
    raw = np.asarray(density)
    if raw.ndim != 2 or raw.shape[0] != raw.shape[1] or not 1 <= len(raw) <= 256:
        raise ValueError('bounded square density required')
    rho = _matrix(density, len(raw))
    if np.linalg.norm(rho-rho.conj().T, 2) > 1e-12 or abs(np.trace(rho)-1) > 1e-10:
        raise ValueError('Hermitian unit-trace density required')
    values = np.linalg.eigvalsh(rho)
    if values[0] <= 0:
        raise ValueError('faithful density required without support deletion')
    return float(-values@np.log(values))


def car_restriction(modes, kept):
    """Complete signed Fock restriction, including every occupation number.

    F(A direct_sum B)=F(A) tensor F(B), with kept modes ordered before
    the environment. The wedge permutation sign is retained on all inputs.
    """
    if type(modes) is not int or not 1 <= modes <= 8:
        raise ValueError('one through eight CAR modes required')
    if (type(kept) is not tuple or not kept
            or any(type(i) is not int or not 0 <= i < modes for i in kept)
            or tuple(sorted(set(kept))) != kept):
        raise ValueError('nonempty increasing tuple of distinct kept modes required')
    environment = tuple(i for i in range(modes) if i not in kept)
    order = kept+environment
    position = {mode: i for i, mode in enumerate(order)}
    kraus = [np.zeros((1 << len(kept), 1 << modes), complex)
             for _ in range(1 << len(environment))]
    for mask in range(1 << modes):
        occupied = [i for i in range(modes) if mask >> i & 1]
        inversions = sum(position[i] > position[j]
                         for index, i in enumerate(occupied) for j in occupied[index+1:])
        local = sum(((mask >> i) & 1) << j for j, i in enumerate(kept))
        outside = sum(((mask >> i) & 1) << j for j, i in enumerate(environment))
        kraus[outside][local, mask] = (-1)**inversions
    return tuple(frozen(k) for k in kraus)


def restrict(density, kraus):
    """Apply the complete restriction, without conditioning on an outcome."""
    if type(kraus) not in (tuple, list) or not kraus:
        raise ValueError('nonempty complete Kraus family required')
    first = np.asarray(kraus[0])
    if first.ndim != 2 or not 1 <= first.shape[1] <= 256:
        raise ValueError('bounded Kraus family required')
    if any(np.ma.isMaskedArray(k) for k in kraus):
        raise ValueError('masked Kraus operators are not admitted')
    matrices = tuple(np.asarray(k) for k in kraus)
    if any(k.shape != first.shape or k.dtype.kind not in 'iufc' or
           not np.isfinite(k).all() for k in matrices):
        raise ValueError('finite common-shaped Kraus operators required')
    if np.linalg.norm(sum(k.conj().T@k for k in matrices)-np.eye(first.shape[1]), 2) > 1e-12:
        raise ValueError('trace-preserving Kraus family required')
    rho = _matrix(density, first.shape[1])
    return sum(k@rho@k.conj().T for k in matrices)


@dataclass(frozen=True)
class ThermalField:
    """One or two full two-component sites, with a declared beta>0.

    Four field modes plus four purifier bits fit the existing complex
    eight-qubit preparation compiler. This is a bounded execution domain,
    not a particle-number cutoff or a thermodynamic-limit preparation.
    """
    geometry: field.FieldGeometry = field.FieldGeometry((2, 1, 1))
    beta: float = 1.

    def __post_init__(self):
        if not isinstance(self.geometry, field.FieldGeometry) or self.geometry.volume not in (1, 2):
            raise ValueError('one or two complete field sites required by the preparation domain')
        beta = real(self.beta, 'declared inverse temperature')
        if beta <= 0:
            raise ValueError('strictly positive inverse temperature required')
        object.__setattr__(self, 'beta', beta)

    @property
    def dimension(self):
        return 1 << self.geometry.modes

    def target(self):
        """Gibbs density and canonical square-root purification in field-low order."""
        h = field.FullFock(self.geometry.modes).hamiltonian(self.geometry).toarray()
        energies, vectors = np.linalg.eigh(h)
        exponents = -self.beta*(energies-energies[0])
        if not np.isfinite(exponents).all() or np.min(exponents) < -600:
            raise ValueError('temperature exceeds the faithful floating preparation domain')
        raw = np.exp(exponents)
        weights = raw/np.sum(raw)
        root = (vectors*np.sqrt(weights))@vectors.conj().T
        rho = (vectors*weights)@vectors.conj().T
        log_partition = math.log(float(np.sum(raw)))-self.beta*float(energies[0])
        state = root.T.reshape(-1).copy()
        # Ill-conditioned numerical densities must fail closed, even when all
        # formal Gibbs weights are positive. No eigenvalue is clipped away.
        minimum = float(np.linalg.eigvalsh(rho)[0])
        if (minimum <= 100*NUMERICAL_LIMIT or abs(np.linalg.norm(state)-1) > 2e-12
                or abs(np.trace(rho)-1) > 2e-12):
            raise ValueError('faithful reference lacks the declared numerical margin')
        return {'hamiltonian': frozen(h), 'reference': frozen(rho),
                'purification': frozen(state), 'minimum_weight': minimum,
                'log_partition': log_partition}

    def word(self):
        target = self.target()['purification']
        return compiler.compile_real_state(np.abs(target))+compiler.compile_diagonal(np.angle(target))

    def prepare(self):
        """Execute every source pulse on initialized work bits; retain the purifier."""
        target = self.target()
        word = self.word()
        blank = np.zeros(self.dimension**2, complex)
        blank[0] = 1
        state, failures = execute_with_exports(word, 2*self.geometry.modes, blank)
        error = math.sqrt(float(np.linalg.norm(state-target['purification']))**2+
                          failure_probability(failures))
        if not math.isfinite(error) or error > NUMERICAL_LIMIT:
            raise ArithmeticError('source-prepared purification exceeds its numerical admission bound')
        columns = state.reshape(self.dimension, self.dimension).T
        rho = columns@columns.conj().T
        return {'state': frozen(state), 'success_field_density': frozen(rho),
                'first_export_failures': failures,
                'purification_error': error, 'resources': {
                    **compiler.resources(word),
                    'initialized_work_bits_required': 2*self.geometry.modes,
                    'field_modes': self.geometry.modes,
                    'retained_purifier_bits': self.geometry.modes,
                    'retained_first_export_failure_slots': len(failures),
                    'retained_export_complement_Kraus_branches': sum(map(len, failures)),
                    'purifier_discard_or_postselection_applied': False,
                    'physical_initialization_cost_calibrated': False}}

    def report(self):
        target, prepared = self.target(), self.prepare()
        rho, h = target['reference'], target['hamiltonian']
        actual = prepared['success_field_density']
        size, modes = self.dimension, self.geometry.modes
        columns = prepared['state'].reshape(size, size).T
        space = field.FullFock(modes)
        duration, steps = 1/16, 4
        exact = space.exact_evolve(self.geometry, columns, duration)
        continuation = list(field.source_word(self.geometry, duration, steps))
        continuation += compiler.cnot(0, 2*modes)
        continued_input = np.zeros(2*size**2, complex)
        continued_input[:size**2] = prepared['state']
        continued, continuation_failures = execute_with_exports(continuation, 2*modes+1, continued_input)
        pointer = tuple(continued[bit*size**2:(bit+1)*size**2].reshape(size, size).T
                        for bit in (0, 1))
        effects = (np.eye(size)-space.number(0).toarray(), space.number(0).toarray())
        pointer_error = math.sqrt(math.fsum(float(np.linalg.norm(x-p@exact))**2
                                          for x, p in zip(pointer, effects, strict=True))+
                                  failure_probability(continuation_failures))
        failed_probability = failure_probability(prepared['first_export_failures']+continuation_failures)
        target_vector = np.zeros(size, complex); target_vector[3] = 1
        direction = np.outer(target_vector, target_vector.conj())-rho
        epsilon = 1/64
        varied = rho+epsilon*direction
        modular = self.beta*h+target['log_partition']*np.eye(size)
        delta_energy = float(np.trace((varied-rho)@h).real)
        delta_entropy = entropy(varied)-entropy(rho)
        relative = float(np.trace(varied@modular).real)-entropy(varied)
        remainder = (epsilon**2*float(np.vdot(direction, direction).real)/
                     (2*(1-epsilon)*target['minimum_weight']))
        n0, n1 = (space.number(i).toarray() for i in (0, 1))
        a0, a1 = (space.annihilator(i).toarray() for i in (0, 1))
        occupations = [float(np.trace(rho@n).real) for n in (n0, n1)]
        coherence = np.trace(rho@a0.conj().T@a1)
        wick_defect = float(np.trace(rho@n0@n1).real)-(occupations[0]*occupations[1]-abs(coherence)**2)
        keep = (0, modes-1)
        first = car_restriction(modes, keep)
        middle = restrict(rho, first)
        second = car_restriction(2, (1,))
        direct = car_restriction(modes, (modes-1,))
        a = ReferenceObservation(rho, first)
        b = ReferenceObservation(middle, second)
        ab = ReferenceObservation(rho, direct)
        # A complete observable basis is tested, not only the reference itself.
        tower = 0.
        for i in range(size):
            for j in range(size):
                unit = np.zeros((size, size), complex); unit[i, j] = 1
                tower = max(tower, float(np.linalg.norm(b.observe(a.observe(unit))-ab.observe(unit), 2)))
        budget = (self.geometry.speed*self.geometry.spatial_edges+
                  self.geometry.interaction*self.geometry.volume)**3*abs(duration)**3/(3*steps**2)
        checks = {'source_reference_trace_norm_error': float(np.linalg.norm(actual-rho, 'nuc')),
                  'source_purification_vector_error': prepared['purification_error'],
                  'complete_pointer_and_abort_probability_sum': float(sum(np.vdot(x, x).real for x in pointer))+failed_probability,
                  'complete_pointer_source_error': pointer_error,
                  'retained_export_failure_probability': failed_probability,
                  'ideal_source_split_operator_bound': budget,
                  'exact_evolution_field_reference_drift': float(np.linalg.norm(exact@exact.conj().T-actual, 'nuc')),
                  'noncommuting_variation_norm': float(np.linalg.norm(direction@rho-rho@direction)),
                  'first_law_identity_residual': abs(delta_entropy+relative-self.beta*delta_energy),
                  'quantum_relative_entropy': relative,
                  'relative_entropy_remainder_upper': remainder,
                  'onsite_non_Gaussian_Wick_defect': wick_defect,
                  'common_reference_observation_tower_max_error': tower}
        if (not all(math.isfinite(value) for value in checks.values())
                or checks['source_reference_trace_norm_error'] > 3*NUMERICAL_LIMIT
                or abs(checks['complete_pointer_and_abort_probability_sum']-1) > 10*NUMERICAL_LIMIT
                or checks['complete_pointer_source_error'] > budget+10*NUMERICAL_LIMIT
                or checks['first_law_identity_residual'] > NUMERICAL_LIMIT
                or not -NUMERICAL_LIMIT <= relative <= remainder+NUMERICAL_LIMIT
                or tower > NUMERICAL_LIMIT):
            raise ArithmeticError('same-reference field checks exceed their admission bounds')
        return {'field_modes': modes, 'full_Fock_dimension': size,
                'geometry': list(self.geometry.sides), 'speed': self.geometry.speed,
                'interaction': self.geometry.interaction, 'declared_beta': self.beta,
                'reference_minimum_weight': target['minimum_weight'],
                'number_sector_projection_applied': False,
                'preparation_resources': prepared['resources'],
                'continuation_resources': {**compiler.resources(continuation),
                    'additional_initialized_pointer_bits_required': 1,
                    'retained_pointer_outcomes': 2,
                    'retained_first_export_failure_slots': len(continuation_failures),
                    'retained_export_complement_Kraus_branches': sum(map(len, continuation_failures))},
                'checks': checks,
                'scope': {'source_preparation_executed': True,
                          'continuation_error_starts_from_same_executed_success_block': True,
                          'purifier_retained': True,
                          'joint_purification_faithful_on_full_algebra': False,
                          'reference_observation_maps_source_compiled': False,
                          'native_reference_selected': False,
                          'physical_temperature_or_clock_calibrated': False,
                          'geometric_G1_G6_or_Einstein_premises_discharged': False,
                          'floating_error_outward_certified': False}}


if __name__ == '__main__':
    import json
    print(json.dumps(ThermalField().report(), indent=2, sort_keys=True))
