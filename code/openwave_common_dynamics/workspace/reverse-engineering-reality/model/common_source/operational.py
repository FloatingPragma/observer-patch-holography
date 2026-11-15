"""GNS and modular readouts of the executed common-source observation state.

The matrix algebra, Born pairing and subsystem access are declared inputs.
This constructs their observer representation; it does not derive them from
integer repair. No thermal state, identity noise or clock calibration is
inserted. A modular parameter is never identified with model-coordinate time.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from ..observer.source_scalar.field import quadratures
from . import experiment, scalar_coupling


RANK_TOLERANCE = 1e-12


def density_spectrum(density):
    """Validate a state without renormalizing it or filling its kernel.

    The rank tolerance is a numerical diagnostic, not a certified eigenvalue
    enclosure. Tiny negative eigenvalues from floating contraction are treated
    as zero only in the quotient; the input density is retained unchanged.
    """
    rho = np.asarray(density, dtype=complex)
    if (rho.ndim != 2 or rho.shape[0] != rho.shape[1] or not rho.size
            or not np.all(np.isfinite(rho))):
        raise ValueError('a finite nonempty square density matrix required')
    if np.linalg.norm(rho-rho.conj().T) > 1e-10 or abs(np.trace(rho)-1) > 1e-9:
        raise ValueError('a Hermitian trace-one state required')
    values, vectors = np.linalg.eigh((rho+rho.conj().T)/2)
    if values[0] < -RANK_TOLERANCE:
        raise ValueError('a positive density matrix required')
    return rho, values, vectors


@dataclass(frozen=True)
class MatrixGNS:
    """M_d / {A: Tr(rho A* A)=0}, realized as d by rank(rho) matrices.

    [A] maps to A V sqrt(p); pi(B) acts by left multiplication. In particular
    a rank-r state gives dimension d*r, not r*r. Restricting the algebra to
    its support corner would be a different observer.
    """
    density: np.ndarray
    root_columns: np.ndarray

    @classmethod
    def from_density(cls, density):
        rho, values, vectors = density_spectrum(density)
        support = values > RANK_TOLERANCE
        return cls(rho.copy(), vectors[:, support]*np.sqrt(values[support]))

    @property
    def dimension(self):
        return int(self.root_columns.size)

    @property
    def rank(self):
        return self.root_columns.shape[1]

    def lift(self, observable):
        return np.asarray(observable)@self.root_columns

    def expectation(self, observable):
        return np.vdot(self.root_columns, self.lift(observable))

    def inner(self, first, second):
        return np.vdot(self.lift(first), self.lift(second))


def modular_hamiltonian(density):
    """K=-log(rho) on the declared full algebra; reject nonfaithful states."""
    _, values, vectors = density_spectrum(density)
    if values[0] <= RANK_TOLERANCE:
        raise ValueError('full-algebra modular flow requires a faithful state')
    return (vectors*(-np.log(values)))@vectors.conj().T


def modular_flow(density, observable, parameter):
    """sigma_s(A)=rho^(is) A rho^(-is), for real modular parameter s."""
    if isinstance(parameter, bool) or not np.isreal(parameter) or not math.isfinite(parameter):
        raise ValueError('a finite real modular parameter required')
    _, values, vectors = density_spectrum(density)
    if values[0] <= RANK_TOLERANCE:
        raise ValueError('full-algebra modular flow requires a faithful state')
    unitary = (vectors*np.exp(1j*float(parameter)*np.log(values)))@vectors.conj().T
    return unitary@np.asarray(observable)@unitary.conj().T


def terminal_vectors(result):
    """Decode all four unnormalized states; never select a successful branch."""
    vectors = {}
    size = 4*result['dimension']**2
    for branch in result['branches']:
        record = tuple(branch['record'])
        values = np.asarray(branch['unnormalized_terminal_state'])
        if (record in vectors or values.shape != (size, 2)
                or values.dtype.kind not in 'fiuc' or not np.all(np.isfinite(values))
                or not math.isfinite(branch['probability'])):
            raise ValueError('distinct complete records and matching working states required')
        vector = values[:, 0]+1j*values[:, 1]
        if abs(np.vdot(vector, vector).real-branch['probability']) > 1e-10:
            raise ValueError('record weight differs from its working-state norm')
        vectors[record] = vector
    if set(vectors) != {(0, 0), (0, 1), (1, 0), (1, 1)}:
        raise ValueError('the complete two-bit record law required')
    if abs(sum(np.vdot(v, v).real for v in vectors.values())-1) > 1e-9:
        raise ValueError('the complete observation state must be normalized')
    return vectors


def scalar_reduction(d, density, site=0):
    """Trace the actual gauge/link and other scalar factor: |n1,n0,gauge>."""
    if type(d) is not int or d not in (2, 4) or type(site) is not int or site not in (0, 1):
        raise ValueError('executed scalar dimension 2 or 4 and site 0 or 1 required')
    matrix = np.asarray(density).reshape(d, d, 4, d, d, 4)
    return (np.einsum('abcadc->bd', matrix) if site == 0
            else np.einsum('abcdbc->ad', matrix))


def scalar_state(result, site=0, record=None):
    """An explicit observer access choice: keep a record or forget all records.

    Conditional states are normalized only after retaining the outcome weight.
    Forgetting records is a stated restriction of access, not postselection.
    """
    vectors = terminal_vectors(result)
    if record is not None:
        record = tuple(record)
        if record not in vectors:
            raise ValueError('an actual complete two-bit record required')
        selected = [vectors[record]]
    else:
        selected = list(vectors.values())
    weight = float(sum(np.vdot(v, v).real for v in selected))
    if weight <= RANK_TOLERANCE:
        raise ValueError('cannot condition on a zero-weight record')
    rho = sum(scalar_reduction(result['dimension'], np.outer(v, v.conj()), site)
              for v in selected)
    return (rho/weight if record is not None else rho), weight


def regional_gns_embedding(result, site=0):
    """Isometry from regional GNS space into the common record/work GNS space.

    J[A]=[embed(A)] follows from equality of the state pairings. This binds
    both regional observers to one representation and one cyclic history
    vector, rather than matching independently chosen probability tables.
    Coordinates use row-major d by d Hilbert--Schmidt matrices locally and
    the concatenated four actual branch vectors globally.
    """
    from .spatial import scalar_operator
    rho, _ = scalar_state(result, site)
    local = MatrixGNS.from_density(rho)
    d = result['dimension']
    if local.rank != d:
        raise ValueError('this coordinate implementation requires a faithful regional state')
    inverse_root = np.linalg.inv(local.root_columns)
    vectors = terminal_vectors(result)
    columns = []
    for index in range(d*d):
        coordinate = np.zeros((d, d))
        coordinate.flat[index] = 1
        observable = coordinate@inverse_root
        embedded = scalar_operator(d, site, observable)
        columns.append(np.concatenate([embedded@v for v in vectors.values()]))
    return np.column_stack(columns), local, np.concatenate(list(vectors.values()))


def record_aware_scalar_view(result):
    """Keep BOTH committed bits and access scalar site 1 in each record block.

    For this executed state those reduced blocks are numerically faithful,
    unlike the just-measured scalar-site-0 blocks. Modular automorphisms act
    within each block and fix the central record projectors exactly. Thus a
    nontrivial finite modular observer need not discard its classical memory.
    """
    d = result['dimension']
    number = np.diag(np.arange(d))
    rows = []
    for record in sorted(terminal_vectors(result)):
        rho, probability = scalar_state(result, site=1, record=record)
        gns = MatrixGNS.from_density(rho)
        moved = modular_flow(rho, number, .3)
        rows.append({'record': list(record), 'probability': probability,
                     'rank': gns.rank,
                     'density_eigenvalues': [float(x) for x in np.linalg.eigvalsh(rho)],
                     'state_preservation_error': float(abs(np.trace(rho@moved)-np.trace(rho@number))),
                     'number_probe_moves': bool(np.linalg.norm(moved-number) > 1e-3)})
    if not all(row['rank'] == d and row['number_probe_moves']
               and row['state_preservation_error'] < 1e-10 for row in rows):
        raise AssertionError('record-aware scalar1 observer failed its declared finite checks')
    return {'access': 'both records and scalar-site-1 algebra within each record block',
            'algebra': f'direct sum of four M_{d} blocks',
            'gns_hilbert_dimension': 4*d*d,
            'modular_flow_fixes_committed_records': True,
            'blocks': rows,
            'physical_clock_calibration_established': False}


def observation_receipt(result):
    """Bind GNS, modular flow and source-time comparison to this same history."""
    d = result['dimension']
    vectors = terminal_vectors(result)
    rho, weight = scalar_state(result)
    gns = MatrixGNS.from_density(rho)
    k = modular_hamiltonian(rho)
    eigenvalues = np.linalg.eigvalsh(rho)
    number = np.diag(np.arange(d))
    q = quadratures(d)[0]
    q_values, q_vectors = np.linalg.eigh(q)
    plus = (q_vectors*(q_values > 0))@q_vectors.conj().T
    s = .7
    moved = modular_flow(rho, number, s)
    probe = number-np.trace(rho@number)*np.eye(d)
    correlation = np.trace(rho@probe@modular_flow(rho, probe, s))
    work = sum(np.outer(v, v.conj()) for v in vectors.values())
    h = scalar_coupling.ideal_hamiltonian(d, result['coupling'])
    rate = scalar_reduction(d, -1j*(h@work-work@h))
    conditional_ranks = []
    for record in sorted(vectors):
        conditional, probability = scalar_state(result, record=record)
        conditional_ranks.append({'record': list(record), 'probability': probability,
                                  'rank': MatrixGNS.from_density(conditional).rank})
    probabilities = [float(np.vdot(v, v).real) for v in vectors.values()]
    gram_error = abs(gns.inner(q, number)-np.trace(rho@q.conj().T@number))
    plus_error = abs(float(np.trace(rho@plus).real)-result['second_record_marginal'][1])
    inclusion, local, common_cyclic = regional_gns_embedding(result)
    from .spatial import scalar_operator
    embedded_q = scalar_operator(d, 0, q)
    common_action = np.concatenate([embedded_q@block for block in
        inclusion.reshape(4, 4*d*d, d*d)], axis=0)
    isometry_error = float(np.linalg.norm(inclusion.conj().T@inclusion-np.eye(d*d)))
    cyclic_error = float(np.linalg.norm(inclusion@local.root_columns.ravel()-common_cyclic))
    intertwiner_error = float(np.linalg.norm(common_action-inclusion@np.kron(q, np.eye(d))))
    checks = {
        'gns_reproduces_observation_pairing': gram_error < 1e-10,
        'regional_gns_embeds_isometrically_into_same_history': isometry_error < 1e-9,
        'same_cyclic_history_and_regional_operator_action': max(cyclic_error, intertwiner_error) < 1e-9,
        'same_history_detector_marginal': plus_error < 1e-10,
        'unconditional_scalar_state_numerically_faithful': eigenvalues[0] > 1e-8,
        'modular_flow_nontrivial_on_number_probe': np.linalg.norm(moved-number) > 1e-3,
        'modular_flow_preserves_state': abs(np.trace(rho@moved)-np.trace(rho@number)) < 1e-10,
        'terminal_sign_effect_fixed_by_modular_flow':
            np.linalg.norm(modular_flow(rho, plus, s)-plus) < 1e-9,
        'source_time_is_not_this_stationary_modular_flow': np.linalg.norm(rate) > 1e-3,
        'conditional_sharp_read_states_are_not_full_rank':
            all(row['rank'] <= d//2 for row in conditional_ranks),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError(checks)
    return {
        'dimension': d, 'steps_per_segment': result['steps_per_segment'],
        'native_current': result['native_current'],
        'source': 'complete terminal vectors of the common interacting source',
        'access': 'site-0 scalar algebra, forgetting both classical records',
        'state_trace': weight,
        'density': [[[float(z.real), float(z.imag)] for z in row] for row in rho],
        'density_eigenvalues': [float(x) for x in eigenvalues],
        'modular_hamiltonian_eigenvalues': [float(x) for x in np.linalg.eigvalsh(k)],
        'gns_hilbert_dimension': gns.dimension,
        'working_hilbert_dimension': 4*d*d,
        'record_and_work_block_algebra_gns_dimension': sum(4*d*d for p in probabilities if p > RANK_TOLERANCE),
        'records_only_gns_dimension': sum(p > RANK_TOLERANCE for p in probabilities),
        'records_only_modular_flow': 'identity',
        'conditional_scalar_ranks': conditional_ranks,
        'record_aware_scalar1_observer': record_aware_scalar_view(result),
        'modular_parameter': s,
        'modular_number_correlation': [float(correlation.real), float(correlation.imag)],
        'number_observable_modular_displacement_frobenius': float(np.linalg.norm(moved-number)),
        'source_model_time_local_state_derivative_frobenius': float(np.linalg.norm(rate)),
        'gns_pairing_error': float(gram_error), 'detector_marginal_error': float(plus_error),
        'regional_gns_isometry_error': isometry_error,
        'regional_gns_cyclic_vector_error': cyclic_error,
        'regional_gns_operator_intertwiner_error': intertwiner_error,
        'rank_tolerance': RANK_TOLERANCE,
        'numerical_spectra_have_outward_error_certificate': False,
        'clock_calibration_established': False,
        'quantum_algebra_and_born_pairing_derived_from_native_repair': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }


def run(d=4, steps=512):
    return observation_receipt(experiment.history(d, steps))
