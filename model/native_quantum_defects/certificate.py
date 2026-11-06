"""Complete native spectrum, compact-support search and recorded-law controls.

The matrix is counted from native repair. Its modes are reconstructed history
fluctuations, not supplied coherent controls or identified physical particles.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp

from tools import strict_json
from ..native_quantum_reconstruction import certificate as reconstruction
from ..observer import native, quantum_mixing
from ..observer.source_discrimination import circuit
from ..repair_scheduler import coins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
STATES = native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}
PAIRS = tuple(tuple(i for i, value in enumerate(state) if value) for state in STATES)
PRIME = 65521


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_core(core):
    require(type(core) is tuple and all(type(i) is int and 0 <= i < 12 for i in core)
            and core == tuple(sorted(set(core))), 'a sorted distinct native port tuple is required')
    return core


@lru_cache(maxsize=1)
def counts():
    matrix, denominator = reconstruction.source_counts()
    require(denominator == 60 and len(matrix) == 66, 'native source dimensions changed')
    return np.array(matrix, dtype=np.int64)


@lru_cache(maxsize=1)
def transfer_spectrum():
    """Exact characteristic factors, independently checked by Newton traces."""
    x = sp.Symbol('x')
    factors = [([1, -60], 1), ([1, -54], 5), ([1, -52], 5),
               ([1, -50], 11), ([1, -48], 5), ([1, -110, 3020], 3),
               ([1, -100, 2488], 1), ([1, -100, 2490], 4),
               ([1, -100, 2498], 4), ([1, -150, 7480, -124024], 5)]
    polynomial = sp.Poly(sp.prod(sp.Poly.from_list(a, x).as_expr()**m for a, m in factors), x)
    exact = sp.Matrix(counts().tolist()).charpoly(x).as_poly()
    require(exact == polynomial, 'native characteristic polynomial changed')
    matrix = counts().astype(object)
    power = np.eye(66, dtype=object)
    traces = []
    coefficients = [1]
    for k in range(1, 67):
        power = power @ matrix
        traces.append(int(np.trace(power)))
        numerator = -sum(coefficients[k-i]*traces[i-1] for i in range(1, k+1))
        require(numerator % k == 0, 'Newton trace coefficient is not integral')
        coefficients.append(numerator//k)
    require(coefficients == [int(a) for a in polynomial.all_coeffs()],
            'independent native trace reconstruction differs')
    intervals = []
    for coefficients_, multiplicity in factors:
        for (lo, hi), root_multiplicity in sp.Poly.from_list(coefficients_, x).intervals(eps=sp.Rational(1, 10000)):
            require(root_multiplicity == 1, 'factor contains a repeated root')
            if lo != 60:
                require(40 < lo <= hi < sp.Rational(229, 4), 'nonconstant spectral enclosure changed')
            intervals.append({'integer_transfer_interval': [str(lo), str(hi)],
                              'multiplicity': multiplicity})
    intervals.sort(key=lambda row: Q(row['integer_transfer_interval'][0]))
    require(sum(row['multiplicity'] for row in intervals) == 66, 'spectrum lost modes')
    return {'characteristic_factors': [{'coefficients': a, 'multiplicity': m} for a, m in factors],
            'all_66_Newton_trace_coefficients_checked': True,
            'eigenvalue_intervals': intervals, 'distinct_eigenvalues': len(intervals),
            'transfer_denominator': 60, 'stationary_eigenvalue_multiplicity': 1,
            'largest_nonconstant_transfer_eigenvalue': '(55+sqrt(5))/60',
            'rational_nonconstant_contraction_upper': '229/240',
            'reconstructed_energy_rule': '-log(integer_transfer_eigenvalue/60)',
            'physical_masses_identified': False}


def port_action(permutation):
    require(tuple(sorted(permutation)) == tuple(range(12)), 'a port permutation is required')
    result = []
    for state in STATES:
        out = [0]*12
        for i, value in enumerate(state):
            out[permutation[i]] = value
        result.append(INDEX[tuple(out)])
    return tuple(result)


@lru_cache(maxsize=1)
def seam_actions():
    result = []
    for a, b in quantum_mixing.seams():
        permutation = list(range(12))
        permutation[a], permutation[b] = b, a
        action = port_action(permutation)
        actual = np.zeros((66, 66), dtype=np.int64)
        expected = np.eye(66, dtype=np.int64)
        for i, state in enumerate(STATES):
            expected[i, action[i]] += 1
            for coin in (0, 1):
                final, _ = native.native(state, (a, b), coin)
                actual[i, INDEX[final]] += 1
        require(np.array_equal(actual, expected), 'native seam is not the asserted hold/exchange transfer')
        result.append(action)
    return tuple(result)


@lru_cache(maxsize=1)
def sector_projectors():
    """Integer numerators of three transfer commutant projections, denominator330."""
    incidence = np.array(STATES, dtype=np.int64)
    overlap = incidence @ incidence.T
    ones = np.ones((66, 66), dtype=np.int64)
    return (5*ones, 33*overlap-11*ones,
            330*np.eye(66, dtype=np.int64)-33*overlap+6*ones)


def transfer_sectors():
    projections = sector_projectors()
    require(np.array_equal(sum(projections), 330*np.eye(66, dtype=np.int64)),
            'transfer sectors fail completeness')
    for i, a in enumerate(projections):
        require(np.array_equal(a, a.T), 'a transfer sector is not self-adjoint')
        for j, b in enumerate(projections):
            require(np.array_equal(a @ b, 330*a if i == j else np.zeros_like(a)),
                    'transfer sectors fail orthogonality or idempotence')
        for action in seam_actions():
            require(np.array_equal(a[np.ix_(action, action)], a),
                    'a transfer sector does not commute with a native seam')
    # Independently solve commutation with seam permutations by exact orbits
    # of matrix entries, rather than infer it from the three displayed ranks.
    unseen = set(range(66*66))
    orbit_sizes = []
    while unseen:
        seed = min(unseen)
        orbit, queue = {seed}, deque([seed])
        while queue:
            left, right = divmod(queue.popleft(), 66)
            for action in seam_actions():
                pair = action[left]*66+action[right]
                if pair not in orbit:
                    orbit.add(pair); queue.append(pair)
        unseen.difference_update(orbit)
        orbit_sizes.append(len(orbit))
    require(sorted(orbit_sizes) == [66, 1320, 2970], 'common transfer commutant changed')
    ranks = [int(np.trace(a))//330 for a in projections]
    require(ranks == [1, 11, 54], 'transfer projector ranks changed')
    occupation = np.diag(np.array(STATES)[:, 0])
    require(np.any(occupation @ projections[2] != projections[2] @ occupation),
            'full occupation algebra unexpectedly preserves a transfer sector')
    return {'projector_denominator': 330, 'ranks': ranks,
            'projector_formulas': ['5J', '33BB^T-11J', '330I-33BB^T+6J'],
            'B_definition': 'B[x,a]=native occupation at port a',
            'common_seam_transfer_commutant_dimension': len(orbit_sizes),
            'matrix_entry_orbit_sizes': sorted(orbit_sizes),
            'all_30_native_seam_transfers_checked': True,
            'nonzero_subspace_invariant_under_every_seam_requires_all_12_ports': True,
            'full_occupation_effect_algebra_preserves_these_sectors': False,
            'particle_or_gauge_group_identification': False}


def modular_row_basis(matrix):
    matrix = np.asarray(matrix, dtype=np.int64).copy() % PRIME
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank+int(candidates[0])
        matrix[[pivot, rank]] = matrix[[rank, pivot]]
        matrix[rank] = matrix[rank]*pow(int(matrix[rank, column]), -1, PRIME) % PRIME
        others = np.flatnonzero(matrix[:, column])
        others = others[others != rank]
        matrix[others] = (matrix[others]-matrix[others, column, None]*matrix[rank]) % PRIME
        rank += 1
        if rank == min(matrix.shape):
            break
    return matrix[:rank]


def core_indices(core):
    checked_core(core)
    inside = [i for i, pair in enumerate(PAIRS) if set(pair) <= set(core)]
    return inside, [i for i in range(66) if i not in inside]


def compact_dimension_upper(core, matrix=None):
    """A rigorous complex dimension upper bound from reduction modulo a prime."""
    inside, outside = core_indices(core)
    if not inside:
        return 0, []
    matrix = counts() if matrix is None else np.asarray(matrix, dtype=np.int64)
    require(matrix.shape == (66, 66), 'the full66 transfer matrix is required')
    internal = matrix[np.ix_(inside, inside)]
    rows = modular_row_basis(matrix[np.ix_(outside, inside)])
    ranks = [len(rows)]
    while len(rows) < len(inside):
        enlarged = modular_row_basis(np.concatenate([rows, rows @ internal % PRIME]))
        if len(enlarged) == len(rows):
            break
        rows = enlarged
        ranks.append(len(rows))
    return len(inside)-len(rows), ranks


@lru_cache(maxsize=1)
def automorphisms():
    graph = nx.Graph()
    graph.add_nodes_from(range(12)); graph.add_edges_from(quantum_mixing.seams())
    result = tuple(sorted(tuple(mapping[i] for i in range(12)) for mapping in
                          nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter()))
    require(len(result) == 120 and len(set(result)) == 120, 'native automorphism enumeration changed')
    for permutation in result:
        action = port_action(permutation)
        require(np.array_equal(counts()[np.ix_(action, action)], counts()),
                'purported automorphism does not preserve the actual transfer')
    return result


def transformed_mask(mask, permutation):
    return sum(1 << permutation[i] for i in range(12) if mask >> i & 1)


def exact_compact_basis(core):
    """Exact rational lower certificate for a nonzero modular candidate."""
    inside, outside = core_indices(core)
    matrix = sp.Matrix(counts().tolist())
    internal = matrix.extract(inside, inside)
    rows = matrix.extract(outside, inside)
    _, independent = rows.T.rref()
    rows = rows[list(independent), :]
    while True:
        enlarged = rows.col_join(rows*internal)
        _, independent = enlarged.T.rref()
        enlarged = enlarged[list(independent), :]
        if enlarged.rows == rows.rows:
            rows = enlarged
            break
        rows = enlarged
    basis = rows.nullspace()
    lifted = sp.zeros(66, len(basis))
    for j, vector in enumerate(basis):
        for i, value in zip(inside, vector):
            lifted[i, j] = value
    require(lifted.rank() == len(basis), 'compact basis lost independence')
    require((matrix*lifted).row_join(lifted).rank() == len(basis),
            'compact basis is not invariant')
    return lifted


@lru_cache(maxsize=1)
def compact_search():
    require(sp.isprime(PRIME), 'modular rank modulus is not prime')
    representatives = {}
    rows = []
    for mask in range(4096):
        if mask in representatives:
            continue
        orbit = {transformed_mask(mask, permutation) for permutation in automorphisms()}
        require(not (orbit & representatives.keys()), 'core orbits overlap')
        for item in orbit:
            representatives[item] = mask
        core = tuple(i for i in range(12) if mask >> i & 1)
        upper, ranks = compact_dimension_upper(core)
        if upper and len(core) < 12:
            require(exact_compact_basis(core).cols == upper,
                    'rational witness does not attain modular dimension bound')
        rows.append({'mask': mask, 'core': list(core), 'orbit_size': len(orbit),
                     'compact_invariant_dimension': upper, 'modular_rank_sequence': ranks})
    require(len(representatives) == 4096 and len(rows) == 82, 'compact search coverage changed')
    dimensions = {row['mask']: row['compact_invariant_dimension'] for row in rows}
    summary = Counter((mask.bit_count(), dimensions[representatives[mask]]) for mask in range(4096))
    expected = {(k, 0): int(sp.binomial(12, k)) for k in range(10)}
    expected.update({(10, 0): 60, (10, 3): 6, (11, 3): 12, (12, 66): 1})
    require(dict(summary) == expected, 'complete compact-support classification changed')
    return {'cores_checked': 4096, 'automorphism_representatives': len(rows),
            'modular_prime': PRIME, 'all_positive_dimensions_attained_by_exact_rational_basis': True,
            'classification': [{'ports': k, 'dimension': d, 'core_count': n}
                               for (k, d), n in sorted(summary.items())],
            'representatives': rows,
            'minimum_nonzero_compact_core_ports': 10,
            'quantified_class': 'all complex invariant subspaces of P supported on configurations with both loads in a fixed port core'}


@lru_cache(maxsize=1)
def compact_modes():
    core = (0, 1, 2, 3, 4, 5, 6, 7, 9, 10)
    basis = exact_compact_basis(core)
    matrix = sp.Matrix(counts().tolist())
    restricted = (basis.T*basis).inv()*basis.T*matrix*basis
    modes = []
    for value, multiplicity, vectors in restricted.eigenvects():
        require(multiplicity == 1 and len(vectors) == 1, 'compact spectrum is not simple')
        vector = basis*vectors[0]
        vector *= sp.ilcm(*[a.q for a in vector])
        vector /= sp.igcd(*[int(a) for a in vector if a])
        vector = np.array([int(a) for a in vector], dtype=np.int64)
        require(np.array_equal(counts() @ vector, int(value)*vector), 'compact mode equation failed')
        require(int(vector.sum()) == 0 and int(vector @ vector) == 20, 'compact fluctuation normalization changed')
        require(np.array_equal(sector_projectors()[2] @ vector, 330*vector),
                'compact mode leaves the correlation transfer sector')
        modes.append((int(value), vector))
    require([value for value, _ in modes] == [48, 50, 52], 'compact eigenvalues changed')
    return core, tuple(modes)


def recorded_maps():
    """Actual endpoint execution, full raw menus, independent coarse-map check."""
    rows = []
    for law in ('C', 'E'):
        for seam in quantum_mixing.seams():
            coarse = {(coin, label): mapping for coin, label, mapping in coins.coin_maps(law, seam)}
            for coin in (0, 1):
                raw = {i: [] for i in range(1 << circuit.RECORD_BITS[law])}
                for i, state in enumerate(STATES):
                    final, (_, current) = native.native(state, seam, coin)
                    label = circuit.CURRENT_CODES[current]
                    if law == 'E':
                        label += 4*((state[seam[0]]+state[seam[1]]) % 2)
                    raw[label].append((i, INDEX[final]))
                for label, mapping in raw.items():
                    decoded = circuit.record_values(law, coin).get(label)
                    require(tuple(mapping) == coarse.get((coin, decoded), ()),
                            'complete raw branch differs from native coarse instrument')
                    rows.append((law, seam, coin, label, tuple(mapping)))
    require(len(rows) == 720, 'raw source menu lost outcomes')
    return rows


def mode_and_instrument_certificate():
    core, modes = compact_modes()
    all_maps = recorded_maps()
    sector = sector_projectors()[2]
    result = []
    for value, vector in modes:
        norm = int(vector @ vector)
        # Search one all-positive schedule perturbation retaining every seam.
        outside = core_indices(core)[1]
        witness = None
        for seam, action in zip(quantum_mixing.seams(), seam_actions()):
            seam_image = vector+vector[list(action)]  # Twice P_seam vector.
            leakage = int(seam_image[outside] @ seam_image[outside])
            if leakage:
                witness = {'seam': list(seam),
                           'schedule': 'one-half uniform plus one-half this seam; every seam has positive weight',
                           'one_attempt_exterior_amplitude_norm_squared': str(Q(leakage, 16*norm))}
                break
        require(witness is not None, 'compact mode survived every individual seam')
        operational = {}
        for law in ('C', 'E'):
            total_norm = sector_numerator = fidelity_numerator = 0
            for name, _, _, _, mapping in all_maps:
                if name != law:
                    continue
                output = np.zeros(66, dtype=np.int64)
                for source, target in mapping:
                    output[target] = vector[source]
                total_norm += int(output @ output)
                sector_numerator += int(output @ sector @ output)
                fidelity_numerator += int(vector @ output)**2
            require(total_norm == 60*norm, 'complete native mode output lost probability')
            weight = Q(sector_numerator, 60*norm*330)
            require(weight < 1, 'recorded native law unexpectedly preserves the transfer sector')
            operational[law] = {'complete_probability': '1',
                                'correlation_sector_weight_after_one_attempt': str(weight),
                                'initial_mode_fidelity_after_one_attempt': str(Q(fidelity_numerator, 60*norm*norm))}
        result.append({'integer_eigenvalue': value, 'transfer_eigenvalue': str(Q(value, 60)),
                       'reconstructed_energy': f'log(60/{value})',
                       'integer_mode': [{'pair': list(PAIRS[i]), 'value': int(a)}
                                        for i, a in enumerate(vector) if a],
                       'squared_norm': norm,
                       'stationary_native_variance': str(Q(norm, 66)),
                       'normalized_native_autocorrelation_at_n': f'({Q(value, 60)})^n',
                       'integrated_normalized_autocorrelation': str(1/(1-Q(value, 60))),
                       'all_positive_nonuniform_schedule_control': witness,
                       'full_recorded_instrument': operational})
    inside = core_indices(core)[0]
    exterior = []
    for length in (1, 2, 4):
        power = reconstruction.path_counts(length)
        trace = sum(power[i][i] for i in range(66))
        distance = 1-Q(sum(power[i][i] for i in inside), trace)
        exterior.append({'closed_length': length, 'full_exterior_trace_distance': str(distance)})
    return {'core': list(core), 'omitted_ports': [8, 11], 'modes': result,
            'full_raw_C_E_slots': len(all_maps),
            'zero_raw_slots_retained': sum(not row[-1] for row in all_maps),
            'full_exterior_reference_comparison': exterior,
            'ground_reference_exterior_trace_distance': '7/22',
            'coherent_mode_preparation_or_native_measurement_implemented': False,
            'native_classical_readout': 'the integer mode is an explicit signed function of the two occupied ports; its equilibrium autocorrelation follows from native counts',
            'small_localized_reference_like_protected_particle': False}


def local_phase_excitation():
    """A faithful, exterior-matched positive-energy state of the count algebra.

    Conjugation here constructs a reconstructed functional. It does not execute
    a native phase operation or prepare the state on a quantum working carrier.
    """
    matrix = counts()
    trace = int(np.trace(matrix))
    states = np.array(STATES, dtype=np.int64)
    for port in range(12):
        signs = 1-2*states[:, port]
        difference = matrix*signs[:, None]*signs[None, :]-matrix
        require(np.all(difference[states[:, port, None] == states[None, :, port]] == 0),
                'local phase changes the exterior restriction')
        require(int(np.sum(difference*difference)) == 400 and int(np.trace(difference)) == 0,
                'one-port phase commutator changed')
    port = 0
    signs = 1-2*states[:, port]
    difference = matrix*signs[:, None]*signs[None, :]-matrix
    commutator_norm = Q(int(np.sum(difference*difference)), 60**2)
    trace_p = Q(trace, 60)
    lower = commutator_norm/(2*trace_p)
    upper = 3*commutator_norm/(4*trace_p)
    require((trace_p, commutator_norm, lower, upper) ==
            (56, Q(1, 9), Q(1, 1008), Q(1, 672)), 'energy bound changed')
    # Singular values of the rectangular cross-occupation block determine
    # the exact initial half trace norm; no floating diagonalization is used.
    occupied = np.flatnonzero(states[:, port])
    empty = np.flatnonzero(1-states[:, port])
    cross = matrix[np.ix_(occupied, empty)]
    x = sp.Symbol('x')
    gram = sp.Matrix((cross @ cross.T).tolist()).charpoly(x).as_expr()
    require(sp.expand(gram-(x-8)*(x-3)**4*(x-5)**6) == 0,
            'local phase trace-distance spectrum changed')
    entries = tuple(zip(*np.nonzero(difference)))
    rows = recorded_maps()
    nonzero_branches = 0
    for _, _, _, _, mapping in rows:
        if mapping:
            displacements = {tuple(states[target]-states[source]) for source, target in mapping}
            require(len(displacements) == 1, 'a native raw branch changes signed input differences')
            nonzero_branches += 1
    laws = {}
    for law, multiplier in (('C', 50), ('E', 42)):
        accumulated = np.zeros((66, 66), dtype=np.int64)
        survival_counts = Counter()
        for name, _, _, _, mapping in rows:
            if name != law or not mapping:
                continue
            inputs, outputs = np.array(mapping, dtype=int).T
            accumulated[np.ix_(outputs, outputs)] += difference[np.ix_(inputs, inputs)]
            mapping_dict = dict(mapping)
            for i, j in entries:
                if i in mapping_dict and j in mapping_dict:
                    require(difference[mapping_dict[i], mapping_dict[j]] == difference[i, j],
                            'phase-difference support is not closed under native history')
                    survival_counts[i, j] += 1
        require(np.array_equal(accumulated, multiplier*difference),
                'phase difference is not the asserted native eigenoperator')
        require(set(survival_counts.values()) == {multiplier} and len(survival_counts) == 100,
                'complete-record entrywise contraction changed')
        laws[law] = {'unread_difference_multiplier': str(Q(multiplier, 60)),
                     'all_history_entrywise_norm_multiplier': str(Q(multiplier, 60)),
                     'all_history_half_trace_distance_upper_prefactor': '5/168',
                     'reference_for_comparison': 'the same native evolution of rho, not a fixed stationary rho'}
    return {'core': [port], 'reference': 'rho=P/Tr(P)',
            'excited_functional': 'rho_prime=D rho D; D[x,x]=(-1)^x_0',
            'faithful_and_same_entropy': True,
            'full_11_port_exterior_exactly_equal': True,
            'all_12_one_port_phase_locality_and_energy_norms_checked': True,
            'configuration_probabilities_exactly_equal': True,
            'hamiltonian': 'H=-log(P)', 'trace_P': str(trace_p),
            'commutator_P_D_Frobenius_squared': str(commutator_norm),
            'strict_positive_energy_excess_interval': [str(lower), str(upper)],
            'exact_initial_half_trace_distance': '(sqrt(8)+4sqrt(3)+6sqrt(5))/1680',
            'native_nonzero_raw_displacement_branches_checked': nonzero_branches,
            'full_native_record_plus_exterior_equal_at_every_finite_history': True,
            'native_decay_comparison': laws,
            'state_prepared_by_native_phase_control': False,
            'reconstructed_phase_functional_has_native_history_count_dictionary': True,
            'physical_localized_particle_or_mass': False}


def source_bindings():
    parent = strict_json.load(ROOT/'model/native_quantum_reconstruction/results.json')
    for name, expected in parent['source_sha256'].items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected,
                'native reconstruction parent source is stale: '+name)
    names = set(parent['source_sha256']) | {
        'model/native_quantum_reconstruction/results.json',
        'model/repair_scheduler/coins.py', 'model/observer/source_discrimination/circuit.py',
        'model/common_source/native_defect_decision.py',
        'model/common_source/native_defect_decision_results.json',
        'model/native_quantum_defects/certificate.py', 'model/native_quantum_defects/PROOF.md',
        'model/native_quantum_defects/test_certificate.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def report():
    return {'schema': 'oph.native_quantum_defects.v1', 'all_pass': True,
            'source': {'states': 66, 'seams': 30, 'placement_coins': 2,
                       'complete_endpoint_branches': 3960,
                       'law': 'native fair endpoint repair; independent uniform seam attempts',
                       'state_and_field_dictionary': 'binary total-two configurations and explicit occupation functions'},
            'spectrum': transfer_spectrum(), 'transfer_sectors': transfer_sectors(),
            'compact_support_search': compact_search(), 'constructed_modes': mode_and_instrument_certificate(),
            'localized_positive_energy_functional': local_phase_excitation(),
            'scope': {'all_complex_compact_subspaces_on_this_fixed_carrier_exhausted': True,
                      'all_refinements_or_history_observers_exhausted': False,
                      'supplied_coherent_service_used': False, 'native_unitary_time_identified': False,
                      'physical_particle_or_mass_identified': False,
                      'existing_collision_and_classical_shell_exclusions_retained': True,
                      'noncompact_local_density_perturbations_excluded': False,
                      'surviving_routes': ['controlled noncompact or approximate localization',
                                           'different admitted scheduler or reference with its own proof',
                                           'relational history observables and explicitly accounted active recovery',
                                           'refining carriers with proved compatible source maps']},
            'source_sha256': source_bindings()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(path) == value, 'native quantum defect receipt differs from exact replay')
    print(json.dumps({'all_pass': True, 'modes_retained': 66, 'cores_checked': 4096,
                      'smallest_compact_core': 10, 'physical_particle_identified': False}))


if __name__ == '__main__':
    main()
