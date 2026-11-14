"""Geometry and joint fields from two visible native incident-event clocks.

The comparison interval ends at a tick of A. Only the ordered union of A/B
incident records is read; no remote attempt counter or dimensional clock is
needed. All remote waiting words are probability-summed, including their tail.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

import sympy as sp

from ..observer import native
from ..observer_reconstruction import native_geometry as geometry
from ..primitive_observer import transfer
from ..repair_identification import core as exact_core

# The retained carrier's package initializers import these modules even though
# the calculation calls only its seam provider. Pin the executable import
# closure without treating those modules' scientific claims as premises.
_ARCHIVE = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
_ARCHIVE_IMPORTS = tuple(_ARCHIVE+name for name in (
    'oph_exact/__init__.py', 'oph_fpe/__init__.py', 'oph_fpe/claims.py',
    'oph_fpe/constants/__init__.py', 'oph_fpe/constants/oph_pixel.py',
    'oph_fpe/core/__init__.py', 'oph_fpe/core/array_geometry.py',
    'oph_fpe/core/echosahedral_dynamics.py', 'oph_fpe/core/echosahedral_federation.py',
    'oph_fpe/core/graph.py', 'oph_fpe/core/icosahedral.py', 'oph_fpe/core/patch_state.py',
    'oph_fpe/core/patchnet.py', 'oph_fpe/core/pixel_scale.py',
    'oph_fpe/core/screen_microphysics.py', 'oph_fpe/core/screen_ports.py',
    'oph_fpe/dynamics/__init__.py', 'oph_fpe/dynamics/anneal.py',
    'oph_fpe/dynamics/consensus_certificate.py', 'oph_fpe/dynamics/kernel_dispatcher.py',
    'oph_fpe/dynamics/modular.py', 'oph_fpe/dynamics/positive_geometry.py',
    'oph_fpe/dynamics/repair.py', 'oph_fpe/evidence/__init__.py',
    'oph_fpe/evidence/bundle.py', 'oph_fpe/evidence/cross_repo_artifacts.py',
    'oph_fpe/evidence/hashes.py', 'oph_fpe/evidence/particle_input_policy.py',
    'oph_fpe/evidence/production_envelope.py', 'oph_fpe/evidence/verifier.py',
    'oph_fpe/finite_groups.py', 'oph_fpe/gauge/__init__.py',
    'oph_fpe/gauge/a5_sm_certificate.py', 'oph_fpe/gauge/covariant_overlap.py',
    'oph_fpe/gauge/higgs_carrier.py', 'oph_fpe/gauge/mar_sieve.py',
    'oph_fpe/gauge/physical_a5_sm_requirements.py', 'oph_fpe/gauge/repair_projection.py',
    'oph_fpe/gauge/source_to_eft.py', 'oph_fpe/gauge/yang_mills_gap.py',
    'oph_fpe/groups/__init__.py', 'oph_fpe/groups/base.py', 'oph_fpe/groups/clock.py',
    'oph_fpe/groups/registry.py', 'oph_fpe/groups/s3.py', 'oph_fpe/groups/z2.py'))
SOURCE_PATHS = tuple(sorted(set((*transfer.SOURCE_PATHS, *exact_core.SOURCE_PATHS,
    *_ARCHIVE_IMPORTS,
    'model/laboratory_clock/event_geometry.py',
    'model/laboratory_clock/__init__.py', 'model/repair_identification/__init__.py',
    'model/observer_reconstruction/__init__.py',
    'model/observer_reconstruction/native_geometry.py',
    'model/__init__.py', 'model/anharmonic_transfer/__init__.py',
    'model/anharmonic_transfer/environment.py', 'model/common_source/__init__.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/native_join.py', 'model/common_source/protocol.py',
    'model/common_source/records.py', 'model/common_source/scalar_coupling.py',
    'model/common_source/smooth_protocol.py', 'model/observer/__init__.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py', 'model/observer/source_scalar/__init__.py',
    'model/observer/source_scalar/bounds.py', 'model/observer/source_scalar/compiler.py',
    'model/observer/source_scalar/cutoff.py', 'model/observer/source_scalar/experiment.py',
    'model/observer/source_scalar/field.py', 'model/predictions/__init__.py',
    'model/primitive_observer/__init__.py', 'model/record_laws/__init__.py',
    'model/relativistic_scalar/__init__.py', 'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/__init__.py', 'model/scalar_observations/lattice.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py'))))
PORTS = 12
RHO = sp.Rational(1, 6)
CHI = sp.Rational(5, 33)
P0 = sp.ImmutableMatrix(sp.eye(PORTS)-sp.ones(PORTS)/PORTS)


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original bounded integer is required')
    return value


def _ports(a, b):
    _integer(a, 0, 11); _integer(b, 0, 11)
    if a == b:
        raise ValueError('two distinct clock ports are required')


def _argument(value):
    if type(value) is not F:
        raise ValueError('an exact canonical Fraction is required')
    value = exact_core._fraction(value)
    if not 0 <= value <= 1:
        raise ValueError('a probability-generating argument in [0,1] is required')
    return value


@lru_cache(maxsize=1)
def _source():
    return transfer.native_counts()


def renewal_parameters(delta):
    """Independent Bernoulli plus geometric-on-zero decomposition of one count.

    delta is the inferred adjacency bit. The geometric success probability
    and terminating-event Bernoulli parameter follow from the same four
    incidence classes; they are not additional stochastic clock premises.
    """
    _integer(delta, 0, 1)
    return {'terminal_B_probability': F(delta, 5),
            'geometric_success_probability': F(5, 10-delta),
            'mean': F(1), 'variance': F(2)-F(2*delta, 5)}


def clock_pgf(a, b, z, count=1):
    """E[z**B(T_A(count))], including a simultaneous terminating B tick."""
    _ports(a, b); _integer(count, 1, 64)
    z = _argument(z)
    probabilities = {(u, v): F(0) for u in (0, 1) for v in (0, 1)}
    for seam in _source().seams:
        probabilities[int(a in seam), int(b in seam)] += F(1, 30)
    p00, p01 = probabilities[0, 0], probabilities[0, 1]
    p10, p11 = probabilities[1, 0], probabilities[1, 1]
    return ((p10+p11*z)/(1-p00-p01*z))**count


def count_law(a, b, count=1, maximum=12):
    """Exact finite count bins plus the complete overflow probability."""
    _ports(a, b); _integer(count, 1, 64); _integer(maximum, 0, 128)
    adjacent = int(tuple(sorted((a, b))) in _source().seams)
    ratio = F(5-adjacent, 10-adjacent)
    single = [ratio]+[(ratio*ratio+F(adjacent, 10-adjacent))*ratio**(k-1)
                       for k in range(1, maximum+1)]
    law = [F(1)]+[F(0)]*maximum
    for _ in range(count):
        law = [sum((law[j]*single[k-j] for j in range(k+1)), F(0))
               for k in range(maximum+1)]
    return tuple(law), 1-sum(law, F(0))


def clock_records():
    """Pair labels and zero-B-count probabilities; no seam labels or fields."""
    return [{'ports': [a, b], 'zero_B_before_next_A': str(clock_pgf(a, b, F(0)))}
            for a, b in combinations(range(PORTS), 2)]


def recover(records):
    """Recover adjacency and its unit rank-three Gram from clock laws only.

    The estimator receives exact ensemble probabilities, not sampled estimates.
    The uniform exogenous scheduler and five incident seams per port are the
    declared class being tested. Malformed or incompatible tables are rejected.
    """
    if type(records) not in (list, tuple) or len(records) != 66:
        raise ValueError('all 66 distinct pair-clock laws are required')
    laplacian, seen = sp.zeros(PORTS), set()
    for row in records:
        if type(row) is not dict or set(row) != {'ports', 'zero_B_before_next_A'}:
            raise ValueError('only pair labels and their declared clock probability are admitted')
        ports = row['ports']
        if type(ports) not in (list, tuple) or len(ports) != 2:
            raise ValueError('two clock port labels are required')
        a, b = ports; _ports(a, b)
        pair = tuple(sorted((a, b)))
        if pair in seen:
            raise ValueError('duplicate clock pair')
        seen.add(pair)
        token = row['zero_B_before_next_A']
        if type(token) is not str:
            raise ValueError('a canonical exact probability string is required')
        try:
            probability = F(token)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError('a canonical exact probability string is required') from error
        if str(probability) != token or probability not in (F(4, 9), F(1, 2)):
            raise ValueError('clock probability is outside the declared source class')
        if probability == F(4, 9):
            laplacian[a, a] += 1; laplacian[b, b] += 1
            laplacian[a, b] -= 1; laplacian[b, a] -= 1
    if any(laplacian[i, i] != 5 for i in range(PORTS)):
        raise ValueError('clock laws do not recover the declared degree-five carrier')
    checks, projector = geometry._band_certificate(
        [[F(value) for value in row] for row in laplacian.tolist()])
    if not all(checks.values()):
        raise ValueError('clock-inferred graph fails the exact spatial band certificate')
    gram = sp.ImmutableMatrix([[4*(sp.Rational(a)+sp.Rational(b)*sp.sqrt(5))
                                for a, b in row] for row in projector])
    return sp.ImmutableMatrix(laplacian), gram


def geometry_error_bound(intervals):
    """Union Chebyshev bound for 66 zero-count classifiers at threshold 17/36.

    Each pair uses independent successive A-to-A comparison intervals. Pair
    datasets may overlap. No independence between different pairs is assumed.
    """
    _integer(intervals, 1, 10**12)
    return min(F(1), F(66*324, intervals))


def recover_from_counts(rows):
    """Apply the fixed zero-count threshold to finite pair-comparison samples.

    The error budget bounds rejection or an incorrect graph jointly. It is
    unconditional under the declared source class, not a confidence level
    conditioned on passing the structural checks or on selected samples.
    """
    if type(rows) not in (list, tuple) or len(rows) != 66:
        raise ValueError('all 66 pair-comparison samples are required')
    records, error = [], F(0)
    for row in rows:
        if type(row) is not dict or set(row) != {'ports', 'intervals', 'zero_counts'}:
            raise ValueError('pair labels, interval count and zero-count total are required')
        n = _integer(row['intervals'], 1, 10**12)
        zero = _integer(row['zero_counts'], 0, n)
        probability = F(4, 9) if F(zero, n) < F(17, 36) else F(1, 2)
        records.append({'ports': row['ports'], 'zero_B_before_next_A': str(probability)})
        error += F(324, n)
    laplacian, gram = recover(records)
    return {'laplacian': laplacian, 'unit_gram': gram,
            'unconditional_rejection_or_wrong_graph_probability_bound': min(F(1), error)}


def predicted_field(a, b, z=F(1), count=2, *, records=None):
    """Joint PGF of initial A bit and final B bit from clock-inferred geometry.

    The stationary binary total-two preparation is declared, not inferred from
    the clocks. Only two initial local bits and their own currents are needed
    to decode these observables. Hidden initial configurations are integrated.
    """
    _ports(a, b); _integer(count, 1, 8)
    z = _argument(z)
    laplacian, _ = recover(clock_records() if records is None else records)
    remote, incident = sp.zeros(PORTS), sp.zeros(PORTS)
    for x, y in combinations(range(PORTS), 2):
        if not laplacian[x, y]:
            continue
        vector = sp.zeros(PORTS, 1); vector[x] = 1; vector[y] = -1
        step = sp.Rational(z)**int(b in (x, y))*(sp.eye(PORTS)-vector*vector.T/2)/30
        if a in (x, y):
            incident += step
        else:
            remote += step
    kernel = incident*(sp.eye(PORTS)-remote).inv()
    mass = F((sp.ones(1, PORTS)*kernel**count)[0])
    covariance = F(CHI*(kernel**count*P0)[b, a])
    p11 = mass/36+covariance
    values = {(1, 1): p11, (1, 0): mass/6-p11,
              (0, 1): mass/6-p11, (0, 0): 2*mass/3+p11}
    if any(value < 0 for value in values.values()) or sum(values.values()) != mass:
        raise ArithmeticError('joint clock/field generating weights lost positivity or mass')
    return values


@lru_cache(maxsize=4)
def _native_groups(a, b):
    """Independent complete 66-state endpoint enumeration and local decoding."""
    _ports(a, b)
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    groups = {(u, v): sp.zeros(66) for u in (0, 1) for v in (0, 1)}
    checks = 0
    for column, before in enumerate(states):
        for seam in _source().seams:
            for coin in (0, 1):
                after, (_, current) = native.native(before, seam, coin)
                for port in (a, b):
                    decoded = before[port]+current*(int(port == seam[1])-int(port == seam[0]))
                    if decoded != after[port]:
                        raise ArithmeticError('local initial bit and own currents lost readback')
                groups[int(a in seam), int(b in seam)][index[after], column] += sp.Rational(1, 60)
                checks += 1
    return states, {key: sp.ImmutableMatrix(value) for key, value in groups.items()}, checks


def native_field(a, b, z=F(1), count=2):
    """Full-state oracle with every unobserved remote waiting word summed."""
    _ports(a, b); _integer(count, 1, 8)
    z = _argument(z)
    states, groups, _ = _native_groups(a, b)
    remote = groups[0, 0]+sp.Rational(z)*groups[0, 1]
    incident = groups[1, 0]+sp.Rational(z)*groups[1, 1]
    # Two initial-bit subensembles; no complete-state data enter the observer.
    state = sp.Matrix(66, 2, lambda i, bit: sp.Rational(int(states[i][a] == bit), 66))
    inverse = (sp.eye(66)-remote).inv(method='DM')
    for _ in range(count):
        state = incident*inverse*state
    return {(first, last): F(sum(state[i, first] for i in range(66) if states[i][b] == last))
            for first in (0, 1) for last in (0, 1)}


def identity_field_control(a, b, z=F(1), count=2):
    """Change repair to identity, retaining the full scheduler/clock process."""
    _ports(a, b); _integer(count, 1, 8)
    z = _argument(z)
    states = tuple(native.states(binary=True))
    mass = clock_pgf(a, b, z, count)
    return {(first, last): mass*F(sum(state[a] == first and state[b] == last
                                    for state in states), len(states))
            for first in (0, 1) for last in (0, 1)}


def run():
    records = clock_records()
    laplacian, gram = recover(records)
    expected_edges = {tuple(sorted(edge)) for edge in _source().seams}
    recovered_edges = {(a, b) for a, b in combinations(range(PORTS), 2) if laplacian[a, b] == -1}
    if recovered_edges != expected_edges:
        raise ArithmeticError('clock-only geometry differs from the executed source')
    distances = {0: 0}; queue = [0]
    for a in queue:
        for b in range(PORTS):
            if laplacian[a, b] == -1 and b not in distances:
                distances[b] = distances[a]+1; queue.append(b)
    rows, oracle_checks, native_11, identity_11 = [], 0, {}, {}
    for distance in (1, 2, 3):
        b = min(port for port in range(1, PORTS) if distances[port] == distance)
        for count, z in ((1, F(0)), (1, F(1, 2)), (2, F(1))):
            predicted = predicted_field(0, b, z, count, records=records)
            if predicted != native_field(0, b, z, count):
                raise ArithmeticError('clock-derived field prediction fails complete native replay')
            oracle_checks += 1
            rows.append({'ports': [0, b], 'hop_distance': distance, 'A_ticks': count,
                'B_count_PGF_argument': str(z),
                'joint_bit_generating_weights': [{'initial_A': a, 'final_B': c, 'weight': str(p)}
                    for (a, c), p in sorted(predicted.items())]})
            if count == 2 and z == 1:
                native_11[str(distance)] = str(predicted[1, 1])
                identity_11[str(distance)] = str(identity_field_control(0, b)[1, 1])
    return {'clock': 'visible incident events, including zero currents',
        'comparison': 'B count through and including the nth A event in their ordered union tape',
        'one_interval_PGF': '(5-delta+delta*z)/(10-delta-(5-delta)*z)',
        'delta': 'one for adjacency, zero otherwise',
        'mean_B_per_A': '1', 'variance_B_per_A': {'adjacent': '8/5', 'nonadjacent': '2'},
        'clock_only_records': records, 'laplacian': [[str(x) for x in row] for row in laplacian.tolist()],
        'unit_gram': [[str(x) for x in row] for row in gram.tolist()],
        'geometry_identified_up_to': 'orthogonal isometry and transported port labels',
        'finite_sample_all_pairs_error_bound': 'min(1,21384/N), threshold 17/36 for zero-count frequencies',
        'joint_native_field_predictions': rows, 'independent_full_state_comparisons': oracle_checks,
        'native_endpoint_branches_per_pair': _native_groups(0, 1)[2],
        'same_clock_identity_repair_control': {'complete_clock_laws_identical': True,
            'two_tick_joint_11_by_shell': identity_11,
            'native_joint_11_by_shell': native_11},
        'scope': {'hidden_attempt_or_movement_counter_read': False,
            'remote_current_tape_read': False, 'joint_incident_order_and_coincidence_read': True,
            'only_initial_local_bits_and_own_currents_decode_fields': True,
            'uniform_exogenous_scheduler': True, 'stationary_total_two_field_preparation': True,
            'all_remote_waiting_tails_summed': True, 'unit_choice_requires_external_calibration': False,
            'selected_classical_readout_marginal_only': True,
            'automatic_identification_with_observed_universe': False,
            'classical_clock_and_field_records_select_quantum_C_over_E': False}}
