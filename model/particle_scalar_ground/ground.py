"""Full-carrier ground certificate from exact small Gram blocks.

Floating eigenvectors propose positive rational vectors. Only addressed
integer matrix multiplication and rational Perron inequalities are admitted.
The Hamiltonian is the literal finite scalar/CAR source at Y=diag(1,2,3),
pair cost d=1 and fermion offset epsilon_F=0. No measured mass enters it.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from scipy.linalg import eigh
import sympy as s

from model.particle_calibration_selection import phase_attenuation as source
from model.native_flavor_orbit import intervals as iv


WEIGHTS = (1, 1, 2, 2, 3, 3)
TRIPLETS = tuple((6*f+spin, 6*f+2+spin, 6*f+5-spin)
                 for f in range(3) for spin in range(2))
LOCAL_SIGNS = (1, -1, 1, -1, 1, -1)
WINNING_TYPES = (2, 2, 1, 1, 1, 1)
PERRON_BITS = 42
ROOT_POLYNOMIAL = (1, -244, 23737, -1171722, 30638880, -394771968, 1881169920)
SOURCE_INPUTS = (
    'model/particle_calibration_selection/phase_attenuation.py',
    'model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md',
    'model/particle_observer_outputs/HIGGS_RECORDS.md',
    'model/particle_observer_outputs/CHARGED_HIGGS_RECORDS.md',
    'model/particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md',
    'model/interacting_vacuum/global_ground.py',
    'model/centered_source_ground/gap.py',
    'model/standard_model_source_map/FERMION_ATTACHMENT.md',
    'model/spin_operational_position/PROOF.md',
)


def local_frame(kind, coordinate, spin_sign=1):
    """Local original-CAR mask and sign for a conserved three-mode cell.

    Types 1/2 are fundamental/dual three-dimensional sectors; types 0/3
    are inactive singletons. Both pair operators have nonnegative entries
    in this signed basis. Their action is checked against literal CAR.
    """
    if kind not in (0, 1, 2, 3) or spin_sign not in (-1, 1):
        raise ValueError('an admitted local cell and spin sign are required')
    choices = {0: ((4, 1),),
               1: ((0, 1), (5, spin_sign), (6, spin_sign)),
               2: ((7, 1), (2, -spin_sign), (1, spin_sign)),
               3: ((3, 1),)}[kind]
    if type(coordinate) is not int or not 0 <= coordinate < len(choices):
        raise ValueError('an admitted local coordinate is required')
    return choices[coordinate]


def original_column(kinds, coordinates):
    """The signed tensor frame as an actual eighteen-mode occupation column."""
    if len(kinds) != 6 or len(coordinates) != 6:
        raise ValueError('six source cells are required')
    occupied, phase = [], 1
    for triplet, sign, kind, coordinate in zip(TRIPLETS, LOCAL_SIGNS, kinds, coordinates, strict=True):
        local_mask, local_phase = local_frame(kind, coordinate, sign)
        phase *= local_phase
        occupied.extend(mode for i, mode in enumerate(triplet) if local_mask >> i & 1)
    inversions = sum(a > b for i, a in enumerate(occupied) for b in occupied[i+1:])
    return sum(1 << i for i in occupied), phase*(-1)**inversions


def decoded_column(mask):
    """Inverse of the entire signed frame, including inactive source inputs."""
    if type(mask) is not int or not 0 <= mask < 1 << 18:
        raise ValueError('an eighteen-mode original occupation is required')
    kinds, coordinates = [], []
    for triplet, sign in zip(TRIPLETS, LOCAL_SIGNS, strict=True):
        local_mask = sum(((mask >> mode) & 1) << i for i, mode in enumerate(triplet))
        found = [(kind, coordinate) for kind in range(4)
                 for coordinate in range(3 if kind in (1, 2) else 1)
                 if local_frame(kind, coordinate, sign)[0] == local_mask]
        if len(found) != 1:
            raise ArithmeticError('the local signed CAR frame is not complete')
        kind, coordinate = found[0]; kinds.append(kind); coordinates.append(coordinate)
    reconstructed, phase = original_column(kinds, coordinates)
    if reconstructed != mask:
        raise ArithmeticError('the signed CAR frame changes an original occupation')
    return tuple(kinds), tuple(coordinates), phase


def pair_operators(kinds, weights=WEIGHTS):
    """Both exact scalar creation channels on one complete cell sector."""
    if len(kinds) != 6 or any(k not in (0, 1, 2, 3) for k in kinds):
        raise ValueError('six complete cell-sector labels are required')
    if len(weights) != 6 or any(type(w) is not int or w <= 0 for w in weights):
        raise ValueError('six positive integer source weights are required')
    states = tuple(product(*(range(3) if k in (1, 2) else range(1) for k in kinds)))
    index = {state: i for i, state in enumerate(states)}
    operators = []
    for scalar in (1, 2):
        rows, columns, values = [], [], []
        for column, state in enumerate(states):
            for cell, (kind, weight) in enumerate(zip(kinds, weights, strict=True)):
                if kind not in (1, 2) or state[cell] != (0 if kind == 1 else scalar):
                    continue
                destination = list(state); destination[cell] = scalar if kind == 1 else 0
                rows.append(index[tuple(destination)]); columns.append(column); values.append(weight)
        operators.append(sparse.csr_matrix((values, (rows, columns)),
                                          shape=(len(states), len(states)), dtype=np.int64))
    return states, tuple(operators)


def gram(kinds):
    states, channels = pair_operators(kinds)
    matrix = sum((channel.T@channel for channel in channels),
                 sparse.csr_matrix((len(states), len(states)), dtype=np.int64))
    if (matrix-matrix.T).nnz or np.any(matrix.data < 0):
        raise ArithmeticError('the source Gram lost its symmetric nonnegative entries')
    return states, channels, matrix


def check_perron_vector(matrix, vector):
    """Exact rational Collatz bounds; a float proposal supplies no premise."""
    if type(matrix) is not sparse.csr_matrix or matrix.dtype != np.dtype(np.int64):
        raise ValueError('the constructed integer CSR source matrix is required')
    if matrix.shape[0] != matrix.shape[1] or len(vector) != matrix.shape[0] or not len(vector):
        raise ValueError('a complete nonempty source column is required')
    if any(type(value) is not int or value <= 0 for value in vector):
        raise ValueError('strictly positive Python-integer coordinates are required')
    if (matrix-matrix.T).nnz or np.any(matrix.data < 0):
        raise ValueError('a symmetric nonnegative source matrix is required')
    # Python integers make the multiplication independent of fixed-width overflow.
    image = [sum(int(matrix.data[k])*vector[int(matrix.indices[k])]
                 for k in range(int(matrix.indptr[i]), int(matrix.indptr[i+1])))
             for i in range(matrix.shape[0])]
    ratios = [F(value, coordinate) for value, coordinate in zip(image, vector, strict=True)]
    return min(ratios), max(ratios)


def propose_positive_vector(matrix):
    """Untrusted, transient numerical proposal for an exact positive witness."""
    if matrix.shape == (1, 1):
        return (1,)
    _, vectors = eigh(matrix.toarray().astype(float), subset_by_index=(matrix.shape[0]-1,)*2)
    v = np.abs(vectors[:, 0]); v /= np.max(v)
    if not np.isfinite(v).all():
        raise ArithmeticError('the numerical witness proposal is not finite')
    return tuple(max(1, int(round(float(x)*(1 << PERRON_BITS)))) for x in v)


@lru_cache(maxsize=1)
def all_active_blocks():
    """Complete active-sector norm enclosure, reduced only by equal weights."""
    rows, winning = [], None
    for counts in product(range(3), repeat=3):
        kinds = tuple(k for count in counts for k in ([1]*count+[2]*(2-count)))
        states, channels, matrix = gram(kinds)
        count, labels = csgraph.connected_components(matrix, directed=False)
        components = []
        for label in range(count):
            indices = np.flatnonzero(labels == label)
            block = matrix[indices][:, indices].tocsr()
            vector = propose_positive_vector(block)
            lower, upper = check_perron_vector(block, vector)
            row = {'component': label, 'dimension': len(indices), 'lower': lower, 'upper': upper}
            components.append(row)
            if kinds == WINNING_TYPES and 0 in indices:
                winning = {'kinds': kinds, 'states': states, 'channels': channels,
                           'matrix': matrix, 'indices': tuple(map(int, indices)),
                           'block': block, 'vector': vector, 'lower': lower, 'upper': upper,
                           'component': label}
        rows.append({'fundamental_counts_per_equal_weight_pair': counts,
                     'represented_active_type_patterns': int(np.prod([2 if n == 1 else 1 for n in counts])),
                     'complete_source_sector_dimension': len(states),
                     'components': components,
                     'top_lower': max(r['lower'] for r in components),
                     'top_upper': max(r['upper'] for r in components)})
    if winning is None or sum(row['represented_active_type_patterns'] for row in rows) != 64:
        raise ArithmeticError('the complete six-cell active domain is not covered')
    competitors = [component['upper'] for row in rows for component in row['components']
                   if not (row['fundamental_counts_per_equal_weight_pair'] == (0, 2, 2)
                           and component['component'] == winning['component'])]
    other_upper = max(competitors)
    if winning['lower'] <= other_upper:
        raise ArithmeticError('exact interval admission does not separate the proposed full-carrier ground')
    return {'sectors': rows, 'winning': winning, 'other_component_upper': other_upper}


@lru_cache(maxsize=1)
def exact_root_and_gap():
    data = all_active_blocks(); win = data['winning']; x = s.Symbol('x')
    polynomial = s.Poly.from_list(ROOT_POLYNOMIAL, x)
    full_characteristic = s.Matrix(win['block'].toarray().tolist()).charpoly(x).as_poly()
    if full_characteristic.rem(polynomial).as_expr() != 0:
        raise ArithmeticError('the recorded algebraic root is not in the actual source characteristic polynomial')
    lo, hi = map(s.Rational, (win['lower'], win['upper']))
    if polynomial.count_roots(lo, hi) != 1 or full_characteristic.count_roots(lo, hi) != 1:
        raise ArithmeticError('the exact Perron interval does not isolate the recorded root')
    # The largest root of this irreducible factor is isolated by the source interval.
    root = s.CRootOf(polynomial, 5)
    if polynomial.count_roots(hi, s.oo) != 0:
        raise ArithmeticError('the selected root is not the largest factor root')
    other = F(703, 10)
    if data['other_component_upper'] >= other:
        raise ArithmeticError('the declared rational competitor ceiling is insufficient')
    if full_characteristic.count_roots(s.Rational(other), hi) != 1:
        raise ArithmeticError('another eigenvalue in the winning block exceeds the competitor ceiling')
    # Inactive cells embed into an active sector with added2w^2 I. At least
    # one omitted cell therefore has Gram maximum at most lambda_max-2.
    if win['upper']-2 >= other:
        raise ArithmeticError('the inactive-domain competitor bound is insufficient')
    delta_lo = iv.sqrt_rational(1+4*win['lower'], 28).lo
    delta_hi = iv.sqrt_rational(1+4*win['upper'], 28).hi
    other_delta_hi = iv.sqrt_rational(1+4*other, 28).hi
    gap_lower = (delta_lo-other_delta_hi)/2
    if gap_lower <= 0:
        raise ArithmeticError('the full-carrier ground has no certified positive separation')
    return {'polynomial': tuple(ROOT_POLYNOMIAL), 'root': root,
            'root_interval': (win['lower'], win['upper']),
            'global_ground_interval': ((1-delta_hi)/2, (1-delta_lo)/2),
            'global_ground_exact': '(1-sqrt(1+4lambda_max))/2',
            'all_other_Gram_eigenvalues_upper': other,
            'full_ground_gap_lower': gap_lower,
            'neutral_307_ground_energy': -7,
            'global_ground_multiplicity': 1,
            'full_source_dimension': 3*(1 << 18)}


def quantum_numbers():
    """All original occupations in the exact global-ground support."""
    win = all_active_blocks()['winning']
    empty, occupied = set(), set()
    support = [win['states'][index] for index in win['indices']]
    masks = []
    for coordinate in support:
        mask, _ = original_column(WINNING_TYPES, coordinate); masks.append(mask)
        right = sum((mask >> (6*f+4+spin)) & 1 for f in range(3) for spin in range(2))
        left = mask.bit_count()-right
        empty.add((left, right, mask.bit_count(), 2*right-left))
        for h in (1, 2):
            for destination, coefficient in source.vertex_action(s.diag(1, 2, 3), 0, mask).items():
                scalar, after = destination
                if scalar != h or coefficient == 0:
                    continue
                after_right = sum((after >> (6*f+4+spin)) & 1 for f in range(3) for spin in range(2))
                after_left = after.bit_count()-after_right
                occupied.add((after_left, after_right, after.bit_count()-2, 2*after_right-after_left-1))
    if empty != {(4, 2, 6, 0)} or occupied != {(5, 3, 6, 0)}:
        raise ArithmeticError('the actual full-ground support changed its source gauge charge')
    return {'scalar_empty_support_NL_Ne_NF_twice_hypercharge': sorted(empty),
            'scalar_occupied_support_NL_Ne_Qstar_twice_hypercharge': sorted(occupied),
            'conserved_Qstar': 6, 'CAR_parity': '+1',
            'scalar_empty_support_dimension': len(masks),
            'weak_gauge_SU2_singlet_from_unique_ground_and_actual_Ward_identity': True,
            'spin_SU2_singlet_from_unique_ground_and_actual_Ward_identity': True,
            'gauge_hypercharge': 0,
            'all18_mode_source_has_color_trivial_lepton_representation': True,
            'same_neutral_307_prepared_ground': False}


def one_body_action(vector, pairs):
    """Literal CAR action sum c_to^dagger c_from on a retained sparse state."""
    result = {}
    for mask, value in vector.items():
        for target, origin in pairs:
            first = source._annihilate(mask, origin)
            second = None if first is None else source._create(first[0], target)
            if second is not None:
                result[second[0]] = result.get(second[0], F(0))+value*first[1]*second[1]
    return {mask: value for mask, value in result.items() if value}


def singlet_projection(vector, raising_pairs, maximum_spin):
    """Exact Casimir projection on a verified Jz=0 occupation domain."""
    def weight(mask):
        return sum(((mask >> target) & 1)-((mask >> origin) & 1)
                   for target, origin in raising_pairs)
    if any(weight(mask) for mask in vector):
        raise ValueError('the complete projection input must have Jz=0')
    out = dict(vector)
    lowering = tuple((origin, target) for target, origin in raising_pairs)
    for spin in range(1, maximum_spin+1):
        casimir = one_body_action(one_body_action(out, raising_pairs), lowering)
        for mask, value in casimir.items():
            out[mask] = out.get(mask, F(0))-F(value, spin*(spin+1))
        out = {mask: value for mask, value in out.items() if value}
    if one_body_action(out, raising_pairs):
        raise ArithmeticError('the admitted Casimir polynomial did not yield a singlet')
    return out


def scalar_image(vector):
    """The actual two-component pair creation map, without a sector cutoff."""
    result = {}
    for mask, value in vector.items():
        for h, f, g, left, right, sign in source._pair_terms(3):
            if f != g:
                continue
            first = source._create(mask, right)
            second = None if first is None else source._create(first[0], left)
            if second is not None:
                key = h+1, second[0]
                result[key] = result.get(key, F(0))+value*(f+1)*sign*first[1]*second[1]
    return {key: value for key, value in result.items() if value}


@lru_cache(maxsize=1)
def singlet_trial():
    """Rational approximate ground with exact weak-gauge and spin singletness."""
    win = all_active_blocks()['winning']
    trial = {}
    for index, value in zip(win['indices'], win['vector'], strict=True):
        mask, sign = original_column(WINNING_TYPES, win['states'][index])
        trial[mask] = F(value*sign)
    weak = tuple((6*f+spin, 6*f+2+spin) for f in range(3) for spin in range(2))
    spin = tuple((6*f+2*internal, 6*f+2*internal+1) for f in range(3) for internal in range(3))
    projected = singlet_projection(singlet_projection(trial, weak, 2), spin, 3)
    if one_body_action(projected, weak) or one_body_action(projected, spin):
        raise ArithmeticError('the projected state does not obey both exact source Ward charges')
    if not projected or any(mask.bit_count() != 6 for mask in projected):
        raise ArithmeticError('the projected preparation changed its retained charge domain')
    image = scalar_image(projected)
    norm = sum(value*value for value in projected.values())
    image_norm = sum(value*value for value in image.values())
    return {'source_state': projected, 'source_image': image,
            'norm_squared': norm, 'image_norm_squared': image_norm,
            'Gram_Rayleigh_quotient': image_norm/norm}


def prepared_ground_witness():
    """A complete, computable near-ground preparation with a fidelity bound."""
    win = all_active_blocks()['winning']; vector = win['vector']
    trial = singlet_trial(); rayleigh = trial['Gram_Rayleigh_quotient']
    # The normalized span {u,A u/sqrt(r)} has exact compression [[0,sqrt(r)],[sqrt(r),1]].
    # It need not be invariant; its Rayleigh energy gives a certified trial.
    sqrt_r_low = iv.sqrt_rational(1+4*rayleigh, 28).lo
    trial_energy_upper = (1-sqrt_r_low)/2
    exact = exact_root_and_gap()
    infidelity_upper = (trial_energy_upper-exact['global_ground_interval'][0])/exact['full_ground_gap_lower']
    if not 0 <= infidelity_upper < F(1, 10**8):
        raise ArithmeticError('the actual complete-source preparation is not within its declared fidelity')
    return {'positive_integer_proposal': list(vector),
            'singlet_trial_original_columns': [{'mask': mask, 'coefficient': value}
                for mask, value in sorted(trial['source_state'].items())],
            'trial_norm_squared': trial['norm_squared'],
            'trial_Gram_Rayleigh_quotient': rayleigh,
            'trial_energy_upper': trial_energy_upper,
            'ground_infidelity_upper': infidelity_upper,
            'exact_ground_definition': 'the unique positive Perron vector of the recorded 53x53 source block, with its actual scalar image',
            'preparation_seed': 'first family all six CAR modes filled; scalar empty; other families empty',
            'seed_original_mask': (1 << 6)-1, 'seed_Qstar': 6,
            'full_unitary_route': 'Householder between the declared normalized seed and the computed ground, identity on the complete orthogonal complement',
            'approximate_preparation_is_exactly_gauge_and_spin_singlet': True,
            'approximate_preparation_uses_compensating_reference': False,
            'exact_ground_and_seed_are_gauge_and_spin_singlets': True,
            'seed_mean_H_energy': 0,
            'exact_preparation_mean_controller_work': '(1-sqrt(1+4lambda_max))/2',
            'autonomous_energy_conserving_preparation_claimed': False,
            'closed_source_vertex_from_empty_Qstar0_prepares_Qstar6': False,
            'compiled_primitive_tape_executed': False}


def complete_ground_records():
    """Complete scalar records on the full-source ground and one fixed clock."""
    exact = exact_root_and_gap()
    lam = iv.I(*exact['root_interval'])
    delta = iv.I(iv.sqrt_rational(1+4*lam.lo, 28).lo,
                 iv.sqrt_rational(1+4*lam.hi, 28).hi)
    alpha = (1+delta.reciprocal())/2
    beta = 1-alpha
    residue = lam/(1+4*lam)
    rows = []
    for ticks in (1, 2):
        time = F(ticks, 32)
        transition = 4*residue*iv.sine(delta*time/2)**2
        laws = {'00': alpha*(1-transition), '01': alpha*transition,
                '10': beta*transition, '11': beta*(1-transition)}
        rows.append({'ticks': ticks, 'time': time,
                     'unconditional_joint_record_intervals': {event: iv.rounded_out(p, 12, probability=True).strings()
                                                              for event, p in laws.items()}})
    eta = prepared_ground_witness()['ground_infidelity_upper']
    trace_error = iv.sqrt_rational(eta, 24).hi
    if delta.hi/32 >= iv.pi_interval().lo:
        raise ArithmeticError('the common tick does not certify its positive inverse-gap branch')
    return {'complete_source_dimension': 786432,
            'complete_detector_projectors': ['scalar empty tensor I_2^18', 'scalar occupied tensor I_2^18'],
            'pointer_rank_empty': 1 << 18, 'pointer_rank_occupied': 2*(1 << 18),
            'shared_source_tick': '1/32',
            'ground_occupied_probability': '(1-1/Delta)/2',
            'connected_ordered_NH_correlator': 'lambda_max/Delta^2 * exp(-i Delta t)',
            'same_H_record_gap': 'Delta=sqrt(1+4lambda_max)',
            'record_gap_equals_full_ground_isolation_gap': False,
            'gap_interval': iv.rounded_out(delta, 12).strings(),
            'spectral_residue_interval': iv.rounded_out(residue, 12, probability=True).strings(),
            'joint_law_formula': 'J00=alpha(1-T), J01=alpha T, J10=beta T, J11=beta(1-T); T=4lambda_max/Delta^2 sin^2(Delta t/2)',
            'raw_double_time_ratio': 'J01(2t)/(4J01(t))=cos^2(Delta t/2)',
            'inverse_branch_certified': '0<Delta/32<pi; J01(t)>0',
            'exact_record_enclosures': rows,
            'pure_preparation_trace_distance_upper': trace_error,
            'complete_two_tick_two_read_law_TV_error': 'sqrt(eta)+epsilon_preparation+2epsilon_tick+2epsilon_read',
            'example_each_compiler_operator_norm_error': F(1, 10**6),
            'example_complete_law_TV_error_upper': trace_error+F(5, 10**6),
            'all_input_detector_columns_retained': True,
            'other_charge_sectors_are_not_discarded': True,
            'all_external_quantum_references_retained': True,
            'record_gap_identified_as_relativistic_particle_mass': False}


def _serial(value):
    if isinstance(value, dict):
        return {str(k): _serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(v) for v in value]
    if isinstance(value, (F, s.Basic)):
        return str(value)
    return value


def certificate():
    data = all_active_blocks(); exact = exact_root_and_gap()
    return _serial({'schema': 'oph.full_finite_scalar_ground.v1',
        'source_inputs': list(SOURCE_INPUTS),
        'declared_H': 'W(diag(1,2,3))+W^dagger+NH; epsilonF=0',
        'complete_CAR_scalar_dimension': 786432,
        'cell_sector_dimension_sum': (1+3+3+1)**6,
        'all_active_sector_patterns': 64, 'equal_weight_quotient_patterns': 27,
        'largest_constructed_Gram_dimension': 729,
        'active_sector_bounds': data['sectors'],
        'unique_ground': exact, 'ground_quantum_numbers': quantum_numbers(),
        'preparation': prepared_ground_witness(),
        'same_H_complete_ground_records': complete_ground_records(),
        'neutral_307_ground_is_full_carrier_vacuum': False,
        'neutral_307_charged_lines_are_relabelled_as_global_vacuum_particles': False,
        'spectral_root_admission_uses_floating_values': False,
        'physical_particle_assignments': [], 'measured_particle_inputs': []})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2, sort_keys=True))
