"""Centered-source marginal enclosures and complete native one-site records."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from .defect_budget import exact_source

ROOT = Path(__file__).resolve().parents[2]
NATIVE_INPUTS = (
    'model/centered_native_preparation/results.json',
    'model/centered_native_preparation/results.npz',
    'model/centered_particle_spectrum/witness.npz',
)


def _rational(x, *, positive=False):
    if type(x) not in (int, F):
        raise ValueError('an exact rational is required')
    x = F(x)
    if positive and x <= 0:
        raise ValueError('a positive exact rational is required')
    return x


def _interval(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError('a complete two-endpoint interval is required')
    lo, hi = map(_rational, value)
    if lo > hi:
        raise ValueError('reversed interval')
    return lo, hi


def energy_density_interval(side, coupling, finite_energy, *,
                            field=F(), observable_norm=F(1), spans=(0, 0, 0)):
    """Two-sided tiling bound, retaining every omitted perturbation translate."""
    if (type(side) is not int or side < 2 or side % 2
            or type(spans) is not tuple or len(spans) != 3
            or any(type(x) is not int or not 0 <= x < side for x in spans)):
        raise ValueError('an even cube side and contained coordinate spans are required')
    c = _rational(coupling)
    if c < 0:
        raise ValueError('a nonnegative hopping coefficient is required')
    h, a = _rational(field), _rational(observable_norm, positive=True)
    lo, hi = _interval(finite_energy)
    contained = F(1)
    for span in spans:
        contained *= F(side-span, side)
    omitted = 1-contained
    surface = 3*c/side
    perturbation = abs(h)*a*omitted
    return {'interval': (lo/side**3-surface-perturbation,
                         hi/side**3+perturbation),
            'missing_bond_norm_per_site': surface,
            'omitted_perturbation_translate_fraction': omitted,
            'onsite_terms_removed': False,
            'finite_ground_substituted_for_infinite_marginal': False}


def expectation_interval(side, coupling, step, energy_zero,
                         energy_plus, energy_minus, *,
                         observable_norm=F(1), spans=(0, 0, 0)):
    """Rigorous pressure secants for every translation-invariant source ground."""
    h = _rational(step, positive=True)
    a = _rational(observable_norm, positive=True)
    zero = energy_density_interval(side, coupling, energy_zero)
    plus = energy_density_interval(side, coupling, energy_plus,
                                   field=h, observable_norm=a, spans=spans)
    minus = energy_density_interval(side, coupling, energy_minus,
                                    field=-h, observable_norm=a, spans=spans)
    upper_zero = zero['interval'][1]
    lo = max(-a, (plus['interval'][0]-upper_zero)/h)
    hi = min(a, (upper_zero-minus['interval'][0])/h)
    if lo > hi:
        raise ArithmeticError('energy enclosures contradict the variational bound')
    return {'expectation_interval': (lo, hi),
            'ground_energy_density_interval': zero['interval'],
            'plus_energy_density_interval': plus['interval'],
            'minus_energy_density_interval': minus['interval'],
            'all_translation_invariant_ground_marginals_covered': True,
            'unique_marginal_or_convergence_modulus_inferred': False}


def pauli_marginal_bound(qubits, reference, intervals):
    """Squared half-trace bound from every Pauli coefficient, identity first."""
    if (type(qubits) is not int or qubits < 1
            or type(reference) is not tuple or type(intervals) is not tuple
            or len(reference) != 4**qubits or len(intervals) != len(reference)):
        raise ValueError('complete Pauli tables are required')
    ref = tuple(_rational(x) for x in reference)
    boxes = tuple(_interval(x) for x in intervals)
    if ref[0] != 1 or boxes[0] != (F(1), F(1)):
        raise ValueError('the identity coefficient must be retained exactly')
    if any(abs(x) > 1 for x in ref) or any(lo < -1 or hi > 1 for lo, hi in boxes):
        raise ValueError('Pauli coefficients must be in [-1,1]')
    squared = sum((max(abs(x-lo), abs(x-hi))**2
                   for x, (lo, hi) in zip(ref[1:], boxes[1:], strict=True)), F())/4
    return {'half_trace_distance_squared_upper': min(F(1), squared),
            'complete_Pauli_coefficients': len(ref),
            'reference_density_positivity_is_a_separate_input': True,
            'conditional_success_normalization_used': False}


def complete_record_interval(ideal_interval, marginal_halftrace,
                             preparation_error, detector_error):
    """Compose complete-state and complete-instrument errors before selection."""
    lo, hi = _interval(ideal_interval)
    errors = tuple(_rational(x) for x in
                   (marginal_halftrace, preparation_error, detector_error))
    if not 0 <= lo <= hi <= 1 or any(x < 0 for x in errors):
        raise ValueError('probabilities and nonnegative errors are required')
    error = sum(errors, F())
    return max(F(), lo-error), min(F(1), hi+error)


def default_faithfulness():
    """Exact rational consequences of the actual cubic dimer and hopping bounds."""
    c, u = F(1, 8), F(1)
    root_lower = F(103, 25)
    kinetic_root_upper = F(49, 40)
    if root_lower**2 >= 17 or kinetic_root_upper**2 <= F(3, 2):
        raise ArithmeticError('rational radical bounds failed')
    delta_lower = (root_lower-4)/16
    lower = delta_lower**2/(324*c*c)
    upper = (c*kinetic_root_upper-delta_lower)/u
    if lower != F(1, 90000) or upper != F(233, 1600):
        raise ArithmeticError('source one-site faithfulness constants changed')
    if F(1, 2)-upper <= lower:
        raise ArithmeticError('single-occupation marginal lost its lower bound')
    return {'coupling': c, 'contact': u,
            'doublon_density_strict_interval': (lower, upper),
            'complete_one_site_density': ['d', '1/2-d', '1/2-d', 'd'],
            'least_one_site_eigenvalue_strict_lower': lower,
            'dimer_energy_gain_strict_lower': delta_lower,
            'thermodynamic_gap_assumed': False,
            'larger_region_faithfulness_inferred': False,
            'unique_thermodynamic_ground_inferred': False}


def algebraic_ground_projector(matrix):
    """Exact finite rational source oracle, including degeneracies and all sectors.

    The executable controls use at most sixteen rows. The polynomial proof
    applies at every finite dimension; no practical runtime claim is made.
    """
    matrix = s.ImmutableMatrix(matrix)
    if matrix.rows != matrix.cols or not 1 <= matrix.rows <= 16 or matrix != matrix.H:
        raise ValueError('a complete Hermitian source of at most sixteen rows is required')
    if any(not x.is_Rational for x in matrix):
        raise ValueError('this finite oracle uses rational signed-frame entries')
    z = s.Symbol('z')
    q = matrix.charpoly(z).as_poly().sqf_part()
    roots = q.real_roots()
    if len(roots) != q.degree():
        raise ArithmeticError('Hermitian source has a nonreal characteristic root')
    ground = roots[0]
    quotient, remainder = s.div(q.as_expr(), z-ground, z)
    if s.simplify(remainder) != 0:
        raise ArithmeticError('algebraic root division failed')
    polynomial = s.Poly(quotient, z)
    value = s.zeros(matrix.rows)
    for coefficient in polynomial.all_coeffs():
        value = value*matrix+coefficient*s.eye(matrix.rows)
    value = (value/q.diff().eval(ground)).applyfunc(s.simplify)
    rank = s.simplify(s.trace(value))
    if (not rank.is_Integer or rank <= 0 or value != value.H
            or (value*value-value).applyfunc(s.simplify) != s.zeros(matrix.rows)
            or (matrix*value-ground*value).applyfunc(s.simplify) != s.zeros(matrix.rows)):
        raise ArithmeticError('full finite ground projector identity failed')
    return ground, s.ImmutableMatrix(value/rank), int(rank)


@lru_cache(maxsize=1)
def _finite_source_controls():
    c, u = F(1, 8), F(1)
    q, hopping, defect = exact_source(2, ((0, 1, 1),), c, u)
    ground, rho, rank = algebraic_ground_projector(q)
    expected = (4-s.sqrt(17))/8
    if s.simplify(ground-expected) != 0 or rank != 1:
        raise ArithmeticError('complete source dimer ground differs')
    single = s.diag(*(int(((mask & 3).bit_count() == 1)
                         and (((mask >> 2) & 3).bit_count() == 1))
                      for mask in range(16)))
    if single*hopping*single != s.zeros(16):
        raise ArithmeticError('exact single-occupation hopping compression failed')
    if max(abs(x) for x in hopping.eigenvals()) != s.Rational(1, 8):
        raise ArithmeticError('whole-Fock bond norm differs')
    # The one-particle block is extracted by occupation addresses, independently
    # of the whole-Fock diagonalization and its selected ground.
    indices = [1 << j for j in range(4)]
    one = hopping.extract(indices, indices)
    if s.trace(one) != 0 or s.trace(one*one) != s.Rational(1, 64):
        raise ArithmeticError('both-mode kinetic trace-square factor changed')
    d = s.simplify(sum(rho[j, j] for j in range(16) if j & 3 == 3))
    if s.simplify(d-(1-4/s.sqrt(17))/4) != 0:
        raise ArithmeticError('actual dimer local density changed')
    if s.simplify(s.trace(rho*defect)/4-d) != 0:
        raise ArithmeticError('onsite energy Ud factor changed')
    return {'source_sites': 2, 'full_Fock_rows': 16,
            'absolute_ground_multiplicity': rank,
            'ground_energy': str(ground), 'doublon_density': str(d),
            'all_Fock_sectors_used': True,
            'whole_bond_norm': c, 'one_particle_trace_squared': F(1, 64),
            'single_occupation_hopping_compression_exactly_zero': True}


def finite_source_controls():
    return {key: F(value.numerator, value.denominator) if type(value) is F else value
            for key, value in _finite_source_controls().items()}


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def _native_cube_reference():
    """Authenticate the existing executed preparation, then replay its exact target."""
    from model.centered_native_preparation import target
    path = ROOT/NATIVE_INPUTS[0]
    receipt = strict_json.load(path)
    if not all(receipt.get(key) is True for key in
               ('all_pass', 'new_checks_pass', 'dependency_custody_pass')):
        raise ArithmeticError('the existing native preparation is not admitted')
    for name, expected in receipt['source_bindings'].items():
        if _hash(ROOT/name) != expected:
            raise ArithmeticError('native preparation source binding differs: '+name)
    archive = path.parent/receipt['complete_record_archive']['filename']
    if (archive != ROOT/NATIVE_INPUTS[1]
            or _hash(archive) != receipt['complete_record_archive']['sha256']):
        raise ArithmeticError('complete native preparation archive differs')
    prep = receipt['complete_error']['preparation']
    if (target.summary() != prep['target'] or prep['prepared_number'] != 8
            or prep['postselection_required'] is not False
            or not F(prep['complete_preparation_error']) < F(1, 10**7)):
        raise ArithmeticError('actual target or complete preparation bound differs')
    values, norm, _ = target.ground()
    density = tuple(F(sum(value*value for mask, value in enumerate(values)
                          if mask & 3 == state), norm) for state in range(4))
    d = density[0]
    if density != (d, F(1, 2)-d, F(1, 2)-d, d):
        raise ArithmeticError('certified full target one-site symmetry differs')
    return receipt, archive, density


@lru_cache(maxsize=1)
def _native_local_records():
    """Execute six native tagging pulses on the archived unnormalized full state."""
    from model.certified_reference import certificate as cert
    from model.native_spectral_detector import circuits, bounds
    from model.protected_particle_readout import execute
    receipt, archive, density = _native_cube_reference()
    with np.load(archive, allow_pickle=False) as data:
        state = np.array(data['prepared_ground'], copy=True)
        earlier_failure = float(data['tree_failure_weights'].sum()
                                +data['frame_failure_weights'].sum())
    if state.shape != (65536, 1):
        raise ArithmeticError('native full-state archive shape changed')
    template = cert.cnot(0, 1)
    word = circuits.emit(template)
    block = execute.compile_block(execute.key(word), 2)
    initial = np.zeros((131072, 1), complex)
    initial[:65536] = state
    first = execute.complete(((block, (0, 16)),), 17, initial)
    # Destructively commit the first tag and supply a second blank register.
    # Both proper first-record columns remain separate; no reset is inferred
    # from measurement. The other four raw pointer values have zero maps.
    columns = first.success.reshape(2, 65536).T
    second_input = np.zeros((131072, 2), complex)
    second_input[:65536] = columns
    second = execute.complete(((block, (1, 16)),), 17, second_input)
    final = second.success.reshape(2, 65536, 2)
    observed = np.sum(np.abs(final)**2, axis=1)
    if abs(float(observed.sum())+first.failure_probability
           +second.failure_probability+earlier_failure-1) >= 10**-8:
        raise ArithmeticError('complete original preparation/record accounting failed')
    detector_error = 2*bounds.local(word, template)
    error = F(receipt['complete_error']['preparation']['complete_preparation_error'])+detector_error
    if error >= F(1, 10**6):
        raise ArithmeticError('complete native local record error missed its enclosure')
    if max(abs(float(observed[b, a])-float(density[a+2*b]))
           for a in range(2) for b in range(2)) >= 10**-6:
        raise ArithmeticError('native local record diagnostic differs from exact target')
    # Stable outward intervals are accepted only after enclosing the complete
    # source-state error; numeric execution is an independent diagnostic.
    low, high = F(15954, 10**6), F(15958, 10**6)
    if not low < density[0]-error <= density[0]+error < high:
        raise ArithmeticError('fixed local occupation interval is not enclosing')
    diagonal = ((low, high), (F(1, 2)-high, F(1, 2)-low),
                (F(1, 2)-high, F(1, 2)-low), (low, high))
    records = []
    for a in range(6):
        for b in range(6):
            interval = diagonal[a+2*b] if a < 2 and b < 2 else (F(), F())
            records.append({'raw_pointer_records': (a, b),
                            'original_probability_interval': interval})
    # Reconstruct every second-stage first-failure branch separately for its
    # earlier classical label; aggregating amplitudes would erase that record.
    second_failures = []
    retained_failures = [
        {'stage': 'first_tag', 'event': event, 'raw_export_slot': slot,
         'original_probability_interval': (F(), F(1, 10**12))}
        for event in range(3) for slot in range(4)]
    for a in range(6):
        for event in range(3):
            for slot in range(4):
                if a < 2 and slot < len(block.failures[event]):
                    value = second.failure(0, event, slot)[:, a]
                    weight = float(np.vdot(value, value).real)
                else:
                    weight = 0.
                second_failures.append(weight)
                retained_failures.append({
                    'stage': 'second_tag', 'first_pointer_record': a,
                    'event': event, 'raw_export_slot': slot,
                    'original_probability_interval':
                    (F(), F(1, 10**12)) if a < 2 else (F(), F())})
    if abs(sum(second_failures)-second.failure_probability) >= 10**-20:
        raise ArithmeticError('fine classical first-failure records were lost')
    return {'existing_preparation_receipt_sha256': _hash(ROOT/'model/centered_native_preparation/results.json'),
            'existing_preparation_archive_sha256': _hash(archive),
            'existing_source_bindings_authenticated': len(receipt['source_bindings']),
            'existing_preparation_tape_reexecuted': False,
            'additional_native_tagging_pulses_executed': len(word)*2,
            'fresh_blank_tag_registers_consumed': 2,
            'pointer_measurement_assumed_to_reset_the_register': False,
            'full_field_rows_retained': 65536,
            'complete_pointer_records': records,
            'complete_additional_first_failure_records': retained_failures,
            'first_tag_failure_slots': 12, 'second_tag_failure_slots_with_first_record': 72,
            'earlier_preparation_failure_archive_retained': True,
            'complete_preparation_and_tagging_error_strict_upper': F(1, 10**6),
            'native_record_probabilities_renormalized': False,
            'all_four_nonzero_local_records_strict_lower': F(3, 200),
            'thermodynamic_reference_prepared': False,
            'growing_preparation_tapes_executed': False,
            'physical_particle_or_mass_assigned': False}


def cube_energy_density_interval():
    receipt, _, _ = _native_cube_reference()
    exact = tuple(F(x)+4 for x in receipt['complete_error']['preparation']['target']
                  ['ground_energy_interval'])
    fixed = (F(-6203, 50000), F(-2481, 20000))
    if not fixed[0] < exact[0] <= exact[1] < fixed[1]:
        raise ArithmeticError('the fixed full-cube energy enclosure failed')
    return {'source_cube_side': 2, 'full_Fock_rows': 65536,
            'scalar_shift_K_to_Q': F(4),
            'finite_Q_ground_energy_interval': fixed,
            **energy_density_interval(2, F(1, 8), fixed)}


def one_site_prepared_reference_bound():
    """A uniform finite-to-infinite local bound, without a convergence assertion."""
    receipt, _, density = _native_cube_reference()
    prep = receipt['complete_error']['preparation']
    target_error = F(prep['target']['normalized_target_to_exact_ground_error'])
    cube = (F(15955, 10**6), F(15957, 10**6))
    if not cube[0] < density[0]-target_error <= density[0]+target_error < cube[1]:
        raise ArithmeticError('exact cube marginal is outside its fixed interval')
    lower, upper = default_faithfulness()['doublon_density_strict_interval']
    ideal_distance = 2*max(abs(cube[0]-upper), abs(cube[1]-lower))
    complete_distance = ideal_distance+F(prep['complete_preparation_error'])
    if complete_distance >= F(13, 50):
        raise ArithmeticError('prepared one-site comparison lost its strict bound')
    return {'ideal_cube_doublon_interval': cube,
            'ideal_cube_to_symmetric_ground_halftrace_upper': ideal_distance,
            'complete_prepared_one_site_halftrace_strict_upper': F(13, 50),
            'every_selected_symmetric_centered_ground_covered': True,
            'full_source_or_multisite_state_bound': False,
            'arbitrarily_small_error_modulus_inferred': False}


def _encode(x):
    if type(x) is F:
        return str(x)
    if type(x) is dict:
        return {k: _encode(v) for k, v in x.items()}
    if type(x) in (list, tuple):
        return [_encode(v) for v in x]
    return x


def certificate():
    return json.loads(json.dumps(_encode({
        'uniform_one_site_faithfulness': default_faithfulness(),
        'whole_Fock_source_controls': finite_source_controls(),
        'executed_native_local_records': _native_local_records(),
        'actual_cube_energy_density_enclosure': cube_energy_density_interval(),
        'prepared_one_site_thermodynamic_comparison': one_site_prepared_reference_bound(),
        'thermodynamic_local_marginal_convergence_modulus_proved': False,
        'cut_ground_substituted_for_infinite_ground': False,
        'physical_state_selection_proved': False})))
