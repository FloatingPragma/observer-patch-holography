"""Gauge-coordinate cutoff bridge and complete spatial-overlap rotor records.

The coordinate completion and Hamiltonian are specified construction choices.
The executed control is one cell and one abelian native charge. The general
cutoff estimate applies to a finite cubic gauge-coordinate/Fock program; it
does not establish a chiral or interacting continuum limit.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import json
import math

import sympy as s

from model.interacting_extension.chiral_lattice import source_data
from model.particle_candidate_atlas.interacting_chiral import creator


SOURCE_INPUTS = (
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/interacting_extension/chiral_lattice/rotor_response.py',
    'model/interacting_extension/chiral_lattice/PROOF.md',
    'model/interacting_extension/chiral_electric/PROOF.md',
    'model/common_source/overlap_regulator.py',
    'model/particle_candidate_atlas/interacting_chiral.py',
)


def _int(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(name+' must be an integer at least '+str(minimum))
    return value


def _frac(value, name, *, positive=False):
    if type(value) not in (F, int, str):
        raise ValueError(name+' must be exact rational data')
    value = F(value)
    if value < 0 or positive and value == 0:
        raise ValueError(name+' has the wrong sign')
    return value


def _clean(matrix):
    matrix = s.ImmutableSparseMatrix(matrix)
    return s.ImmutableSparseMatrix(matrix.rows, matrix.cols,
        {k: y for k, x in matrix.todok().items()
         if (y := s.simplify(s.expand(x))) != 0})


def polynomial_bounds(degree, spacing='1'):
    """Same Legendre polynomial as the retained overlap implementation.

    The polynomial is globally bounded even outside the admitted gap. Its
    accuracy relative to sign(H) requires |H|>=1/2; ||H||<=8 is supplied by
    the finite Wilson assembly, independent of the gauge configuration.
    """
    degree = _int(degree, 'Legendre degree')
    a = _frac(spacing, 'spacing', positive=True)
    epsilon = 16*F(15, 17)**(degree+1)
    return {'degree': degree, 'maximum_coordinate_word_degree': 2*degree+1,
            'sign_error_on_admitted_gap': epsilon,
            'global_sign_polynomial_norm_upper': 1+epsilon,
            'global_one_body_H_norm_upper': (2+epsilon)/a,
            'gap_required_for_accuracy': '1/2', 'Wilson_norm_upper': '8',
            'uniform_global_sign_accuracy_claimed': False}


def sign_polynomial(kernel, degree):
    """Exact finite-matrix evaluation; caller supplies the Wilson kernel."""
    degree = _int(degree, 'Legendre degree')
    h = _clean(kernel)
    if h.rows != h.cols or not _clean(h-h.H).is_zero_matrix:
        raise ValueError('a Hermitian square kernel is required')
    identity = s.eye(h.rows)
    z = (257*identity-8*h*h)/255
    previous = identity
    value = s.Rational(4, 17)*previous
    if degree:
        current = z
        value += s.Rational(4, 17)*s.Rational(15, 17)*current
        for j in range(1, degree):
            following = _clean(((2*j+1)*z*current-j*previous)/(j+1))
            value += s.Rational(4, 17)*s.Rational(15, 17)**(j+1)*following
            previous, current = current, following
    return _clean(h*value)


def exponential_tail_upper(x, first, extra=40):
    """Rational upper enclosure of sum_(n>=first) x^n/n!."""
    x = _frac(x, 'integrated interaction norm')
    first = _int(first, 'first tail index', 1)
    extra = _int(extra, 'tail summation terms', 1)
    if x == 0:
        return F(0)
    last = max(first+extra, int(x)+2)
    term = x**first/math.factorial(first)
    total = term
    for n in range(first+1, last+1):
        term *= x/n
        total += term
    next_term = term*x/(last+1)
    ratio = x/(last+2)
    if ratio >= 1:
        raise ArithmeticError('geometric factorial-tail majorant is not contracting')
    return total+next_term/(1-ratio)


def cutoff_bound(cutoff, initial_degree, coordinate_step, integrated_norm):
    """Complete-history bound with all occupation and pointer outcomes.

    Electric Casimirs and the stated pointer readouts must preserve the
    degree filtration. Bounded interaction and pointer-coupling insertions
    raise its total degree by at most coordinate_step.
    """
    k = _int(cutoff, 'cutoff')
    k0 = _int(initial_degree, 'initial degree')
    step = _int(coordinate_step, 'coordinate step', 1)
    if k0 > k:
        raise ValueError('initial support lies outside the cutoff')
    x = _frac(integrated_norm, 'integrated interaction norm')
    first = (k-k0)//step+1
    tail = exponential_tail_upper(x, first)
    return {'cutoff': k, 'initial_degree': k0, 'coordinate_step': step,
            'first_possible_omitted_Dyson_order': first,
            'integrated_interaction_norm_upper': str(x),
            'state_vector_error_upper': str(min(F(2), 2*tail)),
            'complete_record_TV_error_upper': str(min(F(1), 2*tail)),
            'all_histories_including_zero_outputs_retained': True,
            'valid_for_arbitrary_untouched_references': True,
            'initial_ground_state_approximation_included': False}


def cubic_program_bound(side, degree, cutoff, initial_degree, time,
                        spacing='1', plaquette_coefficient='0', onsite_norm='0',
                        integrated_pointer_norm='0'):
    """Full M15 Dirac Fock bound for a declared finite cubic program.

    The configuration algebra is generated by M15 matrix coefficients and
    their conjugates. Cutoff denotes total representation-word degree.
    Onsite terms are coordinate-independent. A current pointer generated
    from the same polynomial has the same step.
    """
    side = _int(side, 'periodic cubic side', 2)
    bounds = polynomial_bounds(degree, spacing)
    t = _frac(time, 'common duration')
    magnetic = _frac(plaquette_coefficient, 'plaquette coefficient')
    local = _frac(onsite_norm, 'total onsite operator norm')
    pointer = _frac(integrated_pointer_norm, 'integrated pointer norm')
    sites = side**3
    modes = 4*15*sites
    plaquettes = 3*sites
    b = modes*bounds['global_one_body_H_norm_upper']+plaquettes*magnetic+local
    result = cutoff_bound(cutoff, initial_degree,
                          max(bounds['maximum_coordinate_word_degree'], 4), b*t+pointer)
    return {'spatial_dimension': 3, 'sites': sites, 'Dirac_CAR_modes': modes,
            'full_Fock_dimension': str(1 << modes),
            'spatial_link_count': 3*sites, 'plaquette_count': plaquettes,
            'interaction_norm_upper': str(b),
            'one_body_bounds': {k: str(v) if isinstance(v, F) else v for k, v in bounds.items()},
            'complete_history_bound': result,
            'coordinate_completion_is_a_design_choice': True,
            'onsite_terms_are_coordinate_independent': True,
            'both_chiralities_retained': True,
            'physical_chiral_continuum_established': False}


@lru_cache(maxsize=1)
def rotor_coefficients():
    """Actual native singlet charge and source Dirac matrices, all16 CAR states."""
    source = source_data.load()
    if sum(source['charges']) != 0:
        raise ArithmeticError('the complete native matter determinant charge differs')
    singlets = [b for b in source['blocks'] if b['dimension'] == 1]
    if len(singlets) != 1:
        raise ArithmeticError('the native singlet multiplicity differs')
    charge = singlets[0]['q6Y']
    raw_gamma, _ = source['gamma']
    exact = tuple(s.ImmutableMatrix([[s.nsimplify(z, rational=True) for z in row]
                                    for row in g]) for g in raw_gamma)
    chart = s.kronecker_product(s.Matrix([[1, 1], [1, -1]])/s.sqrt(2), s.eye(2))
    beta = _clean(chart.H*exact[3]*chart)
    gamma = _clean(chart.H*exact[0]*chart)
    if beta != s.diag(1, 1, -1, -1) or not _clean(beta*gamma+gamma*beta).is_zero_matrix:
        raise ArithmeticError('the source spatial Dirac convention differs')
    one_hop = _clean((-beta+beta*gamma)/2)
    a = tuple(creator(4, j).H for j in range(4))
    def lift(matrix):
        return _clean(sum((value*a[i].H*a[j] for (i, j), value in matrix.todok().items()), s.zeros(16)))
    onsite, hop = lift(beta), lift(one_hop)
    number = _clean(sum((v.H*v for v in a), s.zeros(16)))
    if not _clean(number*hop-hop*number).is_zero_matrix:
        raise ArithmeticError('the specified rotor field changes total occupation')
    if not _clean(one_hop**2).is_zero_matrix or not _clean(one_hop.H*one_hop+one_hop*one_hop.H-s.eye(4)).is_zero_matrix:
        raise ArithmeticError('the complete one-cell overlap factor changed')
    return {'charge': charge, 'beta': beta, 'gamma': gamma, 'one_hop': one_hop,
            'onsite': onsite, 'hop': hop, 'number': number,
            'source_provenance': source['provenance']}


def exact_rotor(cutoff):
    k = _int(cutoff, 'electric flux cutoff')
    c = rotor_coefficients(); q = c['charge']
    if k < abs(q):
        raise ValueError('the prepared native flux lies outside the cutoff')
    flux = tuple(range(-k, k+1)); n = len(flux)
    shift = s.ImmutableSparseMatrix(n, n,
        {(j+q, j): 1 for j in range(n) if 0 <= j+q < n})
    electric = s.diag(*flux)
    identity = s.eye(n)
    h = _clean(s.kronecker_product(s.eye(16), electric**2/2)
               +s.kronecker_product(c['onsite'], identity)
               +s.kronecker_product(c['hop'], shift)
               +s.kronecker_product(c['hop'].H, shift.H))
    e = _clean(s.kronecker_product(s.eye(16), electric))
    current = _clean(s.I*q*(s.kronecker_product(c['hop'].H, shift.H)
                           -s.kronecker_product(c['hop'], shift)))
    gauss = _clean(s.kronecker_product(q*(c['number']-2*s.eye(16)), identity))
    if h != h.H or not _clean(s.I*(h*e-e*h)-current).is_zero_matrix:
        raise ArithmeticError('same-Hamiltonian electric-current identity differs')
    if not _clean(gauss*h-h*gauss).is_zero_matrix:
        raise ArithmeticError('the neutralized one-cell Gauss sector is not reducing')
    state = s.ImmutableSparseMatrix(16*n, 1,
        {(3*n+k, 0): 1/s.sqrt(2), (3*n+k+q, 0): s.I/s.sqrt(2)})
    if not _clean(gauss*state).is_zero_matrix or (state.H*state)[0] != 1:
        raise ArithmeticError('the complete prepared state is not normalized Gauss zero')
    expectation = lambda op: s.simplify((state.H*op*state)[0])
    if (expectation(e), expectation(current), expectation(h)) != (s.Rational(q, 2), q, 2+s.Rational(q*q, 4)):
        raise ArithmeticError('the explicit state/current/clock readout differs')
    return {'H': h, 'E': e, 'current': current, 'Gauss': gauss, 'state': state,
            'flux': flux, 'shift': shift, 'charge': q,
            'initial_flux_mean': expectation(e), 'initial_current': expectation(current),
            'initial_energy': expectation(h)}


def _gaussian(value):
    value = s.expand(value)
    real, imag = value.as_real_imag()
    if real.is_Integer is not True or imag.is_Integer is not True:
        raise ArithmeticError('the exact rotor matrix is not Gaussian integral after doubling')
    return int(real), int(imag)


def _gaussian_add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _gaussian_multiply(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _taylor_action(rows, vector, duration, order):
    """Gaussian-integer Horner evaluation of exp(-itH), H=rows/2.

    Return an integer vector and its single common rational denominator.
    The supplied vector can itself be an unnormalized record numerator.
    """
    duration = _frac(duration, 'record duration')
    order = _int(order, 'Taylor order')
    if len(rows) != len(vector):
        raise ValueError('the complete operator and input dimensions differ')
    step = 2*duration.denominator
    numerator, denominator = list(vector), 1
    for k in range(order, 0, -1):
        next_denominator = denominator*step*k
        following = []
        for row, entries in enumerate(rows):
            value = (vector[row][0]*next_denominator, vector[row][1]*next_denominator)
            moved = (0, 0)
            for col, entry in entries:
                moved = _gaussian_add(moved, _gaussian_multiply(entry, numerator[col]))
            # (-i * duration.numerator) times the doubled Hamiltonian.
            value = _gaussian_add(value, (duration.numerator*moved[1],
                                          -duration.numerator*moved[0]))
            following.append(value)
        numerator, denominator = following, next_denominator
    return numerator, denominator


def _prepared_rotor_rows(cutoff):
    """Exact 2H on the whole prepared N_F=2, flux mod q=0 sector."""
    c = rotor_coefficients(); q = c['charge']
    occupations = tuple(mask for mask in range(16) if mask.bit_count() == 2)
    flux = tuple(value for value in range(-cutoff, cutoff+1) if value % abs(q) == 0)
    index = {(mask, value): i for i, (mask, value)
             in enumerate((mask, value) for mask in occupations for value in flux)}
    n = len(index); entries = [dict() for _ in range(n)]
    def add(row, col, value):
        entries[row][col] = entries[row].get(col, s.Integer(0))+value
    for (mask, value), col in index.items():
        add(col, col, value*value)
        for other in occupations:
            add(index[other, value], col, 2*c['onsite'][other, mask])
            if (other, value+q) in index:
                add(index[other, value+q], col, 2*c['hop'][other, mask])
            if (other, value-q) in index:
                add(index[other, value-q], col, 2*c['hop'].H[other, mask])
    rows = tuple(tuple((col, _gaussian(value)) for col, value in sorted(row.items()) if value != 0)
                 for row in entries)
    initial = [(0, 0)]*n
    initial[index[3, 0]] = (1, 0)
    initial[index[3, q]] = (0, 1)
    # Its norm is sqrt(2); the exact probability denominator retains the 2.
    return rows, tuple(initial), occupations, flux


def recorded_rotor(cutoff, step_times=(F(1, 64), F(1, 64)), taylor_order=64):
    """Every two-time signed-flux record in the invariant N_F=2 sector.

    Exact coefficient identities define the parent on all16 Fock states.
    This execution uses the complete six-state Gauss-zero occupation block.
    Evolution uses Gaussian-rational Taylor polynomials with analytic
    Hermitian remainder bounds. Every record remains unnormalized; displayed
    decimals are derived from the exact nonnegative rational record weights.
    """
    k = _int(cutoff, 'electric flux cutoff')
    c = rotor_coefficients(); q = c['charge']
    if k < abs(q) or len(step_times) != 2:
        raise ValueError('complete prepared support and exactly two readout times required')
    times = tuple(_frac(t, 'record duration') for t in step_times)
    order = _int(taylor_order, 'Taylor order')
    flux = tuple(range(-k, k+1)); n = len(flux)
    rows, initial, occupations, active_flux = _prepared_rotor_rows(k)
    prepared, first_denominator = _taylor_action(rows, initial, times[0], order)
    records = [[F(0) for _ in flux] for _ in flux]
    small_n = len(active_flux)
    for first_small, first_value in enumerate(active_flux):
        conditional = [v if i % small_n == first_small else (0, 0)
                       for i, v in enumerate(prepared)]
        final, second_denominator = _taylor_action(rows, conditional, times[1], order)
        denominator = 2*(first_denominator*second_denominator)**2
        for last_small, last_value in enumerate(active_flux):
            numerator = sum(final[i][0]**2+final[i][1]**2
                            for i in range(last_small, len(final), small_n))
            records[first_value+k][last_value+k] = F(numerator, denominator)
    norm = sum((sum(row, F(0)) for row in records), F(0))
    hnorm = F(k*k, 2)+10
    errors = tuple((t*hnorm)**(order+1)/math.factorial(order+1) for t in times)
    delta = errors[0]+errors[1]+errors[0]*errors[1]
    half_trace = delta+delta*delta/2
    norm_error = abs(norm-1)
    if norm_error > 2*half_trace:
        raise ArithmeticError('unnormalized complete records exceed their exact Taylor bound')
    first_mean = sum((records[i][j]*flux[i] for i in range(n) for j in range(n)), F(0))
    second_mean = sum((records[i][j]*flux[j] for i in range(n) for j in range(n)), F(0))
    # For four CAR modes, ||dGamma(A)||<=4||A||<=4. Only hopping raises
    # flux, hence V norm<=8 and every insertion changes flux by at most q.
    bound = cutoff_bound(k, abs(q), abs(q), 8*sum(times))
    return {'cutoff': k, 'charge': q, 'Fock_parent_dimension': 16*n,
            'complete_Gauss_zero_dimension': 6*n,
            'exact_prepared_flux_residue_modulo_q': 0,
            'executed_prepared_residue_dimension': len(initial),
            'all_flux_values': list(flux), 'step_times': list(map(str, times)),
            'all_joint_flux_probabilities': [[float(p) for p in row] for row in records],
            'all_joint_flux_probabilities_exact': [[str(p) for p in row] for row in records],
            'record_entry_definition': 'squared unnormalized Taylor amplitudes; unitary probabilities enclosed by the stated bound',
            'record_count_including_zero_records': n*n,
            'normalization_error': float(norm_error), 'normalization_error_exact': str(norm_error),
            'first_flux_mean': float(first_mean), 'second_flux_mean': float(second_mean),
            'first_flux_mean_exact': str(first_mean), 'second_flux_mean_exact': str(second_mean),
            'finite_time_evolution': {'method': 'Gaussian-integer Horner Taylor action',
                'Taylor_order': order, 'Hamiltonian_norm_upper': str(hnorm),
                'single_step_operator_errors': list(map(str, errors)),
                'complete_record_isometry_error_upper': str(delta),
                'complete_record_half_l1_error_upper': str(half_trace),
                'any_single_event_error_upper': str(2*half_trace),
                'flux_mean_error_upper': str(2*k*half_trace),
                'record_probabilities_renormalized': False,
                'floating_evolution_used': False},
            'analytic_cutoff_bound': bound,
            'evolution_error_rationally_bounded': True,
            'spatial_continuum_or_chiral_species_identified': False}


def cutoff_comparison():
    smaller, larger = recorded_rotor(18), recorded_rotor(30)
    a = [[F(v) for v in row] for row in smaller['all_joint_flux_probabilities_exact']]
    b = [[F(v) for v in row] for row in larger['all_joint_flux_probabilities_exact']]
    tv = sum((abs((a[i-12][j-12] if 12 <= i < 49 and 12 <= j < 49 else F(0))-b[i][j])
              for i in range(61) for j in range(61)), F(0))/2
    bound = sum((F(r['analytic_cutoff_bound']['complete_record_TV_error_upper'])
                 for r in (smaller, larger)), F(0))
    time_bound = sum((F(r['finite_time_evolution']['complete_record_half_l1_error_upper'])
                      for r in (smaller, larger)), F(0))
    if tv > bound+time_bound:
        raise ArithmeticError('complete cutoff-comparison records exceed the analytic allowance')
    compact = []
    for r in (smaller, larger):
        compact.append({k: v for k, v in r.items()
                        if k not in ('all_joint_flux_probabilities', 'all_joint_flux_probabilities_exact')})
    return {'runs': compact, 'complete_joint_record_half_l1_difference': float(tv),
            'complete_joint_record_half_l1_difference_exact': str(tv),
            'sum_of_infinite_cutoff_TV_bounds': str(bound),
            'sum_of_time_evolution_half_l1_bounds': str(time_bound),
            'exact_finite_record_TV_difference_upper': str(min(F(1), tv+time_bound)),
            'comparison_is_full_record_law': True}


def certificate():
    exact = exact_rotor(6)
    return {'schema': 'oph.observer_outputs.field_route.v1', 'source_inputs': list(SOURCE_INPUTS),
            'chosen_architecture': 'gauge-coordinate completion with full Dirac Fock occupations',
            'architecture_uniqueness_required': False, 'measured_mass_inputs': [],
            'exact_native_charge': exact['charge'],
            'exact_initial_flux': str(exact['initial_flux_mean']),
            'exact_initial_current': str(exact['initial_current']),
            'exact_initial_energy': str(exact['initial_energy']),
            'initial_state_and_H_use_one_declared_clock': True,
            'native_C_E_pulse_tape_executed_for_this_field_program': False,
            'complete_rotor_cutoff_comparison': cutoff_comparison(),
            'three_dimensional_program_example': cubic_program_bound(2, 8, 200, 1, '1/1024'),
            'physical_chiral_particles_selected': [], 'physical_mass_predictions': [],
            'remaining_output_requirements': [
                'gauge-state control near Wilson-kernel gap closures for the exact-overlap limit',
                'a chiral fermion measure or canonical regulator with controlled complementary modes',
                'uniform interacting volume/cutoff/refinement limits and a suitable positive state',
                'the same interacting detector spectral response, units and experimental particle dictionary']}


if __name__ == '__main__':
    print(json.dumps(certificate(), ensure_ascii=False, indent=2, sort_keys=True))
