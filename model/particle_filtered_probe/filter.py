"""Complete finite instruments for a time-filtered primary boundary probe.

The primary operation is specified by exact finite algebraic apparatus
services. A small separate matrix control is executed; no spatial source
matrix or native tape is allocated.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import isqrt
import hashlib
import json
from pathlib import Path

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_physical_boundary_state import boundary
from model.particle_spatial_ground_state import ground
from model.particle_joint_reconstruction import positive_energy
from model.particle_primary_field_algebra import algebra as field_algebra
from model.particle_positive_window import window
from model.source_family_compatibility.spectral_window_access import cayley_service_bound
from model.spin_operational_position.bounds import generic_instrument_plan
from model.spin_packet_readout import controls as inherited_lcu
from model.spin_packet_readout.bounds import sqrt_upper

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_filtered_probe/__main__.py',
    'model/particle_filtered_probe/__init__.py',
    'model/particle_filtered_probe/filter.py',
    'model/particle_filtered_probe/PROOF.md',
    'model/particle_filtered_probe/test_filter.py',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_positive_window/PROOF.md',
    'model/particle_ground_spectral_records/PROOF.md',
    'model/spin_operational_position/PROOF.md',
    'model/spin_packet_readout/README.md',
    'model/source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md',
)


def exact(value, name, positive=False):
    return ground._q(value, name, positive)


def integer(value, name, minimum=1):
    return ground._i(value, name, minimum)


def clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(lambda x: s.cancel(s.expand_complex(x)))


def exact_matrix(matrix):
    result = s.ImmutableMatrix(matrix)
    if any(x.has(s.Float) or x.free_symbols or x.is_algebraic is not True for x in result):
        raise ValueError('finite exact algebraic matrix entries required')
    return result


def quadrature(half_width='1/16', cells=4):
    """Exact positive cell integrals, with midpoint nodes and total weight one."""
    t = exact(half_width, 'positive filter half-width', True)
    m = integer(cells, 'even quadrature cell count', 2)
    if m % 2 or m > 4096:
        raise ValueError('an even executed quadrature count at most 4096 is required')
    step = 2*t/m
    rows = []
    for j in range(m):
        lo = -t+j*step; hi = lo+step; node = (lo+hi)/2
        # The kink at zero is a cell boundary. Integrate the linear density.
        weight = step*(1-abs(node)/t)/t
        rows.append({'time': node, 'weight': weight, 'interval': [lo, hi]})
    if any(row['weight'] <= 0 for row in rows) or sum(row['weight'] for row in rows) != 1:
        raise ArithmeticError('the complete triangular quadrature is not normalized')
    return {'half_width': t, 'cells': m, 'step': step, 'rows': rows,
            'maximum_node_error': step/2, 'total_weight': F(1)}


def _ceil_log2(value):
    value = F(value); n = value.numerator.bit_length()-value.denominator.bit_length()
    power = F(2)**n
    return n+(power < value)


def factored_compiler_plan(system_bits, outcomes, error='1/1000000'):
    """Symbolic specialization of the inherited universal instrument bound."""
    q = integer(system_bits, 'padded source bits')
    k = integer(outcomes, 'fine output count')
    eps = exact(error, 'compiler half-diamond allowance', True)
    if eps >= 1:
        raise ValueError('compiler allowance below one required')
    r = (k-1).bit_length(); total = q+r
    return {'padded_source_bits': q, 'retained_record_bits': r,
            'padded_unitary_bits': total,
            'pulse_bound_exact_expression': f'20*2^{2*total}+4*2^{total}+10*{q}',
            'pulse_bound_factored_upper': {'coefficient': 34, 'power_of_two': 2*total},
            'port_fraction_bits': 2*total+_ceil_log2(1972/eps),
            'rank_one_tree_angle_bits': total+_ceil_log2(8*(total+2)/eps),
            'complete_halfdiamond_error_upper': eps,
            'coherent_source_code_load_export_service_required': True,
            'gatewise_Gauss_preservation_of_compiled_tape': False,
            'all_compiler_failure_and_code_complement_records_retained': True,
            'native_spatial_tape_executed': False}


def error_budget(spectral_lower, spectral_upper, half_width='1/16', cells=4,
                 ramsey_time='1', steps=16, compiler_error='1/1000000',
                 preparation_error='0'):
    lo = ground._signed(spectral_lower, 'source spectral lower')
    hi = ground._signed(spectral_upper, 'source spectral upper')
    if lo > hi:
        raise ValueError('ordered source spectral bounds required')
    t = exact(half_width, 'filter half-width', True)
    tau = exact(ramsey_time, 'Ramsey source time')
    count = integer(cells, 'even quadrature count', 2)
    if count % 2:
        raise ValueError('an even quadrature count required')
    r = integer(steps, 'Cayley subdivisions')
    sigma = exact(compiler_error, 'compiler allowance')
    prep = exact(preparation_error, 'preparation allowance')
    if max(sigma, prep) > 1:
        raise ValueError('half-trace allowances cannot exceed one')
    width = hi-lo; norm = max(abs(lo), abs(hi))
    # No scalar is subtracted from the executed source Hamiltonian.
    qerr = min(F(2), width*t/count)
    def clock_error(duration):
        if duration == 0 or norm == 0:
            return F(0)
        return cayley_service_bound(duration, duration, norm, r)['unitary_operator_error']
    filt_clock = clock_error(t); read_clock = clock_error(tau)
    amp = min(F(2), qerr+2*filt_clock)
    complete = min(F(1), 2*filt_clock+2*read_clock+sigma+prep)
    f0 = min(F(1), 2*amp+sigma+prep)
    score = min(F(2), 2*amp+2*read_clock+2*(sigma+prep))
    return {'source_spectral_interval': [lo, hi], 'source_norm_upper': norm,
            'source_commutator_norm_upper_for_unit_probe': width,
            'quadrature_operator_error_upper': qerr,
            'filter_clock_operator_error_per_evolution': filt_clock,
            'Ramsey_clock_operator_error_per_evolution': read_clock,
            'filtered_probe_operator_error_upper': amp,
            'complete_output_halftrace_error_to_exact_node_instrument': complete,
            'original_F0_error_to_continuous_triangle': f0,
            'original_ReF_or_ImF_error_to_continuous_triangle': score,
            'original_contrast_error_to_continuous_triangle': f0+score,
            'compiler_halfdiamond_allowance': sigma, 'preparation_halftrace_allowance': prep,
            'full_failure_maps_compared_to_continuous_filter': False,
            'quadrature_error_is_operator_and_original_response_error': True,
            'source_energy_zero_changed': False}


def precision_parameters(member, half_width='1/16', ramsey_time='1', response_error='1/1024'):
    """Choose finite counts without allocating the spatial source or selector."""
    contract.verify_member(member)
    eps = exact(response_error, 'response error target', True)
    if eps >= 1:
        raise ValueError('response target below one required')
    t = exact(half_width, 'filter half-width', True)
    tau = exact(ramsey_time, 'Ramsey time')
    source = member['preparation_plan']['source']
    lo, hi = F(source['spectral_lower']), F(source['spectral_upper'])
    norm = max(abs(lo), abs(hi)); width = hi-lo
    needed = 8*width*t/eps
    cells = 2*max(1, -(-needed.numerator//needed.denominator))
    square = 64*max(t, tau)**3*norm**3/(12*eps)
    ceiling = -(-square.numerator//square.denominator)
    steps = isqrt(ceiling)
    if steps*steps < ceiling:
        steps += 1
    return {'half_width': str(t), 'cells': cells, 'ramsey_time': str(tau),
            'steps': max(1, steps), 'compiler_error': str(eps/32),
            'preparation_error': '0'}


def _relative_sqrt_upper(value):
    """Reuse the inherited dyadic bound without an absolute precision floor."""
    value = F(value)
    if value == 0:
        return F(0)
    scale = F(2)**(_ceil_log2(value)//2)
    return scale*sqrt_upper(value/(scale*scale))


def lueders_error_budget(effect_error, compiler_error='0', history_length=1,
                         preparation_error='0', interleaved_clock_error='0'):
    """Complete histories, including every binary record and untouched reference.

    The clock allowance is the SUM of unitary operator errors over all
    interleaved source evolutions. Exact source evolution has allowance zero.
    """
    delta = min(F(1), exact(effect_error, 'effect operator error'))
    sigma = exact(compiler_error, 'per-instrument compiler allowance')
    zeta = exact(preparation_error, 'initial preparation allowance')
    clocks = exact(interleaved_clock_error, 'total interleaved clock allowance')
    count = integer(history_length, 'finite history length')
    root = _relative_sqrt_upper(2*delta)
    return {'effect_operator_error_upper': delta,
            'two_Kraus_Stinespring_error_formula': 'sqrt(2*delta)',
            'two_Kraus_Stinespring_error_rational_upper': root,
            'one_complete_instrument_halftrace_error_upper': min(F(1),root+sigma),
            'history_length': count,
            'history_complete_halftrace_error_upper': min(F(1),count*(root+sigma)+zeta+clocks),
            'all_binary_histories_and_quantum_references_retained': True,
            'interleaved_source_clock_error_sum': clocks,
            'per_instrument_compiler_halfdiamond_allowance': sigma,
            'preparation_halftrace_allowance': zeta,
            'history_probability_conditioned_on_selected_outcome': False}


def selected_plan(member, half_width='1/16', cells=4, ramsey_time='1', steps=16,
                  compiler_error='1/1000000', preparation_error='0', probe_descriptor=None):
    contract.verify_member(member)
    if probe_descriptor is not None:
        field_algebra.verify_probe(member, probe_descriptor)
        if F(probe_descriptor['half_width']) != F(half_width):
            raise ValueError('apparatus filter width differs from the authenticated probe')
    base = member['preparation_plan']['source']
    budgets = error_budget(base['spectral_lower'], base['spectral_upper'],
                           half_width, cells, ramsey_time, steps,
                           compiler_error, preparation_error)
    degree = member['source']['per_link_word_degree']
    charges = boundary.charge_spectrum(degree)
    spatial = boundary.common_payload(member)
    p = spatial['physical_density_floor_factored']['negative_exponent']
    even = any(row['cut_charges'][0] % 2 == 0 for row in spatial['strictly_positive_selected_state_witnesses'])
    odd = any(row['cut_charges'][0] % 2 != 0 for row in spatial['strictly_positive_selected_state_witnesses'])
    if not even or not odd:
        raise ArithmeticError('the actual boundary witnesses do not support both parity signs')
    count = integer(cells, 'quadrature count', 2)
    t = exact(half_width, 'filter half-width', True)
    def compiler_plan(outcomes):
        if F(compiler_error) == 0:
            return {'route':'exact declared algebraic apparatus services',
                    'compiled_tape_required':False,'retained_fine_outcomes':outcomes,
                    'full_spatial_matrix_or_tape_executed':False}
        return factored_compiler_plan(base['dimension_bits_upper'],outcomes,compiler_error)
    result = {
        'member_id': member['member_id'], **member['ids'],
        'probe': {'name': 'first cut-endpoint hypercharge parity',
                  'formula': 'A=sum_(q in S_ell) (-1)^q E_q(Q6 at e1)',
                  'owned_region': 'origin cell and its three outgoing links',
                  'single_endpoint_spectrum': charges, 'operator_norm': 1,
                  'Hermitian_unitary': True, 'source_Gauss_preserving': True,
                  'all_Fock_and_link_word_states_retained': True},
        'filter': {'half_width': t, 'cells': count, 'step': 2*t/count,
                   'node_formula': '-T+(j+1/2)*2T/M',
                   'weight_formula': '(2/M)*(1-abs(node)/T)',
                   'weights_sum': 1, 'success_operator': 'B_M=sum_j w_j alpha_tj(A)',
                   'all_ancilla_outcomes': count, 'Ramsey_outcomes_per_phase_setting': 2*count,
                   'Ramsey_source_time': exact(ramsey_time, 'Ramsey time'),
                   'Cayley_subdivisions_per_evolution': steps,
                   'success_probability_renormalized': False},
        'apparatus': {'source_primitive': 'controlled exact algebraic Cayley transform of the original finite H',
                      'neutral_ancilla_primitive': 'exact algebraic Householder PREP, its adjoint, and Ramsey phase/Hadamard',
                      'boundary_primitive': 'controlled coherent parity of the complete central-charge PVM',
                      'reference_H': member['reference']['idle_H'],
                      'reference_actively_counter_evolved': False,
                      'negative_times_are_declared_reversible_controls': True,
                      'all_final_quantum_references_retained': True,
                      'local_point_time_detector_claimed': False,
                      'autonomous_uncontrolled_repair_implements_apparatus': False,
                      'full_spatial_matrix_or_tape_executed': False},
        'errors': budgets,
        'compiler': compiler_plan(2*count),
        'finite_survival': {'unfiltered_connected_variance_lower': {'coefficient': 4, 'base': 2, 'negative_exponent': 2*p},
                            'continuous_filter_distance_from_A_upper': min(F(2), F(budgets['source_commutator_norm_upper_for_unit_probe'])*t/3),
                            'uniform_positive_survival_established': False,
                            'nonzero_energy_window_weight_computed': False,
                            'formula': 'sqrt(Var_rho(A_f)) >= max(0,sqrt(Var_rho(A))-||A_f-A||)'},
        'primary_F0_or_Ramsey_contrast_computed': False,
        'calibration_or_learned_mass_inputs': False,
        'joint_refinement_or_physical_particle_identification': False,
    }
    if probe_descriptor is not None:
        result['probe'] = {'authenticated_menu': probe_descriptor,
                           'Hermitian_contraction': True, 'operator_norm_upper': 1,
                           'neutral_dilation': '[[A,sqrt(I-A^2)],[sqrt(I-A^2),-A]]',
                           'source_clock_on_dilation': 'I_probe tensor H',
                           'blank_probe_ancilla_to_blank_block': 'A'}
        result['filter']['probe_ancilla_outcomes'] = 2
        result['filter']['all_ancilla_outcomes'] = 2*count
        result['filter']['Ramsey_outcomes_per_phase_setting'] = 4*count
        result['filter']['success_operator'] = 'B_M=<0_probe,0_selector|W|0_probe,0_selector>'
        result['apparatus']['boundary_primitive'] = 'controlled exact algebraic Hermitian-contraction dilation of the authenticated probe'
        result['compiler'] = compiler_plan(4*count)
        result['finite_survival']['unfiltered_connected_variance_lower'] = None
        result['finite_survival']['parity_witness_floor_transferred_to_this_probe'] = False
        charge = probe_descriptor['kind'] == 'charge'
        effect_error = F(budgets['filtered_probe_operator_error_upper'])/(1 if charge else 2)
        result['lueders_history_apparatus'] = {
            'ideal_effect': 'E=A_f' if charge else 'E=(I+A_f)/2',
            'finite_effect': 'E_M=B_tilde' if charge else 'E_M=(I+B_tilde)/2',
            'complete_Kraus_menu': ['sqrt(E_M)','sqrt(I-E_M)'],
            'exact_primitive': 'algebraic positive square roots and their complete neutral two-outcome dilation',
            'positive_contraction_from_positive_quadrature': True,
            'success_effect_equals_LCU_blank_effect': False,
            'LCU_blank_effect': 'B_tilde^dagger B_tilde',
            'reference_H': member['reference']['idle_H'],
            'errors': lueders_error_budget(effect_error,compiler_error,
                                          preparation_error=preparation_error),
            'compiler': compiler_plan(2),
            'full_spatial_matrix_or_tape_executed': False}
    return contract.serial(result)


def bind_member(member, **kwargs):
    result = selected_plan(member, **kwargs)
    return {'schema': 'oph.filtered_primary_probe.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, receipt, **kwargs):
    if receipt != bind_member(member, **kwargs):
        raise ValueError('filtered probe receipt differs from the actual source and apparatus')


def cayley(matrix, duration, steps=1):
    h = exact_matrix(matrix); value = ground._signed(duration, 'signed source duration')
    t = s.Rational(value.numerator, value.denominator)
    r = integer(steps, 'Cayley subdivisions')
    if h.rows != h.cols or h != h.H:
        raise ValueError('Hermitian original Hamiltonian required')
    eye = s.eye(h.rows); z = s.I*t/(2*r)
    return clean(((eye-z*h)*(eye+z*h).inv())**r)


def lcu_instrument(hamiltonian, probe, half_width='1/16', cells=4, steps=16):
    """Execute the entire finite SELECT and retain every analyzer outcome."""
    h = exact_matrix(hamiltonian); a = exact_matrix(probe)
    if a.shape != h.shape or clean(a.H*a) != s.eye(a.rows):
        raise ValueError('the source probe must be an exact unitary')
    q = quadrature(half_width, cells)
    vector = s.Matrix([s.sqrt(s.Rational(row['weight'].numerator, row['weight'].denominator)) for row in q['rows']])
    blank = s.eye(cells)[:, 0]; difference = blank-vector
    prep = clean(s.eye(cells)-difference*difference.H/(1-vector[0]))
    if clean(prep.H*prep) != s.eye(cells) or clean(prep*blank-vector) != s.zeros(cells, 1):
        raise ArithmeticError('the exact neutral PREP is invalid')
    branches = []
    for row in q['rows']:
        u = cayley(h, row['time'], steps)
        branches.append(clean(u.H*a*u))
    unitary = clean(s.kronecker_product(prep.H, s.eye(a.rows))*s.diag(*branches)
                    *s.kronecker_product(prep, s.eye(a.rows)))
    kraus = tuple(unitary[r*a.rows:(r+1)*a.rows, :a.rows] for r in range(cells))
    complete = clean(sum((k.H*k for k in kraus), s.zeros(a.rows)))
    success = clean(sum((s.Rational(row['weight'].numerator, row['weight'].denominator)*u
                         for row, u in zip(q['rows'], branches, strict=True)), s.zeros(a.rows)))
    if complete != s.eye(a.rows) or clean(kraus[0]-success) != s.zeros(a.rows):
        raise ArithmeticError('LCU analyzer lost original outcomes or changed its success operator')
    return {'unitary': unitary, 'kraus': kraus, 'success': success,
            'prep': prep, 'quadrature': q, 'node_unitaries': tuple(branches)}


def positive_sqrt(matrix):
    """Finite algebraic functional calculus on a positive Hermitian matrix."""
    value = exact_matrix(matrix)
    if value.rows != value.cols or value != value.H:
        raise ValueError('an exact Hermitian matrix required')
    if any(eigenvalue.is_nonnegative is not True for eigenvalue in value.eigenvals()):
        raise ValueError('a positive semidefinite matrix required')
    root = clean(value**s.Rational(1,2))
    if root != root.H or clean(root*root-value) != s.zeros(value.rows):
        raise ArithmeticError('the exact positive square root failed')
    return root


def hermitian_dilation(probe):
    """Exact neutral-ancilla dilation of an algebraic Hermitian contraction."""
    a = exact_matrix(probe)
    if a.rows != a.cols or a != a.H:
        raise ValueError('an exact Hermitian contraction required')
    rest = clean(s.eye(a.rows)-a*a)
    root = positive_sqrt(rest)
    dilation = clean(a.row_join(root).col_join(root.row_join(-a)))
    if clean(root*root-rest) != s.zeros(a.rows) or clean(dilation.H*dilation) != s.eye(2*a.rows):
        raise ArithmeticError('the algebraic contraction dilation failed')
    return dilation


def contraction_instrument(hamiltonian, probe, half_width='1/16', cells=4, steps=16):
    """Retain selector AND neutral probe-ancilla outputs for the same source."""
    h = exact_matrix(hamiltonian); a = exact_matrix(probe)
    if a.shape != h.shape:
        raise ValueError('source and contraction dimensions differ')
    dilation = hermitian_dilation(a)
    extended = s.kronecker_product(s.eye(2),h)
    lifted = lcu_instrument(extended,dilation,half_width,cells,steps)
    dimension = h.rows
    kraus = tuple(clean(k[p*dimension:(p+1)*dimension,:dimension])
                  for k in lifted['kraus'] for p in (0,1))
    if clean(sum((k.H*k for k in kraus),s.zeros(dimension))) != s.eye(dimension):
        raise ArithmeticError('a neutral dilation outcome was omitted')
    expected = s.zeros(dimension)
    for row in lifted['quadrature']['rows']:
        u = cayley(h,row['time'],steps)
        expected += s.Rational(row['weight'].numerator,row['weight'].denominator)*u.H*a*u
    if clean(kraus[0]-expected) != s.zeros(dimension):
        raise ArithmeticError('the blank selector/dilation block is not the filtered probe')
    return {'dilation':dilation,'lifted_unitary':lifted['unitary'],
            'kraus':kraus,'success':kraus[0],
            'fine_labels':[(r,p) for r in range(cells) for p in (0,1)],
            'source_dimension':dimension,'reference_clock_changed':False}


def lueders_instrument(effect):
    """The complete two-outcome instrument, not an LCU success instrument."""
    e = exact_matrix(effect)
    root = positive_sqrt(e); complement = positive_sqrt(s.eye(e.rows)-e)
    kraus = (root,complement)
    dilation = clean(root.row_join(-complement).col_join(complement.row_join(root)))
    if clean(sum((k.H*k for k in kraus),s.zeros(e.rows))) != s.eye(e.rows):
        raise ArithmeticError('Lueders instrument is not complete')
    if clean(dilation.H*dilation) != s.eye(2*e.rows):
        raise ArithmeticError('Lueders neutral dilation is not unitary')
    return {'effect':e,'kraus':kraus,'unitary':dilation,
            'fine_labels':['effect','complement'],
            'success_probability_renormalized':False}


def filtered_lueders_instrument(hamiltonian, probe, kind='charge',
                               half_width='1/16', cells=4, steps=16):
    """Execute positive quadrature before exact two-outcome functional calculus."""
    h = exact_matrix(hamiltonian); a = exact_matrix(probe)
    if a.shape != h.shape or a != a.H or kind not in ('charge','plaquette'):
        raise ValueError('a matching Hermitian charge or plaquette menu required')
    if kind == 'charge' and clean(a*a-a) != s.zeros(a.rows):
        raise ValueError('the unsmeared charge probe must be a projection')
    positive_sqrt(s.eye(a.rows)-a*a)
    q = quadrature(half_width,cells); b = s.zeros(a.rows)
    for row in q['rows']:
        u = cayley(h,row['time'],steps)
        b += s.Rational(row['weight'].numerator,row['weight'].denominator)*u.H*a*u
    b = clean(b)
    e = b if kind == 'charge' else clean((s.eye(a.rows)+b)/2)
    result = lueders_instrument(e)
    result.update({'filtered_probe':b,'kind':kind,'quadrature':q,
                   'LCU_blank_effect':clean(b.H*b),
                   'source_H':h,'reference_clock_changed':False})
    return result


def compose_instruments(first, second):
    """Sequential finite instruments; index zero retains both original blanks."""
    first=tuple(exact_matrix(k) for k in first);second=tuple(exact_matrix(k) for k in second)
    if not first or not second:
        raise ValueError('two complete nonempty finite instruments required')
    dimension=first[0].rows
    for menu in (first,second):
        if any(k.shape != (dimension,dimension) for k in menu):
            raise ValueError('instrument source dimensions differ')
        if clean(sum((k.H*k for k in menu),s.zeros(dimension))) != s.eye(dimension):
            raise ValueError('each input instrument must retain every outcome')
    operators=tuple(clean(right*left) for left in first for right in second)
    dimension=operators[0].rows
    if clean(sum((k.H*k for k in operators),s.zeros(dimension))) != s.eye(dimension):
        raise ArithmeticError('composition omitted a fine outcome')
    return operators


def finite_effect_control():
    """Separate exact two-level contraction and complete-history execution."""
    h=s.diag(0,4);p=s.ones(2)/2
    lcu=contraction_instrument(h,p,'1/2',2,1)
    charge=filtered_lueders_instrument(h,p,'charge','1/2',2,1)
    loop=filtered_lueders_instrument(h,2*p-s.eye(2),'plaquette','1/2',2,1)
    history=compose_instruments(charge['kraus'],loop['kraus'])
    rho=s.diag(s.Rational(7,10),s.Rational(3,10))
    if charge['effect'] != lcu['success'] or loop['effect'] != charge['effect']:
        raise ArithmeticError('the complete contraction control changed its effect')
    return {'source_dimension':2,'reference_dimension':2,'is_selected_primary_source':False,
            'H':str(h),'rho':str(rho),'projection':str(p),
            'filtered_effect':str(charge['effect']),
            'LCU_blank_effect':str(charge['LCU_blank_effect']),
            'LCU_original_blank_probability':str(s.trace(rho*charge['LCU_blank_effect'])),
            'Lueders_original_effect_probability':str(s.trace(rho*charge['effect'])),
            'LCU_selector_and_probe_output_count':len(lcu['kraus']),
            'Lueders_complete_outcome_count':len(charge['kraus']),
            'two_measurement_history_probabilities':[str(s.simplify(s.trace(rho*k.H*k))) for k in history],
            'both_complete_instruments_verified':True,
            'success_effects_identified_with_each_other':False,
            'reference_clock_changed':False,'calibration_or_learned_masses':False}


def ramsey_instrument(hamiltonian, kraus, duration='1', steps=16, imaginary=False):
    if type(imaginary) is not bool:
        raise ValueError('a real or imaginary Ramsey setting required')
    u = cayley(hamiltonian, duration, steps); phase = -s.I if imaginary else s.S.One
    rows = []
    for label, k in enumerate(kraus):
        moved = clean(u.H*k*u)
        for bit in (0, 1):
            operator = clean((k+(-1)**bit*phase*moved)/2)
            rows.append({'filter_label': label, 'Ramsey_bit': bit, 'K': operator,
                         'original_score': (-1)**bit if label == 0 else 0})
    complete = clean(sum((row['K'].H*row['K'] for row in rows), s.zeros(u.rows)))
    if complete != s.eye(u.rows):
        raise ArithmeticError('complete Ramsey instrument lost a failure outcome')
    return rows


def finite_control():
    """Exact small matrix execution on a separate stationary three-level source."""
    data = positive_energy.rotated_three_level_control('1/1024')
    h, rho = data['H'], data['rho']
    rotation = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5), 0],
                         [s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    a = clean(rotation*s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])*rotation.T)
    lcu = lcu_instrument(h, a)
    b = lcu['success']; u = cayley(h, '1', 16)
    f0 = s.simplify(s.trace(rho*b.H*b))
    corr = s.expand_complex(s.trace(rho*b.H*u.H*b*u)).simplify()
    records = []
    for imaginary in (False, True):
        rows = ramsey_instrument(h, lcu['kraus'], imaginary=imaginary)
        for row in rows:
            row['p'] = s.simplify(s.trace(rho*row['K'].H*row['K']))
        if s.simplify(sum(row['p'] for row in rows)-1) != 0:
            raise ArithmeticError('the complete small control is not normalized')
        score = s.simplify(sum(row['original_score']*row['p'] for row in rows))
        if s.simplify(score-(s.im(corr) if imaginary else s.re(corr))) != 0:
            raise ArithmeticError('original Ramsey score does not equal the unnormalized Gram correlation')
        records.append({'imaginary_setting': imaginary,
                        'records': [{'filter_label': r['filter_label'], 'Ramsey_bit': r['Ramsey_bit'],
                                     'original_probability': str(r['p']), 'original_score': r['original_score']}
                                    for r in rows], 'score': str(score), 'probability_sum': '1'})
    errors = error_budget('0', '3', compiler_error='0')
    contrast = F(f0-s.re(corr)); margin = max(F(0), contrast-errors['original_contrast_error_to_continuous_triangle'])
    band = window.positive_window(contrast, errors['original_contrast_error_to_continuous_triangle'],
                                  1, 1, 32, F(1,256), F(1,4), 1024)
    if band['original_positive_band_weight_lower'] <= 0:
        raise ArithmeticError('the separate finite response lost its positive-band certificate')
    return {'H': h, 'rho': rho, 'A': a, 'LCU': lcu,
            'receipt': contract.serial({'source_dimension': 3, 'reference_dimension': 3,
                'is_selected_primary_source': False, 'all_eight_outcomes_per_phase_retained': True,
                'exact_original_F0': str(f0), 'exact_original_ReF': str(s.re(corr)),
                'exact_original_ImF': str(s.im(corr)), 'exact_original_contrast': str(contrast),
                'continuous_triangle_contrast_lower': margin,
                'continuous_triangle_contrast_lower_dyadic20': ground.outward_dyadic(margin,20)[0],
                'mean_ground_excess': F(1, 1024), 'half_width': F(1, 16), 'Ramsey_time': 1,
                'operator_norm_squared': 1, 'first_absolute_filtered_moment_upper': 32,
                'errors': errors, 'phase_settings': records,
                'positive_window_control': band,
                'conditional_success_normalization_applied': False,
                'calibration_or_learned_masses': False})}


def certificate():
    paths = provenance.definition_inputs(SOURCE_INPUTS)
    hashes = provenance.definition_hashes(paths)
    menu_member = contract.member(1)
    menu_parameters = precision_parameters(menu_member, half_width='1/100000')
    menus = [field_algebra.probe(menu_member,kind=kind) for kind in ('charge','plaquette')]
    result = {'schema': 'oph.filtered_probe_certificate.v1',
              'authenticated_field_menu_plans': [bind_member(menu_member,probe_descriptor=probe,**menu_parameters) for probe in menus],
              'selected_member_plan': bind_member(contract.member(0), **precision_parameters(contract.member(0))),
              'exact_small_control': finite_control()['receipt'],
              'exact_contraction_and_Lueders_control':finite_effect_control(),
              'inherited_complete_LCU_control': inherited_lcu.exact_control(),
              'inherited_compiler_control': generic_instrument_plan(2, 8, F(1, 10**6)),
              'source_sha256': hashes,
              'primary_nonzero_energy_survival_or_continuum_proved': False}
    if provenance.definition_hashes(paths) != hashes:
        raise ValueError('filter source inputs changed during calculation')
    return result
