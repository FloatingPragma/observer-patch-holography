"""Exact finite transfer controls for a continuum spectral-gap criterion.

The continuum criterion needs a reconstructed positive-energy source and a
dense local test family. Finite matrices do not establish those hypotheses.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_yang_mills_gap_contract/transfer.py',
    'model/particle_yang_mills_gap_contract/PROOF.md',
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/provenance.py',
    'model/particle_common_contract/PROOF.md',
    'extra/yang_mills_gap_clay_problem.tex',
    'Lean/ObserverPatchHolography/RepairGapChain.lean',
    'Lean/Screen/KogutSusskindFiberRateComparison.lean',
    'model/physical_space_time/occupation_field_results.json',
    'model/observer_time_fields/REFINEMENT.md',
)


def rational(value, name, *, positive=False):
    if type(value) not in (int, str, F):
        raise ValueError(name+' requires exact rational input')
    result = F(value)
    if result < 0 or positive and result == 0:
        raise ValueError(name+' has the wrong sign')
    return result


def clean(matrix):
    return s.ImmutableMatrix(matrix).applyfunc(s.simplify)


def exact_matrix(matrix):
    result = clean(matrix)
    if any(x.has(s.Float) or x.free_symbols or x.is_algebraic is not True
           for x in result):
        raise ValueError('exact algebraic matrix entries are required')
    return result


def finite_transfer(transfer, vacuum, observables, ratio='1/2', time_step='1'):
    """Check original centered weights of a positive finite transfer matrix.

    The bound is C_A(k h) <= Var(A) ratio**k. A complete centered span
    establishes a finite global gap; a smaller span certifies its probes only.
    """
    t = exact_matrix(transfer); omega = exact_matrix(vacuum)
    q = rational(ratio, 'decay ratio', positive=True)
    h = rational(time_step, 'physical time step', positive=True)
    if q >= 1 or t.rows != t.cols or omega.shape != (t.rows, 1):
        raise ValueError('a square transfer, vacuum column and ratio below one are required')
    if t != t.H or clean(omega.H*omega) != s.ones(1):
        raise ValueError('a Hermitian transfer and normalized vacuum are required')
    if clean(t*omega-omega) != s.zeros(t.rows, 1):
        raise ValueError('the supplied vacuum is not fixed by this transfer')
    eigenvalues = t.eigenvals()
    if eigenvalues.get(s.S.One) != 1:
        raise ValueError('a unique fixed vacuum is required')
    if any(x.is_real is not True or x.is_positive is not True
           or (1-x).is_nonnegative is not True for x in eigenvalues):
        raise ValueError('the transfer must have spectrum in (0,1]')
    projectors = {}
    for value in eigenvalues:
        projector = s.eye(t.rows)
        for other in eigenvalues:
            if other != value:
                projector = clean(projector*(t-other*s.eye(t.rows))/(value-other))
        if clean(projector*projector-projector) != s.zeros(t.rows) or projector != projector.H:
            raise ArithmeticError('exact transfer spectral projection failed')
        projectors[value] = projector
    rows = []; columns = []
    for observable in observables:
        a = exact_matrix(observable)
        if a.shape != t.shape:
            raise ValueError('each observable must act on the declared source')
        mean = s.simplify((omega.H*a*omega)[0])
        vector = clean(a*omega-mean*omega)
        variance = s.simplify((vector.H*vector)[0]); columns.append(vector)
        weights = []
        for value in sorted(eigenvalues, key=s.default_sort_key):
            weight = s.simplify((vector.H*projectors[value]*vector)[0])
            if weight.is_nonnegative is not True:
                raise ArithmeticError('a spectral record has a negative weight')
            if weight != 0 and (s.Rational(q.numerator, q.denominator)-value).is_nonnegative is not True:
                raise ValueError('a retained spectral weight violates the proposed decay rate')
            weights.append({'transfer_value': str(value), 'original_weight': str(weight)})
        if s.simplify(sum(s.sympify(row['original_weight']) for row in weights)-variance) != 0:
            raise ArithmeticError('the centered spectral weights are incomplete')
        rows.append({'mean': str(mean), 'variance': str(variance), 'weights': weights})
    if not columns:
        raise ValueError('at least one actual observable is required')
    rank = s.Matrix.hstack(*columns).rank()
    return {'dimension': t.rows, 'time_step': str(h), 'decay_ratio': str(q),
            'physical_gap_lower_expression': '-log('+str(q)+')/('+str(h)+')',
            'centered_span_rank': rank, 'full_excitation_span': rank == t.rows-1,
            'all_original_spectral_weights_retained': True, 'probes': rows,
            'continuum_or_Yang_Mills_identification': False}


def correlation(transfer, vacuum, observable, steps):
    """Compute a complete centered Euclidean matrix element, with no filtering."""
    if type(steps) is not int or steps < 0:
        raise ValueError('a nonnegative integer time count is required')
    t=exact_matrix(transfer); omega=exact_matrix(vacuum); a=exact_matrix(observable)
    mean=(omega.H*a*omega)[0]
    vector=clean(a*omega-mean*omega)
    return s.simplify((vector.H*(t**steps)*vector)[0])


def low_energy_weight_bound(prefactor, ratio, threshold, steps, error=0):
    """Weight at transfer values >= threshold from one diagonal correlator.

    The continuum theorem sends the regulator error to zero at fixed steps
    before sending steps to infinity. This function makes neither limit.
    """
    c=rational(prefactor, 'prefactor'); q=rational(ratio, 'ratio', positive=True)
    r=rational(threshold, 'threshold', positive=True); err=rational(error, 'error')
    if not q < r <= 1 or type(steps) is not int or steps < 0:
        raise ValueError('require 0 < ratio < threshold <= 1 and integer steps')
    return (c*q**steps+err)/r**steps


def lattice_gap_window(spacing, physical_floor, variance_floor, energy_upper):
    """Necessary lattice window under BOTH a uniform gap and a surviving witness."""
    a=rational(spacing, 'spacing', positive=True)
    delta=rational(physical_floor, 'physical gap floor', positive=True)
    v=rational(variance_floor, 'variance floor', positive=True)
    m=rational(energy_upper, 'energy upper', positive=True)
    if delta > m/v:
        raise ValueError('the proposed gap contradicts the finite-energy noncollapse witness')
    return {'dimensionless_gap_lower': str(a*delta),
            'dimensionless_gap_upper': str(a*m/v),
            'requires_uniform_finite_gap_and_witness': True,
            'not_necessary_for_invisible_vanishing_weight_prelimit_modes': True}


def verify_primary_scaling_source(source):
    """The inverse-spacing proof consumes these literal action premises."""
    expected={
        'electric_coefficient':'105', 'magnetic_coefficient':'1/7',
        'regulator':'finite_range_Wilson',
        'onsite_one_body':'3 beta_D/a',
        'interior_forward_hop':'(-beta_D+beta_D gamma_i) R(U_i)/(2a), plus Hermitian adjoint',
        'electric_action':'105 sum_link C_trace/(2a)',
        'magnetic_action':'-sum_interior_plaquettes ReTr_M15 R(U_p)/(105a)',
        'energy_zero':'unchanged negative plaquette trace; no added scalar',
        'compression_order':'form complete parent products, then compress',
        'word_projection':'tensor_product_per_link_Haar_orthogonal',
    }
    if any(source.get(key)!=value for key,value in expected.items()):
        raise ValueError('the declared primary action differs from the scaling theorem premises')


def primary_scaling():
    from model.particle_common_contract import contract
    members=[contract.member(level) for level in (0,1)]
    rows=[]
    for member in members:
        source=member['source']
        verify_primary_scaling_source(source)
        rows.append({'member_id':member['member_id'], **member['ids'],
                     'level':source['level'], 'spacing':source['spacing'],
                     'side':source['side'], 'link_degree':source['per_link_word_degree'],
                     'ground_excess_target':member['state']['target_ground_excess']})
    return {'members': rows,
        'Hamiltonian_identity':'H_n = a_n^(-1) Hhat(L_n,ell_n;105,1/7)',
        'all_terms_share_inverse_spacing':True,
        'dimensionless_operator_changes_with_volume_and_cutoff':True,
        'uniform_gap_or_continuum_vacuum_derived':False,
        'ground_excess_target_is_a_gap_lower_bound':False,
        'bare_couplings_run_toward_asymptotic_freedom':False,
        'pure_compact_simple_Yang_Mills_source':False}


def native_clock_diagnostics():
    """Re-express retained native eigenvalue diagnostics in two declared clocks.

    Exact arithmetic on stored decimal diagnostics does not certify their
    eigenvalues. This is a separate two-dimensional comparison source.
    """
    raw=json.loads((ROOT/'model/physical_space_time/occupation_field_results.json').read_text(),
                   parse_float=str)
    rows=[]
    for row in raw['actual_source_screen_Green']['actual_source_levels']:
        level=row['source_level'];scale=4**level;seams=630*scale
        value=F(row['low_global_modes'][0]['native_L_eigenvalue_diagnostic'])
        attempts=1260*16**level;native=value/(2*seams)
        rows.append({'level':level, 'ports':row['ports'], 'seams':seams,
            'laplacian_eigenvalue_decimal_diagnostic':str(value),
            'per_attempt_gap_diagnostic':str(native),
            'attempts_per_comparison_time':attempts,
            'comparison_time_gap_diagnostic':str(attempts*native),
            'exact_normalization_identity':attempts==2*seams*scale})
    return {'rows':rows, 'certified_eigenvalue_enclosures':False,
        'same_source_as_selected_Wilson_Gauss_program':False,
        'four_dimensional_Yang_Mills_mass_gap':False,
        'clock_was_fitted_to_measured_mass':False}


def certificate():
    omega=s.Matrix([3,4])/5; p=omega*omega.T; t=p+(s.eye(2)-p)/2
    a=s.diag(0,1)
    lap=s.Matrix([[1,-1,0],[-1,2,-1],[0,-1,1]])
    coupled=s.eye(3)-lap/4; uniform=s.ones(3,1)/s.sqrt(3)
    return {'schema':'oph.euclidean_gap_transfer.v1',
        'source_sha256':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                         for path in SOURCE_INPUTS},
        'finite_nonuniform_vacuum':finite_transfer(t,omega,[a]),
        'finite_coupled_source':finite_transfer(coupled,uniform,
            [s.diag(-1,0,1),s.diag(1,-2,1)],ratio='3/4'),
        'same_primary_scaling':primary_scaling(),
        'separate_native_clock_diagnostics':native_clock_diagnostics(),
        'continuum_decay_hypotheses_instantiated_on_primary':False,
        'Clay_Yang_Mills_existence_or_mass_gap_proved':False}
