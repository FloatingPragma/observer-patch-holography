"""Rank-two characteristic polynomials for the actual native charge voltage.

The magnetic compression belongs to the prepared cyclic walk span. The full
event carrier also has fixed phases. Its eigenphases are not identified as
particle rest poles, and the controlled clock schedule has no proved common
stationary physical generator.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
import sympy as s

from model.native_charged_propagation import source as motion
from model.characteristic_constraints.bounds import trig_interval
from model.calibrated_measurement.bounds import root_interval
from model.openwave_native_process.source import graph as native_graph

ROOT = Path(__file__).resolve().parents[2]
X = s.Symbol('x')


def integer_laplacian(character):
    if character not in (0,3):
        raise ValueError('integer zero or parity parent required')
    g = motion.carrier()
    matrix = np.diag(g.degree).astype(np.int64)
    for i,j,voltage in g.edges:
        sign = (-1)**(voltage if character == 3 else 0)
        matrix[i,j] = matrix[j,i] = -sign
    return matrix


def cofactor_polynomial(matrix, coefficients, i, j):
    """One adjugate column from Cayley-Hamilton; exact integer vectors."""
    n = len(matrix)
    vectors = [int(k == j) for k in range(n)]
    entries = []
    rows = [[(k,int(v)) for k,v in enumerate(row) if v] for row in matrix]
    for order in range(n):
        entries.append(vectors[i])
        if order < n-1:
            vectors = [sum(v*vectors[k] for k,v in row) for row in rows]
            vectors[j] += int(coefficients[order+1])
    return s.Poly.from_list(entries, X, domain=s.ZZ)


@lru_cache(None)
def polynomials():
    g = motion.carrier()
    i,j = g.index[(0,1)],g.index[(1,5)]
    parents = {}; result = {}; details = []
    for k in (0,3):
        matrix = integer_laplacian(k)
        poly = s.Poly(s.Matrix(matrix).charpoly(X).as_expr(), X, domain=s.ZZ)
        minor = np.delete(np.delete(matrix,[i,j],axis=0),[i,j],axis=1)
        minor_poly = s.Poly(s.Matrix(minor).charpoly(X).as_expr(), X, domain=s.ZZ)
        cofactor = cofactor_polynomial(matrix,poly.all_coeffs(),i,j)
        # Delta L has diagonal zero and one conjugate off-diagonal pair.
        # 2Re(alpha) = +3 for the zero parent, -3 for the parity parent.
        real_twice = 3 if k == 0 else -3
        changed = poly-real_twice*cofactor-3*minor_poly
        parents[k] = matrix
        result[k] = poly
        if k == 0:
            result[2] = result[4] = changed
        else:
            result[1] = result[5] = changed
        details.append({'parent_character':k,'changed_pair_indices':[i,j],
                        'two_real_alpha_exact':str(real_twice),'alpha_squared_norm_exact':'3',
                        'parent_coefficients':list(map(int,poly.all_coeffs())),
                        'offdiagonal_cofactor_coefficients':list(map(int,cofactor.all_coeffs())),
                        'deleted_pair_minor_coefficients':list(map(int,minor_poly.all_coeffs()))})
    if result[0].eval(0) != 0 or any(result[k].eval(0) == 0 for k in range(1,6)):
        raise ArithmeticError('zero and nonzero character kernels disagree')
    return result,details


def finite_field_control():
    """Independent full-matrix Faddeev coefficients in three prime fields."""
    polys,_ = polynomials(); g = motion.carrier(); checks = []
    for prime in (103,109,127):
        zeta = next(z for z in range(2,prime) if (z*z-z+1)%prime == 0)
        for k in range(6):
            matrix = np.diag(g.degree).astype(np.int64)
            for i,j,voltage in g.edges:
                matrix[j,i] = -pow(zeta,(k*voltage)%6,prime)%prime
                matrix[i,j] = -pow(zeta,(-k*voltage)%6,prime)%prime
            previous = np.eye(g.states,dtype=np.int64); coefficients = [1]
            for order in range(1,g.states+1):
                value = (matrix@previous)%prime
                coefficient = -int(np.trace(value))*pow(order,-1,prime)%prime
                coefficients.append(coefficient)
                previous = value.copy()
                previous[np.diag_indices(g.states)] = (previous.diagonal()+coefficient)%prime
            if coefficients != [int(c)%prime for c in polys[k].all_coeffs()]:
                raise ArithmeticError('independent full magnetic characteristic polynomial differs')
            checks.append({'prime':prime,'sixth_root':zeta,'character':k,'coefficients_checked':67})
    return checks


def smallest_positive_interval(character):
    poly = polynomials()[0][character]
    g = motion.carrier()
    T = motion.apply_T(g,np.eye(g.states),character)
    diagnostic = float(np.linalg.eigvalsh(g.events*(np.eye(g.states)-T))[0])
    if character == 0:
        if poly.eval(0) != 0:
            raise ArithmeticError('constant zero mode absent')
        return F(0),F(0),1,diagnostic
    upper = s.Rational(int(np.ceil((diagnostic+.01)*1000)),1000)
    intervals = poly.intervals(eps=s.Rational(1,10**13),inf=0,sup=upper)
    if len(intervals)!=1 or intervals[0][0][0] <= 0:
        raise ArithmeticError('minimum positive eigenvalue isolation failed')
    (lo,hi),multiplicity = intervals[0]
    if poly.count_roots(0,lo)!=0:
        raise ArithmeticError('an omitted lower positive root exists')
    return F(lo),F(hi),int(multiplicity),diagnostic


def arccos_interval(target):
    """Directed rational cosine bounds on the principal phase in [0,1]."""
    if not F(0) < target <= F(1):
        raise ValueError('positive prepared transfer eigenvalue required')
    if target == 1:return F(0),F(0)
    a,b = F(0),F(1)
    if trig_interval(b,True)[1] >= target:
        raise ValueError('phase outside the stated principal interval')
    for _ in range(52):
        mid = (a+b)/2
        lo,hi = trig_interval(mid,True)
        if lo > target:a = mid
        elif hi < target:b = mid
        else:break
    return a,b


def enlargement_controls(scale):
    """Exact Rayleigh and local-observable bounds on three native carriers.

    For characters 2 and 4 the parity voltage is invisible. The remaining
    voltage changes exactly one configuration edge. A normalized constant
    source vector therefore has magnetic Dirichlet energy 3/N. This proves
    a gap upper bound without enumerating the enlarged two-load Hilbert space.
    """
    local=motion.carrier()
    for k in (2,4):
        cost=sum(motion.norm2(motion.sub((1,0),motion.phase(k,c))) for _,_,c in local.edges)
        if cost!=3:raise ArithmeticError('the source has more than one charged flux edge')
    rows=[]
    for level in (0,1,2):
        graph=native_graph(level);P=graph.ports;N=P*(P-1)//2;D=2*len(graph.seams)
        if graph.seams.count((0,5))!=1 or len(set(graph.seams))!=len(graph.seams):
            raise ArithmeticError('the literal distinguished native seam changed')
        rayleigh=F(3,N)
        # arcsin(x) lies between x and x/sqrt(1-x²), x²=lambda/(2D).
        frequency_square_upper=2*rayleigh/(1-rayleigh/(2*D))
        gev_upper=scale*root_interval(frequency_square_upper)[1]
        regional_probability=F(66,N)
        rows.append({'native_level':level,'ports':P,'two_load_configurations':N,'uniform_events':D,
            'nontrivial_charge_characters':[2,4],
            'constant_source_Rayleigh_quotient_exact':str(rayleigh),
            'prepared_cyclic_frequency_squared_upper_exact':str(frequency_square_upper),
            'same_anchor_phase_energy_GeV_upper_exact':str(gev_upper),
            'uniform_reference_both_loads_in_first_cell_probability_exact':str(regional_probability),
            'bounded_regional_indicator_variance_exact':str(regional_probability*(1-regional_probability)),
            'all_two_load_coordinates_enumerated':False,
            'proof_shortcut':'All even-character edge phases are trivial except the one native(0,5),spectator1 configuration edge; its exact squared phase difference is3.'})
    if not all(F(a['prepared_cyclic_frequency_squared_upper_exact'])>F(b['prepared_cyclic_frequency_squared_upper_exact']) for a,b in zip(rows,rows[1:])):
        raise ArithmeticError('enlargement gap bound does not decrease')
    return {'clock':'one coherent event has duration1/sqrt(2M), with the same source-phase energy unit on every carrier',
        'length_convention':'one native graph edge has length1; this is carrier enlargement, not an identified physical mesh limit',
        'state':'constant source trial vector for the gap bound; uniform two-load population reference for the regional detector',
        'regional_observable':'indicator that both loads occupy the first twelve-port cell',
        'rows':rows,'held_out_level':2,
        'proved_consequence':'The even-character prepared-sector frequency gap has no positive uniform lower bound in this declared enlargement and clock convention.',
        'reference_consequence':'The bounded first-cell indicator has vanishing variance in the fixed-two-load uniform reference. A nontrivial field limit needs a different established state, density or observable scaling.',
        'other_constructions_or_reference_states_excluded':False,
        'four_dimensional_continuum_or_Yang_Mills_limit_proved':False}


def run():
    polys, details = polynomials()
    modular = finite_field_control()
    old = json.loads((ROOT/'model/openwave_native_calibrated_energy/results.json').read_text())['energy_attachment']
    anchor = F(old['electron_anchor_GeV_exact'])
    scale = F(old['one_common_GeV_per_source_unit_exact'])
    if scale*F(1,3072) != anchor:
        raise ArithmeticError('inherited single electron calibration changed')
    rows=[];g=motion.carrier();clock_factor=root_interval(F(g.events))
    for k in range(6):
        lo,hi,multiplicity,diagnostic=smallest_positive_interval(k)
        phase=(arccos_interval(1-lo/g.events)[0],arccos_interval(1-hi/g.events)[1])
        frequency=(phase[0]*clock_factor[0],phase[1]*clock_factor[1])
        rows.append({'character':k,'compression_L_min_interval_exact':list(map(str,(lo,hi))),
                     'multiplicity':multiplicity,'minimum_diagnostic':diagnostic,
                     'prepared_cyclic_principal_phase_interval_exact':list(map(str,phase)),
                     'prepared_cyclic_frequency_interval_exact':list(map(str,frequency)),
                     'conditional_phase_energy_GeV_interval_exact':list(map(str,(scale*frequency[0],scale*frequency[1]))),
                     'characteristic_polynomial_coefficients':list(map(int,polys[k].all_coeffs()))})
    return {'schema':'oph.native_magnetic_phase_spectrum.v1',
        'verdict':'PASS_EXACT_SIX_CHARACTER_PREPARED_PHASE_SPECTRA',
        'same_source':'complete H66 two-load configuration graph and declared active voltage',
        'rank_two_identity':'det(xI-L-deltaL)=det(xI-L)-2Re(alpha)cofactor_ij(xI-L)-abs(alpha)^2 det((xI-L)_delete_ij)',
        'rank_two_parents':details,'independent_finite_field_checks':modular,'rows':rows,
        'enlargement_controls':enlargement_controls(scale),
        'unit_contract':{'electron_anchor_GeV_exact':str(anchor),'GeV_per_declared_source_unit_exact':str(scale),
                         'extra_time_map':'One coherent W event is assigned 1/sqrt(2M) source phase-time units, on each carrier with M seams.',
                         'time_map_is_architecture_choice':True,
                         'source_clock_and_walk_have_proved_common_stationary_generator':False,
                         'physical_rest_mass_assignment':False,
                         'additional_mass_targets_used':False,
                         'electron_anchor_assignment_is_empirical':True},
        'scope':{'full_event_carrier_positive_gap':False,
                 'spectrum_of_flux_dressed_joint_electric_program':False,
                 'fixed_event_modes_retained':True,
                 'prepared_sector_phase_is_physical_mass':False,
                 'spatial_translation_or_physical_momentum_proved':False,
                 'normalizable_continuum_particle_pole_proved':False,
                 'exact_charged_lepton_hierarchy_proved':False,
                 'new_independent_GeV_prediction':False,
                 'literal_OpenWave_mass_row_complete':False,
                 'new_Nature_validation':False}}
