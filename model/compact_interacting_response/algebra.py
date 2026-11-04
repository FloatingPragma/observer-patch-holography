"""Exact compact-band identities and rational spectral/readout bounds."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.common_source.chiral_candidate import action
from model.compact_defect_pair.certification import complete_certificate
from model.compact_defect_pair.topology import spectral_bounds
from model.interacting_defect_response.algebra import annihilators
from model.mcpr_native_memory_bridge import certificate as mass_certificate
from model.mcpr_native_memory_bridge.construction import source
from model.native_flavor_orbit.intervals import I, rounded_out, sqrt_rational


def exact_block(mass, coupling, e1, e2):
    c = annihilators(4)
    ns = [a.H*a for a in c]
    n = sum(ns, s.zeros(16))
    pair = sum((c[i].H*c[j].H+c[j]*c[i] for i, j in ((0, 1), (2, 3))), s.zeros(16))
    kinetic = e1*(ns[0]-ns[1])+e2*(ns[2]-ns[3])
    return kinetic+mass*pair+coupling*(n-2*s.eye(16))**2, n, pair


def exact_identities():
    m, u, e1, e2 = s.symbols('m u e1 e2', real=True)
    z = s.Symbol('z')
    h, n, pair = exact_block(m, u, e1, e2)
    expected = z*(z-4*u)*((z-2*u)**2-4*(u*u+m*m))
    expected *= s.prod(z-a*e1-b*e2 for a in (-1, 1) for b in (-1, 1))
    expected *= s.prod((z-u-a*e)**2-m*m for e in (e1, e2) for a in (-1, 1))
    if s.expand(h.charpoly(z).as_expr()-expected) != 0:
        raise ArithmeticError('complete symbolic compact interaction spectrum failed')
    q4 = s.diag(*(sum((-3, 3, -3, 3)[j] for j in range(4) if i & (1 << j)) for i in range(16)))
    if h*q4 != q4*h or (n-2*s.eye(16))**2*pair == pair*(n-2*s.eye(16))**2:
        raise ArithmeticError('charge preservation or nonquadratic dynamics failed')
    labels = action.ordering()['all_left_components']
    q = {'Q': 1, 'u^c': -1, 'd^c': -1, 'L': -3, 'e^c': 3}
    q15 = s.diag(*(q[row['internal_type']] for row in labels))
    currents = action.generators(families=False)
    if any(q15*g != g*q15 for g in currents):
        raise ArithmeticError('native-current charge commutator failed')
    for channel in action._representation()[4].values():
        for k in channel['components']:
            if q15*k+k*q15.T != s.zeros(15):
                raise ArithmeticError('charged pairing charge identity failed')
    charges = [int(q15[i, i]) for _ in range(12) for i in range(15)]
    return {'compact_spatial_modes': 4, 'native_internal_modes': 45,
            'localized_annihilators': 180, 'localized_Nambu_rows': 360,
            'charged_pairing_blocks': 84, 'neutral_unpaired_modes': 12,
            'four_mode_slice_spectator_charged_blocks': 82,
            'four_mode_slice_occupied_neutral_spectators': 6,
            'four_mode_slice_spectator_number': 88,
            'native_current_count': len(currents),
            'interaction_commutes_with_every_uniform_native_current': True,
            'charge_3B_minus_3L_Fock_range': [sum(x for x in charges if x < 0), sum(x for x in charges if x > 0)],
            'charged_pairing_conserves_3B_minus_3L': True,
            'neutral_Majorana_pairing_included': False,
            'full_Fock_dimension_of_slice': 16,
            'characteristic_polynomial': 'z*(z-4u)*((z-2u)^2-4*(u^2+m^2))*prod_{a,b=+/-1}(z-a*e1-b*e2)*prod_{e=e1,e2;a=+/-1}((z-u-a*e)^2-m^2)',
            'ground_energy': 'min(2u-2sqrt(u^2+m^2), -e1-e2, u-m-max(e1,e2))',
            'blocked_ground_poles_for_first_positive_energy_mode': ['e1+u-m', 'e1+u+m'],
            'blocked_ground_addition_plus_removal_weights': ['1/2', '1/2'],
            'full_core_vacuum_is_four_mode_slice_ground': False}


def cube_root_box(value):
    value = I(value) if not isinstance(value, I) else value
    if not 0 <= value.lo <= value.hi <= 1:
        raise ValueError('cube-root certificate domain is [0,1]')
    brackets = []
    for target in (value.lo, value.hi):
        lo, hi = F(0), F(1)
        for _ in range(90):
            middle = (lo+hi)/2
            if middle**3 <= target:
                lo = middle
            else:
                hi = middle
        brackets.append((lo, hi))
    return I(brackets[0][0], brackets[1][1])


def mass_boxes(branch='forward'):
    squares = mass_certificate.source_shape_boxes(branch)
    alpha = I(source(branch)[1])
    kappa = (alpha/50)/(1+alpha/31+alpha**2/310)
    common_cube = mass_certificate.exponential(3*kappa)/(2*6**14*squares[0]*squares[1]*squares[2])
    return [cube_root_box(common_cube*x**3) for x in squares]


@lru_cache(maxsize=1)
def ground_certificate():
    certificate = complete_certificate()
    positive = [I(*row) for row in certificate['ordered_energy_intervals'] if F(0) < F(row[0]) < F(row[1]) < 2]
    if len(positive) != 2:
        raise ArithmeticError('complete compact certificate has wrong positive low-band count')
    delta = F(spectral_bounds(64, 1, 3)['ordered_band_energy_infinite_volume_error_upper'])
    infinite = [value+I(-delta, delta) for value in positive]
    coupling = F(1, 128)
    masses = mass_boxes()
    inverse = mass_certificate.admissibility('forward')['inverse_branches'][0]
    p_box = I(*inverse['p'])
    gain_box = sum(masses, I(0))/p_box
    if not gain_box.hi < 1 or not p_box.hi < F(1, 12):
        raise ArithmeticError('actual memory source does not satisfy the compact norm bound class')
    rows = []
    for mass in masses:
        r = I(sqrt_rational(coupling*coupling+mass.lo*mass.lo).lo,
              sqrt_rational(coupling*coupling+mass.hi*mass.hi).hi)
        # Compare the actual blocked energy with the lowest energy in every
        # other pair-parity sector, including both spatial signs.
        paired_margin = 2*coupling-2*r+infinite[0]+infinite[1]
        odd_margin = min(x.lo for x in infinite)+coupling-mass.hi
        if paired_margin.lo <= 0 or odd_margin <= 0:
            raise ArithmeticError('compact source slice ground is not certified blocked')
        poles = [infinite[0]+coupling-mass, infinite[0]+coupling+mass]
        if any(x.lo <= 0 for x in poles):
            raise ArithmeticError('single-fermion readout pole positivity failed')
        rows.append({'source_mass_interval': rounded_out(mass, 20).strings(),
                     'paired_sector_energy_margin_lower': str(paired_margin.lo),
                     'single_odd_sector_energy_margin_lower': str(odd_margin),
                     'single_fermion_pole_intervals': [rounded_out(x, 15).strings() for x in poles],
                     'residues': ['1/2', '1/2'], 'unique_slice_ground_regime': 'blocked'})
    return {'finite_source_eigenvalues_exhausted': certificate['dimension'],
            'finite_positive_spatial_energy_intervals': [x.strings() for x in positive],
            'infinite_positive_spatial_energy_intervals': [x.strings() for x in infinite],
            'exterior_error_per_ordered_spatial_energy': str(delta),
            'arbitrary_four_mode_ordered_energy_error_upper': str(2*delta),
            'arbitrary_four_mode_gap_error_upper': str(4*delta),
            'blocked_pole_error_upper': str(delta),
            'actual_memory_inverse_p_interval': p_box.strings(),
            'actual_memory_gain_interval': rounded_out(gain_box, 15).strings(),
            'actual_memory_gain_below_one_and_contraction_below_one_over12': True,
            'declared_coupling': str(coupling), 'source_mass_branch': 'forward',
            'charged_lepton_slices': rows,
            'method': 'complete dyadic source spectrum, rational Taylor mass bounds, cubic and square-root enclosures'}


def compact_memory_bound(gain, coupling, p, time, steps, *, mass=64, half_separation=1):
    if any(isinstance(x, (float, bool)) for x in (gain, coupling, p, time)):
        raise ValueError('exact rational input bounds required')
    g, u, p, t = map(F, (gain, coupling, p, time))
    if g < 0 or u < 0 or not 0 <= p < 1 or t < 0 or type(steps) is not int or steps < 1:
        raise ValueError('invalid compact memory bound domain')
    band = F(spectral_bounds(mass, half_separation, half_separation+2)['band_abs_energy_upper'])
    norm = 90*band+90*(1+g)+8100*u
    transient = 4*norm*t/(steps*(1-p))
    stationary = (8/(1-p)+4)*norm*norm*t*t/steps
    return {'gain_upper': str(g), 'coupling': str(u), 'contraction_upper': str(p),
            'time': str(t), 'steps': steps, 'band_abs_energy_upper': str(band),
            'core_Fock_norm_upper': str(norm),
            'unread_arbitrary_initial_CQ_error_upper': str(min(F(2), transient+stationary)),
            'unread_stationary_product_error_upper': str(min(F(2), stationary)),
            'complete_history_arbitrary_CQ_error_squared_upper': str(min(F(4), 2*(transient+stationary))),
            'complete_history_stationary_product_error_squared_upper': str(min(F(4), 2*stationary)),
            'localized_annihilators': 180, 'uniform_in_exterior_volume': True,
            'uniform_in_lattice_spacing': False, 'rare_conditioned_history_guarantee': False}
