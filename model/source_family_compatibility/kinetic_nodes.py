"""Exact kinetic winding and its conditional free-field interpretation.

The source hopping is fixed. A many-body particle pole or a physical Standard
Model assignment is not inferred from a one-body winding number.
"""
from fractions import Fraction as F
from itertools import product

import sympy as s

from model.native_mass_channels import algebra as masses
from model.weyl_geometric_modular import algebra as geometric
from model.ballistic_free_boundary import packet as free_parent


VALLEYS = tuple(product((0, 1), repeat=3))


def _positive_rational(value, name):
    if type(value) is int:
        value = F(value)
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or F(value.numerator, value.denominator) != value or value <= 0):
        raise ValueError(name+' must be an exact positive rational')
    return F(value.numerator, value.denominator)


def sphere_gap(radius, speed=1):
    """Certified kinetic boundary gap for a node sphere of radius <= 1."""
    radius = _positive_rational(radius, 'radius')
    speed = _positive_rational(speed, 'speed')
    if radius > 1:
        raise ValueError('radius must be at most one radian')
    return speed*radius*(1-radius*radius/6)


def winding_stability(radius, perturbation_norm, speed=1):
    """Strict sufficient condition, for continuous traceless 2x2 perturbations."""
    gap = sphere_gap(radius, speed)
    if type(perturbation_norm) is int:
        perturbation_norm = F(perturbation_norm)
    if (type(perturbation_norm) is not F or type(perturbation_norm.numerator) is not int
            or type(perturbation_norm.denominator) is not int
            or perturbation_norm.denominator <= 0 or perturbation_norm < 0
            or F(perturbation_norm.numerator, perturbation_norm.denominator) != perturbation_norm):
        raise ValueError('perturbation norm must be an exact nonnegative rational')
    return {'boundary_gap_lower_bound': gap,
            'degree_preserved_by_straight_homotopy': perturbation_norm < gap,
            'continuous_traceless_symbol_required': True,
            'interacting_particle_pole_implied': False}


def projector_orientation():
    """Derive the Chern sign from the two actual Pauli bands."""
    theta, phi = s.symbols('theta phi', real=True)
    n = (s.sin(theta)*s.cos(phi), s.sin(theta)*s.sin(phi), s.cos(theta))
    h = sum((x*y for x, y in zip(n, geometric.pauli())), s.zeros(2))
    rows = []
    for band in (-1, 1):
        P = (s.eye(2)+band*h)/2
        commutator = P.diff(theta)*P.diff(phi)-P.diff(phi)*P.diff(theta)
        curvature = s.trigsimp(s.expand(s.trace(P*commutator)))
        expected = band*s.I*s.sin(theta)/2
        if s.trigsimp(curvature-expected) != 0:
            raise ArithmeticError('Pauli projector orientation changed')
        rows.append({'band_sign': band, 'unit_degree_Chern_number': band})
    return {'definition': '(1/(2*pi*i))*integral Tr(P dP wedge dP)',
            'Berry_connection': '-i u^dagger du',
            'outward_oriented_momentum_sphere': True, 'bands': rows}


def node_rows():
    rows = []
    for b in VALLEYS:
        signs, chi, lift = geometric.frame(b)
        jacobian = tuple(-x for x in signs)
        if jacobian[0]*jacobian[1]*jacobian[2] != chi:
            raise ArithmeticError('kinetic winding and Weyl chirality disagree')
        for j, sigma in enumerate(geometric.pauli()):
            if lift.H*(jacobian[j]*sigma)*lift != chi*sigma:
                raise ArithmeticError('kinetic node changed its continuum frame')
        rows.append({'source_node_pi_units': list(b),
                     'source_Jacobian_over_speed': list(jacobian),
                     'kinetic_winding': chi,
                     'positive_band_Chern_number': chi,
                     'negative_band_Chern_number': -chi,
                     'free_field_chirality': chi,
                     'positive_charge_free_helicity': str(F(chi, 2)),
                     'conjugate_free_helicity': str(F(-chi, 2)),
                     'free_positive_energy': 'v*|p|',
                     'free_rest_mass': '0',
                     'physical_particle_assignment': None})
    return rows


def mass_mixing():
    """Exact sixteen-row source channels connect opposite winding sectors."""
    gamma, beta = masses.folded_generators()
    chirality = s.diag(*(geometric.frame(b)[1]*s.eye(2) for b in VALLEYS))
    zero = s.zeros(16)
    for matrix in gamma:
        if chirality*matrix-matrix*chirality != zero:
            raise ArithmeticError('kinetic evolution changes valley chirality')
    rows = []
    for channel, matrix in enumerate(beta):
        if chirality*matrix+matrix*chirality != zero:
            raise ArithmeticError('mass channel does not pair opposite winding')
        for kinetic in gamma:
            if matrix*kinetic+kinetic*matrix != zero:
                raise ArithmeticError('source mass channel fails to open a kinetic gap')
        flip = (1, 1, 1) if channel == 0 else tuple(int(j == channel-1) for j in range(3))
        rows.append({'source_mass_channel': channel, 'valley_flip': list(flip),
                     'opposite_chirality_only': True,
                     'broken_one_step_translations': [j for j, value in enumerate(flip) if value],
                     'two_step_translations_preserved': True})
    return {'physical_spin_valley_rows': 16, 'channels': rows,
            'mass_scales_for_supplied_coefficients': '|m0 +/- sqrt(m1^2+m2^2+m3^2)|',
            'couplings_or_order_parameters_selected': False}


def certificate():
    rows = node_rows()
    if sum(row['kinetic_winding'] for row in rows) != 0:
        raise ArithmeticError('total lattice kinetic winding is nonzero')
    return {'actual_Weyl_symbol': '-c*sum_j sigma_j*sin(k_j)',
            'nodes': rows, 'total_kinetic_winding': 0,
            'each_winding_multiplicity': 4,
            'projector_orientation': projector_orientation(),
            'radius_one_gap_over_speed': sphere_gap(1),
            'mass_mixing': mass_mixing(),
            'replayed_free_parent_theorem': free_parent.build()['theorem'],
            'free_field_scope': {
                'family': 'centered full-ground c=v/a, 0<=U=o(c^(2/7)), fixed v>0',
                'parent': 'model/ballistic_free_boundary/README.md',
                'entire_vacuum_hierarchy_reused': True,
                'Lorentz_spin_representation_parent': 'model/weyl_geometric_modular/README.md',
                'distinct_Weyl_fields': 8, 'Dirac_pairings': 4,
                'negative_band_counted_as_extra_independent_species': False,
                'default_U_over_c_eight_in_this_class': False},
            'scope': {'kinetic_momentum_space_topology': True,
                      'all_native_repair_protected_charge': False,
                      'real_space_defect_classification': False,
                      'default_interacting_Green_function_winding_derived': False,
                      'physical_Standard_Model_assignment': False}}
