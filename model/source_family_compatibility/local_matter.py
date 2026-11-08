"""Complete local-history repair bounds on the declared full CAR source."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.defect_field_bridge.full_domain import finite_branches, RAW_SLOTS
from model.source_family_compatibility.defect_budget import exact_source, rational
from model.field_symmetry import symmetry
from model.full_fock_field.field import FieldGeometry


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('integer outside the admitted domain')
    return value


def _ceil(value):
    return -((-value.numerator)//value.denominator)


def boundary_upper(radius, support_radius, support_sites, duration, coupling):
    """Rational upper bound using exp(x)<=2**ceil(2x), without a Fock array."""
    _integer(radius, 1); _integer(support_radius); _integer(support_sites, 1)
    if radius <= support_radius or support_sites > (2*support_radius+1)**3:
        raise ValueError('nested cubes and valid support size required')
    time, c = rational(duration), rational(coupling)
    if not time or not c:
        return F(0)
    exponent = _ceil(96*c*time)-(radius-support_radius)
    factor = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
    return min(F(2), 6*c*support_sites*(2*radius+1)**2*time*factor)


def local_history_bound(radius, support_radius, support_sites, duration,
                        coupling, interventions, attempts, edges):
    _integer(interventions); _integer(attempts); _integer(edges, 1)
    boundary = boundary_upper(radius, support_radius, support_sites, duration, coupling)
    suffix = min(F(2), (2*interventions+1)*boundary)
    near = min(F(1), F(6*(2*radius+1)**3, edges))
    return {'boundary_operator_error': boundary, 'suffix_operator_error': suffix,
            'near_seam_probability_upper': near,
            'accessible_output_halftrace_upper': min(F(1), attempts*(near+suffix)),
            'global_raw_archive_compared': False,
            'discarded_raw_feedback_allowed': False,
            'initial_reference_entanglement_allowed': True}


def plan_local_history(support_radius, support_sites, duration, coupling,
                       interventions, attempts, tolerance):
    """A terminating exact radius/volume search for a declared local experiment."""
    _integer(support_radius); _integer(support_sites, 1)
    _integer(interventions); _integer(attempts)
    time, c, eta = rational(duration), rational(coupling), rational(tolerance, positive=True)
    if eta > 1 or support_sites > (2*support_radius+1)**3:
        raise ValueError('valid support and tolerance at most one required')
    radius = support_radius+max(1, _ceil(96*c*time))
    while attempts*(2*interventions+1)*boundary_upper(
            radius, support_radius, support_sites, time, c) > eta/2:
        radius += 1
    edge_target = 12*attempts*(2*radius+1)**3/eta
    low = 2*radius+1; high = low
    while 3*(high-1)*high*high < edge_target:
        high *= 2
    while low < high:
        middle = (low+high)//2
        if 3*(middle-1)*middle*middle >= edge_target:
            high = middle
        else:
            low = middle+1
    side = low; edges = 3*(side-1)*side*side
    bound = local_history_bound(radius, support_radius, support_sites, time, c,
                                interventions, attempts, edges)
    if bound['accessible_output_halftrace_upper'] > eta:
        raise ArithmeticError('finite local-history planner missed its requested bound')
    return {'radius': radius, 'outer_cube_side': side, 'outer_cube_sites': side**3,
            'outer_cube_edges': edges, 'attempts': attempts,
            'interventions': interventions, 'tolerance': eta, **bound}


def interaction_witness(coupling, contact, duration):
    c, u, time = rational(coupling, positive=True), rational(contact, positive=True), rational(duration, positive=True)
    remainder = 9*c*c+12*c*u+u*u
    if time > u/(8*remainder):
        raise ValueError('duration exceeds the quantitative interaction interval')
    return {'centered_ground_occupation': F(1, 2),
            'variance_quadratic_coefficient': u*u/4,
            'four_point_variance_lower': u*u*time*time/16,
            'one_connected_four_point_absolute_lower': u*u*time*time/64,
            'double_commutator_norm_upper': remainder,
            'duration_maximum': u/(8*remainder)}


@lru_cache(maxsize=None)
def annihilator(sites, mode):
    _integer(sites, 1); _integer(mode)
    if sites > 3 or mode >= 2*sites:
        raise ValueError('finite oracle admits at most three two-mode sites')
    matrix = s.zeros(1 << (2*sites))
    for mask in range(matrix.cols):
        if mask & (1 << mode):
            matrix[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return s.ImmutableMatrix(matrix)


def repair_dual(observable, *, sites, seam, law, omit_signs=False):
    """All padded actual coin/raw slots, including their fair-coin weights."""
    rows = finite_branches(sites=sites, flavors=2, law=law, seam=seam)
    dimension = 1 << (2*sites)
    if observable.shape != (dimension, dimension) or len(rows) != 2*RAW_SLOTS:
        raise ValueError('complete full-domain observable required')
    result = s.zeros(dimension)
    for entries in rows.values():
        # Sparse partial injections avoid dense 64**3 products.
        for initial, final, sign in entries:
            for other_initial, other_final, other_sign in entries:
                value = observable[final, other_final]
                if value:
                    result[initial, other_initial] += value*(1 if omit_signs else sign*other_sign)/2
    return s.ImmutableMatrix(result)


def finite_locality_controls():
    identity = s.eye(64)
    a, b = annihilator(3, 2), annihilator(3, 3)
    basis = tuple(left*right for left in (identity, a, a.H, a.H*a)
                  for right in (identity, b, b.H, b.H*b))
    for law in ('C', 'E'):
        for observable in basis:
            if repair_dual(observable, sites=3, seam=(0, 2), law=law) != observable:
                raise ArithmeticError('disjoint complete repair changed a local CAR observable')
    wrong = repair_dual(a, sites=3, seam=(0, 2), law='E', omit_signs=True)
    if wrong == a:
        raise ArithmeticError('crossing-seam CAR-sign countercontrol lost sensitivity')

    # Both alternatives have total occupation three and the same parity.
    alternatives = (21, 14)
    density = s.zeros(64)
    for x in alternatives:
        for y in alternatives:
            density[x, y] = s.Rational(1, 2)
    a0, b0 = annihilator(3, 0), annihilator(3, 1)
    for left in (identity, a0, a0.H, a0.H*a0):
        for right in (identity, b0, b0.H, b0.H*b0):
            observable = left*right
            after = repair_dual(observable, sites=3, seam=(1, 2), law='E')
            if s.trace(density*after) != s.trace(density*observable):
                raise ArithmeticError('remote-archive control changed the local output')
    collision = F(0)
    rows = finite_branches(sites=3, flavors=2, law='E', seam=(1, 2))
    for (coin, raw), entries in rows.items():
        if raw == 8:
            collision += sum(F(int(density[x, x].p), int(density[x, x].q))/2 for x, _, _ in entries)
    if collision != F(1, 2):
        raise ArithmeticError('global archive distinguishes the two experiments')
    return {'sites': 3, 'whole_Fock_dimension': 64,
            'local_CAR_basis_elements': len(basis), 'laws': ['C', 'E'],
            'coin_and_raw_slots_per_law': 2*RAW_SLOTS,
            'crossing_seam_local_dual_identity': True,
            'omitted_CAR_sign_countercontrol_detected': True,
            'reference_extension_by_operator_identity': True,
            'archive_countercontrol_masks': list(alternatives),
            'original_remote_collision_probability': str(collision),
            'global_archive_equality_claimed': False}


def centered_interaction_control():
    h, _, _ = exact_source(2, ((0, 1, 1),), F(3, 2), F(4))
    eigenvalues = tuple(h.eigenvals())
    ground_energy = min(eigenvalues)
    projection = s.eye(16)
    for energy in eigenvalues:
        if energy != ground_energy:
            projection = projection*(h-energy*s.eye(16))/(ground_energy-energy)
    projection = projection.applyfunc(s.cancel)
    if projection*projection != projection or s.trace(projection) != 1:
        raise ArithmeticError('full centered source ground projector failed')
    a, b = annihilator(2, 0), annihilator(2, 1)
    n = b.H*b
    derivative = s.I*((h*a-a*h)*a.H+a.H*(h*a-a*h))
    expected = -4*s.I*(n-s.eye(16)/2)
    if derivative != expected or s.trace(projection*n) != s.Rational(1, 2):
        raise ArithmeticError('centered local interaction identity failed')
    mean = s.trace(projection*derivative)
    coefficient = s.trace(projection*derivative.H*derivative)-s.conjugate(mean)*mean
    if coefficient != 4:
        raise ArithmeticError('centered same-ground non-Gaussian coefficient differs')
    return {'sites': 2, 'whole_Fock_dimension': 16, 'coupling': '3/2', 'contact': '4',
            'ground_H_energy': str(ground_energy), 'full_ground_rank': 1,
            'centered_clock_differs_from_H_only_by_scalar': True,
            'one_site_opposite_occupation': '1/2',
            'exact_four_point_variance_quadratic_coefficient': str(coefficient),
            'finite_dimer_is_infinite_ground_identification': False}


def original_weyl_bare_repair_controls():
    """Use the literal Weyl action and the unconjugated full-domain instrument."""
    def exact(matrix, scale=1):
        return (s.Matrix(matrix.real.toarray().tolist())
                +s.I*s.Matrix(matrix.imag.toarray().tolist()))/scale
    identity=s.eye(64);c=s.Rational(3,2);u=s.Integer(4)
    a0,b0=annihilator(3,0),annihilator(3,1)
    basis=tuple(left*right for left in(identity,a0,a0.H,a0.H*a0)
                for right in(identity,b0,b0.H,b0.H*b0))
    for law in('C','E'):
        if any(repair_dual(observable,sites=3,seam=(1,2),law=law)!=observable
               for observable in basis):
            raise ArithmeticError('bare original repair changed the disjoint local CAR algebra')
    rows=[]
    for axis in range(3):
        sides=tuple(3 if j==axis else 1 for j in range(3))
        source=symmetry.operators(FieldGeometry(sides,1.5,4.0))
        clock=exact(source.K,source.scale)
        for bond in source.bonds:
            spectrum=exact(bond,source.scale).eigenvals()
            if spectrum!={-c:4,-c/2:16,s.Integer(0):24,c/2:16,c:4}:
                raise ArithmeticError('literal Weyl bond norm differs from the local theorem')
        a,b=annihilator(3,2),annihilator(3,3)
        eta=b.H*b-identity/2
        derivative=clock*a.H-a.H*clock
        second=derivative.H*derivative+derivative*derivative.H
        expected=(u*u+2*c*c)/4
        if second!=expected*identity:
            raise ArithmeticError('original centered Weyl charged second moment differs')
        anti_derivative=s.I*((clock*a-a*clock)*a.H+a.H*(clock*a-a*clock))
        if anti_derivative!=-s.I*u*eta:
            raise ArithmeticError('original Weyl same-site interaction derivative differs')
        rows.append({'axis':axis,'source_geometry':list(sides),
                     'full_Fock_dimension':64,'coupling':'3/2','contact':'4',
                     'literal_H_and_K_reconstructed_by_source':True,
                     'each_bond_norm':'3/2','probe_site_degree':2,
                     'charged_second_moment_sum_rule':str(expected),
                     'centered_interaction_derivative_exact':True})
    return {'source_cases':rows,'bare_disjoint_physical_seam':[1,2],
            'observed_site':0,'local_CAR_basis_elements':16,
            'both_laws_all_coin_raw_slots_retained':True,
            'repair_Kraus_maps_conjugated':False,
            'hidden_alignment_or_chart_covariance_required':False,
            'Hamiltonian_clock_and_attempt_schedule_selected_by_repair':False}


def certificate():
    def encode(value):
        if isinstance(value, F): return str(value)
        if isinstance(value, dict): return {key: encode(item) for key, item in value.items()}
        if isinstance(value, list): return [encode(item) for item in value]
        return value
    return encode({'finite_locality': finite_locality_controls(),
                   'centered_interaction': centered_interaction_control(),
                   'original_Weyl_bare_repair_local_instance':original_weyl_bare_repair_controls(),
                   'local_history_plans': [plan_local_history(0, 1, F(1, 8), 1, 2, 8, eta)
                                           for eta in (F(1, 10), F(1, 1000))],
                   'interaction_example': interaction_witness(F(3, 2), 4, F(1, 1000)),
                   'chart': 'original Weyl action with bare repair, or separately signed Hubbard action with same-chart repair; whole-protocol transport conjugates both',
                   'original_Weyl_bare_repair_local_history_instantiated':True,
                   'bare_repairs_identified_under_site_dependent_frame':False,
                   'same_centered_thermodynamic_ground_subsequence': True,
                   'effective_ground_subsequence_preparation_modulus': False,
                   'odd_matter_and_even_spin_share_full_CAR_action': True,
                   'global_raw_archive_equality': False,
                   'source_selected_capacity_or_repair_duration': False,
                   'relativistic_continuum_or_physical_particle_assignment': False})
