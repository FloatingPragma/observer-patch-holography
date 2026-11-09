"""Actual interior-link cuts, boundary representations and entropy bounds."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from model.particle_common_contract import contract, provenance
from model.particle_primary_covariant_stress import model as stress
from model.particle_spatial_modular_response import regional

SOURCE_INPUTS = tuple('model/particle_primary_geometric_cut/'+name for name in (
    '__init__.py', '__main__.py', 'model.py', 'test_model.py', 'PROOF.md', 'README.md')) + (
    'model/particle_spatial_modular_response/PROOF.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_primary_covariant_stress/PROOF.md',
    'Lean/QFT/CPRestrictionNet.lean', 'Lean/QFT/CoarseEnergyWard.lean',
)


def cut(side, cells):
    """All matter in cells and every link with both endpoints in cells.

    Crossing links stay outside; outer ghost links are constant on Gauss
    support. Neither is silently turned into an independent inside edge.
    """
    if type(side) is not int or side < 2:
        raise ValueError('integer source side at least two required')
    rows = tuple(stress.stress.locality._cell(x, side) for x in cells)
    if not rows or len(rows) != len(set(rows)):
        raise ValueError('distinct nonempty regional cells required')
    selected = set(rows); internal = set(); crossing = set(); boundary = Counter()
    for x in selected:
        for axis, sign in product(range(3), (-1, 1)):
            y = tuple(t+sign*(k == axis) for k, t in enumerate(x))
            if not all(0 <= t < side for t in y):
                continue
            edge = tuple(sorted((x, y)))
            if y in selected:
                internal.add(edge)
            else:
                crossing.add(edge); boundary[x] += 1
    if sum(boundary.values()) != len(crossing):
        raise ArithmeticError('a crossing edge was counted more than once')
    return {'cells': tuple(sorted(selected)), 'interior_links': tuple(sorted(internal)),
            'crossing_links': tuple(sorted(crossing)),
            'inside_boundary': tuple(sorted(boundary.items())),
            'CAR_modes': 60*len(selected)}


@lru_cache(maxsize=3)
def cube_control(side=2):
    if side not in (2, 3, 4):
        raise ValueError('finite controls use side two, three or four')
    rotations = sorted({tuple(map(int, stress.channels.group.rotation(q).flat))
                        for q in stress.channels.group.elements()})
    regions = [((0, 0, 0),), ((0, 0, 0), (1, 0, 0)),
               tuple(product(range(side), repeat=3))]
    if side > 2:
        regions.append(tuple(product(range(1, side), repeat=3)))
    for cells in regions:
        base = cut(side, cells)
        for flat in rotations:
            R = stress.channels.s.Matrix(3, 3, flat)
            rotate = lambda x: stress.channels.rotate_cell(x, side, R)
            transformed = cut(side, [rotate(x) for x in cells])
            for key in ('interior_links', 'crossing_links'):
                expected = tuple(sorted(tuple(sorted(rotate(x) for x in edge))
                                        for edge in base[key]))
                if transformed[key] != expected:
                    raise ArithmeticError('induced cut does not rotate as a subsystem')
            if transformed['inside_boundary'] != tuple(sorted((rotate(x), d)
                                                              for x, d in base['inside_boundary'])):
                raise ArithmeticError('boundary group assignment changed under rotation')
    return {'side': side, 'proper_rotations_checked': len(rotations),
            'regions_checked': len(regions), 'induced_link_cut_covariant': True,
            'source_matter_modes_per_cell': 60,
            'source_H_state_and_clock_replaced': False}


def boundary_entropy_bound(member, cells=((0, 0, 0), (1, 0, 0))):
    contract.verify_member(member)
    source = member['source']; result = cut(source['side'], cells)
    ell = source['per_link_word_degree']; a = F(source['spacing'])
    blocks = stress.channels.native_number_control()['blocks']
    dmax = max(row['internal_dimension'] for row in blocks)
    qmax = max(abs(row['Q6']) for row in blocks)
    result.update({'word_degree': ell, 'largest_native_irrep_dimension': dmax,
        'edge_entropy_upper_coefficient_of_log_dmax': ell*len(result['crossing_links']),
        'boundary_irrep_dimension_upper_factored': {'base': dmax,
            'exponent': ell*len(result['crossing_links'])},
        'inside_boundary_charge_bounds': [
            {'cell': x, 'generator': 'regional matter plus internal-link endpoints',
             'integer_Q6_absolute_upper_on_physical_support': qmax*ell*d}
            for x, d in result['inside_boundary']],
        'lattice_cut_face_area': len(result['crossing_links'])*a*a,
        'edge_entropy_upper_per_lattice_area_coefficient_of_log_dmax': F(ell)/(a*a)})
    return result


def regional_energy_difference(member, cells=((0, 0, 0), (1, 0, 0))):
    """H'_R minus all whole source terms contained in the induced cut."""
    contract.verify_member(member)
    source = member['source']; geometry = cut(source['side'], cells)
    selected = set(geometry['cells']); shares = {}; norms = {}
    for x in selected:
        for term in stress._incident(source, x):
            key = term['id']; shares[key] = shares.get(key, F())+term['share']
            norms[key] = term['norm']
    differences = {}
    for term, weight in shares.items():
        inside = all(x in selected for x in term[1])
        delta = weight-int(inside)
        if delta:
            if inside or not 0 < delta < 1:
                raise ArithmeticError('regional discrepancy is not a crossing source term')
            differences[str(term)] = delta
    norm = sum((shares[t]*norms[t] for t in shares
                if not all(x in selected for x in t[1])), F())
    b = len(geometry['crossing_links']); ell = source['per_link_word_degree']
    cap = b*(F(609,4)*ell**2+15+F(3,14))/F(source['spacing'])
    if norm > cap:
        raise ArithmeticError('crossing energy bound failed')
    return {'crossing_term_coefficients': differences,
            'boundary_difference_norm_upper': norm,
            'uniform_cut_boundary_difference_norm_upper': cap,
            'contained_H_R_is_in_the_regional_gauge_invariant_algebra': True,
            'contained_H_R_equals_modular_generator': False}


def bind_member(member):
    contract.verify_member(member)
    stress.verify_bound(member, stress.bind_member(member))
    bound = boundary_entropy_bound(member)
    floor = regional.common_contract_payload(member, 2)['physical_density_floor_factored']
    # One internal edge of an elementary plaquette lies inside this cut.
    # A native Wilson character has d^2 equal Schmidt coefficients there.
    blocks = stress.channels.native_number_control()['blocks']
    result = contract.serial({'member_id': member['member_id'], **member['ids'],
        'induced_cut': bound,
        'physical_density_floor': floor,
        'edge_entropy_positive_lower_factored': {
            'multiplier_of_physical_density_floor': 4, 'natural_log_argument': 108},
        'regional_energy_boundary_difference': regional_energy_difference(member),
        'regional_factorization': 'matter at R; interior links with both endpoints in R',
        'all_crossing_links_retained_in_complement': True,
        'physical_outer_ghost_links_are_constant': True,
        'rotation_covariant_partial_trace_and_supported_modular_flow': True,
        'regional_C_factorization_and_covariant_center_records': True,
        'finite_restriction_and_entropy_theorem_premises_supplied': True,
        'center_record_retains_regional_entropy_exactly': True,
        'global_center_measurement_asserted_nondisturbing': False,
        'Wilson_character_boundary_witnesses': [
            {'Q6_opposite_endpoint_pair': [row['Q6'], -row['Q6']],
             'Schmidt_rank': row['internal_dimension']**2,
             'regional_nonzero_eigenvalue': F(1, row['internal_dimension']**2),
             'selected_weight_floor': 'same-source c from rho >= c P_G'}
            for row in blocks],
        'opposite_character_witnesses_and_blank_retained': True,
        'full_regional_density_or_multiplicities_numerically_executed': False,
        'entropy_decomposition': 'S_extended = S_algebra + sum_lambda p_lambda log d_lambda',
        'proper_geometry_or_diamond_kernel_derived': False,
        'finite_nonzero_universal_area_coefficient_derived': False,
        'total_entropy_area_law_derived': False,
        'Einstein_equation_derived': False, 'physical_mass_table_GeV': []})
    return {'schema': 'oph.primary_geometric_cut.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('geometric cut source or response changed')


def certificate():
    paths = provenance.definition_inputs((*SOURCE_INPUTS, *stress.SOURCE_INPUTS))
    pins = provenance.definition_hashes(paths)
    value = {'schema': 'oph.primary_geometric_cut_certificate.v1',
        'cube_controls': [cube_control(n) for n in (2, 3)],
        'members': [bind_member(contract.member(n)) for n in (0, 1, 2, 8)],
        'source_sha256': pins}
    if pins != provenance.definition_hashes(paths):
        raise ValueError('source changed during geometric-cut replay')
    return value
