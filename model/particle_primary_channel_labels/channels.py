"""Exact channel restrictions of the full physical Wilson source.

The cubic average is a distinct finite observable on the same source.
Neither its cubic label nor its positive finite band identifies a mass.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract
from model.particle_primary_field_algebra import algebra
from model.particle_primary_survival import survival
from model.particle_filtered_probe import filter as apparatus
from model.particle_local_wilson_density import density
from model.native_quantum_gauss import algebra as native
from model.interacting_extension.chiral_lattice import source_data
from model.cubic_particle_channels import group

SOURCE_INPUTS = (
    'model/particle_common_contract/PROOF.md',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_compact_mass_probe/PROOF.md',
    'model/particle_local_wilson_density/density.py',
    'model/native_quantum_gauss/algebra.py',
    'model/cubic_particle_channels/group.py',
)


def _sym(value):
    return s.Rational(value.a.numerator, value.a.denominator) + s.sqrt(2)*s.Rational(value.b.numerator, value.b.denominator)


def spin_lift(q):
    """Binary cubic lift in the actual four-component Clifford chart."""
    g = density.gamma_data()
    w, x, y, z = map(_sym, q)
    return apparatus.clean(w*s.eye(4)-x*g[1]*g[2]-y*g[2]*g[0]-z*g[0]*g[1])


def rotate_cell(cell, side, rotation):
    doubled = tuple(2*x-(side-1) for x in cell)
    raw = tuple(sum(int(rotation[i, j])*doubled[j] for j in range(3))+side-1 for i in range(3))
    if any(x % 2 for x in raw):
        raise ArithmeticError('cube rotation left the lattice')
    return tuple(x//2 for x in raw)


def plaquette_vertices(anchor, axes):
    return tuple(sorted(tuple(anchor[k]+int(k == axes[0])*i+int(k == axes[1])*j
                              for k in range(3)) for i, j in product(range(2), repeat=2)))


def rotate_plaquette(vertices, side, rotation):
    return tuple(sorted(rotate_cell(v, side, rotation) for v in vertices))


@lru_cache(maxsize=1)
def native_number_control():
    """Check both actual internal basis charts, without identifying their indices."""
    currents = native.native_currents()
    mask_rows = []
    for center, indices in native.native_blocks():
        charge = int(-2*s.Rational(center))
        projector = s.diag(*(int(i in indices) for i in range(15)))
        if any(apparatus.clean(projector*t-t*projector) != s.zeros(15) for t in currents):
            raise ArithmeticError('a native generator mixes occupation blocks')
        mask_rows.append((charge, len(indices)))
    charges = source_data.load()['charges']
    expected = ((-4, 3), (-3, 2), (1, 6), (2, 3), (6, 1))
    if tuple(sorted(mask_rows)) != expected or tuple((q, charges.count(q)) for q, _ in expected) != expected:
        raise ArithmeticError('the two actual source charts have different native blocks')
    packet_projectors = tuple(s.diag(*(int(x == q) for x in charges)) for q, _ in expected)
    # A non-diagonal genuine group element checks the packet chart as well.
    link = density.native_weak_link(0, s.Rational(3, 5), s.Rational(4, 5))
    if any(apparatus.clean(p*link-link*p) != s.zeros(15) for p in packet_projectors):
        raise ArithmeticError('the native hopping representation mixes blocks')
    return {'blocks': [{'Q6': q, 'internal_dimension': d, 'CAR_modes_per_cell': 4*d}
                       for q, d in expected],
            'native_generators_checked': len(currents),
            'both_mask_and_packet_charts_checked': True,
            'conserved_global_occupation_counts': 5,
            'full_H_has_pair_creation_or_interblock_terms': False,
            'physical_zero_external_charge_constraint': 'sum_q q*N_q=0',
            'baryon_lepton_or_species_identification': None}


@lru_cache(maxsize=1)
def rotation_control():
    g = density.gamma_data(); beta = g[3]
    for q in group.elements():
        lift = spin_lift(q); rotation = group.rotation(q)
        if apparatus.clean(lift.H*lift) != s.eye(4) or apparatus.clean(lift*beta*lift.H) != beta:
            raise ArithmeticError('the native spin lift changed beta or unitarity')
        for j in range(3):
            target = sum((int(rotation[k, j])*g[k] for k in range(3)), s.zeros(4))
            if apparatus.clean(lift*g[j]*lift.H) != target:
                raise ArithmeticError('native Wilson Clifford covariance failed')
    return {'binary_rotations': 48, 'proper_spatial_rotations': 24,
            'actual_Dirac_matrices_checked': True,
            'central_spin_lift': '-I_4; CAR lift is (-1)^N_F',
            'rotation_domain': 'physical interior graph with ghost links fixed to constants',
            'continuous_rotation_or_Lorentz_spin_proved': False}


def orbit_control(side=8, anchor=None, axes=(0, 1)):
    if type(side) is not int or side < 4 or side % 2:
        raise ValueError('an even side at least four is required')
    anchor = (side//2,)*3 if anchor is None else tuple(anchor)
    vertices = plaquette_vertices(anchor, axes)
    if any(not 0 <= x < side for v in vertices for x in v):
        raise ValueError('the plaquette must be inside the source box')
    orbit = tuple(sorted({rotate_plaquette(vertices, side, group.rotation(q)) for q in group.elements()}))
    multiplicities = group.decompose(lambda q: group.Q(sum(
        rotate_plaquette(p, side, group.rotation(q)) == p for p in orbit)))
    return {'vertices': vertices, 'orbit': orbit, 'orbit_size': len(orbit),
            'coefficient_representation_irreps': multiplicities,
            'A1_average_coefficient': str(F(1, len(orbit))),
            'original_single_plaquette_is_cubic_scalar': len(orbit) == 1,
            'representation_injective_on_physical_blank': True}


def boundary_graph_control(side=3):
    """Exact geometric control of the asymmetric allocation and physical graph."""
    cells = tuple(product(range(side), repeat=3))
    all_edges = tuple((x, tuple(x[k]+int(k == j) for k in range(3))) for x in cells for j in range(3))
    ghosts = [y for _, y in all_edges if any(k == side for k in y)]
    if len(ghosts) != len(set(ghosts)) or len(ghosts) != 3*side*side:
        raise ArithmeticError('ghost endpoints are not one-link vertices')
    interior = {tuple(sorted((x, y))) for x, y in all_edges if y not in ghosts}
    for q in group.elements():
        mapped = {tuple(sorted((rotate_cell(x, side, group.rotation(q)), rotate_cell(y, side, group.rotation(q)))))
                  for x, y in interior}
        if mapped != interior:
            raise ArithmeticError('physical interior graph is not cube covariant')
    return {'side': side, 'owned_links': len(all_edges), 'one_link_ghost_vertices': len(ghosts),
            'interior_links': len(interior), 'interior_graph_cube_covariant': True,
            'ambient_owned_link_allocation_claimed_cube_covariant': False}


def car_control():
    """Separate exact two-site/two-block CAR test, including forbidden mutations."""
    size = 16
    destroy = []
    for j in range(4):
        c = s.zeros(size)
        for mask in range(size):
            if mask & (1 << j):
                c[mask ^ (1 << j), mask] = (-1)**((mask & ((1 << j)-1)).bit_count())
        destroy.append(c)
    numbers = [c.H*c for c in destroy]
    blocks = (numbers[0]+numbers[2], numbers[1]+numbers[3])
    h = numbers[0]-numbers[2]+2*numbers[1]-2*numbers[3]
    for j in range(2):
        h += (j+1)*(destroy[j].H*destroy[j+2]+destroy[j+2].H*destroy[j])
    conserved = [h*n-n*h == s.zeros(size) for n in blocks]
    cross = destroy[0].H*destroy[1]+destroy[1].H*destroy[0]
    pair = destroy[0].H*destroy[1].H+destroy[1]*destroy[0]
    total = sum(numbers, s.zeros(size))
    return {'dimension': size, 'number_blocks': blocks, 'H': h,
            'separate_block_commutators_zero': conserved,
            'cross_block_mutation_preserves_total': cross*total == total*cross,
            'cross_block_mutation_breaks_individual': cross*blocks[0] != blocks[0]*cross,
            'pair_mutation_breaks_total': pair*total != total*pair,
            'actual_primary_spatial_matrix': False}


def selected_channel(member):
    contract.verify_member(member)
    witness = survival.selected_witness(member)
    source = member['source']; side = source['side']; a = F(source['spacing'])
    probe = witness['probe']
    orbit = orbit_control(side, probe['nearest_lower_owner_cell'], tuple(probe['plaquette_axes']))
    if orbit['orbit_size'] != 24:
        raise ArithmeticError('the authenticated selected plaquette orbit changed')
    width = F(witness['source_spectral_width_upper'])
    floor = witness['selected_density_floor']['negative_exponent']
    k = F(witness['blank_commutator_norm_squared'])*a*a/orbit['orbit_size']
    multiplier = F(witness['triangle_multiplier_lower_on_full_finite_spectrum'])
    contrast = F(11, 24)*multiplier*k/(a*a*width*width)
    band = (multiplier*k-1)/(2*a*a*width*width)
    finite = band-contrast/4
    if finite <= 0 or width < 1/a:
        raise ArithmeticError('cubic scalar finite-band lower bound is not positive')
    return contract.serial({'schema': 'oph.primary_channel_labels.v1',
        'member_id': member['member_id'], **member['ids'],
        'original_probe': probe, 'native_occupation_control': native_number_control(),
        'delta_N_q': {str(row['Q6']): 0 for row in native_number_control()['blocks']},
        'gauge_invariant': True, 'CAR_even': True,
        'original_probe_cubic_orbit': orbit,
        'cubic_scalar_operator': 'A_A1=(1/24) sum_p_in_orbit A_p',
        'cubic_scalar_is_distinct_observable': True,
        'cubic_scalar_norm_upper': 1,
        'cubic_scalar_blank_commutator_norm_squared': k/(a*a),
        'filter_half_width': 1/width, 'Ramsey_source_time': 1/width,
        'continuous_original_contrast_lower': survival.dyadic(contrast, floor),
        'finite_apparatus_original_contrast_lower': survival.dyadic(contrast/2, floor),
        'original_positive_band': {'lower': survival.dyadic(1/a, (floor+1)//2), 'upper': width},
        'ideal_original_band_weight_lower': survival.dyadic(band, floor),
        'finite_apparatus_original_band_weight_lower': survival.dyadic(finite, floor),
        'apparatus': survival.factored_apparatus(contrast, floor),
        'direct_compact_four_outcome_service_admissible': True,
        'source_marginal_rotation_invariant': True,
        'canonical_reference_rotation_rule': 'R_source tensor conjugate(R)_reference',
        'source_only_rotation_fixes_purification_claimed': False,
        'reference_clock': 'idle H_R=0',
        'charge_conjugation_quantum_number': None,
        'continuum_spin': None, 'particle_species': None, 'rest_mass_or_GeV': None,
        'fixed_positive_band_survival_proved': False,
        'massless_spectrum_excluded': False,
        'individual_initial_occupation_sector_selected': False,
        'full_primary_response_matrix_executed': False})


def bind_member(member):
    value = selected_channel(member)
    return {'result': value, 'result_sha256': contract.digest(value)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('the channel, original source or retained weight differs')


def certificate():
    car = car_control()
    return {'schema': 'oph.primary_channel_family.v1',
            'native_numbers': native_number_control(), 'rotations': rotation_control(),
            'boundary_control': boundary_graph_control(),
            'CAR_control': {k: v for k, v in car.items() if k not in ('H', 'number_blocks')},
            'finite_members': [bind_member(contract.member(n)) for n in (1, 2)],
            'source_derived_physical_mass_table': []}
