"""Exact allocation, boundary improvement and symmetry of source energy."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

from model.particle_common_contract import contract, provenance
from model.particle_primary_stress import stress
from model.particle_primary_channel_labels import channels
from model.particle_primary_charge_conjugation import model as conjugation
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_primary_covariant_stress/__init__.py',
    'model/particle_primary_covariant_stress/__main__.py',
    'model/particle_primary_covariant_stress/model.py',
    'model/particle_primary_covariant_stress/test_model.py',
    'model/particle_primary_covariant_stress/PROOF.md',
    'model/particle_primary_covariant_stress/README.md',
    'model/particle_primary_stress/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/particle_primary_charge_conjugation/PROOF.md',
    'model/particle_spatial_modular_response/PROOF.md',
    'Lean/QFT/CoarseEnergyWard.lean',
)


def _edge(x, axis):
    return (x, tuple(t + (i == axis) for i, t in enumerate(x)))


def _term(kind, vertices):
    return (kind, tuple(sorted(vertices)))


def _incident(source, cell):
    """Use geometric vertices as labels; retain actual outgoing-link owners."""
    L = source['side']; x = stress.locality._cell(cell, L)
    a = F(source['spacing']); ell = source['per_link_word_degree']
    electric = F(source['electric_coefficient']) * F(29, 10) * ell**2 / a
    rows = [{'id': _term('onsite_matter', (x,)), 'vertices': (x,),
             'owner_support': (x,), 'norm': F(90)/a, 'share': F(1)}]
    for axis in range(3):
        for delta in (-1, 0):
            anchor = tuple(t + delta*(i == axis) for i, t in enumerate(x))
            vertices = _edge(anchor, axis)
            if not all(stress._inside(v, L) for v in vertices):
                continue
            rows.append({'id': _term('electric', vertices), 'vertices': vertices,
                         'owner_support': (anchor,), 'norm': electric, 'share': F(1, 2)})
            rows.append({'id': _term('bond', vertices), 'vertices': vertices,
                         'owner_support': vertices, 'norm': F(30)/a, 'share': F(1, 2)})
    for i, j in combinations(range(3), 2):
        for di, dj in product((-1, 0), repeat=2):
            anchor = tuple(t + di*(k == i) + dj*(k == j) for k, t in enumerate(x))
            vertices = tuple(tuple(t + u*(k == i) + v*(k == j)
                                   for k, t in enumerate(anchor))
                             for u, v in product((0, 1), repeat=2))
            if not all(stress._inside(v, L) for v in vertices):
                continue
            owners = (anchor, _edge(anchor, i)[1], _edge(anchor, j)[1])
            rows.append({'id': _term('plaquette', vertices), 'vertices': vertices,
                         'owner_support': owners,
                         'norm': F(source['magnetic_coefficient'])/a, 'share': F(1, 4)})
    return rows


def cell_energy(member, cell=(0, 0, 0)):
    contract.verify_member(member)
    x = stress.locality._cell(cell, member['source']['side'])
    rows = _incident(member['source'], x)
    owners = tuple(sorted({v for row in rows for v in row['owner_support']}))
    if any(sum(abs(u-v) for u, v in zip(x, y)) > 2 for y in owners):
        raise ArithmeticError('covariant energy support exceeded radius two')
    return {'cell': x, 'terms': rows, 'owner_support': owners,
            'norm_upper': sum((row['norm']*row['share'] for row in rows), F())}


def _small_source(side):
    L = ground._i(side, 'small control side', 2)
    if L > 4:
        raise ValueError('enumerated controls use side at most four')
    return {'side': L, 'spacing': '1', 'per_link_word_degree': 1,
            'electric_coefficient': '105', 'magnetic_coefficient': '1/7'}


def _allocations(side):
    source = _small_source(side)
    return {x: {r['id']: r['share'] for r in _incident(source, x)}
            for x in product(range(side), repeat=3)}


def _old_allocations(side):
    source = _small_source(side); result = {}
    for x in product(range(side), repeat=3):
        terms = {}
        for row in stress._cell_terms(source, x):
            kind, anchor, axes = row['id']
            if kind == 'onsite_matter': vertices = (anchor,)
            elif kind in ('electric', 'bond'):
                vertices = _edge(anchor, axes[0])
                if any(not stress._inside(v, side) for v in vertices):
                    continue  # actual ghost Gauss projection makes this electric term zero
            else:
                i, j = axes
                vertices = tuple(tuple(t + u*(k == i) + v*(k == j)
                                       for k, t in enumerate(anchor))
                                 for u, v in product((0, 1), repeat=2))
            terms[_term(kind, vertices)] = row['cell_share']
        result[x] = terms
    return result


def _routes(term):
    kind, vertices = term
    if kind == 'electric':
        yield vertices[0], vertices[1], F(1, 2)
    elif kind == 'plaquette':
        a = tuple(min(v[i] for v in vertices) for i in range(3))
        axes = [i for i in range(3) if len({v[i] for v in vertices}) == 2]
        b = _edge(a, axes[0])[1]; c = _edge(a, axes[1])[1]
        d = _edge(b, axes[1])[1]
        yield a, b, F(1, 24)
        yield a, c, F(1, 24)
        yield b, d, F(1, 8)
        yield c, d, F(1, 8)


@lru_cache(maxsize=3)
def allocation_control(side=3):
    new = _allocations(side); old = _old_allocations(side)
    terms = set().union(*(set(v) for v in new.values()))
    if terms != set().union(*(set(v) for v in old.values())):
        raise ArithmeticError('physical action inventory changed')
    for term in terms:
        if sum((v.get(term, 0) for v in new.values()), F()) != 1:
            raise ArithmeticError('symmetric shares fail to reproduce the action')
        divergence = defaultdict(F)
        for u, v, weight in _routes(term):
            if sum(abs(a-b) for a, b in zip(u, v)) != 1:
                raise ArithmeticError('improvement must use only nearest edges')
            divergence[u] -= weight; divergence[v] += weight
        for x in new:
            if new[x].get(term, 0) - old[x].get(term, 0) != divergence[x]:
                raise ArithmeticError('operator-valued boundary improvement failed')
    rotations = {tuple(map(int, channels.group.rotation(q).flat)) for q in channels.group.elements()}
    for flat in rotations:
        R = channels.s.Matrix(3, 3, flat)
        for x, weights in new.items():
            y = channels.rotate_cell(x, side, R)
            rotated = {_term(kind, [channels.rotate_cell(v, side, R) for v in vertices]): w
                       for (kind, vertices), w in weights.items()}
            if rotated != new[y]:
                raise ArithmeticError('geometric source-energy allocation is not cubic covariant')
    return {'side': side, 'cells_checked': side**3, 'proper_rotations_checked': len(rotations),
            'physical_terms_checked': len(terms),
            'term_counts': {k: sum(t[0] == k for t in terms)
                            for k in ('onsite_matter', 'electric', 'bond', 'plaquette')},
            'exact_partition_and_nearest_edge_improvement': True,
            'dangling_electric_terms_zero_only_after_ghost_Gauss_projection': True}


def regional_improvement_control(side=3, region=((0, 0, 0),)):
    new = _allocations(side); old = _old_allocations(side)
    cells = tuple(stress.locality._cell(x, side) for x in region)
    if not cells or len(set(cells)) != len(cells):
        raise ValueError('nonempty distinct regional cells required')
    selected = set(cells); terms = set().union(*(set(v) for v in new.values()))
    differences = {}; boundary = {}; crossing = 0
    for term in sorted(terms):
        difference = sum((new[x].get(term, 0)-old[x].get(term, 0) for x in selected), F())
        flow = F()
        for u, v, weight in _routes(term):
            if (u in selected) != (v in selected):
                flow += weight*((v in selected)-(u in selected)); crossing += 1
        if difference != flow:
            raise ArithmeticError('regional improvement contains a volume contribution')
        if difference:
            differences[str(term)] = difference; boundary[str(term)] = flow
    return {'region_cells': len(cells), 'crossing_routes': crossing,
            'energy_difference_coefficients': differences,
            'inward_boundary_flow_coefficients': boundary}


def bind_member(member):
    contract.verify_member(member)
    conjugation.verify_bound(member, conjugation.bind_member(member))
    x = cell_energy(member)
    a = F(member['source']['spacing']); ell = member['source']['per_link_word_degree']
    cap = (F(180)+F(1827, 2)*ell**2+F(3, 7))/a
    if x['norm_upper'] > cap:
        raise ArithmeticError('source cap does not bound the physical corner energy')
    result = contract.serial({'member_id': member['member_id'], **member['ids'],
        'partition': 'onsite matter + half endpoint electric and bonds + quarter vertex plaquettes',
        'physical_source_H_and_energy_zero_unchanged': True, 'source_clock_unchanged': True,
        'ghost_Gauss_support_required': True, 'source_corner_energy': x,
        'energy_norm_upper': cap, 'pair_current_norm_upper': 2*cap**2,
        'support_radius_graph': 2, 'current_zero_beyond_graph_distance': 4,
        'proper_cube_energy_covariance': 'R h_x R^-1 = h_(R x)',
        'proper_cube_current_covariance': 'R J_xy R^-1 = J_(R x,R y)',
        'declared_scalar_control_covariance': 'R H_f R^-1 = H_(f composed with R^-1)',
        'C_even_energy_and_current': True,
        'exact_boundary_improvement_from_owner_allocation': True,
        'finite_Ward_identity_premises_supplied': True,
        'complete_square_root_energy_and_current_records': True,
        'source_and_idle_quantum_reference_retained': True,
        'rotation_invariant_regional_spectra_and_complete_record_laws': True,
        'stationary_source_net_regional_flux': 0,
        'stationary_source_each_pair_flux_claimed_zero': False,
        'C_odd_source_tangent_has_zero_first_order_energy_and_regional_entropy': True,
        'regional_entropy_support_uses_actual_partial_trace': True,
        'owner_cell_partial_trace_is_automatically_rotation_covariant': False,
        'proper_metric_or_Lorentz_stress_derived': False,
        'scalar_control_is_a_derived_proper_lapse': False,
        'Einstein_modular_kernel_and_area_law_derived': False,
        'refinement_uniform_stress_norm': False, 'physical_mass_table_GeV': []})
    return {'schema': 'oph.primary_covariant_stress.v1', 'result': result,
            'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('covariant stress source or result changed')


def certificate():
    paths = provenance.definition_inputs((*SOURCE_INPUTS, *stress.SOURCE_INPUTS,
                                         *conjugation.SOURCE_INPUTS))
    pins = provenance.definition_hashes(paths)
    value = {'schema': 'oph.primary_covariant_stress_certificate.v1',
             'allocation_controls': [allocation_control(L) for L in (2, 3)],
             'regional_control': contract.serial(regional_improvement_control()),
             'members': [bind_member(contract.member(n)) for n in (0, 1, 2, 8)],
             'source_sha256': pins}
    if provenance.definition_hashes(paths) != pins:
        raise ValueError('source changed during replay')
    return value
