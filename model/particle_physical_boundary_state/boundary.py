"""Actual cut-charge sectors of the common Gauss-constrained source state.

The complete central-charge menu is finite. Spatial matrices and all
nonabelian multiplicities need not be allocated to prove supported sectors.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.interacting_extension.chiral_lattice import source_data
from model.native_quantum_gauss import algebra as native
from model.particle_common_contract import contract
from model.particle_spatial_modular_response import regional
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_common_contract/contract.py',
    'model/particle_common_contract/PROOF.md',
    'model/particle_gauss_ground_state/gauss.py',
    'model/particle_gauss_ground_state/PROOF.md',
    'model/particle_spatial_ground_state/ground.py',
    'model/particle_spatial_ground_state/PROOF.md',
    'model/particle_spatial_modular_response/regional.py',
    'model/particle_spatial_modular_response/PROOF.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/native_quantum_gauss/algebra.py',
    'model/particle_field_preparation/PROOF.md',
)


def charge_spectrum(degree):
    """Exact central right-endpoint spectrum of the native word space V_ell.

    Every sum occurs: multiply diagonal coefficient entries, each nonzero
    at the identity. Identity padding includes every smaller word degree.
    """
    degree = ground._i(degree, 'per-link word degree', 0)
    charges = {b['q6Y'] for b in source_data.load()['blocks']}
    alphabet = {0} | charges | {-q for q in charges}
    result = {0}
    for _ in range(degree):
        result = {left + right for left in result for right in alphabet}
    return tuple(sorted(result))


def charge_projectors(matrix, spectrum):
    """Exact complete spectral projectors for a finite central charge matrix.

    Exhaustiveness is checked, not inferred from a proposed outcome list.
    This routine accepts exact matrices and does not allocate a spatial one.
    """
    q = s.ImmutableSparseMatrix(matrix)
    values = tuple(spectrum)
    if (q.rows != q.cols or q != q.H or not values
            or len(set(values)) != len(values)
            or any(type(x) is not int for x in values)
            or any(x.has(s.Float) or x.free_symbols for x in q)):
        raise ValueError('exact Hermitian matrix and distinct integer spectrum required')
    ident = s.eye(q.rows)
    rows = []
    diagonal = all(i == j for i, j in q.todok())
    for value in values:
        if diagonal:
            entries = []
            for i in range(q.rows):
                entry = s.prod((q[i, i]-other)/s.Integer(value-other)
                               for other in values if other != value)
                entries.append(s.cancel(entry))
            p = s.ImmutableSparseMatrix(s.diag(*entries))
        else:
            p = s.eye(q.rows)
            for other in values:
                if other != value:
                    p = p*(q-other*ident)/s.Integer(value-other)
            p = s.ImmutableSparseMatrix(p.applyfunc(s.simplify))
        if p != p.H or p*p != p or q*p != value*p:
            raise ValueError('the proposed menu does not resolve the supplied charge')
        rows.append((value, p))
    if sum((p for _, p in rows), s.zeros(q.rows)) != ident:
        raise ValueError('the proposed charge spectrum omits an outcome')
    if any(p*r != s.zeros(q.rows) for i, (_, p) in enumerate(rows)
           for _, r in rows[:i]):
        raise ArithmeticError('charge outcomes are not mutually exclusive')
    return tuple(rows)


@lru_cache(maxsize=1)
def native_degree_one_control():
    """The complete119-dimensional one-link V1, including its constant word."""
    diagonal = [0]
    labels = [('constant', 0, 0)]
    for block in source_data.load()['blocks']:
        q = block['q6Y']; d = block['dimension']
        for conjugate in (False, True):
            endpoint = q if conjugate else -q
            for row, column in product(range(d), repeat=2):
                diagonal.append(endpoint)
                labels.append((('conjugate_' if conjugate else '')+str(q), row, column))
    matrix = s.ImmutableSparseMatrix(s.diag(*diagonal))
    pvm = charge_projectors(matrix, charge_spectrum(1))
    expected = 1+2*sum(b['dimension']**2 for b in source_data.load()['blocks'])
    if len(diagonal) != expected or expected != 119:
        raise ArithmeticError('the complete native degree-one Haar basis changed')
    return {'matrix': matrix, 'labels': tuple(labels), 'projectors': pvm,
            'dimension': expected,
            'outcome_ranks': {str(q): int(s.trace(p)) for q, p in pvm},
            'all_native_coefficients_and_conjugates_retained': True}


@lru_cache(maxsize=1)
def source_irreducibility_control():
    """Exact commutants of all five actual native irreducible blocks."""
    currents = native.native_currents()
    expected = {row['q6Y']: row['dimension'] for row in source_data.load()['blocks']}
    results = []
    # This generator packet uses exterior-mask order, while source_data uses
    # the identified operator packet order. Its central current is -i Q6/2.
    for center, indices in native.native_blocks():
        charge = -2*s.Rational(center)
        if charge.q != 1 or expected.get(int(charge)) != len(indices):
            raise ArithmeticError('native central normalization or multiplet changed')
        d = len(indices); eye = s.eye(d)
        equations = []
        for current in currents:
            t = current.extract(indices, indices)
            if t+t.H != s.zeros(d):
                raise ArithmeticError('native current is not anti-Hermitian')
            equations.append(s.kronecker_product(t, eye)-s.kronecker_product(eye, t.T))
        commutant_dimension = d*d-s.Matrix.vstack(*equations).rank()
        if commutant_dimension != 1:
            raise ArithmeticError('a native character block is reducible')
        results.append({'q6Y': int(charge), 'dimension': d,
                        'commutant_dimension': commutant_dimension})
    return tuple(sorted(results, key=lambda row: row['q6Y']))


def loop_witnesses():
    """Actual physical plaquette characters across the origin-cell cut.

    Haar-normalized matrix coefficients give d² equal Schmidt weights.
    Each vector is in the declared spatial carrier at every selected level.
    """
    result = [{'label': 'vacuum', 'cut_charges': (0, 0, 0),
               'irrep_dimension': 1, 'Schmidt_rank': 1,
               'regional_Schmidt_eigenvalue': F(1), 'word_degree_required': 0}]
    for block in source_data.load()['blocks']:
        for conjugate in (False, True):
            q = (-1 if conjugate else 1)*block['q6Y']
            d = block['dimension']
            result.append({'label': ('conjugate_' if conjugate else 'native_')+str(block['q6Y']),
                'cut_charges': (-q, q, 0), 'irrep_dimension': d,
                'Schmidt_rank': d*d, 'regional_Schmidt_eigenvalue': F(1, d*d),
                'word_degree_required': 1})
    if len({row['cut_charges'] for row in result}) != 11:
        raise ArithmeticError('the eleven physical witness sectors overlap')
    return result


def common_payload(member):
    """Recompute the finite physical-state consequence of the selected member."""
    contract.verify_member(member)
    source = member['source']; state = member['state']
    side = source['side']; degree = source['per_link_word_degree']
    if side < 2 or degree < 1 or not state['source_only_Gauss_support']:
        raise ValueError('the selected source must contain an admitted physical plaquette')
    support = regional.common_contract_payload(member, region_cells=1)
    exponent = support['physical_density_floor_factored']['negative_exponent']
    spectrum = charge_spectrum(degree)
    witnesses = loop_witnesses()
    for row in witnesses:
        if any(q not in spectrum for q in row['cut_charges']):
            raise ArithmeticError('an admitted physical witness left the complete charge menu')
        row['selected_state_sector_probability_lower'] = {'base': 2, 'negative_exponent': exponent}
        row['regional_Schmidt_span_floor'] = {
            'base': 2, 'negative_exponent': exponent,
            'coefficient': row['regional_Schmidt_eigenvalue']}
    result = {
        'region': {'cells': [(0, 0, 0)], 'owned_link_directions': [0, 1, 2],
                   'cut_vertices': [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                   'plaquette_vertices': [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]},
        'source_only_Gauss_support': True,
        'outer_boundary': 'zero external charge at every owned-link ghost endpoint',
        'regional_boundary': 'inside charges paired with opposite outside charges',
        'state': {'density': state['density'], 'amplitude': state['amplitude'],
                  'beta': state['beta'], 'Taylor_order': state['Taylor_order'],
                  'dyadic_scale': state['dyadic_scale'], 'source_stationary': True,
                  'conjugate_reference_idle_H': member['reference']['idle_H'],
                  'direct_deterministic_preparation': True, 'postselection': False},
        'physical_density_floor_factored': support['physical_density_floor_factored'],
        'regional_support': 'support of Tr_out(P_G)',
        'complete_central_charge_measurement': {
            'single_endpoint_spectrum': spectrum, 'endpoint_count': 3,
            'outcome_menu': 'Cartesian product of the three endpoint spectra',
            'outcome_count': len(spectrum)**3,
            'projector_formula': 'E_q=product_j product_(r!=q_j) (Q_j-r I)/(q_j-r)',
            'isometry': 'V=sum_q E_q tensor |q>_record',
            'unitary_extension': 'U=sum_q E_q tensor S_q; S_q shifts the neutral record blank to q',
            'effect_sum': 'identity on the entire regional carrier',
            'commutes_with_full_source_Gauss_constraints': True,
            'all_zero_probability_outcomes_retained': True,
            'vacuum_and_other_outcomes_retained': True,
            'all_nonabelian_center_labels_resolved': False},
        'strictly_positive_selected_state_witnesses': witnesses,
        'witness_sectors_exhaust_physical_boundary_support': False,
        'full_spatial_multiplicities_computed': False,
        'full_spatial_density_or_gate_tape_executed': False,
        'efficient_volume_scaling_claimed': False,
        'finite_physical_preparation_including_boundary_charges_constructed': True,
        'refinement_compatible_state_family_proved': False,
        'interacting_relativistic_limit_or_particle_identification': False,
    }
    return contract.serial(result)


def bind_member(member):
    """Typed consumer; verify its result by rerunning this producer."""
    payload = common_payload(member)
    return {'schema': 'oph.physical_boundary_state.v1', 'member_id': member['member_id'],
            **member['ids'], 'result': payload, 'result_sha256': contract.digest(payload)}


def common_contract_payload(member):
    """Producer spelling used by the primary contract's claim registry."""
    return common_payload(member)


def common_contract_bound(member):
    """Bind after the primary contract registers this recomputed producer."""
    return contract.bind_result(member,
        'physical boundary-charge preparation and complete central readout',
        common_contract_payload(member))


def verify_bound(member, receipt):
    if receipt != bind_member(member):
        raise ValueError('physical boundary receipt differs from the selected source or result')


def certificate():
    member = contract.member(0); bound = bind_member(member)
    control = native_degree_one_control()
    return {'schema': 'oph.physical_boundary_certificate.v1',
            'selected_state_boundary': bound,
            'degree_one_link_control': {key: control[key] for key in
                ('dimension', 'outcome_ranks', 'all_native_coefficients_and_conjugates_retained')},
            'native_irreducibility_controls': list(source_irreducibility_control()),
            'selected_full_spatial_run_executed': False,
            'boundary_preparation_requirement_closed_at_finite_member': True,
            'joint_physical_refinement_proved': False}
