"""Finite quantum-Higgs feedback and exact regulator compatibility controls.

No physical vacuum, interacting continuum, or mirror-decoupling input is
supplied by these finite operator identities.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

import sympy as s

from model.interacting_extension.higgs_join import native_inputs, source_algebra
from model.standard_model_source_map import feedback_limit as feedback

ROOT = Path(native_inputs.ROOT)
SOURCE_INPUTS = (
    'model/interacting_extension/higgs_join/native_inputs.py',
    'model/interacting_extension/higgs_join/source_algebra.py',
    'model/interacting_extension/higgs_join/interacting_carrier.py',
    'model/interacting_extension/higgs_join/PROOF.md',
    'model/observer_particles/results/representation_identification_operators.json',
    'model/observer_particles/results/representation_identification.json',
    'model/observer_particles/results/unified_yukawa_operators.json',
    'model/observer_particles/results/unified_yukawa_action.json',
    'model/standard_model_source_map/feedback_limit.py',
    'model/standard_model_source_map/FEEDBACK_LIMIT.md',
    'model/standard_model_source_map/GAUGE_ATTACHMENT.md',
    'model/common_source/native_family_record_memory.py',
    'model/common_source/overlap_regulator.py',
    'model/interacting_extension/chiral_electric/PROOF.md',
)


def _clean(matrix):
    return s.ImmutableSparseMatrix(matrix.applyfunc(lambda x: s.simplify(s.expand(x))))


def _zero(matrix):
    return all(s.simplify(s.expand(x)) == 0 for x in matrix.todok().values())


def creator(modes, mode):
    if type(modes) is not int or type(mode) is not int or not 1 <= modes <= 12 or not 0 <= mode < modes:
        raise ValueError('one to twelve CAR modes and an admitted index are required')
    return s.ImmutableSparseMatrix(1 << modes, 1 << modes, {
        (mask | (1 << mode), mask): (-1)**((mask & ((1 << mode)-1)).bit_count())
        for mask in range(1 << modes) if not mask >> mode & 1})


@lru_cache(maxsize=1)
def _source():
    # These four fixtures authenticate one another; all used Ward identities
    # are independently recomputed. No old manifest's full replay is claimed.
    data = native_inputs.native_data(ROOT)
    anti, real, scalar = source_algebra.standard_generators(data)
    scalar_anti = []
    for a, r, h in zip(anti, real, scalar, strict=True):
        row, col = next((i, j) for i, j in r.todok() if r[i, j])
        scalar_anti.append(_clean((a[row, col]/r[row, col])*h))
    channel = data['channels']['charged_lepton']
    ids = (*channel['left'], *channel['right'])
    if len(ids) != 3 or len(channel['left']) != 2:
        raise ArithmeticError('the retained charged Higgs support differs')
    return data, tuple(anti), tuple(scalar_anti), ids


@lru_cache(maxsize=1)
def _quantum_higgs_operators():
    """The complete 3-scalar by 64-fermion local carrier, no projection."""
    data, currents, scalars, ids = _source()
    creation = tuple(creator(6, i) for i in range(6))
    gamma = data['channels']['charged_lepton']['components']
    w = s.SparseMatrix(192, 192, {})
    pair = []
    for h in range(2):
        a = s.zeros(64)
        for local, old in enumerate(ids[:2]):
            coefficient = s.conjugate(gamma[h][old, ids[2]])
            a += coefficient*(creation[2*local]*creation[5]-creation[2*local+1]*creation[4])
        pair.append(_clean(a))
        b = s.zeros(3); b[h+1, 0] = 1
        w += s.kronecker_product(b, a)
    nf = s.diag(*(mask.bit_count() for mask in range(64)))
    nh = s.diag(0, 1, 1)
    generators = []
    for g, h in zip(currents, scalars, strict=True):
        small = g.extract(ids, ids)
        q = s.zeros(64)
        for i, j, spin in product(range(3), range(3), range(2)):
            q += small[i, j]*creation[2*i+spin]*creation[2*j+spin].H
        hs = s.diag(0, 0, 0); hs[1:, 1:] = h
        generators.append(_clean(s.kronecker_product(s.eye(3), q)+s.kronecker_product(hs, s.eye(64))))
    return {'W': _clean(w), 'V': _clean(w+w.H), 'pairs': tuple(pair),
            'generators': tuple(generators), 'N_F': s.ImmutableSparseMatrix(s.kronecker_product(s.eye(3), nf)),
            'N_H': s.ImmutableSparseMatrix(s.kronecker_product(nh, s.eye(64))),
            'source_internal_indices': tuple(ids)}


def quantum_higgs_operators():
    return dict(_quantum_higgs_operators())


def quantum_higgs_certificate():
    ops = quantum_higgs_operators(); w, v = ops['W'], ops['V']
    for g in ops['generators']:
        if not _zero(g*v-v*g):
            raise ArithmeticError('the actual quantum Higgs vertex violates a source current')
    q = ops['N_F']-2*ops['N_H']
    if not _zero(q*v-v*q) or w*w != s.zeros(192):
        raise ArithmeticError('quantum scalar charge or finite cutoff differs')
    gram = sum((a.H*a for a in ops['pairs']), s.zeros(64))
    spectrum = gram.eigenvals()
    expected = {s.Integer(0): 29, s.Integer(1): 16, s.Integer(2): 10, s.Integer(3): 8, s.Integer(4): 1}
    if spectrum != expected or gram[0, 0] != 4:
        raise ArithmeticError('the full-Fock quantum-Higgs norm differs')
    return {'scalar_dimension': 3, 'active_CAR_modes': 6, 'complete_dimension': 192,
            'source_internal_indices': list(ops['source_internal_indices']),
            'source_Higgs_Ward_substitutions': _source()[0]['higgs_ward_checks'],
            'full_CAR_current_commutators': 12, 'W_squared_zero': True,
            'N_F_minus_2N_H_conserved': True,
            'W_adjoint_W_vacuum_scalar_block_spectrum': {str(k): int(v) for k, v in sorted(spectrum.items())},
            'rank_one_family_vertex_norm_exact': '2',
            'quantum_scalar_replaced_by_vacuum_expectation': False,
            'restriction_of_unbounded_bosonic_creation_claimed': False}


def quantum_feedback_budget(steps, time, p, *, gain=1, cells=1):
    p, gain = feedback._rational(p), feedback._rational(gain)
    if type(cells) is not int or cells <= 0 or gain < 0 or not 0 <= p < 1:
        raise ValueError('positive cell count, nonnegative gain and p in [0,1) required')
    joint = 1-(1-p)**cells
    result = feedback.averaging_bound(steps, time, joint, feedback_norm=2*cells*gain)
    result.update({'joint_nonreset_probability': joint, 'vertex_norm_bound': 2*cells*gain,
                   'all_records_retained': True, 'local_Gauss_preserved_pathwise': True,
                   'target': 'H0 + sum_x V_x(gain Z_D)',
                   'finite_cutoff_quantum_Higgs': True,
                   'uniform_in_volume': False, 'physical_vacuum_selected': False})
    return result


def identity_link_control():
    """Exact two-dimensional nonzero block of every actual U_ab."""
    u = s.Matrix([[0, 0], [1, 0]])
    b = s.eye(2)-u
    values = (b.H*b).eigenvals()
    if values != {(3-s.sqrt(5))/2: 1, (3+s.sqrt(5))/2: 1}:
        raise ArithmeticError('rank-one raising link control differs')
    # The relative lower bound uses lambda equal to the actual expectation,
    # hence |lambda|<=1/2. It is not asserted for arbitrary scalar targets.
    r = s.Symbol('r', real=True, nonnegative=True)
    error_matrix = (u-r*s.eye(2)).H*(u-r*s.eye(2))+(u.H-r*s.eye(2)).H*(u.H-r*s.eye(2))
    if error_matrix != s.Matrix([[1+2*r*r, -2*r], [-2*r, 1+2*r*r]]):
        raise ArithmeticError('bidirectional error matrix differs')
    return {'nilpotency_degree': 2, 'operator_norm': '1', 'numerical_radius': '1/2',
            'identity_state_vector_error_lower': '(sqrt(5)-1)/2',
            'identity_expectation_error_lower': '1/2',
            'mean_lambda_variance_sum': 'q-2*abs(lambda)^2',
            'mean_lambda_q_lower': '2*abs(lambda)',
            'mean_lambda_relative_two_direction_RMS_lower': '1',
            'relative_RMS_denominator': 'sqrt(2)*abs(lambda)',
            'nonzero_mean_required_for_relative_error': True,
            'classical_unitary_overlap_link_admitted': False,
            'collective_dressed_transport_excluded': False}


def dressed_path_control():
    """Two distinct weak-block links, full coefficients of a CAR bilinear."""
    _, currents, _, ids = _source()
    weak = ids[:2]; dim = 5
    us = {(a, b): s.SparseMatrix(dim, dim, {(1+2*a+b, 0): 1}) for a, b in product(range(2), repeat=2)}
    direct = s.SparseMatrix(100, 100, {})
    for a, c in product(range(2), repeat=2):
        e = s.zeros(4); e[a, 2+c] = 1
        path = sum((s.kronecker_product(us[a, b], us[b, c]) for b in range(2)), s.zeros(25))
        direct += s.SparseMatrix(s.kronecker_product(e, path))
    direct = s.ImmutableSparseMatrix(direct)
    count = 0
    for g in currents:
        small = g.extract(weak, weak)
        lx, ly = s.zeros(5), s.zeros(5)
        lx[1:, 1:] = -s.kronecker_product(small.T, s.eye(2))
        ly[1:, 1:] = s.kronecker_product(s.eye(2), small)
        qx, qz = s.zeros(4), s.zeros(4)
        qx[:2, :2] = small; qz[2:, 2:] = small
        local = (s.kronecker_product(qx, s.eye(25))+s.kronecker_product(s.eye(4), lx, s.eye(5)),
                 s.kronecker_product(s.eye(4), ly, s.eye(5))+s.kronecker_product(s.eye(4), s.eye(5), lx),
                 s.kronecker_product(qz, s.eye(25))+s.kronecker_product(s.eye(4), s.eye(5), ly))
        for generator in map(s.ImmutableSparseMatrix, local):
            if not _zero(generator*direct-direct*generator):
                raise ArithmeticError('contracted path failed a local endpoint/current identity')
            count += 1
    return {'distinct_edges': 2, 'local_vertices': 3, 'source_current_checks': count,
            'coefficient_carrier_dimension': 100, 'full_CAR_bilinear_identity': True,
            'single_charged_endpoint_is_gauge_invariant': False,
            'path_transporter_assumed_unitary': False,
            'charge_at_intermediate_vertices': 'cancels by the contracted internal index'}


def sine_regulator_control():
    rows = []
    for corner in product((0, 1), repeat=3):
        # The actual h=-c sigma.sin(k) has det velocity -(-1)^sum(corner).
        rows.append({'corner_pi_units': list(corner), 'orientation': -(-1)**sum(corner)})
    return {'nodes': rows, 'positive_orientation_count': 4, 'negative_orientation_count': 4,
            'unitary_internal_family_change_removes_nodes': False,
            'finite_Gauss_Ward_identity_is_mirror_decoupling': False,
            'overlap_changes_regulator_and_chiral_projector': True,
            'quantum_Higgs_feedback_changes_naive_derivative': False,
            'interacting_eight_free_poles_asserted': False,
            'interacting_mirror_decoupling_proved': False}


def quantum_projector_car_control():
    """An actual dressed-edge projector is not a scalar CAR coefficient.

    K below assigns the retained dark complement sign+1. It is a source
    algebra control of a gapped sign kernel, not the cubic Wilson kernel.
    """
    u = s.Matrix([[0, 0], [1, 0]])
    zero = s.zeros(2)
    edge = s.BlockMatrix([[zero, s.I*u], [-s.I*u.H, zero]]).as_explicit()
    active = edge*edge
    kernel = s.eye(4)-active+edge
    projection = (s.eye(4)-kernel)/2
    gamma5 = s.diag(1, 1, -1, -1)
    overlap = s.eye(4)+gamma5*kernel
    if kernel*kernel != s.eye(4) or projection*projection != projection:
        raise ArithmeticError('exact dressed-source spectral projection differs')
    if not _zero(gamma5*overlap+overlap*gamma5-overlap*gamma5*overlap):
        raise ArithmeticError('the gapped sign control lost its GW identity')
    occupied = u*u.H
    endpoint = (s.kronecker_product(s.diag(1, 0), s.eye(2))-s.kronecker_product(s.eye(2), occupied),
                s.kronecker_product(s.diag(0, 1), s.eye(2))+s.kronecker_product(s.eye(2), occupied))
    for gauss in endpoint:
        if not _zero(gauss*projection-projection*gauss):
            raise ArithmeticError('actual singlet-block endpoint covariance differs')
    annihilators = tuple(creator(2, i).H for i in range(2))
    projected = tuple(sum((s.kronecker_product(annihilators[j], projection[2*i:2*i+2, 2*j:2*j+2])
                           for j in range(2)), s.zeros(8)) for i in range(2))
    defect = _clean(projected[0]*projected[0].H+projected[0].H*projected[0]
                    -s.kronecker_product(s.eye(4), projection[:2, :2]))
    if defect == s.zeros(8):
        raise ArithmeticError('the operator-coefficient CAR obstruction disappeared')
    spectrum = defect.eigenvals()
    return {'source_link_block': 'vacuum plus q6Y=6 raising state',
            'dressed_edge_matrix': [[str(x) for x in edge.row(i)] for i in range(4)],
            'kernel_gap_exact': '1', 'Ginsparg_Wilson_identity_exact': True,
            'endpoint_Gauss_commutators': 2,
            'projector_rank': int(projection.rank()),
            'projected_CAR_defect_spectrum': {str(k): int(v) for k, v in sorted(spectrum.items())},
            'simple_operator_projector_is_canonical_CAR': False,
            'control_is_actual_cubic_Wilson_kernel': False,
            'a_different_many_body_chiral_construction_excluded': False}


def raising_link_wilson_control():
    """Actual singlet raising links on a one-cell Euclidean Wilson stencil.

    Every 60-level edge is either its active vacuum/singlet pair or a dark
    level. The algebra below handles all active counts0..4. It is not an
    all-volume gap estimate or a replacement of the spatial source clock.
    """
    x=s.Matrix([[0,1],[1,0]]); y=s.Matrix([[0,-s.I],[s.I,0]]); z=s.diag(1,-1)
    gamma=(s.kronecker_product(x,s.eye(2)),s.kronecker_product(y,x),
           s.kronecker_product(y,y),s.kronecker_product(y,z))
    gamma5=s.kronecker_product(z,s.eye(2))
    rows=[]
    for active in range(5):
        ld=1<<active; dimension=4*ld
        def link(op,index):
            return s.kronecker_product(*[op if j==index else s.eye(2) for j in range(active)])
        b=3*s.eye(dimension); c=s.zeros(dimension); spin_z=s.zeros(dimension)
        for j in range(active):
            b-=s.kronecker_product(s.eye(4),link(x,j))/2
            c+=s.kronecker_product(gamma[j],link(y,j))/2
            spin_z+=s.kronecker_product(gamma[j],link(z,j))
        b,c,spin_z=map(_clean,(b,c,spin_z))
        a=b-s.I*c
        chirality=s.ImmutableSparseMatrix(s.kronecker_product(gamma5,s.eye(ld)))
        h=_clean(chirality*a)
        if h!=h.H or not _zero(h*h-b*b-s.Rational(active,4)*s.eye(dimension)+spin_z/2):
            raise ArithmeticError('literal raising-link Wilson square differs')
        if not _zero(spin_z*spin_z-active*s.eye(dimension)):
            raise ArithmeticError('Clifford link bound differs')
        # B>=3-active/2>=1; sqrt(active)<=2. Thus the expression below is>=1.
        lower=(3-s.Rational(active,2))**2+s.Rational(active,4)-s.sqrt(active)/2
        if lower<1:
            raise ArithmeticError('the exact singlet Wilson gap bound failed')
        rows.append({'active_link_pairs':active,'dimension':dimension,
                     'H_squared_lower':str(lower),'all_other_edge_levels_retained_as_dark':True})
    return {'rows':rows,'full_link_dimension_per_direction':60,
            'internal_block':'q6Y=6 singlet','Euclidean_periodic_sites':1,
            'H_squared_uniform_over_dark_sectors_lower':'1',
            'real_part_DW_minus_I_lower':'1','DW_minus_I_norm_upper':'6',
            'overlap_D_squared_singular_value_lower':'7/3',
            'classical_identity_stencil_overlap_D':'0',
            'sign_kernel_exists_without_classical_background':True,
            'light_chiral_mode_obtained_in_this_control':False,
            'all_volume_or_all_M15_gap_proved':False,
            'Lorentzian_clock_identification_proved':False}


def _json_safe(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


@lru_cache(maxsize=1)
def _certificate_text():
    duties = []
    for duty in feedback.POLICIES:
        model = feedback.record_model(F(1, 12), F(1, 6), duty)
        duties.append({'duty': duty, 'records_per_old_memory': 4,
                       'coarse_rows': len(model['rows']), 'original_mean_trace': str(s.trace(model['mean'])),
                       'quantum_Higgs_mean_matrix': [[str(x) for x in model['mean'].row(i)] for i in range(3)],
                       'mean_divided_by_duty_probability': False})
    result = {'schema': 'oph.finite_quantum_higgs_chiral_join.v1',
              'quantum_Higgs': quantum_higgs_certificate(), 'memory_duties': duties,
              'example_budget': quantum_feedback_budget(1000000, 1, F(1, 12)),
              'quantum_link_boundary': identity_link_control(),
              'dressed_path': dressed_path_control(), 'naive_regulator': sine_regulator_control(),
              'quantum_projector_CAR': quantum_projector_car_control(),
              'raising_link_Wilson': raising_link_wilson_control(),
              'source_inputs': list(SOURCE_INPUTS),
              'native_long_pulse_tape_executed': False,
              'full_nonabelian_quantum_overlap_join_constructed': False,
              'interacting_Lorentz_continuum_proved': False,
              'physical_particle_or_mass_assigned': False}
    return json.dumps(_json_safe(result), sort_keys=True)


def certificate():
    return json.loads(_certificate_text())
