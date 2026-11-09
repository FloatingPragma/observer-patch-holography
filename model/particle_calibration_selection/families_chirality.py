"""Native event cost and a canonical full-Fock quantum-link dressing.

The carrier cost is computed from the existing twelve-rotation service.
The dressing uses the existing singlet component of the complete 60-level
link. Neither construction selects a physical generation count, an
interacting chiral continuum, a vacuum, or a laboratory mass.
"""
from __future__ import annotations

from functools import lru_cache
import json

import sympy as s

from model.common_source import native_family_selection as family
from model.particle_candidate_atlas.interacting_chiral import creator
from model.standard_model_source_map import gauge_attachment as gauge


SOURCE_INPUTS = (
    'model/common_source/native_family_selection.py',
    'model/common_source/native_gauge_object.py',
    'model/standard_model_source_map/gauge_attachment.py',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/common_source/overlap_regulator.py',
    'model/interacting_extension/chiral_electric/PROOF.md',
    'model/interacting_extension/chiral_lattice/PROOF.md',
    'code/a5_closure/family_band_attachment_certificate.py',
    'Lean/Screen/A5FamilyBand.lean',
)


def _clean(value):
    value = s.ImmutableSparseMatrix(value)
    return s.ImmutableSparseMatrix(value.rows, value.cols,
        {key: v for key, x in value.todok().items()
         if (v := s.simplify(s.expand(x))) != 0})


def _zero(value):
    return not _clean(value).todok()


@lru_cache(maxsize=1)
def event_cost_data():
    """Mean squared displacement per near-rotation event, on all36 states."""
    data = family.geometry()
    identity = s.eye(36)
    rotations = tuple(family.objects.permutation_matrix(data['pair_actions'][j])
                      for j in data['near_indices'])
    direct = _clean(sum(((r-identity).H*(r-identity) for r in rotations),
                        s.zeros(36))/len(rotations))
    cost = _clean(2*identity-data['C']/6)
    if direct != cost:
        raise ArithmeticError('the complete near-rotation event cost differs')
    coefficient = _clean(data['W'].H*cost*data['W'])
    adjacency = (data['C12']-2*s.eye(12))/2
    if not _zero(coefficient-(5*s.eye(12)-adjacency)/3):
        raise ArithmeticError('the operational cost is not the carrier Laplacian/3')
    spectrum = {s.simplify(2-root/6): multiplicity
                for root, multiplicity in data['roots'].items()}
    positive = sorted((x for x in spectrum if x != 0), key=lambda x: float(x))
    least = positive[0]
    if not all((x-least).is_positive is True for x in positive[1:]):
        raise ArithmeticError('the smallest nonzero event cost is not strict')
    if least != (5-s.sqrt(5))/3 or spectrum[least] != 3 or spectrum[0] != 2:
        raise ArithmeticError('native full-carrier multiplicities changed')
    if not _zero(cost*data['P']-least*data['P']):
        raise ArithmeticError('the event-cost band differs from the retained cooling band')
    all_rotations = tuple(family.objects.permutation_matrix(p) for p in data['pair_actions'])
    all_cost = _clean(2*identity-sum((r+r.H for r in all_rotations), s.zeros(36))/60)
    if not _zero(all_cost-2*(identity-data['fixed'])):
        raise ArithmeticError('uniform whole-group cost lost its complete non-invariant degeneracy')
    return {'cost': cost, 'coefficient_cost': coefficient,
            'spectrum': spectrum, 'selected': least, 'projector': data['P'],
            'invariant': data['fixed'], 'rotations': rotations,
            'all_rotation_cost': all_cost,
            'alternative_projector': data['projectors'][s.Integer(0)]}


def family_certificate():
    data = event_cost_data()
    rank = int(s.trace(data['projector']))
    alternate = int(s.trace(data['alternative_projector']))
    # The corrected C/E branches are scalar multiples of identity on this
    # protected object. Their complete unread map fixes the maximally mixed
    # state. The added cooling service instead maps it to P/rank(P).
    full_mixed = s.eye(36)/36
    p = data['projector']
    cooled = _clean(p*full_mixed*p+s.trace((s.eye(36)-p)*full_mixed)*p/rank)
    if not _zero(cooled-p/rank) or _zero(cooled-full_mixed):
        raise ArithmeticError('passive repair and added cooling were conflated')
    return {
        'carrier_dimension': 36,
        'near_rotation_event_count': len(data['rotations']),
        'cost_formula': 'E_g ||(R_g-I)v||^2 = <v,(2I-C/6)v>',
        'coefficient_cost': '(5I-A)/3',
        'cost_spectrum': [{'cost': str(x), 'rank': int(n)}
                          for x, n in sorted(data['spectrum'].items(), key=lambda r: float(r[0]))],
        'smallest_positive_cost': str(data['selected']),
        'smallest_positive_cost_rank': rank,
        'dimension_window_used_in_this_derivation': False,
        'invariant_zero_cost_rank': int(s.trace(data['invariant'])),
        'same_existing_native_cooling_projector': True,
        'uniform_all60_rotation_positive_cost_rank': 34,
        'near_rotation_rule_selected_by_symmetry_alone': False,
        'passive_corrected_C_E_selects_the_cooling_service': False,
        'alternative_covariant_cooling_rank': alternate,
        'all36_carrier_levels_and_complete_records_retained': True,
        'three_physical_generations_selected': False,
        'missing_source_operator': 'a selected complete family controller and a current-preserving matter attachment to its output algebra',
        'physical_family_attachment_is_the_cost_rank_identity': False,
    }


def replica_trace_control(copies):
    """Independent finite witness of external-copy blindness.

    General proof: Tr((rho tensor I_n/n)(A tensor I_n))=Tr(rho A).
    Every product of source operators has the same tensor factorization.
    """
    if type(copies) is not int or not 1 <= copies <= 8:
        raise ValueError('an integer copy count from one through eight required')
    rho = s.Matrix([[s.Rational(2, 3), s.I/6], [-s.I/6, s.Rational(1, 3)]])
    x = s.Matrix([[0, 1], [1, 0]])
    z = s.diag(1, -1)
    word = (s.eye(2)+x/3)*(s.eye(2)+z/5)*(s.eye(2)-x/7)
    return s.simplify(s.trace(s.kronecker_product(rho, s.eye(copies)/copies)
                             *s.kronecker_product(word, s.eye(copies))))


@lru_cache(maxsize=1)
def dressed_data():
    """Canonical dressing on two CAR modes and every native link level.

    The two modes carry the retained charged singlet at opposite endpoints.
    The link components and all twelve endpoint currents are read from the
    owning quantum-link implementation. This is a finite local operator
    construction, not a chosen chiral lattice action.
    """
    labels = gauge.link_labels()
    singleton = [a for a, b in labels[1:] if a == b
                 and sum(left == a for left, _ in labels[1:]) == 1]
    if len(labels) != 60 or len(singleton) != 1:
        raise ArithmeticError('the unique native charged singlet link changed')
    internal = singleton[0]
    active = labels.index((internal, internal))
    link = s.ImmutableSparseMatrix(60, 60, {(active, 0): 1})
    a = tuple(creator(2, j).H for j in range(2))
    number = tuple(_clean(op.H*op) for op in a)
    identity = s.eye(240)
    forward = s.kronecker_product(a[0].H*a[1], link)
    hopping = _clean(forward+forward.H)
    if not _zero(hopping**3-hopping):
        raise ArithmeticError('the native active link is not an exact two-level exchange')
    # exp(i*pi*S/4) on the complete domain, including vacuum, double
    # occupation and every dark link state. No projected linear CAR ansatz.
    unitary = _clean(identity+(1/s.sqrt(2)-1)*hopping**2+s.I*hopping/s.sqrt(2))
    if not _zero(unitary.H*unitary-identity):
        raise ArithmeticError('the complete dressing is not unitary')
    bare = tuple(s.kronecker_product(op, s.eye(60)) for op in a)
    dressed = tuple(_clean(unitary*op*unitary.H) for op in bare)
    for i in range(2):
        for j in range(2):
            if not _zero(dressed[i]*dressed[j]+dressed[j]*dressed[i]):
                raise ArithmeticError('dressed annihilators violate CAR')
            if not _zero(dressed[i]*dressed[j].H+dressed[j].H*dressed[i]
                         -(identity if i == j else s.zeros(240))):
                raise ArithmeticError('dressed creation/annihilation CAR differs')
    gauss = []
    for generator in gauge.currents():
        if any(generator[j, internal] or generator[internal, j]
               for j in range(15) if j != internal):
            raise ArithmeticError('the selected internal line is not a current singlet')
        for endpoint in (0, 1):
            current = _clean(s.kronecker_product(generator[internal, internal]*number[endpoint], s.eye(60))
                             +s.kronecker_product(s.eye(4), gauge.link_generator(generator, endpoint)))
            if not _zero(current*unitary-unitary*current):
                raise ArithmeticError('full-Fock unitary violates an endpoint Gauss identity')
            gauss.append(current)
    total_number = s.kronecker_product(number[0]+number[1], s.eye(60))
    if not _zero(unitary*total_number-total_number*unitary):
        raise ArithmeticError('dressing changed total occupation')
    one_indices = [60*mask+j for mask in (1, 2) for j in range(60)]
    one_unitary = unitary.extract(one_indices, one_indices)
    p0 = s.diag(*([1]*60+[0]*60))
    p = _clean(one_unitary*p0*one_unitary.H)
    edge = s.BlockMatrix([[s.zeros(60), s.I*link], [-s.I*link.H, s.zeros(60)]]).as_explicit()
    old_p = _clean((edge**2-edge)/2)
    dark = _clean(p-old_p)
    if s.trace(old_p) != 1 or s.trace(p) != 60 or s.trace(dark) != 59 or not _zero(dark*dark-dark):
        raise ArithmeticError('the exact rank obstruction or retained dark completion changed')
    kernel = _clean(s.eye(120)-2*p)
    gamma5 = p0*2-s.eye(120)
    overlap = _clean(s.eye(120)+gamma5*kernel)
    if not _zero(kernel*kernel-s.eye(120)) or not _zero(gamma5*overlap+overlap*gamma5-overlap*gamma5*overlap):
        raise ArithmeticError('the balanced sign control lost its GW identity')
    heavy = _clean(dressed[1].H*dressed[1])
    h0 = _clean(identity-hopping)
    return {'unitary': unitary, 'hopping': hopping, 'bare': bare,
            'dressed': dressed, 'gauss': tuple(gauss), 'total_number': total_number,
            'one_particle_projector': p, 'old_projector': old_p,
            'dark_completion': dark, 'kernel': kernel, 'overlap': overlap,
            'heavy': heavy, 'parent': h0, 'active_link_index': active,
            'internal_index': internal, 'one_indices': one_indices}


def penalty_certificate(strength):
    """Exact finite heavy-complement bound, with its coefficient supplied."""
    strength = s.sympify(strength)
    if strength.is_Rational is not True or strength <= 1:
        raise ValueError('an exact rational penalty greater than one is required')
    # The only mixed block of H0 is a Pauli X. A common unitary preserves its
    # spectrum. All other occupation/link states have energy1 or1+strength.
    low = s.simplify(1+strength/2-s.sqrt(strength**2/4+1))
    high = s.simplify(1+strength/2+s.sqrt(strength**2/4+1))
    error = s.simplify(1-low)
    bound = 1/(strength-1)
    if not (bound-error).is_positive:
        raise ArithmeticError('the finite Schur-complement bound failed')
    return {'supplied_strength': str(strength), 'full_dimension': 240,
            'low_subspace_dimension': 120, 'retained_complement_dimension': 120,
            'mixed_low_energy': str(low), 'mixed_high_energy': str(high),
            'other_low_energy': '1', 'other_low_multiplicity': 119,
            'other_high_energy': str(1+strength), 'other_high_multiplicity': 119,
            'maximum_low_energy_error': str(error), 'schur_bound': str(bound),
            'bound_is_physical_mirror_decoupling': False}


def chirality_certificate():
    data = dressed_data()
    # The canonical completion changes dark-kernel signs. It is not a
    # finite-volume relativistic kinetic operator:118 dark one-particle
    # columns of this sign control have D=0.
    dark_zero = sum(all(data['overlap'][j, i] == 0 for j in range(120))
                    for i in range(120))
    return {
        'source_link_dimension': 60, 'full_Fock_dimension': 240,
        'fermion_modes': 2, 'retained_occupation_sectors': [0, 1, 2],
        'native_singlet_internal_index': data['internal_index'],
        'native_singlet_link_index': data['active_link_index'],
        'active_matter_subrepresentation': 'the charged e^c-type internal singlet',
        'full_nonabelian_matter_regulator_constructed': False,
        'endpoint_Gauss_identities_exact': len(data['gauss']),
        'all_full_Fock_CAR_identities_exact': True,
        'dressed_operators': 'b_i=W c_i W^dagger; W=exp(i*pi*S/4)',
        'source_operator': 'S=c_x^dagger U c_y+c_y^dagger U^dagger c_x',
        'total_occupation_preserved': True,
        'linear_operator_projector_replaced_by_many_body_unitary': True,
        'old_one_particle_projector_rank': int(s.trace(data['old_projector'])),
        'canonical_one_particle_projector_rank': int(s.trace(data['one_particle_projector'])),
        'retained_dark_completion_rank': int(s.trace(data['dark_completion'])),
        'old_rank_one_projector_is_one_mode_tensor_full_link': False,
        'balanced_sign_Ginsparg_Wilson_identity_exact': True,
        'endpoint_grading_identified_with_spacetime_chirality': False,
        'balanced_sign_control_dark_zero_columns': dark_zero,
        'balanced_sign_control_is_a_physical_kinetic_regulator': False,
        'complete_link_register_replaced_by_classical_unitaries': False,
        'source_pulse_tape_executed_for_this_dressing': False,
        'native_state_or_physical_clock_selected': False,
        'interacting_mirror_decoupling_or_Lorentz_limit_proved': False,
        'penalty_examples': [penalty_certificate(n) for n in (2, 4, 8)],
        'minimal_source_data': [
            'a selected gauge-covariant local projector path with compatible charge-sector ranks',
            'its many-body gauge-invariant transport and induced electric connection',
            'a source-selected kinetic/complement Hamiltonian with uniform gap and locality control',
            'anomaly-compatible chiral measure and a positive physical state on the same clock',
        ],
    }


def certificate():
    return {'schema': 'oph.particle_calibration.families_chirality.v1',
            'source_inputs': list(SOURCE_INPUTS),
            'empirical_mass_or_mixing_inputs': [],
            'family_cost': family_certificate(),
            'canonical_quantum_link_dressing': chirality_certificate(),
            'external_copy_trace_witness': {str(n): str(replica_trace_control(n)) for n in (1, 2, 3, 4)},
            'physical_particles_selected': [], 'physical_mass_predictions': []}


if __name__ == '__main__':
    print(json.dumps(certificate(), ensure_ascii=False, indent=2, sort_keys=True))
