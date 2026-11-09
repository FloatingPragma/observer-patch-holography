"""Same-current free Weyl frames and exact mirror/antiparticle distinctions."""
from fractions import Fraction as F
from itertools import product

import sympy as s

from model.common_source import native_matter_attachment as matter
from model.common_source.chiral_candidate import action
from model.weyl_geometric_modular import algebra as geometry


def _zero(matrix):
    return all(s.expand(x) == 0 for x in matrix)


def _internal():
    data = matter.representation()
    p, _ = matter.action_intertwiner()
    currents = data['matter_currents']
    old = action.generators(families=False)
    ports = data['native']['native_to_fixture_port']
    if p.T*p != s.eye(15) or len(currents) != 12:
        raise ArithmeticError('the existing source current dictionary differs')
    for index, g in enumerate(currents):
        if not _zero(old[ports[index]]*p-p*g):
            raise ArithmeticError('the actual source current permutation differs')
    return data, s.ImmutableSparseMatrix(p), currents, tuple(old[ports[i]] for i in range(12))


def valley_frame(valley):
    """Map native internal/source spin coordinates to canonical old-action ones."""
    _, p, _, _ = _internal()
    signs, chi, spin = geometry.frame(valley)
    change = s.ImmutableSparseMatrix(s.kronecker_product(spin.H, p, s.eye(3)))
    return {'map': change, 'chirality': chi, 'source_axis_signs': tuple(-x for x in signs),
            'direction': 'native source spin/internal/family to canonical Weyl/action/family'}


def current_frames():
    """Exact ninety-row same-current identities at every source Weyl corner."""
    data, p, native, old = _internal()
    pauli = geometry.pauli()
    identities = 0
    rows = []
    for valley in geometry.VALLEYS:
        signs, chi, spin = geometry.frame(valley)
        change = s.ImmutableSparseMatrix(s.kronecker_product(spin.H, p, s.eye(3)))
        if change.H*change != s.eye(90):
            raise ArithmeticError('the source frame does not preserve the CAR metric')
        for axis in range(3):
            source_axis = s.ImmutableSparseMatrix(s.kronecker_product(-signs[axis]*pauli[axis], s.eye(45)))
            canonical_axis = s.ImmutableSparseMatrix(s.kronecker_product(chi*pauli[axis], s.eye(45)))
            if not _zero(canonical_axis*change-change*source_axis):
                raise ArithmeticError('the source kinetic axis lost its current frame')
        for native_g, old_g in zip(native, old, strict=True):
            ng = s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), native_g, s.eye(3)))
            og = s.ImmutableSparseMatrix(s.kronecker_product(s.eye(2), old_g, s.eye(3)))
            if not _zero(og*change-change*ng):
                raise ArithmeticError('a gauge current changes under the Weyl frame')
            identities += 1
        rows.append({'node_pi_units': list(valley), 'chirality': chi,
                     'physical_Cauchy_rows_per_node': 90, 'internal_family_components': 45})
    return {'rows': rows, 'same_current_intertwining_identities': identities,
            'principal_symbol_axis_identities': 24,
            'field_order': 'spin*45+3*matter_index+family',
            'native_matter_masks': list(data['selected']),
            'inherited_configured_coordinate_speed': str(action.source_coefficients()['Omega']/2),
            'kinetic_or_Yukawa_coefficients_refitted': False,
            'interacting_vacuum_or_continuum_attachment_inferred': False}


def lorentz_currents():
    """Free spinor current identities; a pointwise boost is not a CAR unitary."""
    pauli = geometry.pauli()
    epsilon = s.Matrix([[0, 1], [-1, 0]])
    checked = 0
    # A rational nontrivial boost: cosh(u)=5/3, sinh(u)=4/3.
    ch, sh = s.Rational(5, 3), s.Rational(4, 3)
    for chi in (-1, 1):
        current = (s.eye(2), *(chi*x for x in pauli))
        rotations = tuple(-s.I*x/2 for x in pauli)
        boosts = tuple(chi*x/2 for x in pauli)
        for i, j in product(range(3), repeat=2):
            bracket = lambda x, y: x*y-y*x
            r = sum((s.LeviCivita(i, j, k)*rotations[k] for k in range(3)), s.zeros(2))
            b = sum((s.LeviCivita(i, j, k)*boosts[k] for k in range(3)), s.zeros(2))
            if bracket(rotations[i], rotations[j]) != r or bracket(rotations[i], boosts[j]) != b or bracket(boosts[i], boosts[j]) != -r:
                raise ArithmeticError('Lorentz spin generators differ')
        for generator in (*rotations, *boosts):
            if generator.T*epsilon+epsilon*generator != s.zeros(2):
                raise ArithmeticError('the Weyl antisymmetric invariant differs')
        for axis in range(3):
            transform = ch*s.eye(2)+chi*sh*pauli[axis]
            lorentz = s.eye(4)
            lorentz[0, 0] = lorentz[axis+1, axis+1] = ch*ch+sh*sh
            lorentz[0, axis+1] = lorentz[axis+1, 0] = 2*ch*sh
            metric = s.diag(1, -1, -1, -1)
            if transform.det() != 1 or lorentz.T*metric*lorentz != metric:
                raise ArithmeticError('the finite boost lost its Lorentz norm')
            for mu in range(4):
                target = sum((lorentz[mu, nu]*current[nu] for nu in range(4)), s.zeros(2))
                if not _zero(transform.H*current[mu]*transform-target):
                    raise ArithmeticError('the same-current Lorentz vector identity differs')
                checked += 1
    return {'finite_spinor_current_identities': checked,
            'Hermitian_internal_current': 'i*T for the actual anti-Hermitian T',
            'vector_current': 'psi^dagger ((I,chi*sigma) tensor i*T) psi',
            'time_coordinate': 'x0=v*t with the configured coordinate speed v',
            'pointwise_boost_is_unitary_on_equal_time_CAR': False,
            'finite_momentum_menu_closed_under_all_Lorentz_boosts': False,
            'free_full_field_representation_required_for_unitary_boost': True}


def mirror_invariant():
    """A fifth charge moment separates M15 from its conjugate, without anomaly claims."""
    data = matter.representation()
    selected = list(data['selected'])
    y = data['operators']['q'].extract(selected, selected)/6
    conjugate = -y.T
    weights = {}
    for value in y.diagonal():
        weights[str(value)] = weights.get(str(value), 0)+1
    traces = {str(power): str(s.trace(y**power)) for power in (1, 3, 5)}
    if traces != {'1': '0', '3': '0', '5': '5/9'} or s.trace(conjugate**5) != -s.Rational(5, 9):
        raise ArithmeticError('the exact selected hypercharge spectrum differs')
    return {'normalization': 'Y=q6Y/6, electric charge Q=T3+Y',
            'one_generation_internal_dimension': 15, 'weights': weights,
            'hypercharge_power_traces': traces, 'conjugate_fifth_trace': '-5/9',
            'M15_equivalent_to_conjugate_preserving_hypercharge': False,
            'fifth_trace_is_a_four_dimensional_gauge_anomaly': False,
            'spatial_spinor_rows_counted_as_additional_multiplets': False,
            'Nambu_holes_counted_as_additional_multiplets': False}


def local_valley_weights(side=4):
    """Original complete local addition/removal law of the free sine stencil."""
    if type(side) is not int or side < 4 or side % 4:
        raise ValueError('periodic side divisible by four and at least four required')
    volume = side**3
    rows = []
    for valley in geometry.VALLEYS:
        chi = geometry.frame(valley)[1]
        rows.append({'node_pi_units': list(valley), 'chirality': chi,
                     'original_combined_weight': F(1, 8),
                     'original_strictly_positive_energy_weight': F(1, 8)-F(1, volume)})
    return {'side': side, 'spatial_sites': volume, 'rows': rows,
            'local_instrument': 'a_(x,s)/sqrt(2), a_(x,s)^dagger/sqrt(2), s=0,1',
            'local_records_before_energy_and_momentum': 4,
            'complete_charge_spin_momentum_records': 4*volume,
            'one_chirality_combined_weight': F(1, 2),
            'one_chirality_strictly_positive_energy_weight': F(1, 2)-F(4, volume),
            'success_probability_divided_out': False,
            'uses_free_filled_band_reference': True,
            'quantum_link_vacuum_replaced_by_identity': False,
            'additional_native_tape_executed': False}


def local_gauge_obstruction():
    """Valley restriction is not reducing for an admitted local gauge current."""
    data = matter.representation(); indices = list(data['selected'])
    y = data['operators']['q'].extract(indices, indices)/6
    trace2 = s.trace(y*y)
    if trace2 != s.Rational(10, 3):
        raise ArithmeticError('source hypercharge trace normalization differs')
    # A valley Fourier projector has every spatial diagonal entry 1/8.
    scalar = 2*s.Rational(1, 8)*s.Rational(7, 8)
    return {'one_spin_one_generation_commutator_HS_squared': str(scalar*trace2),
            'full_two_spin_three_family_commutator_HS_squared': str(6*scalar*trace2),
            'generator': 'site projector tensor i*Y',
            'projector': 'one pi-shift valley cell tensor I15',
            'valley_projection_commutes_with_all_local_source_gauges': False,
            'failure_complement_may_be_deleted': False}


def certificate():
    value = {'same_current_valley_frames': current_frames(),
            'free_Lorentz_current_algebra': lorentz_currents(),
            'mirror_versus_hole_charge_invariant': mirror_invariant(),
            'complete_free_local_valley_weights': local_valley_weights(),
            'local_gauge_valley_obstruction': local_gauge_obstruction(),
            'configured_finite_menu': {'left_internal_family_fields': 45,
                'Cauchy_rows': 90, 'Nambu_rows': 180,
                'Nambu_doubling_adds_independent_fields': False},
            'local_sine_regulator': {'Weyl_fields_counting_internal_and_family': 360,
                'each_chirality_counting_internal_and_family': 180,
                'Dirac_pairings_counting_internal_and_family': 180,
                'mirror_removal_or_weight_suppression_proved': False},
            'scope': {'source_current_and_free_kinetic_intertwiner_exact': True,
                'same_dynamical_gauge_vacuum_attached': False,
                'interacting_chiral_Lorentz_limit_proved': False,
                'physical_particle_or_mass_assignment': False}}
    def serialize(item):
        if isinstance(item, F):
            return str(item)
        if isinstance(item, dict):
            return {key: serialize(child) for key, child in item.items()}
        if isinstance(item, (list, tuple)):
            return [serialize(child) for child in item]
        return item
    return serialize(value)
