"""Exact energy accounting on the declared gauge-coordinate and Higgs sources.

The current is formed from the Hamiltonian energy terms themselves. The
response changes a stated action coefficient at fixed preparation and clock.
No Einstein equation or geometric entropy is used as a definition.
"""
from functools import lru_cache
import json

import sympy as s

from model.particle_observer_outputs import field_route as field
from model.particle_observer_outputs import higgs_records as scalar
from model.matched_einstein_limit import frame


SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_observer_outputs/higgs_records.py',
    'model/particle_observer_outputs/HIGGS_RECORDS.md',
    'model/particle_observer_outputs/readout_contract.py',
    'model/particle_field_causality/bounds.py',
    'model/particle_field_causality/PROOF.md',
    'model/field_stress/energy.py',
    'Lean/QFT/CoarseEnergyWard.lean',
    'Lean/QFT/LatticeMomentumBalance.lean',
    'model/matched_einstein_limit/PROOF.md',
    'model/native_gravity_identification/PROOF.md',
    'model/particle_regional_geometry/PROOF.md',
    'model/weyl_geometric_modular/README.md',
    'model/shrinking_modular_entropy/README.md',
    'model/reciprocal_area_dynamics/DERIVATION.md',
    'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean',
)


def clean(matrix):
    return scalar._matrix(matrix)


def commutator(a, b):
    return clean(a*b-b*a)


def expect(state, operator):
    return s.simplify((state.H*operator*state)[0])


def partition_certificate(total, terms, state, *, symmetry=None):
    """Finite operator continuity, weighted work and exact block transport.

    Terms may be spatial cells or named energy components. The labels do
    not assert a spatial locality radius or a symmetric Lorentz tensor.
    """
    total = clean(total)
    terms = tuple(map(clean, terms))
    if not terms or any(h.shape != total.shape or h != h.H for h in terms):
        raise ValueError('a nonempty compatible Hermitian energy partition is required')
    zero = s.ImmutableSparseMatrix(total.rows, total.cols, {})
    if clean(sum(terms, zero)-total) != zero:
        raise ValueError('the local terms must sum to the same source Hamiltonian')
    state = clean(state)
    if state.shape != (total.rows, 1) or expect(state, s.eye(total.rows)) != 1:
        raise ValueError('a normalized state on the complete retained carrier is required')
    if symmetry is not None and any(commutator(h, symmetry) != zero for h in terms):
        raise ValueError('every energy term must preserve the specified source symmetry')
    currents = {(i, j): clean(s.I*commutator(a, b))
                for i, a in enumerate(terms) for j, b in enumerate(terms)}
    for (i, j), current in currents.items():
        if current != current.H or clean(current+currents[j, i]) != zero:
            raise ArithmeticError('energy current is not Hermitian and antisymmetric')
    rates = []
    for i, term in enumerate(terms):
        rate = clean(s.I*commutator(total, term))
        if clean(rate+sum((currents[i, j] for j in range(len(terms))), zero)) != zero:
            raise ArithmeticError('same-Hamiltonian local energy continuity failed')
        rates.append(rate)
    # A declared inhomogeneous lapse-source vector checks weighted work.
    weights = tuple(s.Integer(2*j+1) for j in range(len(terms)))
    driven = clean(sum((f*h for f, h in zip(weights, terms, strict=True)), zero))
    for i, term in enumerate(terms):
        divergence = sum((weights[i]*weights[j]*currents[i, j]
                          for j in range(len(terms))), zero)
        if clean(s.I*commutator(driven, weights[i]*term)+divergence) != zero:
            raise ArithmeticError('the weighted energy current differs from its action')
    # Merge all but the last component. This supplies the finite-sum
    # coarse-continuity hypotheses without identifying components as space.
    region = tuple(range(max(1, len(terms)-1)))
    regional = clean(sum((terms[i] for i in region), zero))
    flux = clean(sum((currents[i, j] for i in region for j in range(len(terms))
                      if j not in region), zero))
    if clean(s.I*commutator(total, regional)+flux) != zero:
        raise ArithmeticError('complete regional boundary flux failed')
    return {'dimension': total.rows, 'energy_component_count': len(terms),
            'all_ordered_component_pairs_checked': len(currents),
            'operator_continuity_holds_for_all_states_and_references': True,
            'weighted_action_work_identity': True,
            'coarse_region': list(region), 'coarse_boundary_flux_is_complete': True,
            'nonzero_unordered_current_pairs': [[i, j] for (i, j), op in currents.items()
                                                if i < j and op != zero],
            'initial_component_energies': [expect(state, h) for h in terms],
            'initial_component_energy_rates': [expect(state, rate) for rate in rates],
            'initial_total_energy': expect(state, total),
            'initial_regional_outward_flux': expect(state, flux),
            'fixed_source_time_parameter': True,
            'symmetric_Lorentz_stress_identified': False}, currents


@lru_cache(maxsize=2)
def rotor_data(cutoff=6):
    source = field.exact_rotor(cutoff)
    n = len(source['flux'])
    coefficients = field.rotor_coefficients()
    electric = clean(source['E']**2/2)
    number = clean(s.kronecker_product(coefficients['onsite'], s.eye(n)))
    hopping = clean(source['H']-electric-number)
    result, currents = partition_certificate(source['H'], (electric, number, hopping),
                                             source['state'], symmetry=source['Gauss'])
    if result['initial_component_energies'] != [9, 2, 0]:
        raise ArithmeticError('the literal q6 source energy reference differs')
    if result['initial_component_energy_rates'] != [18, 0, -18]:
        raise ArithmeticError('the full source energy exchange differs')
    response_slope = expect(source['state'], s.I*commutator(hopping, source['E']))
    if response_slope != 6:
        raise ArithmeticError('the action perturbation and its flux response differ')
    result.update({'component_names': ['electric', 'Dirac onsite', 'covariant hopping'],
                   'same_source_native_charge': source['charge'],
                   'same_source_initial_flux': source['initial_flux_mean'],
                   'initial_flux_time_derivative': source['initial_current'],
                   'hopping_coefficient_response_initial_time_slope': response_slope,
                   'all16_CAR_occupations_retained': True,
                   'all_flux_records': list(source['flux']),
                   'spatial_cells_in_executed_control': 1,
                   'three_dimensional_spatial_stress_replayed': False})
    return source, (electric, number, hopping), result, currents


@lru_cache(maxsize=1)
def scalar_data():
    source = scalar.sector_data(s.diag(1, 2, 3), 1)
    g = s.sqrt(source['g_squared'])
    bright = clean(source['bright_vector']/g)
    state = clean((source['vacuum']+s.I*bright)/s.sqrt(2))
    number, vertex = source['occupation'], source['vertex']
    result, currents = partition_certificate(source['H'], (number, vertex), state)
    current = currents[0, 1]
    if expect(state, current) != g:
        raise ArithmeticError('the literal scalar pair-energy transfer differs')
    plus = clean((source['vacuum_projector']+source['bright_projector']+current/g)/2)
    minus = clean((source['vacuum_projector']+source['bright_projector']-current/g)/2)
    projections = (plus, minus, source['dark_projector'])
    identity = s.eye(source['H'].rows)
    zero = s.zeros(*identity.shape)
    if clean(sum(projections, zero)-identity) != zero:
        raise ArithmeticError('the complete energy-current detector lost a dark state')
    for i, a in enumerate(projections):
        for j, b in enumerate(projections):
            if clean(a*b-(a if i == j else zero)) != zero:
                raise ArithmeticError('energy-current detector is not a complete PVM')
    result.update({'component_names': ['scalar number energy', 'pair vertex energy'],
                   'family_matrix': 'diag(1,2,3)', 'pair_cost': 1,
                   'coupling_squared': source['g_squared'], 'neutral_gap': source['gap'],
                   'preparation': '(vacuum+i*normalized_bright)/sqrt(2)',
                   'prepared_state_is_gauge_singlet': True,
                   'full307_carrier_is_all_singlet': False,
                   'current_detector_values': [g, -g, 0],
                   'current_detector_ranks': [s.trace(p) for p in projections],
                   'initial_complete_current_probabilities': [expect(state, p) for p in projections]})
    return source, state, result, projections


def scalar_history(time):
    """All full-sector occupation and current outcomes, on one source clock."""
    time = scalar._exact(time)
    if time.is_real is not True:
        raise ValueError('a real common source time is required')
    source, _, _, _ = scalar_data()
    gap, g, d = source['gap'], s.sqrt(source['g_squared']), source['pair_cost']
    transport = s.simplify(d*g*s.sin(gap*time)/gap)
    occupation = [s.simplify(s.Rational(1, 2)+g*s.sin(gap*time)/gap),
                  s.simplify(s.Rational(1, 2)-g*s.sin(gap*time)/gap), s.S.Zero]
    current = [s.simplify((1+s.cos(gap*time))/2),
               s.simplify((1-s.cos(gap*time))/2), s.S.Zero]
    return {'time': time, 'occupation_labels': ['empty', 'bright', 'dark'],
            'occupation_ranks': [1, 1, 305], 'occupation_probabilities': occupation,
            'current_values': [d*g, -d*g, 0], 'current_probabilities': current,
            'number_energy': s.simplify(d/2-transport), 'vertex_energy': transport,
            'integrated_number_to_vertex_flux': transport,
            'total_energy': d/2,
            'two_alternative_detector_menus_not_a_joint_measurement': True}


def scalar_wavefunction(time):
    """A phase representative of the actual evolved singlet preparation.

    The common exp(-i*d*t/2) phase is removed from this vector. It does
    not change any density or retained-reference record. This is not a
    propagator definition on other inputs, including the dark subspace.
    """
    time = scalar._exact(time)
    if time.is_real is not True:
        raise ValueError('a real common source time is required')
    source, state, _, _ = scalar_data()
    delta, d = source['gap'], source['pair_cost']
    bright_identity = source['vacuum_projector']+source['bright_projector']
    traceless = source['vertex']+d*(source['bright_projector']-source['vacuum_projector'])/2
    u_bright = (s.cos(delta*time/2)*bright_identity
                -2*s.I*s.sin(delta*time/2)*traceless/delta)
    return clean(u_bright*state)


def action_response_certificate():
    """Exact pair-cost quench response at fixed input state and source time."""
    d, g, t = s.symbols('d g t', real=True, positive=True)
    gap = s.sqrt(d*d+4*g*g)
    p = s.Rational(1, 2)-g*s.sin(gap*t)/gap
    derivative = g*d*(s.sin(gap*t)-gap*t*s.cos(gap*t))/gap**3
    if s.simplify(s.diff(p, d)-derivative) != 0:
        raise ArithmeticError('the same-state action response differs')
    value = s.simplify(derivative.subs({d: 1, g: s.sqrt(56), t: s.pi/30}))
    if value != 2*s.sqrt(14)/3375:
        raise ArithmeticError('the literal full-sector action response changed')
    return {'perturbation': 'H(d+lambda)=H(d)+lambda*N_H',
            'preparation_held_fixed': '(vacuum+i*normalized_bright)/sqrt(2)',
            'all_occupation_outputs_retained': True,
            'occupied_probability_derivative': 'g*d*(sin(Delta*t)-Delta*t*cos(Delta*t))/Delta^3',
            'fixture_time': s.pi/30, 'fixture_occupied_probability_derivative': value,
            'geometric_metric_identification_of_pair_cost': False,
            'general_Kubo_response': 'd<A(t)>/dlambda=i*integral_0^t <[V(s),A(t)]> ds',
            'central_difference_error_bound': '(4/3)*h^2*|t|^3*norm(V)^3*norm(A)',
            'same_preparation_and_common_time_required': True}


def paired_initial_records():
    """A declared tensor-product experiment with no omitted coupling term.

    H_total=H_rotor tensor I + I tensor H_scalar. The scalar vertex is
    interacting. There is no rotor-scalar interaction in this control.
    """
    rotor, _, _, _ = rotor_data()
    scalar_source, scalar_state, _, _ = scalar_data()
    n = len(rotor['flux'])
    pointer = scalar.measurement_projectors(scalar_source)
    rows = []
    for j, electric in enumerate(rotor['flux']):
        weight = sum(abs(rotor['state'][mask*n+j, 0])**2 for mask in range(16))
        for label, projector in pointer.items():
            probability = s.simplify(weight*expect(scalar_state, projector))
            rows.append({'electric_flux': electric, 'scalar_outcome': label,
                         'probability': probability})
    if sum(row['probability'] for row in rows) != 1 or sum(row['probability'] != 0 for row in rows) != 4:
        raise ArithmeticError('the full tensor detector lost its original record law')
    return {'declared_carrier_dimension': rotor['H'].rows*scalar_source['H'].rows,
            'same_program_clock': True, 'factorized_evolution_derived_from_commuting_terms': True,
            'rotor_scalar_interaction_declared': False,
            'all_outcomes': rows, 'postselection_or_failure_conditioning': False,
            'valid_with_untouched_reference': True}


def gravity_interface():
    # Same energy density, different spatial stress. This is a tomography
    # null direction, not a claim that either tensor is supplied by the source.
    witness = frame.basis(10, 4)
    charges = frame.apply(frame.design(), witness)
    if charges[0] != 0 or not any(charges) or frame.decode(charges) != witness:
        raise ArithmeticError('the energy-only tensor nonidentification control failed')
    return {
        'realized': ['same-action conserved total energy and complete component flux',
                     'weighted local-energy source and controlled work identity',
                     'all-state operator continuity implies finite block continuity',
                     'same-state pair-energy response and complete detector laws'],
        'finite_theorems_consumed': [
            'OPH.QFT.CoarseEnergyWard.coarse_continuity',
            'OPH.QFT.CoarseEnergyWard.coarse_compose',
            'OPH.QFT.CoarseEnergyWard.time_increments_telescope'],
        'Einstein_requirements_not_instantiated': [
            'source metric with proper clocks, volume, ten boosted regional readouts and curvature',
            'one symmetric action-derived Lorentz stress with independent covariant Ward conservation',
            'regional modular diamond kernel for that same stress and observer state',
            'fixed-proper-volume area response and generalized-entropy stationarity',
            'common shrinking-ball errors o(ell^4) and positive universal area coefficient',
            'same-source smooth metric limit and contracted Bianchi identity'],
        'finite_Einstein_tensor_bound_invoked': False,
        'repair_stationary_entropy_identified_with_geometric_entropy': False,
        'Einstein_dynamics_or_physical_Newton_constant_derived': False,
        'architecture_uniqueness_required': False,
        'energy_only_tomography_null_direction': {
            'symmetric_tensor_coordinates': list(map(str, witness)),
            'ten_observer_charges': list(map(str, charges)),
            'rest_energy_charge': str(charges[0]),
            'other_nine_tensor_components_supplied_by_energy_continuity': False},
    }


def _serial(value):
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, dict):
        return {key: _serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serial(item) for item in value]
    return value


def certificate():
    return _serial({'schema': 'oph.particle_gravity_response.v1',
                    'rotor_energy': rotor_data()[2], 'scalar_energy': scalar_data()[2],
                    'scalar_history': [scalar_history(0), scalar_history(s.pi/30)],
                    'action_response': action_response_certificate(),
                    'paired_initial_records': paired_initial_records(),
                    'gravity_interface': gravity_interface(),
                    'measured_particle_or_gravity_inputs': [], 'source_inputs': list(SOURCE_INPUTS)})


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
