"""Exact illustrative four-bin instruments for the horizontal particle atlas.

The state and Fourier analyser are supplied pedagogical data. No native
physical position, particle identity or universal effect implementation is
inferred by this calculation.
"""
import sympy as s


SOURCE_INPUTS = (
    'Lean/EventAlgebra/FiniteBuschGleason.lean',
    'Lean/EventAlgebra/OperationalBornBridge.lean',
    'Lean/EventAlgebra/Robertson.lean',
    'Lean/EventAlgebra/InstrumentSelectionByCertainStates.lean',
    'model/native_quantum_reconstruction/PROOF.md',
    'model/native_born_bridge/PROOF.md',
    'model/carrier_position_translation/PROOF.md',
    'model/spin_covariant_localization/CONTINUUM_POSITION.md',
    'model/spin_dynamical_localization/PROOF.md',
    'model/spin_operational_position/PROOF.md',
)


def clean(value):
    if isinstance(value, s.MatrixBase):
        return s.ImmutableMatrix(value.applyfunc(s.simplify))
    return s.simplify(value)


def construction():
    fourier = s.ImmutableMatrix(4, 4, lambda k, j: (-s.I)**(k*j)/2)
    psi = s.ImmutableMatrix([1, 1, 0, 0])/s.sqrt(2)
    position = tuple(s.ImmutableMatrix(4, 4, lambda j, k: int(j == k == r))
                     for r in range(4))
    momentum = tuple(clean(fourier.H*p*fourier) for p in position)
    return fourier, psi, position, momentum


def outcome_probabilities(rho, instrument):
    return tuple(clean(s.trace(k*rho*k.H)) for k in instrument)


def conditional_state(rho, kraus):
    unnormalized = clean(kraus*rho*kraus.H)
    probability = clean(s.trace(unnormalized))
    if probability == 0:
        return probability, None
    return probability, clean(unnormalized/probability)


def unread_state(rho, instrument):
    return clean(sum((k*rho*k.H for k in instrument), s.zeros(4)))


def sequential_records(rho, first, second):
    """Every ordered pair is retained, including zero-probability records."""
    return [{'first': i, 'second': j,
             'probability': str(clean(s.trace(b*a*rho*a.H*b.H)))}
            for i, a in enumerate(first) for j, b in enumerate(second)]


def finite_robertson(rho, x, p):
    mx, mp = clean(s.trace(rho*x)), clean(s.trace(rho*p))
    vx = clean(s.trace(rho*x*x)-mx*mx)
    vp = clean(s.trace(rho*p*p)-mp*mp)
    comm = clean(s.trace(rho*(x*p-p*x)))
    bound = clean(comm*s.conjugate(comm)/4)
    if clean(vx*vp-bound).is_nonnegative is not True:
        raise ArithmeticError('finite Robertson inequality failed')
    return {'variance_position_bins': str(vx), 'variance_momentum_bins': str(vp),
            'variance_product': str(clean(vx*vp)),
            'commutator_expectation': str(comm),
            'Robertson_rhs_squared': str(bound),
            'trace_commutator': str(clean(s.trace(x*p-p*x))),
            'canonical_CCR_claimed': False}


def double_path_construction(phase=0, overlap=1):
    """Exact supplied two-path input, optional marker, and four-port screen.

    The marker convention is gamma=<d0|d1>. Its normalized states are
    d0=(1,0), d1=(gamma,sqrt(1-|gamma|^2)). The unused input paths retain
    complete isometric action, although this preparation gives them zero
    weight. A screen is a Fourier-output detector in this finite illustration.
    """
    phase, overlap = s.sympify(phase), s.sympify(overlap)
    if phase.has(s.Float) or overlap.has(s.Float):
        raise ValueError('exact phase and marker overlap are required')
    if phase.is_real is not True:
        raise ValueError('a real relative phase is required')
    modulus_squared = clean(overlap*s.conjugate(overlap))
    if clean(1-modulus_squared).is_nonnegative is not True:
        raise ValueError('marker overlap must lie in the closed unit disk')
    f, _, _, _ = construction()
    d0 = s.ImmutableMatrix([1, 0])
    d1 = s.ImmutableMatrix([overlap, s.sqrt(1-modulus_squared)])
    marker = s.zeros(8, 4)
    for path, tag in enumerate((d0, d1, d0, d0)):
        marker[2*path:2*path+2, path] = tag
    marker = clean(marker)
    if clean(marker.H*marker) != s.eye(4):
        raise ArithmeticError('the full marker coupling is not an isometry')
    psi = s.ImmutableMatrix([1, s.exp(s.I*phase), 0, 0])/s.sqrt(2)
    state_before_screen = clean(marker*psi)
    state_at_screen = clean(s.kronecker_product(f, s.eye(2))*state_before_screen)
    density = clean(state_before_screen*state_before_screen.H)
    path_density = clean(s.Matrix(4, 4, lambda j, k:
        sum(density[2*j+r, 2*k+r] for r in range(2))))
    screen = tuple(clean(sum(state_at_screen[2*k+r]*s.conjugate(state_at_screen[2*k+r])
                             for r in range(2))) for k in range(4))
    formula = tuple(clean((1+s.re(overlap*s.exp(s.I*phase)*(-s.I)**k))/4)
                    for k in range(4))
    if any(clean(actual-expected) != 0 for actual, expected in zip(screen, formula)):
        raise ArithmeticError('the marker-overlap interference formula differs')
    if clean(sum(screen)) != 1:
        raise ArithmeticError('screen outcomes are incomplete')
    joint = tuple({'marker_basis_outcome': r, 'screen_outcome': k,
                   'probability': clean(state_at_screen[2*k+r]*s.conjugate(state_at_screen[2*k+r]))}
                  for r in range(2) for k in range(4))
    return {'phase': phase, 'overlap': overlap, 'marker_isometry': marker,
            'input_state': psi, 'path_density_after_unread_marker': path_density,
            'joint_pure_state_at_screen': state_at_screen,
            'screen_probabilities': screen, 'joint_marker_screen_records': joint,
            'phase_scan_visibility': s.sqrt(modulus_squared)}


def double_path_certificate():
    """Complete graphical records for alternative marked/untagged experiments."""
    unmarked = double_path_construction(0, 1)
    marked = double_path_construction(0, 0)
    shifted = double_path_construction(s.pi, 1)
    partial = double_path_construction(0, s.Rational(1, 2))
    return {
        'schema': 'oph.horizontal_atlas.double_path.v1',
        'scope': 'supplied finite two-path interference apparatus; not an emergent laboratory double slit',
        'input_state': '(|0>+exp(i phi)|1>)/sqrt(2)',
        'same_preparation_in_both_apparatus_branches': True,
        'fourier_convention': 'F4[k,j]=(-i)^(k*j)/2',
        'screen_labels': ['D0', 'D1', 'D2', 'D3'],
        'marker_overlap_convention': 'gamma=<d0|d1>',
        'screen_formula': 'p_k(phi,gamma)=[1+Re(gamma exp(i phi)(-i)^k)]/4',
        'visibility': '(max_phi p_k-min_phi p_k)/(max_phi p_k+min_phi p_k)=|gamma|',
        'visibility_requires_phase_scan_at_fixed_screen_bin': True,
        'unmarked': {
            'phase': '0', 'marker_overlap': '1',
            'screen_probabilities': list(map(str, unmarked['screen_probabilities'])),
            'path_record_exists': False,
            'zero_weight_screen_outcomes_retained': [2]},
        'orthogonal_path_tags': {
            'phase': '0', 'marker_overlap': '0',
            'screen_probabilities': list(map(str, marked['screen_probabilities'])),
            'path_screen_records': [{'path_record': row['marker_basis_outcome'],
                                     'screen_record': row['screen_outcome'],
                                     'probability': str(row['probability'])}
                                    for row in marked['joint_marker_screen_records']],
            'screen_given_each_path_record': ['1/4']*4,
            'each_path_record_probability': '1/2',
            'unread_path_record_suffices': True,
            'human_observation_required': False,
            'all_joint_records_retained': True},
        'unmarked_phase_pi_screen_probabilities': list(map(str, shifted['screen_probabilities'])),
        'partial_overlap_one_half_screen_probabilities': list(map(str, partial['screen_probabilities'])),
        'partial_overlap_one_half_phase_scan_visibility': '1/2',
        'endpoint_record_fields': ['run', 'preparation_phase', 'apparatus_setting',
                                   'path_record_if_available', 'screen_record'],
        'consensus_action': 'copy and compare the realized record; preserve its original probability',
        'source_selected_apparatus_or_particle_identity': False,
        'postselection_used': False,
    }


def certificate():
    f, psi, px, pk = construction()
    rho = clean(psi*psi.H)
    if clean(f.H*f) != s.eye(4):
        raise ArithmeticError('Fourier analyser is not unitary')
    for instrument in (px, pk):
        if clean(sum((k.H*k for k in instrument), s.zeros(4))) != s.eye(4):
            raise ArithmeticError('incomplete measurement instrument')
    before_x = outcome_probabilities(rho, px)
    before_k = outcome_probabilities(rho, pk)
    _, selected_x0 = conditional_state(rho, px[0])
    _, selected_k0 = conditional_state(rho, pk[0])
    # A pointer coupling is an isometry; copying its orthogonal record labels
    # is classical record sharing, with no cloning of an arbitrary state.
    pointer = s.zeros(16, 4)
    for r, k in enumerate(px):
        pointer[4*r:4*r+4, :] = k
    if clean(pointer.H*pointer) != s.eye(4):
        raise ArithmeticError('complete record coupling failed')
    x = s.diag(0, 1, 2, 3)
    p = clean(f.H*x*f)
    branches = []
    for setting, instrument, other in [('position', px, pk), ('momentum', pk, px)]:
        for outcome, kraus in enumerate(instrument):
            probability, state = conditional_state(rho, kraus)
            branches.append({'setting': setting, 'outcome': outcome,
                             'probability': str(probability),
                             'conditional_other_distribution': None if state is None
                             else list(map(str, outcome_probabilities(state, other))),
                             'zero_weight_record_retained': probability == 0,
                             'normalized_conditional_state_exists': state is not None})
    return {
        'schema': 'oph.horizontal_atlas.four_bin_measurement.v1',
        'scope': 'supplied pedagogical state and Fourier instruments',
        'native_physical_detector_derived': False,
        'all_particle_families_share_this_position_PVM': False,
        'preparation': {'amplitudes': list(map(str, psi)),
                        'state_formula': '(1,1,0,0)/sqrt(2)',
                        'relative_phase_between_occupied_bins': '0',
                        'source_selection_proved': False},
        'position_labels': ['x0', 'x1', 'x2', 'x3'],
        'momentum_labels': ['k0', 'k1', 'k2', 'k3'],
        'position_probabilities': list(map(str, before_x)),
        'momentum_probabilities': list(map(str, before_k)),
        'momentum_amplitudes': list(map(str, clean(f*psi))),
        'after_position_x0': {
            'branch_probability': str(before_x[0]),
            'position_probabilities': list(map(str, outcome_probabilities(selected_x0, px))),
            'momentum_probabilities': list(map(str, outcome_probabilities(selected_x0, pk)))},
        'after_momentum_k0': {
            'branch_probability': str(before_k[0]),
            'position_probabilities': list(map(str, outcome_probabilities(selected_k0, px))),
            'momentum_probabilities': list(map(str, outcome_probabilities(selected_k0, pk)))},
        'all_conditional_branches': branches,
        'unread_position_then_momentum': list(map(str, outcome_probabilities(unread_state(rho, px), pk))),
        'position_then_momentum_records': sequential_records(rho, px, pk),
        'momentum_then_position_records': sequential_records(rho, pk, px),
        'public_copies': [{'observer_A': i, 'observer_B': j,
                           'probability': str(before_x[i] if i == j else 0)}
                          for i in range(4) for j in range(4)],
        'record_sharing_changes_Born_weights': False,
        'finite_uncertainty': finite_robertson(rho, x, p),
        'finite_position_eigenstate_uncertainty': finite_robertson(selected_x0, x, p),
        'Born_probability_formula': 'p(r)=Tr(rho K_r^dagger K_r)',
        'conditional_state_formula': 'rho_r=K_r rho K_r^dagger/p(r), only when p(r)>0',
        'complete_outcomes_retained': True,
        'all_results_are_exact': True,
        'physical_lengths_momenta_or_hbar_calibrated': False,
        'double_path': double_path_certificate(),
    }
