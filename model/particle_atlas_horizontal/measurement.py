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
    }
