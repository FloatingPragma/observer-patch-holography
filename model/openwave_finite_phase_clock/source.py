"""Exact finite phase transfer on the native charged rest pair."""
from fractions import Fraction as F
import sympy as s
from model.openwave_clock_phase_resource.source import exact_source
from model.openwave_passive_clock_records.clock import native_attachment


def clean(x):
    return x.applyfunc(s.simplify) if isinstance(x, s.MatrixBase) else s.simplify(x)


def exchange(c=F(3, 5), sine=F(4, 5)):
    """All eight source/reference coordinates, including inactive occupations."""
    c, sine = s.Rational(c), s.Rational(sine)
    if c*c+sine*sine != 1:
        raise ValueError('an exact unit-circle pair is required')
    N, A, _ = exact_source()
    lower = s.Matrix([[0, 1], [0, 0]])
    R = s.diag(0, 1)
    B = s.kronecker_product(A, lower)+s.kronecker_product(A.H, lower.H)
    H = s.kronecker_product(N, s.eye(2))+s.kronecker_product(s.eye(4), R)
    U = s.eye(8)+(c-1)*(B*B)-s.I*sine*B
    assert B**3 == B
    assert U.H*U == s.eye(8)
    assert H*U == U*H
    return N, A, R, B, H, U


def partial_reference(rho):
    return s.Matrix(4, 4, lambda i, j: sum(rho[2*i+k, 2*j+k] for k in range(2)))


def prepared(c=F(3, 5), sine=F(4, 5), *, twirl=False):
    N, A, R, B, H, U = exchange(c, sine)
    ref = s.eye(2)/2 if twirl else s.ones(2)/2
    initial = s.kronecker_product(s.diag(1, 0, 0, 0), ref)
    rho = clean(U*initial*U.H)
    source = partial_reference(rho)
    assert s.trace(rho) == 1
    assert s.trace(H*rho) == s.Rational(1, 2)
    assert s.trace(B*initial) == s.trace(B*rho) == 0
    branches = []
    for b in range(2):
        block = s.Matrix(4, 4, lambda i, j: rho[2*i+b, 2*j+b])
        branches.append({'reference_energy_outcome': b,
            'probability_exact': str(s.trace(block)),
            'unnormalized_source_matrix_exact': [[str(z) for z in row] for row in block.tolist()],
            'unnormalized_source_energy_over_mass_exact': str(s.trace(N*block))})
    expected = s.Rational(sine)**2/2
    assert s.trace(N*source) == expected
    assert sum(s.sympify(x['probability_exact']) for x in branches) == 1
    return rho, source, {'cosine_exact': str(c), 'sine_exact': str(sine),
        'resource_energy_twirled': twirl, 'complete_reference_records': branches,
        'source_mean_energy_over_mass_exact': str(expected),
        'reference_mean_energy_over_mass_exact': str(s.Rational(1, 2)-expected),
        'reference_initial_energy_over_mass_exact': '1/2',
        'source_plus_reference_energy_balance_exact': '0',
        'coupling_on_off_mean_work_for_stated_input_exact': '0',
        'switch_or_collision_duration_implemented_autonomously': False,
        'source_A_expectation_exact': str(s.trace(A*source))}


def free_source_law(c=F(3, 5), sine=F(4, 5), *, twirl=False):
    _, rho, record = prepared(c, sine, twirl=twirl)
    N, A, *_ = exchange(c, sine)
    phase = s.symbols('phi', real=True)
    D = s.diag(1, 1, s.exp(-s.I*phase), s.exp(-s.I*phase))
    out = clean(D*rho*D.H)
    coherence = clean(s.trace(A*out))
    assert clean(s.diff(coherence, phase)-s.I*coherence) == 0
    if twirl:
        assert coherence == 0
    else:
        assert clean(coherence-s.I*s.Rational(sine)*s.exp(s.I*phase)/2) == 0
    record.update({'phase_definition': 'phi=m*t in the declared native source parameter',
        'free_source_A_expectation_exact': str(coherence),
        'frequency_over_native_rest_energy_exact': '1',
        'isolated_source_linewidth_exact': '0',
        'isolated_coherence_lifetime': 'unbounded for this exact isolated reducing sector',
        'physical_laboratory_time_identified': False,
        'fixed_source_quadrature_is_a_complete_autonomous_detector': False})
    return record


def complete_terminal_records(phase_power=0, *, fresh_terminal_phase=False,
                              twirl=False, c=F(3, 5), sine=F(4, 5)):
    """Every native and two-reference occupation, with no success selection.

    Phase powers are exact multiples of pi/2. A fresh terminal phase has its
    preparation time supplied separately; it is a counterfactual resource.
    """
    if type(phase_power) is not int:
        raise ValueError('an integral quarter-cycle phase is required')
    rho, _, _ = prepared(c, sine, twirl=twirl)
    N, A, R, _, _, _ = exchange(c, sine)
    z = (-s.I)**phase_power
    Ds = s.diag(1, 1, z, z)
    Dr = s.diag(1, z)
    D = s.kronecker_product(Ds, Dr)
    rho = clean(D*rho*D.H)
    terminal = s.ones(2)/2
    if not fresh_terminal_phase:
        terminal = Dr*terminal*Dr.H
    joint = s.kronecker_product(rho, terminal)
    lower = s.Matrix([[0, 1], [0, 0]])
    B = s.kronecker_product(A, s.eye(2), lower)+s.kronecker_product(A.H, s.eye(2), lower.H)
    H = (s.kronecker_product(N, s.eye(4))+
         s.kronecker_product(s.eye(4), R, s.eye(2))+
         s.kronecker_product(s.eye(8), R))
    cr, sr = s.Rational(3, 5), s.Rational(4, 5)
    V = s.eye(16)+(cr-1)*(B*B)-s.I*sr*B
    assert V.H*V == s.eye(16) and V*H == H*V
    final = clean(V*joint*V.H)
    assert s.trace(final) == 1
    assert s.trace(H*final) == s.trace(H*joint) == 1
    rows = [{'source_occupation': i, 'preparation_reference': j, 'readout_reference': k,
             'probability_exact': str(final[4*i+2*j+k, 4*i+2*j+k])}
            for i in range(4) for j in range(2) for k in range(2)]
    assert all(s.sympify(x['probability_exact']) >= 0 for x in rows)
    assert sum(s.sympify(x['probability_exact']) for x in rows) == 1
    p = sum(s.sympify(x['probability_exact']) for x in rows if x['readout_reference'])
    u = s.Rational(sine)
    wave = 0 if twirl else (s.re(s.I**phase_power) if fresh_terminal_phase else 1)
    expected = cr**2/2+sr**2*u*u/2-cr*sr*u*wave/2
    assert p == expected
    return {'phase_quarter_cycles': phase_power,
        'terminal_resource_phase_prepared_at_readout_time': fresh_terminal_phase,
        'initial_preparation_resource_twirled': twirl,
        'complete_sixteen_records': rows,
        'readout_excited_probability_exact': str(p),
        'total_mean_energy_over_mass_exact': '1',
        'terminal_exchange_total_energy_balance_exact': '0',
        'terminal_occupation_measurement_total_energy_work_exact': '0'}


def certificate():
    from .controller import certificate as controller_certificate
    normal = [complete_terminal_records(k) for k in range(4)]
    assert all(x['complete_sixteen_records'] == normal[0]['complete_sixteen_records'] for x in normal)
    fresh = [complete_terminal_records(k, fresh_terminal_phase=True) for k in range(4)]
    assert len({x['readout_excited_probability_exact'] for x in fresh}) == 3
    twirled = [complete_terminal_records(k, fresh_terminal_phase=True, twirl=True) for k in range(4)]
    assert len({x['readout_excited_probability_exact'] for x in twirled}) == 1
    strengths = [(F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(1), F(0))]
    masses = native_attachment()
    return {'schema': 'oph.finite_phase_native_clock.v1',
        'verdict': 'PASS_NATIVE_PHASE_RESOURCE_AND_AUTONOMOUS_FINITE_CONTROLLER',
        'native_source': 'Actual charged pair, complete four-occupation CAR memory and original quantum Gauss dressing; every orthogonal native Fock block retained by its bounded partial isometry',
        'native_mass_and_probe_attachments': masses,
        'preparation': 'H=m*N_source+m*N_reference+g*(A*b+A^dagger*b^dagger); the exact exchange collision has cos(g*tau)=c and sin(g*tau)=s',
        'resource': 'One initial two-level phase state (|0>+|1>)/sqrt(2); a second two-level phase state for the complete terminal experiment',
        'auxiliary_strength_controls': [free_source_law(c, u) for c, u in strengths],
        'twirled_preparation': free_source_law(twirl=True),
        'complete_coevolving_reference_readouts': normal,
        'fresh_terminal_phase_counterfactual': fresh,
        'twirled_fresh_phase_controls': twirled,
        'autonomous_finite_clock_controller': controller_certificate(masses),
        'frequency_statement': 'For every native charged m, source coherence is i*s*exp(i*m*t)/2; g and collision duration affect s, not the free frequency',
        'operational_no_go': 'Any fixed total-energy-compatible terminal instrument after joint free evolution has time-independent complete outcome probabilities, for arbitrary phaseful initial states',
        'linewidth_and_lifetime': 'The isolated source pole has zero linewidth and no intrinsic decay. The finite references are used once. No repeated-resource lifetime, physical seconds, or autonomous switching schedule is established.',
        'scope': {'phaseful_finite_preparation_executed': True,
            'complete_energy_compatible_terminal_instrument_executed': True,
            'co_evolving_reference_exhibits_mass_frequency': False,
            'fresh_terminal_phase_is_autonomously_supplied': False,
            'native_mass_ratios_source_derived_without_target_exposure': False,
            'literal_openwave_row3_validated': False,
            'additional_nature_validations': 0}}
