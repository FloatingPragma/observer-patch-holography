"""Bounded reducing-defect tests on a chosen source response space."""
from fractions import Fraction as F
import sympy as s
from model.native_ordered_ground import source
from . import asymptotic


def _rational(value, name):
    if type(value) not in (int, F):
        raise TypeError(name+' must be an exact rational')
    return F(value)


def finite_defect(clock, embedding, resolution=F(1)):
    """Exact small-matrix certificate; no continuum embedding is inferred."""
    z = _rational(resolution, 'resolution')
    if z <= 0:
        raise ValueError('positive resolution required')
    k, v = s.Matrix(clock), s.Matrix(embedding)
    if k.rows == 0 or k.rows != k.cols or v.rows != k.rows or not v.cols:
        raise ValueError('nonempty compatible source and response spaces required')
    if k.has(s.Float) or v.has(s.Float):
        raise TypeError('exact matrix entries required')
    if k != k.H or any(e.is_nonnegative is not True for e in k.eigenvals()):
        raise ValueError('a positive selfadjoint source clock is required')
    if s.simplify(v.H*v) != s.eye(v.cols):
        raise ValueError('the response embedding must be an isometry')
    bounded = s.Rational(z.numerator, z.denominator)*(k+z*s.eye(k.rows)).inv()
    compression = s.simplify(v.H*bounded*v)
    residual = s.simplify((s.eye(k.rows)-v*v.H)*bounded*v)
    defect = s.simplify(v.H*bounded**2*v-compression**2)
    if s.simplify(defect-residual.H*residual) != s.zeros(v.cols):
        raise AssertionError('bounded reducing-defect identity failed')
    return {'resolvent': bounded, 'compression': compression,
            'residual': residual, 'defect': defect}


def polynomial_coefficient(coefficients):
    """The exact multiplier sum j|c_j| in the telescoping error estimate."""
    values = tuple(_rational(value, 'coefficient') for value in coefficients)
    if not values:
        raise ValueError('a nonempty coefficient sequence is required')
    return sum((j*abs(value) for j, value in enumerate(values)), F(0))


def scalar_score_error(original_weight, first_error, second_error):
    """Known-weight error, with the estimated first score clipped to [0,w]."""
    w, e1, e2 = (_rational(x, name) for x, name in zip(
        (original_weight, first_error, second_error),
        ('original weight', 'first error', 'second error'), strict=True))
    if not 0 < w <= 1 or min(e1, e2) < 0:
        raise ValueError('positive original probability and nonnegative errors required')
    return e2/w+2*e1/w


def high_energy_control(scale=1, resolution=F(1)):
    m, z = _rational(scale, 'scale'), _rational(resolution, 'resolution')
    if min(m, z) <= 0:
        raise ValueError('positive energy scale and resolution required')
    v = s.Matrix([1, 1])/s.sqrt(2)
    k = s.diag(s.Rational(m), 3*s.Rational(m))
    result = finite_defect(k, v, z)
    t = s.pi/(2*s.Rational(m))
    evolved = s.diag(s.exp(-s.I*m*t), s.exp(-3*s.I*m*t))*v
    assert s.simplify((v.H*evolved)[0]) == 0
    assert s.simplify((evolved.H*evolved)[0]) == 1
    expected = s.Rational(z*z*m*m/((m+z)**2*(3*m+z)**2))
    assert result['defect'][0] == expected
    return {'energy_scale': str(m), 'resolution': str(z),
            'defect': str(expected), 'escape_time': str(t),
            'original_probability_outside_response': '1',
            'uniform_time_control_from_small_resolvent_defect': False}


def _native_dimer():
    raw, scale = source.matrix(2, ((0, 1),), F(1), F(3, 8))
    h = s.Matrix(raw.tolist())/scale
    omega = s.Matrix([1, 3, 3, 1])/s.sqrt(20)
    upper = s.Matrix([3, -1, -1, 3])/s.sqrt(20)
    x = s.zeros(4)
    for mask in range(4):
        x[mask ^ 1, mask] = 1
    z = s.diag(1, -1, 1, -1)
    e0 = -s.Rational(9, 8)
    assert h.eigenvals() == {e0: 1, s.Integer(0): 2, s.Rational(1, 8): 1}
    assert h*omega == e0*omega and h*upper == upper/8
    assert s.simplify((omega.H*omega)[0]) == 1
    assert s.simplify((upper.H*upper)[0]) == 1
    assert x*omega == 3*omega/5+4*upper/5
    return h, omega, upper, x, z


def native_dimer_control():
    """Actual source, original outcome weights and complete finite operators."""
    h, omega, upper, x, z = _native_dimer()
    ident = s.eye(4)
    k, resolution = h+9*ident/8, F(1, 8)
    sharp = z*omega
    assert k*sharp == 9*sharp/8
    assert (sharp.H*upper)[0] == 0 and (omega.H*sharp)[0] == 0
    sharp_control = finite_defect(k, sharp, resolution)
    assert sharp_control['defect'] == s.zeros(1)
    bounded = sharp_control['resolvent']
    filt = ident-bounded
    a, prep_fail = (x+z)/2, ident/s.sqrt(2)
    assert a.H*a+prep_fail.H*prep_fail == ident
    p0, p_upper = omega*omega.H, upper*upper.H
    p_middle = ident-p0-p_upper
    filter_fail = p0+s.sqrt(19)*p_middle/10+s.sqrt(21)*p_upper/11
    assert s.simplify(filter_fail**2+filt**2-ident) == s.zeros(4)
    kraus = (prep_fail, filter_fail*a, filt*a)
    assert s.simplify(sum((op.H*op for op in kraus), s.zeros(4))-ident) == s.zeros(4)
    probabilities = [s.simplify((omega.H*op.H*op*omega)[0]) for op in kraus]
    assert sum(probabilities) == 1
    psi = filt*a*omega
    assert psi == 9*sharp/20+4*upper/11
    original_weights = [s.Rational(81, 400), s.Rational(16, 121)]
    weight = sum(original_weights)
    assert weight == s.Rational(16201, 48400) == probabilities[-1]
    scalar = finite_defect(k, psi/s.sqrt(weight), resolution)
    bounded_scores = [s.Rational(1, 10), s.Rational(1, 11)]
    first = sum(w*score for w, score in zip(original_weights, bounded_scores, strict=True))
    second = sum(w*score**2 for w, score in zip(original_weights, bounded_scores, strict=True))
    defect = s.simplify(second/weight-(first/weight)**2)
    assert defect == scalar['defect'][0] and defect > 0
    spectral_variance = s.prod(original_weights)/weight**2*(bounded_scores[0]-bounded_scores[1])**2
    assert spectral_variance == defect
    reducing = finite_defect(k, sharp.row_join(upper), resolution)
    assert reducing['defect'] == s.zeros(2)
    assert reducing['compression'] == s.diag(*bounded_scores)
    return {'source_producer': 'model.native_ordered_ground.source.matrix',
            'sites': 2, 'edges': [[0, 1]], 'coupling': '1', 'field': '3/8',
            'ground_energy': '-9/8', 'source_spectrum': ['-9/8', '0', '0', '1/8'],
            'positive_clock': 'K=H_B+9I/8', 'resolution': '1/8',
            'sharp_local_Z_positive_energy': '9/8', 'sharp_local_Z_defect': '0',
            'original_success_energies': ['9/8', '5/4'],
            'original_success_spectral_weights': list(map(str, original_weights)),
            'complete_record_names': ['preparation_failure', 'filter_failure', 'success'],
            'complete_record_probabilities': list(map(str, probabilities)),
            'original_success_probability': str(weight),
            'original_first_bounded_score': str(first),
            'original_second_bounded_score': str(second),
            'normalized_rank_one_defect': str(defect),
            'two_component_response_defect': [['0', '0'], ['0', '0']],
            'finite_complete_instrument_specified': True,
            'native_pulse_tape_executed': False,
            'continuum_Gram_projection_implemented': False,
            'physical_particle_mass': None}



def general_two_speed_coefficients(order=F(1, 4), bond=F(1, 4)):
    """Match any positive m,kappa; no measured source coefficients are chosen."""
    m, kappa = _rational(order, 'order'), _rational(bond, 'bond coefficient')
    if min(m, kappa) <= 0:
        raise ValueError('positive canonical order and bond coefficient required')
    # Symbolic coefficient identities apply to every positive pair.
    ms, ks = s.symbols('m kappa', positive=True)
    aa, bb1, bb2 = 12*ks/ms, ks/(2*ms), 3*ks/(2*ms)
    assert s.simplify(ms*aa/2-6*ks) == 0
    assert s.simplify(ms*(bb1+bb2)/4-ks/2) == 0
    assert s.simplify(aa*bb1-6*ks**2/ms**2) == 0
    assert s.simplify(aa*bb2-18*ks**2/ms**2) == 0
    assert s.simplify(aa*(bb1+bb2)/2-12*ks**2/ms**2) == 0
    assert s.simplify(aa*(bb1-bb2)/2+6*ks**2/ms**2) == 0
    actual_source_range = F(1, 10) < m <= F(1, 2) and m*m/2 <= kappa <= F(1, 2)
    q_norm = s.Rational(m)*(s.sqrt(6)+s.sqrt(2))/2
    p_norm = s.Rational(m)*(s.sqrt(6)+3*s.sqrt(2))/48
    q_inverse, p_inverse = 2*m*m/(3*kappa), m*m/(24*kappa)
    cone_squared = 54*kappa*kappa/(m*m)
    if actual_source_range:
        assert q_norm < 3 and p_norm < F(1, 4)
        assert q_inverse <= 4 and p_inverse <= F(1, 12)
        assert cone_squared < 72**2
    return {'canonical_coefficient': str(m), 'bond_coefficient': str(kappa),
            'kinetic_coefficient': str(12*kappa/m),
            'gradient_coefficients': [str(kappa/(2*m)), str(3*kappa/(2*m))],
            'Q_first_energy_moment': str(6*kappa),
            'P_first_energy_coefficient_over_momentum_squared': str(kappa/2),
            'two_speeds_squared': [str(6*kappa*kappa/(m*m)), str(18*kappa*kappa/(m*m))],
            'visible_Delta_coefficient': str(12*kappa*kappa/(m*m)),
            'hidden_Delta_coefficient': str(-6*kappa*kappa/(m*m)),
            'Q_norm_coefficient_over_inverse_momentum': str(q_norm),
            'P_norm_coefficient_over_momentum': str(p_norm),
            'Q_inverse_energy_coefficient_over_inverse_momentum_squared': str(q_inverse),
            'P_inverse_energy_moment': str(p_inverse),
            'sufficient_l1_cone_speed_squared': str(cone_squared),
            'within_stated_joint_source_parameter_range': actual_source_range,
            'source_coefficients_selected_by_this_control': False,
            'actual_OPH_source_limit_claimed': False}


def two_pair_circle_control():
    """Double the inherited field control, without asserting a source limit."""
    parent = asymptotic.two_speed_control()
    assert parent['order'] == parent['bond_coefficient'] == '1/4'
    assert parent['kinetic_coefficient'] == '12'
    assert parent['gradient_coefficients'] == ['1/2', '3/2']
    assert parent['two_speeds_squared'] == ['6', '18']
    assert parent['Q_first_energy_moment'] == '3/2'
    assert parent['P_first_energy_coefficient_over_momentum_squared'] == '1/8'
    mixing = s.Matrix([[1, 1], [1, -1]])/s.sqrt(2)
    diagonal_gradient = s.diag(*(s.Rational(x) for x in parent['gradient_coefficients']))
    gradient = s.simplify(mixing*diagonal_gradient*mixing.T)
    assert gradient == s.Matrix([[1, -s.Rational(1, 2)], [-s.Rational(1, 2), 1]])
    # Coordinates (q1,x1,q2,x2) and the corresponding canonical momenta.
    doubled_gradient = s.kronecker_product(s.eye(2), gradient)
    a = s.Rational(parent['kinetic_coefficient'])
    hessian = s.diag(doubled_gradient, a*s.eye(4))
    j2 = s.Matrix([[0, 1], [-1, 0]])
    internal = s.kronecker_product(j2, s.eye(2))
    rotation_generator = s.diag(internal, internal)
    symplectic = s.BlockMatrix([[s.zeros(4), s.eye(4)], [-s.eye(4), s.zeros(4)]]).as_explicit()
    assert rotation_generator.T*hessian+hessian*rotation_generator == s.zeros(8)
    assert rotation_generator.T*symplectic+symplectic*rotation_generator == s.zeros(8)
    assert rotation_generator**2 == -s.eye(8)
    # Exponentiation therefore gives a full continuous circle, not only its generator.
    c, t = s.symbols('c t', real=True)
    rotation = c*s.eye(8)+t*rotation_generator
    assert rotation.T*hessian*rotation == (c*c+t*t)*hessian
    assert rotation.T*symplectic*rotation == (c*c+t*t)*symplectic
    m = s.Rational(1, 4)
    # Projection onto visible q1,q2 and their momenta, with the source CCR scaling.
    visible = s.zeros(4, 8)
    for row, column in enumerate((0, 2, 4, 6)):
        visible[row, column] = s.sqrt(m)
    visible_symplectic = m*s.BlockMatrix([[s.zeros(2), s.eye(2)], [-s.eye(2), s.zeros(2)]]).as_explicit()
    assert visible*symplectic*visible.T == visible_symplectic
    frequency_squared = a*gradient
    assert frequency_squared.eigenvals() == {s.Integer(6): 1, s.Integer(18): 1}
    assert frequency_squared[0, 1] == -6
    return {'inherited_single_pair_control': parent,
            'visible_regular_CCR_pairs': 2,
            'full_internal_U1_preserves_clock_and_vacuum': True,
            'gradient_matrix': [[str(value) for value in gradient.row(row)] for row in range(2)],
            'two_speeds_squared': ['6', '18'],
            'visible_wave_equation': 'd_t^2 Q=12 Delta Q-6 Delta X',
            'hidden_force_coefficient': '-6',
            'scalar_pure_direction_twist_Hessians': True,
            'projected_memory_proves_scattering': False,
            'actual_OPH_source_limit_claimed': False,
            'source_specific_force_estimate_replaced': False}


def certificate():
    return {'bounded_resolvent_reducing_identity': 'D=R*R',
            'energy_moment_domain_required': False,
            'arbitrary_response_projection_implemented': False,
            'native_dimer': native_dimer_control(),
            'high_energy_control': [high_energy_control(n) for n in (1, 10, 100)],
            'two_pair_internal_circle_control': two_pair_circle_control(),
            'general_positive_coefficient_control': general_two_speed_coefficients(),
            'actual_continuum_response_defect_proved_zero': False,
            'actual_source_dispersion_identified': False}
