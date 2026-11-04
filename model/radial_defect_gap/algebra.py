"""Exact source symmetries and the two-row quaternion frame identities."""
import sympy as s

from ..native_mass_channels import algebra as source


def certificate():
    sigma = source.SIGMA
    J = s.I*sigma[1]
    if J*J.conjugate() != -s.eye(2):
        raise ArithmeticError('source spin antiunitary lost its square minus one')
    for item in sigma:
        if J*item.conjugate()*J.H != -item or J*item.T*J.H != -item:
            raise ArithmeticError('source spin conjugation convention changed')
        link = s.I*item/2
        if -J*link.conjugate()*J.H != -link:
            raise ArithmeticError('staggered antiunitary does not reverse source hopping')
        if -J*link.H.T*J.H != -link:
            raise ArithmeticError('staggered particle-hole map does not preserve the quadratic source')
    origin_coupling = sum((2*(s.I*item/2).H*(s.I*item/2) for item in sigma), s.zeros(2))
    if origin_coupling != s.Rational(3, 2)*s.eye(2):
        raise ArithmeticError('the complete six-neighbor origin coupling changed')
    z1, z2 = s.symbols('z1 z2', complex=True)
    first = s.Matrix([z1, z2])
    density = s.conjugate(z1)*z1+s.conjugate(z2)*z2
    for parity in (-1, 1):
        frame = s.Matrix.hstack(first, -parity*J*first.conjugate())
        if ((frame.H*frame-density*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
                or (frame*frame.H-density*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
                or s.expand(s.conjugate(frame.det())*frame.det()-density**2) != 0):
            raise ArithmeticError('a local Kramers frame is not a scalar multiple of a unitary')
    rho = s.symbols('rho', real=True)
    charged = 1-(1-rho)+2*((1-rho)/2)**2
    singly = 1-1+2*((1+rho)/2)*((1-rho)/2)
    if s.expand(charged-(1+rho*rho)/2) != 0 or s.expand(charged-singly-rho*rho) != 0:
        raise ArithmeticError('the exact origin number-deviation contrast changed')
    return {'original_source_spin_rows': 2, 'folded_source_valleys': 8,
            'origin_coupling_Gram_in_units_v2_over_a2': '3/2 * I2',
            'stagger': 'D(n)=(-1)^(n1+n2+n3)',
            'anticommuting_antiunitary': 'C=D*(i*sigma2)*complex_conjugation',
            'commuting_Kramers_antiunitary': 'Theta=R_inv*C',
            'Kramers_square': -1,
            'particle_hole': 'c_n -> D(n)*(i*sigma2)*c_n^dagger',
            'local_normalized_zero_frame': 'W_n^dagger*W_n=W_n*W_n^dagger=rho_n*I2',
            'local_zero_frame_determinant_squared': 'rho_n^2',
            'local_weights_sum': 'sum_n rho_n=1',
            'local_weights_inversion_even': True,
            'origin_negative_sea_density': '(1-rho_0)*I2/2',
            'origin_D0_observable': '(n_origin-1)^2',
            'origin_D0_empty_or_full_core_expectation': '(1+rho_0^2)/2',
            'origin_D0_singly_occupied_core_expectation': '(1-rho_0^2)/2',
            'origin_D0_charged_minus_single_contrast': 'rho_0^2',
            'origin_contrast_is_noninteracting_reference_identity': True,
            'interacting_ground_or_charge_gap_proved_in_this_package': False}
