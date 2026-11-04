"""Full-CAR local norm from the actual origin hopping block."""
import sympy as s

from model.radial_interacting_core.algebra import source_matrix


def annihilator(mode, modes=4):
    if type(mode) is not int or type(modes) is not int or not 0 <= mode < modes <= 6:
        raise ValueError('a valid small full-CAR carrier is required')
    result = s.zeros(1 << modes)
    for state in range(1 << modes):
        if state & (1 << mode):
            result[state ^ (1 << mode), state] = (-1)**(state & ((1 << mode)-1)).bit_count()
    return s.ImmutableSparseMatrix(result)


def local_normal_form():
    identity = s.eye(16)
    a = [annihilator(j) for j in range(4)]
    n = [x.H*x for x in a]
    D = (n[0]+n[1]-identity)**2
    parity = 2*D-identity
    star = sum((a[j].H*a[j+2]+a[j+2].H*a[j] for j in range(2)), s.zeros(16))
    if D*D != D or parity*parity != identity or star*parity+parity*star != s.zeros(16):
        raise ArithmeticError('origin projector or incident-hopping parity failed')
    commutator = star*D-D*star
    if commutator.H*commutator != star*star:
        raise ArithmeticError('projector and sign-parity commutator factors were confused')
    eigenvalues = star.eigenvals()
    if eigenvalues != {-2: 1, -1: 4, 0: 6, 1: 4, 2: 1}:
        raise ArithmeticError('full-Fock star spectrum differs from its two-channel CAR sum')
    b, U = s.symbols('b U', real=True)
    local = b*star+U*parity/4
    if (local*local-b*b*star*star-U*U*identity/16).applyfunc(s.expand) != s.zeros(16):
        raise ArithmeticError('center contact and star did not anticommute')
    return {
        'full_CAR_modes_in_reduced_star': 4,
        'complete_reduced_star_Fock_rows': 16,
        'star_over_singular_value_spectrum': {str(key): value for key, value in sorted(eigenvalues.items())},
        'origin_projector': 'D0=(n_up+n_down-1)^2',
        'origin_sign_parity': 'P0=2D0-I=(2n_up-I)(2n_down-I)',
        'center_contact': 'W0=(U/4)P0',
        'star_Fock_norm_exact': '2b',
        'projector_commutator_norm_exact': '2b',
        'sign_parity_commutator_norm_exact': '4b',
        'star_plus_center_contact_norm_squared_exact': '4b^2+U^2/16',
        'full_CAR_normal_form_replaces_full_source_Fock_carrier': False,
    }


def source_origin_block():
    h, sites = source_matrix()
    index = sites.index((0, 0, 0))
    center = [2*index, 2*index+1]
    remote = [j for j in range(h.rows) if j not in center]
    V = h.extract(remote, center)
    gram = s.simplify(V.H*V)
    if h.extract(center, center) != s.zeros(2) or gram != s.Rational(3, 2)*s.eye(2):
        raise ArithmeticError('the actual bounded radial source has a different origin block')
    return {
        'actual_parent_source_one_body_rows': h.rows,
        'origin_block_zero_exact': True,
        'origin_coupling_Gram_at_a_v_one': [[str(x) for x in gram.row(j)] for j in range(2)],
        'two_equal_remote_channel_singular_values_squared': '3/2',
        'general_dimension_d_singular_value_squared': 'd*v^2/(2*a^2)',
        'origin_star_support_sites_in_three_dimensions': 7,
        'far_source_is_identity_on_origin_factor': True,
        'radial_mass_and_remote_terms_deleted_from_actual_dynamics': False,
    }


def certificate():
    return {'actual_origin_block': source_origin_block(), 'exact_complete_CAR_star': local_normal_form()}
