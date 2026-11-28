"""Execute the link, CAR and rotor parts of the relative weak-current meter."""

from collections import defaultdict
from itertools import product

import sympy as sp

from ..charged_endpoints.cube import A, BETA_SEA, CELLS, CENTER, MASKS, c_action


DOUBLET = (15, 23)
N = sp.symbols("nx ny nz", real=True)
LINK = sp.symbols("u00 u01 u10 u11")
U = sp.Matrix(2, 2, LINK)
SIGMA = (
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    sp.diag(1, -1),
)
N_SIGMA = sum((a * b for a, b in zip(N, SIGMA, strict=True)), sp.zeros(2))
HALF = sp.Rational(1, 2)


def clean(matrix):
    return matrix.applyfunc(sp.expand)


def gauge_checks():
    assert all(s.H == s and s.trace() == 0 for s in SIGMA)
    gauge_actions = 0
    for r in (1, -1):
        frame = (sp.eye(2) + r * N_SIGMA) / 2
        for axis, sigma in enumerate(SIGMA):
            # Under exp(i sigma/2), the spherical coordinate rotates with
            # the opposite signed spatial cross product in this chart.
            direction = sp.eye(3)[:, axis]
            delta_n = sp.Matrix(N).cross(direction)
            delta_p = sum(
                (frame.diff(n) * dn for n, dn in zip(N, delta_n, strict=True)),
                sp.zeros(2),
            )
            generator = sp.I * sigma / 2
            assert clean(delta_p - (generator * frame - frame * generator)) == sp.zeros(
                2
            )
            delta_x = delta_p * U + frame * generator * U
            delta_y = -frame * U * generator
            assert clean(delta_x - generator * frame * U) == sp.zeros(2)
            assert clean(delta_y + frame * U * generator) == sp.zeros(2)
            adjoint = (frame * U).H
            assert clean(delta_x.H + adjoint * generator) == sp.zeros(2)
            assert clean(delta_y.H - generator * adjoint) == sp.zeros(2)
            gauge_actions += 4
    return gauge_actions


def mode(cell, beta, spin, component):
    return (
        60 * CELLS.index(cell)
        + 15 * (spin + (2 if beta == -1 else 0))
        + MASKS.index(component)
    )


def one_body(word, create_mode, drop_mode):
    if not (word >> drop_mode) & 1:
        return None
    dropped, parity1 = c_action(word, drop_mode, False)
    if (dropped >> create_mode) & 1:
        return None
    result, parity2 = c_action(dropped, create_mode, True)
    return result, parity1 * parity2


def pair_column(x, y, spin, sign, oriented_coefficient):
    particle, hole = (x, y) if sign == 1 else (y, x)
    coeff = oriented_coefficient if sign == 1 else oriented_coefficient.H
    vector = {}
    for i, j in product(range(2), repeat=2):
        outcome = one_body(
            BETA_SEA,
            mode(particle, 1, spin, DOUBLET[i]),
            mode(hole, -1, spin, DOUBLET[j]),
        )
        assert outcome is not None
        new, parity = outcome
        vector[new] = sp.expand(parity * coeff[i, j])
    assert len(vector) == 4
    return vector


def rotor_projection(expr):
    """Exact projection to span(1,nx,ny,nz) for polynomials through degree two."""
    out = 0
    for powers, coefficient in sp.Poly(sp.expand(expr), *N).terms():
        if sum(powers) <= 1:
            out += coefficient * sp.prod(n**p for n, p in zip(N, powers, strict=True))
        elif 2 in powers:
            assert sum(powers) == 2
            out += coefficient / 3
        else:
            assert (
                sum(powers) == 2
            )  # mixed second moments and all odd third moments vanish
    return sp.expand(out)


def projected_current(vector, x, spin, sign):
    output = defaultdict(lambda: 0)
    for word, amplitude in vector.items():
        for i, j in product(range(2), repeat=2):
            result = one_body(
                word, mode(x, sign, spin, DOUBLET[i]), mode(x, sign, spin, DOUBLET[j])
            )
            if result is not None:
                new, parity = result
                output[new] += parity * N_SIGMA[i, j] * amplitude / 2
    return {word: rotor_projection(value) for word, value in output.items()}


def gauge_transfer(word, weights):
    values = []
    for cell in CELLS:
        begin = 60 * CELLS.index(cell)
        values.append(
            sum(
                weights[mask]
                * (((word >> (begin + 15 * beta_spin + i)) & 1) - int(beta_spin >= 2))
                for beta_spin in range(4)
                for i, mask in enumerate(MASKS)
            )
        )
    return values


def haar_coefficients(vector, sign):
    """Use the exact Haar Gram matrices of the degree-one link and frame."""
    variables = LINK if sign > 0 else tuple(sp.conjugate(u) for u in LINK)
    components = []
    for expr in vector.values():
        for u in variables:
            coefficient = sp.expand(expr).coeff(u)
            assert coefficient == rotor_projection(coefficient)
            components.append(
                [
                    coefficient.subs(dict.fromkeys(N, 0)),
                    *(coefficient.coeff(n) for n in N),
                ]
            )
    return sp.Matrix(components)


def exact_sector(x, y, spin, sign, weights):
    """Run the two independent upper/lower trials through the same pulse."""
    columns = [
        pair_column(x, y, spin, sign, clean((sp.eye(2) + r * N_SIGMA) * U / 2))
        for r in (1, -1)
    ]
    expected_shift = [-3 * sign * (int(z == x) - int(z == y)) for z in CELLS]
    for r, column in zip((1, -1), columns, strict=True):
        output = projected_current(column, x, spin, sign)
        assert set(output) <= set(column)
        assert all(
            sp.expand(output.get(word, 0) - sign * r * HALF * coefficient) == 0
            for word, coefficient in column.items()
        )
        assert all(gauge_transfer(word, weights) == expected_shift for word in column)
    a, b = (haar_coefficients(c, sign) for c in columns)
    # Link matrix coefficient inner products are delta/d. Rotor Gram is
    # diag(1,1/3,1/3,1/3). U(1) phases have unit Haar norm.
    frame_metric = sp.diag(1, sp.Rational(1, 3), sp.Rational(1, 3), sp.Rational(1, 3))
    for left, right, wanted in ((a, a, 1), (b, b, 1), (a, b, 0)):
        assert (
            sp.expand(sp.trace(left.conjugate() * frame_metric * right.T) / 2) == wanted
        )
    source_a = a * frame_metric * a.H / 2
    source_b = b * frame_metric * b.H / 2
    assert source_a == source_b
    assert (a[:, 0].H * a[:, 0])[0] / 2 == HALF
    return {
        "direction": [list(y), list(x)],
        "receiver_beta_sign": sign,
        "spin_component": spin,
        "signed_native_number_transfers": expected_shift,
        "joint_frame_overlap_between_preparations": 0,
        "source_density_after_unreading_frame_equal": True,
        "source_gauss_projection_both_keep_weight": "1/2",
        "source_gauss_projected_upper_lower_columns_equal": True,
        "upper_lower_current_generator_eigenvalues": [
            [str(sign * HALF), str(-sign * HALF)],
            [str(-sign * HALF), str(-sign * HALF)],
        ],
    }


def singlet_car_commutators():
    """Check full local occupations for one singlet hop and a weak receiver."""
    singular = (0, 1)  # addressed singlet modes at two cells
    weak = (2, 3)  # either weak doublet at the addressed cell
    for receiver, donor in (singular, singular[::-1]):
        for word in range(16):
            hopped = one_body(word, receiver, donor)
            if hopped is not None:
                new, _ = hopped
                assert ((new >> 0) & 1) - ((word >> 0) & 1) == (
                    1 if receiver == 0 else -1
                )
            for i, j in product(weak, repeat=2):
                turned = one_body(word, i, j)
                after_hop = one_body(hopped[0], i, j) if hopped else None
                before_hop = one_body(turned[0], receiver, donor) if turned else None
                assert (after_hop is None) == (before_hop is None)
                if after_hop is not None:
                    assert before_hop is not None
                    assert after_hop[0] == before_hop[0]
                    assert after_hop[1] * hopped[1] == before_hop[1] * turned[1]
    return 2 * 16 * 4


def rotor_pulse_domain(local_doublets):
    """Check a finite compression that keeps all weak Fock spectators."""
    doubled_weights = {0}
    for _ in range(local_doublets):
        # Each fundamental fermion pair has two singlets (empty and filled)
        # and two singly occupied states of doubled weight +1 and -1.
        doubled_weights = {weight + w for weight in doubled_weights for w in (-1, 0, 1)}
    assert sorted(doubled_weights) == list(range(-local_doublets, local_doublets + 1))
    # Pointwise rotation to the third axis fixes this list of eigenvalues.
    # Spectral interpolation needs degree at most one less than its length,
    # independently of the weak spectator occupations or rotor coordinates.
    pulse_degree = len(doubled_weights) - 1
    rotation_degree = 2 * pulse_degree + 1
    assert local_doublets == 16 and rotation_degree == 65
    # Admission and read columns use only the smaller space. Its three
    # compressed coordinates are not mutually commuting as full pulses.
    basis = (sp.S(1), *N)
    coordinate = [
        sp.Matrix.hstack(
            *(
                sp.Matrix(
                    [
                        p.coeff(n) if n in N else p.subs(dict.fromkeys(N, 0))
                        for n in basis
                    ]
                )
                for element in basis
                for p in (rotor_projection(n_axis * element),)
            )
        )
        for n_axis in N
    ]
    assert coordinate[0] * coordinate[1] != coordinate[1] * coordinate[0]
    assert clean(N_SIGMA**2 - sum(n**2 for n in N) * sp.eye(2)) == sp.zeros(2)
    return {
        "total_weak_fundamental_copies_in_full_receiver_Fock_space": local_doublets,
        "doubled_weight_spectral_set": sorted(doubled_weights),
        "maximum_multiplication_degree_in_each_full_spectator_pulse": pulse_degree,
        "maximum_degree_after_two_pulses_and_one_preparation": rotation_degree,
        "declared_rotor_harmonic_degree": rotation_degree,
        "declared_rotor_finite_dimension": (rotation_degree + 1) ** 2,
        "admitted_and_emitted_control_column_rotor_dimension": 4,
        "four_mode_compressed_coordinates_noncommuting": True,
    }


def verify_link_and_frame(weights, local_doublets):
    assert sum(weights[x] for x in MASKS) == 0
    assert all(weights[c] == -3 for c in DOUBLET)
    # Symbolic SU(2) infinitesimal covariance implies exact connected group
    # covariance. The colour factor is identity and U(1) compensates -3.
    verified = gauge_checks()
    examples = [
        exact_sector(x, y, spin, sign, weights)
        for (x, y), spin, sign in product(((A, CENTER), (CENTER, A)), (0, 1), (1, -1))
    ]
    return {
        "local_weak_covariance_identities": verified,
        "full_fock_rotor_pulse_domain": rotor_pulse_domain(local_doublets),
        "symbolic_car_rotor_controls": examples,
        "weak_rotor_moment_gram": ["1", "1/3", "1/3", "1/3"],
        "wilson_coefficient_haar_squared_norm": "1/2",
        "link_polynomial_degree": 1,
        "singlet_transfer_weak_current_CAR_commutators_checked": singlet_car_commutators(),
    }
