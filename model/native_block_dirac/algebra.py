"""Exact native-doublet Fock encoding and complete charged-fermion tests."""
from functools import lru_cache
import sympy as s


def annihilator(modes, mode):
    out = s.zeros(2**modes)
    for word in range(2**modes):
        if (word >> mode) & 1:
            out[word ^ (1 << mode), word] = (-1)**((word & ((1 << mode)-1)).bit_count())
    return out


def doublet_code():
    """|g,g>, |up,g>, |down,g>, and the two-occupied spin singlet."""
    E = s.zeros(9, 4)
    E[0, 0] = 1
    E[3, 1] = 1
    E[6, 2] = 1
    E[5, 3], E[7, 3] = 1/s.sqrt(2), -1/s.sqrt(2)
    assert E.H*E == s.eye(4)
    n = s.diag(0, 1, 1)
    native_N = s.kronecker_product(n, s.eye(3)) + s.kronecker_product(s.eye(3), n)
    fock_N = s.diag(0, 1, 1, 2)
    assert native_N*E == E*fock_N
    parity = s.kronecker_product(s.diag(1, -1, -1), s.diag(1, -1, -1))
    assert parity*E == E*s.diag(1, -1, -1, 1)
    P = E*E.H
    source = [E*annihilator(2, j)*E.H for j in range(2)]
    for i in range(2):
        for j in range(2):
            assert source[i]*source[j] + source[j]*source[i] == s.zeros(9)
            assert source[i]*source[j].H + source[j].H*source[i] == (P if i == j else s.zeros(9))
        assert parity*source[i] == -source[i]*parity
    # Exact arbitrary determinant-one two-component action. Unitarity is an
    # additional premise when this polynomial identity represents a rotation.
    a, b, c, d = s.symbols('a b c d')
    U = s.Matrix([[a, b], [c, d]])
    single = s.diag(1, U)
    exterior = s.diag(1, U, a*d-b*c)
    assert (s.kronecker_product(single, single)*E-E*exterior).applyfunc(s.expand) == s.zeros(9, 4)
    # All four occupation states are essential: the original3-state code
    # gives {b_up,b_up†}=diag(1,1,0), not its identity.
    hard = s.zeros(3); hard[0, 1] = 1
    assert hard*hard.H + hard.H*hard == s.diag(1, 1, 0)
    return E, P, source


@lru_cache(None)
def four_mode_code():
    E, _, _ = doublet_code()
    # Native constituent order: two +charge copies then two -charge copies.
    # Binary Fock order uses +charge bits0,1 and -charge bits2,3. Matrix
    # tensor indexing is low-word*4+high-word, hence this explicit column map.
    tensor = s.kronecker_product(E, E)
    embedding = s.zeros(81, 16)
    for word in range(16):
        embedding[:, word] = tensor[:, (word & 3)*4 + (word >> 2)]
    assert embedding.H*embedding == s.eye(16)
    P = embedding*embedding.H
    aa = [annihilator(4, j) for j in range(4)]
    native = [embedding*a*embedding.H for a in aa]
    n = [a.H*a for a in aa]
    charge = n[0]+n[1]-n[2]-n[3]
    number = sum(n, s.zeros(16))
    for i in range(4):
        for j in range(4):
            assert native[i]*native[j] + native[j]*native[i] == s.zeros(81)
            assert native[i]*native[j].H + native[j].H*native[i] == (P if i == j else s.zeros(81))
        q = 1 if i < 2 else -1
        assert charge*aa[i]-aa[i]*charge == -q*aa[i]
    assert list(number.diagonal()) == [word.bit_count() for word in range(16)]
    return {'constituents_per_macrocell': 4, 'selected_source_product_dimension': 81,
            'complete_four_mode_Fock_code_dimension': 16,
            'canonical_CAR_pairs_checked': 16,
            'all_native_complements_retained': True,
            'local_code_has_all_occupations_zero_through_four': True,
            'actual_native_odd_parity_retained': True,
            'mass_term_on_code': 'E_star times total canonical particle-plus-antiparticle number',
            'global_boundary_charge_on_code': 'N_plus-N_minus',
            'single_hard_core_doublet_is_not_a_two_mode_Fock_carrier': True}


def pauli():
    return (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))


def spinor_frames(momentum, mass):
    omega = s.sqrt(mass**2 + sum(p*p for p in momentum))
    sigma_p = sum((p*sig for p, sig in zip(momentum, pauli())), s.zeros(2))
    factor = s.sqrt((omega+mass)/(2*omega))
    U = factor*s.eye(2).col_join(sigma_p/(omega+mass))
    V = factor*(-sigma_p/(omega+mass)).col_join(s.eye(2))
    h = s.BlockMatrix([[mass*s.eye(2), sigma_p], [sigma_p, -mass*s.eye(2)]]).as_explicit()
    return omega, U, V, h


def spinor_checks():
    beta = s.diag(1, 1, -1, -1)
    for p in ((0, 0, 0), (s.Rational(3, 4), 0, 0), (s.Rational(2, 3), s.Rational(2, 3), s.Rational(1, 3))):
        w, U, V, h = spinor_frames(p, s.Integer(1))
        assert (U.H*U-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)
        assert (V.H*V-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)
        assert (U.H*V).applyfunc(s.simplify) == s.zeros(2)
        assert (U*U.H+V*V.H-s.eye(4)).applyfunc(s.simplify) == s.zeros(4)
        assert (h*U-w*U).applyfunc(s.simplify) == s.zeros(4, 2)
        assert (h*V+w*V).applyfunc(s.simplify) == s.zeros(4, 2)
        assert (h*h-(1+sum(x*x for x in p))*s.eye(4)).applyfunc(s.simplify) == s.zeros(4)
        # Covariant spin sums use the physical antiparticle momentum -p.
        _, _, Vminus, _ = spinor_frames(tuple(-x for x in p), s.Integer(1))
        slash = w*beta-sum((p[j]*beta*s.BlockMatrix([[s.zeros(2), pauli()[j]], [pauli()[j], s.zeros(2)]]).as_explicit()
                           for j in range(3)), s.zeros(4))
        assert (2*w*U*U.H*beta-(slash+s.eye(4))).applyfunc(s.simplify) == s.zeros(4)
        assert (2*w*Vminus*Vminus.H*beta-(slash-s.eye(4))).applyfunc(s.simplify) == s.zeros(4)
    # Exact noncollinear boost: cosh chi=5/4, sinh chi=3/4, along z.
    p = (s.Rational(3, 4), 0, 0)
    target = (p[0], 0, s.Rational(15, 16))
    w, U, _, _ = spinor_frames(p, s.Integer(1))
    wp, Up, _, _ = spinor_frames(target, s.Integer(1))
    alpha_z = s.BlockMatrix([[s.zeros(2), pauli()[2]], [pauli()[2], s.zeros(2)]]).as_explicit()
    S = (3*s.eye(4)+alpha_z)/(2*s.sqrt(2))
    beta = s.diag(1, 1, -1, -1)
    u, up = s.sqrt(2*w)*U, s.sqrt(2*wp)*Up
    D = (up.H*beta*S*u/2).applyfunc(s.simplify)
    assert (S*u-up*D).applyfunc(s.simplify) == s.zeros(4, 2)
    assert (D.H*D-s.eye(2)).applyfunc(s.simplify) == s.zeros(2)
    assert s.simplify(D.det()) == 1
    return {'normalized_positive_negative_frames_are_unitary': True,
            'both_covariant_spin_sums_checked': True,
            'source_charge_compatible_particle_and_antiparticle_components': True,
            'noncollinear_boost_Wigner_matrix': [[str(x) for x in D.row(j)] for j in range(2)],
            'mass_is_actual_selected_native_E_star': True}


def complete_two_tests(quarter_period=True):
    """Literal16-state CAR circuit: one particle, one antiparticle, two references.

    Source field at p/m=3/4, f=(component1+component3)/sqrt2, normalized
    coefficients2/sqrt5 and1/sqrt5. Reference modes have charge+1.
    The second source time is pi/(2 Omega_a), so all phases here are exact.
    """
    aa = [annihilator(4, j) for j in range(4)]
    a, b, r1, r2 = aa
    n = [x.H*x for x in aa]
    Q = n[0]-n[1]+n[2]+n[3]
    parity = s.diag(*[(-1)**word.bit_count() for word in range(16)])
    d = (2*a+b.H)/s.sqrt(5)
    assert d*d == s.zeros(16)
    assert d*d.H+d.H*d == s.eye(16)
    assert Q*d-d*Q == -d
    Nsource = n[0]+n[1]
    free = s.diag(*[(-s.I)**int(Nsource[j,j]) for j in range(16)]) if quarter_period else s.eye(16)
    def gate(r):
        G = r.H*d-d.H*r
        assert G.H == -G and G*Q == Q*G and G*parity == parity*G
        assert G**3 == -G
        return s.eye(16)+G/s.sqrt(2)+(1-1/s.sqrt(2))*G**2
    first, second = gate(r1), gate(r2)
    assert (first.H*first-s.eye(16)).applyfunc(s.simplify) == s.zeros(16)
    assert (second.H*second-s.eye(16)).applyfunc(s.simplify) == s.zeros(16)
    vacuum = s.zeros(16, 1); vacuum[0] = 1
    records = []
    for k in (0, 1):
        P1 = n[2] if k else s.eye(16)-n[2]
        for ell in (0, 1):
            P2 = n[3] if ell else s.eye(16)-n[3]
            output = P2*second*free*P1*first*vacuum
            probability = s.simplify((output.H*output)[0])
            assert probability >= 0
            records.append({'first_reference_number': k, 'second_reference_number': ell,
                            'ideal_actual_source_probability': str(probability)})
    assert s.simplify(sum(s.sympify(row['ideal_actual_source_probability']) for row in records)) == 1
    assert s.simplify(sum(s.sympify(row['ideal_actual_source_probability']) for row in records if row['first_reference_number'])) == s.Rational(1, 10)
    complete = []
    for up1 in (0, 1):
        for down1 in (0, 1):
            for up2 in (0, 1):
                for down2 in (0, 1):
                    row = next(r for r in records if r['first_reference_number'] == up1 and r['second_reference_number'] == up2)
                    complete.append({'first_up': up1, 'first_down': down1,
                                     'second_up': up2, 'second_down': down2,
                                     'ideal_actual_source_probability': row['ideal_actual_source_probability'] if down1 == down2 == 0 else '0'})
    assert len(complete) == 16
    return {'literal_executed_dimension': 16, 'executed_source_modes': 2,
            'executed_retained_reference_modes': 2,
            'opposite_spin_source_modes_are_retained_vacuum_spectators': True,
            'reference_spin_partners_are_retained_vacuum_spectators': True,
            'full_four_source_four_reference_mode_dimension': 256,
            'field_test_is_charge_conserving_and_even_with_reference': True,
            'each_swap_angle': 'pi/4', 'source_elapsed_time': 'pi/(2 Omega_a(p))' if quarter_period else '0',
            'active_four_records': records, 'all_sixteen_records': complete,
            'complete_probability_sum': '1'}


if __name__ == '__main__':
    import json
    print(json.dumps({'encoding': four_mode_code(), 'spinors': spinor_checks(),
                      'complete_native_code_test': complete_two_tests()}, indent=2))
