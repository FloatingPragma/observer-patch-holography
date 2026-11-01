"""Independent exact minor, eta intertwiner and actual repair checks."""
from itertools import combinations
import math

import numpy as np
import pytest
from scipy import sparse

from ..charged_defect import dictionary
from ..common_spin_protection import spectrum, spin
from . import algebra as a


def determinant_mod(matrix, prime):
    """Row-elimination determinant, independent of streaming column pivots."""
    value = [[int(x) % prime for x in row] for row in matrix]
    result = 1
    for j in range(len(value)):
        pivot = next((i for i in range(j, len(value)) if value[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            value[pivot], value[j] = value[j], value[pivot]; result = -result
        coefficient = value[j][j]; result = result*coefficient % prime
        inverse = pow(coefficient, -1, prime)
        for i in range(j+1, len(value)):
            factor = value[i][j]*inverse % prime
            for k in range(j, len(value)):
                value[i][k] = (value[i][k]-factor*value[j][k]) % prime
    return result % prime


def test_exact_ninety_by_ninety_minor_and_characteristic_zero_upper_bound():
    proof = a.cyclic_certificate()
    h, raising, lowering, c, d = a.exact_data()
    values = (c+44*d) % 101
    columns = []
    for _ in range(43):
        columns.append(values); values = h@values % 101
    krylov = np.column_stack(columns)
    rows, cols, _ = zip(*proof['pivot_witness'], strict=True)
    assert len(rows) == len(set(rows)) == len(set(cols)) == 90
    assert determinant_mod(krylov[np.ix_(rows, cols)], 101) != 0
    assert np.count_nonzero(raising@c) == np.count_nonzero(raising@d) == 0
    assert np.count_nonzero(lowering@c) == np.count_nonzero(lowering@d) == 0
    # Independently count the exact joint multiplicity from binomial weights.
    def w(s, e):
        up, down = (6+e+s)//2, (6+e-s)//2
        return math.comb(6, up)*math.comb(6, down)
    assert w(2, 2)-w(4, 2)-w(2, 4)+w(4, 4) == 90
    # The second distinct root is a useful arithmetic control, not an upper
    # bound inferred from rank deficiency.
    values = (c+57*d) % 101
    columns = []
    for _ in range(43):
        columns.append(values); values = h@values % 101
    assert a.modular_rank(np.column_stack(columns))[0] == 90
    normalized = (c+np.sqrt(17)*d)/np.sqrt(272-64*np.sqrt(17))
    actual = spin.operators().frame.matrix()@dictionary.embedded_columns(3, 'hole')[:, [0, 4, 8]]
    np.testing.assert_allclose(actual[a.WHERE], normalized*np.array([-1, 1, -1]), atol=2e-15, rtol=0)
    assert proof['actual_repair_seed_binding_residual'] < 2e-12


def test_decoupling_removes_active_cyclicity_and_no_larger_rank_is_invented():
    h, _, _, c, d = a.exact_data()
    frame = spin.operators().frame.matrix()
    decoupled = frame@dictionary.tensor_generator(3)@frame.conj().T
    numeric = (16*decoupled[a.WHERE][:, a.WHERE]).toarray()
    np.testing.assert_array_equal(numeric.imag, 0)
    np.testing.assert_array_equal(numeric.real, np.rint(numeric.real))
    decoupled = numeric.real.astype(np.int64)
    values = (c+44*d) % 101
    assert a.modular_rank(np.column_stack([values, decoupled@values % 101]))[0] == 3
    assert np.any(h != decoupled)


def test_dark_fifteen_orthonormal_columns_and_every_complete_repair_outcome():
    numerator, _ = a.dark_columns()
    d = numerator/2
    frame = spin.operators().frame.adjoint().matrix()
    d = frame@d
    # Complete all three spin weights by the actual hidden lowering operator.
    spinops = spin.operators().spin_twice
    lower = (spinops[0].matrix()-1j*spinops[1].matrix())/2
    columns = [d, lower@d/np.sqrt(2), lower@lower@d/2]
    complete = sparse.hstack(columns, format='csr')
    assert sparse.linalg.norm(complete.conj().T@complete-sparse.eye(45)) < 1e-14
    for law in ('C', 'E'):
        for seam in combinations(range(3), 2):
            for label, op in spin.repair_outcomes('hole', law, seam):
                expected = complete if label == ('outside',) else sparse.csr_matrix(complete.shape)
                assert sparse.linalg.norm(op@complete-expected) < 1e-13
    proof = a.dark_certificate()
    assert proof['distinct_energies'] == 13
    assert proof['dark_joint_commutant_complex_dimension'] == 12*3**2+9**2


def test_dark_energies_from_independent_slater_sine_functions():
    d, pairs = a.dark_columns()
    h = (spin.operators().frame@spin.operators().K@spin.operators().frame.adjoint()).matrix(16)
    for p, q in combinations(range(1, 7), 2):
        first = [np.sqrt(2/7)*np.sin(p*np.pi*(x+1)/7) for x in range(6)]
        second = [np.sqrt(2/7)*np.sin(q*np.pi*(x+1)/7) for x in range(6)]
        state = d@np.array([first[x]*second[y]-first[y]*second[x] for x, y in pairs])/2
        eigenvalue = (np.cos(p*np.pi/7)+np.cos(q*np.pi/7))/8-1
        assert abs(np.vdot(state, state)-1) < 2e-15
        assert np.linalg.norm(h@state-eigenvalue*state) < 2e-15


def test_actual_complete_repair_generates_the_source_bound_three_state_corner():
    report = a.repair_corner_certificate()
    assert report['raw_slots_including_zero'] == 72
    assert len(set(map(tuple, report['effect_signatures']))) == 3
    assert report['corner_star_algebra_dimension'] == 9
    assert report['full_Fock_dictionary_max_residual'] < 2e-12


@pytest.mark.parametrize('value', [np.zeros((3, 3)), [1, 2], np.array([['x']])])
def test_invalid_exact_rank_matrix_rejects(value):
    with pytest.raises(ValueError):
        a.modular_rank(value)
