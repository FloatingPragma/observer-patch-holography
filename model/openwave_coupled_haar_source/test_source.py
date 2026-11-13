from fractions import Fraction as F
import sympy as sy
import pytest
from .source import admissible, labels, electric, wilson_entry, cg_integral_entry, matrices, hamiltonian, run


def test_exact_Haar_oracle_and_complete_graph():
    receipt = run()
    assert len(receipt['graph']['edges']) == 7
    assert len(receipt['graph']['vertices']) == 6
    assert receipt['full_successor_commutators_checked'] == 30
    assert receipt['diagnostic_caps'][-1]['complete_admissible_states'] == 285


def test_one_loop_reduction_and_shared_edge_interaction():
    for a in range(6):
        assert wilson_entry((a, 0, a), (a+1, 0, a+1), 0) == F(1, 2)
    # Identical loop spins carry different shared-edge electric energies.
    assert electric((1, 1, 0)) == F(9, 2)
    assert electric((1, 1, 2)) == F(13, 2)
    assert electric((1, 1, 0), False) == electric((1, 1, 2), False) == 6
    assert hamiltonian(2) != hamiltonian(2, shared=False)


def test_no_selected_J_or_invalid_spin_network():
    assert labels(1) == ((0,0,0),(0,1,1),(1,0,1),(1,1,0),(1,1,2))
    assert not admissible(1, 1, 1)
    with pytest.raises(ValueError):
        wilson_entry((1, 1, 1), (0, 0, 0), 0)
    with pytest.raises(ValueError):
        labels(True)
    for cap in (1, 2, 3):
        _, operators = matrices(cap)
        assert all(matrix == matrix.T for matrix in operators)
