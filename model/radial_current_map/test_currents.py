"""Independent full-support current, multiplet and leakage controls."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from ..native_mass_channels.algebra import SIGMA
from ..radial_defect_gap.finite import kernel_frame
from ..radial_interacting_core.algebra import source_matrix
from . import algebra, controls, representation


def test_full_even_algebra_has_quartic_and_pairing_directions():
    row = algebra.local_completeness()
    assert row['all_even_Hermitian_dimension'] == 8
    assert row['quartic_local_term_retained']
    assert row['all_eta_Nambu_ladder_equations_exact']
    c = algebra.annihilators(2)
    q = algebra.onsite_basis(c, 0)
    # The quartic direction cannot be represented by a scalar and quadratic currents.
    assert s.Matrix.hstack(*(s.Matrix(x).reshape(16, 1) for x in q[:-1])).rank() == 7
    assert q[7] not in q[:7]


@pytest.mark.parametrize('axis', range(3))
def test_each_full_edge_kernel_has_no_quartic_escape(axis):
    row = algebra.edge_classification(axis)
    assert row['equation_rank'] == row['kernel_dimension'] == 8
    assert all(v[7] == v[15] == 0 for v in row['kernel_columns'])
    c = algebra.annihilators(4)
    z0, z1 = algebra.onsite_basis(c, 0)[-1], algebra.onsite_basis(c, 1)[-1]
    hopping = sum((s.I*SIGMA[axis][i, j]*c[i].H*c[2+j]/2
                   for i in range(2) for j in range(2)), s.zeros(16))
    hopping += hopping.H
    assert (z0+z1)*hopping != hopping*(z0+z1)
    assert (z0-z1)*hopping != hopping*(z0-z1)


def test_transport_is_path_independent_and_literal_mass_kills_every_spin_axis():
    h, points = source_matrix()
    mass = s.diag(*(h[2*x:2*x+2, 2*x:2*x+2] for x in range(len(points))))
    hopping = h-mass
    for p in SIGMA:
        blocks = []
        for point in points:
            A = p
            for axis, n in enumerate(point):
                if n % 2: A = SIGMA[axis]*A*SIGMA[axis]
            blocks.append(A)
        current = s.diag(*blocks)
        assert current*hopping == hopping*current
        assert current*mass != mass*current
    for i, j in product(range(3), repeat=2):
        for A in SIGMA:
            assert SIGMA[i]*SIGMA[j]*A*SIGMA[j]*SIGMA[i] == SIGMA[j]*SIGMA[i]*A*SIGMA[i]*SIGMA[j]


def test_all_source_rows_obey_pairing_and_scalar_mass_breaks_it():
    assert algebra.literal_source_pairing()['h Delta + Delta h^T_exact_zero']
    h, points = source_matrix()
    delta = s.kronecker_product(s.diag(*((-1)**(sum(x) % 2) for x in points)), s.I*SIGMA[1])
    changed = h+s.eye(54)
    assert changed*delta+delta*changed.T == 2*delta


def test_ground_and_charged_current_representations_are_distinct():
    row = representation.low_representation()
    assert row['eta_Casimir'] == ['3/4', '0', '0', '3/4']
    assert row['eta_current_leakage_from_complete_low_quartet'] == '0'
    assert row['charged_ladder_norm_squared'] == '1'
    assert row['particle_numbers'] == [26, 27, 27, 28]
    assert not row['charge_doublet_is_a_selected_Standard_Model_multiplet']


def test_exact_actual_interacting_leakage_strictly_exceeds_one_over_4096():
    row = representation.leakage_bounds()
    alpha = F(2904, 2905); t = F(27, 16252874)
    bound = F(row['each_actual_charged_center_Pauli_leakage_squared_lower'])
    assert bound == (1-alpha*alpha)/2-t > F(1, 4096)
    assert row['actual_charged_leakage_is_outside_entire_low_quartet']


def test_literal_zero_frame_has_nonzero_compression_bracket_defect_and_leakage():
    value = kernel_frame(1, mass=64)
    W = value['frame']; P = W@W.conj().T; origin = value['P_rows']
    alpha = value['density'][13]
    compressed = []
    for p in SIGMA:
        S = np.zeros((54, 54), complex); S[np.ix_(origin, origin)] = np.array(p, complex)
        small = W.conj().T@S@W; compressed.append(small)
        leakage = (np.eye(54)-P)@S@W
        assert np.linalg.norm(small-alpha*np.array(p, complex)) < 1e-11
        assert np.linalg.norm(leakage.conj().T@leakage-alpha*(1-alpha)*np.eye(2)) < 1e-11
    assert np.linalg.norm(compressed[0]@compressed[1]-compressed[1]@compressed[0]-2j*compressed[2]) > 1e-5


@pytest.mark.parametrize('mass,coupling', [(4, F(1, 256)), (64, F(1, 65536))])
def test_full_interacting_Fock_controls_and_wrong_law_countercontrols(mass, coupling):
    row = controls.chain_control(mass, coupling)
    assert row['all_Fock_states'] == 64 and row['all_number_sectors_checked'] == 7
    assert row['all_three_center_spin_currents_leak_outside_actual_quartet']
    assert row['wrong_unstaggered_pairing_rejected']
    assert row['uncentered_contact_has_nonzero_pair_commutator']
    assert row['all_sites_and_ground_columns_addition_removal_amplitudes_intertwine']


def test_actual_native_target_has_nontrivial_color_and_weak_actions():
    row = representation.native_current_comparison()
    assert (row['native_ports'], row['native_color_current_rank'], row['native_weak_current_rank']) == (12, 8, 3)
    assert sum(v['dimension'] for v in row['target_multiplets_from_source_Casimirs']) == 15
    assert not row['faithful_full_native_current_action_in_source_local_conserved_algebra']
    assert not row['nonlocal_or_emergent_current_algebras_excluded_by_this_theorem']
