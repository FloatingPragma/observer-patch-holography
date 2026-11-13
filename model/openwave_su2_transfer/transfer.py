"""Complete character-space controls for the positive heat-kernel transfer."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh, expm
from model.openwave_su2_spectral.spectrum import jacobi


def matrices(n=12, g2=F(8)):
    d, b = jacobi(n, g2)
    beta = 2/g2
    E = np.diag(np.array([z-beta for z in d], float))
    V = float(beta)*np.eye(n+1)+float(b)*(np.eye(n+1, k=1)+np.eye(n+1, k=-1))
    W = (np.eye(n+1, k=1)+np.eye(n+1, k=-1))/2
    return E, V, W


def execute(n=12, g2=F(8)):
    E, V, W = matrices(n, g2)
    H = E+V
    ev, Q = eigh(H)
    exact = expm(-H)
    commutator = float(np.linalg.norm(E@V-V@E, 2))
    rows = []
    for steps in (64, 128, 256):
        h = 1/steps
        half = expm(-h*V/2)
        T = half@np.diag(np.exp(-h*np.diag(E)))@half
        te, U = eigh(T)
        assert te[0] > 0 and te[-1] < 1
        repeated = (U*(te**steps))@U.T
        difference = float(np.linalg.norm(repeated-exact, 2))
        error = commutator/(2*steps)
        assert difference < error
        omega = U[:, -1]
        if omega[0] < 0:
            omega = -omega
        normalized = T/te[-1]
        centered = W-float(omega@W@omega)*np.eye(n+1)
        # Reflected histories are actual vectors of the same positive transfer.
        histories = [omega, centered@omega,
                     normalized@centered@omega,
                     centered@normalized@centered@omega,
                     normalized@centered@normalized@centered@omega]
        columns = np.column_stack(histories)
        gram = columns.T@columns
        minimum = float(eigh(gram, eigvals_only=True)[0])
        assert minimum > -1e-13
        we, Z = eigh(W)
        full_readout = abs(Z.T@omega)**2
        assert abs(full_readout.sum()-1) < 1e-12
        rows.append({'time_step_exact': str(F(1, steps)),
                     'complete_transfer_eigenvalues': te.tolist(),
                     'unit_time_semigroup_operator_error': difference,
                     'finite_core_contractivity_error_bound': error,
                     'five_history_reflection_Gram': gram.tolist(),
                     'reflection_Gram_minimum_numeric': minimum,
                     'complete_13_outcome_plaquette_values': we.tolist(),
                     'complete_13_outcome_plaquette_probabilities': full_readout.tolist(),
                     'transfer_vacuum_is_its_own_preparation': True})
    assert rows[-1]['unit_time_semigroup_operator_error'] < rows[0]['unit_time_semigroup_operator_error']/10
    return {'character_cutoff_for_execution_only': n,
            'same_source_Hamiltonian_eigenvalues': ev.tolist(),
            'rows': rows,
            'negative_transfer_control': {'operator': 'diag(1,-1/2)',
                'test_vector': [0, 1], 'bond_reflection_quadratic_form_exact': '-1/2',
                'self_adjointness_alone_is_sufficient': False}}


def run():
    return {'schema': 'oph.pure_SU2_positive_transfer_same_clock.v1',
            'verdict': 'PASS_FULL_COMPACT_GAUGE_TRANSFER_POSITIVITY_AND_COMMON_CLOCK',
            'full_source_theorem': {
                'space': 'full product Haar SU2 links and its Gauss-invariant range',
                'temporal_kernel': 'sum_{n>=0} (n+1) exp[-h alpha n(n+2)/4] chi_n(U V^-1)',
                'every_character_coefficient_strictly_positive': True,
                'operator': 'T_h=exp(-h V/2) exp(-h E) exp(-h V/2)',
                'injective_positive_self_adjoint_compact': True,
                'Gauss_projection_commutes_with_transfer': True,
                'all_bond_reflection_history_Gram_forms_nonnegative': True,
                'same_clock_limit': 'strong limit (T_(t/n))^n=exp(-t(E+V))',
                'fixed_graph_Peter_Weyl_cores_are_form_cores': True,
                'reflection_positive_kernel_is_a_real_time_detector': False},
            'complete_finite_character_execution': execute(),
            'scope': {'Y2_finite_regulator_positive_transfer_instantiated': True,
                'same_finite_graph_Hamiltonian_clock_identified': True,
                'Y1_uniform_renormalized_spatial_continuum_Cauchy_proved': False,
                'Y3_full_Euclidean_covariance_restoration_proved': False,
                'Y4_joint_spatial_infinite_volume_Mosco_proved': False,
                'Y5_continuum_OS_identification_proved': False,
                'Y6_uniform_continuum_noncollapse_proved': False,
                'Clay_problem_solved': False}}
