"""Actual binary cubic subgroup and retained regular-reference controls."""
from functools import lru_cache
import sympy as s

from .source import ROOT
from .algebra import doublet_code, pauli
from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import cubic_rotations


def clean(matrix):
    return matrix.applyfunc(lambda value: s.expand(s.simplify(value)))


@lru_cache(None)
def spatial_group():
    lifts = [(R, sign*S) for R, S in cubic_rotations() for sign in (1, -1)]
    keys = {tuple(clean(S)): i for i, (_, S) in enumerate(lifts)}
    assert len(keys) == 48
    generators = [(s.eye(2)-s.I*pauli()[axis])/s.sqrt(2) for axis in (0, 2)]
    rotations = [s.Matrix(3, 3, lambda j, i: s.simplify(s.trace(pauli()[j]*G*pauli()[i]*G.H)/2))
                 for G in generators]
    adjacency = []
    for R, S in lifts:
        row = []
        assert clean(S.H*S) == s.eye(2) and s.simplify(S.det()) == 1
        for G, rotation in zip(generators, rotations):
            index = keys[tuple(clean(S*G))]
            assert lifts[index][0] == R*rotation
            row.append(index)
        adjacency.append(row)
    reached = {keys[tuple(s.eye(2))]}; pending = list(reached)
    while pending:
        for target in adjacency[pending.pop()]:
            if target not in reached:
                reached.add(target); pending.append(target)
    assert len(reached) == 48
    E, _, _ = doublet_code()
    for _, S in lifts:
        assert clean(s.kronecker_product(s.diag(1, S), s.diag(1, S))*E
                     - E*s.diag(1, S, 1)) == s.zeros(9, 4)
    return {'actual_native_spin_lifts': 48, 'proper_signed_permutation_rotations': 24,
            'exact_generator_products_checked': 96,
            'all_two_copy_spinor_and_occupied_singlet_intertwiners_checked': 48,
            'central_element_is_actual_fermion_parity': True,
            'periodic_macro_group': '(Z/NZ)^3 semidirect 2O',
            'regular_reference_dimension': '48 N^3',
            'macro_translations_are_declared_periodic_not_native_internal_AP_translations': True,
            'native_internal_translation_covariance_claimed': False}


def complete_lift_control():
    """Exact all-Kraus identity for a noninvariant two-record instrument."""
    I = s.eye(2); X = pauli()[0]; Z = pauli()[2]
    source_group = [I, X]
    J = s.Matrix.vstack(I, X)/s.sqrt(2)
    W = s.diag(1, s.I)
    projectors = [s.diag(1, 0), s.diag(0, 1)]
    Wlift = s.diag(*[U*W*U.H for U in source_group])
    Hlift = s.diag(X, X)
    assert J.H*J == I and Hlift*J == J*X
    assert Wlift*J == J*W
    total = s.zeros(4)
    for P in projectors:
        Plift = s.diag(*[U*P*U.H for U in source_group])
        K = P*W; Klift = Plift*Wlift
        assert Klift*J == J*K
        total += Klift.H*Klift
    assert total == s.eye(4)
    # The reference is retained: its source-only marginal is the group twirl.
    rho = s.diag(1, 0)
    assert (rho+X*rho*X)/2 == I/2
    return {'noninvariant_instrument_all_Kraus_intertwiners_exact': True,
            'complete_lift_effect_sum_is_full_identity': True,
            'source_only_marginal_is_twirl': True,
            'joint_records_are_reference_relative': True,
            'reference_Hamiltonian': '0'}
