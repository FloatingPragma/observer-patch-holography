"""Independent exterior and native-branch oracles; capacity counterexamples."""
from itertools import product
import json
from pathlib import Path

import pytest
import sympy as sp

from . import identity as d


@pytest.fixture(scope='module')
def receipt():
    return d.certificate()


def test_existing_car_representation_against_independent_minor_derivative():
    # det((I+tA)[target,source]) has derivative dGamma(A). This oracle
    # uses determinants, not the creation/annihilation implementation.
    t = sp.symbols('t')
    basis = d.matter_basis()
    for one_body, actual in d.generators():
        transform = sp.eye(5)+t*one_body
        for target, a in enumerate(basis):
            for source, b in enumerate(basis):
                expected = 0
                if len(a) == len(b):
                    expected = sp.expand(transform.extract(a, b).det()).coeff(t, 1)
                assert actual[target, source] == expected


def test_all_unitary_involution_spectra_exhausted_by_six_classes():
    # Every unitary involution is unitarily diagonalizable with +/-1
    # eigenvalues. Exhaust all determinant-one five-mode sign assignments,
    # including ones not chosen as the canonical representative.
    classes = {(row['negative_color'], row['negative_weak']): row for row in d.involution_classes()}
    assert len(classes) == 6
    for signs in product((-1, 1), repeat=5):
        if sp.prod(signs) != 1:
            continue
        a, b = signs[:3].count(-1), signs[3:].count(-1)
        row = classes[a, b]
        t = sp.symbols('t')
        generating = sp.Poly(sp.prod(1+t*s for s in signs), t)
        trace = generating.nth(2)+generating.nth(4)
        assert row['plus']-row['minus'] == trace
        assert row['plus']+row['minus'] == 15
        assert (row['plus'], row['minus']) == ((15, 0) if (a, b) == (0, 0) else (7, 8))


@pytest.mark.parametrize('a,b', [(1, 0), (0, 1), (3, 2), (-1, 1), (4, 0), (True, 1)])
def test_non_group_or_nonexact_involution_rejected(a, b):
    with pytest.raises(ValueError):
        d.involution(a, b)


def test_only_weak_parity_supports_full_commuting_gauge_action():
    compatible = []
    for row in d.involution_classes():
        J = sp.diag(*row['matter_signs'])
        commutes = all(J*g == g*J for _, g in d.generators())
        assert commutes == row['central']
        if commutes and row['minus']:
            compatible.append((row['negative_color'], row['negative_weak']))
    assert compatible == [(0, 2)]
    blocks = d.involution(0, 2)['blocks']
    assert blocks == {'Q': {'plus': 0, 'minus': 6}, 'u_c': {'plus': 3, 'minus': 0},
                      'e_c': {'plus': 1, 'minus': 0}, 'd_c': {'plus': 3, 'minus': 0},
                      'L': {'plus': 0, 'minus': 2}}


def test_global_faithfulness_uses_generic_polynomial_not_lie_rank():
    value = d.global_faithfulness_certificate()
    assert value['generic_matrix_variables'] == 25
    assert value['exact_polynomial_entries_checked'] == 25
    H = sp.Matrix(value['H_volume_pairing'])
    assert H*H.T == sp.eye(5)
    # The nontrivial SU(5) center has zero Lie-level diagnostic but is not
    # in the group kernel: Lambda^4(zI)=z^4 I != I for primitive fifth z.
    z = sp.symbols('z')
    assert sp.rem(z**4-1, sum(z**k for k in range(5)), domain=sp.QQ) != 0


@pytest.mark.parametrize('multiplicity', range(1, 8))
def test_exact_weak_parity_capacity_obstruction(multiplicity):
    with pytest.raises(ValueError, match='insufficient capacity'):
        d.embedding(d.involution(0, 2)['matter_signs'], multiplicity)


def test_capacity_is_not_a_choice_of_actual_memory_size():
    signs = d.involution(0, 2)['matter_signs']
    for multiplicity in (8, 9, 17):
        d.validate_embedding(d.embedding(signs, multiplicity), signs, multiplicity)
    # A trivial gauge image is a different hypothesis; it needs fifteen
    # even slots, and does not identify the nontrivial deck character.
    with pytest.raises(ValueError, match='insufficient capacity'):
        d.embedding([1]*15, 8)
    d.validate_embedding(d.embedding([1]*15, 15), [1]*15, 15)


def test_loss_of_norm_or_wrong_charge_rejected():
    signs = d.involution(0, 2)['matter_signs']
    matrix = sp.MutableDenseMatrix(d.embedding(signs, 8))
    matrix[:, 1] = matrix[:, 0]
    with pytest.raises(ValueError, match='loses norm'):
        d.validate_embedding(matrix, signs, 8)
    matrix = d.embedding(signs, 8)
    wrong = signs.copy()
    wrong[0] *= -1
    with pytest.raises(ValueError, match='intertwine'):
        d.validate_embedding(matrix, wrong, 8)


def test_noncentral_signature_match_is_not_full_gauge_compatibility():
    row = d.involution(1, 1)
    matrix = d.embedding(row['matter_signs'], 8)
    d.validate_embedding(matrix, row['matter_signs'], 8)
    Q = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(8))
    embedded = [matrix*g*matrix.H/2 for _, g in d.generators()]
    assert any(Q*g != g*Q for g in embedded)


def test_weak_grading_cannot_give_fifteen_uniform_negative_labels():
    weak = d.involution(0, 2)['matter_signs']
    W = d.embedding(weak, 8)
    Q = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(8))
    assert Q*W != -W
    assert weak.count(1) == 7 and weak.count(-1) == 8
    # Global faithfulness promotes rho(g)^2=I to g^2=I; these are all
    # involution classes, none with the required all-negative signature.
    assert all(row['minus'] != 15 for row in d.involution_classes())


@pytest.mark.parametrize('multiplicity', range(1, 15))
def test_uniform_negative_collective_capacity_obstruction(multiplicity):
    with pytest.raises(ValueError, match='insufficient capacity'):
        d.embedding([-1]*15, multiplicity)


def test_uniform_negative_initialization_keeps_full_positive_sector():
    W = d.embedding([-1]*15, 15)
    assert W == sp.Matrix.vstack(sp.eye(15), -sp.eye(15))
    H = sp.Matrix([[1, 1], [1, -1]])/sp.sqrt(2)
    Z = sp.diag(1, -1)
    initial = sp.Matrix.vstack(sp.eye(15), sp.zeros(15))
    preparation = sp.kronecker_product(Z*H, sp.eye(15))
    assert preparation.H*preparation == sp.eye(30)
    assert preparation*initial == W/sp.sqrt(2)
    negative = W*W.H/2
    positive = sp.eye(30)-negative
    assert positive.rank() == 15
    Q = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(15))
    assert Q*negative == -negative and Q*positive == positive
    for _, g in d.generators():
        full = sp.kronecker_product(sp.eye(2), g)
        assert full*W == W*g
        assert full*positive == positive*full


def test_group_extension_not_only_lie_algebra():
    signs = d.involution(0, 2)['matter_signs']
    W = d.embedding(signs, 8)
    P = W*W.H/2
    # Real rotations in each block have exact rational entries. Verify the
    # full group law from exterior minors, independently of dGamma.
    color = sp.eye(5)
    color[0, 0], color[1, 1] = sp.Rational(3, 5), sp.Rational(3, 5)
    color[0, 1], color[1, 0] = sp.Rational(4, 5), sp.Rational(-4, 5)
    weak = sp.eye(5)
    weak[3, 3], weak[4, 4] = sp.Rational(5, 13), sp.Rational(5, 13)
    weak[3, 4], weak[4, 3] = sp.Rational(12, 13), sp.Rational(-12, 13)
    def exterior(g):
        basis = d.matter_basis()
        return sp.Matrix([[g.extract(a, b).det() if len(a) == len(b) else 0
                           for b in basis] for a in basis])
    def extension(g):
        return W*exterior(g)*W.H/2+sp.eye(16)-P
    assert extension(color)*extension(weak) == extension(color*weak)
    assert extension(color).H*extension(color) == sp.eye(16)
    assert extension(weak).H*extension(weak) == sp.eye(16)
    assert extension(sp.diag(1, 1, 1, -1, -1)) == sp.kronecker_product(
        sp.Matrix([[0, 1], [1, 0]]), sp.eye(8))


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_intertwiner_against_named_token_oracle(law):
    # Check the source transport without using row.signs or lifted_pairs.
    from ..exchange_sector import statistics as exchange
    signs = d.involution(0, 2)['matter_signs']
    W = d.embedding(signs, 8)
    # The two exact memory identities apply to every column below; compute
    # them once so this exhaustive transport oracle stays inexpensive.
    images = [W, sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(8))*W]
    assert images == [W, W*sp.diag(*signs)]
    for row in d.protection.branches(law):
        a, b = row.seam
        for i, pair in enumerate(exchange.PAIRS):
            n = int(a in pair)+int(b in pair)
            current = int(a in pair)-(n//2+row.coin*(n % 2))
            raw = {0: 0, 1: 1, -1: 2}[current] + (4*(n % 2) if law == 'E' else 0)
            if raw != row.raw:
                assert i not in row.inputs
                continue
            tokens = list(pair)
            if current:
                origin, target = (a, b) if current == 1 else (b, a)
                tokens[tokens.index(origin)] = target
            final_pair = tuple(sorted(tokens))
            parity = int(tuple(tokens) != final_pair)
            at = row.inputs.index(i)
            assert row.outputs[at] == exchange.INDEX[final_pair]
            assert row.signs[at] == (-1)**parity


def test_native_branch_mutation_is_not_hidden_by_signature(monkeypatch):
    original = d.protection.lifted_pairs
    def erase_memory_transport(row):
        return tuple((r*66+i, r*66+j) for i, j in zip(row.inputs, row.outputs) for r in (0, 1))
    monkeypatch.setattr(d.protection, 'lifted_pairs', erase_memory_transport)
    signs = d.involution(0, 2)['matter_signs']
    with pytest.raises(ArithmeticError, match='intertwiner failed'):
        d.native_intertwining('C', d.embedding(signs, 8), signs, 8)
    monkeypatch.setattr(d.protection, 'lifted_pairs', original)


def test_receipt_scope_and_replay(receipt):
    assert receipt == json.loads(Path(d.__file__).with_name('results.json').read_text())
    assert receipt['physical_identification'] is False
    assert 'explicit permitted design' in receipt['spectator_extension']['fiber_and_gauge_action_status']
    assert receipt['spectator_extension']['minimal_feasible_m_for_this_declared_extension'] == 8
    for row in receipt['complete_branch_checks'].values():
        assert row['complete_native_branch_intertwining']
        assert row['all_source_and_matter_columns'] == row['raw_slots']*66*15
        assert row['abort_slots_retained'] > 0
    fermions = receipt['uniform_negative_character_extension']
    assert fermions['chosen_m'] == 15
    assert fermions['matter_exchange_signature'] == [0, 15]
    assert fermions['native_extended_dimension'] == 1980
    assert receipt['internal_gauge_uniform_fermion_parity']['exists'] is False
    assert 'not a uniformly negative' in receipt['spectator_extension']['scope']
    for law, row in fermions['complete_branch_checks'].items():
        assert row == receipt['complete_branch_checks'][law]
