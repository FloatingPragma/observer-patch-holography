"""Full-carrier rank, spectator transitions and scope-mutation controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from ..coarse_energy_ward import ward as w
from ..field_symmetry.symmetry import operators
from ..full_fock_field.field import FieldGeometry
from . import algebra as a
from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet.produce()


def test_complete_cube_ranks_and_independent_full_fock_entries(receipt):
    results = receipt['certificate']
    cube = [row for row in results['classifications'] if row['sides'] == [2, 2, 2]]
    assert [(r['quadratic_dimension_excluding_identity'],
             r['rational_and_real_commutator_rank'],
             r['independent_exact_kernel_dimension']) for r in cube] == [(496, 492, 4), (496, 490, 6)]
    assert [r['exact_rank'] for r in results['independent_full_fock']] == [116, 114]
    assert all(r['all_input_output_matrix_entries'] == 256**2
               for r in results['independent_full_fock'])
    assert results['spectator_witnesses']['isolated_complex_offsite_pair_coefficients'] == 112


@pytest.mark.parametrize('sides', [(2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1)])
@pytest.mark.parametrize('centered', [False, True])
def test_source_matches_independent_integer_fock_hamiltonian(sides, centered):
    geometry = FieldGeometry(sides)
    exact = operators(geometry)
    actual = (exact.K if centered else exact.H).matrix(exact.scale).toarray()
    assert np.array_equal(a.source(geometry, centered=centered).matrix(geometry.modes), actual)


@pytest.mark.parametrize('speed,interaction', [(F(1, 4), F(3, 2)), (F(3, 8), F(1, 16))])
@pytest.mark.parametrize('centered', [False, True])
def test_no_special_coupling_assumption_and_rational_rank(speed, interaction, centered):
    geometry = FieldGeometry((2, 1, 1), speed=float(speed), interaction=float(interaction))
    h = a.source(geometry, centered=centered)
    rank = a.rational_rank(tuple(w.commutator(h, q) for q in a.quadratic_basis(4)))
    result = a.classification(geometry, centered=centered)
    assert result['rational_and_real_commutator_rank'] == rank == (22 if centered else 24)


def test_spectator_transition_detects_pair_invisible_on_vacuum():
    geometry = FieldGeometry((2, 1, 1), speed=0.)
    contact = a.source(geometry, centered=False)
    pair = w.annihilator(0).adjoint()@w.annihilator(2).adjoint()
    residual = w.commutator(contact, pair).matrix(4)
    assert np.array_equal(residual[:, 0], np.zeros(16))
    assert abs(residual[0b111, 0b010]) == 1
    assert a.spectator_witnesses(3)['isolated_complex_offsite_pair_coefficients'] == 12


def test_commuting_internal_algebras_and_scalar_dimension():
    geometry = FieldGeometry((2, 2, 1))
    nz, sx, sy, sz, ex, ey = a.generators(geometry, centered=True)
    assert w.commutator(sx, sy) == sz.scale((0, 2))
    assert w.commutator(ex, ey) == nz.scale((0, 2))
    for spin in (sx, sy, sz):
        for eta in (nz, ex, ey):
            assert w.commutator(spin, eta).is_zero()
    h = a.source(geometry, centered=True)
    assert w.commutator(h, sx+w.Clifford.scalar(F(5, 7))).is_zero()
    assert a.modular_rank(a.generators(geometry, centered=True)+(w.Clifford.scalar(1),)) == 7


def test_hypothesis_countercontrols(receipt):
    assert receipt['certificate']['boundary_controls'] == {
        'U_zero_centered_dimer_quadratic_kernel_dimension': 16,
        'c_zero_centered_two_site_kernel_dimension': 12,
        'detuned_centered_dimer_quadratic_kernel_dimension': 4,
        'separate_real_hopping_triangle_kernel_dimension': 4,
        'proper_site_eta_pair_conserved': False,
    }


def test_source_mutation_with_same_nullspace_is_rejected(monkeypatch):
    original = FieldGeometry.groups.fget
    monkeypatch.setattr(FieldGeometry, 'groups', property(
        lambda geometry: tuple(tuple((p, q, 2*z) for p, q, z in group)
                               for group in original(geometry))))
    with pytest.raises(ArithmeticError, match='actual source'):
        a.classification(FieldGeometry((2, 1, 1)))


def test_pairing_generator_omission_is_rejected(monkeypatch):
    original = a.generators
    monkeypatch.setattr(a, 'generators',
                        lambda geometry, centered: original(geometry, centered=centered)[:-1])
    with pytest.raises(ArithmeticError, match='exhaust'):
        a.classification(FieldGeometry((2, 1, 1)))


def test_nonconserved_local_pair_cannot_replace_global_eta(monkeypatch):
    original = a.generators
    pair = w.annihilator(0).adjoint()@w.annihilator(1).adjoint()
    monkeypatch.setattr(a, 'generators', lambda geometry, centered:
                        original(geometry, centered=centered)[:-1]+(pair+pair.adjoint(),))
    with pytest.raises(ArithmeticError, match='conservation'):
        a.classification(FieldGeometry((2, 1, 1)))


@pytest.mark.parametrize('field,value', [
    ('integrated_gate_accepted', True), ('physical_comparison_count_added', False),
    ('all_pass', 1), ('canonical_physical_premises_discharged', ['PR-24']),
    ('source_law_changed', True),
])
def test_receipt_scope_changes_fail(receipt, field, value):
    altered = deepcopy(receipt)
    altered[field] = value
    with pytest.raises(ValueError, match='reconstruction'):
        packet.validate(altered, receipt)


def test_nested_claim_and_manifest_mutations_fail(receipt):
    for mutate in (
        lambda r: r['certificate']['classifications'][-1].__setitem__('rational_and_real_commutator_rank', 491),
        lambda r: r['certificate'].__setitem__('physical_species_identified', True),
        lambda r: r['certificate'].__setitem__('anomalous_search_is_complete_on_all_actions', True),
        lambda r: r['source_bindings'].pop('model/quadratic_conservation/PROOF.md'),
    ):
        altered = deepcopy(receipt)
        mutate(altered)
        with pytest.raises(ValueError, match='reconstruction'):
            packet.validate(altered, receipt)


def test_retained_receipt_and_string_path_cold_replay():
    assert packet.verify(str(packet.OUTPUT)) is True


def test_duplicate_json_keys_fail_before_replay(tmp_path):
    path = tmp_path/'duplicate.json'
    path.write_text('{"all_pass":true,"all_pass":false}', encoding='utf-8')
    with pytest.raises(ValueError):
        packet.verify(path)


@pytest.mark.parametrize('geometry', [
    FieldGeometry((2, 1, 1), speed=0.), FieldGeometry((2, 1, 1), interaction=0.),
    FieldGeometry((3, 3, 1)), FieldGeometry((2, 1, 1), speed=0.1),
])
def test_classification_domain_is_explicit(geometry):
    with pytest.raises(ValueError):
        a.classification(geometry)


def test_prime_and_denominator_checks(monkeypatch):
    with pytest.raises(ValueError, match='singular'):
        a.modular_rank((w.Clifford.scalar(F(1, a.PRIME)),))
    monkeypatch.setattr(a, 'PRIME', 15)
    with pytest.raises(ArithmeticError, match='prime'):
        a.modular_rank((w.Clifford.scalar(1),))
