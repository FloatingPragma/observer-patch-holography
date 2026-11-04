"""Native-current derivation and independent invariant-tensor controls."""
import sympy as s

from ..common_source import calibrated_matter as calibrated
from ..common_source.chiral_candidate import action
from . import native


def test_fresh_native_LH_dictionary_and_six_complex_coefficient_transport():
    value = native.certificate()
    assert (value['LH_equation_rank'], value['LH_invariant_line_dimension']) == (3, 1)
    assert value['native_current_count'] == 12
    assert value['complete_K45_transport_exact']
    data = native.dictionary()
    assert data['permutation'].T*data['permutation'] == s.eye(15)
    assert data['line_scale']**2 == 1


def test_neutral_slot_and_wrong_LH_sign_are_checked_by_actual_currents():
    tensor = calibrated._LH_tensor()
    assert tensor[:, 1].dot(tensor[:, 1]) == 1
    neutral = tensor[:, 1]*tensor[:, 1].T
    assert neutral.rank() == 1
    q = action.electric_charge(families=False)
    assert q*neutral == s.zeros(15)
    wrong = tensor.copy().as_mutable()
    nonzero = list(wrong.todok())
    wrong[nonzero[0]] *= -1
    residuals = [g.T*wrong+wrong*h[3:, 3:] for g, h in zip(action.generators(families=False), action._representation()[1], strict=True)]
    assert any(not action._zero(row) for row in residuals)
