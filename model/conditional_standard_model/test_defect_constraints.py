"""Interaction controls for proposed protected defect charges."""
import pytest

from .defect_constraints import SPECIES, allowed_assignment, constraint_certificate


def assignment(d=1, e=3):
    return dict(zip(SPECIES, [d-e]*3+[d]*3+[e]*4+[0]*2, strict=True))


def test_integer_and_torsion_groups_keep_the_instanton_factor_three():
    value = constraint_certificate()
    p = value["perturbative_vertices"]
    assert p["universal_additive_charge_group_free_rank"] == 2
    assert p["universal_additive_charge_group_torsion_orders"] == []
    i = value["optional_electroweak_instanton"]
    assert i["universal_additive_charge_group_free_rank"] == 1
    assert i["universal_additive_charge_group_torsion_orders"] == [3]
    assert allowed_assignment(assignment(), electroweak_instanton=True)


@pytest.mark.parametrize("name", ["muon", "tau", "charm", "strange", "Z", "Higgs"])
def test_independent_family_winding_breaks_a_required_interaction(name):
    charges = assignment()
    charges[name] += 1
    assert not allowed_assignment(charges)


def test_baryon_direction_and_finite_group_control():
    baryon = assignment(d=1, e=0)
    assert allowed_assignment(baryon)
    assert not allowed_assignment(baryon, electroweak_instanton=True)
    # A torsion target cannot be treated by division by three.
    torsion = assignment(d=0, e=1)
    assert not allowed_assignment(torsion, electroweak_instanton=True)
    assert allowed_assignment(torsion, modulus=3, electroweak_instanton=True)


def test_malformed_or_incomplete_dictionaries_are_rejected():
    with pytest.raises(ValueError, match="complete"):
        allowed_assignment({"electron": 1})
    with pytest.raises(ValueError, match="integer"):
        allowed_assignment({**assignment(), "electron": 3.0})
    with pytest.raises(ValueError, match="modulus"):
        allowed_assignment(assignment(), modulus=True)
