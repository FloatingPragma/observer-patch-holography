"""Independent label, restricted-domain and incompatible-character controls."""
from fractions import Fraction as F

import pytest

from . import many_body_compatibility as c
from . import many_body_topology as t
from . import many_body


def test_both_complete_families_and_extra_sector_phase_rejection():
    for sign in (1, -1):
        assert c.validate_family(c.character_family(sign)) == sign
    family = list(c.character_family(1))
    family[3] = (F(1, 4),)  # An allowed Z4 character, but not spectator stable.
    with pytest.raises(ValueError, match="spectators"):
        c.validate_family(tuple(family))
    family[3] = (F(1, 3),)
    with pytest.raises(ValueError, match="torsion"):
        c.validate_family(tuple(family))
    family = list(c.character_family(-1))
    family[11] = (F(0),)*19
    with pytest.raises(ValueError, match="one-hole"):
        c.validate_family(tuple(family))


@pytest.mark.parametrize("bad", (True, -1.0, F(-1), 0, 2, "-1"))
def test_character_type_aliases_rejected(bad):
    c.character_family(-1)
    with pytest.raises(ValueError):
        c.character_family(bad)
    with pytest.raises(ValueError):
        c.insertion_phase((1,), 0, bad)


@pytest.mark.parametrize("bad", ((False,), (1.0,), (1, 1), (2, 1), [1], (12,)))
def test_state_type_and_domain_aliases_rejected(bad):
    c.movement_sign((0,), (5,))
    with pytest.raises(ValueError):
        c.movement_sign(bad, bad)
    with pytest.raises(ValueError):
        c.insertion_phase(bad, 3)


def test_exact_family_schema_rejects_missing_and_numeric_aliases():
    family = c.character_family(1)
    for bad in (family[:-1], list(family)):
        with pytest.raises(ValueError):
            c.validate_family(bad)
    for bad in (0, False, 0.0, F(-1), F(1)):
        altered = list(family)
        altered[2] = (bad,)
        with pytest.raises(ValueError):
            c.validate_family(tuple(altered))


@pytest.mark.parametrize("numerator,denominator", (
    (False, 1), (0.0, 1), (1, True), (1, 20.0), (1, 0), (-1, -20),
    (0, 20), (2, 40), (20, 20), (-1, 20), (1, -20)))
def test_mutated_fraction_storage_is_rejected_before_arithmetic(numerator, denominator):
    bad = F(0)
    object.__setattr__(bad, "_numerator", numerator)
    object.__setattr__(bad, "_denominator", denominator)
    family = list(c.character_family(1))
    family[2] = (bad,)
    with pytest.raises(ValueError, match="canonical Fraction"):
        c.validate_family(tuple(family))
    with pytest.raises(ValueError, match="canonical Fraction"):
        c.uniform_hole_character(bad)
    # Canonical exact arithmetic continues to admit both cross-sector signs.
    assert c.validate_family(c.character_family(1)) == 1
    assert c.validate_family(c.character_family(-1)) == -1
    assert c.uniform_hole_character(F(1, 20))


def test_all_movements_agree_with_independent_between_port_count():
    seams, _ = t._source()
    for mask in range(4096):
        state = tuple(p for p in range(12) if mask >> p & 1)
        for a, b in seams:
            if bool(mask >> a & 1) == bool(mask >> b & 1):
                continue
            target = tuple(p for p in range(12) if (mask ^ (1 << a) ^ (1 << b)) >> p & 1)
            expected = (-1)**sum(mask >> p & 1 for p in range(a+1, b))
            assert c._ordered_parity(state, target) == (expected == -1)
    assert c.movement_sign((0,), (5,), 1) == 1
    with pytest.raises(ValueError):
        c.movement_sign((0,), (3,))  # Antipodal, not a source seam.


def test_every_permitted_spectator_cell_and_edge_embeds():
    result = c.spectator_certificate()
    assert result["vertex_inclusions"] == 24576
    assert result["edge_inclusions"] == 153600
    assert result["particle_triangle_inclusions"] == 46080
    assert result["disjoint_square_inclusions"] == 322560
    assert result["whole_complex_creation_map_claimed"] is False


def test_complete_native_columns_intertwine_on_the_declared_domain():
    # All configurations, both laws/coins/signs; a seam avoiding the spectator.
    seam, p = t._source()[0][0], 11
    if p in seam:
        p = next(i for i in range(12) if i not in seam)
    for mask in range(4096):
        if mask >> p & 1:
            continue
        state = tuple(i for i in range(12) if mask >> i & 1)
        for law in ("C", "E"):
            for coin in (0, 1):
                for character in (1, -1):
                    raw, out, phase, current = many_body.column(mask, law, seam, coin, character)
                    raw2, out2, phase2, current2 = many_body.column(mask | (1 << p), law, seam, coin, character)
                    target = tuple(i for i in range(12) if out >> i & 1)
                    assert (raw2, out2, current2) == (raw, out | (1 << p), current)
                    assert phase*c.insertion_phase(target, p, character) == c.insertion_phase(state, p, character)*phase2


def test_addition_maps_only_the_avoiding_subcomplex_not_all_null_homotopies():
    g = t.geometry(1)
    a, b, p = g.faces[0]
    index = {state[0]: i for i, state in enumerate(g.states)}
    backtrack = (index[a], index[b], index[a])
    image = c.insert_path(1, backtrack, p)
    assert image[0] == image[-1]
    face_boundary = (index[a], index[b], index[p], index[a])
    with pytest.raises(ValueError, match="entire path"):
        c.insert_path(1, face_boundary, p)
    with pytest.raises(ValueError):
        c.insert_path(12, (0, 0), 0)
    for invalid in (True, 1.0, -1, 12):
        with pytest.raises(ValueError):
            c.insert_path(1, backtrack, invalid)


def test_one_hole_relation_is_integral_and_complement_does_not_preserve_cells():
    g = t.geometry(11)
    vectors = tuple(t.loop_class(11, path) for _, _, path in t.exchange_loops(11))
    assert len(vectors) == 20 and all(len(v) == 19 for v in vectors)
    assert all(sum(v[j] for v in vectors) == 0 for j in range(19))
    assert all(abs(t._determinant(vectors[:i]+vectors[i+1:])) == 1 for i in range(20))
    obstruction = c.particle_hole_obstruction()
    assert obstruction["one_particle_face_is_filled"]
    assert not obstruction["one_hole_image_is_filled"]
    assert any(obstruction["one_hole_homology_coordinates"])


def test_twenty_uniform_one_hole_phases_exist_but_only_two_extend():
    vectors = tuple(t.loop_class(11, path) for _, _, path in t.exchange_loops(11))
    observed = set()
    for k in range(20):
        turn = F(k, 20)
        coordinates = c.uniform_hole_character(turn)
        assert all(sum((a*b for a, b in zip(row, coordinates)), F(0)) % 1 == turn
                   for row in vectors)
        observed.add(coordinates)
        family = list(c.character_family(-1 if k == 10 else 1))
        family[11] = coordinates
        if k in (0, 10):
            assert c.validate_family(tuple(family)) == (1 if k == 0 else -1)
        else:
            with pytest.raises(ValueError, match="one-hole"):
                c.validate_family(tuple(family))
    assert len(observed) == 20
    for invalid in (0, False, 0.0, F(1, 21), F(1), F(-1, 20)):
        with pytest.raises(ValueError):
            c.uniform_hole_character(invalid)


def test_reversal_and_vertex_gauge_cannot_change_exchange_holonomy():
    # Deterministic nonconstant U(1) endpoint gauge, in exact rational turns.
    for n in (2, 6, 11):
        g = t.geometry(n)
        _, _, path = t.exchange_loops(n)[0]
        gauge = tuple(F((i*i+3*i+2) % 17, 17) for i in range(len(g.states)))
        turns = sum((F(c._ordered_parity(g.states[a], g.states[b]), 2)
                     + gauge[b]-gauge[a] for a, b in zip(path, path[1:])), F(0)) % 1
        assert turns == F(1, 2)


def test_certificate_replays_exact_groups_and_rejects_forged_ordered_lift(monkeypatch):
    result = c.certificate()
    assert result["compatible_scalar_families_up_to_vertex_gauge"] == 2
    assert result["ordered_lift_flat_cells"] == 90880
    assert result["ordered_lift_odd_exchange_loops"] == 10240
    assert result["physical_statistics_or_particle_identity_selected"] is False
    monkeypatch.setattr(c, "_ordered_parity", lambda source, target: 0)
    with pytest.raises(ArithmeticError, match="swap ordered labels"):
        c.certificate()
