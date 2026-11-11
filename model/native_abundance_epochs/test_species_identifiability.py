"""Independent linear algebra, extremal populations and physical-domain controls."""
from fractions import Fraction as F
import itertools
import json

import pytest

from . import species_identifiability as s


def rational_rank(matrix):
    rows = [[F(x) for x in row] for row in matrix]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        leading = rows[rank][col]
        rows[rank] = [x / leading for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank:
                factor = rows[i][col]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


@pytest.mark.parametrize("grammar,bosons,size,multiplicity", [
    ("full32", False, 32, (4, 6, 6, 4, 6, 6)),
    ("restricted30", False, 30, (2, 6, 6, 4, 6, 6)),
    ("full32", True, 48, (16, 6, 6, 8, 6, 6)),
    ("restricted30", True, 46, (14, 6, 6, 8, 6, 6)),
])
def test_partition_rank_and_complete_kernel_independently(grammar, bosons, size, multiplicity):
    rows = s.components(grammar, bosons=bosons)
    matrix, basis = s.matrix(rows), s.kernel_basis(rows)
    assert len(rows) == size
    assert tuple(map(sum, matrix)) == multiplicity
    assert rational_rank(matrix) == 6
    assert rational_rank(basis) == size - 6
    for vector in basis:
        assert all(sum(x * y for x, y in zip(row, vector)) == 0 for row in matrix)


def test_mask_classes_and_restricted_grammar_are_independent_bit_calculations():
    rows = s.components()
    assert [r["mask"] for r in rows] == list(range(32))
    for row in rows:
        mask = row["mask"]
        assert row["k"] == (3 * (mask >> 3).bit_count() - 2 * (mask & 7).bit_count()) % 6
    restricted = s.components("restricted30")
    assert {r["mask"] for r in restricted} == set(range(1, 31))
    assert all(r["species"] != "neutrino" or not r["right_neutrino"] for r in restricted)
    assert {r["name"] for r in rows if r["mask"] in (0, 31)} == {"nu_R", "anti-nu_R"}


def test_sharp_signed_bounds_equal_exhaustive_integer_polytope():
    # All integer populations with two units in class0 and three in class3.
    rows = tuple(dict(key=str(j), k=k) for j, k in enumerate((0, 0, 0, 3, 3)))
    census = (2, 0, 0, 3, 0, 0)
    weights = (-2, 4, 1, 7, -1)
    result = s.sharp_bounds(rows, census, weights)
    values = [sum(x * y for x, y in zip(population, weights))
              for population in itertools.product(range(4), repeat=5)
              if s.pushforward(rows, population) == census]
    assert (result["lower"], result["upper"]) == (min(values), max(values)) == (-7, 29)
    for side in ("lower", "upper"):
        assert s.pushforward(rows, result[side + "_witness"]) == census
    with pytest.raises(ValueError):
        s.sharp_bounds(rows, (2, 1, 0, 3, 0, 0), weights)


@pytest.mark.parametrize("grammar", s.GRAMMARS)
def test_nonuniform_census_identifies_quarks_but_not_up_vs_down(grammar):
    rows = s.components(grammar)
    census = (1, 2, 3, 4, 5, 6)
    results = s.infer_aggregates(rows, census, observable_kind="code_quantum_component_counts")
    assert (results["quark"]["lower"], results["quark"]["upper"]) == (16, 16)
    assert (results["lepton"]["lower"], results["lepton"]["upper"]) == (5, 5)
    assert (results["up"]["lower"], results["up"]["upper"]) == (0, 16)
    assert (results["neutrino"]["lower"], results["neutrino"]["upper"]) == (0, 5 if grammar == "full32" else 4)
    assert results["charged_lepton"]["lower"] == (0 if grammar == "full32" else 1)


def test_boson_extension_has_sharp_unidentifiable_fermion_content():
    rows = s.components("restricted30", bosons=True)
    result = s.infer_aggregates(rows, (2, 1, 1, 3, 1, 1), observable_kind="code_quantum_component_counts")
    for name, bounds in {"quark": (4, 4), "lepton": (0, 5), "fermion": (4, 9),
                         "boson": (0, 5), "gauge": (0, 2), "higgs": (0, 3)}.items():
        assert (result[name]["lower"], result[name]["upper"]) == bounds


def test_zero_populated_classes_do_not_create_false_local_ambiguity():
    rows = s.components()
    result = s.sharp_bounds(rows, (0, 1, 0, 0, 0, 0), s.coefficients(rows, "right_neutrino"))
    assert result["identifiable_for_this_census"]
    assert not result["identifiable_for_every_census"]
    assert s.sharp_bounds(rows, (0,) * 6, s.coefficients(rows, "neutrino"))["upper"] == 0


def test_representation_counts_do_not_select_a_population():
    rows = s.components()
    quark = next(j for j, r in enumerate(rows) if r["family"] == "quark")
    lepton = next(j for j, r in enumerate(rows) if r["family"] == "lepton")
    for selected, expected in ((quark, 1), (lepton, 0)):
        population = [int(j == selected) for j in range(len(rows))]
        census = s.pushforward(rows, population)
        q = s.sharp_bounds(rows, census, s.coefficients(rows, "quark"))
        assert q["lower"] == q["upper"] == expected


@pytest.mark.parametrize("grammar", s.GRAMMARS)
def test_minimal_additional_label_and_pigeonhole_obstruction(grammar):
    rows = s.components(grammar, bosons=True)
    families = s.minimal_extra_label(rows, ("quark", "lepton", "higgs", "gauge"))
    assert families["minimum_alphabet_size"] == 2
    assert families["minimum_fixed_binary_bits"] == 1
    assert families["joint_class_label_rank"] == 8
    all_species = s.minimal_extra_label(rows, s.AGGREGATES)
    assert all_species["minimum_alphabet_size"] == 3
    assert all_species["minimum_fixed_binary_bits"] == 2
    assert all_species["joint_class_label_rank"] == (14 if grammar == "full32" else 13)
    # Three different aggregate signatures share class3. Two labels cannot separate them.
    indices = [next(i for i, row in enumerate(rows) if row["name"] == name)
               for name in ("e_L", "nu_L", "H")]
    assert all(rows[i]["k"] == 3 for i in indices)
    signatures = [tuple(s.coefficients(rows, a)[i] for a in s.AGGREGATES) for i in indices]
    assert len(set(signatures)) == 3
    for two_labels in itertools.product(range(2), repeat=3):
        assert any(two_labels[i] == two_labels[j] and signatures[i] != signatures[j]
                   for i in range(3) for j in range(i))
    assert not all_species["physical_instrument_implemented"]


@pytest.mark.parametrize("invalid", ["native_flux_faces", "thermal_number", "thermal_energy", ""])
def test_face_and_thermal_records_cannot_silently_enter_species_decoder(invalid):
    with pytest.raises(ValueError, match="decoder"):
        s.infer_aggregates(s.components(), (1, 0, 0, 0, 0, 0), observable_kind=invalid)


@pytest.mark.parametrize("invalid", [[0] * 31, [-1] + [0] * 31, [0.0] * 32,
                                     [True] + [0] * 31, [float("nan")] * 32])
def test_invalid_population_is_rejected(invalid):
    with pytest.raises(ValueError):
        s.pushforward(s.components(), invalid)


def test_thermal_helicity_count_number_and_energy_normalization():
    result = s.thermal_benchmark()
    assert result["degree_counts"] == {"fermion": 90, "boson": 28}
    assert len(result["physical_components"]) == 118
    assert not any(row["right_neutrino"] for row in result["physical_components"])
    number, energy = (result["quantities"][key] for key in ("number", "energy"))
    assert number["effective_degrees"] == F(191, 2)
    assert energy["effective_degrees"] == F(427, 4)
    assert number["class_fractions"] == tuple(F(x, 191) for x in (57, 27, 27, 26, 27, 27))
    assert energy["class_fractions"] == tuple(F(x, 427) for x in (117, 63, 63, 58, 63, 63))
    for quantity, denominator, expected in ((number, 191, (108, 27, 8, 48)),
                                           (energy, 427, (252, 63, 16, 96))):
        assert sum(quantity["class_fractions"]) == 1
        assert [quantity["aggregate_fractions"][key] for key in ("quark", "lepton", "higgs", "gauge")] == [F(x, denominator) for x in expected]
        assert quantity["quark_fraction_among_fermions"] == F(4, 5)
    assert number["class_fractions"] != energy["class_fractions"]
    # Extra spin doubling of Weyl components or 7/8 for number density fails.
    assert 28 + F(7, 8) * 180 != energy["effective_degrees"]
    assert 28 + F(7, 8) * 90 != number["effective_degrees"]


def test_receipt_replays_and_contains_attained_indistinguishable_populations():
    result = s.receipt()
    assert result == json.loads(s.RESULT.read_text())
    assert result["all_checks_pass"] and not result["physical_species_claim"]
    for example in result["indistinguishable_populations"]:
        rows = s.components(example["grammar"], bosons=example["bosons"])
        populations = [[F(x) for x in p] for p in example["populations"]]
        assert s.pushforward(rows, populations[0]) == s.pushforward(rows, populations[1])
        assert example["aggregate_values"][0] != example["aggregate_values"][1]
