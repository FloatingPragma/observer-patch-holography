from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

import pytest

from . import restriction as r, packet


def _distance(first, second):
    return sum(abs(first.get(key, F())-second.get(key, F()))
               for key in first.keys() | second.keys())/2


def test_exact_formula_nesting_distances_and_bounds_against_native_partial_trace():
    result = r.exact_restriction_certificate()
    assert result['actual_partial_trace_cases'] == 284
    assert result['nested_restriction_cases'] > 1000
    assert result['all_number_sectors_including_zeros_retained']


def test_every_retained_subset_has_the_same_formula():
    population = 6
    for occupied in range(population+1):
        density = r._canonical_density(population, occupied)
        for retained in range(population+1):
            expected = r._window_density(population, occupied, retained)
            for keep in combinations(range(population), retained):
                assert r.native.occupation_partial_trace(density, population, keep) == expected


def test_owning_two_port_formula_is_reused_exactly(monkeypatch):
    original = r.owning.initial_weight
    calls = []
    def wrapped(*args):
        calls.append(args)
        return original(*args)
    monkeypatch.setattr(r.owning, 'initial_weight', wrapped)
    for population in (2, 240, 960, 1000000):
        for occupied in (0, population//3, population):
            for ones in range(3):
                assert r.pattern_weight(population, occupied, 2, ones) == original(
                    population, occupied, (int(ones >= 1), int(ones == 2)))
    assert len(calls) == 36


def test_number_law_retains_impossible_sectors_and_large_local_number():
    assert r.number_law(240, 0, 12) == (F(1),)+(F(),)*12
    assert r.number_law(240, 240, 12) == (F(),)*12+(F(1),)
    law = r.number_law(240, 239, 12)
    assert len(law) == 13 and law[:11] == (F(),)*11
    assert law[11:] == (F(1, 20), F(19, 20))
    assert r.number_law(5, 2, 5)[2] == 1


def test_canonical_distance_is_only_a_lower_bound_for_other_sector_states():
    reduced = r._window_density(4, 2, 3)
    pure_one = {(1, 0, 1, 0): F(1)}
    canonical_one = r._canonical_density(3, 1)
    lower = r.canonical_sector_distance(4, 2, 3, 1)
    assert lower == F(1, 2)
    assert _distance(reduced, canonical_one) == lower
    assert _distance(reduced, pure_one) == F(5, 6) > lower


def test_nontrivial_proper_restriction_has_positive_number_variance():
    for population in range(2, 16):
        for occupied in range(1, population):
            for retained in range(1, population):
                assert r.number_variance(population, occupied, retained) > 0
    assert r.number_variance(1, 1, 1) == 0
    assert r.number_variance(12, 0, 6) == 0
    assert r.number_variance(12, 12, 6) == 0
    assert r.number_variance(12, 3, 12) == 0


def test_independent_coherent_reference_extension_is_preserved():
    population, occupied, keep = 5, 2, (0, 2, 4)
    density = r._canonical_density(population, occupied)
    reference = {(a,b): F(1, 2) for a,b in product(range(2), repeat=2)}
    extended = {(i,a,j,b): value*coefficient for (i,_,j,_),value in density.items()
                for (a,b),coefficient in reference.items()}
    actual = r.native.occupation_partial_trace(extended, population, keep)
    expected = {(i,a,j,b): value*coefficient
                for (i,_,j,_),value in r._window_density(population,occupied,len(keep)).items()
                for (a,b),coefficient in reference.items()}
    assert actual == expected
    assert r.native.reference_marginal([actual]) == reference


def test_equal_source_marginal_does_not_identify_correlated_reference():
    density = {(1,0,1,0): F(1,2), (2,1,2,1): F(1,2)}
    actual = r.native.occupation_partial_trace(density, 2, (0,))
    independent = {(i,a,i,a): F(1,4) for i,a in product(range(2), repeat=2)}
    assert _distance(actual, independent) == F(1,2)
    marginal = defaultdict(F)
    for (i,a,j,b),value in actual.items():
        if a == b:
            marginal[i,0,j,0] += value
    assert dict(marginal) == r._window_density(2,1,1)


@pytest.mark.parametrize('law', ('C','E'))
def test_actual_complete_native_records_preserve_each_internal_canonical_sector(law):
    edge = r.native.federation_source.federation(0).seams[0]
    extra = tuple(i for i in range(12) if i not in edge)[:4]
    keep = tuple(sorted((*edge,*extra)))
    for occupied in range(len(keep)+1):
        value = r.native.program(federation={'level':0,'total':occupied}, commands=[{
            'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':law,'repeat':1}])
        density = {(mask,a,mask,b): F(1,2*comb(len(keep),occupied))
                   for chosen in combinations(keep,occupied)
                   for mask in (sum(1 << i for i in chosen),)
                   for a,b in product(range(2),repeat=2)}
        output = r.native.execute(value,density)
        assert output.status == 'COMPLETE' and len(output.branches) == 25
        assert any(not branch.density for branch in output.branches)
        total = defaultdict(F)
        for branch in output.branches:
            for key, coefficient in branch.density.items():
                total[key] += coefficient
        assert dict(total) == density
        assert sum(r.native.trace(branch.density) for branch in output.branches) == 1


def test_big_integer_native_population_is_evaluated_without_global_carrier():
    population = r.federation_population(200)
    assert population.bit_length() == 408
    law = r.number_law(population,population//3,12)
    assert len(law) == 13 and sum(law) == 1
    comparison = r.bernoulli_law(12,F(1,3))
    distance = sum(abs(a-b) for a,b in zip(law,comparison,strict=True))/2
    assert 0 < distance <= r.bernoulli_bound(population,12) == F(66,population)
    assert r.nested_law(population,population//3,12,4) == r.number_law(population,population//3,4)


def test_fixed_total_dilution_bound_does_not_require_full_state_storage():
    distances = []
    for level in range(5):
        population = r.federation_population(level)
        distance, bound = r.vacuum_distance(population,2,12)
        assert distance <= bound == F(24,population)
        distances.append(distance)
    assert all(a > b for a,b in zip(distances,distances[1:]))


@pytest.mark.parametrize('args', [
    (0,0,0,0), (True,0,0,0), (2,True,1,0), (2,3,1,0), (2,-1,1,0),
    (2,1,3,0), (2,1,1,2), (2,1,1,-1), (2,1,F(1),0), (2.0,1,1,0)])
def test_malformed_pattern_domains_rejected(args):
    with pytest.raises(ValueError):
        r.pattern_weight(*args)


def test_empty_window_and_no_silent_conditioning():
    assert r.number_law(240,2,0) == (F(1),)
    assert r.nested_law(240,2,0,0) == (F(1),)
    assert r.vacuum_distance(240,2,0) == (F(),F())
    assert r.bernoulli_bound(240,0) == 0
    with pytest.raises(ValueError):
        r.nested_law(240,2,4,5)
    with pytest.raises(ValueError):
        r.canonical_sector_distance(240,2,3,4)
    with pytest.raises(ValueError):
        r.bernoulli_law(2,.5)


def test_receipt_rejects_semantic_mutation(tmp_path, monkeypatch):
    expected = {'source_bound':True,'all_number_sectors':[0,1,2]}
    monkeypatch.setattr(packet,'produce',lambda:expected)
    destination = tmp_path/'result.json'
    destination.write_text(packet.canonical(expected))
    assert packet.verify(destination)
    destination.write_text(packet.canonical({'source_bound':True,'all_number_sectors':[1]}))
    with pytest.raises(ValueError):
        packet.verify(destination)


def test_source_bindings_include_proof_tests_and_owning_canonical_interface():
    bindings = packet.source_bindings()
    for name in ('__init__.py','restriction.py','packet.py','test_restriction.py',
                 'product.py','test_product.py','PRODUCT.md','README.md','PROOF.md'):
        assert 'model/native_capacity_restriction/'+name in bindings
    assert 'model/observer_time_fields/canonical_sector.py' in bindings
    assert 'model/native_program/interpreter.py' in bindings
    assert 'model/completion/program_class_v2.json' in bindings
    assert 'Lean/ObserverPatchHolography/Locality/NoSignalling.lean' in bindings
    assert not any(name.startswith('docs/') for name in bindings)
