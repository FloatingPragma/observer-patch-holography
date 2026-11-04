"""Integer presentations, finite cyclic controls and mixing counterexamples."""
from itertools import product

import pytest

from . import topology as module


def assignment(d, e, n):
    return dict(zip(module.SPECIES, [d-e+n]*3+[d]*3+[e]*3+[e-n, 0, 0]+[n]*3, strict=True))


def test_exact_groups_retain_majorana_torsion_and_massless_exclusion():
    value = module.certificate()
    expected = {'all_three_positive_Majorana_masses': (2, [2]),
                'all_three_positive_with_optional_instanton': (1, [6]),
                'calibrated_one_massless_excluded': (2, []),
                'calibrated_one_massless_with_optional_instanton': (1, [3])}
    for key, pair in expected.items():
        assert (value[key]['free_rank'], value[key]['torsion_invariant_factors']) == pair
    assert value['instanton_relation_before_reduction'] == '9*D-3*E+6*N=0'


@pytest.mark.parametrize('modulus', [2, 3, 4, 6])
@pytest.mark.parametrize('instanton', [False, True])
def test_independent_cyclic_assignment_counts_match_presented_group(modulus, instanton):
    from math import gcd
    group = module.universal_group((True, True, True), electroweak_instanton=instanton)
    count = sum(module.allowed_assignment(assignment(d, e, n), (True, True, True), modulus=modulus,
                                         electroweak_instanton=instanton)
                for d, e, n in product(range(modulus), repeat=3))
    expected = modulus**group['free_rank']
    for divisor in group['torsion_invariant_factors']: expected *= gcd(modulus, divisor)
    assert count == expected


def test_full_instanton_retains_fermion_parity_and_excluded_zero_neutrino_removes_it():
    parity = assignment(1, 1, 1)
    assert module.allowed_assignment(parity, (True, True, True), modulus=2, electroweak_instanton=True)
    assert not module.allowed_assignment(parity, (False, True, True), modulus=2, electroweak_instanton=True)
    assert module.allowed_assignment(assignment(1, 3, 0), electroweak_instanton=True)
    assert not module.allowed_assignment(assignment(0, 0, 1), (True, True, True))


def test_disconnected_PMNS_does_not_imply_generation_charge_identification():
    diagonal = ((0, 0), (1, 1), (2, 2))
    group = module.universal_group((True, True, True), pmns_edges=diagonal)
    assert (group['free_rank'], group['torsion_invariant_factors']) == (2, [2, 2, 2])
    values = assignment(0, 0, 0)
    values['electron'] = values['nu_1'] = 1
    assert module.allowed_assignment(values, (True, True, True), modulus=2, pmns_edges=diagonal)
    assert not module.allowed_assignment(values, (True, True, True), modulus=2)
    # Even a connected non-complete graph gives the same universal relation.
    tree = ((0, 0), (0, 1), (1, 1), (1, 2), (2, 2))
    connected = module.universal_group((True, True, True), pmns_edges=tree)
    assert (connected['free_rank'], connected['torsion_invariant_factors']) == (2, [2])


def test_no_positive_majorana_mass_and_no_massless_exclusion_leave_neutrino_integer_charge_free():
    group = module.universal_group((False, False, False), exclude_massless=False)
    assert (group['free_rank'], group['torsion_invariant_factors']) == (3, [])
    assert module.allowed_assignment(assignment(0, 0, 1), (False, False, False), exclude_massless=False)


@pytest.mark.parametrize('options', [{'massive': (0, True, True)}, {'exclude_massless': 1},
                                     {'pmns_edges': ([0, 1],)}, {'pmns_edges': ((0, 0), (0, 0))},
                                     {'pmns_edges': ((True, 0),)}])
def test_topological_hypotheses_require_exact_types(options):
    with pytest.raises(ValueError): module.relation_matrix(**options)
