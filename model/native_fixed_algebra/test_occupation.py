"""Independent endpoint, graph, reference and locality controls for full occupation.

The oracle below uses bit occupations and ordered-token permutations.  It
does not call the certificate producer's signature or transport routines.
"""
from collections import defaultdict
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import json
from math import comb, prod
from pathlib import Path

import numpy as np
import pytest

from . import occupation as o
from ..exchange_sector import many_body as native


def _endpoint(mask, law, seam, coin, character):
    a, b = seam
    occupied = [p for p in range(12) if mask & (1 << p)]
    na, nb = int(a in occupied), int(b in occupied)
    new_a = (na + nb) // 2 + coin * ((na + nb) % 2)
    current = na - new_a
    moved = occupied.copy()
    if current:
        source, target = (a, b) if current == 1 else (b, a)
        moved[moved.index(source)] = target
    inversions = sum(x > y for i, x in enumerate(moved) for y in moved[i + 1:])
    sign = (-1) ** inversions if character == -1 else 1
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == 'E':
        raw += 4 * ((na + nb) % 2)
    return raw, sum(1 << p for p in moved), sign, current


@pytest.fixture(scope='module', params=tuple(product(('C', 'E'), (1, -1))))
def native_case(request):
    law, character = request.param
    return law, character, native.branches(law, character)


def test_actual_full_endpoint_instrument_matches_independent_ordered_tokens(native_case):
    law, character, branches = native_case
    assert len(branches) == (240 if law == 'C' else 480)
    expected = defaultdict(list)
    for seam, coin in product(native.SEAMS, (0, 1)):
        for mask in range(4096):
            raw, output, sign, _ = _endpoint(mask, law, seam, coin, character)
            expected[seam, coin, raw].append((mask, output, sign))
    for branch in branches:
        entries = expected[branch.seam, branch.coin, branch.raw]
        assert list(zip(branch.inputs, branch.outputs, branch.signs, strict=True)) == entries
        assert not branch.abort or not entries


def test_joint_endpoint_difference_signature_separates_every_nonextremal_state():
    # A connected graph of n_a-n_b differences determines the occupation
    # string up to a constant.  Only all-empty/all-full share a signature.
    groups = defaultdict(list)
    for mask in range(4096):
        key = tuple(((mask >> a) & 1) - ((mask >> b) & 1)
                    for a, b in native.SEAMS)
        groups[key].append(mask)
    assert len(groups) == 4095
    assert [row for row in groups.values() if len(row) != 1] == [[0, 4095]]
    for number in range(1, 12):
        start = (1 << number) - 1
        seen, todo = {start}, [start]
        while todo:
            mask = todo.pop()
            for a, b in native.SEAMS:
                if bool(mask & (1 << a)) != bool(mask & (1 << b)):
                    target = mask ^ (1 << a) ^ (1 << b)
                    if target not in seen:
                        seen.add(target)
                        todo.append(target)
        assert len(seen) == comb(12, number)
    assert sum(comb(12, n) - 1 for n in range(1, 12)) == 4083


def test_extremal_qubit_is_preserved_in_each_complete_raw_record(native_case):
    law, character, branches = native_case
    rng = np.random.default_rng(71103)
    logical_reference = rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3))
    logical_reference /= np.linalg.norm(logical_reference)
    # Unweighted native branches have a common fair amplitude 1/sqrt(2).
    # Both endpoints have raw zero, including E's endpoint-parity read.
    for seam in native.SEAMS:
        mass = 0.
        for branch in branches:
            if branch.seam != seam:
                continue
            mapping = dict(zip(branch.inputs, zip(branch.outputs, branch.signs), strict=True))
            expected = branch.raw == 0
            assert (0 in mapping) == (4095 in mapping) == expected
            if expected:
                assert mapping[0] == (0, 1) and mapping[4095] == (4095, 1)
                out = logical_reference / np.sqrt(2)
                np.testing.assert_allclose(out.conj().T @ out,
                                           logical_reference.conj().T @ logical_reference / 2)
                mass += np.linalg.norm(out) ** 2
        assert mass == pytest.approx(1, abs=4e-16)


def test_removing_all_seams_at_one_port_restores_local_conservation():
    retained = tuple(edge for edge in native.SEAMS if 0 not in edge)
    assert retained and len(retained) < len(native.SEAMS)
    for mask, seam, coin in product(range(4096), retained, (0, 1)):
        output = _endpoint(mask, 'E', seam, coin, -1)[1]
        assert (output & 1) == (mask & 1)
    assert any((_endpoint(1, 'C', seam, 0, 1)[1] & 1) == 0
               for seam in native.SEAMS if 0 in seam)


def test_global_extremal_memory_has_no_proper_region_privacy_or_local_repair_protection():
    # Every proper nonempty exterior can distinguish the two logical basis
    # states perfectly; the phase itself needs all ports to be read.
    for region in range(1, 4095):
        exterior = 4095 ^ region
        assert (0 & exterior) != (4095 & exterior)
        assert (0 & region) != (4095 & region)
    plus = np.array([[1, 1], [1, 1]], dtype=complex) / 2
    local_number_on_code = np.diag([0, 1])
    complement = np.eye(2) - local_number_on_code
    read_and_forget = (local_number_on_code @ plus @ local_number_on_code
                       + complement @ plus @ complement)
    np.testing.assert_array_equal(read_and_forget, np.eye(2) / 2)
    assert np.trace(plus @ read_and_forget).real == .5


def test_stationary_cat_references_allow_unfixed_modular_frequencies():
    # Both references are faithful on the SAME full carrier and stationary
    # under every native branch; their extremal modular frequencies differ.
    # Native records cannot select their relative sector normalization.
    for p0, p12 in ((F(1, 8), F(1, 4)), (F(1, 16), F(5, 16))):
        rest = (1 - p0 - p12) / 4094
        probabilities = [rest] * 4096
        probabilities[0], probabilities[-1] = p0, p12
        assert min(probabilities) > 0 and sum(probabilities) == 1
        assert p12 / p0 in (F(2), F(5))
        for law in ('C', 'E'):
            for seam in native.SEAMS:
                output = [F(0)] * 4096
                for coin, mask in product((0, 1), range(4096)):
                    target = _endpoint(mask, law, seam, coin, -1)[1]
                    output[target] += probabilities[mask] / 2
                assert output == probabilities
    assert np.log(2.) != np.log(5.)


def test_fixed_number_admission_removes_extremal_coherence_even_though_parity_allows_it():
    path = Path(__file__).resolve().parents[2] / 'model/completion/program_class_v2.json'
    program = json.loads(path.read_text(encoding='utf-8'))
    carriers = {row['id']: row for row in program['carriers']}
    assert 'cross-global-number coherence is not admitted' in carriers['binary_native_federation']['definition']
    assert 'fixed global total' in carriers['binary_native_federation']['definition']
    assert 'no arbitrary total or new species' in carriers['local_integer']['admissible_loads']
    cat_transfer = np.array([[0, 1], [0, 0]], dtype=int)
    number = np.diag([0, 12])
    parity = np.eye(2, dtype=int)
    np.testing.assert_array_equal(number @ cat_transfer - cat_transfer @ number,
                                  -12 * cat_transfer)
    np.testing.assert_array_equal(parity @ cat_transfer - cat_transfer @ parity,
                                  np.zeros((2, 2), dtype=int))


@pytest.mark.parametrize('weights', ((1,) * 12, (25, 25) + (8,) * 10, (1, 4) * 6))
def test_full_product_reference_matches_independent_partition_and_every_seam(weights):
    reference = o.weighted_reference(weights)
    denominator = prod(1 + value for value in weights)
    expected = tuple(F(prod(weights[p] for p in range(12) if mask & (1 << p)), denominator)
                     for mask in range(4096))
    assert reference == expected and sum(reference) == 1 and min(reference) > 0
    for seam in native.SEAMS:
        a, b = seam
        output = [F(0)] * 4096
        for coin, mask in product((0, 1), range(4096)):
            probability = F(weights[a] if coin else weights[b], weights[a] + weights[b])
            target = _endpoint(mask, 'C', seam, coin, 1)[1]
            output[target] += probability * expected[mask]
        assert tuple(output) == expected


@pytest.fixture(scope='module')
def algebra_cases():
    return {(law, sign): o.algebra_certificate(law, sign)
            for law, sign in product(('C', 'E'), (1, -1))}


@pytest.mark.parametrize('law,character', tuple(product(('C', 'E'), (1, -1))))
def test_certified_matrix_algebra_has_all_sectors_and_complete_signed_witnesses(algebra_cases, law, character):
    receipt = algebra_cases[law, character]
    assert o.verify_algebra(law, character, deepcopy(receipt)) is True
    assert receipt['sector_dimensions'] == [comb(12, n) for n in range(13)]
    assert receipt['algebra_dimension'] == sum(comb(12, n) ** 2 for n in range(1, 12)) + 1
    assert receipt['commutant_dimension'] == 15
    assert receipt['inventory']['raw_slots'] == (240 if law == 'C' else 480)
    assert receipt['inventory']['columns'] == 60 * 4096
    signatures = [0] * 4096
    for bit, (index, coin, raw) in enumerate(receipt['effect_labels']):
        for mask in range(4096):
            if _endpoint(mask, law, native.SEAMS[index], coin, character)[0] == raw:
                signatures[mask] |= 1 << bit
    assert signatures == receipt['effect_signatures']
    rows = native.branches(law, character)
    for source, target, index, phase in receipt['transport_forest']:
        branch = rows[index]
        raw, output, sign, _ = _endpoint(source, law, branch.seam, branch.coin, character)
        assert (raw, output, sign) == (branch.raw, target, phase)


@pytest.mark.parametrize('mutation', (
    'empty', 'missing', 'extra', 'float_character', 'float_dimension',
    'bool_signature', 'wrong_signature', 'duplicate_effect', 'bad_effect',
    'missing_root', 'wrong_sign', 'bad_branch', 'missing_edge',
    'duplicate_edge', 'zero_slot_removed', 'fake_sector',
))
def test_certificate_replay_rejects_malformed_and_forged_witnesses(algebra_cases, mutation):
    value = deepcopy(algebra_cases['E', -1])
    if mutation == 'empty':
        value = {}
    elif mutation == 'missing':
        del value['inventory']
    elif mutation == 'extra':
        value['unchecked_promotion'] = True
    elif mutation == 'float_character':
        value['character'] = -1.
    elif mutation == 'float_dimension':
        value['algebra_dimension'] = float(value['algebra_dimension'])
    elif mutation == 'bool_signature':
        value['effect_signatures'][0] = False
    elif mutation == 'wrong_signature':
        value['effect_signatures'][7] ^= 1
    elif mutation == 'duplicate_effect':
        value['effect_labels'][0] = value['effect_labels'][1].copy()
    elif mutation == 'bad_effect':
        value['effect_labels'][0][0] = 30
    elif mutation == 'missing_root':
        value['roots'].pop()
    elif mutation == 'wrong_sign':
        value['transport_forest'][0][-1] *= -1
    elif mutation == 'bad_branch':
        value['transport_forest'][0][-2] = 480
    elif mutation == 'missing_edge':
        value['transport_forest'].pop()
    elif mutation == 'duplicate_edge':
        value['transport_forest'][-1] = value['transport_forest'][0].copy()
    elif mutation == 'zero_slot_removed':
        value['inventory']['zero_slots'] -= 1
    elif mutation == 'fake_sector':
        value['sector_dimensions'][2] += 1
    with pytest.raises((ValueError, ArithmeticError)):
        o.verify_algebra('E', -1, value)


@pytest.mark.parametrize('law,character', (('', 1), ('c', 1), (None, 1), ('C', True), ('E', -1.), ('C', 0)))
def test_certificate_domain_rejects_implicit_coercion(law, character):
    with pytest.raises(ValueError):
        o.algebra_certificate(law, character)


def test_returned_certificate_mutation_does_not_change_next_proof(algebra_cases):
    first = o.algebra_certificate('C', 1)
    first['effect_signatures'][0] = 987
    first['transport_forest'][0][-1] = 987
    first['inventory']['columns'] = 0
    assert o.algebra_certificate('C', 1) == algebra_cases['C', 1]


@pytest.mark.parametrize('weights', ((), [1] * 12, (1,) * 11, (True,) * 12,
                                     (1.,) * 12, (0,) * 12, (-1,) * 12))
def test_reference_rejects_unadmitted_weights(weights):
    with pytest.raises(ValueError):
        o.weighted_reference(weights)


def test_reference_rejects_forged_fraction_internals():
    value = F(1, 2)
    object.__setattr__(value, '_denominator', 0)
    with pytest.raises(ValueError):
        o.weighted_reference((value,) + (F(1),) * 11)


def test_zero_raw_slots_are_part_of_the_mathematical_source_inventory(monkeypatch):
    actual = o._rows('E', -1)
    index = next(i for i, row in enumerate(actual) if not row.inputs)
    incomplete = actual[:index] + actual[index + 1:]
    monkeypatch.setattr(o, '_rows', lambda law, character: incomplete)
    with pytest.raises((ValueError, ArithmeticError)):
        o.algebra_certificate('E', -1)


def test_deleted_seam_cannot_be_silently_reclassified_as_the_complete_source(monkeypatch):
    actual = o._rows('C', 1)
    incomplete = tuple(row for row in actual if row.seam != native.SEAMS[0])
    monkeypatch.setattr(o, '_rows', lambda law, character: incomplete)
    with pytest.raises((ValueError, ArithmeticError)):
        o.algebra_certificate('C', 1)


def test_extremal_relative_kraus_phase_breaks_noiseless_multiplicity(monkeypatch):
    actual = list(o._rows('C', 1))
    index = next(i for i, row in enumerate(actual) if 4095 in row.inputs)
    row = actual[index]
    signs = list(row.signs)
    signs[row.inputs.index(4095)] *= -1
    actual[index] = replace(row, signs=tuple(signs))
    monkeypatch.setattr(o, '_rows', lambda law, character: tuple(actual))
    with pytest.raises((ValueError, ArithmeticError)):
        o.algebra_certificate('C', 1)


def test_public_reference_return_cannot_poison_later_reference_values():
    first = o.weighted_reference((1,) * 12)
    object.__setattr__(first[0], '_numerator', 2)
    later = o.weighted_reference((1,) * 12)
    assert later == (F(1, 4096),) * 4096


@pytest.mark.parametrize('law,character', tuple(product(('C', 'E'), (1, -1))))
def test_faithful_certificate_keeps_unselected_sector_weights_explicit(law, character):
    weights = (1, 2) * 6
    report = o.faithful_certificate(weights, law, character)
    coefficients = [0] * 13
    for mask in range(4096):
        coefficients[mask.bit_count()] += prod(weights[p] for p in range(12) if mask & (1 << p))
    assert list(map(F, report['sector_normalizers'])) == coefficients
    assert F(report['reference_normalizer']) == prod(1 + w for w in weights)
    assert report['columns_checked'] == 60 * 4096
    assert report['all_seams_complete_and_faithfully_invariant']
    assert not report['sector_populations_or_chemical_potential_selected']
    assert not report['weighted_source_service_executed']


@pytest.mark.parametrize('law,character', tuple(product(('C', 'E'), (1, -1))))
def test_peripheral_witness_is_actual_strictly_positive_native_diagonal(law, character):
    result = o.peripheral_certificate(law, character)
    diagonal = [0] * 4096
    rows = native.branches(law, character)
    for index in result['holding_branch_indices']:
        row = rows[index]
        for source, target, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            assert source == target and sign == 1
            diagonal[source] += 1
    assert diagonal == result['positive_diagonal']
    assert min(diagonal) == result['smallest_diagonal'] > 0
    assert max(diagonal) == result['largest_diagonal']
    assert result['only_peripheral_eigenvalue'] == 1
    assert not result['convergence_rate_claimed']
    assert 'no adaptive-history mixing assertion' in result['scope']


def test_locality_and_extremal_receipts_preserve_the_component_and_number_boundaries():
    local, extreme = o.locality_certificate(), o.extremal_certificate()
    assert local['proper_nonempty_core_sizes'] == list(range(1, 12))
    assert local['proper_port_local_fixed_algebra'] == 'C I'
    assert 'dimension13' in local['number_superselected_fixed_algebra']
    assert local['fixed_number_source_algebra'] == 'C I_N for each N=0,...,12'
    for row in local['overlap_witnesses']:
        for left, outside_left, right, outside_right in row['adjacent_count_witnesses']:
            assert left + outside_left == right + outside_right
            assert right == left + 1
            assert 0 <= outside_left <= row['outside_size']
    assert extreme['quantum_dimension'] == 2
    assert extreme['local_error_code_distance'] == 1
    assert extreme['cross_global_number_coherence_required']
    assert not extreme['class_V2_fixed_global_number_admission_satisfied']
    assert not extreme['full_code_locally_indistinguishable']
    assert not extreme['localized_topological_defect_constructed']
