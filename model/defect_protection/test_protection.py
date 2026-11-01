"""Independent endpoint/order tests and adversarial memory-law checks."""
from copy import deepcopy
import math

import numpy as np
import pytest

from ..exchange_sector import statistics as exchange, statistics_service as source
from ..observer.source_scalar import compiler
from . import protection as p, service


def oracle(row):
    """Track named first/second tokens, using endpoint balancing directly."""
    result = {}
    a, b = row.seam
    for i, pair in enumerate(exchange.PAIRS):
        occupied = set(pair)
        total = int(a in occupied)+int(b in occupied)
        final_left = total//2+row.coin*(total % 2)
        current = int(a in occupied)-final_left
        raw = {0: 0, 1: 1, -1: 2}[current]
        if row.law == 'E':
            raw += 4*(total % 2)
        if raw != row.raw:
            continue
        for r in (0, 1):
            tokens = list(pair if r == 0 else pair[::-1])
            if current:
                origin, target = (a, b) if current == 1 else (b, a)
                tokens[tokens.index(origin)] = target
            ordered = tuple(sorted(tokens))
            final_r = int(tuple(tokens) != ordered)
            result[r*66+i] = final_r*66+exchange.INDEX[ordered]
    return result


@pytest.mark.parametrize('law', ['C', 'E'])
def test_every_lifted_column_against_named_token_oracle(law):
    for row in p.branches(law):
        assert dict(p.lifted_pairs(row)) == oracle(row)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_exact_charge_decomposition_on_all_columns(law):
    for row in p.branches(law):
        lifted = dict(p.lifted_pairs(row))
        for character in (-1, 1):
            for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
                result = {lifted[i]: 1, lifted[66+i]: character}
                phase = sign if character == -1 else 1
                assert result == {j: phase, 66+j: character*phase}


@pytest.mark.parametrize('law', ['C', 'E'])
def test_each_chart_generator_conjugates_every_complete_map(law):
    for vertex in range(66):
        chart = tuple(int(i == vertex) for i in range(66))
        def chart_map(index):
            r, i = divmod(index, 66)
            return (r ^ chart[i])*66+i
        for row in p.branches(law):
            wanted = {(chart_map(a), chart_map(b)) for a, b in p.lifted_pairs(row)}
            assert set(p.lifted_pairs(row, chart)) == wanted


@pytest.mark.parametrize('law', ['C', 'E'])
def test_algebra_constructive_certificate_and_nontrivial_closed_loop(law):
    value = p.algebra_certificate(law)
    assert value['commutant_complex_dimension'] == 2
    assert value['generated_star_algebra_complex_dimension'] == 8712
    assert len(set(value['source_effect_signatures'])) == 66
    assert len(value['source_spanning_tree']) == 65
    assert value['ordered_cover_reachable_states'] == 132
    loop = p.exchange_loop(law)
    maps = {(r.law, r.seam, r.coin, r.raw): oracle(r) for r in p.branches(law)}
    at = exchange.INDEX[tuple(loop['initial_pair'])]
    image = {i: i for i in range(132)}
    for label in loop['history']:
        key = (label[0], tuple(label[1]), label[2], label[3])
        image = {i: maps[key][j] for i, j in image.items() if j in maps[key]}
    assert image == {at: at+66, at+66: at}
    assert loop['probability_on_initial_pair'] == '1/8'


def test_memory_loss_and_reference_localization_exact_matrices():
    X = np.array([[0, 1], [1, 0]])
    Z = np.diag([1, -1])
    states = [(np.eye(2)+sign*X)/2 for sign in (-1, 1)]
    assert np.array_equal(X@states[0]@X, states[0])
    assert np.array_equal(Z@states[0]@Z, states[1])
    for rho in states:
        assert np.trace(X@np.diag(np.diag(rho))) == 0
        assert np.array_equal((rho+Z@rho@Z)/2, np.eye(2)/2)
    for probability in (0, .125, .5, 1):
        channel = lambda rho: (1-probability)*rho+probability*Z@rho@Z
        rho = states[1]
        for count in range(5):
            assert np.trace(X@rho) == (1-2*probability)**count
            rho = channel(rho)
    value = p.memory_and_localization_certificate()
    assert value['source_trace_distance'] == '0'
    assert not value['source_localized_defect']
    assert not value['sector_coherence_protected']


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_full_register_feedback_padding_reference_and_charge(current):
    rng = np.random.default_rng(1337+current)
    state = rng.normal(size=(1024, 2))+1j*rng.normal(size=(1024, 2))
    state /= np.linalg.norm(state)
    target = state.copy()
    seam = (0, 5)
    for pair in exchange.PAIRS:
        destination = 5 if current == 1 else 0
        origin = 0 if current == 1 else 5
        # Independent named-token ordering crossing criterion.
        if (current and destination in pair and origin not in pair
                and any(0 < spectator < 5 for spectator in pair if spectator != destination)):
            at = source._address(pair)
            target[[at, at+512]] = state[[at+512, at]]
    actual = service.apply_memory(seam, current, state)
    np.testing.assert_allclose(actual, target, atol=2e-10, rtol=0)
    flipped = np.concatenate((state[512:], state[:512]))
    moved = service.apply_memory(seam, current, flipped)
    np.testing.assert_allclose(moved, np.concatenate((actual[512:], actual[:512])), atol=2e-10, rtol=0)


def test_exact_phase_audit_rejects_missing_scalar_term():
    phases = service.controlled_parity_quarters((0, 5), 1)
    word = source.diagonal_word(phases)
    assert word[0]['kind'] == 'I'
    with pytest.raises(ArithmeticError):
        source._audit_diagonal_word(phases, word[1:])


def test_all_seam_feedback_join_rejects_nonwitness_parity_mutation(monkeypatch):
    original = source.phase_quarters
    changed_seam = (0, 11)
    address = source._address((5, 11))
    def wrong(seam, current, character):
        values = list(original(seam, current, character))
        if seam == changed_seam and current == 1:
            values[address] = (values[address]+2) % 4
        return tuple(values)
    monkeypatch.setattr(source, 'phase_quarters', wrong)
    with pytest.raises(ArithmeticError, match='token order'):
        service.audit_feedback()


@pytest.fixture(scope='module')
def executed():
    return service.certificate()


def test_complete_source_execution_and_resources(executed):
    assert executed['all_seam_current_feedback_words_audited'] == 90
    assert executed['native_source_columns_joined_to_exact_feedback'] == 7920
    for law, count, aborts in [('C', 8, 2), ('E', 16, 6)]:
        row = executed['laws'][law]
        assert row['raw_slots'] == count
        assert row['abort_slots'] == aborts
        assert row['all_joint_input_columns'] == 132
        assert row['terminal_addresses_per_slot'] == 1024
        resources = row['native_and_feedback_resource_rows']
        assert len(resources) == count
        assert all(r['proper_code_service_events'] == 2*r['native_pulses'] for r in resources)
        assert all(r['proper_code_complement_outcome_slots'] ==
                   4*r['single_code_pulses']+2*r['pair_code_pulses'] for r in resources)
        assert any(r['orientation_feedback_pulses'] for r in resources)


@pytest.mark.parametrize('character', [-1, 1])
def test_actual_memory_preparation_native_continuation_and_complete_readout(character):
    state = np.zeros(1024, complex)
    state[source._address((0, 1))] = 1
    prepare = [compiler.pulse('Y', -character*math.pi/4, [9])]
    state = compiler.execute(prepare, 10, state)
    for law in ('C', 'E'):
        probabilities = np.zeros(2)
        for row in service.repair_branches(law, (0, 5), state):
            output = compiler.execute([compiler.pulse('Y', math.pi/4, [9])], 10, row['amplitude'])
            probabilities += np.sum(abs(output.reshape(2, 512))**2, axis=1)
        np.testing.assert_allclose(probabilities, [1, 0] if character == 1 else [0, 1], atol=2e-9, rtol=0)


def test_receipt_comparison_rejects_omitted_outputs_and_changed_memory_law():
    from .packet import compare
    value = {'raw': [0, 1, 2], 'memory': 'coherent', 'charge': True}
    compare(value, deepcopy(value))
    for field, replacement in [('raw', [0, 1]), ('memory', 'reset'), ('charge', 1)]:
        wrong = deepcopy(value)
        wrong[field] = replacement
        with pytest.raises(ValueError):
            compare(value, wrong)


@pytest.mark.parametrize('chart', [(0,)*65, (False,)*66, [0]*66, (2,)*66])
def test_chart_rejects_ambiguous_or_incomplete_bits(chart):
    with pytest.raises(ValueError):
        p.lifted_pairs(p.branches('C')[0], chart)
