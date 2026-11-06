from fractions import Fraction as Q
from itertools import repeat
from unittest.mock import patch

import pytest

from . import certificate as c


def test_source_prefix_law_on_all_78_inputs_is_nonuniform_and_input_independent():
    assert c.source_law() == (c.WEIGHTS,)*78
    assert c.COSTS == tuple(map(len, c.PREFIXES))
    assert c.WEIGHTS == tuple(Q(1, 2**n) for n in c.COSTS)


def test_every_native_prefix_stops_and_retains_its_actual_records():
    for state in (c.born.BASIS[0], c.born.BASIS[20], c.born.BASIS[-1]):
        for label, word in enumerate(c.PREFIXES):
            result = c.prefix(state, word)
            assert result.label == label and len(result.records) == len(word)
            output = state
            for coin in word:
                output, _ = c.born.endpoint(output, c.born.native.E, coin)
            assert result.output == output


def test_replica_agreement_is_counted_before_normalizing():
    for beta in (1, 2, 3, 4):
        accepted, rejected, mean = c.independent_proposal(beta)
        pi, z = c.reference(beta)
        assert accepted == tuple(w**beta for w in c.WEIGHTS)
        assert rejected == 1-z and mean == Q(3*beta, 2)
        assert tuple(a/z for a in accepted) == pi


def test_sequential_replica_controller_needs_no_source_reset():
    initial = c.born.BASIS[-1]
    # First batch labels 0,1 disagree; second batch labels 2,2 agree.
    result = c.execute(initial, (0, 1, 0, 1, 1, 1, 1), 0, beta=2, cap=3)
    assert result.flag == 'accept' and result.register == 2
    assert len(result.batches) == 2 and result.attempts == 7
    assert [tuple(item.label for item in batch) for batch in result.batches] == [(0, 1), (2, 2)]
    output = initial
    for coin in (0, 1, 0, 1, 1, 1, 1):
        output, _ = c.born.endpoint(output, c.born.native.E, coin)
    assert result.output == output


def test_failure_keeps_old_register_and_all_native_records():
    result = c.execute(c.born.BASIS[0], (0, 1, 0, 0, 1, 0), 2, beta=2, cap=2)
    assert result.flag == 'timeout' and result.register == 2
    assert len(result.batches) == 2 and result.attempts == 6
    law = c.complete_register_law(2, 2, 2)
    assert sum(law.values()) == 1
    assert law['timeout', 2] == Q(5, 8)**2


def test_full_bounded_executor_counts_match_the_matrix_and_resource_law():
    for state in (c.born.BASIS[0], c.born.BASIS[-1]):
        for register in range(3):
            assert c.enumerated_controller(state, register) == c.complete_register_law(register, 2, 2)


def test_beta_one_is_an_actual_bounded_exact_heat_bath():
    for initial_register in range(3):
        law = c.complete_register_law(initial_register, 1)
        assert law['timeout', initial_register] == 0
        assert tuple(law['accept', h] for h in range(3)) == c.WEIGHTS
    assert c.kernel(1) == (c.WEIGHTS,)*3
    assert c.expected_attempts(1) == Q(3, 2)


def test_finite_beta_two_sampler_is_not_mistaken_for_a_projector():
    k = c.kernel(2, 3)
    pi, z = c.reference(2)
    assert pi == (Q(2, 3), Q(1, 6), Q(1, 6))
    assert c.entropy.multiply(k, k) != k
    assert pi[0].denominator == 3
    assert c.entropy.push(pi, k) == pi
    for i in range(3):
        assert sum(abs(k[i][j]-pi[j]) for j in range(3))/2 <= (1-z)**3


def test_exact_renewal_cost_includes_failed_comparisons():
    for beta in (1, 2, 3, 8):
        _, z = c.reference(beta)
        for cap in (1, 2, 5):
            assert c.expected_attempts(beta, cap) == sum(Q(3*beta, 2)*(1-z)**k for k in range(cap))
            assert c.expected_attempts(beta, cap) <= 2*beta*cap
        assert c.expected_attempts(beta) == Q(3*beta, 2)/z


def test_same_nondegenerate_energy_identifies_integer_beta():
    references = [c.reference(beta)[0] for beta in range(1, 9)]
    assert len(set(references)) == 8
    for beta, pi in enumerate(references, 1):
        assert pi[0]/pi[1] == 2**beta
        assert pi[1]+pi[2] == Q(2, 2**beta+2)
        assert pi[1]+pi[2] <= Q(2, 2**beta)
        assert all(x > 0 for x in pi)


def test_copied_prefixes_do_not_emulate_independent_replicas():
    assert c.reference(2)[0] != c.WEIGHTS
    assert sum(c.WEIGHTS) == 1


def test_stationarity_fluctuation_and_entropy_split_are_exact():
    for beta in (1, 2, 3, 8):
        for cap in (1, 3):
            result = c.check_kernel(beta, cap)
            assert result['positive_support_fluctuation_average'] == 1
            assert result['projector'] == (beta == 1)


def test_moving_a_source_prefix_outcome_breaks_the_source_certificate():
    original = c.prefix
    def changed(state, word, seam=c.born.native.E):
        item = original(state, word, seam)
        return c.Prefix((item.label+1)%3, item.records, item.output)
    with patch.object(c, 'prefix', side_effect=changed):
        with pytest.raises(ValueError, match='stopping rule'):
            c.source_law()


@pytest.mark.parametrize('bad', [0, True, 2.0, -1])
def test_replica_parameter_is_a_positive_integer(bad):
    with pytest.raises(ValueError):c.reference(bad)
    with pytest.raises(ValueError):c.kernel(1, bad)


def test_no_coin_or_seam_coercion_in_native_execution():
    with pytest.raises(ValueError):c.prefix(c.born.BASIS[0], (True,))
    with pytest.raises(ValueError):c.prefix(c.born.BASIS[0], (1,))
    with pytest.raises(ValueError):c.prefix(c.born.BASIS[0], (0,), (0, 2))
