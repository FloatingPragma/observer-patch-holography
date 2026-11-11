from fractions import Fraction as F

import pytest

from .bridge import casimir_kernel, certificate, lump_kernel, minorization, native_quotient, stochastic, su2_controls


@pytest.mark.parametrize("kernel", [[], [[1, 0]], [[2, -1], [0, 1]], [[F(1, 2), 0], [0, 1]],
                                    [[1.0, 0], [0, 1]], [[True, 0], [0, 1]]])
def test_invalid_kernel_cannot_receive_exact_certificate(kernel):
    with pytest.raises(ValueError):
        certificate(kernel, [1, 1])


@pytest.mark.parametrize("y", [0, -1, F(11, 10), 0.9, True])
def test_acceptance_rejects_invalid_temperature_parameter(y):
    with pytest.raises(ValueError):
        casimir_kernel([[0, 1], [1, 0]], [0, 1], y)


def test_invalid_proposal_is_rejected_before_acceptance_hides_its_excess():
    # Multiplying the upward probability by y would hide an invalid row.
    with pytest.raises(ValueError):
        casimir_kernel([[0, 2], [1, 0]], [0, 1], F(1, 4))


def test_rejection_record_is_required_to_recover_complete_kernel():
    kernel, records = casimir_kernel([[0, 1], [1, 0]], [0, 2], F(1, 2))
    assert kernel == ((F(3, 4), F(1, 4)), (F(1), F(0)))
    assert (0, 1, 0, "rejected", F(3, 4)) in records
    assert sum(row[4] for row in records if row[0] == 0 and row[3] != "rejected") == F(1, 4)
    assert certificate(kernel, [1, F(1, 4)])["unique_stationary_law"]


def test_stationarity_does_not_imply_unique_selection():
    kernel = [[1, 0], [0, 1]]
    for weights in ([1, 1], [1, 9]):
        row = certificate(kernel, weights)
        assert row["stationary"] and row["detailed_balance"]
        assert not row["unique_stationary_law"]
        assert not row["convergence_certified"]


def test_irreducibility_requires_return_paths():
    row = certificate([[0, 1], [0, 1]], [1, 1])
    assert not row["irreducible"]


def test_periodic_chain_has_no_claimed_convergence():
    row = certificate([[0, 1], [1, 0]], [1, 1])
    assert row["unique_stationary_law"]
    assert not row["convergence_certified"]


def test_asymmetric_flux_changes_the_stationary_dimension_power():
    rows = su2_controls()
    for name in ("label_walk", "one_sided", "two_sided_fusion"):
        row = rows[name]
        assert row["stationary"] and row["detailed_balance"] and row["convergence_certified"]
        assert row["record_rows_complete"] and row["record_pushforward_matches_kernel"]
    assert rows["one_sided"]["one_sided_gibbs_is_stationary"]
    assert not rows["two_sided_fusion"]["one_sided_gibbs_is_stationary"]
    assert not rows["label_walk"]["one_sided_gibbs_is_stationary"]


def test_nonlumpable_observer_partition_is_rejected():
    # Representatives 0 and 1 have incompatible probabilities of outcome b.
    with pytest.raises(ValueError, match="not lumpable"):
        lump_kernel([[1, 0, 0], [0, 0, 1], [0, 0, 1]], ["a", "a", "b"])


def test_mixing_bound_requires_actual_stationary_target():
    with pytest.raises(ValueError, match="stationary"):
        minorization([[F(1, 2), F(1, 2)], [1, 0]], [1, 1], 2)
    assert minorization([[0, 1], [1, 0]], [1, 1], 2) == 0


def test_actual_native_events_discharge_finite_quotient_premises():
    row = native_quotient()
    assert row["native_states"] == 66
    assert row["complete_native_records"] == 3960
    assert dict(zip(row["distance_labels"], row["orbit_sizes"])) == {1: 30, 2: 30, 3: 6}
    assert row["native_certificate"]["convergence_certified"]
    assert row["quotient_certificate"]["convergence_certified"]
    probabilities = dict(zip(row["distance_labels"], row["quotient_certificate"]["probabilities"]))
    assert probabilities == {1: "5/11", 2: "5/11", 3: "1/11"}
    assert row["two_attempt_minorization"] == "11/450"
    kernel = stochastic([[F(x) for x in values] for values in row["kernel"]])
    assert not certificate(kernel, [1, 1, 1])["stationary"]
