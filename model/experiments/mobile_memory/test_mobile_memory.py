import numpy as np
from scipy.sparse.csgraph import connected_components

from model.executable_model import carrier
from model.experiments.radiative_repair import source
from .__main__ import run
from .cage_transfer import transfer
from .current import one_charge_read, pair_current_read
from .law import charged_operator, pair_generator


def test_current_read_is_the_oriented_continuity_instrument():
    car = carrier.Carrier("one addressed seam", 3, [(0, 1)])
    law = charged_operator(car).toarray()
    target = np.diag([0, 1, 0])
    current = 1j * (law @ target - target @ law)
    state = np.array([1, 2j, -1], dtype=complex) / np.sqrt(6)
    for charge in (-1, 1):
        records = one_charge_read(state, (0, 1), charge)
        np.testing.assert_allclose(sum(records.values()), state, atol=1e-15)
        expectation = sum(r * np.vdot(branch, branch) for r, branch in records.items())
        np.testing.assert_allclose(expectation, charge * np.vdot(state, current @ state), atol=1e-15)
        reversed_records = one_charge_read(state, (1, 0), charge)
        for record, branch in records.items():
            np.testing.assert_allclose(reversed_records[-record], branch, atol=1e-15)


def test_pair_current_keeps_the_whole_zero_current_sector():
    nports = 3
    charged = np.arange(nports ** 2).reshape(nports, nports) * (1 + 2j)
    neutral = charged + charged.T
    state = np.stack((charged, neutral)).ravel()
    branches = pair_current_read(state, nports, (0, 1))
    assert set(branches) == {-2, -1, 0, 1, 2}
    np.testing.assert_allclose(sum(branches.values()), state, atol=1e-13)
    np.testing.assert_allclose(
        [np.vdot(branch, branch) for branch in branches.values()],
        [np.vdot(branch, state) for branch in branches.values()],
        atol=1e-10,
    )
    for record, branch in branches.items():
        np.testing.assert_array_equal(branch.reshape(2, nports, nports)[1], neutral if record == 0 else 0)


def test_charged_change_keeps_the_native_neutral_block_and_converter():
    car = carrier.Carrier("one addressed seam", 3, [(0, 1), (1, 2)])
    charged_size = car.n ** 2
    before = source.generator(car, "exchange").toarray()
    changed = pair_generator(car).toarray()
    parked_converter = pair_generator(car, blocked_conversion=True).toarray()
    np.testing.assert_array_equal(changed[charged_size:, :], before[charged_size:, :])
    np.testing.assert_array_equal(changed[:, charged_size:], before[:, charged_size:])
    np.testing.assert_array_equal(parked_converter.diagonal(), changed.diagonal())
    np.testing.assert_array_equal(parked_converter[:charged_size, charged_size:], 0)


def test_mobile_memory_records_replay():
    record = run()
    assert record["mobility"]["receiver_gain"] > 0.02
    assert min(
        row["feedback_shuffled_total_variation"]
        for row in record["current_controlled_cage"]["joint_probabilities"]
    ) > 1e-4
    contrast = record["contrast_cage"]
    coupled = np.asarray(contrast["all_link_joint_record"])
    additive = np.asarray(contrast["all_link_additive_record"])
    assert coupled.shape == additive.shape == (2, 20)
    assert coupled[0].sum() <= contrast["all_time_no_pair_bound"]
    assert coupled[1, 1:].sum() > 0.2
    np.testing.assert_array_equal(additive[:, 1:], 0)


def test_carrier_hypotheses_hold_through_all_three_federations():
    for level in range(3):
        car = carrier.federation(level)
        facts = carrier.federation_checks(car)
        assert facts["every_cell_complete"]
        assert facts["external_seams_form_a_matching"]
        assert connected_components(car.adjacency(), directed=False, return_labels=False) == 1


def test_whole_cage_link_controls_on_larger_federations():
    for level in (1, 2):
        enabled = transfer(level=level)
        blocked = transfer(level=level, memory=False)
        np.testing.assert_allclose(enabled["exciton_record"], (0, 0, 0, 1), atol=1e-12)
        np.testing.assert_allclose(blocked["exciton_record"], (1, 0, 0, 0), atol=1e-12)
        np.testing.assert_allclose(enabled["reference_exciton_record"], (1, 0, 0, 0), atol=1e-12)
