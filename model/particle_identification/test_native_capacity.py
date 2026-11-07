"""Capacity changes are checked on full native instruments, not marginals."""
from copy import deepcopy
from fractions import Fraction

import pytest

from . import native_capacity as bridge

p = bridge.source


@pytest.mark.parametrize('shell', ['H66', 'H78'])
@pytest.mark.parametrize('law', ['C', 'E', 'mixture', 'alternating'])
def test_complete_adaptive_instrument_with_reference(shell, law):
    value, _ = bridge.fixture()
    value = p.program(shell=shell, commands=[
        dict(value['commands'][0], law=law),
        dict(value['commands'][-1], law=law)])
    # Retain a nonzero coherence with an external reference.
    rho = {(i, i, j, j): Fraction(1, 2) for i in (0, 1) for j in (0, 1)}
    baseline = p.execute(value, rho)
    for capacity in (2, 8):
        actual = p.execute(bridge.with_capacity(value, capacity), rho)
        assert bridge.complete_output(actual) == bridge.complete_output(baseline)
        assert actual.branches == tuple(bridge.resegment(b, capacity) for b in baseline.branches)
    assert p.reference_marginal(b.density for b in baseline.branches) == p.reference_marginal([rho])


@pytest.mark.parametrize('selector', [
    {'kind': 'uniform'},
    {'kind': 'fixed', 'seam': [0, 1]},
    {'kind': 'last_current', 'positive': [0, 5], 'otherwise': [2, 3]},
    {'kind': 'observer', 'positive': [0, 5], 'otherwise': [2, 3]},
])
def test_old_tape_does_not_change_complete_one_step_kernel(selector):
    value, rho = bridge.fixture()
    command = dict(value['commands'][0], selector=selector)
    positive = p.Record('repair', current=1)
    old = (p.Record('repair', current=-1), p.Record('hold'), positive)
    small = p._outcomes(value, command, rho, old[-1:], old[:-1])
    large_value = bridge.with_capacity(value, 8)
    large = p._outcomes(large_value, command, rho, old, ())
    assert tuple(bridge.resegment(b, 8) for b in small) == large


@pytest.mark.parametrize('limits', [{'max_steps': 1}, {'max_branches': 25}])
def test_resource_stop_at_same_prefix(limits):
    value, rho = bridge.fixture()
    small = p.execute(value, rho, **limits)
    large = p.execute(bridge.with_capacity(value, 8), rho, **limits)
    assert small.status == large.status == 'RESOURCE_STOP'
    assert bridge.complete_output(small) == bridge.complete_output(large)


def test_hold_resets_latest_current_even_with_long_archive():
    value, rho = bridge.fixture()
    command = value['commands'][-1]
    large_value = bridge.with_capacity(value, 8)
    old = (p.Record('repair', current=1), p.Record('hold'))
    branches = p._outcomes(large_value, command, rho, old, ())
    assert all(b.record.seam == (2, 3) for b in branches if b.record.kind == 'repair')
    positive = p._outcomes(large_value, command, rho, old[:1], ())
    assert all(b.record.seam == (0, 5) for b in positive if b.record.kind == 'repair')


def test_halt_policy_is_separated_by_observable_status():
    _, rho = bridge.fixture()
    small = p.program(commands=[{'op': 'hold', 'repeat': 2}], overflow='halt')
    large = p.program(commands=small['commands'], capacity=2, overflow='halt')
    assert p.execute(small, rho).status == 'FIFO_TIMEOUT'
    assert p.execute(large, rho).status == 'COMPLETE'
    with pytest.raises(ValueError, match='archive'):
        bridge.with_capacity(small, 2)


@pytest.mark.parametrize('capacity', [0, 9, True, 1.0, '2'])
def test_capacity_change_rejects_unadmitted_values(capacity):
    value, _ = bridge.fixture()
    original = deepcopy(value)
    with pytest.raises(ValueError):
        bridge.with_capacity(value, capacity)
    assert value == original


def test_certificate_includes_all_blocks_and_distinct_costs():
    report = bridge.certificate()
    assert report['fixture_complete_raw_blocks'] == 625
    assert report['fixture_total_probability'] == '1'
    assert report['physical_global_capacity_law_derived'] is False
    value, _ = bridge.fixture()
    original = deepcopy(value)
    changed = bridge.with_capacity(value, 8)
    assert value == original
    assert value['declared_resources']['FIFO_slots'] == 1
    assert changed['declared_resources']['FIFO_slots'] == 8
    assert p.digest(value) != p.digest(changed)
