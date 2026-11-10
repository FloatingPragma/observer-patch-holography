"""Operational checks of the finite coherent observer and its boundaries."""
from fractions import Fraction as F
import json
from pathlib import Path
import pytest
import sympy as s

from .common import local_pump_kraus, local_repair_unitary
from .program import Observer, State
from .session import Session


def test_repair_keeps_the_removed_information_in_the_record():
    U = local_repair_unitary()
    assert U.H * U == s.eye(24)
    for record, expected in enumerate(local_pump_kraus()):
        block = s.Matrix(6, 6, lambda i, j: U[4*i+record, 4*j])
        assert s.simplify(block-expected) == s.zeros(6)
    assert sum((K.H*K for K in local_pump_kraus()), s.zeros(6)) == s.eye(6)
    # An actual double load splits with a readable origin, and the inverse
    # restores its input rather than irreversibly erasing that origin.
    for source, record in ((4, 2), (5, 3)):
        incoming = s.eye(24)[:, 4*source]
        outgoing = s.eye(24)[:, 4*3+record]
        assert U*incoming == outgoing and U*outgoing == incoming


def test_idle_walk_is_reversible_but_changes_two_event_populations():
    observer = Observer()
    initial = observer.embedded_port(0)
    first = observer.step(initial)
    second = observer.step(first)
    assert observer.inverse_step(first) == initial
    assert first.norm_squared() == second.norm_squared() == 1
    assert observer.source_populations(first)[0] == s.Rational(11, 12)
    assert observer.source_populations(second)[0] == s.Rational(643, 900)
    assert observer.source_populations(second)[0] != s.Rational(101, 120)


def test_optional_complex_probe_and_filter_keep_all_outcomes():
    observer = Observer()
    assert observer.basis_state(0, 0, 0).norm_squared() == 1
    initial = observer.embedded_port(0)
    phased, label = observer.quarter_turn(initial, observer.index(0, 0, 0))
    assert label == {'kind': 'quarter_turn', 'address': 0}
    assert phased.norm_squared() == 1
    readout = observer.register_measurement(phased, 0)
    assert [record['outcome'] for record, _ in readout] == [True, False]
    assert [state.norm_squared() for _, state in readout] == [s.Rational(1, 60), s.Rational(59, 60)]
    filtered = observer.preparation_filter(phased)
    assert [record['outcome'] for record, _ in filtered] == ['success', 'failure']
    assert s.simplify(sum(state.norm_squared() for _, state in filtered)) == 1
    assert all(state.norm_squared() > 0 for _, state in filtered)


def test_central_holonomy_changes_a_complete_executed_return_detector():
    baseline = Observer(fibre=5)
    edge = baseline.seams[0]
    flux = Observer(fibre=5, central_edge=edge)
    probabilities = []
    for observer in (baseline, flux):
        state = observer.embedded_port(edge[0], internal=3)
        for _ in range(3):
            state = observer.step(state)
        branches = observer.embedded_return_measurement(state, edge[0], internal=3)
        assert [record['outcome'] for record, _ in branches] == [True, False]
        assert sum(branch.norm_squared() for _, branch in branches) == 1
        probabilities.append(branches[0][1].norm_squared())
        # The idle program keeps this state in its declared internal charge.
        for port in range(12):
            for parity in (0, 1):
                for address in range(60):
                    assert all(state.amplitudes[observer.index(port, parity, address, j)] == 0 for j in range(3))
    assert probabilities[1]-probabilities[0] == -s.Rational(9353, 91125000)
    # A deliberately removed background recovers the flat experiment.
    assert Observer(fibre=5).central_edge is None


def test_declared_domain_and_exact_arithmetic_are_enforced():
    with pytest.raises(ValueError):
        Observer(fibre=2)
    with pytest.raises(ValueError):
        Observer(central_edge=(0, 1))
    with pytest.raises(ValueError):
        State((0.5,))
    with pytest.raises(ValueError):
        Observer().step(State((F(1),)))
    with pytest.raises(ValueError):
        Observer().index(0, 0, True)


def test_receipt_does_not_promote_static_flux_to_particle_mass():
    receipt = json.loads(Path(__file__).with_name('results.json').read_text())
    assert receipt['scope']['kinematic_gauge_coupling']
    assert not receipt['scope']['physical_particle_identification']
    assert not receipt['scope']['GeV_mass_prediction']
    assert not receipt['scope']['dynamical_gauge_field']
    assert not receipt['scope']['full_7200_dimensional_tomography']
    assert receipt['gauged_walk']['parallel_sections']['general_single_negative_edge_bound'] == 'lambda_min(L_signed) <= 4/N'


def test_bounded_observer_retains_both_histories_and_stops_at_capacity():
    observer = Observer()
    session = Session(observer, observer.embedded_port(0), max_events=2)
    advanced = session.advance()
    assert advanced.wave_steps == 1 and len(advanced.history) == 1
    outcomes = advanced.detect(observer.index(0, 0, 0))
    assert len(outcomes) == 2
    assert sum(branch.state.norm_squared() for branch in outcomes) == 1
    assert {dict(branch.history[-1])['outcome'] for branch in outcomes} == {True, False}
    assert all(branch.wave_steps == 1 and len(branch.history) == 2 for branch in outcomes)
    for branch in outcomes:
        with pytest.raises(OverflowError):
            branch.advance()
    assert session.history == ()
