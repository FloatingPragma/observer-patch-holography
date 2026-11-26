from fractions import Fraction as F
from itertools import permutations, product
import json

import numpy as np
import pytest

from ...memories import Pair, THETA
from ...native_pair import NativePair
from .experiment import HERE, correlated_source, deposit, forgotten_address_control, trace_distance
from .service import SWEEP, bath_branch, native_gate, output_without_conditioning, replacement_horizon, stock_read


@pytest.mark.parametrize("addresses", list(permutations(SWEEP)) + [(0,0,1,2)])
def test_arbitrary_reference_factorizes_on_each_complete_displacement_record(addresses):
    code = NativePair(1, "left")
    w = stock_read("held", "xy", 32)["work"]
    target = code.state(THETA+Pair.from_weighted_work(w).angle_increment)
    source = correlated_source()
    marginal = np.trace(source.reshape(4,8,4,8), axis1=0, axis2=2)
    total = np.zeros_like(source)
    for record in product((0,1), repeat=len(addresses)):
        block = bath_branch(source, addresses, record, w, code)
        reference = np.trace(block.reshape(4,8,4,8), axis1=0, axis2=2)
        assert trace_distance(block, np.kron(target, reference)) < 8e-14
        total += block
    assert trace_distance(total, np.kron(target, marginal)) < 3e-13


def test_every_joint_record_is_archived_and_no_postselection_is_needed():
    archived = json.loads((HERE/"results.json").read_text())
    assert deposit() == archived
    table = np.array(archived["native_banks"]["joint_record_probabilities"])
    assert table.shape == (8,8) and (table > 0).sum() == 16
    assert abs(table.sum()-1) < 5e-12


def test_control_retains_one_incoming_address():
    success, broken = forgotten_address_control()
    assert success["high_low_word_probabilities"] == [[5/8,1/8,1/8,1/8]]*2
    assert broken["high_low_word_probabilities"] == [[1.,0,0,0], [0.25]*4]
    assert broken["distance_between_preparations"] == 0.75


def test_exact_schedule_erases_the_previously_excluded_signed_transfer():
    w = stock_read("blank", "rest", 32)["work"]
    assert w == 0
    code = NativePair(0, "right")
    gate = native_gate(code, w)
    vacuum = code.state(THETA)
    tangent = np.outer(gate[:,3],gate[:,3].conj()) - np.outer(gate[:,0],gate[:,0].conj())
    for transfer in (-1/16, 0, 1/16):
        returned = output_without_conditioning(vacuum + transfer*tangent, SWEEP, w, code)
        assert trace_distance(returned, vacuum) < 4e-13
    assert trace_distance(vacuum+tangent/16, vacuum) == pytest.approx(1/16)


def test_random_scheduler_chosen_horizon_is_sufficient_and_minimal():
    for levels in (1,2,3,6,32):
        m = 7*levels
        r = F(1,32*4**(levels-1))
        goal = r**8/(1+m)**2
        s = replacement_horizon(levels)
        assert 3*m*F(2,3)**s <= goal < 3*m*F(2,3)**(s-1)
    # Direct finite scheduler enumeration verifies its event bound on one memory.
    for ticks in range(2,6):
        all_schedules = list(product(SWEEP, repeat=ticks))
        missed = sum(len(set(s)) < 3 for s in all_schedules)
        assert F(missed, len(all_schedules)) <= 3*F(2,3)**ticks


def test_bath_rejects_invalid_address_record_and_native_code():
    q = np.eye(4)/4
    for addresses, outcomes in [((3,), (0,)), ((1,), (2,)), (SWEEP, (0,0))]:
        with pytest.raises(ValueError):
            bath_branch(q, addresses, outcomes, 0)
    with pytest.raises(ValueError):
        bath_branch(np.eye(7)/7, SWEEP, (0,0,0), 0)
