from fractions import Fraction as F
from math import log
from random import Random

import numpy as np
import pytest

from .channel import ClockPair, LagGrid, PairTape, VACUUM_SLOPE
from .experiment import local_word_counts, verify


def test_integer_source_matches_scalar_and_batch_meters():
    grid = LagGrid(3,6,8)
    for w in (F(1,256),-F(1,256),F(0)):
        source = ClockPair(grid,w)
        times = (-93,-53,-5,0,11,16,59,103)
        probe = np.asarray(times)
        batch = source.entropy_from_tick_differences(probe-11, F(1,1), 0.25)
        for i,t in enumerate(times):
            counts = local_word_counts(source,11,t,1,F(1,4))
            total = sum(counts)
            d = F(counts[0]-counts[1]-counts[2]+counts[3],total)-F(1,2)
            assert d == source.shift(F(11,64),F(t,64),1,1,F(1,4))
            assert min(counts)>0 and counts[0]+counts[1] == counts[2]+counts[3]
            assert counts[0]+counts[2] == counts[1]+counts[3]
            control = local_word_counts(source,11,t,1,F(1,4),shared_lag=False)
            assert control[0]-control[1]-control[2]+control[3] == sum(control)//2
            direct = source.entropy_change(F(11,64),F(t,64),1,1,F(1,4))
            assert batch[i] == pytest.approx(direct, abs=1e-17)
            kl = source.divergence(float(d))
            assert -1e-19 <= kl <= float(d)**2*4/3 + 1e-18
            assert -VACUUM_SLOPE*float(d)-direct == pytest.approx(kl,abs=1e-17)
            H = -sum((q/total)*log(q/total) for q in counts)
            Hvac = -2*((3/8)*log(3/8)+(1/8)*log(1/8))
            assert direct == pytest.approx(H-Hvac, abs=5e-16)


def test_local_tape_uses_own_time_and_replays_each_endpoint():
    source = ClockPair(LagGrid(3,6,8), F(1,256))
    # Offered offset 2 is fixed before choosing either appointment.
    for k in range(16):
        tape = PairTape(source,1,"lag",2,0,Random(k))
        for t in (0, F(1,8), F(-1,8)):
            a = tape.read("a",t)
            assert tape.read("b",t+F(2,8)) == -a
            assert tape.read("a",t) == a
    rng = Random(197)
    seen = set()
    for _ in range(512):
        tape = source.issue(rng,1)
        seen.add(tape.kind)
        assert tape.read("a",0) == tape.read("a",0)
    assert seen == {"base","white","independent","lag"}


def test_held_boost_and_cross_record_controls():
    data = verify()
    read = data["sign_readings"]
    b = read["banks"]
    assert b["ordinary"]["correlation"] < b["attenuated"]["correlation"] < b["vacuum"]["correlation"]
    assert b["independent_tables"]["correlation"] == b["vacuum"]["correlation"]
    trials = read["starts_per_setting"]
    assert (b["independent_tables"]["sample"][0]-b["ordinary"]["sample"][0])/trials > 2*read["frequency_radius"]
    for case in data["tangent_ellipsoid_holdouts"].values():
        assert case["levels"][-1]["absolute_error"] < 4e-5
    fast = data["tangent_ellipsoid_holdouts"]["fast"]
    assert fast["target"]-fast["flat_lag_control"] > 0.9


def test_a_fixed_cohort_can_cover_larger_finite_boosts_by_declaring_a_finer_grid():
    grid = LagGrid(8,254,32)
    source = ClockPair(grid,F(1,1 << 27))
    # |beta|=sqrt(8/9) is the fastest of the ten original timelike observers.
    # Also check a faster, independent pair of directions on the same source.
    for x in (F(47,50), F(19,20), -F(97,100)):
        signal = source.entropy_change(0,x,1)/float(source.weight)
        assert signal == pytest.approx(float(1/(1-x*x)), rel=0.012)
    with pytest.raises(ValueError,match="weight"):
        ClockPair(grid,F(1,256))
