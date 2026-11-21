import math
from fractions import Fraction

import numpy as np
import pytest

from .clock_inference import _lower, decode_elapsed, infer_rate


def _word(phase, starts, local=.7, visibility=.91):
    hand = visibility*np.exp(1j*phase)
    p = local*(1+np.array([hand.real,hand.imag,-hand.real,-hand.imag]))/4
    clicks = np.rint(starts*p).astype(np.int64)
    return [*clicks, int(starts-clicks.sum())]


def test_one_sided_limits_are_conservative_even_at_endpoints():
    tail = Fraction(1, 100)
    for n in (4, 12, 23):
        for k in range(1,n+1):
            p = Fraction(_lower(k,n,tail))
            area = sum(math.comb(n,h)*p**h*(1-p)**(n-h) for h in range(k,n+1))
            assert area <= tail


def test_decode_unlabelled_even_lapse_with_wrapped_station_phases():
    rate = -math.pi/4096
    first = [_word(phase,n) for phase,n in ((-2.9,2**32),(.6,2**31))]
    later = [_word(phase+rate*6172,n) for phase,n in ((-2.9,2**31),(.6,2**33))]
    result = decode_elapsed(first,later,rate=rate,cap=8191,family_size=48,cell_error=2e-12)
    assert result["even_runs_step_two"] == [[6172,6172]]
    assert result["original_start_count_total"] == 2**32+2**33+2*2**31
    assert result["informative_stations"] == 2
    # Reversing chronological count words reverses the increment modulo one turn.
    reverse = decode_elapsed(later,first,rate=rate,cap=8191,family_size=48)
    assert reverse["even_runs_step_two"] == [[8192-6172,8192-6172]]


def test_retains_every_alias_and_untested_odd_lapse():
    a = [_word(.2,10**8)]
    b = [_word(.2+math.pi*.75,10**8)]
    result = decode_elapsed(a,b,rate=math.pi/8,cap=50,even=False)
    assert result["even_runs_step_two"] == [[6,6],[22,22],[38,38]]
    assert result["untested_odd_run_step_two"] == [1,49]
    flat = [_word(0,10**8,visibility=0)]
    for other in ([0,0,0,0,24], flat[0]):
        unknown = decode_elapsed(a,[other],rate=-.04,cap=61)
        assert unknown["even_runs_step_two"] == [[0,60]]
        assert unknown["informative_stations"] == 0
    wrong = decode_elapsed([a[0],a[0]], [b[0],a[0]],rate=math.pi/8,cap=60)
    assert wrong["even_runs_step_two"] == []


def test_rate_needs_only_local_spacings_and_retains_parity_aliases():
    ticks = np.array([0, 1, 2, 6, 9, 33, 66])
    rate = -.2
    records = [[_word(a+(1.1 if t%2 else 0)+t*rate,2**28) for t in ticks]
               for a in (2.6, -2.8)]
    decoded = infer_rate(records,ticks,family_size=24)
    bands = decoded["per_service_intervals"]
    assert all(any(lo <= x <= hi for lo,hi in bands) for x in (rate,rate+math.pi))
    assert not any(lo <= .2 <= hi for lo,hi in bands)
    assert len(bands) == 2
    assert infer_rate(records,ticks+340,family_size=24)["per_service_intervals"] == bands
    assert infer_rate(np.asarray(records)[:,:2],ticks[:2])["per_service_intervals"] == [[-math.pi,math.pi]]


def test_invalid_hand_and_original_start_inputs():
    impossible = [[10**7,10**7,0,0,0]]
    result = decode_elapsed(impossible,impossible,rate=.01)
    assert result["station_incompatible_with_unit_disk"] == 1
    assert not result["even_runs_step_two"]
    with pytest.raises(ValueError, match="miss count"):
        decode_elapsed([[20,10,5,5]],[[20,10,5,5]],rate=.01,cell_error=1e-8)
    with pytest.raises(ValueError, match="enrolled"):
        decode_elapsed(impossible,impossible,rate=.01,family_size=1)
