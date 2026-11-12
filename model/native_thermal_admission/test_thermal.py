from fractions import Fraction as F

import pytest

from .thermal import classification, drift_obstruction, stationary_push
from model.native_abundance_epochs import model as native


def test_all_recurrent_classes_and_their_exact_native_stationarity():
    result=classification();closed=result['closed']
    assert len(closed)==46 and sum(map(len,closed))==66
    assert set(result['class_of'])=={s for s in native.STATES if native.extinct_or_frozen(s)}
    for part in closed:
        p={s:F(1,len(part)) for s in part}
        assert stationary_push(p)==p
        for s in part:
            for t in part:
                assert p[s]*native.kernel()[s].get(t,0)==p[t]*native.kernel()[t].get(s,0)


def test_all_positive_temperature_obstruction_sector_by_sector():
    for load in range(16):
        states=[s for s in native.STATES if sum(s[:3])==load]
        matrix={s:native.kernel()[s] for s in states}
        costs={s:native.energy(s) for s in states}
        assert drift_obstruction(matrix,costs)==(0<load<15)
        for s in states:
            for t in matrix[s]:
                if costs[t]<costs[s]:assert matrix[t].get(s,0)==0


def test_support_and_reverse_moves_are_real_premises():
    # The absorbing high-cost state is stationary; full support cannot be omitted.
    state=(2,2,2,3)
    assert stationary_push({state:F(1)})=={state:F(1)}
    assert native.energy(state)>native.energy((2,2,2,0))
    # A true reversible comparison has positive drift at the lower-energy state.
    matrix={0:{0:F(1,2),1:F(1,2)},1:{0:F(1,2),1:F(1,2)}}
    assert not drift_obstruction(matrix,{0:0,1:1})
    matrix[0]={0:F(1)}
    assert drift_obstruction(matrix,{0:0,1:1})
    with pytest.raises(ValueError):drift_obstruction({0:{0:F(1,2)}},{0:0})


def test_complete_attempt_law_includes_all_six_records_on_absorbing_states():
    state=(1,1,1,5)
    records=[native.event(state,s,c) for s in range(3) for c in range(2)]
    assert len(records)==6 and all(t==state and r['no_operation'] for t,r in records)
    assert len({(r['seam'],r['coin']) for _,r in records})==6
