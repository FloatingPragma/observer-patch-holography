from copy import deepcopy
from fractions import Fraction as F
import pytest
from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.build()


def test_actual_finite_source_has_positive_complete_readout(result):
    rows=result['finite_source_readout']['rows']
    assert len(rows)==3
    for row in rows:
        lo,hi=map(F,row['ideal_LCU_success_probability_interval'])
        rlo,rhi=map(F,row['all_rejected_outcomes_probability_interval'])
        assert 0<lo<=hi<1
        assert lo+rlo<=1<=hi+rhi
        elo,ehi=map(F,row['conditional_excitation_energy_interval'])
        assert F(3,4)<elo<=ehi<F(9,4)
    assert F(result['finite_native_resource_plan']['conditional']['conditional_trace_error_upper'])<F(1,100)


@pytest.mark.parametrize('damage',['scope','bool','source','energy','extra'])
def test_mutations_are_rejected(result,damage):
    value=deepcopy(result)
    if damage=='scope':value['theorem']['new_native_execution_performed']=True
    elif damage=='bool':value['new_checks_pass']=1
    elif damage=='source':value['source_bindings'][next(iter(value['source_bindings']))]='0'*64
    elif damage=='energy':value['finite_source_readout']['rows'][0]['conditional_excitation_energy_interval'][0]='0'
    else:value['unreviewed_premise']=True
    with pytest.raises(ValueError):packet.validate(value,result)


def test_mid_build_source_change_rejected(monkeypatch):
    before=packet.bindings();calls=0
    def changed():
        nonlocal calls
        calls+=1
        return before if calls==1 else {**before,'unexpected':'0'*64}
    monkeypatch.setattr(packet,'bindings',changed)
    with pytest.raises(RuntimeError,match='source changed'):packet.build()
