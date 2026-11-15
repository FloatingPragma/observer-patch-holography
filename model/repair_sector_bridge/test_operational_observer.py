from dataclasses import replace
import pytest
from .operational_observer import Pair,report

def test_actual_source_and_resources():
    r=report();assert all(r['checks'].values())
    assert r['resource_contract']['per_observer_registers']['retained_register_bits_upper']==298
    assert r['resource_contract']['two_observer_retained_bits_upper']==596

def test_halt_missing_address_and_mismatch():
    p=Pair.prepare(0,0);assert p.step(0,0)==p
    with pytest.raises(ValueError):Pair.prepare(0,1).step(None,0)
    p=Pair.prepare(0,2);p=replace(p,right=replace(p.right,initial_index=1))
    q=p.step(0,0);assert q.status=='observer_mismatch' and q.source.halted

def test_durable_readback_and_sequence_failure():
    p=Pair.prepare(0,2).step(0,0)
    assert len(p.readback()['shared'])==1
    assert p.left.consume(4,p.left.tape[0]).status=='sequence_error'

def test_Lean_literal_tables_bound_to_fresh_actual_source_and_control():
    from pathlib import Path
    import re
    from .operational_observer import transition_table
    from .record_observer import Observer
    source=(Path(__file__).resolve().parents[2]/'Lean/Tower/RepairRecordObserver.lean').read_text()
    def read(name):
        body=re.search(r'def '+name+r' : Array \(Array Nat\) := #\[(.*?)\n\]',source,re.S).group(1)
        return [int(x) for x in re.findall(r'\d+',body)]
    rows=transition_table()
    assert read('nextTable')==[r['next'] for r in rows]
    assert read('letterTable')==[r['letter'] for r in rows]
    expected=[]
    for state in range(66):
        for letter in range(120):
            out=Observer.prepare(state,1).consume(0,letter)
            expected.append(out.state_index if out.status=='full' else 66)
    assert read('controlTable')==expected


def test_hidden_resource_and_shared_budget_mutations_rejected():
    from copy import deepcopy
    from .operational_observer import resources,validate_resource_contract
    contract=resources();assert validate_resource_contract(contract)
    for key,value in [('excluded_resources',[]),('two_observer_retained_bits_upper',298),('readback_output_workspace_bits_upper',0)]:
        changed=deepcopy(contract);changed[key]=value
        with pytest.raises(ValueError):validate_resource_contract(changed)


def test_missing_shared_event_or_premature_full_controller_stops_source():
    p=Pair.prepare(0,2)
    for q in (replace(p,source=replace(p.source,attempts=1)),replace(p,left=replace(p.left,status='full'))):
        result=q.step(0,0)
        assert result.status=='observer_mismatch' and result.source.halted
        assert result.source.state_index==q.source.state_index
        assert result.source.attempts==q.source.attempts
