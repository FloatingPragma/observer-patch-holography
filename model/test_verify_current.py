from copy import deepcopy
import pytest
from . import verify,verify_current as current


def test_exact_three_tuples_replaced_and_all_other_positions_preserved():
    rows=current.current_replays();old=verify.REPLAYS
    assert len(rows)==len(old)
    assert [(a,b)for a,b in zip(old,rows)if a!=b]==[(a,current.ROUTES[a])for a in old if a in current.ROUTES]
    assert len(current.ROUTES)==3
    assert current.LEGACY in verify.REPLAYS
    assert current.CURRENT not in verify.REPLAYS


def test_unlisted_route_deletion_addition_reorder_and_mutation_fail():
    base=list(verify.REPLAYS)
    variants=[base[:-1],base+[('model.extra',None,'extra.json')],base[1:]+base[:1]]
    changed=base.copy();changed[0]=(base[0][0],base[0][1],'other.json');variants.append(changed)
    for rows in variants:
        with pytest.raises(ValueError,match='unreviewed'):current.current_replays(rows)


def test_every_shard_keeps_original_group_order_and_empirical_batch():
    for count in (1,2,6):
        for index in range(count):
            old=verify.selected_replay_groups(index,count)
            expected=tuple(tuple(current.ROUTES.get(x,x) for x in group)for group in old)
            assert current.selected_groups(index,count)==expected


def test_original_hash_guard_cannot_be_bypassed(monkeypatch):
    calls=[]
    def fail():calls.append('original');raise ValueError('original source pin rejected')
    monkeypatch.setattr(verify,'verify_artifacts',fail)
    monkeypatch.setattr(current,'verify_successor_source_bindings',lambda:calls.append('successor'))
    with pytest.raises(ValueError,match='original source'):current.verify_all_hashes()
    assert calls==['original']


def test_successor_hash_guard_is_additional_and_failure_propagates(monkeypatch):
    calls=[]
    monkeypatch.setattr(verify,'verify_artifacts',lambda:calls.append('original')or 17)
    def fail():calls.append('successor');raise ValueError('successor source changed')
    monkeypatch.setattr(current,'verify_successor_source_bindings',fail)
    with pytest.raises(ValueError,match='successor source'):current.verify_all_hashes()
    assert calls==['original','successor']
