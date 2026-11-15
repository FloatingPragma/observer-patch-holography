import json
from pathlib import Path
import numpy as np
from interacting_modular_observer import compare_scientific

def test_actual_interacting_source_and_complete_continuations():
    here=Path(__file__).parent;r=json.loads((here/'INTERACTING_MODULAR_OBSERVER.json').read_text())
    assert all(type(x) is bool and x for group in ('checks','mutations') for x in r[group].values())
    assert r['construction_id']=='triangle-CAR-allN-Z6-periodic-half-record-native-MIX-v1'
    assert len(r['complete_observer_records'])==5
    assert abs(sum(x['probability'] for x in r['complete_observer_records'])-1)<1e-12
    with np.load(here/r['state_archive']['filename'],allow_pickle=False) as z:
        for record in r['complete_observer_records']:
            post=z[record['poststate_archive_key']];out=z[record['continued_archive_key']]
            assert abs(np.trace(post)-record['probability'])<1e-11
            assert abs(np.trace(out)-record['probability'])<1e-11
            assert abs(sum(x['joint_probability'] for x in record['complete_native_continuation'])-record['probability'])<1e-11
    assert not r['full_criterion_met']

def test_receipt_replay_rejects_lost_complement_and_bool_promotion():
    assert not compare_scientific({'probability':.988},{'probability':0.})
    assert not compare_scientific({'full_criterion_met':False},{'full_criterion_met':0})
    assert not compare_scientific({'indices':[43,66]},{'indices':[43,True]})


def test_declared_modular_and_native_parameters_replayed_exactly():
    assert not compare_scientific({'s':.37},{'s':.37+1e-12})
    assert not compare_scientific({'supplied_native_duration':.37},{'supplied_native_duration':.37+1e-12})
