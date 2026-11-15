from pathlib import Path
import json,hashlib
import numpy as np


def test_actual_canonical_source_contract_and_compact_branches():
    here=Path(__file__).resolve().parent;base=here.parents[1]
    receipt=json.loads((here/'CANONICAL_CONSTRUCTION.json').read_text());contract=json.loads((here/'CANONICAL_SOURCE_CONTRACT.json').read_text())
    assert all(receipt['checks'].values()) and not receipt['full_criterion_met']
    for relative,pin in contract['source_sha256'].items():assert hashlib.sha256((base/relative).read_bytes()).hexdigest()==pin
    archive=np.load(here/'canonical_states.npz');spectral=receipt['spectral_readout_instrument']
    for label,rows in spectral['complete_branch_records'].items():
        for row in rows:
            coefficient=archive[row['coefficient_archive_key']]
            assert coefficient.shape==(row['degeneracy'],row['degeneracy'])
            assert abs(np.trace(coefficient).real-row['probability'])<1e-12
    branch=spectral['actual_populated_degenerate_continuation'];post=archive[branch['poststate_archive_key']]
    assert abs(np.trace(post).real-branch['branch']['probability'])<1e-12


def test_declared_metadata_exact_and_measured_values_tolerant():
    import copy
    from canonical_construction import compare
    good={"selected":{"time":.37,"event_probability":.29},"checks":{"ok":True},"readouts":{"mean_selected_native_work":.1},"spectral_readout_instrument":{"energy_group_tolerance":1e-9}}
    bad=copy.deepcopy(good);bad["selected"]["time"]+=1e-11;assert not compare(good,bad)
    bad=copy.deepcopy(good);bad["checks"]["ok"]=1;assert not compare(good,bad)
    bad=copy.deepcopy(good);bad["spectral_readout_instrument"]["energy_group_tolerance"]+=1e-12;assert not compare(good,bad)
    close=copy.deepcopy(good);close["readouts"]["mean_selected_native_work"]+=1e-12;assert compare(good,close)
