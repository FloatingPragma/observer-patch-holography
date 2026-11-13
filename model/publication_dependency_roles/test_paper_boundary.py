import shutil
import pytest
from pathlib import Path
from model.particle_yang_mills_gap_contract import transfer


def test_presentation_edit_does_not_change_science_and_operative_edit_invalidates_pin(tmp_path,monkeypatch):
    for name in transfer.SOURCE_INPUTS:
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(transfer.ROOT/name,target)
    paper=tmp_path/transfer.PRESENTATION_INPUTS[0]
    paper.parent.mkdir(parents=True,exist_ok=True);paper.write_text('presentation A')
    monkeypatch.setattr(transfer,'ROOT',tmp_path)
    first=transfer.certificate()
    paper.write_text('different presentation B')
    assert transfer.certificate()==first
    operative=tmp_path/'model/particle_yang_mills_gap_contract/transfer.py'
    operative.write_bytes(operative.read_bytes()+b'\n# source mutation\n')
    changed=transfer.certificate()
    assert changed['source_sha256']!=first['source_sha256']
    assert first['source_sha256'].keys()==changed['source_sha256'].keys()


def test_exact_presentation_tuple_does_not_drop_owner_proof_or_Lean():
    assert transfer.PRESENTATION_INPUTS==('extra/yang_mills_gap_clay_problem.tex',)
    assert 'model/particle_yang_mills_gap_contract/PROOF.md' in transfer.SOURCE_INPUTS
    assert 'Lean/ObserverPatchHolography/RepairGapChain.lean' in transfer.SOURCE_INPUTS
    assert 'Lean/Screen/KogutSusskindFiberRateComparison.lean' in transfer.SOURCE_INPUTS


def test_complete_binding_union_hashes_once_and_detects_mutations(tmp_path,monkeypatch):
    """Repeated references never weaken exact operative-byte or conflict checks."""
    import hashlib,json
    from . import paper_check as p
    sources=['model/particle_yang_mills_gap_contract/PROOF.md',
             'Lean/ObserverPatchHolography/RepairGapChain.lean',
             'Lean/Screen/KogutSusskindFiberRateComparison.lean','operative.py']
    for name in sources+[p.PAPER]:
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_text(name+'\n')
    pins={name:hashlib.sha256((tmp_path/name).read_bytes()).hexdigest()for name in sources}
    owners=('one.json','two.json');monkeypatch.setattr(p,'OWNERS',owners)
    for name in owners:(tmp_path/name).write_text(json.dumps({'source_sha256':pins,'nested':[{'input_sha256':pins}]*8}))
    calls=[];original=p.sha
    def counted(path):calls.append(path.relative_to(tmp_path).as_posix());return original(path)
    monkeypatch.setattr(p,'sha',counted)
    result=p.check(tmp_path)
    assert result['scientific_source_pin_count']==4
    assert all(calls.count(name)==1 for name in sources)
    assert sum(r['authenticated_binding_edges']for r in result['accepted_scientific_receipts'])==72
    (tmp_path/'operative.py').write_text('changed operative input\n')
    with pytest.raises(AssertionError,match='scientific source changed'):p.check(tmp_path)
    conflict=dict(pins);conflict['operative.py']='0'*64
    (tmp_path/owners[1]).write_text(json.dumps({'source_sha256':conflict}))
    with pytest.raises(AssertionError,match='conflicting source binding'):p.check(tmp_path)
