"""The official navigation renderer writes prose, without premise promotion."""
from tools import renew_finite_attachment_navigation as owner


def test_only_current_navigation_outputs(monkeypatch,tmp_path):
    outputs={'model/PREMISES.md':'current theorem domains\n','docs/COMPLETION_PROGRAM.md':'current source acceptance\n'}
    monkeypatch.setattr(owner,'build',lambda root:outputs)
    for path in outputs:
        target=tmp_path/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_text('stale\n')
    preview=owner.execute(tmp_path,False)
    assert set(preview['changed_paths'])==set(outputs)
    assert (tmp_path/'model/PREMISES.md').read_text()=='stale\n'
    written=owner.execute(tmp_path,True)
    assert written['physical_premises_discharged']==[]
    assert written['gate_statuses_changed'] is False
    assert all((tmp_path/p).read_text()==value for p,value in outputs.items())
    assert owner.execute(tmp_path,False)['changed_paths']==[]
