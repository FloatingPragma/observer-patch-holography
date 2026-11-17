"""Source contracts cannot be promoted or substituted by planning metadata."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
from tools import build_completion_program as program


@pytest.fixture
def fixture_root(tmp_path):
    names=[program.PROGRAM,program.HYPOTHESIS,program.JOINT_DEPENDENCY_PATH,'model/SPEC.md','claims/ccf_lt2.json']
    graph=json.loads((program.ROOT/program.JOINT_DEPENDENCY_PATH).read_text())
    names+= [row['path'] for n in graph['nodes'] for row in n['evidence']]
    for name in set(names):
        target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((program.ROOT/name).read_bytes())
    return tmp_path


def mutate(root,name,fn):
    path=root/name;data=json.loads(path.read_text());fn(data);path.write_text(json.dumps(data))


def test_current_program_and_generated_view():
    p,h=program.load_program()
    assert p['official_model']=='CCF-LT-2'
    assert (program.ROOT/program.OUTPUT).read_text()==program.render(p,h)
    assert all(s['status']!='physical_establishment' for s in p['stages'])


@pytest.mark.parametrize('damage',['model','status','domain','acceptance','evidence','cycle','history','verdict'])
def test_rejects_contract_corruption(fixture_root,damage):
    def change(p):
        if damage=='model':p['official_model']='different-law'
        if damage=='status':p['stages'][9]['status']='finite_constructed'
        if damage=='domain':p['stages'][0]['domain']='physical universe'
        if damage=='acceptance':p['stages'][9]['acceptance']=[]
        if damage=='evidence':p['stages'][0]['evidence']=['model/SPEC.md']
        if damage=='cycle':p['stages'][0]['requires']=['geometric_gravity']
        if damage=='history':p['predecessor_program']={}
        if damage=='verdict':p['decision']='COMPLETED_DOMAIN'
    mutate(fixture_root,program.PROGRAM,change)
    with pytest.raises(ValueError):program.load_program(fixture_root)


def test_rejects_physical_hypothesis_promotion(fixture_root):
    mutate(fixture_root,program.HYPOTHESIS,lambda h:h.update(einstein_premises_realized=True))
    with pytest.raises(ValueError,match='promotion'):program.load_program(fixture_root)


def test_ledger_binds_actual_law_and_retains_physical_boundary(fixture_root):
    view=program.ledger_view(fixture_root)
    assert 'model/SPEC.md' in view['parents']
    assert view['program']['physical_comparison_count_added']==0
    assert view['hypothesis']['einstein_premises_realized'] is False
