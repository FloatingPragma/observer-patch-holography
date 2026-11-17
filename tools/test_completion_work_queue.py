"""Acceptance tasks retain the official source and physical boundary."""
from copy import deepcopy
import pytest
from tools import build_completion_work_queue as queue


def test_official_queue():
    data=queue.load_queue()
    assert data['official_model']=='CCF-LT-2'
    assert all(t['status']=='required' for t in data['tasks'])
    assert 'Same-source calculations' in queue.render_section(data)


@pytest.mark.parametrize('damage',['physics','gate','source','scope','completion','cycle'])
def test_rejects_task_corruption(damage):
    data=deepcopy(queue.load_queue())
    if damage=='physics':data['gate_acceptance_changed']=True
    if damage=='gate':data['tasks'][0]['gate']='unproved'
    if damage=='source':data['tasks'][0]['source_paths']=['../other.py']
    if damage=='scope':data['tasks'][0]['scope_boundary']='all physical states'
    if damage=='completion':data['tasks'][0]['status']='completed'
    if damage=='cycle':data['tasks'][0]['depends_on']=[data['tasks'][0]['id']]
    with pytest.raises(ValueError):queue.validate_queue(data)
