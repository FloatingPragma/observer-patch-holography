"""The live maintenance projection cannot mask scientific changes."""
import copy
import pytest
from tools.renew_source_workflow_receipts import conserved_payload, replay_owner


def test_only_the_named_top_level_custody_fields_are_renewable():
    old = {'schema': 1, 'pins': {'source.py': 'old'},
           'science': {'accepted': False, 'count': 2, 'sha256': 'scientific'}}
    fresh = copy.deepcopy(old)
    fresh['pins']['source.py'] = 'fresh'
    assert conserved_payload(old, fresh, {'pins'})
    for key, value in [('accepted', True), ('count', 2.0), ('sha256', 'changed')]:
        damaged = copy.deepcopy(fresh)
        damaged['science'][key] = value
        with pytest.raises(ValueError, match='scientific payload differs'):
            conserved_payload(old, damaged, {'pins'})


def test_unknown_owner_cannot_choose_an_arbitrary_producer():
    with pytest.raises(ValueError, match='unknown live workflow owner'):
        replay_owner('unreviewed_producer')
