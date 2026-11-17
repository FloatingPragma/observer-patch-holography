"""Navigation refresh retains scientific payloads and rejects changed physics."""
from copy import deepcopy
import pytest
from tools import refresh_native_projection as owner


def fixture():
    ledger={'sections':{'mass':{'value':'3/7','status':'conditional'}},
            'native_records':[0,.2,.8],owner.METADATA:{'program':{'official_model':'separate-component'}}}
    completion={'program':{'official_model':'CCF-LT-2'},'verification':'metadata only'}
    return ledger,completion


def test_projection_changes_only_programme_metadata():
    ledger,completion=fixture();before=deepcopy(ledger)
    result=owner.project(ledger,completion)
    assert ledger==before
    assert result[owner.METADATA]==completion
    assert owner.scientific_projection(result)==owner.scientific_projection(ledger)


@pytest.mark.parametrize('damage',['value','status','records','numeric_type','omission','model'])
def test_scientific_mutation_is_rejected(damage):
    ledger,completion=fixture();result=owner.project(ledger,completion)
    if damage=='value':result['sections']['mass']['value']='3/8'
    if damage=='status':result['sections']['mass']['status']='physical'
    if damage=='records':result['native_records']=[0,0,1]
    if damage=='numeric_type':result['native_records'][0]=False
    if damage=='omission':del result['sections']
    if damage=='model':result[owner.METADATA]['program']['official_model']='another-law'
    with pytest.raises(ValueError):owner.validate_projection(ledger,result)
