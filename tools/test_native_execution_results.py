from copy import deepcopy
from pathlib import Path
import json,pytest
from tools.native_execution_results import NAMES,project,render

@pytest.fixture(scope='module')
def data():
    root=Path(__file__).resolve().parents[1]
    return {n:json.loads((root/f'model/openwave_{n}/results.json').read_text())for n in NAMES}

def test_native_does_not_mean_physical(data):
    x=project(data)
    assert x['literal_repair_process_executed'] and not x['all_effective_fields_implemented_by_that_process']
    assert not x['new_physical_GeV_rest_masses'] and x['official_openwave_untested_rows']==31
    text='\n'.join(render(x));assert 'author-side reproduction' in text and 'zero validated Nature' in text

@pytest.mark.parametrize('mutation',('coherent','cosmology','official'))
def test_unsupported_promotions_rejected(data,mutation):
    d=deepcopy(data)
    if mutation=='coherent':d['native_geometry_readout']['coherent_extension']['native_repair_alone_implements_coherent_service']=True
    if mutation=='cosmology':d['native_capacity_history']['normalization_symmetries']['H0_or_Lambda_or_primordial_amplitude_selected']=True
    if mutation=='official':d['native_platform_gate']['scope']['new_official_awards']=1
    with pytest.raises(ValueError):project(d)
