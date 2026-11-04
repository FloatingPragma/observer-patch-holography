"""The finite radial join preserves physical scope and complete source custody."""
from copy import deepcopy
import pytest
from tools import strict_json
from . import radial_source_campaign as campaign


def packet():
    return strict_json.load(campaign.ROOT/campaign.PATH)


def test_five_issues_and_five_source_results():
    result=campaign.validate(packet())
    assert result['selected_issues']==[1139,1142,1143,1144,1145]
    assert result['reviewed_components']==5
    assert result['canonical_premises_discharged']==0


@pytest.mark.parametrize('key',list(campaign.INTERFACES))
def test_each_source_domain_is_load_bearing(key):
    value=deepcopy(packet())
    value['conditional_interfaces'][key]=not value['conditional_interfaces'][key]
    with pytest.raises(ValueError,match='domain promoted or erased'):
        campaign.validate(value)


@pytest.mark.parametrize('key',[key for key,value in campaign.JOIN.items() if type(value) is bool])
def test_typed_finite_source_join_rejects_changed_identifications(key):
    value=deepcopy(packet())
    value['supported_source_join'][key]=not value['supported_source_join'][key]
    with pytest.raises(ValueError,match='source join broadened or erased'):
        campaign.validate(value)


@pytest.mark.parametrize('damage',['missing','clock','reference','evidence','state','gate',
                                  'proof','mean','temperature','area','physical_count'])
def test_source_and_reused_readout_contract_fail_closed(damage):
    value=deepcopy(packet())
    if damage=='missing':value['component_results'].pop()
    elif damage=='clock':value['supported_source_join']['clock']='An uncalibrated common physical clock'
    elif damage=='reference':value['composition_boundary']['common_continuum_QFT_Einstein_reference_and_refinement_constructed']=True
    elif damage=='evidence':
        evidence=value['component_results'][0]['evidence'];evidence[next(iter(evidence))]='0'*64
    elif damage=='state':value['component_results'][2]['interface']['state']='An unrelated diagonal four-state reference'
    elif damage=='gate':value['integrated_gate_accepted']=True
    elif damage=='proof':value['component_results'][2]['evidence'].pop('model/parity_entropy_response/PROOF.md')
    elif damage=='mean':value['supported_source_join']['finite_reused_readout']='The reused field follows the power of its averaged native mean.'
    elif damage=='temperature':value['supported_source_join']['interaction_to_entropy']='The affine state is the full Gibbs state at the regional contact temperature.'
    elif damage=='area':value['supported_source_join']['regional_entropy_is_geometric_area']=True
    else:value['physical_comparison_count_added']=False
    with pytest.raises(ValueError):campaign.validate(value)


@pytest.mark.parametrize('domain',['particle_domain','geometry_domain','extra_unreviewed_domain'])
def test_domain_prose_and_key_set_are_reviewed(domain):
    value=deepcopy(packet())
    value['conditional_interfaces'][domain]='A complete identified relativistic universe and Einstein metric.'
    with pytest.raises(ValueError,match='source domain changed'):
        campaign.validate(value)
