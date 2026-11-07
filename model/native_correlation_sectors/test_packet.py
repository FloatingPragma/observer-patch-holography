"""Receipt replay, typed mutations and source-bound native controls."""
from copy import deepcopy
import json
from fractions import Fraction as F

import pytest

from . import packet as p


@pytest.fixture(scope='module')
def candidate():return p.produce()


def test_owner_packet_replays_and_keeps_physical_boundaries(candidate):
    assert p.verify()
    assert candidate['scope']['physical_particle_or_gauge_identification'] is False
    assert candidate['scope']['raw_quantum_branches_intertwine_all_occupations'] is False
    assert candidate['native_attachment']['H66_source_columns']==66
    for row in candidate['evaluated_federations']:
        values=[F(x['normalized_covariance']) for x in row['rows']]
        assert values[0]==1 and all(0<b<=a for a,b in zip(values,values[1:]))


@pytest.mark.parametrize('mutation',('bool','float','norm','branch','hash'))
def test_receipt_rejects_typed_or_scientific_mutations(candidate,tmp_path,mutation):
    value=deepcopy(candidate)
    if mutation=='bool':value['native_attachment']['H66_fair_raw_maps']=True
    elif mutation=='float':value['native_attachment']['H66_source_columns']=66.0
    elif mutation=='norm':value['evaluated_federations'][0]['rows'][0]['normalized_covariance']='2'
    elif mutation=='branch':value['complete_reflected_record_identity']['zero_two_E_histories']-=1
    else:value['source_bindings']['model/native_correlation_sectors/sectors.py']='0'*64
    path=tmp_path/'changed.json';path.write_text(p.canonical(value))
    with pytest.raises(ValueError):p.verify(path)


def test_duplicate_json_keys_rejected_before_replay(candidate,tmp_path):
    path=tmp_path/'duplicate.json';path.write_text('{"schema":"wrong",'+p.canonical(candidate)[1:])
    with pytest.raises(ValueError):p.verify(path)
