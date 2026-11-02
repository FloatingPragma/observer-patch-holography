import copy
import pytest
from . import packet


@pytest.fixture(scope='module')
def certificate():
    return packet.build()


def test_current_certificate(certificate):
    packet.validate(certificate)
    assert certificate['full_repository_replay_performed'] is False


@pytest.mark.parametrize('path,value',[
    (('theorem','source_U'),0.0),
    (('theorem','source_valleys_retained'),1),
    (('theorem','actual_cell_norm_used'),1),
    (('theorem','uniform_covariance_gap_assumed'),True),
    (('theorem','Einstein_endpoint_closed'),True),
    (('theorem','single_lattice_family_of_boosted_sharp_algebras_constructed'),True),
    (('bounds','modulus_example','mesh_is_strictly_positive'),0),
    (('algebra','exact_controls','block_Douglas_identity'),False),
    (('custody','predecessor_full_repository_receipts_replayed'),True),
])
def test_scope_and_exact_type_mutations_rejected(certificate,path,value):
    changed=copy.deepcopy(certificate);cursor=changed
    for key in path[:-1]:cursor=cursor[key]
    cursor[path[-1]]=value
    with pytest.raises(ArithmeticError):packet.validate(changed)


def test_current_source_mutation_rejected(certificate):
    changed=copy.deepcopy(certificate)
    key=next(iter(changed['source_bindings']))
    changed['source_bindings'][key]='0'*64
    with pytest.raises(ArithmeticError):packet.validate(changed)


def test_extra_unproved_top_level_claim_rejected(certificate):
    changed=copy.deepcopy(certificate)
    changed['Einstein_dynamics_proved']=True
    with pytest.raises(ArithmeticError,match='root fields'):packet.validate(changed)


def test_source_change_during_build_rejected(monkeypatch):
    original=packet.source_bindings
    calls=0
    def changing():
        nonlocal calls
        calls+=1
        result=original()
        if calls>1:result[next(iter(result))]='0'*64
        return result
    monkeypatch.setattr(packet,'source_bindings',changing)
    with pytest.raises(ArithmeticError,match='source changed'):
        packet.build()
