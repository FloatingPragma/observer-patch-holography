import sympy as s
import pytest

from .model import certificate,reservoir_swap


@pytest.fixture(scope='module')
def value():
    return certificate()


def test_complete_reservoir_energy_charge_and_work_records(value):
    for experiment in value['experiments']:
        records=experiment['complete_records'];e=s.sympify(experiment['one_particle_energy'])
        assert len(records)==36 and experiment['nonzero_record_count']==6
        assert s.simplify(sum(s.sympify(row['probability']) for row in records))==1
        for row in records:
            assert s.simplify(s.sympify(row['energy'])+row['battery_final_level']*e-4*e)==0
            if row['probability']!='0':
                initial=int(row['setting']=='add')
                assert row['relative_pair_charge']+row['reservoir_out']==initial
        assert any(not row['success'] and row['relative_pair_charge']==0
                   and row['energy']==str(2*e) and row['probability']!='0' for row in records)


def test_source_energies_and_unit_ancestry(value):
    assert [x['one_particle_energy'] for x in value['experiments']]==['2','sqrt(10)']
    assert value['experiments'][0]['negative_band_probe_weight']=='1/2'
    assert all(row['inverted_source_parameter']=='2' for row in value['exact_regulator_relations'])
    assert value['unit_map']['source_energy_unit_in_GeV'] is None
    assert value['scope']['physical_Standard_Model_identification'] is False
    assert value['scope']['interacting_relativistic_mass'] is False


def test_reservoir_is_fermionic_and_complete():
    swap=reservoir_swap()
    assert swap.shape==(32,32)
    # Two exchanged fermions retain the canonical exchange sign.
    assert swap[17,17]==-1
    assert swap.H*swap==s.eye(32)
