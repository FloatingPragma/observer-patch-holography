import numpy as np
from character_minima import run


def test_Z6_all_characters_and_complete_native_energy_support():
    r=run();assert all(r['checks'].values())
    assert r['CAR_minimizing_characters']==[3]
    assert r['one_load_lowest_signed_energy_minimizers']==[3]
    assert r['one_load_positive_band_edge_minimizers']==[0]
    for law,data in r['native_ground_energy_readout'].items():
        assert abs(sum(x['probability'] for x in data['all_energy_groups'])-1)<1e-12
        for event in data['complete_joint_native_energy_records']:
            assert abs(sum(event['all_energy_group_probabilities'])-event['event_probability'])<1e-12
        assert data['first_numerically_supported_excitation']>r['global_fixed_number_gap']+.1
        assert data['mean_excitation_work']>0
