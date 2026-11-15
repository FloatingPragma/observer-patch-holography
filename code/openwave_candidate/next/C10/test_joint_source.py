import numpy as np
from joint_source import execute,JointNativeSource


def test_one_interacting_source_all_responses_and_complete_archive():
    r,states=execute();assert all(r['checks'].values())
    assert r['parameters']['source_dimension']==384
    h=states['hamiltonian']
    for label,ref in r['references'].items():
        assert abs(sum(ref['reference_number_probabilities'].values())-1)<1e-12
        assert abs(ref['field_energy_after_reservoir']-ref['field_energy_before_reservoir'])<1e-12
        for row in ref['complete_reservoir_response']:
            branch=states[row['poststate_archive_key']]
            assert abs(np.trace(branch)-row['probability'])<1e-12
            assert abs(sum(x['probability'] for x in row['complete_source_number_energy_records'])-row['probability'])<1e-12
            assert abs(np.trace(h@branch)-row['unnormalized_final_source_energy'])<1e-12
        for law,readout in ref['native_instruments'].items():
            assert abs(sum(x['probability'] for x in readout['complete_outcomes'])-1)<1e-12
            assert abs(sum(x['energy']*x['probability'] for x in readout['final_source_number_energy_records'])-np.trace(h@states[readout['poststate_archive_key']]))<1e-12
    # The interacting reference is recomputed, rather than importing bareC9 energy.
    assert abs(r['references']['global_allN']['reference_energy']+4.010659107)>1
    assert not r['full_criterion_met']
