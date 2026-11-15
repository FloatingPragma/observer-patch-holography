import numpy as np
import pytest
from oph_openwave.regional_charge import Construction,JointRecords,RegionalCharge,run


def test_complete_operational_charge_variants():
    result=run()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    assert result['physical_species_mapping'] is None
    assert result['full_criterion_met'] is False
    assert all(v['continuity']['rejected_outcomes_retained']>0 for v in result['variants'])


def test_luders_measurement_complement_transit_and_postselection():
    j=JointRecords(Construction(3,((0,1),(1,2),(2,0)),.23))
    rho,_=j.step(j.initial_density(),.37,.29)
    region=RegionalCharge(j,(0,),(0,));output,branches=region.instrument(rho)
    assert len(branches)==2
    assert np.isclose(np.trace(output),1)
    for branch in branches:
        if branch['conditional_density'] is not None:
            r=region.readout(branch['conditional_density'])
            assert np.isclose(r['mean_charge'],branch['sharp_charge'])
            assert np.isclose(r['variance'],0)
        for row in region.readout(rho)['joint_charge_boundary_gauss_probabilities']:
            assert row['charge']+row['boundary_record']==row['gauss']
    complement=RegionalCharge(j,ports=(1,2),transit_seams=(1,2))
    assert np.array_equal(region.charge+complement.charge,np.ones(len(j.states)))


def test_impossible_outcomes_retained_and_no_half_transit_assignment():
    j=JointRecords(Construction(2,((0,1),)))
    empty=RegionalCharge(j)
    output,branches=empty.instrument(j.initial_density(),True)
    assert len(branches)==4
    assert sum(row['probability']==0 for row in branches)==3
    region=RegionalCharge(j,transit_seams=(0,))
    assert set(region.charge)=={0,1}
    with pytest.raises(ValueError):RegionalCharge(j,ports=(.5,))
    with pytest.raises(ValueError):RegionalCharge(j,transit_seams=(1,))
