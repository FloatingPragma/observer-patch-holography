import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.quantum_records import JointRecords
from oph_openwave.record_preparation import run as preparations
from oph_openwave.record_energy import RecordEnergy,run as energies


def test_modular_intertwiner_does_not_transfer_to_bounded_integer_state():
    result=preparations()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    from oph_openwave.runner import validate_result
    validate_result("record_operator",result)
    assert result['measurements']['modular_sector_states']==18


def test_spectral_energy_identity_and_integer_vs_real_control():
    result=energies()
    assert all(result['checks'].values()) and all(result['mutations'].values())
    from oph_openwave.runner import validate_result
    validate_result("record_operator",result)


def test_committed_current_fields_use_actual_record_energy_identity():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    joint=JointRecords(c)
    for kind in ('first_order','quadratic','massive'):
        energy=RecordEnergy(c,kind,.23)
        for pos,f in joint.states:
            if pos<c.ports:
                e=np.array(f[::2]);assert np.array_equal(e,-np.array(f[1::2]))
                assert energy.components(e)['identity_error']<1e-12
        h=energy.augment(joint,.29)
        assert np.allclose(h,h.T)
    with pytest.raises(ValueError):RecordEnergy(c,'unsupported')
    with pytest.raises(ValueError):RecordEnergy(c,'massive',float('nan'))
    with pytest.raises(ValueError):RecordEnergy(c).augment(joint,-.1)


def test_modular_control_has_exact_gauss_sector_and_integer_cutoff_does_not_wrap():
    from oph_openwave.record_preparation import modular_control
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    states,h,v=modular_control(c,3)
    subdivided=np.zeros((6,6),dtype=int)
    for edge,(a,b) in enumerate(c.seams):
        subdivided[a,2*edge]=1;subdivided[b,2*edge+1]=1
        subdivided[3+edge,2*edge:2*edge+2]=-1
    initial=np.array([1,0,0,0,0,0])
    for pos,f in states:
        occupancy=np.eye(6,dtype=int)[pos]
        assert np.array_equal((occupancy+subdivided@f)%3,initial)
    assert np.allclose(h,h.T)
    integer=JointRecords(c,cutoff=1)
    assert len(integer.states)==13 and len(states)==18
