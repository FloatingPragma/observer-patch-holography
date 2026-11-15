import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.pinned_charge import PinnedCharge,charge_energy,run
from oph_openwave.runner import validate_result

def test_charge_functional_source_and_uniform_intertwiner():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    u,charges,kernel=charge_energy(c)
    assert u[0]==0
    assert np.allclose(charges.sum(axis=1),0)
    assert np.allclose(kernel@np.ones(3),0)
    m=PinnedCharge(c)
    assert np.allclose(m.hamiltonian-c.hamiltonian,.29*np.diag(u))

def test_complete_canonical_component_and_rejecting_controls():
    result=run()
    validate_result('pinned_charge',result)
    assert all(result['checks'].values())
    assert all(result['mutations'].values())
    assert result['full_criterion_met'] is False
    assert result['local_literal_outcome']=='NOTTESTABLE'
    b=result['branches'][0]
    assert b['dimension']==87 and len(b['all_eigenbranches'])==87
    assert len(b['selected_complete_position_amplitudes'])==87
    assert abs(sum(x['probability'] for x in b['monitor_complete_outcomes'])-1)<1e-12
    assert b['counts']['negative']==63

def test_larger_flag_is_explicit_boolean():
    with pytest.raises(ValueError):run(include_larger=1)
