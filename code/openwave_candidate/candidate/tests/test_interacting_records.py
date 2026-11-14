import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.quantum_records import JointRecords
from oph_openwave.interacting_records import InteractingRecords,run
from oph_openwave.runner import validate_result


def test_selected_channel_energy_ledger_and_free_limit():
    result=run();validate_result('interacting_records',result)
    assert all(result['checks'].values()) and all(result['mutations'].values())


def test_energy_accounting_reconstructs_complete_instrument_independently():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    j=InteractingRecords(c,energy_kind='massive',energy_coupling=.29,energy_mass=.23)
    psi=np.arange(1,len(j.states)+1,dtype=complex)*(1+1j);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    output,ledger,account=j.step_with_energy(rho,.37,.29)
    h=j.hamiltonian
    assert np.isclose(sum(row['weighted_output_energy'] for row in ledger),np.trace(output@h).real)
    assert abs(account['coherent_energy_change'])<1e-12
    for row in ledger:
        if row['status']!='no_commit':
            assert np.isclose(row['weighted_output_energy'],row['probability']*h[row['target'],row['target']])
    assert all(j.gauss(s)==j.gauss(j.states[0]) for s in j.states)
    assert np.linalg.eigvalsh(output).min()>-1e-12


def test_selected_generator_validates_energy_configuration():
    c=Construction(2,((0,1),),.23)
    for coupling in (-.1,float('nan'),float('inf')):
        with pytest.raises(ValueError):InteractingRecords(c,energy_coupling=coupling)
    with pytest.raises(ValueError):InteractingRecords(c,energy_kind='missing')
    with pytest.raises(ValueError):InteractingRecords(c,energy_mass=-.1)


def test_observer_sparse_evolution_uses_selected_interaction_once():
    from oph_openwave.observer_protocol import ObserverProtocol
    from oph_openwave.record_preparation import finite_evolve
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    selected=InteractingRecords(c,energy_kind='massive',energy_coupling=.29,energy_mass=.23)
    observer=ObserverProtocol(selected)
    assert np.allclose(observer.h.toarray(),selected.hamiltonian)
    initial=np.zeros(len(selected.states),complex);initial[0]=1
    assert np.allclose(observer.evolve(initial,.37),finite_evolve(selected.hamiltonian,initial,.37))
    bare=JointRecords.hamiltonian.fget(selected)
    assert np.allclose(np.diag(observer.h.toarray()-bare),selected.coherent_diagonal_addition)
    assert not np.allclose(observer.h.toarray(),bare+2*np.diag(selected.coherent_diagonal_addition))
