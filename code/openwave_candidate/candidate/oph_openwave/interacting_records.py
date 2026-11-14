"""Selected record-energy variants of the existing JointRecords instrument.

The coherent generator includes the supplied record energy. The native monitor
is unchanged and generally exchanges energy with the system. Work below denotes
change in expected selected Hamiltonian energy, without a heat/work split.
"""
from dataclasses import dataclass
import math
import numpy as np
from .quantum_records import JointRecords
from .record_energy import RecordEnergy
from .construction import Construction


@dataclass
class InteractingRecords(JointRecords):
    energy_kind: str = 'first_order'
    energy_coupling: float = 0.
    energy_mass: float = 0.

    def __post_init__(self):
        if not math.isfinite(self.energy_coupling) or self.energy_coupling<0:
            raise ValueError('finite nonnegative energy coupling required')
        super().__post_init__()
        self.record_energy=RecordEnergy(self.construction,self.energy_kind,self.energy_mass)
        self.energy_diagonal=self.record_energy.joint_diagonal(self)

    @property
    def coherent_diagonal_addition(self):
        return self.energy_coupling*self.energy_diagonal

    @property
    def hamiltonian(self):
        bare=JointRecords.hamiltonian.fget(self)
        return bare+np.diag(self.coherent_diagonal_addition)

    def energy_expectation(self,rho):
        rho=self._density(rho)
        return float(np.trace(rho@self.hamiltonian).real)

    def step_with_energy(self,rho,time,commit_probability):
        """Use the existing step and retain exact finite expected energy changes.

        Each commit outcome contributes p_i <target|H|target>. The uncommitted
        branch retains its state. Coherent evolution preserves this H's mean.
        Monitoring work includes measurement backaction and native placement.
        """
        rho=self._density(rho)
        initial=self.energy_expectation(rho)
        coherent,_=JointRecords.step(self,rho,time,0.)
        coherent_energy=self.energy_expectation(coherent)
        output,ledger=JointRecords.step(self,rho,time,commit_probability)
        final=self.energy_expectation(output)
        diagonal=np.diag(self.hamiltonian)
        for outcome in ledger:
            if outcome['status']=='no_commit':
                outcome['weighted_output_energy']=outcome['probability']*coherent_energy
            else:
                outcome['weighted_output_energy']=outcome['probability']*float(diagonal[outcome['target']])
        accounted=sum(row['weighted_output_energy'] for row in ledger)
        return output,ledger,{'initial_energy':initial,'after_coherent_energy':coherent_energy,
                              'final_energy':final,'coherent_energy_change':coherent_energy-initial,
                              'monitoring_energy_change':final-coherent_energy,
                              'total_energy_change':final-initial,
                              'ledger_output_energy':accounted,
                              'ledger_energy_error':abs(accounted-final),
                              'scope':'finite expected Hamiltonian energy; monitor supplies or removes energy, with no conservation assertion'}


def run():
    from .observer_protocol import ObserverProtocol
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    bare=JointRecords(c)
    free=InteractingRecords(c,energy_kind='massive',energy_coupling=0.,energy_mass=.23)
    rho=bare.initial_density();time=.37;probability=.29
    expected,_=bare.step(rho,time,probability)
    zero,_,_=free.step_with_energy(rho,time,probability)
    variants=[]
    for kind in ('first_order','massive'):
        selected=InteractingRecords(c,energy_kind=kind,energy_coupling=.29,energy_mass=.23)
        output,ledger,energy=selected.step_with_energy(selected.initial_density(),time,probability)
        observer=ObserverProtocol(selected)
        actual_h=observer.h.toarray()
        readout=observer.readout(output)
        free_readout=ObserverProtocol(bare).readout(expected)
        variants.append({'energy_kind':kind,'energy_mass':.23,'energy_coupling':.29,
                         'energy_accounting':energy,'readout':readout,
                         'position_change_from_free':float(np.max(np.abs(np.array(readout['position_probabilities'])-
                                                                            free_readout['position_probabilities']))),
                         'observer_generator_max_error':float(np.max(np.abs(actual_h-selected.hamiltonian))),
                         'trace':float(np.trace(output).real),
                         'ledger_probability':sum(row['probability'] for row in ledger),
                         'exact_gauss_sector':all(selected.gauss(s)==selected.gauss(selected.states[0]) for s in selected.states)})
    return {'scope':'Explicit record-energy variants using JointRecords.step; monitored energy exchange retained. Default candidate energy remains unselected.',
            'checks':{'zero_coupling_exact_channel':bool(np.array_equal(zero,expected)),
                      'gauss_sector':all(v['exact_gauss_sector'] for v in variants),
                      'trace_and_ledger':all(abs(v['trace']-1)<1e-12 and abs(v['ledger_probability']-1)<1e-12 for v in variants),
                      'coherent_energy_conservation':all(abs(v['energy_accounting']['coherent_energy_change'])<1e-12 for v in variants),
                      'instrument_energy_accounting':all(v['energy_accounting']['ledger_energy_error']<1e-12 for v in variants),
                      'observer_uses_selected_generator':all(v['observer_generator_max_error']<1e-12 for v in variants)},
            'mutations':{'monitor_energy_conservation_rejected':all(abs(v['energy_accounting']['monitoring_energy_change'])>1e-6 for v in variants),
                         'interaction_affects_positions':all(v['position_change_from_free']>1e-7 for v in variants)},
            'full_criterion_met':False,'variants':variants,
            'supplied':{'graph':'oriented triangle','mass':.23,'cutoff':1,'time':time,
                        'commit_probability':probability,'preparation':'port0, zero records',
                        'transit_record_extension':'(F_ea-F_eb)/2'}}
