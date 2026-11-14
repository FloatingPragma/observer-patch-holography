"""Finite preparation controls for J02; modular control is a distinct regulator.

Uniform bounded integer record fibers are not assumed to intertwine with the
bare generator. The positive intertwiner control uses cyclic Z_3 registers,
where uniform fibers are normalizable and shifts are unitary.
"""
from collections import deque
import hashlib
from pathlib import Path
import numpy as np
from .construction import Construction
from .quantum_records import JointRecords


def uniform_fibers(states,positions):
    v=np.zeros((len(states),positions))
    for pos in range(positions):
        indices=[i for i,s in enumerate(states) if s[0]==pos]
        if not indices: raise ValueError('preparation missing a position fiber')
        v[indices,pos]=1/np.sqrt(len(indices))
    return v


def modular_control(c,modulus=3):
    """Complete connected modular sector; no integer-cutoff transfer implied."""
    if type(modulus) is not int or modulus<2:raise ValueError('modulus at least two required')
    start=(0,(0,)*(2*len(c.seams)))
    states=[start]; seen={start}; queue=deque([start])
    def hops(state):
        pos,f=state
        for edge,ends in enumerate(c.seams):
            for side,x in enumerate(ends):
                if pos==x: target,delta=c.ports+edge,1
                elif pos==c.ports+edge:target,delta=x,-1
                else:continue
                record=list(f); h=2*edge+side; record[h]=(record[h]+delta)%modulus
                yield (target,tuple(record)),1 if side==0 else -1
    while queue:
        for target,_ in hops(queue.popleft()):
            if target not in seen:
                if len(states)>=4096:raise ValueError('modular control exceeds state budget')
                seen.add(target);states.append(target);queue.append(target)
    index={s:i for i,s in enumerate(states)}
    h=np.zeros((len(states),)*2)
    for col,state in enumerate(states):
        h[col,col]=c.mass*(1 if state[0]<c.ports else -1)
        for target,sign in hops(state):h[index[target],col]+=sign
    return tuple(states),h,uniform_fibers(states,c.ports+len(c.seams))


def finite_evolve(h,psi,time):
    w,v=np.linalg.eigh(h)
    return v@(np.exp(-1j*time*w)*(v.T@psi))


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    j=JointRecords(c,cutoff=1)
    vb=uniform_fibers(j.states,6)
    states,h,vm=modular_control(c,3)
    bounded_defect=float(np.max(np.abs(j.hamiltonian@vb-vb@c.hamiltonian)))
    modular_defect=float(np.max(np.abs(h@vm-vm@c.hamiltonian)))
    bare=np.eye(6)[0]; coherent=finite_evolve(h,vm@bare,.37)
    base=c.evolve(bare,.37)
    evolved_defect=float(np.max(np.abs(coherent-vm@base)))
    definite=finite_evolve(h,np.eye(len(states))[0],.37)
    def positions(psi):
        return [float(sum(abs(psi[i])**2 for i,s in enumerate(states) if s[0]==p)) for p in range(6)]
    fluxfree=positions(coherent);recorded=positions(definite)
    return {'scope':'J02 finite preparation audit: exact modular control versus bounded integer candidate; no common-preparation transfer.',
            'checks':{'modular_intertwines':modular_defect<1e-12,
                      'modular_preparation_isometry':bool(np.allclose(vm.T@vm,np.eye(6))),
                      'modular_evolution_matches_bare':evolved_defect<1e-12,
                      'bounded_preparation_isometry':bool(np.allclose(vb.T@vb,np.eye(6)))},
            'mutations':{'bounded_uniform_transfer_rejected':bounded_defect>.01,
                         'definite_records_differ_from_fluxfree':max(abs(a-b) for a,b in zip(fluxfree,recorded))>1e-6},
            'full_criterion_met':False,
            'supplied':{'graph':'oriented triangle','mass':.23,'time':.37,
                        'candidate_integer_cutoff':1,'control_modulus':3,
                        'fluxfree_control':'uniform modular Gauss-sector fiber at each position'},
            'measurements':{'bounded_intertwining_max_defect':bounded_defect,
                            'modular_intertwining_max_defect':modular_defect,
                            'modular_evolution_max_defect':evolved_defect,
                            'modular_sector_states':len(states),
                            'modular_fluxfree_position_probabilities':fluxfree,
                            'modular_definite_position_probabilities':recorded},
            'physical_rows':{'3':'preparation and clock transfer prerequisite only',
                             '4':'definite-record preparation cannot inherit bare stability',
                             '30':'bare massive-square identification needs a state-space intertwiner'}}
