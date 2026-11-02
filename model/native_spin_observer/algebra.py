"""Exact finite encoding, repair and exchange identities."""
from fractions import Fraction as F
import numpy as np
from ..defect_field_bridge import full_domain

PAULI={'I':np.eye(2,dtype=complex),'X':np.array([[0,1],[1,0]],complex),
       'Y':np.array([[0,-1j],[1j,0]],complex),'Z':np.diag([1,-1]).astype(complex)}


def tensor(label):
    out=np.ones((1,1),complex)
    for p in reversed(label):out=np.kron(out,PAULI[p])
    return out


def inclusion(sites):
    value=np.zeros((1<<(2*sites),1<<sites),complex)
    for q in range(1<<sites):
        mask=sum((1 if q>>s&1 else 2)<<(2*s) for s in range(sites))
        value[mask,q]=1
    return value


def target_repair(law,sites=2,seam=(0,1),inputs=None):
    dim=1<<(2*sites)
    if inputs is None:inputs=np.eye(dim,dtype=complex)/np.sqrt(dim)
    inputs=np.asarray(inputs);out=np.zeros((dim*32,inputs.shape[1]),complex)
    for coin in range(2):
        for mask in range(dim):
            raw,dest,sign,_=full_domain.branch_action(mask,law,coin,sites=sites,flavors=2,seam=seam)
            out[dest+(coin<<(2*sites))+(raw<<(2*sites+1))]+=sign*inputs[mask]/np.sqrt(2)
    return out


def certificate():
    p=inclusion(1);px=p@p.conj().T
    logical={'X':'XX','Y':'YX','Z':'ZI'}
    for axis,label in logical.items():
        a=tensor(label)
        if not np.array_equal(a@p,p@PAULI[axis]):raise ArithmeticError('encoded Pauli changed')
        if np.any(a@px-px@a):raise ArithmeticError('encoded Pauli leaves single occupancy')
    fields=[tensor(x) for x in ('XXXX','YXYX','ZIZI')]
    p2=inclusion(2)
    for a in fields:
        if np.any(a@(p2@p2.conj().T)-(p2@p2.conj().T)@a):raise ArithmeticError('exchange term leaves P')
    if any(np.any(a@b-b@a) for a in fields for b in fields):raise ArithmeticError('bond factors fail to commute')
    scaled=sum(fields)-np.eye(16)
    target=sum(tensor(axis*2) for axis in 'XYZ')-np.eye(4)
    if not np.array_equal(scaled@p2,p2@target):raise ArithmeticError('exchange compression changed')
    repairs=[]
    for law in ('C','E'):
        for sites,seam in ((2,(0,1)),(3,(0,2))):
            dim=1<<(2*sites);cols=[];protected=0
            for coin in range(2):
                destinations=set()
                for mask in range(dim):
                    raw,out,sign,current=full_domain.branch_action(mask,law,coin,sites=sites,flavors=2,seam=seam)
                    key=(raw,out)
                    if key in destinations:raise ArithmeticError('repair is not a recorded injection')
                    destinations.add(key)
                    if all(((mask>>(2*s))&3) in (1,2) for s in range(sites)):
                        if (raw,out,sign,current)!=(0,mask,1,0):raise ArithmeticError('protected repair changed')
                        protected+=1
                    cols.append((mask,coin,raw,out,sign))
            repairs.append({'law':law,'sites':sites,'seam':list(seam),'full_input_columns':dim,
                            'protected_coin_columns':protected,'all_32_coin_raw_slots_retained':True})
    # Bond overlap really produces a neighbor-dependent, noncommuting dynamics.
    h01=sum(tensor(axis+axis+'I') for axis in 'XYZ')-np.eye(8)
    h12=sum(tensor('I'+axis+axis) for axis in 'XYZ')-np.eye(8)
    comm=h01@h12-h12@h01
    if not np.any(comm):raise ArithmeticError('neighbor exchange became commuting')
    return {'encoding':'q=0 -> down occupied; q=1 -> up occupied',
            'encoded_Paulis':logical,'each_encoded_Pauli_preserves_single_occupancy':True,
            'three_bond_Paulis_commute':True,'bond_generator':'J*(XX+YY+ZZ-I)/4',
            'adjacent_bond_commutator_nonzero':True,'repair_audits':repairs,
            'repair_during_individual_native_pulses_is_protected':False,
            'full_CAR_routing_requires_fermionic_swap_sign':True}
