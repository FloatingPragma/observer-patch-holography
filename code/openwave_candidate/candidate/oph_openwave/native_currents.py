"""Same-H position continuity and complete native outcome/current readouts.

Spectral current readout is declared apparatus; C/E alone does not implement it.
"""
import numpy as np
from .construction import Construction
from .uniform_native import UniformFiberReduction

def regional_current(hamiltonian,positions):
    h=np.asarray(hamiltonian,dtype=complex)
    if h.ndim!=2 or h.shape[0]!=h.shape[1] or not np.isfinite(h).all() or not np.allclose(h,h.conj().T,atol=1e-12):raise ValueError('finite Hermitian H required')
    ids=tuple(positions)
    if len(set(ids))!=len(ids) or any(type(i) is not int or i<0 or i>=len(h) for i in ids):raise ValueError('distinct valid integer position indices required')
    p=np.diag([float(i in ids) for i in range(len(h))])
    return p,1j*(h@p-p@h)

def complete_readout(operator,density):
    values,vectors=np.linalg.eigh(operator)
    records=[];out=np.zeros_like(density,dtype=complex)
    for i,value in enumerate(values):
        v=vectors[:,i];effect=np.outer(v,v.conj());branch=effect@density@effect
        records.append({'basis_index':i,'value':float(value),'probability':float(np.trace(branch).real)})
        out+=branch
    return records,out

def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23);model=UniformFiberReduction(c,6)
    psi=np.array([1,2j,0,1-1j,0,0],dtype=complex);psi/=np.linalg.norm(psi);rho=np.outer(psi,psi.conj())
    p,current=regional_current(model.hamiltonian,[0])
    derivative=-1j*(model.hamiltonian@rho-rho@model.hamiltonian)
    continuity=float(abs(np.trace(p@derivative)-np.trace(current@rho)))
    records,post=complete_readout(current,rho)
    native=[];total=np.zeros_like(rho)
    for meta,k in model.kraus_for(0):
        branch=k@rho@k.conj().T;total+=branch
        native.append({'native_record':meta,'probability':float(np.trace(branch).real),
                       'unnormalized_regional_charge':float(np.trace(p@branch).real)})
    ledger=float(sum(r['unnormalized_regional_charge'] for r in native)-np.trace(p@rho).real)
    # Reject the opposite continuity convention and incomplete event records.
    wrong_sign=float(abs(np.trace(p@derivative)+np.trace(current@rho)))
    omitted=max(r['probability'] for r in native)
    return {'schema':'oph.native.same_source_currents.v1','checks':{'same_H_continuity':continuity<1e-12,
              'current_PVM_complete':abs(sum(r['probability'] for r in records)-1)<1e-12,
              'current_expectation':bool(abs(sum(r['probability']*r['value'] for r in records)-np.trace(current@rho).real)<1e-12),
              'native_records_complete':abs(sum(r['probability'] for r in native)-1)<1e-12,
              'native_charge_ledger':bool(abs(ledger-np.trace(p@(total-rho)).real)<1e-12)},
            'mutations':{'opposite_current_sign_rejected':wrong_sign>1e-3,
                         'omitted_native_outcome_rejected':omitted>1e-3},
            'broken_control_residuals':{'opposite_sign':wrong_sign,'omitted_probability':omitted},
            'continuity_error':continuity,'current_complete_records':records,'native_complete_records':native,
            'native_commit_regional_charge_change':ledger,'coherent_charge_derivative':float(np.trace(current@rho).real),
            'scope':'Whole port/transit position charge and its same-H continuity current. Current eigenbasis PVM is supplied complete apparatus; no electric/spin/SM current identity or C/E synthesis inferred. Native commits add separately retained regional charge changes. Uniform-fiber theorem transfers H and each Kraus, not arbitrary postmeasurement continuation.',
            'full_criterion_met':False,'official_rows_changed':[]}
