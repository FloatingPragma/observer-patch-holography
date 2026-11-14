"""Canonical C10 composition: one interacting all-number native CAR source."""
from pathlib import Path
import json,hashlib
import numpy as np
from charged_reservoir import BoundarySource,BASE
from modular_field_adapter import ModularField


class JointNativeSource:
    def __init__(self):
        self.source=BoundarySource(6)
        self.field=ModularField(self.source,kappa=.37,face_coefficient=.61)
        self.hamiltonian=self.field.hamiltonian;self.states=self.source.states;self.n=self.source.n
        self.numbers=np.array([mask.bit_count() for mask,_ in self.states])
        self.spectral={}
        for number in range(7):
            ids=np.flatnonzero(self.numbers==number);values,vectors=np.linalg.eigh(self.hamiltonian[np.ix_(ids,ids)])
            self.spectral[number]=(ids,values,vectors)
        self.ground=min(float(e[0]) for _,e,_ in self.spectral.values())

    def reference(self,number=None):
        selected=range(7) if number is None else (number,)
        ground=min(float(self.spectral[q][1][0]) for q in selected);rho=np.zeros((self.n,self.n));multiplicity=0
        for q in selected:
            ids,e,v=self.spectral[q];flags=abs(e-ground)<1e-10
            if np.any(flags):rho[np.ix_(ids,ids)]=v[:,flags]@v[:,flags].T;multiplicity+=int(sum(flags))
        return rho/multiplicity,ground,multiplicity

    def energy_number_readout(self,rho):
        records=[]
        for q,(ids,values,vectors) in self.spectral.items():
            probabilities=np.diag(vectors.T@rho[np.ix_(ids,ids)]@vectors).real;groups=[]
            for energy,prob in zip(values,probabilities):
                if groups and abs(energy-groups[-1]['energy'])<1e-9:groups[-1]['probability']+=float(prob);groups[-1]['degeneracy']+=1
                else:groups.append({'energy':float(energy),'probability':float(prob),'degeneracy':1})
            records.extend({'source_number':q,**g} for g in groups)
        return records

    def response(self,rho,reference_energy):
        result=[];states={};reservoir_energy=.23
        for meta,k in self.source.reservoir_instrument():
            state=k@rho@k.T;prob=float(np.trace(state));energy=float(np.trace(self.hamiltonian@state))
            work=energy-reference_energy*prob+reservoir_energy*(meta['reservoir_final']-meta['reservoir_initial'])*prob
            label=f"reservoir_{meta['reservoir_initial']}_{meta['reservoir_final']}";states[label]=state
            result.append({**meta,'reservoir_preparation_probability':.5,'probability':prob,'initial_source_energy':reference_energy,'unnormalized_final_source_energy':energy,'unnormalized_field_energy':float(np.trace(self.field.field@state)),'unnormalized_Hm_energy':float(np.trace(self.field.source@state)),'unnormalized_apparatus_work':work,'complete_source_number_energy_records':self.energy_number_readout(state),'poststate_archive_key':label})
        return result,states


def execute():
    model=JointNativeSource();rho,e,deg=model.reference();rho3,e3,deg3=model.reference(3)
    references={'global_allN':(rho,e,deg),'fixed_N3':(rho3,e3,deg3)};arrays={'hamiltonian':model.hamiltonian,'occupation_masks':np.array([mask for mask,_ in model.states]),'half_records':np.array([f for _,f in model.states]),'boundary_B':3-model.numbers}
    reports={};complete=True;energy_error=0.;reference_error=0.
    for label,(state,energy,multiplicity) in references.items():
        arrays[label+'_reference']=state;reference_error=max(reference_error,float(np.linalg.norm(model.hamiltonian@state-energy*state)))
        response,poststates=model.response(state,energy)
        for key,value in poststates.items():arrays[label+'_'+key]=value
        for row in response:row['poststate_archive_key']=label+'_'+row['poststate_archive_key']
        native={}
        for law in ('C','E','MIX'):
            post,ledger=model.field.commit(state,law);arrays[label+'_native_'+law]=post
            native[law]={'complete_outcomes':ledger,'work':float(np.trace(model.hamiltonian@(post-state)).real),'final_source_number_energy_records':model.energy_number_readout(post),'poststate_archive_key':label+'_native_'+law}
        complete &= abs(sum(r['probability'] for r in response)-1)<1e-12
        for row in response:
            energy_error=max(energy_error,abs(sum(x['energy']*x['probability'] for x in row['complete_source_number_energy_records'])-row['unnormalized_final_source_energy']))
            complete &= abs(sum(x['probability'] for x in row['complete_source_number_energy_records'])-row['probability'])<1e-12
        reports[label]={'reference_energy':energy,'multiplicity':multiplicity,'reference_number_probabilities':{str(q):float(np.trace(state[np.ix_(ids,ids)])) for q,(ids,_,_) in model.spectral.items()},'reference_archive_key':label+'_reference','complete_reservoir_response':response,'field_energy_before_reservoir':float(np.trace(model.field.field@state)),'field_energy_after_reservoir':sum(row['unnormalized_field_energy'] for row in response),'mean_reservoir_apparatus_work':sum(row['unnormalized_apparatus_work'] for row in response),'native_instruments':native}
    # Charge-conjugation n->1−n,F->−F,B->−B changes fixed G to1−G.
    conjugate_images=sum(((63-mask,tuple((-np.array(f))%6)) in model.source.index) for mask,f in model.states)
    native_complete=max(float(np.max(abs(sum(k.T@k for _,k in ops)-np.eye(model.n)))) for law in model.field.operators_by_law.values() for ops in law)
    reservoir_complete=float(np.max(abs(sum(k.T@k for _,k in model.source.reservoir_instrument())-np.eye(model.n))))
    field_commutator=max(float(np.max(abs(model.field.field@k-k@model.field.field))) for _,k in model.source.reservoir_instrument())
    files=[Path(__file__),Path(__file__).with_name('charged_reservoir.py'),Path(__file__).with_name('modular_field_adapter.py'),BASE/'next/C9/fermionic_shared_records.py']
    receipt={'schema':'oph.c10.joint_native_source.v1','checks':{'source384_reservoir768_predeclared':model.n==384,'same_composed_interacting_H':bool(np.array_equal(model.hamiltonian,model.field.source+model.field.field)),'actual_interacting_references_stationary':bool(reference_error<1e-12),'complete_reservoir_responses':bool(complete),'complete_energy_work_readout':bool(energy_error<1e-12),'all_native_C_E_MIX_complete':bool(native_complete<1e-12),'full_reservoir_instrument_complete':bool(reservoir_complete<1e-12),'fixed_background_naive_charge_conjugation_rejected':bool(conjugate_images==0),'reservoir_pulse_commutes_with_field':bool(field_commutator<1e-12)},'parameters':{'source_dimension':384,'joint_source_reservoir_dimension':768,'one_development_case_only':True,'modulus':6,'source_numbers':list(range(7)),'mass':.23,'kappa':.37,'face_coefficient':.61,'source_H':'CARHm + periodic_half_record_energy + g(2I−W−Wdagger), identical supplied coefficients in allnumber sectors','boundary_B':'exact integer3−Nsource; sourceG=(111000)mod6, reservoirGr=nr−B','boundary_extra_energy':0,'reservoir_energy':.23,'reservoir_initial_settings':'fair0/1, setting recorded','clock':'same native source hopparameter; supplied reservoir apparatuspulse has no derivedduration','energy_group_tolerance':1e-9,'reference_degeneracy_tolerance':1e-10},'actual_all_number_spectra':{str(q):{'dimension':len(ids),'all_levels':values.tolist(),'ground':float(values[0])} for q,(ids,values,_) in model.spectral.items()},'references':reports,'errors':{'reference_eigenstate':reference_error,'energy_readout':energy_error,'native_completeness':native_complete,'reservoir_completeness':reservoir_complete,'reservoir_field_commutator':field_commutator},'conjugation_image_states_in_same_background':conjugate_images,'scope':'One actual interacting allnumber source, common nativeC/E/MIX and boundarydressed chargedreservoir instrument. Every reference and response recomputed from THIS H; no bareC9/fieldfree chargedgap transferred. FixedGaussbackground is not automatically chargeconjugation invariant. Zero boundaryextraenergy and pulsework budget are supplied apparatuschoices, not fullphysicalcost. No species, spin, calibratedmass, radiation, physicalvacuum or continuumclaim.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    return receipt,arrays

if __name__=='__main__':
 receipt,arrays=execute();archive=Path(__file__).with_name('joint_source_states.npz');np.savez_compressed(archive,**arrays)
 receipt['state_archive']={'filename':archive.name,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}
 Path(__file__).with_name('JOINT_SOURCE.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'checks':receipt['checks'],'references':{k:{j:v[j] for j in ('reference_energy','multiplicity','reference_number_probabilities')} for k,v in receipt['references'].items()}}))
