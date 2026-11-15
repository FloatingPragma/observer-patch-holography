"""External Peierls intervention and same-H current on the full C10 source."""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'next/C10'));sys.path.insert(0,str(BASE/'candidate'))
from joint_source import JointNativeSource
from oph_openwave.native_currents import regional_current


class ExternalPhaseProbe:
    def __init__(self,model):
        self.model=model;self.halfhops=model.field.halfhops
        self.occupations=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in model.states])
        self.edges=[(port,3+e) for e,ends in enumerate(model.source.construction.seams) for port in ends]

    def phases(self,a):
        if np.iscomplexobj(a):raise ValueError('real external phases required')
        a=np.asarray(a,float)
        if a.shape!=(6,) or not np.isfinite(a).all():raise ValueError('six finite external link phases required')
        return a

    def hamiltonian(self,a):
        a=self.phases(a);h=self.model.hamiltonian.astype(complex).copy()
        for phase,t in zip(a,self.halfhops):h+=(np.exp(-1j*phase)-1)*t+(np.exp(1j*phase)-1)*t.T
        return h

    def currents(self,a):
        a=self.phases(a)
        return [1j*(np.exp(-1j*phase)*t-np.exp(1j*phase)*t.T) for phase,t in zip(a,self.halfhops)]

    def gauge(self,a,lambdas):
        if np.iscomplexobj(lambdas):raise ValueError('real orbital gauge phases required')
        lambdas=np.asarray(lambdas,float)
        if lambdas.shape!=(6,) or not np.isfinite(lambdas).all():raise ValueError('six finite orbital gauge phases required')
        transformed=self.phases(a)+np.array([lambdas[x]-lambdas[y] for x,y in self.edges])
        return transformed,np.exp(1j*self.occupations@lambdas)

    def native_operators(self,a,law='MIX'):
        if law not in self.model.field.operators_by_law:raise ValueError('C/E/MIX required')
        a=self.phases(a);out=[]
        for seam,ops in enumerate(self.model.field.operators_by_law[law]):
            result=[]
            for meta,k in ops:
                d=meta['record_value'];coin=meta['coin'];j=d if meta['law']=='C' else (d-(2*coin-1)*(abs(d)%2))//2
                result.append(({**meta,'applied_current':j,'proposed_current':j,'half_record_delta':[j,-j]},k*np.exp(-1j*j*(a[2*seam]-a[2*seam+1]))))
            out.append(tuple(result))
        return out

    def evolve(self,rho,a,time=.37):
        if not math.isfinite(time):raise ValueError('finite native time required')
        h=csr_matrix(self.hamiltonian(a));left=expm_multiply(-1j*time*h,rho,traceA=-1j*time*h.diagonal().sum())
        return expm_multiply(1j*time*h.T,left.T,traceA=1j*time*h.diagonal().sum()).T

    def native_step(self,rho,a,law='MIX',probability=.29,region=(0,3)):
        rho=np.asarray(rho,complex)
        if rho.shape!=(self.model.n,self.model.n) or not np.isfinite(rho).all():raise ValueError('finite source matrix required')
        if isinstance(probability,(bool,np.bool_)) or not math.isfinite(probability) or not 0<=probability<=1:raise ValueError('native probability in[0,1] required')
        region=tuple(region)
        if len(set(region))!=len(region) or any(type(z) is not int or not 0<=z<6 for z in region):raise ValueError('distinct whole orbital region required')
        h=self.hamiltonian(a);qvalues=self.occupations[:,list(region)].sum(axis=1);q=np.diag(qvalues);out=(1-probability)*rho;records=[]
        for ops in self.native_operators(a,law):
            for meta,k in ops:
                branch=probability*(k@rho@k.conj().T)/3;out+=branch
                branchprob=float(np.trace(branch).real)
                beforecharge=float(probability*np.dot(qvalues*np.sum(abs(k)**2,axis=0),np.diag(rho).real)/3)
                first,second=self.model.source.construction.seams[meta['seam']]
                expected_delta=(int(second in region)-int(first in region))*meta['applied_current']*branchprob
                records.append({**meta,'conditional_input_regional_charge_weight':beforecharge,'applied_regional_charge_delta':expected_delta,'probability':branchprob,'unnormalized_regional_charge':float(np.trace(q@branch).real),'unnormalized_energy':float(np.trace(h@branch).real)})
        records.append({'status':'no_commit','probability':float((1-probability)*np.trace(rho).real),'unnormalized_regional_charge':float((1-probability)*np.trace(q@rho).real),'unnormalized_energy':float((1-probability)*np.trace(h@rho).real)})
        return out,{'complete_outcomes':records,'charge_before':float(np.trace(q@rho).real),'charge_after':float(np.trace(q@out).real),'energy_before':float(np.trace(h@rho).real),'energy_after':float(np.trace(h@out).real),'mean_native_work':float(np.trace(h@(out-rho)).real)}


def execute():
    model=JointNativeSource();probe=ExternalPhaseProbe(model);rho,e,deg=model.reference()
    a=np.array([0,0,0,0,.17,0.]);h=probe.hamiltonian(a);currents=probe.currents(a);coherent=probe.evolve(rho,a)
    region=(0,3);q=sum((np.diag(probe.occupations[:,z]) for z in region),np.zeros_like(h));regional=sum((regional_current(h,[int(i) for i in np.flatnonzero(probe.occupations[:,z])])[1] for z in region),np.zeros_like(h))
    link=sum(((int(x in region)-int(y in region))*j for (x,y),j in zip(probe.edges,currents)),np.zeros_like(h))
    continuity=float(np.max(abs(regional-link)));derivative=-1j*(h@coherent-coherent@h)
    continuity_expect=float(abs(np.trace(q@derivative)-np.trace(link@coherent)))
    epsilon=1e-6;ap=a.copy();am=a.copy();ap[4]+=epsilon;am[4]-=epsilon
    current_error=float(np.max(abs(currents[4]+(probe.hamiltonian(ap)-probe.hamiltonian(am))/(2*epsilon))))
    gauge_a,phase=probe.gauge(a,[.2,-.37,.11,.43,-.19,.07]);trans=lambda x:phase[:,None]*x*phase.conj()[None,:]
    hg=probe.hamiltonian(gauge_a);rho_g=trans(rho);coherent_g=probe.evolve(rho_g,gauge_a)
    hcov=float(np.max(abs(hg-trans(h))));statecov=float(np.max(abs(coherent_g-trans(coherent))))
    kcov=0.;records_error=0.;postcov=0.;native={};arrays={'initial_actual_reference':rho,'intervened_H':h,'coherent_after_probe':coherent}
    for law in ('C','E','MIX'):
        out,ledger=probe.native_step(coherent,a,law);gout,gledger=probe.native_step(coherent_g,gauge_a,law)
        native[law]=ledger;arrays['native_'+law]=out
        postcov=max(postcov,float(np.max(abs(gout-trans(out)))))
        for old,new in zip(ledger['complete_outcomes'],gledger['complete_outcomes']):
            for key in ('probability','unnormalized_regional_charge','unnormalized_energy'):records_error=max(records_error,abs(old[key]-new[key]))
        for left,right in zip(probe.native_operators(a,law),probe.native_operators(gauge_a,law)):
            for (_,k),(_,kg) in zip(left,right):kcov=max(kcov,float(np.max(abs(kg-trans(k)))))
    wrong=probe.evolve(rho,gauge_a);wrong_prep=float(np.linalg.norm(wrong-coherent_g))
    # Static tiny probe: recompute the actual all-number ground, not a bare gap.
    def ground(a):
        h=probe.hamiltonian(a);lowest=[]
        for number,(ids,_,_) in model.spectral.items():
            values,vectors=np.linalg.eigh(h[np.ix_(ids,ids)]);lowest.append((float(values[0]),number,ids,vectors[:,0]))
        energy,number,ids,vec=min(lowest,key=lambda x:x[0]);state=np.zeros_like(h);state[np.ix_(ids,ids)]=np.outer(vec,vec.conj());return energy,number,state
    eg,ng,rg=ground(a);ep,_,_=ground(ap);em,_,_=ground(am)
    persistent=float(np.trace(currents[4]@rg).real);hferror=abs(persistent+(ep-em)/(2*epsilon))
    sourcefiles=[Path(__file__),BASE/'next/C10/joint_source.py',BASE/'next/C10/charged_reservoir.py',BASE/'next/C10/modular_field_adapter.py',BASE/'candidate/oph_openwave/native_currents.py']
    r={'schema':'oph.c11.external_link_phase_probe.v1','checks':{'zero_external_phase_exact_baseline':bool(np.array_equal(probe.hamiltonian(np.zeros(6)),model.hamiltonian)),'J_is_minus_actual_H_derivative':bool(current_error<1e-9),'regional_continuity_operator':bool(continuity<1e-12),'regional_continuity_expectation':bool(continuity_expect<1e-12),'every_native_applied_current_matches_regional_charge_delta':bool(all(abs(row['unnormalized_regional_charge']-row['conditional_input_regional_charge_weight']-row['applied_regional_charge_delta'])<1e-12 for x in native.values() for row in x['complete_outcomes'] if row.get('status')!='no_commit')),'Hamiltonian_gauge_covariance':bool(hcov<1e-12),'transformed_reference_evolution_covariance':bool(statecov<1e-12),'every_C_E_MIX_Kraus_gauge_covariance':bool(kcov<1e-12),'all_native_outcome_charge_energy_covariance':bool(max(records_error,postcov)<1e-12),'untransformed_preparation_control_rejected':bool(wrong_prep>1e-3),'static_ground_current_Hellmann_Feynman':bool(hferror<1e-8),'nonzero_external_current_response':bool(abs(persistent)>1e-5),'complete_all_native_outcomes':bool(all(abs(sum(row['probability'] for row in x['complete_outcomes'])-1)<1e-12 for x in native.values()))},'parameters':{'construction':'same full384state C10 interacting allN source','external_half_link_phases':a.tolist(),'cycle_flux':float(a@np.array([1,-1]*3)),'time':.37,'native_probability':.29,'selected_native_law':'MIX; C/E separate controls','initial_reference_energy':e,'initial_reference_multiplicity':deg,'region_whole_orbitals':list(region),'reservoir_used':False,'gauge_rule':'G=exp(i sum lambda_z n_z); A_hprime=A_h+lambda_port−lambda_transit','native_phase_rule':'K_record(A)=exp(−i j(A_ea−A_eb)) K_record(0)'},'currents':{'after_quench_coherent': [float(np.trace(j@coherent).real) for j in currents],'static_perturbed_ground_energy':eg,'static_perturbed_ground_number':ng,'static_persistent_current':persistent,'quench_work':float(np.trace((h-model.hamiltonian)@rho).real),'regional_coherent_charge_derivative':float(np.trace(q@derivative).real)},'native':native,'errors':{'current_finite_difference':current_error,'continuity_operator':continuity,'continuity_expectation':continuity_expect,'H_gauge':hcov,'state_gauge':statecov,'Kraus_gauge':kcov,'outcome_gauge':records_error,'native_poststate_gauge':postcov,'untransformed_preparation_distance':wrong_prep,'Hellmann_Feynman':hferror},'scope':'Declared external Peierls phase/current apparatus on the same interacting source. This continuous external gauge covariance does not assert an internal integerF commutator at cyclic wrap or identify a physical magnetic field, magnetostatics, spin or g=2. Per-record native dressing is a scalar phase and cancels in a classical-record branch map; it is retained for exact operator covariance, not misreported as a different classical channel.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles}}
    return r,arrays
if __name__=='__main__':
 r,arrays=execute();path=Path(__file__).with_name('link_phase_states.npz');np.savez_compressed(path,**arrays);r['state_archive']={'filename':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()};Path(__file__).with_name('LINK_PHASE_PROBE.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'currents':r['currents']}))
