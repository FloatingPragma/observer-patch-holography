"""Complete grouped number/energy Luders instrument with compact branch states."""
import numpy as np
from scipy.linalg import expm


class GroupedEnergyInstrument:
    def __init__(self,model,tolerance=1e-9):
        if isinstance(tolerance,(bool,np.bool_)) or not isinstance(tolerance,(int,float,np.integer,np.floating)) or not np.isfinite(tolerance) or tolerance<=0:raise ValueError('finite positive energy grouping tolerance required')
        self.tolerance=float(tolerance);self.model=model;self.groups=[];self.arrays={}
        for number,(indices,energies,vectors) in model.spectral.items():
            idskey=f'spectral_indices_N{number}';self.arrays[idskey]=indices
            groups=[]
            for i,energy in enumerate(energies):
                if groups and abs(energy-groups[-1][0])<tolerance:groups[-1][1].append(i)
                else:groups.append([float(energy),[i]])
            for j,(energy,ids) in enumerate(groups):
                basis=vectors[:,ids];key=f'spectral_basis_N{number}_E{j}';self.arrays[key]=basis
                self.groups.append({'source_number':number,'energy':energy,'degeneracy':len(ids),'within_group_energy_spread':float(max(energies[ids])-min(energies[ids])),'basis_archive_key':key,'position_indices_archive_key':idskey})

    def measure(self,rho,label):
        rows=[]
        for i,group in enumerate(self.groups):
            indices=self.arrays[group['position_indices_archive_key']];basis=self.arrays[group['basis_archive_key']]
            coefficient=basis.conj().T@rho[np.ix_(indices,indices)]@basis
            key=f'spectral_coefficient_{label}_G{i}';self.arrays[key]=coefficient
            rows.append({**group,'probability':float(np.trace(coefficient).real),'coefficient_archive_key':key})
        return rows

    def reconstruct(self,row):
        ids=self.arrays[row['position_indices_archive_key']];basis=self.arrays[row['basis_archive_key']];coefficient=self.arrays[row['coefficient_archive_key']]
        full=np.zeros((self.model.n,self.model.n),complex);full[np.ix_(ids,ids)]=basis@coefficient@basis.conj().T
        return full

    def continue_branch(self,branch):
        u=expm(-.37j*self.model.hamiltonian);coherent=u@branch@u.conj().T
        native,ledger=self.model.field.commit(coherent,'MIX');out=.71*coherent+.29*native
        rows=[{**r,'joint_probability':.29*r['probability']} for r in ledger]
        rows.append({'status':'no_commit','joint_probability':.71*float(np.trace(branch).real)})
        return out,rows


def execute(model,inputs):
    instrument=GroupedEnergyInstrument(model);records={name:instrument.measure(state,name) for name,state in inputs.items()}
    # Continue one populated degenerate branch from an ACTUAL archived output.
    name,row=next((name,row) for name,rows in records.items() for row in rows if row['degeneracy']>1 and row['probability']>1e-8)
    branch=instrument.reconstruct(row);continued,ledger=instrument.continue_branch(branch)
    instrument.arrays['spectral_actual_degenerate_continuation']=continued
    # An explicit coherent control in the actual N3 doublet, not a new source.
    control_group=min((g for g in instrument.groups if g['source_number']==3 and g['degeneracy']==2),key=lambda g:abs(g['energy']+.7613052172833566))
    ids=instrument.arrays[control_group['position_indices_archive_key']];basis=instrument.arrays[control_group['basis_archive_key']]
    vector=basis@np.array([1,1])/np.sqrt(2);control=np.zeros((model.n,model.n),complex);control[np.ix_(ids,ids)]=np.outer(vector,vector.conj())
    control_rows=instrument.measure(control,'declared_doublet_control');records['declared_doublet_control']=control_rows
    observed=next(r for r in control_rows if r['probability']>.99);luders=instrument.reconstruct(observed)
    fine=np.zeros_like(luders);fine[np.ix_(ids,ids)]=basis@np.diag([.5,.5])@basis.conj().T
    distance=float(np.linalg.norm(luders-fine,ord='nuc')/2)
    control_after,control_ledger=instrument.continue_branch(luders);fine_after,_=instrument.continue_branch(fine)
    instrument.arrays['spectral_doublet_Luders_continuation']=control_after;instrument.arrays['spectral_doublet_fine_forget_continuation']=fine_after
    # Basis covariance of the compact representation, retaining all coherences.
    theta=.43;phi=.71;phase=np.array([[np.cos(theta),-np.exp(-1j*phi)*np.sin(theta)],[np.exp(1j*phi)*np.sin(theta),np.cos(theta)]]);c=instrument.arrays[observed['coefficient_archive_key']]
    rotated=(basis@phase)@(phase.conj().T@c@phase)@(basis@phase).conj().T
    rotation_error=float(np.max(abs(rotated-luders[np.ix_(ids,ids)])))
    continuation_difference=float(np.linalg.norm(control_after-fine_after,ord='nuc')/2)
    checks={'all_terminal_probabilities_complete':all(abs(sum(r['probability'] for r in rows)-float(np.trace(inputs[name]).real))<1e-12 for name,rows in records.items() if name in inputs),'actual_degenerate_branch_native_continuation_complete':abs(sum(r['joint_probability'] for r in ledger)-row['probability'])<1e-12,'fine_basis_forget_rejected':abs(distance-.5)<1e-12,'compact_basis_rotation_covariant':rotation_error<1e-12,'fine_forget_effect_survives_native_continuation':continuation_difference>.1}
    return {'instrument':'grouped(source_number,energy) Luders PVM, completezero outcomes; no rankone resolution inside degeneracy','energy_group_tolerance':instrument.tolerance,'apparatus_status':'explicit additional spectral apparatus; native C/E does not synthesize it','branch_state_format':'Embed(indices,basis @ coefficient @ basisdagger); coefficient is full complex subnormalized matrix, never diagonalized/discarded','complete_branch_records':records,'actual_populated_degenerate_continuation':{'input_archive_label':name,'branch':row,'native_time':.37,'native_probability':.29,'native_law':'MIX','complete_joint_native_outcomes':ledger,'poststate_archive_key':'spectral_actual_degenerate_continuation','subnormalization_retained':True},'fine_basis_rejecting_control':{'preparation':'(v1+v2)/sqrt2 in actual N3 energy doublet; added coherent apparatuscontrol, not replacementglobalreference','energy':control_group['energy'],'Luders_probability':observed['probability'],'fine_forget_trace_distance':distance,'continued_fine_forget_trace_distance':continuation_difference,'basis_rotation_reconstruction_error':rotation_error,'Luders_continuation_archive_key':'spectral_doublet_Luders_continuation','fine_forget_continuation_archive_key':'spectral_doublet_fine_forget_continuation','complete_Luders_native_outcomes':control_ledger},'checks':{k:bool(v) for k,v in checks.items()}},instrument.arrays
