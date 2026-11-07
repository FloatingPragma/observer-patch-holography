"""Complete finite-carrier addition data and Higgs self-energy audit.

No measured mass or dimensional energy unit is read. The native-prepared
family matrices are evaluated from the declared exact source expressions.
The state is the ground state of an explicitly declared finite Hamiltonian;
neither that selection nor the scalar controls are claimed source-derived.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import splu
from .interacting_carrier import HERE,build_join,primary_family_matrices


def maxabs(a):
    if sparse.issparse(a):return float(np.max(np.abs(a.data),initial=0))
    return float(np.max(np.abs(a),initial=0))


def frob(a):
    if sparse.issparse(a):return float(np.linalg.norm(a.data))
    return float(np.linalg.norm(a))


def block_eigh(model):
    """All spin-weight blocks, with every eigenvector retained."""
    h=model['H'];n=h.shape[0]
    if h.shape!=(n,n) or not np.all(np.isfinite(h.data)) or maxabs(h-h.getH())>1e-12:
        raise ValueError('a finite Hermitian Hamiltonian is required')
    weights=model['observables']['twice_Sz'].diagonal().real.astype(int)
    if maxabs(h@model['observables']['twice_Sz']-model['observables']['twice_Sz']@h)>1e-13:
        raise ArithmeticError('spin weight is not an exact numerical block')
    values=[];vectors=[];block_labels=[];blocks=[]
    for weight in np.unique(weights):
        indices=np.flatnonzero(weights==weight)
        block=h[indices][:,indices].toarray()
        e,u=eigh(block,driver='evd')
        full=np.zeros((n,len(e)),complex);full[indices,:]=u
        values.extend(e);vectors.append(full);block_labels.extend([int(weight)]*len(e))
        blocks.append({'twice_Sz':int(weight),'dimension':len(indices),
                      'lowest':float(e[0]),'highest':float(e[-1])})
    e=np.asarray(values);u=np.column_stack(vectors);label=np.asarray(block_labels)
    order=np.argsort(e,kind='stable');e=e[order];u=u[:,order];label=label[order]
    hu=h@u;res=maxabs(hu-u*e)
    orth=maxabs(u.conj().T@u-np.eye(n))
    # Numerical degeneracy is a reporting diagnostic; the exact spin labels
    # remain available, and no eigenstates or residues are combined/dropped.
    scale=max(1.,maxabs(h));tol=64*np.finfo(float).eps*n*scale
    ground_count=int(np.sum(e-e[0]<=tol))
    receipt={'all_eigenvectors_retained':True,'dimension':n,'blocks':blocks,
             'eigen_residual_max':res,'orthogonality_max':orth,
             'hermiticity_max':maxabs(h-h.getH()),
             'ground_numerical_multiplicity':ground_count,
             'ground_grouping_tolerance':tol,'ground_energy':float(e[0]),
             'first_gap_outside_reported_ground_group':float(e[ground_count]-e[0]) if ground_count<n else None}
    if res>2e-11 or orth>2e-11:raise ArithmeticError('complete eigensystem residual')
    return e,u,label,receipt


def operator_checks(model):
    h=model['H'];o=model['observables']
    comm=lambda a,b:a@b-b@a
    r={'site_continuity_max':maxabs(1j*comm(h,o['N_at1'])-o['J1']+o['J0']),
       'Higgs_occupancy_conversion_max':maxabs(1j*comm(h,o['N_H'])-o['JH']),
       'fermion_pair_conversion_max':maxabs(1j*comm(h,o['N_F'])-2*o['JH']),
       'Qstar_conservation_max':maxabs(comm(h,o['N_F']-2*o['N_H'])),
       'lepton_Higgs_number_identity_max':maxabs(o['N_L_total']-o['N_H']),
       'bulk_charge_identity_max':maxabs(6*o['N_e_total']-3*o['N_L_total']-3*o['N_H']-o['matter_scalar_charge6']),
       'complete_reference_charge_neutrality_max':maxabs(o['matter_scalar_charge6']+o['reference_charge6']),
       'spin_weight_conservation_max':maxabs(comm(h,o['twice_Sz'])),
       'Higgs_hopping_commutator_frobenius':frob(comm(model['H_Higgs'],model['H_hop'])),
       'Higgs_field_commutator_frobenius':frob(comm(model['H_Higgs'],model['H_field'])),
       'Higgs_term_frobenius':frob(model['H_Higgs']),
       'hopping_term_frobenius':frob(model['H_hop']),
       'field_term_frobenius':frob(model['H_field'])}
    if max(v for k,v in r.items() if k.endswith('_max'))>1e-12:raise ArithmeticError('finite operator identities')
    if not r['Higgs_hopping_commutator_frobenius']>0:raise ArithmeticError('missing dynamical interaction join')
    return r


def addition_data(joined,eigensystems,prefix,arrays):
    e0,u0,_,r0=eigensystems[0];e1,u1,_,r1=eigensystems[1]
    ng=r0['ground_numerical_multiplicity']
    # All columns of the normalized ground projector are used if degeneracy
    # is present. No arbitrary pure state is hidden in that case.
    ground=u0[:,:ng]
    labels=list(sorted(joined['addition']))
    created=np.stack([joined['addition'][k]@ground for k in labels],axis=1)
    amplitudes=np.einsum('ni,nag->iag',u1.conj(),created,optimize=True)
    # A numerical low-energy cluster can contain distinct energies. Retain
    # each preparation component rather than replacing them by one energy.
    residues=np.einsum('nag,nbg->ngab',amplitudes.conj(),amplitudes,optimize=True)/ng
    residues=residues.reshape(len(e1)*ng,len(labels),len(labels))
    weights=np.einsum('naa->na',residues).real
    gram=np.einsum('nag,nbg->ab',created.conj(),created,optimize=True)/ng
    closure=maxabs(np.sum(residues,axis=0)-gram)
    energies=(e1[:,None]-e0[None,:ng]).reshape(-1)
    # Each row is retained even when exactly dark. A tolerance is used only
    # to identify lines worth printing in a short human-readable receipt.
    per_level=weights.sum(axis=1);visible=np.flatnonzero(per_level>1e-15)
    arrays.update({prefix+'_addition_frequencies':energies,
                   prefix+'_addition_amplitudes':amplitudes,
                   prefix+'_addition_residues':residues,
                   prefix+'_addition_weights':weights,
                   prefix+'_addition_created_vectors':created,
                   prefix+'_addition_probe_gram':gram,
                   prefix+'_addition_final_state_indices':np.repeat(np.arange(len(e1)),ng),
                   prefix+'_addition_initial_state_indices':np.tile(np.arange(ng),len(e1)),
                   prefix+'_addition_initial_energies':e0[:ng],
                   prefix+'_addition_probe_labels':np.asarray(labels,dtype=int)})
    if closure>1e-11:raise ArithmeticError('complete addition spectral sum rule')
    return {'ground_state_preparation':'normalized numerical low-energy projector with each component energy retained',
            'ground_rank':ng,'all_Q1_eigenstate_rows_retained':len(e1),
            'all_initial_final_pairs_retained':len(energies),
            'spectral_sum_rule_max':closure,
            'total_addition_weight':float(np.trace(gram).real),
            'individual_probe_total_weights':np.diag(gram).real.tolist(),
            'lowest_visible_rows':[{'state':int(i//ng),'initial_state':int(i%ng),'frequency':float(energies[i]),
                                   'sum_six_probe_weights':float(per_level[i]),
                                   'per_probe_weights':weights[i].tolist()} for i in visible[:12]],
            'negative_frequency_visible_rows':int(np.sum((energies<0)&(per_level>1e-15))),
            'print_visibility_threshold':1e-15,
            'threshold_does_not_truncate_arrays':True}


def schur_checks(model,arrays,prefix):
    h=model['H'];nh=model['observables']['N_H'].diagonal().real
    p=np.flatnonzero(nh==0);q=np.flatnonzero(nh==1)
    a=h[p][:,p].toarray();d=h[q][:,q].tocsc();c=h[q][:,p].toarray()
    b=c.conj().T
    coupling_gram=b@c
    # All off-diagonal occupancy blocks must literally be the Higgs term.
    mismatch=maxabs(h[q][:,p]-model['H_Higgs'][q][:,p])
    if mismatch>1e-12:
        raise ArithmeticError('offdiagonal block is not the declared Higgs interaction')
    receipts=[]
    for j,z in enumerate((1.+.5j,-4.+.5j)):
        invq=splu(z*sparse.eye(len(q),format='csc')-d)
        sigma=b@invq.solve(c)
        schur=z*np.eye(len(p))-a-sigma
        reduced=np.linalg.inv(schur)
        rhs=np.zeros((len(nh),len(p)),complex);rhs[p,:]=np.eye(len(p))
        full=splu(z*sparse.eye(len(nh),format='csc')-h.tocsc()).solve(rhs)[p,:]
        imaginary=(sigma-sigma.conj().T)/(2j)
        causality=eigh(-imaginary,eigvals_only=True)
        arrays[prefix+f'_self_energy_{j}']=sigma
        receipts.append({'z':[z.real,z.imag], 'P_dimension':len(p),'Q_dimension':len(q),
            'self_energy_frobenius':frob(sigma),'self_energy_operator_norm':float(np.linalg.norm(sigma,2)),
            'Schur_full_resolvent_max':maxabs(reduced-full),
            'minus_imaginary_self_energy_min_eigenvalue':float(causality[0]),
            'minus_imaginary_self_energy_max_eigenvalue':float(causality[-1])})
        if maxabs(reduced-full)>1e-10 or causality[0]<-1e-11:raise ArithmeticError('Schur/casuality check')
    arrays[prefix+'_self_energy_P_indices']=p;arrays[prefix+'_self_energy_Q_indices']=q
    return {'offdiagonal_is_literal_Higgs_max':mismatch,
            'coupling_gram_min_eigenvalue':float(eigh(coupling_gram,eigvals_only=True)[0]),
            'coupling_gram_max_eigenvalue':float(eigh(coupling_gram,eigvals_only=True)[-1]),
            'coupling_frobenius':frob(c),'resolvent_checks':receipts,
            'Q_resolvent_contains_native_electric_color_loop_and_matter_hopping':True}


def finite_difference_audit(joined_builder,joined):
    # A separate nonzero phase check guards current derivative signs.
    step=1e-5;records={}
    shifted=joined_builder(phase0=.17,phase1=-.23,higgs_phase=.11)
    for phase,key in (('phase0','J0'),('phase1','J1'),('higgs_phase','JH')):
        a={'phase0':.17,'phase1':-.23,'higgs_phase':.11};b=dict(a)
        a[phase]+=step;b[phase]-=step
        plus=joined_builder(**a);minus=joined_builder(**b)
        for q in (0,1):
            derivative=(plus['sectors'][q]['H']-minus['sectors'][q]['H'])/(2*step)
            err=maxabs(derivative-shifted['sectors'][q]['observables'][key])
            records[f'q{q}_{key}']=err
            if err>2e-10:raise ArithmeticError('phase current derivative')
    return {'step':step,'maximum_entry_residuals':records}


def run(root=HERE,prior=None,output=None):
    start=time.monotonic();root=Path(root)
    output=Path(output) if output is not None else root
    output.mkdir(parents=True,exist_ok=True)
    if prior is None:
        prior=HERE.parent/'dependencies'/'exact_family_operators.json'
        y,sigma=primary_family_matrices(prior)
    else:
        prior=Path(prior)
        if prior.suffix=='.json':
            y,sigma=primary_family_matrices(prior)
        else:
            with np.load(prior,allow_pickle=False) as source:
                y=source['source_Y'].copy();sigma=source['source_sigma'].copy()
    parameters=dict(family0=y,family1=sigma,yukawa=y,electric=1.,magnetic=1.,hopping=1.,
                    fermion_energy=0.,scalar_energy=1.,reference_energy=0.)
    builder=lambda **kw:build_join(root,**(parameters|kw))
    dressed=builder(higgs_coupling=1.)
    arrays={'source_Y':y,'source_sigma':sigma}
    results={};all_eigensystems={}
    for name,g in (('dressed',1.),('no_higgs',0.)):
        joined=builder(higgs_coupling=g);eigensystems={};sector_receipts={}
        for q in (0,1):
            model=joined['sectors'][q]
            eigen=block_eigh(model);eigensystems[q]=eigen
            e,u,sz,receipt=eigen
            arrays.update({f'{name}_q{q}_energies':e,f'{name}_q{q}_vectors':u,
                           f'{name}_q{q}_twice_Sz':sz})
            if name=='dressed':receipt['operator_identities']=operator_checks(model)
            sector_receipts[str(q)]=receipt
        all_eigensystems[name]=eigensystems
        results[name]={'higgs_coupling':g,'sectors':sector_receipts,
                      'charged_addition':addition_data(joined,eigensystems,name,arrays)}
        print(json.dumps({'stage':'complete_spectrum','case':name,
                          'ground_Q0':eigensystems[0][0][0],
                          'ground_Q1':eigensystems[1][0][0],
                          'seconds':time.monotonic()-start}),flush=True)
    results['dressed']['Q1_Higgs_self_energy']=schur_checks(dressed['sectors'][1],arrays,'dressed_q1')
    results['dressed']['Q0_Higgs_self_energy']=schur_checks(dressed['sectors'][0],arrays,'dressed_q0')
    # The complete addition correlation can be compared without assigning
    # a particle interpretation to an individual eigenvalue or dark line.
    correlation_changes=[]
    for t in (0.,.5,1.,2.):
        cs={}
        for name in ('dressed','no_higgs'):
            cs[name]=np.einsum('n,nab->ab',np.exp(-1j*t*arrays[name+'_addition_frequencies']),
                               arrays[name+'_addition_residues'])
        arrays[f'addition_correlation_dressed_t{str(t).replace(".","p")}']=cs['dressed']
        correlation_changes.append({'time':t,'complete_correlation_change_frobenius':frob(cs['dressed']-cs['no_higgs'])})
    currents=finite_difference_audit(builder,dressed)
    mu=.375  # exact binary rational, so this calibration check has no roundoff ambiguity
    common_shift=builder(fermion_energy=mu,scalar_energy=1.-2*mu,reference_energy=-mu)
    fixed_reference=builder(fermion_energy=mu,scalar_energy=1.-2*mu,reference_energy=0.)
    calibration={}
    for q in (0,1):
        h=dressed['sectors'][q]['H'];unit=sparse.eye(h.shape[0],format='csr')
        calibration[str(q)]={
            'consistent_common_gauge_shift_H_difference_max':maxabs(common_shift['sectors'][q]['H']-h),
            'fixed_reference_shift_minus_mu_Qstar_max':maxabs(fixed_reference['sectors'][q]['H']-h-mu*q*unit)}
        if max(calibration[str(q)].values())>1e-13:raise ArithmeticError('reference calibration covariance')
    np.savez_compressed(output/'joined_response_arrays.npz',**arrays)
    receipt={'schema':'oph.interacting.lepton_higgs_response.v1',
             'source':dressed['source'],'scope':dressed['scope'],
             'reference_grading':dressed['reference_grading'],
             'controls':dressed['parameters'],
             'family_assignment':{'edge0':'source_Y','edge1':'source_sigma','charged_lepton_vertex':'source_Y'},
             'input_family_packet':str(prior.resolve()),
             'input_family_packet_sha256':hashlib.sha256(prior.read_bytes()).hexdigest(),
             'source_selected_controls_or_ground_state_claimed':False,
             'empirical_particle_data_inputs':[],
             'family_commutator_frobenius':frob(y@sigma-sigma@y),
             'cases':results,'phase_derivative_audit':currents,
             'complete_addition_correlation_changes':correlation_changes,
             'energy_offset_identity':{
                 'exact':'epsilon_F*N_F+epsilon_H*N_H = epsilon_F*Qstar+(2*epsilon_F+epsilon_H)*N_H',
                 'exact_bulk_charge_identity':'6*N_e-3*N_L-3*N_H=6*Qstar; N_L=N_H; N_e=Qstar+N_H',
                 'reference_charge_identity':'q6_reference=-6*Qstar',
                 'addition_frequency_shift_under_mu_Qstar_with_reference_fixed':'mu',
                 'addition_frequency_shift_under_reference_energy':'reference_energy',
                 'within_each_Qstar_sector_density_and_record_dynamics_unchanged_by_constant_shift':True,
                 'consistent_common_gauge_shift':'H_bulk -> H_bulk+mu*Qstar; E_reference -> E_reference-mu; total physical H unchanged',
                 'common_gauge_potential_cancels_in_complete_reference_neutral_detector':True,
                 'direct_calibration_covariance_check_mu':mu,
                 'direct_calibration_covariance_checks':calibration,
                 'relative_reference_calibration_not_source_selected':True,
                 'intrinsic_particle_mass_arbitrariness_or_new_physical_chemical_potential_proved':False},
             'response_arrays_sha256':hashlib.sha256((output/'joined_response_arrays.npz').read_bytes()).hexdigest(),
             'elapsed_seconds':time.monotonic()-start}
    (output/'joined_response_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'stage':'complete_response_audit','seconds':time.monotonic()-start,
                      'output':str(output/'joined_response_receipt.json')}),flush=True)
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=HERE)
    parser.add_argument('--prior',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.root,args.prior,args.output)
