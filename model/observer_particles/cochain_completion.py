#!/usr/bin/env python3
"""Complete the actual Whitney cone through degree three, retaining all modes.

This is an alternative conditional geometry, not a modification of the
240-port radial-cylinder experiment or a selected physical Weyl regulator.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import hashlib,importlib.util,json,math,sys
from pathlib import Path
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
sys.dont_write_bytecode=True
import numpy as np
import sympy as s
from scipy.linalg import block_diag,eigh,expm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PUBLIC=ROOT

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def maximum(x):return float(np.max(abs(x)))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
def roots(a):
    values,basis=eigh((a+a.T)/2)
    assert min(values)>0
    return (basis*np.sqrt(values))@basis.T,(basis/np.sqrt(values))@basis.T,values

def top_incidence(mesh):
    lookup={tuple(sorted(face)):i for i,face in enumerate(mesh['faces'])}
    out=np.zeros((len(mesh['tetrahedra']),len(mesh['faces'])),dtype=np.int64)
    for row,tet in enumerate(mesh['tetrahedra']):
        for omitted in range(4):
            face=[x for i,x in enumerate(tet) if i!=omitted]
            col=lookup[tuple(sorted(face))];registered=mesh['faces'][col]
            permutation=[registered.index(x) for x in face]
            inversions=sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3))
            out[row,col]=(-1)**(omitted+inversions)
    return out

def main():
    upstream=HERE/'results/geometric_matter_join.json'
    upstream_arrays=HERE/'results/geometric_matter_join_arrays.npz'
    previous=json.loads(upstream.read_text())
    assert previous['producer_sha256']==sha(HERE/'geometric_matter_join.py')
    assert previous['artifacts'][upstream_arrays.name]==sha(upstream_arrays)
    public=load('cochain_existing_whitney',PUBLIC/'code/equation_of_state/maxwell/maxwell_eos.py')
    mesh,(_,_,base_local,_)=public.mesh_and_forms()
    d0=np.array(mesh['D'],dtype=np.int64);d1=np.array(mesh['C'],dtype=np.int64);d2=top_incidence(mesh)
    sizes=[13,42,50,20];offsets=np.cumsum([0,*sizes]);total=sum(sizes)
    assert d0.shape==(42,13) and d1.shape==(50,42) and d2.shape==(20,50)
    assert np.array_equal(d1@d0,np.zeros((50,13),dtype=np.int64))
    assert np.array_equal(d2@d1,np.zeros((20,42),dtype=np.int64))
    ranks=[s.Matrix(d).rank() for d in (d0,d1,d2)]
    assert ranks==[12,30,20]
    betti=[sizes[k]-(ranks[k-1] if k else 0)-(ranks[k] if k<3 else 0) for k in range(4)]
    assert betti==[1,0,0,0]
    with np.load(upstream_arrays) as a:
        metric=a['metric'];coframe=a['coframe'];lapse=math.sqrt(-metric[0,0]);h=metric[1:,1:]
        masses=[lapse*a['scalar_kinetic'],lapse*a['Maxwell_electric_kinetic'],a['Maxwell_magnetic_metric']/lapse]
        inherited_kinetic=a['scalar_kinetic'];inherited_stiffness=a['scalar_stiffness']
        assert np.array_equal(a['Whitney_D'],d0) and np.array_equal(a['Whitney_C'],d1)
    mm,nn,physical_local,volume=public.forms(mesh,coframe[1:,1:])
    assert maximum(mm-masses[1])<2e-11 and maximum(nn-masses[2])<2e-11
    # Whitney top forms have unit integral on their own tetrahedron, so
    # their squared norm is the reciprocal physical tetrahedral volume.
    rootdet=math.sqrt(np.linalg.det(h));m3=np.diag([1/(rootdet*row['volume']) for row in base_local])
    independent_m3=np.diag([1/row['volume'] for row in physical_local])
    assert maximum(m3-independent_m3)<2e-11
    independent_m0=np.zeros((13,13))
    for row in physical_local:
        ids=row['vertices'];independent_m0[np.ix_(ids,ids)]+=row['volume']*(np.ones((4,4))+np.eye(4))/20
    assert maximum(independent_m0-masses[0])<2e-11
    masses.append(m3);sqrtm=[];invsqrtm=[];mass_spectra=[]
    for m in masses:
        a,b,values=roots(m);sqrtm.append(a);invsqrtm.append(b);mass_spectra.append(values)
    differentials=[sqrtm[k+1]@d@invsqrtm[k] for k,d in enumerate((d0,d1,d2))]
    full_d=np.zeros((total,total))
    for k,d in enumerate(differentials):full_d[offsets[k+1]:offsets[k+2],offsets[k]:offsets[k+1]]=d
    nilpotence=maximum(full_d@full_d)
    hodge=lapse*(full_d+full_d.T);squared=hodge@hodge
    assert nilpotence<3e-11 and maximum(hodge-hodge.T)==0
    # Compare the scalar block with the independently retained scalar action,
    # not a second expression constructed from the same incidence identity.
    _,inverse_kinetic_root,_=roots(inherited_kinetic)
    scalar_target=inverse_kinetic_root@inherited_stiffness@inverse_kinetic_root
    scalar_error=maximum(squared[:13,:13]-scalar_target)
    assert scalar_error<3e-11
    off_diagonal=squared.copy()
    for k in range(4):off_diagonal[offsets[k]:offsets[k+1],offsets[k]:offsets[k+1]]=0
    assert maximum(off_diagonal)<3e-11
    energies=eigh(hodge,eigvals_only=True)
    assert sum(abs(energies)<1e-9)==1 and sum(energies>1e-9)==sum(energies<-1e-9)==62
    parity=np.diag(np.concatenate([np.full(n,(-1)**k) for k,n in enumerate(sizes)]))
    assert maximum(parity@hodge+hodge@parity)==0
    # A formerly zero coexact edge direction becomes dynamical when faces
    # are retained. No postselection or new potential lifts it.
    seed=np.zeros(50);seed[0]=1
    edge=np.linalg.solve(masses[1],d1.T@seed)
    edge_white=sqrtm[1]@edge;edge_white/=np.linalg.norm(edge_white)
    truncated=lapse*np.block([[np.zeros((13,13)),differentials[0].T],
                              [differentials[0],np.zeros((42,42))]])
    old=np.concatenate([np.zeros(13),edge_white])
    old_residual=float(np.linalg.norm(truncated@old))
    start=np.zeros(total,complex);start[13:55]=edge_white
    full_speed=float(np.linalg.norm(hodge@start));assert old_residual<3e-12 and full_speed>1e-3
    rows=[];laws=[]
    for ticks in (1,2,25):
        t=ticks/math.sqrt(630);state=expm(-1j*t*hodge)@start;law=abs(state)**2
        assert abs(law.sum()-1)<3e-12
        rows.append({'comparison_incident_count':ticks,'comparison_coordinate_time':t,
            'complete_number_probability_sum':float(law.sum()),
            'probabilities_by_form_degree':[float(law[offsets[k]:offsets[k+1]].sum()) for k in range(4)],
            'state_change_norm':float(np.linalg.norm(state-start))})
        laws.append(law)
    assert rows[-1]['probabilities_by_form_degree'][2]>1e-5
    # This is a genuine finite all-sector CAR continuation if CAR statistics
    # are declared: its covariance is faithful and its full one-body U is
    # unitary. No2^125 Fock enumeration is needed for the algebraic lift.
    occupation=1/(1+np.exp(energies));assert np.all(occupation>0) and np.all(occupation<1)
    artifact=HERE/'results/cochain_completion_arrays.npz'
    np.savez_compressed(artifact,d0=d0,d1=d1,d2=d2,M0=masses[0],M1=masses[1],M2=masses[2],M3=masses[3],
        whitened_d=full_d,coordinate_time_Hodge_Dirac=hodge,scalar_action_target=scalar_target,
        full_cochain_energies=energies,CAR_energy_occupations=occupation,
        former_edge_zero_mode=start,complete_cochain_number_laws=np.array(laws),form_degree_offsets=offsets)
    paths=[upstream,upstream_arrays,HERE/'geometric_matter_join.py',
        PUBLIC/'code/equation_of_state/maxwell/maxwell_eos.py',
        PUBLIC/'code/electromagnetism/cone_whitney_bridge.py',
        PUBLIC/'docs/research/FERMIONIC_SOURCE_CLOCKS.md',
        PUBLIC/'model/relative_observer_geometry/README.md']
    for module in tuple(sys.modules.values()):
        p=getattr(module,'__file__',None)
        if p and Path(p).is_file() and Path(p).is_relative_to(PUBLIC):paths.append(Path(p))
    result={'status':'PASS_FULL_WHITNEY_COCHAIN_COMPLETION_WITH_SHARED_CONE_SCALAR_BLOCK',
        'producer_sha256':sha(__file__),'source_sha256':bound_pins({str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}),
        'artifacts':{artifact.name:sha(artifact)},'empirical_inputs':[],
        'existing_domain':{'name':'existing native-boundary tetrahedral cone','cochain_dimensions':sizes,
            'total_one_particle_dimension':total,'metric_and_lapse_consumed_from':'geometric_matter_join',
            'boundary':'full finite cochain complex and its metric adjoints; inherited natural variational/absolute cohomology domain',
            'same_as_240_port_radial_cylinder':False,'new_topology_invented':False},
        'exact_complex':{'d1_d0_zero':True,'d2_d1_zero':True,'differential_ranks':ranks,'Betti_numbers':betti,
            'truncated_node_edge_zero_modes':31,'complete_all_degree_zero_modes':1,
            'zero_modes_by_degree':betti,'Euler_characteristic':sum((-1)**k*n for k,n in enumerate(sizes))},
        'same_metric_action':{'Hodge_mass_eigenvalue_ranges':[[float(min(v)),float(max(v))] for v in mass_spectra],
            'top_form_mass_formula':'M3_tt=1/(sqrt(det h)*coordinate_volume_t)',
            'top_form_independent_pullback_error':maximum(m3-independent_m3),
            'coordinate_generator':'D=N*(d_hat+d_hat^dagger), d_hat_k=M_(k+1)^(1/2)*d_k*M_k^(-1/2)',
            'nilpotence_residual_after_whitening':nilpotence,'square_degree_offdiagonal_error':maximum(off_diagonal),
            'inherited_scalar_action_block_error':scalar_error,
            'scalar_comparison':'D^2|degree0 equals M_scalar^(-1/2) K_scalar M_scalar^(-1/2) from the independent retained action'},
        'full_spectrum':{'negative':62,'zero':1,'positive':62,
            'positive_energy_range':[float(min(energies[energies>1e-9])),float(max(energies))],
            'form_parity_even_dimension':63,'form_parity_odd_dimension':62,
            'form_parity_is_physical_Weyl_chirality':False,
            'finite_CAR_reference_at_beta1_faithful':True,
            'zero_mode_CAR_entropy':'log(2), not an extensive edge-cycle entropy',
            'scalar_constant_mode_thermal_caveat':'The natural scalar constant mode still needs an explicit zero-mode treatment for a trace-class scalar Gibbs reference.'},
        'lifted_edge_mode':{'truncated_generator_action_norm':old_residual,
            'full_generator_action_norm':full_speed,'complete_125_outcome_readouts':rows,
            'new_mass_potential_or_fitted_coefficient_inserted':False,
            'new_native_125_mode_pulse_history_executed':False},
        'conclusion':{'edge_cycle_obstruction_of_truncated_complex_removed_in_this_alternative':True,
            'same_scalar_block_preserved_on_this_cone':True,
            'complete_CAR_extension_available_by_existing_source_theorem':True,
            'minimal_physical_Weyl_field_or_SM_tastes_selected':False,
            'entire_existing_cylinder_experiment_replaced_or_reinterpreted':False,
            'this_cochain_completion_selected_by_native_repair':False,
            'continuum_chiral_regulator_or_full_physical_fermion_net_proved':False}}
    (HERE/'results/cochain_completion.json').write_text(json.dumps(public_receipt(result),indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status':result['status'],'ranks':ranks,'Betti':betti,
        'scalar_block_error':scalar_error,'former_edge_mode_full_speed':full_speed},indent=2))

if __name__=='__main__':main()
