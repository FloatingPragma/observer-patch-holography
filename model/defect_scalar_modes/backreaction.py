"""All-mode quantum scalar response of the native unit interaction.

The reference fills every negative spatial mode with all fifteen native
internal labels. H_Y applied once gives a complete gauge-singlet scalar
response and a variational bound for the full joint Hamiltonian.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import linalg as la, sparse as sp

from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import tensor_checks
from model.defect_mode_overlap.overlap import profile_response
from model.compact_defect_pair.certification import complete_certificate
from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from .modes import build_scalar_source

ROOT=Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def native_coefficients():
    """Actual native tensors, orthogonality and global sea charge."""
    import sympy as s
    _,_,currents,vertices=representation()
    if any(s.trace(g)!=0 for g in currents):
        raise ArithmeticError('uniform occupation of all native labels is not a singlet')
    norms=[]
    for component in range(2):
        tensors=[vertices[name]['components'][component] for name in ('up','down','charged_lepton')]
        gram=s.Matrix(3,3,lambda i,j:s.trace(tensors[i].H*tensors[j]))
        if gram != s.diag(6,6,2):
            raise ArithmeticError('native unit tensor Gram changed')
        norms.append([[str(x) for x in row] for row in gram.tolist()])
    ward=tensor_checks()
    parts,coefficients=car.quantum_yukawa()
    literal=car.clean({key:value.subs({y:1 for y in coefficients.values()})
                      for key,value in car.add(*parts.values()).items()})
    branch_l1=[sum(abs(value) for (_,row,col),value in literal.items()
                   if row==branch and col==0) for branch in range(1,5)]
    if branch_l1 != [14]*4 or literal != car.adjoint(literal):
        raise ArithmeticError('the literal local scalar star changed')
    # The vacuum-to-occupied column has four blocks; each is bounded by 14.
    # The full self-adjoint scalar star has the same norm as that column.
    local_norm_bound=28
    return {'native_internal_labels':15,'all_twelve_generator_traces_zero':True,
            'all_labels_occupied_identically_in_each_negative_spatial_mode':True,
            'global_native_gauge_singlet':True,'local_Gauss_with_quantum_links_claimed':False,
            'component_tensor_Grams':norms,'complete_Ward_checks':ward,
            'pair_kernel_total_factor':56,'primitive_interaction_coefficient':'1',
            'literal_local_matrix_unit_words':len(literal),
            'vacuum_to_scalar_branch_coefficient_absolute_sums':list(map(int,branch_l1)),
            'local_interaction_norm_upper_exact':local_norm_bound}


def pair_kernel(h,energies,vectors,points):
    occupied=energies<0
    g=vectors[:,occupied]@vectors[:,occupied].conj().T
    hg=(vectors[:,occupied]*energies[occupied])@vectors[:,occupied].conj().T
    n=len(points);blocks=g.reshape(n,2,n,2).transpose(0,2,1,3)
    signs=np.array([(-1.)**sum(p) for p in points])
    determinant=np.linalg.det(blocks)*signs[:,None]*signs[None,:]
    determinant=(determinant+determinant.conj().T)/2
    moment=0.
    for x in range(n):
        gx=g[2*x:2*x+2,2*x:2*x+2]
        hx=hg[2*x:2*x+2,2*x:2*x+2]
        adjugate=np.trace(gx)*np.eye(2)-gx
        moment-=float(np.trace(adjugate@hx).real)
        if abs(np.trace(gx)-1)>2e-10:
            raise ArithmeticError('source particle-hole local trace changed')
    if la.eigvalsh(determinant)[0]<-2e-10:
        raise ArithmeticError('complete pair Gram is not positive')
    return determinant,moment,g


def strict_variational_bound():
    """A rational margin, using the same parent's complete h gap certificate."""
    cert=complete_certificate()
    if (cert['dimension']!=686 or cert['mass']!=64 or cert['half_separation']!=1
            or cert['all_finite_source_eigenvalues_certified'] is not True
            or F(cert['certified_core_abs_energy_lower'])<=F(1,1000)):
        raise ArithmeticError('the matching complete h gap certificate is required')
    delta=F(1,1000);norm=F(66);a=delta/(norm+delta)
    coupling_squared=112*a*(1-a)
    energy_upper=F(197)  # Two fermion energies <=132; Omega<=sqrt(4109)<65.
    lowering=coupling_squared*energy_upper/(energy_upper**2+coupling_squared)
    pointer_norm=native_coefficients()['local_interaction_norm_upper_exact']*343
    return {'complete_spatial_gap_certificate_regenerated':True,
            'certified_parent_core_abs_energy_lower_exact':cert['certified_core_abs_energy_lower'],
            'spatial_gap_lower_exact':str(delta),'spatial_norm_upper_exact':str(norm),
            'core_occupation_eigenvalue_lower_exact':str(a),
            'interaction_image_norm_squared_lower_exact':str(coupling_squared),
            'first_image_H0_excitation_upper_exact':str(energy_upper),
            'trial_coefficient_of_HY_Omega_exact':'-1/197',
            'ground_energy_lowering_at_least_exact':str(lowering),
            'ground_energy_lowering_at_least_numeric':float(lowering),
            'complete_pointer_norm_upper':pointer_norm,
            'original_success_probability_lower_exact':str(coupling_squared/pointer_norm**2),
            'statement':'E_ground(H0+H_Y) <= E_sea - positive_margin in the global native singlet domain',
            'exact_interacting_ground_or_pole_reconstructed':False}


@lru_cache(maxsize=1)
def calculate_backreaction():
    native=native_coefficients();data=build_scalar_source();spatial=data['spatial']
    kernel,fermion_moment,g=pair_kernel(spatial['h'],spatial['energies'],spatial['vectors'],data['points'])
    factor=native['pair_kernel_total_factor']
    trace=float(np.trace(kernel).real)
    c2=factor*trace
    scalar_energy=factor*float(np.trace(data['Omega']@kernel).real)
    fermion_energy=factor*fermion_moment
    cost=(scalar_energy+fermion_energy)/c2
    ritz=-2*c2/(cost+np.sqrt(cost*cost+4*c2))
    # Both components of particle and antiparticle sectors are retained.
    particle=8*kernel+6*kernel.conj()
    antiparticle=particle.conj()
    modal_particle=np.diag(data['modes'].T@particle@data['modes']).real
    modal_anti=np.diag(data['modes'].T@antiparticle@data['modes']).real
    probabilities=2*(modal_particle+modal_anti)/c2
    if probabilities.min() < -2e-10 or abs(probabilities.sum()-1)>2e-10:
        raise ArithmeticError('complete scalar first-image record lost weight')
    if la.norm(kernel.imag)>2e-10:
        raise ArithmeticError('particle-hole pair kernel must be real')
    pair_values,pair_vectors=la.eigh(kernel.real)
    dominant=pair_vectors[:,-1]
    if dominant[data['core_sites'][-1]]<0:
        dominant=-dominant
    response=profile_response(dominant,value=spatial)
    bounds=strict_variational_bound();pointer=bounds['complete_pointer_norm_upper']
    success=c2/pointer**2
    return {'native':native,'kernel':kernel,'particle_Gram':particle,'antiparticle_Gram':antiparticle,
        'pair_Gram_spectrum':pair_values,'dominant_pair_profile':dominant,
        'reference':{'negative_spatial_modes':int(np.count_nonzero(spatial['energies']<0)),
            'native_fermion_number':15*int(np.count_nonzero(spatial['energies']<0)),
            'scalar_blank':True,'joint_interacting_ground':False,
            'source_global_native_singlet':True,
            'idle_conjugate_reference_retained':True,
            'combined_source_reference_fermion_parity':'even'},
        'interaction_image_norm_squared':c2,'scalar_energy_numerator':scalar_energy,
        'fermion_energy_numerator':fermion_energy,'normalized_first_image_H0_excitation':cost,
        'two_vector_Ritz_energy_relative_to_sea':float(ritz),
        'strict_variational_bound':bounds,
        'complete_scalar_mode_records':[{'mode':j,'frequency':float(w),
                'particle_weight_per_component':float(modal_particle[j]),
                'antiparticle_weight_per_component':float(modal_anti[j]),
                'conditional_first_image_probability':float(probabilities[j]),
                'original_pointer_and_mode_probability':float(success*probabilities[j])}
            for j,w in enumerate(data['frequencies'])],
        'complete_pointer':{'norm_bound':pointer,'success_probability':success,
            'failure_probability':1-success,'all_four_scalar_charge_component_branches_retained':True,
            'mode_projectors_sum_over_all_four_internal_scalar_branches':True,
            'success_Kraus':'H_Y/B',
            'failure_Kraus':'positive_sqrt(I-H_Y^2/B^2)',
            'complete_outcome_count':len(probabilities)+1,
            'probability_total':float(1-success+sum(success*probabilities))},
        'dominant_pair_mode':{'largest_Gram_eigenvalue':float(pair_values[-1]),
            'next_Gram_eigenvalue':float(pair_values[-2]),'core_profile_values':dominant[data['core_sites']].tolist(),
            'unit_ladder_strengths':response['singulars'].tolist(),
            'strong_to_weak_ratio':float(response['singulars'][0]/response['singulars'][-1]),
            'mode_is_dressed_Higgs_pole':False},
        'scalar_vacuum_normalization_transferred_to_interacting_ground':False}


def small_car_control():
    """Independent full eight-mode CAR first image, all five scalar states."""
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.]);j=np.array([[0,1],[-1,0.]])
    h=np.block([[2*sz,.5j*sx],[-.5j*sx,2*sx]])
    c=la.block_diag(j,-j)
    if la.norm(h@c+c@h.conj())>1e-12:
        raise ArithmeticError('small control must have the same particle-hole law')
    energies,vectors=la.eigh(h)
    points=((0,0,0),(1,0,0))
    kernel,moment,_=pair_kernel(h,energies,vectors,points)
    one_body=np.kron(h,np.eye(2));ef,uf=la.eigh(one_body)
    occ=uf[:,ef<0];modes=8;dimension=1<<modes
    annihilators=[]
    for index in range(modes):
        rows=[];cols=[];entries=[]
        for mask in range(dimension):
            if mask & (1<<index):
                rows.append(mask^(1<<index));cols.append(mask)
                entries.append((-1.)**((mask&((1<<index)-1)).bit_count()))
        annihilators.append(sp.csr_matrix((entries,(rows,cols)),shape=(dimension,dimension)))
    hf=sp.csr_matrix((dimension,dimension),dtype=complex)
    for a,b in zip(*np.nonzero(one_body)):
        hf+=one_body[a,b]*(annihilators[a].conj().T@annihilators[b])
    sea=np.zeros(dimension,complex)
    for labels in combinations(range(modes),occ.shape[1]):
        sea[sum(1<<x for x in labels)]=la.det(occ[list(labels)])
    ground=float(ef[ef<0].sum())
    if la.norm(hf@sea-ground*sea)>1e-10 or abs(np.vdot(sea,sea)-1)>1e-10:
        raise ArithmeticError('complete small Fock sea failed')
    states=np.zeros((5,dimension),complex)
    for x in range(2):
        base=4*x
        pair=(-1.)**x*(annihilators[base]@annihilators[base+3]+annihilators[base+1]@annihilators[base+2])
        states[1+x]=pair@sea
        states[3+x]=pair.conj().T@sea
    direct=states@states.conj().T
    error=max(la.norm(direct[1:3,1:3]-2*kernel),la.norm(direct[3:5,3:5]-2*kernel.conj()))
    c2=float(np.vdot(states,states).real)
    fermi=float(sum(np.vdot(row,hf@row-ground*row).real for row in states))
    if error>1e-10 or abs(c2-4*np.trace(kernel).real)>1e-10 or abs(fermi-4*moment)>1e-10:
        raise ArithmeticError('literal CAR image disagrees with complete Wick kernel')
    return {'full_Fock_dimension':dimension,'scalar_states':5,'joint_dimension':dimension*5,
            'complete_particle_antiparticle_Gram_error':float(error),
            'interaction_image_norm_squared':c2,'fermion_excitation_numerator':fermi,
            'every_state_and_both_scalar_branches_retained':True}
