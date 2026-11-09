"""Boundary interaction bound with the actual polynomial and Gauss sectors."""
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np

from model.particle_common_contract import contract, provenance
from model.particle_primary_boundary_information import model as information
from model.particle_primary_geometric_cut import model as cut
from model.particle_primary_gravity_premises import model as modular

SOURCE_INPUTS=tuple('model/particle_primary_boundary_correlation/'+name for name in
    ('__init__.py','__main__.py','model.py','test_model.py','PROOF.md','README.md'))+(
    'model/particle_primary_boundary_information/PROOF.md',
    'model/particle_primary_gravity_premises/PROOF.md',
    'model/particle_primary_geometric_cut/PROOF.md',
    'model/particle_primary_capacity/PROOF.md')


def interaction_cut(side,cells):
    geometry=cut.cut(side,cells)
    inside=set(geometry['interior_links']);candidates=set()
    # Every mixed plaquette contains an inside link; enumerate only its stars.
    for edge in inside:
        x,y=edge;axis=next(i for i in range(3) if x[i]!=y[i])
        for other in range(3):
            if axis==other:continue
            for shift in (-1,0):
                base=tuple(t+shift*(i==other) for i,t in enumerate(x))
                vertices=tuple(tuple(t+u*(k==axis)+v*(k==other)
                    for k,t in enumerate(base)) for u,v in product((0,1),repeat=2))
                if all(all(0<=t<side for t in z) for z in vertices):
                    candidates.add(tuple(sorted(vertices)))
    mixed=[]
    for vertices in sorted(candidates):
        edges={tuple(sorted((x,y))) for x,y in combinations(vertices,2)
               if sum(abs(a-b) for a,b in zip(x,y))==1}
        n=len(edges&inside)
        if 0<n<4:mixed.append(vertices)
    b=len(geometry['crossing_links'])
    if len(mixed)>2*b:raise ArithmeticError('mixed plaquette incidence bound failed')
    return {'induced_cut':geometry,'crossing_hopping_bonds':b,
            'mixed_plaquettes':tuple(mixed),'mixed_plaquette_count':len(mixed),
            'all_electric_terms_belong_to_one_regional_factor':True,
            'boundary_gauge_sectors_retained':True}


def decorrelate(rho,dimensions,multiplicities):
    """Comparison state, not a linear channel on unknown inputs."""
    W,nr,no,blocks=information.singlet_embedding(dimensions,multiplicities)
    if rho.shape!=(W.shape[1],W.shape[1]) or not np.allclose(np.trace(rho),1):
        raise ValueError('normalized physical Schur density required')
    if min(np.linalg.eigvalsh(rho))<=0:
        raise ValueError('faithful physical density required')
    sigma=np.zeros_like(rho,dtype=complex);within=0.
    for start,end,m,n,d in blocks:
        block=rho[start:end,start:end];p=float(np.trace(block).real)
        tau=block/p;tensor=tau.reshape(m,n,m,n)
        left=np.einsum('abcb->ac',tensor);right=np.einsum('abad->bd',tensor)
        sigma[start:end,start:end]=p*np.kron(left,right)
        within+=p*(information.entropy(left)+information.entropy(right)-information.entropy(tau))
    stats=information.information_control(rho,dimensions,multiplicities)
    gain=information.entropy(sigma)-information.entropy(rho)
    if abs(gain-within-stats['center_coherence'])>2e-11:
        raise ArithmeticError('physical decorrelation entropy identity failed')
    return sigma,{'entropy_gain':gain,'conditional_mutual_information':within,
                  'center_coherence':stats['center_coherence'],
                  'source_state_replaced_by_comparison':False}


def logarithm(matrix):
    values,vectors=np.linalg.eigh(matrix)
    if min(values)<=0:raise ValueError('strict positivity required')
    return (vectors*np.log(values))@vectors.conj().T


def polynomial_density(H,beta,order,scale):
    energy,vectors=np.linalg.eigh(H);U=float(max(energy));L=float(min(energy))
    y=beta*(U-energy)/(2*scale)
    if min(y)<-1e-12 or max(y)>.5+1e-12:
        raise ValueError('control polynomial outside source approximation domain')
    from math import factorial
    logw=2*scale*np.log(sum(y**k/factorial(k) for k in range(order+1)))
    logw-=max(logw);weights=np.exp(logw);weights/=sum(weights)
    return (vectors*weights)@vectors.conj().T,L,U


def control(H,H0,beta,order,scale,dimensions,multiplicities):
    rho,L,U=polynomial_density(H,beta,order,scale)
    sigma,stats=decorrelate(rho,dimensions,multiplicities)
    delta=sigma-rho;K=-logarithm(rho)
    energy=float(np.trace(delta@H).real)
    local_difference=float(np.trace(delta@H0).real)
    residual=K-beta*H;osc=float(np.ptp(np.linalg.eigvalsh(residual)))
    rel=float(np.trace(sigma@(logarithm(sigma)-logarithm(rho))).real)
    remainder=float(np.trace(delta@residual).real)
    gain=stats['entropy_gain'];boundary=float(np.linalg.norm(H-H0,2))
    if abs(local_difference)>2e-11 or abs(rel-(beta*energy+remainder-gain))>2e-11:
        raise ArithmeticError('source free-energy comparison failed')
    if gain>2*beta*boundary+osc+2e-11:
        raise ArithmeticError('boundary information upper bound failed')
    return dict(stats,energy_change=energy,inside_outside_energy_change=local_difference,
                relative_entropy=rel,modular_affine_remainder=remainder,
                actual_modular_oscillation=osc,boundary_operator_norm=boundary,
                comparison_identity_residual=rel-beta*energy-remainder+gain)


def noninteracting_control():
    # Actual positive degree-one polynomial form, on a separate two-qubit control.
    weights=(1296,625,625,256);den=sum(weights)
    determinant=F(weights[0]*weights[3]-weights[1]*weights[2],den**2)
    H=np.diag([0.,1.,1.,2.])
    result=control(H,H,1.,1,2,[1],[(2,2)])
    if determinant==0 or result['entropy_gain']<=0:
        raise ArithmeticError('the polynomial nonfactorization control vanished')
    return dict(result,domain='two-qubit polynomial control; not the primary carrier',
                exact_probability_numerators=weights,probability_denominator=den,
                exact_product_determinant=contract.serial(determinant),
                boundary_interaction_zero=True,exact_Gibbs_state_substituted=False)


def bind_member(member):
    contract.verify_member(member)
    inherited=information.bind_member(member);information.verify_bound(member,inherited)
    source=member['source'];a=F(source['spacing']);ell=source['per_link_word_degree']
    geometry=interaction_cut(source['side'],((0,0,0),(1,0,0)))
    b=geometry['crossing_hopping_bonds'];p=geometry['mixed_plaquette_count']
    bound=(30*b+F(source['magnetic_coefficient'])*p)/a
    uniform=(30+2*F(source['magnetic_coefficient']))*b/a
    slopes=modular.selected_modular_slope(member)
    beta=slopes['beta'];r=slopes['global_modular_affine_defect_range_upper']
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'cut_interactions':geometry,'actual_crossing_interaction_norm_upper':bound,
        'uniform_crossing_interaction_norm_upper':uniform,
        'selected_modular_slope':slopes,
        'conditional_information_plus_center_coherence_upper':2*beta*bound+r,
        'center_entropy_upper_coefficient_of_log_451':ell*b,
        'twice_representation_entropy_upper_coefficient_of_log_6':2*ell*b,
        'ambient_mutual_information_upper':{'log_451_coefficient':ell*b,
            'log_6_coefficient':2*ell*b,'rational_addend':2*beta*bound+r},
        'boundary_bound_per_lattice_area':{'log_451_coefficient':F(ell)/(a*a),
            'log_6_coefficient':F(2*ell)/(a*a),
            'rational_addend':2*beta*(30+2*F(source['magnetic_coefficient']))/a**3,
            'separate_polynomial_remainder':r},
        'comparison_entropy_gain':'sum p_lambda I(tau_lambda) + C_center',
        'comparison_state':'direct sum p_lambda tau_R,lambda tensor tau_out,lambda tensor singlet',
        'comparison_preserves_original_regional_marginals_and_center_weights':True,
        'comparison_preserves_complete_Gauss_support':True,
        'comparison_asserted_linear_channel_on_unknown_states':False,
        'comparison_replaces_original_source_preparation':False,
        'finite_boundary_information_bound_proved':True,
        'regional_modular_generator_identified_with_internal_energy':False,
        'positive_finite_continuum_area_coefficient_proved':False,
        'Einstein_dynamics_proved':False,'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_boundary_correlation.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('boundary correlation or source binding changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*information.SOURCE_INPUTS,*cut.SOURCE_INPUTS,*modular.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    value={'schema':'oph.primary_boundary_correlation_certificate.v1',
        'noninteracting_polynomial_control':noninteracting_control(),
        'members':[bind_member(contract.member(n)) for n in (0,1,2,8)],
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('boundary correlation source changed')
    return value
