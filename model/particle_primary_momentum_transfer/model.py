"""Source-bound finite Fourier filters and their joint spectral dictionary."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import isqrt

import sympy as s

from model.particle_common_contract import contract,provenance
from model.particle_primary_origin_momentum import model as origin
from model.particle_primary_refinement_response import model as response
from model.particle_filtered_probe import filter as apparatus
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS=(
 'model/particle_primary_momentum_transfer/__init__.py',
 'model/particle_primary_momentum_transfer/__main__.py',
 'model/particle_primary_momentum_transfer/model.py',
 'model/particle_primary_momentum_transfer/test_model.py',
 'model/particle_primary_momentum_transfer/PROOF.md',
 'model/particle_primary_momentum_transfer/README.md',
 'model/particle_primary_origin_momentum/PROOF.md',
 'model/particle_primary_spatial_algebra/PROOF.md',
 'model/particle_primary_field_algebra/PROOF.md',
 'model/particle_primary_refinement_response/PROOF.md',
 'model/particle_filtered_probe/PROOF.md',
 'model/particle_primary_momentum/PROOF.md',
 'model/particle_primary_mass_readout/PROOF.md',
 'Lean/QFT/PeriodicMomentum.lean',
)


def ceil_sqrt(value):
    value=F(value)
    if value<0:raise ValueError('nonnegative square required')
    k=isqrt(value.numerator//value.denominator)
    return k if k*k>=value else k+1


def _frequency(value):
    if not isinstance(value,(tuple,list)) or len(value)!=3 or any(type(k) is not int for k in value):
        raise ValueError('three signed integer spatial Fourier indices required')
    return tuple(value)


def kernel_leakage(spatial_aperture='1',time_aperture='1',energy_radius='2',momentum_radii=('2','2','2')):
    Q=ground._q(spatial_aperture,'spatial half aperture',True)
    S=ground._q(time_aperture,'time half aperture',True)
    dE=ground._q(energy_radius,'positive energy radius',True)
    if len(momentum_radii)!=3:raise ValueError('three positive momentum radii required')
    dp=tuple(ground._q(x,'positive momentum radius',True) for x in momentum_radii)
    return min(F(1),max(1/(S*dE)**2,*(1/(Q*x)**2 for x in dp)))


def original_rectangle_lower(success_lower,leakage,probability_error='0'):
    """Use a limiting/admitted spectral success bound, not an isolated finite histogram."""
    p=ground._q(success_lower,'original success lower')
    leak=ground._q(leakage,'kernel leakage');error=ground._q(probability_error,'probability allowance')
    if max(p,leak)>1:raise ValueError('original probability and leakage cannot exceed one')
    return max(F(),p-leak-error)


def weighted_instrument(probes,weights,phases=None):
    """Small exact non-Hermitian coherent selector; every fine outcome survives."""
    probes=tuple(apparatus.exact_matrix(b) for b in probes)
    if not probes or len(probes)>12 or probes[0].rows>4:
        raise ValueError('small executed control requires at most12 probes of dimension at most4')
    d=probes[0].rows;M=len(probes)
    if any(b.shape!=(d,d) for b in probes):raise ValueError('one common square source carrier required')
    weights=tuple(ground._q(p,'positive address weight',True) for p in weights)
    if len(weights)!=M or sum(weights)!=1:raise ValueError('address weights must sum to one')
    phases=tuple(s.S(1) for _ in probes) if phases is None else tuple(map(s.sympify,phases))
    if len(phases)!=M or any(z.is_algebraic is not True or s.simplify(z*s.conjugate(z))!=1 for z in phases):
        raise ValueError('one exact algebraic unit phase per address required')
    v=s.Matrix([s.sqrt(s.Rational(p.numerator,p.denominator)) for p in weights])
    if M==1:W=s.eye(1)
    else:
        e=s.eye(M)[:,0];diff=e-v
        W=apparatus.clean(s.eye(M)-diff*diff.H/(1-v[0]))
    if apparatus.clean(W.H*W)!=s.eye(M) or apparatus.clean(W[0,:]-v.H)!=s.zeros(1,M):
        raise ArithmeticError('coherent weight analyzer failed')
    columns=[(z*b,z*apparatus.positive_sqrt(s.eye(d)-b.H*b)) for b,z in zip(probes,phases)]
    rows=[]
    for branch in (0,1):
        for label in range(M):
            K=apparatus.clean(sum((W[label,j]*v[j]*columns[j][branch] for j in range(M)),s.zeros(d)))
            rows.append({'branch':branch,'analyzer':label,'K':K})
    complete=apparatus.clean(sum((r['K'].H*r['K'] for r in rows),s.zeros(d)))
    expected=apparatus.clean(sum((s.Rational(p.numerator,p.denominator)*z*b
        for p,z,b in zip(weights,phases,probes)),s.zeros(d)))
    if complete!=s.eye(d) or rows[0]['K']!=expected:raise ArithmeticError('complete coherent selector failed')
    return {'rows':rows,'success':expected,'complete_effect':complete,
            'analyzer':W,'address_count':M,'fine_outcomes':2*M}


def binary_instrument(success):
    B=apparatus.exact_matrix(success)
    return (B,apparatus.positive_sqrt(s.eye(B.rows)-B.H*B))


def factored_fourier_instrument(probe_groups,inner_weights,outer_shape):
    """Complete Cartesian outer Fourier rows after each inner weighted analyzer."""
    shape=tuple(ground._i(x,'outer Fourier dimension',1) for x in outer_shape)
    labels=tuple(product(*(range(x) for x in shape)))
    if len(labels)!=len(probe_groups) or len(labels)>8:
        raise ValueError('small control needs the full Cartesian outer grid of at most8 addresses')
    inner=[weighted_instrument(group,inner_weights) for group in probe_groups]
    dimension=inner[0]['success'].rows;count=len(labels);inner_count=len(inner_weights)
    if any(v['success'].shape!=(dimension,dimension) for v in inner):
        raise ValueError('one common source carrier required')
    rows=[]
    for b in (0,1):
        for r in range(inner_count):
            for k in labels:
                K=s.zeros(dimension)
                for j,group in zip(labels,inner):
                    angle=sum((s.Rational(x*y,n) for x,y,n in zip(k,j,shape)),s.S(0))
                    phase=s.expand_complex(s.exp(-2*s.pi*s.I*angle))
                    K+=phase*group['rows'][b*inner_count+r]['K']/count
                rows.append({'branch':b,'inner_analyzer':r,'Fourier_index':k,'K':apparatus.clean(K)})
    complete=apparatus.clean(sum((row['K'].H*row['K'] for row in rows),s.zeros(dimension)))
    if complete!=s.eye(dimension):raise ArithmeticError('Cartesian Fourier analyzer lost records')
    return {'rows':rows,'complete_effect':complete,'fine_outcomes':len(rows),
            'outer_shape':shape,'source_translation_unitary_used':False}


@lru_cache(maxsize=1)
def exact_control():
    B=s.Matrix([[0,s.Rational(1,2)],[0,0]])
    probes=(B,B.T,s.I*B,-B.T);weights=(F(1,2),F(1,4),F(1,8),F(1,8))
    instrument=weighted_instrument(probes,weights,(1,s.I,-1,-s.I))
    rho=s.diag(s.Rational(3,4),s.Rational(1,4));psi=s.Matrix([s.sqrt(3)/2,0,0,s.Rational(1,2)])
    joint=s.zeros(4);probabilities=[]
    for row in instrument['rows']:
        vector=s.kronecker_product(row['K'],s.eye(2))*psi
        joint+=vector*vector.H;probabilities.append(s.simplify((vector.H*vector)[0]))
    reference=s.Matrix(2,2,lambda i,j:sum(joint[2*k+i,2*k+j] for k in range(2)))
    coherent=s.simplify(s.trace(rho*instrument['success'].H*instrument['success']))
    mixture=sum(s.Rational(p.numerator,p.denominator)*s.trace(rho*b.H*b)
        for p,b in zip(weights,probes))
    binary=binary_instrument(instrument['success'])
    if apparatus.clean(reference-rho)!=s.zeros(2) or sum(probabilities)!=1:
        raise ArithmeticError('complete selector lost probability or reference')
    return {'instrument':instrument,'original_success':coherent,'incoherent_success':mixture,
        'probabilities':probabilities,'reference_difference':apparatus.clean(reference-rho),
        'binary_complete':apparatus.clean(sum((k.H*k for k in binary),s.zeros(2))),
        'is_primary_source_execution':False}


def selected_plan(member,spatial_width='1',time_width='1/1000',spatial_aperture='1',
                  time_aperture='1/1000',spatial_index=(1,0,0),time_index=0):
    contract.verify_member(member);n=member['source']['level'];a=F(member['source']['spacing'])
    R=ground._q(spatial_width,'seed spatial width',True);T=ground._q(time_width,'seed time width',True)
    Q=ground._q(spatial_aperture,'spatial half aperture',True);S=ground._q(time_aperture,'time half aperture',True)
    k=_frequency(spatial_index)
    if type(time_index) is not int:raise ValueError('signed integer temporal Fourier index required')
    if (R/a).denominator!=1 or (Q/a).denominator!=1:
        raise ValueError('this finite aperture plan requires grid-aligned spatial widths')
    native=response.native_control()
    if not any(d==6 and charge==1 for d,_,charge in native):raise ArithmeticError('native collective seed changed')
    geometry=origin.origin_geometry(n);bulk=origin.bulk_bound(member,R+Q,T+S,1)
    active=bulk['all_origins_bulk_admitted']
    lo=F(member['state']['spectral_lower']);hi=F(member['state']['spectral_upper'])
    width=hi-lo;norm=max(abs(lo),abs(hi));target=F(1,2**(n+12))
    need=width*T/(2*target);inner_cells=2*max(1,need.__ceil__())
    cayley_steps=max(1,ceil_sqrt(norm**3*(T+S)**3/(12*target)))
    e_U=norm**3*(T+S)**3/(12*cayley_steps**2)
    inner_error=width*T/inner_cells
    outer_cells=max(2,2**n);h=2*S/outer_cells
    grid=a*(6/R+F(22,7)*sum(abs(x) for x in k)/Q)+h*(2/T+F(22,7)*abs(time_index)/S)
    bits=bulk['operator_error']['negative_exponent']
    if bulk['operator_error']['exact_zero']:bulk_usable=F()
    elif bits<0:bulk_usable=F(2)
    else:bulk_usable=F(1,2**min(bits,n+12))
    amp=min(F(2),inner_error+2*e_U+grid+bulk_usable) if active else None
    spatial_nodes=int(2*R/a)**3;offset_nodes=int(2*Q/a)**3
    inner_count=spatial_nodes*inner_cells;outer_count=offset_nodes*outer_cells
    count=inner_count*outer_count
    params={'spatial_width':R,'time_width':T,'spatial_aperture':Q,'time_aperture':S,
        'spatial_index':k,'time_index':time_index}
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'parameters':params,'origin_geometry':geometry,
        'active_source_aperture':active,'inactive_convention':'zero success probe and all failure retained',
        'seed':'native6 Q6=1 collective beta pair divided by sqrt(42)',
        'seed_gauge_charge':0,'seed_full_Fock_norm_upper':1,
        'source_spectral_interval':[lo,hi],'source_norm_upper':norm,
        'inner_spatial_nodes':spatial_nodes,'inner_triangle_cells':inner_cells,
        'outer_spatial_nodes':offset_nodes,'outer_time_nodes':outer_cells,
        'inner_analyzer_dimension':inner_count,'outer_analyzer_dimension':outer_count,
        'selector_dimension':count,'fine_outcomes_per_origin':2*count,
        'declared_analyzer':'inner exact weighted Householder and outer Cartesian spatial/time Fourier analyzers',
        'outer_Fourier_shape':[int(2*Q/a)]*3+[outer_cells],
        'distinguished_Fourier_bin':[x%int(2*Q/a) for x in k]+[time_index%outer_cells],
        'distinguished_row_phase':1 if (sum(k)+time_index)%2==0 else -1,
        'retained_fine_labels':['origin X','probe branch b','inner weighted analyzer r','Cartesian spatial/time Fourier bins'],
        'alternative_weighted_analyzer_has_same_distinguished_success':True,
        'alternative_analyzer_complete_law_claimed_equal':False,
        'complete_fine_outcomes_with_retained_origin':2*count*geometry['origin_count'],
        'Cayley_subdivisions':cayley_steps,'source_energy_zero_changed':False,
        'inner_quadrature_operator_error_upper':inner_error,
        'one_source_evolution_operator_error_upper':e_U,
        'complete_node_instrument_half_trace_error_upper':min(F(1),2*e_U),
        'regular_outer_grid_operator_error_upper':min(F(2),grid),
        'same_source_bulk_bound':bulk,'bulk_error_usable_upper':bulk_usable,
        'fixed_regular_target_comparison_applicable':active,
        'success_amplitude_comparison_upper':amp,
        'original_probability_comparison_upper':min(F(1),2*amp) if active else None,
        'additional_compiler_or_preparation_error_consumed':False,
        'finite_state_to_limit_convergence_modulus_supplied':False,
        'source_reference_and_origin_retained':True,
        'one_origin_across_complete_history':True,
        'joint_limiting_success_kernel':'sinc^2(S*(E-nu))*product_i sinc^2(Q*(p_i-q_i))',
        'Fourier_convention':'fhat(p)=integral f(y)*exp(+i*y.p)dy',
        'frequency_centers':{'q_over_pi':[F(x)/Q for x in k],'nu_over_pi':F(time_index)/S},
        'same_relative_observation_canonical_H_P_used':True,
        'finite_source_translation_symmetry_assumed':False,
        'fixed_origin_state_replaced_or_assumed_invariant':False,
        'finite_Fourier_label_is_exact_momentum_projection':False,
        'binary_failure_instrument_equals_fine_LCU_failure':False,
        'large_source_or_selector_executed':False,
        'nonzero_joint_spectral_weight_proved':False,
        'dispersion_species_or_mass_proved':False,'physical_mass_table_GeV':[]})


def bind_member(member,**kwargs):
    result=selected_plan(member,**kwargs)
    return {'schema':'oph.primary_momentum_transfer.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,receipt,**kwargs):
    if receipt!=bind_member(member,**kwargs):raise ValueError('momentum-transfer plan or source differs')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*origin.SOURCE_INPUTS,*response.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths);control=exact_control()
    B=s.Matrix([[0,s.Rational(1,2)],[0,0]])
    fourier=factored_fourier_instrument(((B,B.T),(B.T,s.I*B),(-B,B.T),(s.I*B,-B.T)),
        (F(3,4),F(1,4)),(2,2))
    result={'schema':'oph.primary_momentum_transfer_certificate.v1',
        'selected_member_bound':bind_member(contract.member(8)),
        'other_members':[bind_member(contract.member(n)) for n in (0,4,12)],
        'exact_control':{'original_success':str(control['original_success']),
            'incoherent_success':str(control['incoherent_success']),
            'complete_probabilities':[str(x) for x in control['probabilities']],
            'arbitrary_reference_completeness_checked':True,'is_primary_source':False},
        'Cartesian_Fourier_control':{'outer_shape':[2,2],'fine_outcomes':fourier['fine_outcomes'],
            'complete_effect_identity':fourier['complete_effect']==s.eye(2),
            'original_probabilities':[str(s.simplify(s.trace(s.diag(s.Rational(2,3),s.Rational(1,3))
                *row['K'].H*row['K']))) for row in fourier['rows']]},
        'rectangle_control':{'leakage':str(kernel_leakage()),
            'original_weight_lower_for_success_half':str(original_rectangle_lower('1/2',kernel_leakage()))},
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('momentum-transfer source changed during replay')
    return result
