"""Complete pair-frame detector and connected finite source responses."""
from fractions import Fraction as F
from math import factorial
import numpy as np
from model.centered_particle_spectrum.certificate import sqrt_bounds,decimal_bounds
from . import algebra


def squared_modulus(re,im,error):
    lower,upper=sqrt_bounds(re*re+im*im)
    return max(F(),lower-error)**2,min(F(1),(upper+error)**2)


def pair_readout(q,t=F(1,4),order=32):
    if q not in(-1,1) or type(t) is not F or t<0 or type(order) is not int or not 1<=order<=64:
        raise ValueError('signed pair charge, rational time and bounded order required')
    op,g,norm,error,checked,creation,sectors=algebra.data()
    sector=sectors[q];c=sector['coordinates'];den=sector['denominator']
    incoming=algebra.PAIRS.index((0,7));vector=c[:,incoming].copy()
    hc=sector['block']+6*op['scale']*np.eye(28,dtype=object)
    re=[F() for _ in range(28)];im=[F() for _ in range(28)]
    for n in range(order+1):
        coefficients=c.T@vector
        factor=t**n/F(den*op['scale']**n*factorial(n))
        for j in range(28):
            term=int(coefficients[j])*factor
            if n%2==0:re[j]+=(-1)**(n//2)*term
            else:im[j]-=(-1)**((n-1)//2)*term
        vector=hc@vector
    eps=(7*t)**(order+1)/factorial(order+1)+2*error
    probabilities=[tuple(v/28 for v in squared_modulus(r,i,eps)) for r,i in zip(re,im)]
    wl=F(int(sector['gram_numerator'][incoming,incoming]),den)-2*error
    wh=wl+4*error
    dark=(max(F(),wl-sum(p[1] for p in probabilities)),min(F(1),wh-sum(p[0] for p in probabilities)))
    if dark[0]<0 or dark[0]>dark[1]:raise ArithmeticError('complete pair detector lost its complement')
    rows=[]
    for j,(p,r,i) in enumerate(zip(probabilities,re,im)):
        rows.append({'detector_pair':list(algebra.PAIRS[j]),
            'original_input_herald_and_detector_probability':decimal_bounds(*p,14),
            'unnormalized_kernel_real_interval':decimal_bounds(r-eps,r+eps,12),
            'unnormalized_kernel_imaginary_interval':decimal_bounds(i-eps,i+eps,12)})
    return {'total_charge':2*q,'source_time':str(t),'input_sites':[0,7],
        'same_charge_input_probability_interval':decimal_bounds(wl,wh,12),
        'all28_detector_outcomes':rows,
        'same_charge_input_and_dark_probability':decimal_bounds(*dark,14),
        'each_other_preparation_branch_has_zero_pair_click_probability':True,
        'all_preparation_and_detector_coarse_records':116,
        'conditional_total_pair_click_probability_lower':str(sector['true_gram_floor']/28),
        'detector_pair_label_claimed_as_two_particle_positions':False,
        'physical_cross_section':None}


def unitary_polynomial_vector(matrix,vector,scale,t,order):
    """Integer real and imaginary numerators over one exact denominator."""
    if type(t) is not F or t<0 or type(order) is not int or not 1<=order<=64:
        raise ValueError('nonnegative rational time and bounded degree required')
    denominator=(scale*t.denominator)**order*factorial(order)
    re=np.zeros(len(vector),dtype=object);im=np.zeros(len(vector),dtype=object)
    power=np.array(vector,dtype=object)
    for n in range(order+1):
        coefficient=t.numerator**n*(scale*t.denominator)**(order-n)*(factorial(order)//factorial(n))
        if n%2==0:re+=(-1)**(n//2)*coefficient*power
        else:im-=(-1)**((n-1)//2)*coefficient*power
        power=matrix@power
    return re,im,denominator


def connected_response(y,t=F(1,8),order=24):
    if type(y) is not int or not 1<=y<8:raise ValueError('distinct actual source sites required')
    op,g,norm,error,checked,creation,sectors=algebra.data()
    scale=op['scale'];indices=np.arange(256);identity=np.eye(256,dtype=object)
    zx=np.diag(1-2*(indices&1)).astype(object)
    zy=np.diag(1-2*((indices>>y)&1)).astype(object)
    h=op['H'].astype(object);hc=h+6*scale*identity
    product=zx@zy
    mixed=product@h@product-zx@h@zx-zy@h@zy+h
    edge=(0,y) in algebra.source.cube()[0]
    xx=identity[indices^1];yy=identity[indices^(1<<y)]
    expected=-scale*xx@yy@(identity-zx@zy) if edge else np.zeros_like(h)
    if not np.array_equal(mixed,expected):raise ArithmeticError('mixed-kick interaction differs from actual edge term')
    double=(h@zx-zx@h)@zy-zy@(h@zx-zx@h)
    if edge and not np.array_equal(double,-mixed):raise ArithmeticError('connected-response sign or bond normalization differs')
    mean=F(int(g@mixed@g),scale*norm);mean_error=4*error if edge else F()
    distance=y.bit_count();nested=zx.copy()
    for n in range(distance):
        if np.any(nested@zy-zy@nested):raise ArithmeticError('source causal support arrived before graph distance')
        nested=h@nested-nested@h
    ur,ui,den=unitary_polynomial_vector(hc,g,scale,t,order)
    vr,vi,den2=unitary_polynomial_vector(hc,zy@g,scale,t,order)
    assert den==den2
    imaginary=F(int(ur@zx@vi-ui@zx@vr),norm*den**2)
    rem=(7*t)**(order+1)/factorial(order+1)
    eps=2*rem+rem**2+2*error
    chi=(2*(imaginary-eps),2*(imaginary+eps))
    pointer=(F(1,2)+(imaginary-eps)/2,F(1,2)+(imaginary+eps)/2)
    return {'sites':[0,y],'graph_distance':distance,'actual_source_edge':edge,
        'mixed_Z_kick_energy_interval':decimal_bounds(mean-mean_error,mean+mean_error,12),
        'mixed_kick_operator_norm_bound':'2' if edge else '0',
        'commutator_orders_before_distance_exactly_zero':list(range(distance)),
        'source_time':str(t),'retarded_Z_response_interval':decimal_bounds(*chi,12),
        'complete_interference_Y_plus_probability':decimal_bounds(*pointer,12),
        'complete_interference_other_probability':decimal_bounds(1-pointer[1],1-pointer[0],12),
        'finite_connected_response_claimed_as_asymptotic_scattering':False}


def report():
    responses=[connected_response(y) for y in(1,3,7)]
    if F(responses[0]['mixed_Z_kick_energy_interval'][1])>=-1:
        raise ArithmeticError('actual nonzero bond interaction not certified')
    if any(F(r['retarded_Z_response_interval'][0])<=0<=F(r['retarded_Z_response_interval'][1]) for r in responses):
        raise ArithmeticError('source propagation response not separated from zero')
    return {'all_pair_frame_readouts':[pair_readout(q,t) for q in(-1,1) for t in(F(),F(1,4))],
        'complete_connected_response_controls':responses,
        'scope':{'same_selected_ground_and_complete_source_H':True,
            'full_two_charge_sector_dynamics_and_observables_constructed':True,
            'two_charge_creation_has_all_four_outcomes':True,
            'complete_pair_frame_POVM_has_29_outcomes':True,
            'positive_source_interaction_and_finite_propagation_response':True,
            'native_arbitrary_reference_instrument_compiler_reused':True,
            'pair_detector_and_connected_response_native_tapes_executed':False,
            'scalable_two_particle_Fock_identification':False,
            'asymptotic_scattering_on_this_selected_source_proved':False,
            'other_Weyl_contact_source_S_matrix_transferred':False,
            'physical_species_mass_or_cross_section':None}}
