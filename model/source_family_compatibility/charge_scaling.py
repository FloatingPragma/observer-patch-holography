"""Actual centered-source charged bands, uniform in source volume.

Bounds retain original response weights. They imply no hard many-body gap,
particle species, or arbitrary-input channel approximation.
"""
from fractions import Fraction as F
from functools import lru_cache

import numpy as np

from . import composite


def response_budget(hopping, interaction, degree, energy_window, band_width=None,
                    scaled_time=F(0)):
    """Squared rational bounds on the same source response and clock."""
    c,u,e,t=map(F,(hopping,interaction,energy_window,scaled_time))
    if c<=0 or u<=0 or e<0 or type(degree)is not int or degree<1:
        raise ValueError('positive source coefficients, integer degree, and nonnegative window required')
    alpha=u/2
    beta2=c*c*degree*degree*(F(101,4)+25*degree*degree)
    gamma2=25*c*c*degree*degree*(1+degree*degree)
    low=F(1) if e>=alpha else min(F(1),beta2/(alpha-e)**2)
    result={'alpha':alpha,'residual_Gram_upper':beta2,
            'original_total_fixed_window_weight_upper':low,
            'original_local_each_charge_window_weight_upper':min(F(1,2),low),
            'interaction_clock_vector_error_squared_upper':min(F(4),t*t*gamma2/alpha**2),
            'hard_many_body_gap_derived':False,'physical_mass':None,
            'arbitrary_input_channel_bound':False}
    if band_width is not None:
        w=F(band_width)
        if w<=0:raise ValueError('positive band radius required')
        result['original_total_outside_two_bands_weight_upper']=min(F(1),beta2/(w*w))
        result['original_local_each_charge_inside_band_weight_lower']=max(F(0),F(1,2)-beta2/(w*w))
    return result




def charge_clock_shift_budget(hopping, interaction, degree, energy_window, chemical_shift):
    """At least one conjugate local branch increases under K-mu N."""
    e,mu=map(F,(energy_window,chemical_shift))
    # Validate the unshifted source even if the selected interval is empty.
    response_budget(hopping,interaction,degree,e)
    original_upper=e-abs(mu)
    weight=F(0) if original_upper<0 else response_budget(
        hopping,interaction,degree,original_upper)['original_local_each_charge_window_weight_upper']
    return {'increased_branch':'removal' if mu>=0 else 'addition',
            'chemical_shift':mu,'original_window_upper':original_upper,
            'original_local_selected_branch_window_weight_upper':weight,
            'shifted_reference_is_assumed_ground':False}


def cell_contact_budget(spacing, speed, contact, degree, energy_window):
    """Exact cell normalization g=a^3 U, retaining the centered clock."""
    a,v,g,e=map(F,(spacing,speed,contact,energy_window))
    if a<=0 or v<=0 or g<=0:
        raise ValueError('positive spacing, speed, and cell contact required')
    row=response_budget(v/a,g/a**3,degree,e)
    denominator=g/2-e*a**3
    explicit=F(1) if denominator<=0 else min(F(1),
        v*v*degree*degree*(F(101,4)+25*degree*degree)*a**4/denominator**2)
    if explicit!=row['original_total_fixed_window_weight_upper']:
        raise ArithmeticError('cell-normalized contact rate differs from source bound')
    return {'spacing':a,'speed':v,'cell_quartic_coefficient':g,
            'c':v/a,'U':g/a**3,'centering_coefficient':g/(2*a**3),
            'original_total_fixed_window_weight_upper':explicit}


def _matrix_fraction(matrix,scale):
    return tuple(tuple(F(int(x),scale)for x in row)for row in matrix)


def _add(a,b):
    return tuple(tuple(x+y for x,y in zip(ar,br,strict=True))for ar,br in zip(a,b,strict=True))


def _gram(rho,operators,scale):
    return tuple(tuple(composite._trace(rho,composite._anti(x.T,y))/(scale*scale)
                       for y in operators)for x in operators)


@lru_cache(maxsize=16)
def dimer(n,scale=F(1,2)):
    """Exact full-source projectors and response residuals for a rational family."""
    if type(n)is not int or n<3:raise ValueError('integer n>=3 required')
    r=F(scale)
    if r<=0:raise ValueError('positive physical source scale required')
    c,u,e0=2*n*r,(n*n-4)*r,-n*n*r
    data=composite.source((2,1,1),c,u)
    numerator,denominator,rank=composite.ground(data,e0)
    rho=(numerator,denominator);h=data['H'];s=data['scale'];a,b=data['a'],data['b']
    identity=np.eye(len(h),dtype=np.int64)
    alpha=u/2;L=h-int(s*e0)*identity
    energies=(r*(F(n*n,2)-n+2),r*(F(n*n,2)+n+2))
    weights=((1+F(4*n,n*n+4))/4,(1-F(4*n,n*n+4))/4)
    # Complete projectors within the actual N=1 and N=3 occupation sectors.
    # Preparation failure is retained separately, without conditioning.
    outcomes={}
    for charge,operator,number in (('addition',a[0].T,3),('removal',a[0],1)):
        sector=np.diag([int(mask.bit_count()==number)for mask in range(len(h))])
        projectors=[];probabilities=[]
        for index,energy in enumerate(energies):
            other=energies[1-index]
            # Rational numerator/denominator, checked before using the probability.
            pn=composite._mul(sector,L-int(s*other)*identity)
            pd=F(s*(energy-other))
            if pd.denominator!=1:raise ArithmeticError('energy projector scale mismatch')
            pd=pd.numerator
            if np.any(composite._mul(pn,pn)-pd*pn):
                raise ArithmeticError('actual charge energy projection is not idempotent')
            if np.any(composite._mul(L,pn)-int(s*energy)*pn):
                raise ArithmeticError('actual source energy projection has the wrong energy')
            value=composite._trace(rho,composite._mul(composite._mul(operator.T,pn),operator))/pd
            if value!=weights[index]:raise ArithmeticError('original charge probability changed')
            projectors.append((pn,pd));probabilities.append(value)
        if np.any(projectors[0][0]*projectors[1][1]+projectors[1][0]*projectors[0][1]
                  -sector*projectors[0][1]*projectors[1][1]):
            raise ArithmeticError('complete charge-sector energy outcomes do not sum to identity')
        if sum(probabilities,F())!=F(1,2):raise ArithmeticError('local charge weight is not one half')
        # The full complement comprises rejected occupation sectors and the
        # original source preparation failure. No zero branch is deleted.
        outcomes[charge]={'energy_atoms':energies,'original_probabilities':tuple(probabilities),
                          'preparation_failure_probability':F(1,2),
                          'other_energy_outcomes_probability':F(0),
                          'total_probability':F(1)}
    one=_matrix_fraction(data['one'],s)
    # Retarded derivative columns are the actual odd source commutators.
    # (Lret-sign alpha)(V0+sign V1) is represented by the corresponding
    # adjoint odd source operator; divide the resulting Gram by two.
    Y=tuple(composite._comm(h,x.T)-int(s*alpha)*y.T for x,y in zip(b,a,strict=True))
    Ygram=_gram(rho,Y,s)
    target=_add(composite._square(one),Ygram)
    for sign in (-1,1):
        operators=tuple(composite._comm(h,x.T+sign*y.T)
                        -int(s*sign*alpha)*(x.T+sign*y.T)
                        for x,y in zip(a,b,strict=True))
        actual=_gram(rho,operators,s)
        if actual!=target:
            raise ArithmeticError('source positive/negative residual Gram identity failed')
    exact_beta=max(sum(abs(x)for x in row)for row in target)
    if exact_beta!=F(5,2)*c*c:
        raise ArithmeticError('dimer source residual normalization failed')
    universal=response_budget(c,u,1,F(1))
    if exact_beta>universal['residual_Gram_upper']:
        raise ArithmeticError('universal bound misses actual source Gram')
    return {'n':n,'source_scale':r,'c':c,'U':u,'ground_energy':e0,
            'whole_Fock_dimension':len(h),'positive_shift_exact_rank':rank,
            'ground_degeneracy':len(h)-rank,'ground_density_denominator':denominator,
            'alpha':alpha,'h2_plus_YGram':target,'exact_residual_Gram_norm_upper':exact_beta,
            'complete_local_charge_records':outcomes,
            'fixed_energy_window_1_bound':universal['original_total_fixed_window_weight_upper'],
            'hard_gap_inferred_from_bound':False,'physical_particle_assignment':None}



def homogeneous_dimer(n,first_scale=F(1,2),second_scale=F(3,2)):
    """Exact fixed-ratio source scaling, including its full ground projector."""
    first=dimer(n,F(first_scale));second=dimer(n,F(second_scale))
    factor=F(second_scale)/F(first_scale)
    sources=[composite.source((2,1,1),row['c'],row['U'])for row in(first,second)]
    one,two=sources
    if np.any(factor.denominator*two['H'].astype(object)*one['scale']
              -factor.numerator*one['H'].astype(object)*two['scale']):
        raise ArithmeticError('full centered source fails homogeneous scaling')
    grounds=[composite.ground(data,row['ground_energy'])for data,row in zip(sources,(first,second),strict=True)]
    g1,g2=grounds
    if np.any(g1[0].astype(object)*g2[1]-g2[0].astype(object)*g1[1])or g1[2]!=g2[2]:
        raise ArithmeticError('homogeneous source scaling changed the full ground projector')
    for charge in ('addition','removal'):
        a=first['complete_local_charge_records'][charge];b=second['complete_local_charge_records'][charge]
        if b['original_probabilities']!=a['original_probabilities']:
            raise ArithmeticError('homogeneous source scaling changed original response weights')
        if b['energy_atoms']!=tuple(factor*energy for energy in a['energy_atoms']):
            raise ArithmeticError('original charge energies fail the exact pushforward')
    if second['exact_residual_Gram_norm_upper']!=factor**2*first['exact_residual_Gram_norm_upper']:
        raise ArithmeticError('source response residual has the wrong homogeneous degree')
    return {'n':n,'lambda':first['U']/first['c'],'positive_clock_factor':factor,
            'first_c':first['c'],'second_c':second['c'],
            'same_complete_ground_projector':True,'same_original_charge_weights':True,
            'charge_energy_pushforward_checked':True,'residual_Gram_scales_quadratically':True,
            'graph_or_ground_subsequence_changed':False,
            'fixed_ratio_volume_uniform_gap_inferred':False}


def certificate():
    ballistic=[dimer(n)for n in (3,4,8)]
    atomic=[dimer(n,F(1,n*n-4))for n in (3,4,8)]
    # This exact family checks the asymptotic error constant independently
    # of the small finite oracle; no floating underflow certifies a limit.
    degree=6
    family=[]
    for n in (1024,4096,16384):
        row=response_budget(n,F(n*n-4,2),degree,1,scaled_time=1)
        family.append({'n':n,'spacing':F(1,n),**row})
    if not all(family[i+1]['original_total_fixed_window_weight_upper']
               <family[i]['original_total_fixed_window_weight_upper']for i in range(len(family)-1)):
        raise ArithmeticError('ballistic strong-contact fixed-window bound fails to decrease')
    return {'source':'same centered full invariant ground and original CAR probes',
            'volume_uniform':True,'requires_total_edge_norm_projection':False,
            'homogeneous_source_controls':[homogeneous_dimer(n)for n in(3,4,8)],
            'ballistic_strong_contact_controls':ballistic,'atomic_controls':atomic,
            'degree_six_uniform_bounds':family,
            'fixed_cell_contact_bounds':[cell_contact_budget(F(1,n),1,1,6,1)for n in (32,64,128)],
            'charge_clock_shift_controls':[charge_clock_shift_budget(1,1000,6,1,mu)for mu in (-2,F(-1,2),0,F(1,2),2)],
            'all_smearings_have_exact_separate_half_weights':False,
            'physical_species_assigned':0,'physical_masses_calculated':False,
            'hard_many_body_gap_derived':False,'arbitrary_input_channel_bound':False}
