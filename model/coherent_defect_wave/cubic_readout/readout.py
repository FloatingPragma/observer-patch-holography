"""Complete eight-momentum and two-quadrature phase records."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from .lattice import State,ZERO,plus,scale,packet,tick,translate,overlap,observer,chi,coordinates,mode
from .fourier import compress,phase_branch


def multiply(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def quarter_phase(n):
    return ((F(1),F(0)),(F(0),F(1)),(-F(1),F(0)),(F(0),-F(1)))[n%4]


def fixed_link_wave_control(initial,length):
    """Remove the five repair controls, keeping the same native wave and hops."""
    from ..model import SparseState
    from .lattice import STAR,hopping_layer
    out={}
    for (position,index,field),value in initial.amplitudes:
        if field:
            raise ValueError('this control specifies the face-zero reference input')
        w=coordinates(index)
        local=SparseState({(1,*w):F(1),(1^STAR,*w):F(chi(index))},2)
        advanced=observer().wave_step(local)
        for (g,p,j,r,a),coefficient in advanced.amplitudes:
            if g==1:
                key=position,mode(p,j,r,a),()
                out[key]=plus(out.get(key,ZERO),scale(value,F(coefficient)))
    state=State(out,initial.norm_denominator)
    for axis,parity in product(range(3),range(2)):
        state=hopping_layer(state,length,axis,parity)
    return state


@lru_cache(maxsize=1)
def higher_momentum_experiment():
    """Actual side-eight tick, all64 coarse momenta, complete exact kernels."""
    initial=packet();advanced=tick(initial,8)
    if advanced.norm_squared()!=1:
        raise ArithmeticError('side-eight actual tick lost probability')
    displacements=tuple(product(range(4),repeat=3))
    correlations={d:overlap(translate(initial,8,tuple(2*x for x in d)),advanced)
                  for d in displacements}
    momenta=[];records=[];cosines=(1,0,-1,0)
    for k in displacements:
        actual=ZERO
        for d,value in correlations.items():
            phase=quarter_phase(-sum(a*b for a,b in zip(k,d,strict=True)))
            actual=plus(actual,multiply(phase,value))
        predicted=(-F(29,30),F(0))
        for label in k:
            cosine=cosines[label]
            factor=(F(1+cosine,2),F(cosine-1,2))
            predicted=multiply(predicted,factor)
        if actual!=predicted:
            raise ArithmeticError('actual autocorrelations disagree with the derived Bloch first moment')
        variance=1-actual[0]**2-actual[1]**2
        if variance<F(59,900):
            raise ArithmeticError('source phase variance bound changed')
        momenta.append({'momentum_in_pi_over_two_units':list(k),
                        'probability_given_positive_Gauss':'1/64',
                        'complex_first_moment':[str(x) for x in actual],
                        'unitary_variance':str(variance)})
        for q in (0,1):
            for sign in (1,-1):
                component=actual[0] if q==0 else -actual[1]
                weight=(1+sign*component)/256
                records.append({'momentum_in_pi_over_two_units':list(k),
                                'quadrature':'cos' if q==0 else 'minus_sin','sign':sign,
                                'probability_given_positive_Gauss':str(weight),
                                'probability_per_original_ambient_trial':str(weight/2**512)})
    if sum(F(row['probability_given_positive_Gauss']) for row in records)!=1:
        raise ArithmeticError('complete finer-momentum records lost normalization')
    moved=sum((a*a+b*b for (_,_,field),(a,b) in advanced.amplitudes if field),F(0))/advanced.norm_denominator
    if moved!=F(1,450):
        raise ArithmeticError('larger-volume defect response differs')
    return {'periodic_side_cells':8,'spatial_cell_count':512,
            'source_wave_occupation':1,'full_tick_amplitude_support':len(advanced.amplitudes),
            'full_tick_norm_squared':'1','coarse_translation_stride_cells':2,
            'coarse_momentum_count':64,'complete_record_count':256,
            'evaluation':'exact source autocorrelations and proved complete polynomial kernels',
            'complete_branch_vectors_materialized':False,
            'source_Bloch_first_moment':'(-29/30) product_i [(1-i)/2 + (1+i) cos(q_i)/2]',
            'this_first_moment_is_not_an_eigenphase_dispersion':True,
            'momenta':momenta,'complete_joint_records':records,
            'probability_any_defect_changes_face':'1/450',
            'Gauss_preparation':{'retained_flag_count':2**512,
                                 'probability_each_flag_and_all_positive_success':str(F(1,2**512)),
                                 'retained_other_flags_total_probability':str(1-F(1,2**512))}}


@lru_cache(maxsize=1)
def experiment():
    initial=packet(); advanced=tick(initial,4)
    if initial.norm_squared()!=1 or advanced.norm_squared()!=1:
        raise ArithmeticError('full cubic tick lost probability')
    if tick(advanced,4,True)!=initial:
        raise ArithmeticError('full cubic inverse failed')
    correlations={d:overlap(translate(initial,4,tuple(2*x for x in d)),advanced)
                  for d in product(range(2),repeat=3)}
    rows=[]; momenta=[]
    for k in product(range(2),repeat=3):
        first=ZERO
        for d,c in correlations.items():
            first=plus(first,scale(c,(-1)**sum(a*b for a,b in zip(k,d,strict=True))))
        projected=compress(initial,k)
        evolved=compress(advanced,k)
        weight=projected.norm_squared()
        if weight!=F(1,8) or evolved.norm_squared()!=weight:
            raise ArithmeticError('momentum projector or covariance lost probability')
        if overlap(projected,evolved)!=scale(first,weight):
            raise ArithmeticError('actual projected vectors disagree with autocorrelation moments')
        # A complete source execution supplies the moment, rather than a fitted dispersion.
        r=sum(k)
        expected=((-F(29,30),F(0)),(F(0),F(29,30)),
                  (F(29,30),F(0)),(F(0),-F(29,30)))[r]
        if first!=expected:
            raise ArithmeticError('momentum-dependent complex response changed')
        for q in (0,1):
            for sign in (1,-1):
                branch=phase_branch(projected,evolved,q,sign)
                probability=branch.norm_squared()
                component=first[0] if q==0 else -first[1]
                expected_probability=weight*(1+sign*component)/4
                if probability!=expected_probability:
                    raise ArithmeticError('executed branch does not equal its phase kernel')
                rows.append({'momentum_pi_units':list(k),'quadrature':'cos' if q==0 else 'minus_sin',
                             'sign':sign,'probability_given_positive_Gauss':str(probability),
                             'probability_given_momentum_and_positive_Gauss':str(probability/weight),
                             'orthonormal_orbit_basis_support':len(branch.amplitudes),
                             'full_cell_basis_support':8*len(branch.amplitudes)})
        momenta.append({'momentum_pi_units':list(k),'probability_given_positive_Gauss':str(weight),
                        'complex_first_moment':[str(x) for x in first],
                        'unitary_variance':'59/900','is_exact_frequency_eigenstate':False})
    if sum(F(row['probability_given_positive_Gauss']) for row in rows)!=1:
        raise ArithmeticError('complete joint instrument is not normalized')
    defect_weight=sum((a*a+b*b for (_,_,field),(a,b) in advanced.amplitudes if field),F(0))/advanced.norm_denominator
    if defect_weight!=F(1,450):
        raise ArithmeticError('actual cubic experiment did not retain its moving-field response')
    # A distinct no-onsite control must remove the measured source disturbance.
    from .lattice import hopping_layer
    no_source=initial
    for axis,parity in product(range(3),range(2)):
        no_source=hopping_layer(no_source,4,axis,parity)
    if any(field for (_,_,field),_ in no_source.amplitudes):
        raise ArithmeticError('hopping-only control moved a link orbit')
    no_source_correlations={d:overlap(translate(initial,4,tuple(2*x for x in d)),no_source)
                            for d in product(range(2),repeat=3)}
    for k in product(range(2),repeat=3):
        response=ZERO
        for d,c in no_source_correlations.items():
            response=plus(response,scale(c,(-1)**sum(a*b for a,b in zip(k,d,strict=True))))
        if response[0]**2+response[1]**2!=1:
            raise ArithmeticError('hopping-only countercontrol is not its exact phase state')
    frozen=fixed_link_wave_control(initial,4)
    if frozen.norm_squared()!=1 or any(field for (_,_,field),_ in frozen.amplitudes):
        raise ArithmeticError('fixed-link comparison did not preserve its source domain')
    frozen_correlations={d:overlap(translate(initial,4,tuple(2*x for x in d)),frozen)
                         for d in product(range(2),repeat=3)}
    if frozen_correlations!=correlations:
        raise ArithmeticError('fixed-link comparison no longer matches these first-moment records')
    return {'periodic_side_cells':4,'spatial_cell_count':64,
            'source_occupation':1,'source_internal_wave_mode':[7,3,0,0],
            'initial_packet_cells':[[*x] for x in product(range(2),repeat=3)],
            'initial_amplitude_support':len(initial.amplitudes),
            'full_tick_amplitude_support':len(advanced.amplitudes),
            'full_tick_norm_squared':'1','full_tick_inverse_exact':True,
            'all_momentum_outputs_use_the_proved_translation_eigenbasis':True,
            'complete_joint_records':rows,'momenta':momenta,
            'translation_autocorrelations':[
                {'supercell_displacement':list(d),'value':[str(x) for x in c]}
                for d,c in correlations.items()],
            'probability_any_defect_changes_face':'1/450',
            'phase_concentration':{'half_width_radians':'pi/3',
                                   'conditional_weight_lower_bound_each_momentum':'14/15',
                                   'center':'argument of the executed complex first moment',
                                   'source_spectral_measure_not_clock_bin_frequency':True},
            'hopping_only_control':{'defect_change_probability':'0',
                                   'each_tested_momentum_state_is_phase_eigenstate':True},
            'fixed_link_wave_control':{'same_first_moment_autocorrelations':True,
                                       'same_complete32_momentum_quadrature_record_probabilities':True,
                                       'defect_change_probability':'0',
                                       'the_first_moment_instrument_does_not_by_itself_certify_interaction':True},
            'preparation':{'individual_local_Gauss_outcome_probabilities':['1/2','1/2'],
                           'retained_joint_flag_count':2**64,
                           'every_joint_flag_probability':str(F(1,2**64)),
                           'all_positive_success_probability':str(F(1,2**64)),
                           'retained_other_flags_total_probability':str(1-F(1,2**64)),
                           'original_trial_joint_record_weights':[
                               str(F(row['probability_given_positive_Gauss'])/2**64) for row in rows]},
            'physical_scope':{'three_dimensional_finite_momentum_phase_readout':True,
                              'same_spatial_CAR_clock_executed':True,
                              'physical_energy_selected':False,'rest_mass_extracted':False,
                              'phase_bins_are_energy_eigenprojectors':False,
                              'relativistic_dispersion_proved':False,'Standard_Model_species_identified':False}}
