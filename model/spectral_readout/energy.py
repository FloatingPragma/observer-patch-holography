"""Exact spectral sum rules and a finite-to-infinite readout error contract."""
from fractions import Fraction as F

from model.thermodynamic_ground.volume import averaged_dynamics_error
from model.thermodynamic_ground.retarded import exp_negative_upper
from .kernels import exact,count,bin_fourier_error


def charge_tag_certificate():
    """Complete addition/removal instrument and positive excitation countercontrol."""
    from model.field_symmetry import symmetry as s
    a=s._annihilator(2,0)+s._annihilator(2,1)
    if not (a.adjoint()@a+a@a.adjoint()-s._scale(s._identity(4),2)).is_zero():
        raise ArithmeticError('normalized addition/removal Kraus pair is incomplete')
    # H=-2n0+3n1 has ground |10>, energy -2. Applying the normalized
    # mode (c0+c1)/sqrt2 removes to vacuum or adds to the doubly occupied
    # state, each with probability 1/2 and positive cost 2 or 3.
    n=[s._annihilator(2,j).adjoint()@s._annihilator(2,j) for j in range(2)]
    h=s._scale(n[0],-2)+s._scale(n[1],3)
    shifted=h+s._scale(s._identity(4),2)
    diagonal=[int(shifted.real[j,j]) for j in range(4)]
    if diagonal!=[2,0,5,3]:raise ArithmeticError('positive excitation control changed')
    removed=[int(a.real[j,1]) for j in range(4)]
    added=[int(a.adjoint().real[j,1]) for j in range(4)]
    if removed!=[1,0,0,0] or added!=[0,0,0,-1]:
        raise ArithmeticError('charge-tagged control lost a branch or its CAR sign')
    return {'addition_Kraus':'c(f)^dagger','removal_Kraus':'c(f)',
            'completeness':'c(f)c(f)^dagger+c(f)^dagger c(f)=I for ||f||=1',
            'physical_readout':'nonnegative excitation energy above the ground plus an addition/removal label',
            'signed_response_readout':'addition:+energy; removal:-energy',
            'two_mode_control':{'positive_excitation_costs':['3','2'],
                                'signed_response_frequencies':['3','-2'],
                                'probabilities':['1/2','1/2']},
            'negative_ground_excitation_energy_inferred':False,
            'native_charge_tagged_full_measurement_executed':False}


def local_moments(speed=F(1,8),interaction=F(1)):
    c=exact(speed,'speed');u=exact(interaction,'interaction')
    kinetic=F(3,2)*c*c
    return {'first_moment':[u/8,u/2],
            'second_moment':[kinetic+u*u/8,kinetic+u*u/2],
            'variance':[kinetic+7*u*u/64,kinetic+u*u/4]}


def energy_tail(radius,second_moment):
    radius=exact(radius,'tail radius',True);moment=exact(second_moment,'second moment')
    return min(F(1),moment/(radius*radius))


def kinetic_controlled_second_moment(kinetic_norm,interaction=F(1)):
    k=exact(kinetic_norm,'kinetic norm');u=exact(interaction,'interaction')
    # r<=1/2 in the selected hidden-SU2 ground and <f,hf><=||hf||.
    return k*k+u*k+u*u/2


def comparison_contract(count_,tau,selected_count,*,half_side,buffer_radius,
                        support_radius,support_sites,speed=F(1,8),
                        native_characteristic_error=F(0),
                        statistical_characteristic_error=F(0),
                        state_discrepancy=None):
    """Complete explicit known errors plus the indispensable state discrepancy.

    State discrepancy is a bound for the characteristic-function difference
    between the averaged finite state and selected infinite state, both with
    infinite dynamics, at all sampled times. None leaves the total unknown.
    Supplying a number states a conditional premise, not its certification.
    """
    m=count(count_);tau=exact(tau,'time step',True)
    time=(m-1)*tau
    dynamics=2*averaged_dynamics_error(half_side,buffer_radius,support_radius,
                                      support_sites,time,speed)
    native=exact(native_characteristic_error,'native characteristic error')
    statistical=exact(statistical_characteristic_error,'statistical characteristic error')
    state=None if state_discrepancy is None else exact(state_discrepancy,'state discrepancy')
    known=min(F(2),dynamics+native+statistical)
    total=None if state is None else bin_fourier_error(m,selected_count,min(F(2),known+state))
    return {'largest_sample_time':str(time),'dynamics_characteristic_error':str(dynamics),
            'native_characteristic_error':str(native),'statistical_characteristic_error':str(statistical),
            'state_discrepancy':None if state is None else str(state),
            'known_bin_error':str(bin_fourier_error(m,selected_count,known)),
            'total_infinite_ground_bin_error':None if total is None else str(total),
            'statistical_event_estimator_clipped_to_unit_interval':True,
            'state_discrepancy_is_supplied_conditional_premise':state is not None,
            'ground_state_convergence_rate_derived':False}


def sampling_failure(count_,shots_per_quadrature,probability_error):
    """Hoeffding union bound for both binary quadratures at M-1 times."""
    m=count(count_)
    if type(shots_per_quadrature) is not int or shots_per_quadrature<=0:
        raise ValueError('positive integer shot count required')
    e=exact(probability_error,'binary probability error',True)
    if e>=1:raise ValueError('binary probability error must be below one')
    return min(F(1),4*(m-1)*exp_negative_upper(2*shots_per_quadrature*e*e))


def certificate():
    moments=local_moments()
    return {'measure':'mu_f(B)=<Vf,1_B(L)Vf>, ||f||=1; signed addition/removal energy variable',
            'charge_tagged_interpretation':charge_tag_certificate(),
            'zeroth_operator_moment':'I',
            'first_operator_moment':'h+U r I',
            'second_operator_moment':'h^2+2U r h+U^2 r I',
            'density_per_spin_interval':['1/8','1/2'],
            'local_probe_moments':{k:[str(x) for x in v] for k,v in moments.items()},
            'local_probe_tail_examples':[{'energy_radius':str(r),'probability_upper':str(energy_tail(r,moments['second_moment'][1]))}
                                         for r in (1,2,4,8)],
            'kinetic_controlled_second_moment':'||h f||^2+U ||h f||+U^2/2',
            'uniform_kinetic_example':str(kinetic_controlled_second_moment(1,1)),
            'native_characteristic_rule':'Binary real/imaginary Ramsey probability errors <=e give complex characteristic error <=3e, since 2sqrt(2)<3.',
            'sampling_rule':'With N independent shots for each of 2(M-1) binary quadratures, failure<=4(M-1) exp(-2N e^2). No such shots are generated here.',
            'sampling_example_failure_upper':str(sampling_failure(16,100000,F(1,50))),
            'weak_spectral_convergence_same_ground_subsequence':True,
            'strong_compressed_resolvent_convergence_at_fixed_spacing':True,
            'common_Hilbert_space_strong_resolvent_convergence_of_dilations_claimed':False,
            'infinite_ground_spectral_values_computed':False,
            'sharp_interval_convergence_requires_zero_endpoint_mass':True,
            'real_axis_density_or_isolated_pole_inferred':False,
            'relativistic_rest_mass_inferred':False}
