"""Explicit finite resources for local histories in the infinite spin dynamics."""
from fractions import Fraction as F
from . import circuits,precision


def rational(value,positive=False):
    if type(value) not in (int,F) or value<0 or (positive and value==0):raise ValueError('nonnegative exact rational required')
    return F(value)


def ceil(value):return -(-value.numerator//value.denominator)


def cutoff(radius,support_radius,support_sites,coupling,duration):
    if type(radius) is not int or type(support_radius) is not int or not 0<=support_radius<radius:
        raise ValueError('strictly nested nonnegative integer cube radii required')
    if type(support_sites) is not int or not 1<=support_sites<=(2*support_radius+1)**3:raise ValueError('valid support count required')
    j,t=rational(coupling,True),rational(duration)
    if not t:return F()
    # e^x<=2^ceil(2x), since log 2>1/2. There is no fixed-time cap.
    exp_upper=1<<ceil(96*j*t)
    return min(F(2),6*j*support_sites*(2*radius+1)**2*t*F(1,1<<(radius-support_radius))*(exp_upper-1))


def history_cutoff(radius,support_radius,support_sites,coupling,duration,interventions):
    if type(interventions) is not int or interventions<0:raise ValueError('nonnegative intervention count required')
    return min(F(2),(2*interventions+1)*cutoff(radius,support_radius,support_sites,coupling,duration))


def resource_plan(*,support_radius=0,support_sites=1,coupling=F(1),duration=F(1,32),
                  interventions=2,repairs=1,probes=1,rotations=0,error=F(1,10)):
    j,t,eps=rational(coupling,True),rational(duration),rational(error,True)
    if eps>=1:raise ValueError('a requested error strictly below one is required')
    counts=(interventions,repairs,probes,rotations)
    if any(type(v) is not int or v<0 for v in counts) or repairs+probes+rotations!=interventions:
        raise ValueError('nonnegative intervention counts with complete budget required')
    radius=support_radius+1
    while history_cutoff(radius,support_radius,support_sites,j,t,interventions)>eps/3:radius+=1
    side=2*radius+1;vertices=side**3;edges=3*(side-1)*side**2
    adjacent_pairs=15*vertices
    steps=ceil(6*adjacent_pairs*j*j*t*t/eps) if t else 0
    total_steps=steps+interventions+1 if t else 0
    count=circuits.pulse_counts()
    # A specified Neel product is prepared from blanks, not queried as a ground.
    preparation=7*vertices
    operation=repairs*(count['E_repair_adjacent_sites']+88*max(vertices-2,0))+probes*max(count['single_'+a+'_probe'] for a in 'XYZ')+11*rotations
    pulses=preparation+count['exchange_bond']*edges*total_steps+operation
    bits=precision.bits_for(pulses,eps/3)
    product_error=F(2)*adjacent_pairs*j*j*t*t/steps if steps else F()
    native_error=2*pulses*sum(precision.generator_norms())*F(1,1<<bits)
    locality=history_cutoff(radius,support_radius,support_sites,j,t,interventions)
    if locality+product_error+native_error>eps:raise ArithmeticError('resource budget failed')
    return {'support_radius':support_radius,'support_sites':support_sites,'duration':t,'exchange':j,
            'interventions':interventions,'repair_operations':repairs,'Pauli_probes':probes,
            'spin_rotations':rotations,'requested_trace_error':eps,'radius':radius,'sites':vertices,
            'edges':edges,'unordered_overlapping_bond_pairs_upper':adjacent_pairs,
            'base_time_steps':steps,'total_substeps_upper':total_steps,
            'native_pulses_upper':pulses,'port_fraction_bits':bits,
            'initial_product_preparation_pulses_upper':preparation,
            'local_history_cutoff_error':locality,'product_formula_trace_error':product_error,
            'native_complete_trace_error':native_error,'total_trace_error':locality+product_error+native_error,
            'input_state_trace_error_is_added_separately':True,
            'schedule':'fixed listed native words; adaptive control-circuit overhead must be added separately',
            'arbitrary_adaptive_control_synthesis_budgeted':False,
            'initial_state':'specified Neel product; arbitrary consistent supplied marginal also allowed with added trace error',
            'infinite_ground_prepared':False,'resource_plan_native_tape_executed':False}


def encode(value):
    if type(value) is F:return str(value)
    if type(value) is dict:return {k:encode(v) for k,v in value.items()}
    if type(value) in (list,tuple):return [encode(v) for v in value]
    return value


def certificate():
    return encode({'local_history_bound':'(2*m+1)*B_R(T), uniform in reference dimension',
                   'proof':'interaction-picture product of m local Stinespring unitaries and one final local observable',
                   'first_order_product_formula_trace_bound':'2*A*J^2*T*delta, A<=15*|V|',
                   'gate_pulses':circuits.pulse_counts(),
                   'examples':[resource_plan(),resource_plan(duration=F(1,16),error=F(1,20))],
                   'fixed_binary64_word_yields_arbitrary_accuracy':False,
                   'refinement_uses_explicit_dyadic_port_commands':True,
                   'global_ground_oracle_used':False})
