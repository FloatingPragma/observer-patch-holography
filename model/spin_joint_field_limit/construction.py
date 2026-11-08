"""Complete finite positive-field Pauli histories with a source and clock budget."""
from fractions import Fraction as F
from math import ceil

from model.native_spin_observer import circuits
from . import bounds, clock_comparison


def history_cutoff(radius, support_radius, support_sites, coupling, duration, interventions):
    if (type(radius) is not int or type(support_radius) is not int
        or not 0<=support_radius<radius or type(support_sites) is not int
        or not 1<=support_sites<=(2*support_radius+1)**3
        or type(interventions) is not int or interventions<0):
        raise ValueError('valid nested local-history support required')
    j=bounds.rational(coupling,'coupling',True)
    t=bounds.rational(duration,'duration')
    if not t:return F()
    exponent=radius-support_radius-ceil(96*j*t)
    if exponent<0:return F(2)
    prefactor=(2*interventions+1)*6*j*support_sites*(2*radius+1)**2*t
    return min(F(2),prefactor/F(1<<exponent))


def probe_history(n, probes, scaled_durations, accuracy, coupling=F(1)):
    """Specified local encoded-Pauli probes; every pointer/failure outcome retained."""
    family=bounds.family(n,coupling)
    eps=bounds.rational(accuracy,'accuracy',True)
    if eps>=1:raise ValueError('accuracy must be less than one')
    probes=tuple(probes); durations=tuple(scaled_durations)
    if not probes or len(durations)!=len(probes)+1:
        raise ValueError('at least one probe and one more clock slot required')
    if any(type(t) not in (int,F) for t in durations):
        raise ValueError('exact rational scaled durations required')
    for item in probes:
        if (type(item) is not tuple or len(item)!=2 or item[0] not in ('X','Y','Z')
            or type(item[1]) is not tuple or len(item[1])!=3
            or any(type(x) is not int for x in item[1])):
            raise ValueError('a Pauli axis and integer source site are required')
    support_radius=max(abs(x) for _,site in probes for x in site)
    support_sites=len({site for _,site in probes})
    time=sum((abs(F(t)) for t in durations),F())/(family['coupling']*family['p'])
    offset=max(1,ceil(96*family['coupling']*time)+1)
    while True:
        radius=support_radius+offset
        cutoff=history_cutoff(radius,support_radius,support_sites,family['coupling'],time,len(probes))
        if cutoff<=eps/4:break
        offset*=2
    # Both the target infinite clock and the finite source clock are compared
    # with this inner clock. Two unhalved trace bounds cost one cutoff in half trace.
    prepared=bounds.local_preparation(n,radius,(2*radius+1)**3,eps/4,coupling)
    source_radius=prepared['finite_source']['radius']
    clock=clock_comparison.cube_clock_plan(source_radius,n,durations,eps/4,coupling)
    pulses=sum(len(circuits.probe(axis)) for axis,_ in probes)
    bits=1
    while 29*pulses*F(1,1<<bits)>eps/4:bits+=1
    instrument_error=29*pulses*F(1,1<<bits)
    total=(cutoff+prepared['prepared_marginal_half_trace_upper']
           +clock['complete_clock_half_trace_error_upper']+instrument_error)
    if total>eps:raise ArithmeticError('complete source history error exceeds target')
    return {**family,'probes':[{'axis':a,'site':list(x)} for a,x in probes],
            'probe_observable':'encoded local Pauli; physical-spin signs use the exact source frame',
            'scaled_durations':list(map(F,durations)),
            'comparison_inner_radius':radius,'comparison_support_sites':support_sites,
            'absolute_physical_duration':time,
            'two_clock_cutoffs_combined_half_trace_upper':cutoff,
            'prepared_source':prepared,'retained_field_clock':clock,
            'probe_native_pulses':pulses,'probe_dyadic_port_bits':bits,
            'probe_complete_isometry_error_upper':instrument_error,
            'complete_output_half_trace_upper':total,'requested_accuracy':eps,
            'ideal_pointer_record_count':1<<len(probes),
            'all_native_failure_and_pointer_records_retained':True,
            'fresh_pointer_registers':len(probes),'large_history_executed':False,
            'initially_independent_reference_allowed':True,
            'target':'the specified complete local H_B history in omega_B',
            'effective_continuum_history_accuracy_included':False,
            'coherent_time_control_included':False}
