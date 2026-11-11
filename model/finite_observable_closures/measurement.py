"""Exact arithmetic enclosure of freshly executed interacting histories."""
from fractions import Fraction as F
from functools import lru_cache
import math

from model.common_source import experiment, finite_signal


def sqrt_upper(value, bits=100):
    if value < 0:
        raise ValueError('nonnegative radicand required')
    scale = 1 << bits
    numerator = math.isqrt(value.numerator*scale*scale//value.denominator)
    return F(numerator if numerator*numerator*value.denominator == value.numerator*scale*scale else numerator+1,scale)


def outward_interval(lower,upper=None):
    if upper is None: upper=lower
    grid=10**12
    return [str(F((lower*grid).numerator//(lower*grid).denominator,grid)),
            str(F(-(-(upper*grid).numerator//(upper*grid).denominator),grid))]


def density_distance_upper(binary, polynomial, remainder):
    """Half trace norm: exact rank-two Frobenius bound plus Taylor error."""
    actual = tuple((F(float(real)),F(float(imag))) for real,imag in binary)
    def outer(v,i,j):
        a,b=v[i]; c,d=v[j]
        return a*c+b*d,b*c-a*d
    square = F(0)
    for i in range(len(actual)):
        for j in range(len(actual)):
            a,b=outer(actual,i,j); c,d=outer(polynomial,i,j)
            square += (a-c)**2+(b-d)**2
    # Each exact ideal branch norm is <=1. Taylor branch norm is <=1+eta.
    error = sqrt_upper(square/2)+remainder*(2+remainder)/2
    return error, sum(real*real+imag*imag for real,imag in actual)


@lru_cache(maxsize=1)
def calculate():
    certificate = finite_signal.certificate(degree=64,steps=512)
    h = finite_signal.hamiltonian_matrix()
    if any(h[i][j] != (h[j][i][0],-h[j][i][1]) for i in range(16) for j in range(16)):
        raise ArithmeticError('actual finite Hamiltonian is not exactly selfadjoint')
    configurations=[]
    event_probabilities=[]
    complete_errors=[]
    for enabled in (True,False):
        execution=experiment.history(2,512,coupling=.25,feedback_enabled=enabled,current=0)
        polynomial,_,_,remainder=finite_signal.branch_vectors(enabled,64)
        expected=dict(polynomial)
        rows=[]
        error=mass=positive=F(0)
        for branch in execution['branches']:
            record=tuple(branch['record'])
            distance,probability=density_distance_upper(branch['unnormalized_terminal_state'],expected[record],remainder)
            error+=distance; mass+=probability
            if record[1] == 1: positive+=probability
            rows.append({'record':list(record),'executed_binary_probability_interval':outward_interval(probability),
                         'complete_branch_half_trace_error_upper':outward_interval(distance)[1]})
        if len(rows) != 4 or abs(mass-1)>F(1,10**9) or error>F(3,10**7):
            raise ArithmeticError('fresh finite source output escaped the exact enclosure')
        configurations.append({'feedback':enabled,'records':rows,
            'source_output_mass_interval':outward_interval(1-F(1,10**9),1+F(1,10**9)),
            'source_mass_defect_upper':'1/1000000000',
            'complete_output_half_trace_error_upper':'3/10000000',
            'recomputed_half_trace_error_upper_diagnostic':outward_interval(error)[1],
            'positive_scalar_executed_probability_interval':outward_interval(positive),
            'source_pulse_resources':execution['resources_per_realized_history']})
        # A computed complete output need not have exactly unit trace.
        # An event needs half trace norm plus half the total trace defect.
        event_probabilities.append(positive);complete_errors.append(F(3,10**7)+F(1,2*10**9))
    contrast=event_probabilities[0]-event_probabilities[1]
    event_error=sum(complete_errors,F(0))
    interval=(contrast-event_error,contrast+event_error)
    if interval[1]>=0 or not certificate['ideal_source_word_nonzero_response_certified']:
        raise ArithmeticError('the retained unconditional feedback contrast is not separated from zero')
    return {'dashboard_id':'measurement','recommended_status':'executed_theorem',
        'recommended_name':'Complete measurement and feedback in a finite interacting source',
        'source_dimension':16,'onsite_scalar_dimension':2,'source_steps_per_segment':512,
        'exact_continuous_Hermitian_generator':True,
        'exact_operator_certificate':certificate,
        'fresh_source_executions':configurations,
        'executed_contrast_interval':outward_interval(contrast),
        'continuous_contrast_from_fresh_execution_interval':outward_interval(*interval),
        'all_original_records_retained':True,'branch_renormalization':False,
        'host_output_error_certified_by_exact_binary_components':True,
        'each_event_probability_error_upper':'601/2000000000',
        'finite_premises_closed':['same16-state native source and rational Hamiltonian','both complete pointers and record-only feedback','fresh pulse execution with every terminal vector','exact source-output arithmetic comparison with independent Taylor target'],
        'source_scope':'fixed projected scalar vacuum and gauge blank, native current zero, exact ideal primitive services; the distinct infinite-occupation comparison keeps its own theorem and detector',
        'physical_hardware_error_or_continuum_claimed':False}
