"""Finite displacement bounds for one native unit in its stationary sector."""
from fractions import Fraction as F
from . import geometry as g


def _count(value):
    if type(value) is not int or value<0:raise ValueError('a nonnegative exact attempt count required')
    return value


def displacement_bound(level,attempts,radius_squared=F(1)):
    level=g.level(level);attempts=_count(attempts)
    if type(radius_squared) is not F or radius_squared<=0:
        raise ValueError('a positive exact squared comparison distance required')
    mean=F(attempts,12600*16**level)
    return {'mean_squared_displacement_upper':mean,
            'displacement_event_probability_upper':min(F(1),mean/radius_squared),
            'radius_squared':radius_squared}


def ballistic_bound(level,time=F(1),radius_squared=F(1)):
    level=g.level(level)
    if type(time) is not F or time<0:raise ValueError('a nonnegative rational count-scale parameter required')
    attempts=(1260*8**level*time.numerator)//time.denominator
    return {'attempts':attempts,**displacement_bound(level,attempts,radius_squared),
            'simpler_event_upper':min(F(1),time/(10*2**level*radius_squared))}


def exact_moment(level,attempts):
    """Uniform one-unit endpoint moment from integer port-vector recurrences.

    No configuration carrier, Markov logarithm, Gaussian limit or numerical
    diagonalization is used. This is a kernel calculation, not host execution.
    """
    attempts=_count(attempts)
    if attempts>8:raise ValueError('this exact finite control caps eight iterations')
    graph=g.graph(level);rows,denominator=g.coordinates(level);factor=2*len(graph.seams)
    columns=[list(row) for row in rows]
    for _ in range(attempts):
        following=[[factor*x for x in row] for row in columns]
        for a,b in graph.seams:
            for j in range(6):
                delta=columns[a][j]-columns[b][j]
                following[a][j]-=delta;following[b][j]+=delta
        columns=following
    numerator=(0,0);power=factor**attempts
    for left,right in zip(rows,columns,strict=True):
        for j in (0,2,4):
            delta=(power*left[j]-right[j],power*left[j+1]-right[j+1])
            numerator=g.add(numerator,g.multiply((left[j],left[j+1]),delta))
    scale=F(2,graph.ports*power)
    result=g.divide((scale*numerator[0],scale*numerator[1]),
                    (F(5*denominator**2,2),F(denominator**2,2)))
    bound=displacement_bound(level,attempts)['mean_squared_displacement_upper']
    if g.sign(result)<0 or not g.at_most(result,bound):
        raise ArithmeticError('the exact native moment violates the proven bound')
    return result


def clock_counts(level,attempts,degree):
    level=g.level(level);attempts=_count(attempts)
    if type(degree) is not int or degree not in (5,6):raise ValueError('an actual native port degree required')
    M=g.census(level)['seams']
    return {'receiver_incident_mean':F(attempts*degree,M),
            'cell_internal_mean':F(attempts*30,M),
            'common_global_attempts':attempts}


def incident_bound(level,ticks,incident_seams,radius_squared=F(1)):
    """Stationary stopped-clock bound; not substitution of a mean stopping time."""
    level=g.level(level);ticks=_count(ticks)
    M=g.census(level)['seams']
    if type(incident_seams) is not int or not 1<=incident_seams<=M:
        raise ValueError('a nonempty literal fixed tick seam set required')
    if type(radius_squared) is not F or radius_squared<=0:
        raise ValueError('a positive exact squared comparison distance required')
    mean=F(ticks*ticks,10*incident_seams*4**level)
    return {'mean_squared_displacement_upper':mean,
            'displacement_event_probability_upper':min(F(1),mean/radius_squared),
            'expected_elapsed_global_attempts':F(ticks*M,incident_seams),
            'tick_probability_per_global_attempt':F(incident_seams,M),
            'mean_stopping_time_substituted_into_deterministic_bound':False}


def timeout(level,ticks,incident_seams,cap):
    from math import comb
    level=g.level(level);ticks=_count(ticks);cap=_count(cap)
    if cap>256:raise ValueError('exact finite timeout evaluation capped at256')
    a=incident_bound(level,ticks,incident_seams)['tick_probability_per_global_attempt']
    return sum((F(comb(cap,j))*a**j*(1-a)**(cap-j) for j in range(min(ticks,cap+1))),F())


def certificate():
    values=[{'level':ell,'attempts':t,'exact_squared_displacement':g.display(exact_moment(ell,t)),
             'upper':str(displacement_bound(ell,t)['mean_squared_displacement_upper'])}
            for ell in (0,1,2) for t in (0,1,2,4)]
    examples=[]
    for ell in (0,1,2,20,200):
        item=ballistic_bound(ell)
        examples.append({'level':ell,**{k:str(v) if type(v) is F else v for k,v in item.items()}})
    return {'finite_exact_moments':values,'ballistic_count_bounds':examples,
            'stationary_internal_clock_bounds':[
                {'level':ell,'ticks':k,**{key:str(v) if type(v) is F else v for key,v in
                   incident_bound(ell,k,30).items()}}
                for ell in (0,1,2,200) for k in (1,2)],
            'source_step':'P=I-L/(2M)',
            'native_record_displacement_bound':'T/(12600*16^level)',
            'stationary_uniform_one_unit_preparation_required':True,
            'Poisson_clock_Gaussian_limit_or_physical_position_assumed':False,
            'new_diffusive_continuum_theorem_claimed':False}
