"""Primitive quadratic projections and controlled collective-field moments."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial
from model.native_correlation_sectors import sectors as s


def profiles(first,second):
    if (type(first) is not tuple or type(second) is not tuple or len(first)!=len(second)
            or len(first)<4 or len(first)>256
            or any(type(x) not in (int,F) for x in first+second)
            or sum(first) or sum(second)):
        raise ValueError('two exact mean-zero profiles on four through256 ports required')
    return tuple(map(F,first)),tuple(map(F,second))


def projection(first,second):
    first,second=profiles(first,second);V=len(first)
    inner=sum(a*b for a,b in zip(first,second,strict=True))
    u=tuple((-2*a*b+inner/F(V-1))/F(V-2) for a,b in zip(first,second,strict=True))
    h={(a,b):first[a]*second[b]+first[b]*second[a]-u[a]-u[b] for a,b in combinations(range(V),2)}
    h={pair:value for pair,value in h.items() if value}
    if s.lower(h,V,2):raise ArithmeticError('projected quadratic field is not primitive')
    return h


def field(first,second,occupied,subset):
    first,second=profiles(first,second);V=len(first);s.degree_data(V,occupied)
    if (not 2<=occupied<=V-2 or type(subset) is not tuple or len(subset)!=occupied
            or tuple(sorted(set(subset)))!=subset or any(type(i) is not int or not 0<=i<V for i in subset)):
        raise ValueError('a literal configuration in a nonzero degree-two sector required')
    rho=F(occupied,V);chi=F(occupied*(V-occupied),V*(V-1))
    inner=sum(a*b for a,b in zip(first,second,strict=True))
    linear=sum(first[i]*second[i]*(int(i in subset)-rho) for i in range(V))
    a=F(V-2*occupied,V-2)
    return (sum(first[i] for i in subset)*sum(second[i] for i in subset)-chi*inner-a*linear)/(V*chi)


def correction_bound(first,second,occupied):
    first,second=profiles(first,second);V=len(first);s.norm_factor(V,2,occupied)
    chi=F(occupied*(V-occupied),V*(V-1));a=F(V-2*occupied,V-2)
    product=tuple(x*y for x,y in zip(first,second,strict=True));mean=sum(product)/V
    exact=a*a*sum((x-mean)**2 for x in product)/(V*V*chi)
    upper=a*a*max(abs(x) for x in first)**2*max(abs(x) for x in second)**2/(V*chi)
    if exact>upper:raise ArithmeticError('primitive correction exceeds uniform profile bound')
    return {'correction_mean_square':exact,'correction_mean_square_upper':upper,
            'chi':chi,'projection_coefficient':a}


def even_moment_bound(ports,occupied,supremum,order):
    s.degree_data(ports,occupied);s.integer(order,1)
    if (not 0<occupied<ports or type(supremum) is not F or supremum<0 or order>8):
        raise ValueError('nontrivial fixed number, rational profile bound and order1..8 required')
    chi=F(occupied*(ports-occupied),ports*(ports-1))
    sigma2=occupied*supremum**2/(ports*chi)
    return 2*(2*sigma2)**order*factorial(order)


def localized_profile(ports,weights):
    s.integer(ports,4)
    if (ports>256 or type(weights) is not dict or any(type(i) is not int or not 0<=i<ports
            or type(w) not in (int,F) for i,w in weights.items())):
        raise ValueError('bounded exact local profile dictionary required')
    mean=sum(weights.values(),F())/ports
    return tuple(F(weights.get(i,0))-mean for i in range(ports))
