"""Exact finite collar bounds for the native regional modular construction."""
from fractions import Fraction as F


def integer(x,low=0):
 if type(x)is not int or x<low:raise ValueError('exact integer in range required')
 return x

def tail_upper(x,r):
 """Upper bound for sum_(n>=r) x^n/n!, with a rational geometric tail."""
 if type(x)is not F or x<0:raise ValueError('nonnegative exact rational required')
 integer(r)
 if r>4096 or x>256:raise ValueError('finite rational tail arithmetic budget exceeded')
 N=max(r,int(x)+1)+32;term=F(1);total=F(1)if r==0 else F()
 for n in range(1,N+1):
  term*=x/n
  if n>=r:total+=term
 first=term*x/(N+1);ratio=x/(N+2)
 if ratio>=1:raise ArithmeticError('the positive tail ratio is not contractive')
 return total+first/(1-ratio)

def locality_bound(s,degree,crossings,r,time,p):
 for x in(s,degree):integer(x,1)
 for x in(crossings,r):integer(x)
 if type(time)is not F or type(p)is not F or not 0<=p<=1:raise ValueError('exact time and density required')
 J=1/(1-p*(1-p));tail=tail_upper(4*degree*J*abs(time),r)
 return min(F(2),2*s*J*crossings*(1+abs(1-2*p))*abs(time)*tail)

def uniform_modular_error(V,M,q,v,internal,crossings,p):
 for x in(V,M,v):integer(x,1)
 for x in(q,internal,crossings):integer(x)
 if not(q<=V and v<V-1 and internal+crossings<=M)or type(p)is not F or not 0<=p<=1:raise ValueError('valid exact source counts and density required')
 b=V-v;L=2*internal+crossings
 if L>=M:raise ValueError('uniform logarithm gate L<M required')
 dr=F(V,b)*abs(F(q,V)-p)+F(v,b)
 da=dr+F(internal+3*crossings,4*M)+F(1,4*(b-1))
 return F(L*L,M)/(1-F(L,M))+4*crossings*dr+F(8,3)*L*da

def combined_bound(V,M,q,v,internal,crossings,s,degree,r,time,p):
 modular=uniform_modular_error(V,M,q,v,internal,crossings,p)
 collar=locality_bound(s,degree,crossings,r,time,p)
 return {'regional_generator_error_upper':modular,'collar_flow_error_upper':collar,
         'regional_modular_to_larger_exchange_flow_upper':min(F(2),2*abs(time)*modular+collar)}
