"""Exact real quadratic field Q(sqrt(5))."""
from dataclasses import dataclass
from fractions import Fraction as F
@dataclass(frozen=True)
class Q:
 a:F=F(0)
 b:F=F(0)
 def __post_init__(self):object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
 def __add__(self,x):
  x=x if isinstance(x,Q)else Q(x);return Q(self.a+x.a,self.b+x.b)
 __radd__=__add__
 def __neg__(self):return Q(-self.a,-self.b)
 def __sub__(self,x):return self+-asq(x)
 def __rsub__(self,x):return asq(x)+-self
 def __mul__(self,x):
  x=asq(x);return Q(self.a*x.a+5*self.b*x.b,self.a*x.b+self.b*x.a)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=asq(x);d=x.a*x.a-5*x.b*x.b
  if not d:raise ZeroDivisionError('zero quadratic-field element')
  return self*Q(x.a/d,-x.b/d)
 def __pow__(self,n):
  if type(n)is not int or n<0:raise ValueError('nonnegative integral exponent')
  out=Q(1)
  for _ in range(n):out=out*self
  return out
 def pair(self):return[str(self.a),str(self.b)]
def asq(x):return x if isinstance(x,Q)else Q(x)
