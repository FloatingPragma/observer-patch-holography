"""Exact radical-field certificates for the finite source phase compiler.

All comparisons are rational. Decimal or binary floating square roots cannot
decide the endpoint, preboost branch or accepted rounding cell in this audit.
"""
from fractions import Fraction as F


def require(ok,message):
    if not ok:raise AssertionError(message)


class Radical:
    """a+b sqrt(d), with one nonnegative rational d shared by a calculation."""
    def __init__(self,a,b=0,d=0):self.a,self.b,self.d=F(a),F(b),F(d)
    def other(self,x):
        if isinstance(x,Radical):
            require(x.d==self.d,'Mixed radical fields');return x
        return Radical(x,0,self.d)
    def __add__(self,x):
        x=self.other(x);return Radical(self.a+x.a,self.b+x.b,self.d)
    __radd__=__add__
    def __neg__(self):return Radical(-self.a,-self.b,self.d)
    def __sub__(self,x):return self+-self.other(x)
    def __rsub__(self,x):return self.other(x)+-self
    def __mul__(self,x):
        x=self.other(x)
        return Radical(self.a*x.a+self.b*x.b*self.d,self.a*x.b+self.b*x.a,self.d)
    __rmul__=__mul__
    def __pow__(self,n):
        require(type(n)is int and n>=0,'Unsupported radical power')
        result=Radical(1,0,self.d)
        for _ in range(n):result=result*self
        return result
    def sign(self):
        a,b,d=self.a,self.b,self.d
        require(d>=0,'Negative discriminant')
        sign=lambda x:(x>0)-(x<0)
        if not b or not d:return sign(a)
        if not a:return sign(b)
        if sign(a)==sign(b):return sign(a)
        difference=a*a-b*b*d
        return sign(difference)*sign(a)
    def serialized(self):return [str(self.a),str(self.b),str(self.d)]


def phase_coordinates(record):
    re,im=F(record['real']),F(record['imag'])
    require(re*re+im*im==1,'Phase is not exactly unitary')
    return re,im


def certify_source(counts,record,preboost):
    q=[[F(x)for x in row]for row in counts]
    require(len(q)==3 and all(len(row)==3 for row in q),'Not a three-position source block')
    constant=sum(q[0])
    require(all(sum(row)==constant for row in q),'Uniform sector is not invariant')
    require(all(q[i][j]==q[j][i]for i in range(3)for j in range(3)),'Nonsymmetric source block')
    trace=sum(q[i][i]for i in range(3))-constant
    square=sum(q[i][j]*q[j][i]for i in range(3)for j in range(3))-constant**2
    disc=2*square-trace**2
    require(disc>=0,'Invalid optical characteristic discriminant')
    s=Radical(abs(trace)/2,F(1,2),disc)
    require((s-F(1,5)).sign()>=0 and (s-1).sign()<=0,'Outside declared one-boost norm domain')
    expected=(s-F(1,2)).sign()<0
    require(type(preboost)is bool and preboost==expected,'Incorrect exact preboost branch')
    boosted=3*s-4*s**3 if expected else s
    h=4*boosted**2-1
    require(h.sign()>=0,'Boosted signal lies below half')
    re,im=phase_coordinates(record)
    if h.sign()==0:
        require(re==-1 and im==0,'Exact half endpoint does not use minus one')
        return {'endpoint':True,'phase_chord_error_bound':'0','exact_h':h.serialized()}
    require('numerator'in record and record.get('endpoint')is None,'False half endpoint')
    n,bits=record['numerator'],record['bits']
    require(type(n)is int and type(bits)is int and bits==40,'Undeclared rounding grid')
    t=F(n,2**bits)
    require(re==(1-t*t)/(1+t*t) and im==2*t/(1+t*t),'Phase does not encode the stated tangent')
    lower=F(2*n-1,2**(bits+1));upper=F(2*n+1,2**(bits+1))
    require(lower>0,'Unexpected tangent cell at zero')
    # t_ideal=1/sqrt(h), so this is the exact reciprocal-square cell.
    left=(h-1/(upper*upper)).sign();right=(h-1/(lower*lower)).sign()
    require(left>=0 and right<=0,'Tangent lies outside claimed nearest rounding cell')
    require((left!=0 and right!=0)or n%2==0,'Half-even tie has odd numerator')
    return {'endpoint':False,'phase_chord_error_bound':str(F(1,2**bits)),
            'exact_h':h.serialized(),'rounding_cell':[str(lower),str(upper)],
            'rounding_inequality_signs':[left,right]}


def certify_five_selector(record):
    re,im=phase_coordinates(record);n,bits=record['numerator'],record['bits']
    require(type(n)is int and bits==40,'Undeclared five-selector grid')
    t=F(n,2**bits);lo=F(2*n-1,2**(bits+1));hi=F(2*n+1,2**(bits+1))
    require(lo>=0 and lo*lo<=F(2,3)<=hi*hi,'Five-selector tangent cell misses sqrt(2/3)')
    require(re==(1-t*t)/(1+t*t) and im==2*t/(1+t*t),'Five-selector phase encoding')
    return {'phase_chord_error_bound':str(F(1,2**bits)),'rounding_cell':[str(lo),str(hi)],
            'exact_target_tangent_square':'2/3'}
