"""Exact native-algebra helpers for the independent interacting join."""
from fractions import Fraction
import math
import sympy as s
from .native_inputs import exterior_derivative, assert_zero

PIN='3cf3c14810671a16e61b7915f69079a1773a12a2'

def add(vector,key,value):
    if not value:return
    value=value+vector.get(key,0)
    if value:vector[key]=value
    else:vector.pop(key,None)

def rational(x):
    q=s.cancel(x)
    if not q.is_Rational:raise ArithmeticError('nonrational exact coefficient')
    return Fraction(int(q.p),int(q.q))

def norm(vector):return sum((x*x for x in vector.values()),Fraction(0))

def dot(a,b):return sum((v*b.get(k,0) for k,v in a.items()),Fraction(0))

def primitive(vector):
    if not vector:raise ArithmeticError('zero kernel vector')
    vector={k:Fraction(v) for k,v in vector.items() if v}
    denominator=math.lcm(*(x.denominator for x in vector.values()))
    ints={k:x.numerator*(denominator//x.denominator) for k,x in vector.items()}
    divisor=math.gcd(*ints.values())
    if ints[min(ints)]<0:divisor=-divisor
    return {k:Fraction(v//divisor) for k,v in ints.items()}

def standard_generators(data):
    basis=[]
    for offset,width in ((0,3),(3,2)):
        for a in range(offset,offset+width):
            for b in range(a+1,offset+width):
                g=s.zeros(5);g[a,b]=1;g[b,a]=-1;basis.append(g)
                g=s.zeros(5);g[a,b]=s.I;g[b,a]=s.I;basis.append(g)
        if width==3:basis.extend((s.I*s.diag(1,-1,0,0,0),s.I*s.diag(1,1,-2,0,0)))
        else:basis.append(s.I*s.diag(0,0,0,1,-1))
    basis.append(s.I*s.diag(-2,-2,-2,3,3))
    anti=[];real=[];scalar=[]
    for g in basis:
        lift=exterior_derivative(g,data['subsets'])
        phase=1 if all(s.im(x)==0 for x in g) else s.I
        anti.append(lift)
        real.append((lift/phase).applyfunc(lambda x:s.Rational(x)))
        scalar.append((-g[3:,3:].T/phase).applyfunc(lambda x:s.Rational(x)))
    assert_zero(anti[-1]/s.I-data['operators']['q'],'source central charge')
    c3=-sum((w*g*g for w,g in zip([s.Rational(1,4)]*7+[s.Rational(1,12)],anti[:8])),s.zeros(15))
    c2=-sum((g*g/4 for g in anti[8:11]),s.zeros(15))
    assert_zero(c3-data['operators']['color_C2'],'source color metric')
    assert_zero(c2-data['operators']['weak_C2'],'source weak metric')
    # Matching Casimirs alone does not certify the same embedded current
    # algebra.  Verify equality of the full exact 12-dimensional spans.
    joined=s.Matrix.hstack(*(g.reshape(225,1) for g in (*data['currents'],*anti)))
    if joined.rank()!=12:raise ArithmeticError('source current algebra span')
    return tuple(anti),tuple(real),tuple(scalar)

def rank_mod_prime(rows,columns,target,prime=1000003):
    pivots={};used=0
    for source in rows:
        if len(pivots)==target:break
        row={int(j):int(v)%prime for j,v in source.items() if int(v)%prime};used+=1
        while row:
            p=min(row)
            if p not in pivots:
                inv=pow(row[p],-1,prime)
                pivots[p]={j:v*inv%prime for j,v in row.items()};break
            amount=row[p]
            for j,v in pivots[p].items():
                x=(row.get(j,0)-amount*v)%prime
                if x:row[j]=x
                else:row.pop(j,None)
    if len(pivots)!=target:raise ArithmeticError('modular lower rank certificate failed')
    return {'prime':prime,'columns':columns,'rank_lower_bound':len(pivots),'rows_consumed':used}

def create(mask,mode):
    bit=1<<mode
    if mask&bit:return None
    return mask|bit,(-1)**((mask&(bit-1)).bit_count())

def annihilate(mask,mode):
    bit=1<<mode
    if not mask&bit:return None
    return mask^bit,(-1)**((mask&(bit-1)).bit_count())

def occupied(mask):
    while mask:
        low=mask&-mask;yield low.bit_length()-1;mask^=low

def mode(site,internal,family,spin,families):
    return ((site*15+internal)*families+family)*2+spin

def decode(m,families):
    u,spin=divmod(m,2);u,family=divmod(u,families);site,internal=divmod(u,15)
    return site,internal,family,spin

def links(blocks):
    records=[None];reverse={};by_block={-1:[0]}
    for r,b in enumerate(blocks):
        by_block[r]=[]
        for i in b['indices']:
            for j in b['indices']:
                reverse[r,i,j]=len(records);by_block[r].append(len(records));records.append((r,i,j))
    if len(records)!=60:raise ArithmeticError('native link size')
    return records,reverse,by_block
