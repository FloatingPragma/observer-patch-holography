from fractions import Fraction as F
ZERO=(F(0),F(0))
def q(a=0,b=0): return F(a),F(b)

def add(x,y): return x[0]+y[0],x[1]+y[1]

def neg(x): return -x[0],-x[1]

def sub(x,y): return add(x,neg(y))

def mul(x,y):
    a,b=x; c,d=y
    return a*c+b*d,a*d+b*c+b*d

def div(x,y):
    a,b=y; denominator=a*a+a*b-b*b
    if not denominator: raise ZeroDivisionError
    return mul(x,((a+b)/denominator,-b/denominator))

def total(xs):
    out=ZERO
    for x in xs: out=add(out,x)
    return out

def dot(x,y): return total(mul(a,b) for a,b in zip(x,y))

def geometry(seams):
    z=ZERO
    points=[(q(-1),q(0,1),z),(q(1),q(0,1),z),(q(-1),q(0,-1),z),(q(1),q(0,-1),z),
            (z,q(-1),q(0,1)),(z,q(1),q(0,1)),(z,q(-1),q(0,-1)),(z,q(1),q(0,-1)),
            (q(0,1),z,q(-1)),(q(0,1),z,q(1)),(q(0,-1),z,q(-1)),(q(0,-1),z,q(1))]
    actual=tuple((a,b) for a in range(12) for b in range(a+1,12) if dot(points[a],points[b])==q(0,1))
    assert actual==seams, 'geometry must use actual row1 seam labels'
    antipode=[points.index(tuple(neg(v) for v in p)) for p in points]
    rotated=[tuple(mul(q(s),v) for s,v in zip((-1,-1,1),p)) for p in points]
    permutation=[points.index(p) for p in rotated]
    assert permutation[0]==antipode[0]==3
    assert all(permutation[permutation[p]]==p for p in range(12))
    assert all(tuple(sorted((permutation[a],permutation[b]))) in seams for a,b in seams)
    # Same explicit port vectors as the Lean finite current quotient.
    lean=[(z,q(-1),q(0,-1)),(q(-1),q(0,-1),z),(q(0,-1),z,q(-1)),(q(1),q(0,-1),z),
          (z,q(1),q(0,-1)),(q(0,-1),z,q(1)),(q(0,1),z,q(-1)),(z,q(-1),q(0,1)),
          (q(-1),q(0,1),z),(q(0,1),z,q(1)),(q(1),q(0,1),z),(z,q(1),q(0,1))]
    return points,antipode,permutation,[lean.index(p) for p in points]
