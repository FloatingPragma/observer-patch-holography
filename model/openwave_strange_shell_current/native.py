"""Exact native epsilon and weak CAR words on all compact uds shapes."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations,product
import sympy as sy
from model.openwave_cluster_bandwidth import lambda_source as lam
from model.openwave_strange_composites.source import INDEX,native
from model.openwave_strange_composites.nonleptonic import word


@lru_cache(None)
def frame():
    d=lam.source(4);labels=d['labels'];positions=d['positions']
    def shape(row):
        origin=min(d['orbs'][j][0]for j in row)
        return tuple((d['orbs'][j][0]-origin,d['orbs'][j][1])for j in row)
    shapes=sorted({shape(row)for row,pos in zip(labels,positions)if max(pos)-min(pos)<=1})
    ids=[next(i for i,row in enumerate(labels)if shape(row)==s and min(positions[i])==1)for s in shapes]
    S=lam.rational(d['spin'](tuple(range(4)))[ids][:,ids].toarray())
    I=lam.rational(d['iso'](tuple(range(4)))[ids][:,ids].toarray())
    Sz=sy.diag(*[sy.Rational(3-2*sum(s for x,s in sh),2)for sh in shapes])
    selector=(S-sy.Rational(3,4)*sy.eye(56))**2+I**2+(Sz-sy.Rational(1,2)*sy.eye(56))**2
    P=sy.Matrix.hstack(*selector.nullspace());assert P.shape==(56,7)
    assert P.T*P==lam.compact()[1]
    return shapes,P


def epsilon_column(shape):
    out=defaultdict(int)
    for cols in permutations(range(3)):
        sign=(-1)**sum(cols[i]>cols[j]for i in range(3)for j in range(i+1,3))
        ops=[]
        for f,((x,s),c)in enumerate(zip(shape,cols)):
            family,w=((0,0),(0,1),(1,1))[f]
            ops.append((True,3*INDEX[x,family,s,w]+c))
        mask,sgn=word(0,ops);assert sgn
        out[mask]+=sign*sgn
    assert sum(x*x for x in out.values())==6
    return dict(out)


@lru_cache(None)
def vertex_words():
    pauli=[sy.Matrix([[0,1],[1,0]]),sy.Matrix([[0,-sy.I],[sy.I,0]]),sy.diag(1,-1)]
    bar=[sy.eye(2),*(-x for x in pauli)];eta=(1,-1,-1,-1);epsilon=sy.Matrix([[0,1],[-1,0]])
    join=native()[2][0]
    exact_join=sy.Matrix([[int(round(z.real)) for z in row]for row in join])
    assert all(complex(exact_join[i,j])==join[i,j]for i in range(3)for j in range(3))
    out=[]
    for x,su,ss,sd,sa,c,d,dd in product(range(2),range(2),range(2),range(2),range(2),range(3),range(3),range(3)):
        spin=sum(eta[mu]*bar[mu][su,ss]*(bar[mu]*epsilon)[sd,sa]for mu in range(4))/2
        assert spin in(0,1,-1)
        coefficient=spin*exact_join[d,dd]
        if not coefficient:continue
        up=3*INDEX[x,0,su,0]+c;strange=3*INDEX[x,1,ss,1]+c
        down=3*INDEX[x,0,sd,1]+d;antiup=48+3*INDEX[x,0,sa,0]+dd
        out.append((int(coefficient),((True,up),(False,strange),(True,down),(True,antiup))))
    return out


@lru_cache(None)
def contraction():
    shapes,P=frame();cols=[]
    for shape in shapes:
        out=defaultdict(int)
        for coefficient,ops in vertex_words():
            for mask,amplitude in epsilon_column(shape).items():
                final,sign=word(mask,ops)
                if sign:out[final]+=coefficient*sign*amplitude
        cols.append({k:v for k,v in out.items()if v})
    masks=sorted(set().union(*(set(col)for col in cols)));index={m:i for i,m in enumerate(masks)}
    V=sy.zeros(len(masks),56)
    for j,col in enumerate(cols):
        for mask,value in col.items():V[index[mask],j]=value
    for mask in masks:
        q=[i for i in range(48)if mask>>i&1];a=[i-48 for i in range(48,96)if mask>>i&1]
        assert len(q)==4 and len(a)==1
        # 16 orbitals have ordered (x,family,spin,weak) binary coordinates.
        assert sum((i//3)%2==0 for i in q)==2
        assert all(((i//3)//4)%2==0 for i in q+a)
    s=F('0.22431');cab2=s*s*(1-s*s)
    M=V.T*V/6;compact=P.T*M*P
    bound=max(sum(abs(x)for x in M.row(i))for i in range(56))
    return {'matrix':V,'raw_output_masks':masks,'Gram_without_Cabibbo':M,'seven_Gram_without_Cabibbo':compact,
        'Gram_norm_squared_upper_without_Cabibbo':F(bound),'Cabibbo_squared':cab2,
        'compact_shape_count':56,'raw_N5_output_count':len(masks),'literal_CAR_word_count':len(vertex_words())}
