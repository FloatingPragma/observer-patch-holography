"""Literal CAR proton-pion channel projector; no fitted mass vector."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement,permutations,product
import sympy as sy
from model.openwave_baryon_dispersion import exact_native as base
from model.openwave_pole_weak_current.native import contraction as weak
from model.openwave_strange_shell_current.native import contraction,frame
from model.openwave_strange_composites.nonleptonic import word
from model.openwave_strange_composites.source import INDEX,native


def epsilon(shape):
    out=defaultdict(int)
    for colours in permutations(range(3)):
        sign=(-1)**sum(colours[i]>colours[j]for i in range(3)for j in range(i+1,3))
        ops=[(True,3*INDEX[x,0,s,w]+c)for (x,s,w),c in zip(shape,colours)]
        mask,sgn=word(0,ops)
        if sgn:out[mask]+=sign*sgn
    out={k:v for k,v in out.items()if v};n=sum(x*x for x in out.values())
    assert n>0
    return {k:sy.simplify(v/sy.sqrt(n))for k,v in out.items()}


def normalized_rational_column(col):
    col={k:sy.simplify(v)for k,v in col.items()if sy.simplify(v)!=0}
    if not col:return {}
    factor=col[min(col)]
    result={k:sy.simplify(v/factor)for k,v in col.items()}
    assert all(x.is_Rational for x in result.values())
    return result


@lru_cache(None)
def proton_fields():
    labels=tuple(combinations_with_replacement(range(16),3))
    shapes=sorted({base.canonical(lab)[1]for lab in labels if max(base.ORBITALS[i][0]for i in lab)-min(base.ORBITALS[i][0]for i in lab)<=1})
    proton=[s for s in shapes if sum(w==0 for x,spin,w in s)==2 and sum(spin==0 for x,spin,w in s)==2]
    B=weak()['native_proton_embedding'];assert B.shape==(12,5)
    columns=[]
    for j in range(5):
        out=defaultdict(lambda:sy.S.Zero)
        for i,shape in enumerate(proton):
            if B[i,j]:
                for mask,v in epsilon(shape).items():out[mask]+=B[i,j]*v
        columns.append(normalized_rational_column(out))
    # The first field is onsite. Its translate supplies the other local origin.
    assert all(mask<(1<<24)for mask in columns[0])
    columns.append({mask<<24:v for mask,v in columns[0].items()})
    return columns


def pion_field(x,y):
    join=native()[2][0];out={}
    for s,c,cc in product(range(2),range(3),range(3)):
        coefficient=(1-2*s)*int(round(join[c,cc].real))
        assert complex(int(round(join[c,cc].real)))==join[c,cc]
        if coefficient:
            q=3*INDEX[x,0,s,1]+c;a=48+3*INDEX[y,0,1-s,0]+cc
            mask,sign=word(0,[(True,q),(True,a)])
            assert sign
            out[mask]=coefficient*sign
    assert sum(v*v for v in out.values())==6
    return out


def product_column(left,right):
    out=defaultdict(lambda:sy.S.Zero)
    for a,u in left.items():
        ops=[(True,i)for i in range(48)if a>>i&1]
        for b,v in right.items():
            mask,sgn=word(b,ops)
            if sgn:out[mask]+=sgn*u*v
    return normalized_rational_column(out)


@lru_cache(None)
def channel():
    fields=proton_fields();columns=[];labels=[]
    for j,(x,y)in product(range(6),product(range(2),repeat=2)):
        if j==5 and x==y==1:continue
        col=product_column(fields[j],pion_field(x,y))
        assert col
        columns.append(col);labels.append((j,x,y))
    vertex=contraction();masks=sorted(set(vertex['raw_output_masks']).union(*(set(c)for c in columns)));index={m:i for i,m in enumerate(masks)}
    X=sy.zeros(len(masks),len(columns));V=sy.zeros(len(masks),56)
    for j,col in enumerate(columns):
        for mask,v in col.items():X[index[mask],j]=v
    for i,mask in enumerate(vertex['raw_output_masks']):V[index[mask],:]=vertex['matrix'][i,:]
    gram=X.T*X;_,pivots=gram.rref();Y=X[:,list(pivots)];G=Y.T*Y
    assert G.det()>0
    _,P=frame();B=Y.T*V*P
    projected=sy.simplify(B.T*G.inv()*B/6)
    total=vertex['seven_Gram_without_Cabibbo']
    assert projected==projected.T
    return {'basis':Y,'basis_Gram':G,'vertex':V,'overlap':B,'projected_Gram':projected,'total_Gram':total,
            'source_Gram':P.T*P,'all_columns_Gram':gram,'independent_columns':pivots,
            'composite_labels':labels,'raw_union_dimension':len(masks),'raw_masks':masks}


@lru_cache(None)
def separated_hopping():
    """Actual native kinetic component, without pretending it is all Hstrong."""
    fields=proton_fields()
    targets=[product_column(fields[0],pion_field(1,1)),
             product_column(fields[5],pion_field(0,0))]
    vertex=contraction();_,P=frame()
    # All 48 quark plus antiquark modes and both directions are retained.
    words=[]
    for family,spin,weak_axis,colour,species in product(range(2),range(2),range(2),range(3),range(2)):
        a=48*species+3*INDEX[0,family,spin,weak_axis]+colour
        b=48*species+3*INDEX[1,family,spin,weak_axis]+colour
        words.extend((((True,a),(False,b)),((True,b),(False,a))))
    assert len(words)==96
    target_norms=[sum(v*v for v in target.values())for target in targets]
    direct=sy.zeros(2,56);hopping=sy.zeros(2,56)
    for row,mask in enumerate(vertex['raw_output_masks']):
        v=vertex['matrix'][row,:]
        for j,target in enumerate(targets):
            if mask in target:direct[j,:]+=target[mask]*v
        for ops in words:
            final,sgn=word(mask,ops)
            if sgn:
                for j,target in enumerate(targets):
                    if final in target:hopping[j,:]-=sgn*target[final]*v
    direct=direct*P;amplitude=hopping*P
    assert direct==sy.zeros(2,7)
    G=sy.diag(*target_norms)
    gram=sy.simplify(amplitude.T*G.inv()*amplitude/6)
    assert gram!=sy.zeros(7)
    return {'amplitude':amplitude,'target_Gram':G,'current_then_hopping_Gram_without_t_Cabibbo':gram,
            'direct_amplitude':direct,'literal_hopping_words':len(words)}
