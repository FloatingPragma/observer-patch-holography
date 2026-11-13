"""Exact CAR current contraction on the complete compact triquark frame."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
import sympy as sy
from model.openwave_baryon_dispersion import exact_native as native


def canonical_basis(P,rank):
    cols=[]
    for j in range(P.cols):
        v=P[:,j]
        for q in cols:v=sy.simplify(v-q*(q.T*v)[0])
        n=sy.simplify((v.T*v)[0])
        if n==0:continue
        assert n>0
        v=sy.simplify(v/sy.sqrt(n))
        first=next(x for x in v if x!=0)
        if first<0:v=-v
        cols.append(v)
        if len(cols)==rank:break
    assert len(cols)==rank
    B=sy.Matrix.hstack(*cols)
    assert sy.simplify(B.T*B)==sy.eye(rank)
    return B


@lru_cache(None)
def contraction():
    A,G,(S,I)=native.complete()
    labels=tuple(combinations_with_replacement(range(16),3))
    selected=[lab for lab in labels if max(native.ORBITALS[i][0]for i in lab)-min(native.ORBITALS[i][0]for i in lab)<=1]
    shapes=sorted({native.canonical(lab)[1]for lab in selected})
    rows={native.canonical(lab)[1]:lab for lab in selected if native.canonical(lab)[0]==1}
    indices={shape:j for j,shape in enumerate(shapes)}
    neutron=[i for i,shape in enumerate(shapes)if sum(w==0 for x,s,w in shape)==1 and sum(s==0 for x,s,w in shape)==2]
    proton=[i for i,shape in enumerate(shapes)if sum(w==0 for x,s,w in shape)==2 and sum(s==0 for x,s,w in shape)==2]
    Dn=sy.diag(*[sy.sqrt(G[i,i])for i in neutron]);Dp=sy.diag(*[sy.sqrt(G[i,i])for i in proton])
    SN=Dn.inv()*S.extract(neutron,neutron)*Dn.inv();IN=Dn.inv()*I.extract(neutron,neutron)*Dn.inv()
    projector=sy.simplify((sy.Rational(15,4)*sy.eye(len(neutron))-SN)*(sy.Rational(15,4)*sy.eye(len(neutron))-IN)/9)
    assert projector==projector.T and sy.simplify(projector*projector)==projector
    cols=[]
    for count,rank in ((3,1),(2,2),(1,2)):
        group=[j for j,i in enumerate(neutron)if sum(x==0 for x,s,w in shapes[i])==count]
        C=canonical_basis(projector.extract(group,group),rank)
        B=sy.zeros(len(neutron),rank)
        for j,i in enumerate(group):B[i,:]=C[j,:]
        cols.append(B)
    Bn=sy.Matrix.hstack(*cols)
    assert sy.simplify(Bn.T*Bn)==sy.eye(5)
    currents=[]
    for axial in(False,True):
        C=sy.zeros(len(proton),len(neutron));lookup={i:j for j,i in enumerate(proton)}
        for j,i in enumerate(neutron):
            label=rows[shapes[i]]
            for x in(1,2):
                for spin in(0,1):
                    coefficient=1-2*spin if axial else 1
                    for target,value in native.one_body(label,native.INDEX[x,spin,0],native.INDEX[x,spin,1]).items():
                        k=indices[native.canonical(target)[1]]
                        assert k in lookup
                        C[lookup[k],j]+=coefficient*value*G[k,k]
        currents.append(sy.simplify(Dp.inv()*C*Dn.inv()))
    V,Acurrent=currents;Bp=sy.simplify(V*Bn)
    assert sy.simplify(Bp.T*Bp)==sy.eye(5)
    vector=sy.simplify(Bp.T*V*Bn);axial=sy.simplify(Bp.T*Acurrent*Bn)
    from model.openwave_baryon_dispersion.source import matrix
    H5,G5=matrix(F(1),F(0))
    d5=sy.diag(*[sy.sqrt(x)for x in G5.diagonal()])
    hn=Dn.inv()*A.extract(neutron,neutron)*Dn.inv()
    assert sy.simplify(Bn.T*hn*Bn-d5.inv()*H5*d5.inv())==sy.zeros(5)
    expected=sy.diag(sy.Rational(5,3),1,1,1,1)
    for i,j in((1,2),(3,4)):expected[i,j]=expected[j,i]=sy.Rational(2,3)
    assert vector==sy.eye(5) and axial==expected
    M=sy.simplify(Acurrent.T*Acurrent)
    assert sy.simplify((M-sy.eye(M.rows))*(M-4*sy.eye(M.rows)))==sy.zeros(M.rows)
    return {'axial_full_fixed_weight_square':M,'native_neutron_embedding':Bn,'native_proton_embedding':Bp,
        'vector_numerator':G5,'axial_numerator':sy.simplify(d5*axial*d5),
        'complete_native_compact_shapes':len(shapes),'neutron_fixed_weight_shapes':len(neutron),
        'proton_fixed_weight_shapes':len(proton),'actual_complete_native_contraction':True,
        'vector_isospin_half_norm_exact':1,'axial_full_three_quark_norm_upper':3}
