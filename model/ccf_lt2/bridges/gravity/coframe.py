"""Actual CCF2 background variations; mechanical response, not Lorentz stress."""
from collections import defaultdict
from pathlib import Path
import sys
import sympy as sy
from ... import core, dirac, source, spatial


def edge_lift(poly,k,nedges):
    return {(w,tuple(sh[0] if j==k else 0 for j in range(nedges))):c for (w,sh),c in poly.items()}


def bases(src):
    if 'dirac' not in src['components']:
        raise ValueError('the declared default Dirac component is required')
    graph=src['graph'];directions,spacings=graph.edge_geometry();rows=[]
    for k,(edge,d,l) in enumerate(zip(graph.edges,directions,spacings)):
        derivative=[]
        for a in range(3):
            unit=tuple(int(i==a) for i in range(3))
            p=dirac.dirac_terms(graph.nvertices,[edge],[unit],[l],c=1,wilson_r=0)
            derivative.append(edge_lift(p,k,len(graph.edges)))
        both=dirac.dirac_terms(graph.nvertices,[edge],[(1,0,0)],[l],c=1,wilson_r=1)
        wilson=core.gadd(edge_lift(both,k,len(graph.edges)),core.scale(derivative[0],-1))
        u,v=edge;momentum={}
        for _,left,right,q in dirac.KINDS:
            shift=tuple(-q if j==k else 0 for j in range(len(graph.edges)))
            for family in range(3):
                x,y=dirac._spinor(u,left,right,family),dirac._spinor(v,left,right,family)
                p=dirac._dress(core.scale(dirac._bilinear(y,sy.eye(4),x),sy.I/(2*l)),shift)
                momentum=dirac._gadd(momentum,p,dirac._gadjoint(p))
        momentum=dirac.to_native(momentum,graph.nvertices)
        rows.append(dict(d=sy.Matrix(d),D=tuple(derivative),W=wilson,P=momentum))
    normals=[]
    for face in graph.faces:
        x,y,z=(sy.Matrix(graph.positions[i]) for i in face[:3])
        normal=(y-x).cross(z-x)
        norm=sy.sqrt(normal.dot(normal))
        if norm==0:raise ValueError('nondegenerate face plane required')
        normals.append(normal/norm)
    static=core.gadd(src['matter'],core.scale(src['components']['dirac'],-1))
    result=dict(edges=tuple(rows),normals=tuple(normals),static=static)
    validate_origin(src,result)
    return result


def validate_origin(src,basis):
    represented=core.gadd(*(core.scale(p,w) for b in basis['edges'] for p,w in zip(b['D'],b['d'])),
                          *(b['W'] for b in basis['edges']))
    if represented!=src['components'].get('dirac'):
        raise ValueError('this adapter requires the exact default c=r=1 Dirac component')


def weights(frame,edges,normals,shift=(0,0,0)):
    V=sy.Matrix(frame);beta=sy.Matrix(shift)
    if V.shape!=(3,3) or beta.shape!=(3,1):raise ValueError('three-dimensional coframe and shift required')
    if any(x.is_real is not True for x in (*V,*beta)):
        raise ValueError('real coframe entries and shift required')
    G=V.T*V;volume=V.det()
    if volume==0:raise ValueError('singular coframe')
    inverse=V.inv().T
    result=[]
    for edge in edges:
        d=edge['d'];kin=inverse*d
        result.append(dict(D=tuple(kin),W=sy.simplify(kin.dot(kin)),
                           E=sy.simplify((d.T*G*d)[0]/volume),P=sy.simplify(beta.dot(kin))))
    magnetic=tuple(sy.simplify((m.T*G*m)[0]/volume) for m in normals)
    return tuple(result),magnetic


def expression(src,basis,frame=sy.eye(3),shift=(0,0,0)):
    # Numeric coframes are physical members; symbolic frames are formal probes near I.
    V=sy.Matrix(frame)
    if not V.free_symbols and V.det().is_positive is not True:
        raise ValueError('positive determinant required')
    coeff,magnetic=weights(V,basis['edges'],basis['normals'],shift)
    terms=[basis['static']]
    for b,c in zip(basis['edges'],coeff):
        terms.extend(core.scale(p,w) for p,w in zip(b['D'],c['D']))
        terms.extend((core.scale(b['W'],c['W']),core.scale(b['P'],c['P'])))
    return dict(matter=core.gadd(*terms),
                electric=tuple(src['field']['electric']*c['E'] for c in coeff),
                magnetic=tuple(src['field']['magnetic']*c for c in magnetic))


def strain(src,basis,S):
    S=sy.Matrix(S)
    if S.shape!=(3,3) or S!=S.T or any(x.is_real is not True for x in S):
        raise ValueError('real symmetric three-dimensional strain required')
    terms=[];electric=[];tr=sy.trace(S)
    for b in basis['edges']:
        d=b['d'];vector=-S*d;quad=(d.T*S*d)[0]
        terms.extend(core.scale(p,w) for p,w in zip(b['D'],vector))
        terms.append(core.scale(b['W'],-2*quad))
        electric.append(src['field']['electric']*(2*quad-tr))
    magnetic=tuple(src['field']['magnetic']*(2*(m.T*S*m)[0]-tr) for m in basis['normals'])
    return dict(matter=core.gadd(*terms),electric=tuple(electric),magnetic=magnetic)


def shift_response(src,basis,axis):
    if axis not in range(3):raise ValueError('axis must be 0, 1 or 2')
    return dict(matter=core.gadd(*(core.scale(b['P'],b['d'][axis]) for b in basis['edges'])),
                electric=(sy.Integer(0),)*len(basis['edges']),magnetic=(sy.Integer(0),)*len(basis['normals']))


def apply_column(expr,src,key,cap=6):
    occupation,flux=key;graph=src['graph'];out=defaultdict(lambda:sy.Integer(0))
    if occupation<0 or occupation.bit_length()>96*graph.nvertices:raise ValueError('native occupation outside carrier')
    if len(flux)!=len(graph.edges) or any(abs(x)>cap for x in flux):raise ValueError('flux outside hard walls')
    for (word,shift),coefficient in expr['matter'].items():
        applied=core.apply_word(word,occupation)
        if applied is None:continue
        target=tuple(e+s for e,s in zip(flux,shift))
        if any(abs(e)>cap for e in target):continue
        mask,sign=applied;out[(mask,target)]+=sign*coefficient
    out[key]+=sum(c*e*e/2 for c,e in zip(expr['electric'],flux))
    for coefficient,row in zip(expr['magnetic'],src['field']['C']):
        out[key]+=coefficient
        for sign in (-1,1):
            target=tuple(e+sign*int(x) for e,x in zip(flux,row))
            if all(abs(e)<=cap for e in target):out[(occupation,target)]-=coefficient/2
    return core.clean(out)


def coefficient_audit():
    s,t=sy.symbols('s t',real=True)
    strains=[sy.diag(*(int(i==a) for i in range(3))) for a in range(3)]
    for a,b in ((0,1),(0,2),(1,2)):
        S=sy.zeros(3);S[a,b]=S[b,a]=sy.Rational(1,2);strains.append(S)
    d=sy.Matrix([sy.Rational(1,3),sy.Rational(2,3),sy.Rational(2,3)])
    edges=[dict(d=d)];normals=[sy.Matrix([0,0,1])]
    first=[];mixed=0
    for S in strains:
        coeff,mag=weights(sy.eye(3)+s*S,edges,normals)
        c=coeff[0];q=(d.T*S*d)[0]
        assert all(sy.simplify(sy.diff(c['D'][i],s).subs(s,0)+(S*d)[i])==0 for i in range(3))
        assert sy.simplify(sy.diff(c['W'],s).subs(s,0)+2*q)==0
        assert sy.simplify(sy.diff(c['E'],s).subs(s,0)-(2*q-sy.trace(S)))==0
        first.append([str(sy.diff(x,s).subs(s,0)) for x in (*c['D'],c['W'],c['E'],mag[0])])
        for T in strains:
            cc,mm=weights(sy.eye(3)+s*S+t*T,edges,normals)
            for f in (*cc[0]['D'],cc[0]['W'],cc[0]['E'],mm[0]):
                assert sy.simplify(sy.diff(f,s,t).subs({s:0,t:0})-sy.diff(f,t,s).subs({s:0,t:0}))==0
            mixed+=1
    return dict(six_strain_coefficients=first,mixed_strain_pairs=mixed,exact_first_and_second_derivatives=True)


def audit():
    result=coefficient_audit()
    src=source.assemble(spatial.edge());basis=bases(src);identity=expression(src,basis)
    assert identity['matter']==src['matter']
    assert identity['electric']==(src['field']['electric'],) and identity['magnetic']==()
    rejected=[]
    bad=dict(src);bad['components']=dict(src['components'])
    bad['components']['dirac']=core.scale(src['components']['dirac'],2)
    for name,call in (
        ('nondefault_Dirac_origin',lambda:validate_origin(bad,basis)),
        ('complex_coframe',lambda:expression(src,basis,sy.diag(sy.I,-sy.I,1))),
        ('complex_shift',lambda:expression(src,basis,shift=(sy.I,0,0))),
        ('complex_strain',lambda:strain(src,basis,sy.I*sy.eye(3)))):
        try:call()
        except ValueError:rejected.append(name)
        else:raise AssertionError('invalid input accepted: '+name)
    keys=[(0,(0,)),(1<<core.index(0,0,10,0),(0,)),((1<<192)-1,(6,))]
    for key in keys:assert apply_column(identity,src,key)==source.apply_source_column(src,*key,cap=6)
    strains=[sy.diag(*(int(i==a) for i in range(3))) for a in range(3)]
    for a,b in ((0,1),(0,2),(1,2)):
        S=sy.zeros(3);S[a,b]=S[b,a]=sy.Rational(1,2);strains.append(S)
    derivatives=[]
    for S in strains:
        e=strain(src,basis,S);check=core.check_gauged(e['matter'],2,src['graph'].edges)
        derivatives.append(dict(**check,digest=core.polynomial_digest(e['matter'])))
    shifts=[]
    for a in range(3):
        e=shift_response(src,basis,a);check=core.check_gauged(e['matter'],2,src['graph'].edges)
        shifts.append(dict(axis=a,**check,digest=core.polynomial_digest(e['matter'])))
    isotropic=strain(src,basis,sy.eye(3))
    D=core.gadd(*(core.scale(p,w) for b in basis['edges'] for p,w in zip(b['D'],b['d'])))
    W=core.gadd(*(b['W'] for b in basis['edges']))
    assert isotropic['matter']==core.gadd(core.scale(D,-1),core.scale(W,-2))
    assert isotropic['matter']!=core.scale(D,-1) and W
    # Compare the actual complete operator-column response at an unselected native input.
    h=sy.Rational(1,4096);plus=expression(src,basis,(1+h)*sy.eye(3));minus=expression(src,basis,(1-h)*sy.eye(3))
    key=keys[0];p=apply_column(plus,src,key);m=apply_column(minus,src,key);analytic=apply_column(isotropic,src,key)
    finite=core.scale(core.gadd(p,core.scale(m,-1)),1/(2*h))
    error=core.gadd(finite,core.scale(analytic,-1))
    maxerr=max((abs(complex(c)) for c in error.values()),default=0)
    assert maxerr<1e-5
    result.update(dict(parent_native_modes=192,parent_source_digest=core.polynomial_digest(src['matter']),
                       identity_member_exact=True,identity_complete_native_columns=len(keys),
                       strain_derivatives=derivatives,shift_derivatives=shifts,
                       isotropic_wilson_factor_exact=True,omitted_wilson_factor_detected=True,
                       central_difference_step=str(h),complete_column_derivative_error=maxerr,
                       explicit_domain_rejections=rejected,
                       scope='declared background mechanical responses; not an identified symmetric Lorentz stress'))
    return result
