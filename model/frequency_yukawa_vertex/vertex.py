"""Exact ordered Euclidean responses and full-leg amputation.

Every intermediate state belongs to the unchanged 320-state source. Imaginary
frequencies are source-clock parameters, not physical momentum or GeV.
"""
from functools import lru_cache
from itertools import permutations
import hashlib
import json
from pathlib import Path

import sympy as s
import numpy as np
from scipy import linalg as la
from sympy.polys.domains import QQ_I
from sympy.polys.matrices import DomainMatrix

from model.yukawa_defect_vertex import vertex as native
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
FREQUENCIES=((1,2),(2,3),(4,8))
PAIRS=((0,4,0),(1,3,0),(0,5,1),(2,3,1))


def clean(matrix):
    return s.ImmutableSparseMatrix(matrix.rows,matrix.cols,
        {key:value for key,x in s.ImmutableSparseMatrix(matrix).todok().items()
         if (value:=s.expand(x))!=0})


def tr(matrix):
    return s.expand(sum(matrix[i,i] for i in range(matrix.rows)))


@lru_cache(maxsize=1)
def data():
    d=native.source_data();p=d['ground'];columns=[]
    for j in range(p.cols):
        if p[j,j]!=0:
            candidate=p[:,j]
            if not columns or s.Matrix.hstack(*columns,candidate).rank()>len(columns):
                columns.append(candidate)
            if len(columns)==2:break
    c=clean(s.Matrix.hstack(*columns));gram=clean(c.H*c)
    weight=clean(gram.inv()/2)
    if clean(c*weight*c.H-d['rho']).todok():raise ArithmeticError('ground factorization changed')
    operators=tuple(d['c'])+tuple(x.H for x in d['phi'])
    return dict(**d,columns=c,weight=weight,operators=operators)


@lru_cache(maxsize=24)
def resolvent(omega):
    """(i omega-(H-E0))^-1, including the 198-dimensional dark space."""
    if type(omega) is not int or omega==0:raise ValueError('nonzero integer source frequency required')
    forms,occupied=native.spectral_trace_forms();z=s.I*omega-2
    dark=occupied-sum((forms[r][2]/r for r in range(1,7)),s.zeros(320))
    result=forms[0][0]*s.expand_complex(1/z)+dark*s.expand_complex(1/(z-1))
    for r in range(1,7):
        a,b,c,_=forms[r]
        inverse=s.expand_complex(1/(z*(z-1)-r))
        result+=(a*(z-1)+b+c*z/r)*inverse
    result=clean(result)
    if clean((s.I*omega*native.eye(320)-data()['H']-2*native.eye(320))*result-native.eye(320)).todok():
        raise ArithmeticError('full source resolvent identity failed')
    return result


@lru_cache(maxsize=256)
def right_leg(operator,omega):
    d=data();return clean(resolvent(omega)*d['operators'][operator]*d['columns'])


@lru_cache(maxsize=2048)
def double_leg(middle,last,first_frequency,last_frequency):
    return clean(resolvent(first_frequency)*data()['operators'][middle]*right_leg(last,-last_frequency))


def ordered_three_point(i,j,k,omega1,omega2):
    """Sum all six time orderings, with the actual fermion exchange sign."""
    d=data();ops=(i,j,6+k);freq=(omega1,omega2,-omega1-omega2);value=0;pieces=[]
    for order in permutations(range(3)):
        a,b,c=order;sign=-1 if order.index(0)>order.index(1) else 1
        matrix=d['columns'].H*d['operators'][ops[a]]*double_leg(ops[b],ops[c],freq[a],freq[c])
        term=sign*tr(d['weight']*matrix)
        value+=term;pieces.append(term)
    return s.expand(value),tuple(pieces)


@lru_cache(maxsize=16)
def propagators(omega):
    d=data();r=resolvent(omega);rm=resolvent(-omega)
    def expectation(a):return native.trace_product(d['rho'],a)
    g=s.zeros(6);boson=s.zeros(2)
    for i,ci in enumerate(d['c']):
        for j,cj in enumerate(d['c']):
            g[i,j]=s.expand(expectation(ci*r*cj.H-cj.H*rm*ci))
    for i,phi in enumerate(d['phi']):
        for j,psi in enumerate(d['phi']):
            boson[i,j]=s.expand(expectation(phi*r*psi.H+psi.H*rm*phi))
    return clean(g),clean(boson)


def encoded(value):
    value=s.expand(value)
    return {'exact':str(value),'real':float(s.re(value)),'imaginary':float(s.im(value))}


def frequency_row(w1,w2):
    g1,_=propagators(w1);g2,_=propagators(w2);_,dh=propagators(-w1-w2)
    if dh!=propagators(w1+w2)[1]:raise ArithmeticError('conjugate scalar leg requires a different orientation')
    if not g1.det() or not g2.det() or not dh.det():raise ArithmeticError('full-leg amputation is singular')
    if g1!=s.diag(*g1.diagonal()) or g2!=s.diag(*g2.diagonal()) or dh!=dh[0,0]*s.eye(2):
        raise ArithmeticError('source representation does not have the verified diagonal leg frame')
    rows=[]
    for i,j,k in PAIRS:
        value,pieces=ordered_three_point(i,j,k,w1,w2)
        exchanged,_=ordered_three_point(j,i,k,w2,w1)
        if s.expand(value+exchanged)!=0:raise ArithmeticError('frequency-exchanged fermion antisymmetry failed')
        gamma=s.expand(value/(g1[i,i]*g2[j,j]*dh[k,k]))
        rows.append({'left':i,'right':j,'scalar':k,'complete_six_ordered_terms':[encoded(x) for x in pieces],
            'connected_three_point':encoded(value),'amputated_coefficient':encoded(gamma),
            'canonically_normalized_amputated_coefficient':encoded(s.sqrt(s.Rational(14,5))*gamma)})
    if any(row['connected_three_point']['exact']=='0' for row in rows):raise ArithmeticError('native finite-frequency response vanished')
    return {'fermion_frequencies':[w1,w2],'scalar_frequency':-w1-w2,
        'fermion_leg1_diagonal':[encoded(x) for x in g1.diagonal()],
        'fermion_leg2_diagonal':[encoded(x) for x in g2.diagonal()],
        'scalar_leg':encoded(dh[0,0]),'components':rows,
        'all_six_time_orderings_retained':True,'full_320_state_resolvent_identities':True,
        'all_external_legs_invertible':True,'scalar_canonical_moment':'14/5'}


def complete_tensor():
    """Every component at the predeclared first frequency pair, including zeros."""
    rows=[];allowed={(i,j,k) for i,j,k in PAIRS}|{(j,i,k) for i,j,k in PAIRS}
    for i in range(6):
        for j in range(6):
            for k in range(2):
                value,_=ordered_three_point(i,j,k,1,2)
                if (value!=0)!=((i,j,k) in allowed):raise ArithmeticError('full native tensor support differs')
                rows.append({'left':i,'right':j,'scalar':k,'value':encoded(value)})
    return rows


def complete_pointer():
    """Complete bounded quadrature instrument for the invariant three-point sum."""
    d=data();freq=(1,2,-3);normalizer=s.S.Zero
    # Exact sparse Gaussian-rational multiplication avoids expression-tree expansion.
    domain=lambda a:DomainMatrix.from_Matrix(a).convert_to(QQ_I)
    rs={omega:domain(resolvent(omega)) for omega in (*freq,*(-w for w in freq))}
    m=DomainMatrix.zeros((320,320),QQ_I)
    for i,j,k in PAIRS:
        ops=tuple(domain(a) for a in (d['c'][i],d['c'][j],d['phi'][k].H))
        for order in permutations(range(3)):
            a,b,c=order;sign=-1 if order.index(0)>order.index(1) else 1
            term=ops[a].matmul(rs[freq[a]]).matmul(ops[b]).matmul(rs[-freq[c]]).matmul(ops[c])
            m=m.add(term.mul(QQ_I(sign)))
            normalizer+=s.Rational(1,abs(freq[a]*freq[c]))
    m=clean(m.to_Matrix()/normalizer)
    if normalizer!=8:raise ArithmeticError('predeclared contraction bound changed')
    for operator in d['operators']:
        effect=clean(operator.H*operator)
        if clean(effect*effect-effect).todok():raise ArithmeticError('field is not a norm-one partial isometry')
    ward=[not clean(q*m-m*q).todok() for q in d['currents']]
    if not all(ward):raise ArithmeticError('invariant three-point pointer changed gauge charge')
    real=clean((m+m.H)/2);imaginary=clean((m-m.H)/(2*s.I))
    effects=[clean((native.eye(320)+sign*a)/4) for a in (real,imaginary) for sign in (1,-1)]
    if clean(sum(effects,s.zeros(320))-native.eye(320)).todok():raise ArithmeticError('quadrature effects are incomplete')
    kraus=[];minimum=[]
    for effect in effects:
        numerical=np.asarray(effect,complex);eig,vectors=la.eigh(numerical)
        minimum.append(float(eig[0]))
        if eig[0]<0:raise ArithmeticError('positive contraction effects failed')
        kraus.append((vectors*np.sqrt(eig))@vectors.conj().T)
    gram=sum((k.conj().T@k for k in kraus),np.zeros((320,320),complex))
    numerical_error=float(la.norm(gram-np.eye(320),2))
    if numerical_error>2e-12:raise ArithmeticError('full quantum-output instrument failed')
    records=[]
    for energy,_,rank in native.spectral_weights(d['rho']):
        if energy!=-2:records.append({'stage':'preparation_other_energy','energy':str(energy),'probability':str(s.Rational(rank,320))})
    for label,effect in zip(('real_plus','real_minus','imaginary_plus','imaginary_minus'),effects,strict=True):
        records.append({'stage':'ground_then_pointer','label':label,'probability':str(native.trace_product(d['rho'],effect)/160)})
    if sum(s.sympify(row['probability']) for row in records)!=1:raise ArithmeticError('original preparation weight was lost')
    expected=sum((ordered_three_point(i,j,k,1,2)[0] for i,j,k in PAIRS),s.S.Zero)/normalizer
    if s.expand(native.trace_product(d['rho'],m)-expected)!=0:raise ArithmeticError('ordered trace and full pointer disagree')
    operator_rows=[[i,j,str(value)] for (i,j),value in sorted(m.todok().items())]
    return {'source_dimension':320,'contraction_bound_normalizer_exact':str(normalizer),
        'invariant_complex_pointer_expectation':encoded(expected),'total_current_commutators_exact':ward,
        'four_positive_effects_sum_identity_exact':True,'all_source_and_reference_outputs_retained':True,
        'numerical_sqrt_instrument_identity_error':numerical_error,'minimum_effect_eigenvalues':minimum,
        'complete_original_records':records,'record_count':len(records),
        'pointer_matrix_sha256':hashlib.sha256(json.dumps(operator_rows,separators=(',',':')).encode()).hexdigest(),
        'apparatus':'declared finite square-root effects, not a compiled native pulse sequence',
        'zero_operator_control_four_probabilities':['1/4']*4}


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in {'.py','.md'}}
    paths.update(native.SOURCE_FILES)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def certificate():
    pins=bindings();rows=[frequency_row(*pair) for pair in FREQUENCIES]
    d=data()
    if any(native.trace_product(d['rho'],o)!=0 for o in (*d['c'],*d['phi'])):
        raise ArithmeticError('one-point disconnected terms do not vanish')
    if any(clean(d['ground']*o*d['ground']).todok() for o in d['operators']):
        raise ArithmeticError('charged external legs act inside the singlet ground space')
    result={'schema':'oph.frequency_yukawa_vertex.v1','source_bindings':pins,'source_dimension':320,
        'ground_rank':2,'ground_energy':-2,'ground_preparation_probability':'1/160','frequency_rows':rows,
        'complete_first_frequency_tensor':complete_tensor(),'complete_invariant_pointer':complete_pointer(),
        'all_one_point_functions_zero_exact':True,
        'primitive_interaction_coefficient':1,'family_matrix_supplied':False,
        'physical_Yukawa_or_on_shell_mass_identification':False,
        'source_clock':'exp(-itH), with the literal unit native lepton/Higgs law',
        'scope':'Finite Euclidean connected three-point function and explicit full-leg amputation; no spatial momentum, physical renormalization scheme or selected Higgs vacuum is inferred.'}
    if bindings()!=pins:raise ArithmeticError('source changed during finite-frequency calculation')
    return result
