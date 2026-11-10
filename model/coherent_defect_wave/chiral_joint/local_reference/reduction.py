"""Exact 140-state maximal quark module and actual lepton/scalar currents."""
from collections import defaultdict
from functools import lru_cache
from math import factorial
import sympy as s

from .algebra import (QUARKS,LEPTONS,transitions,quark_clifford,one_body_fock,
                      transformed_one_body,lepton_data)
from ...physical_attachment.chiral import representation
from ...chiral_higgs.model import PAULI


def coordinate_kernel(matrix):
    vectors=s.Matrix(matrix).nullspace()
    if not vectors:return s.ImmutableSparseMatrix(matrix.cols,0,{}),()
    basis=s.Matrix.hstack(*vectors)
    rows=tuple(basis.T.rref()[1])
    basis=s.ImmutableSparseMatrix(basis*basis.extract(rows,range(basis.cols)).inv())
    assert basis.extract(rows,range(basis.cols))==s.eye(basis.cols)
    assert matrix*basis==s.zeros(matrix.rows,basis.cols)
    return basis,rows


def restrict(operator,basis,coordinates):
    image=operator*basis
    result=s.ImmutableSparseMatrix(image.extract(coordinates,range(basis.cols)))
    assert basis*result==image
    return result


@lru_cache(maxsize=1)
def quark_current_matrices():
    rows,f,gamma,rotations,_,_=quark_clifford()
    generators=gamma+rotations
    full=tuple(m.extract(QUARKS,QUARKS) for m in transitions())
    full_gamma=tuple(m for i in (0,1) for m in
                     ((full[i]+full[i+2])/2,(full[i]-full[i+2])/(2*s.I)))
    full_generators=full_gamma+tuple(s.I*(a*b-b*a)
        for i,a in enumerate(full_gamma) for b in full_gamma[i+1:])
    actual=[]
    for current in (8,9,10,11):
        native=s.kronecker_product(s.eye(2),representation()[2][current])
        projected=s.ImmutableSparseMatrix(30,30,{(i,j):value for (i,j),value in native.todok().items()
                                                if i in QUARKS and j in QUARKS})
        constant,matrix=transformed_one_body(projected)
        if current==11:
            assert constant==-2
            matrix=matrix-s.diag(*(s.Rational(1,6) if i in QUARKS else 0 for i in range(30)))
        else:assert constant==0
        local=s.ImmutableSparseMatrix(matrix.extract(rows,rows))
        coefficients=[s.trace(t*local) for t in generators]
        assert sum((c*t for c,t in zip(coefficients,full_generators)),s.zeros(24))==matrix.extract(QUARKS,QUARKS)
        actual.append(local)
    return (*actual[:3],actual[3])


@lru_cache(maxsize=1)
def vector_five():
    rows,f,gamma,rotations,_,_=quark_clifford()
    generators=gamma+rotations
    fock_generators=tuple(one_body_fock(t) for t in generators)
    casimir=sum((t*t for t in fock_generators),s.zeros(16))
    coordinates=tuple(i for i in range(16) if i.bit_count()==2)
    small=casimir.extract(coordinates,coordinates)
    basis,pivots=coordinate_kernel(small-4*s.eye(6))
    assert basis.cols==5
    metric=s.ImmutableSparseMatrix(basis.T*basis)
    assert metric==s.diag(*metric.diagonal()) and all(v>0 for v in metric.diagonal())
    def lift(matrix):
        exterior=one_body_fock(matrix).extract(coordinates,coordinates)
        return restrict(exterior,basis,pivots)
    return {'transitions':tuple(lift(t) for t in f),
            'generators':tuple(lift(t) for t in generators),
            'currents':tuple(lift(t) for t in quark_current_matrices()),'metric':metric}


def occupations(total,dimension):
    if dimension==1:return ((total,),)
    return tuple((first,*rest) for first in range(total+1)
                 for rest in occupations(total-first,dimension-1))


@lru_cache(maxsize=1)
def symmetric_six():
    vector=vector_five();keys=occupations(6,5);positions={key:i for i,key in enumerate(keys)}
    assert len(keys)==210
    def lift(matrix):
        entries=defaultdict(lambda:s.Integer(0))
        for column,key in enumerate(keys):
            for (i,j),value in matrix.todok().items():
                if not key[j]:continue
                target=list(key);target[j]-=1;target[i]+=1
                entries[positions[tuple(target)],column]+=key[j]*value
        return s.ImmutableSparseMatrix(210,210,entries)
    diagonal=[]
    for key in keys:
        norm=s.Integer(1)
        for i,n in enumerate(key):norm*=factorial(n)*vector['metric'][i,i]**n
        diagonal.append(norm)
    metric=s.ImmutableSparseMatrix(s.diag(*diagonal))
    generators=tuple(lift(t) for t in vector['generators'])
    assert all(metric*t==t.H*metric for t in generators)
    casimir=s.ImmutableSparseMatrix(sum((t*t for t in generators),s.zeros(210)))
    return {'keys':keys,'metric':metric,'Casimir':casimir,
            'transitions':tuple(lift(t) for t in vector['transitions']),
            'currents':tuple(lift(t) for t in vector['currents'])}


@lru_cache(maxsize=1)
def quark_module():
    source=symmetric_six();casimir=source['Casimir'];weak=source['currents'][2];hyper=source['currents'][3]
    assert weak==s.diag(*weak.diagonal()) and hyper==s.diag(*hyper.diagonal())
    groups=defaultdict(list)
    for i in range(210):groups[(weak[i,i],hyper[i,i])].append(i)
    entries={};coordinates=[];weights=[]
    for weight,rows in sorted(groups.items()):
        basis,pivots=coordinate_kernel(casimir.extract(rows,rows)-54*s.eye(len(rows)))
        start=len(coordinates)
        for (i,j),value in basis.todok().items():entries[rows[i],start+j]=value
        coordinates.extend(rows[i] for i in pivots);weights.extend([weight]*basis.cols)
    basis=s.ImmutableSparseMatrix(210,len(coordinates),entries)
    assert basis.cols==140 and basis.extract(coordinates,range(140))==s.eye(140)
    assert casimir*basis==54*basis
    metric=s.ImmutableSparseMatrix(basis.T*source['metric']*basis)
    transitions=tuple(restrict(t,basis,coordinates) for t in source['transitions'])
    currents=tuple(restrict(t,basis,coordinates) for t in source['currents'])
    assert all(metric*t==t.H*metric for t in currents)
    assert metric*transitions[2]==transitions[0].T*metric
    assert metric*transitions[3]==transitions[1].T*metric
    assert currents[2]==s.diag(*(weight[0] for weight in weights))
    assert currents[3]==s.diag(*(weight[1] for weight in weights))
    return {'metric':metric,'transitions':transitions,'currents':currents,
            'weights':tuple((int(2*a),int(2*b)) for a,b in weights),
            'embedding':basis,'coordinates':tuple(coordinates)}


@lru_cache(maxsize=1)
def lepton_operators():
    currents=[]
    for current in (8,9,10,11):
        native=s.kronecker_product(s.eye(2),representation()[2][current])
        projected=s.ImmutableSparseMatrix(30,30,{(i,j):value for (i,j),value in native.todok().items()
                                                if i in LEPTONS and j in LEPTONS})
        constant,matrix=transformed_one_body(projected)
        assert constant==(2 if current==11 else 0)
        currents.append(s.ImmutableSparseMatrix(one_body_fock(matrix.extract(LEPTONS,LEPTONS))+constant*s.eye(64)))
    spin=[]
    for axis in PAULI:
        native=s.kronecker_product(axis/2,s.eye(15))
        projected=s.ImmutableSparseMatrix(30,30,{(i,j):value for (i,j),value in native.todok().items()
                                                if i in LEPTONS and j in LEPTONS})
        constant,matrix=transformed_one_body(projected);assert constant==0
        spin.append(one_body_fock(matrix.extract(LEPTONS,LEPTONS)))
    assert all(t==t.H for t in (*currents,*spin))
    assert all(t==s.diag(*t.diagonal()) for t in (currents[2],currents[3],spin[2]))
    return {'transitions':lepton_data()[0],'currents':tuple(currents),'spin':tuple(spin)}
