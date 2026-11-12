"""Independent source, Gaussian-energy and finite-minimum review.

No implementation from openwave_native_force is imported. Native scalar
matrices, sparse graph solves and high-precision trial states are rebuilt
here. Floating solves only propose exact dyadic or rational witnesses.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib
import json

import mpmath as mp
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg
import sympy as s

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / 'model/openwave_native_force'


def rayleigh(matrix, vector):
    return F(s.cancel((vector.T*matrix*vector)[0]/(vector.T*vector)[0]))


def root_upper(value):
    scale=10**30
    numerator=isqrt(value.numerator*scale*scale//value.denominator)+1
    result=F(numerator,scale)
    assert result*result>=value
    return result


def graph(n):
    axis=sp.diags((-np.ones(n-1),2*np.ones(n),-np.ones(n-1)),(-1,0,1),format='csr')
    eye=sp.eye(n,format='csr')
    return (sp.kron(sp.kron(axis,eye),eye)+sp.kron(sp.kron(eye,axis),eye)
            +sp.kron(sp.kron(eye,eye),axis)).tocsr()


def independent_minimum(n, sites, charge=F(9,140000)):
    matrix=graph(n)
    density=np.zeros(n**3)
    indices=[np.ravel_multi_index(site,(n,n,n)) for site in sites]
    for index in indices:density[index]+=float(charge)
    trial,info=cg(matrix,-density,rtol=1e-13,atol=0,maxiter=1000)
    assert info==0
    scale=2**48
    integers=np.rint(trial*scale).astype(np.int64)
    maximum=max(abs(int(integers.min())),abs(int(integers.max())))
    assert F(maximum,scale)<F(1,50000)
    # Matrix entries are cast to integers before exact residual evaluation.
    action=matrix.astype(np.int64)@integers
    residual=charge.denominator*action
    for index in indices:residual[index]+=charge.numerator*scale
    bound=max(abs(int(residual.min())),abs(int(residual.max())))
    norm_upper=F(n**3*bound*bound,charge.denominator**2*scale**2)
    energy=F(sum(int(a)*int(b) for a,b in zip(integers,action)),2*scale**2)
    energy+=charge*sum(F(int(integers[index]),scale) for index in indices)
    lower_eigenvalue=F(12,(n+1)**2)
    return energy-norm_upper/(2*lower_eigenvalue),energy


def difference(pair,one,two):
    return pair[0]-one[1]-two[1],pair[1]-one[0]-two[0]


def run():
    saved=json.loads((PACKAGE/'results.json').read_text())
    for path,expected in saved['source_bindings'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path

    from model.native_equilibration import vacuum_repair as source
    from model.finite_area_response.response import spectral_enclosures
    r=s.symbols('r',positive=True)
    terms=[s.Matrix(term.astype(int)) for term in source._terms(2)]
    actual=sum((s.simplify(c/s.sqrt(14))*term for c,term in
                zip(source._coefficient_expressions(2,r),terms)),s.zeros(8))
    m=actual.subs(r,1);j=actual.diff(r).subs(r,1)
    contrast=s.factor(actual[1,1]-actual[0,0])
    assert contrast.subs(r,1)==1
    assert s.diff(contrast,r).subs(r,1)==s.Rational(9,14)
    assert F(saved['force']['native_dimensionless_geometry_charge'])==F(9,14)
    assert j[0,3]==s.Rational(1,7)  # Sign reverses the energy off-diagonal.
    assert m[0,3]==-s.Rational(1,7)
    assert j!=s.Rational(9,14)*m  # Energy and response are not proportional.
    assert max(sum(abs(j[a,b]) for b in range(8)) for a in range(8))==s.Rational(89,28)

    q=s.symbols('q',real=True)
    remainder=s.factor(contrast.subs(r,1+q)-1-s.Rational(9,14)*q)
    coefficient=s.cancel(remainder/q**2)
    assert s.simplify(coefficient-(7+3*q)/(14*(1+q)))==0
    assert s.simplify(s.diff(coefficient,q)+s.Rational(2,7)/(1+q)**2)==0
    window=F(1,50000);mixture=F(1,10000)
    bound=F(coefficient.subs(q,-s.Rational(window)))
    error=2*mixture*bound*window**2
    assert error==F(saved['force']['uniform_nonlinear_cross_energy_error_upper'])
    assert all(remainder.subs(q,x)>=0 for x in (-s.Rational(window),0,s.Rational(window)))

    # Exact 27-coordinate graph elimination gives an independent sign and
    # factor-of-two check, including the entire boundary diagonal.
    small=s.Matrix(graph(3).toarray().astype(int))
    inverse=small.inv()
    a,b=13,14
    rho=s.eye(27)[:,a]+s.eye(27)[:,b]
    solution=-inverse*rho
    energy=(solution.T*small*solution)[0]/2+(rho.T*solution)[0]
    single_a=-inverse[a,a]/2;single_b=-inverse[b,b]/2
    cross=s.factor(energy-single_a-single_b)
    assert cross==-inverse[a,b]<0
    assert inverse[a,b]>0
    assert cross!=-inverse[a,b]/2  # Dropping one cross term fails.
    opposite=s.eye(27)[:,a]-s.eye(27)[:,b]
    opposite_min=-(opposite.T*inverse*opposite)[0]/2
    assert opposite_min-single_a-single_b==inverse[a,b]>0
    assert (2*small).inv()==inverse/2
    screened=(small+s.eye(27)).inv()
    assert 0<screened[a,b]<inverse[a,b]
    test=s.Matrix([s.Rational((i*7)%11-5,17) for i in range(27)])
    residual=small*test+rho
    trial_energy=(test.T*small*test)[0]/2+(rho.T*test)[0]
    assert s.factor(trial_energy-energy-(residual.T*inverse*residual)[0]/2)==0

    # Alternate sparse-CG witnesses independently enclose the same finite
    # energy minima; each witness is converted to exact integer arithmetic.
    n=31;center=(15,15,15)
    one=independent_minimum(n,(center,))
    crosses={}
    for distance in (5,6,7):
        position=(15+distance,15,15)
        pair=independent_minimum(n,(center,position))
        second=independent_minimum(n,(position,))
        linear=difference(pair,one,second)
        crosses[distance]=(linear[0]-error,linear[1]+error)
        if distance==6:
            recorded=[F(x) for x in saved['force']['finite_carrier_cross_energy_interval']]
            assert max(recorded[0],linear[0])<=min(recorded[1],linear[1])
            assert crosses[distance][1]<0
    attraction=((crosses[7][0]-crosses[5][1])/2,
                (crosses[7][1]-crosses[5][0])/2)
    recorded=[F(x) for x in saved['force']['full_nonlinear_attractive_force_interval_at_six']]
    assert attraction[0]>0
    assert max(recorded[0],attraction[0])<=min(recorded[1],attraction[1])

    # Ground-relative obstruction: propose a new high-precision reference,
    # then certify its expectation using rational Rayleigh/gap estimates.
    spectrum=spectral_enclosures(m)
    mp.mp.dps=75
    matrix_mp=mp.matrix([[mp.mpf(str(m[i,k].p))/int(m[i,k].q) for k in range(8)] for i in range(8)])
    _,vectors=mp.eigsy(matrix_mp)
    denominator=10**65
    v=s.Matrix([s.Rational(int(mp.nint(vectors[i,0]*denominator)),denominator) for i in range(8)])
    energy=rayleigh(m,v);expectation=rayleigh(j,v)
    gap=spectrum[1].lo-spectrum[0].hi
    infidelity=(energy-spectrum[0].lo)/gap
    assert 0<=infidelity<F(1,10**18)
    charge_error=2*F(89,28)*root_upper(infidelity)
    checked=[]
    for row in saved['actual_ground_nonuniversality']['opposite_charge_positive_energy_states']:
        state=s.Matrix([s.Rational(value) for value in row['binary_rational_state']])
        excess=rayleigh(m,state)-spectrum[0].hi
        charge=rayleigh(j,state)-expectation
        assert excess>0
        assert (charge+charge_error<0) if row['sign']<0 else (charge-charge_error>0)
        checked.append({'sign':row['sign'],'energy_lower':str(excess),
                        'geometry_charge_interval':[str(charge-charge_error),str(charge+charge_error)]})

    scope=saved['scope']
    assert scope['spatial_elastic_geometry_law_is_declared_extension']
    assert not scope['universal_Newton_gravity_or_G_derived']
    assert not scope['native_controller_implements_elastic_minimum_preparation']
    paths=[PACKAGE/name for name in ('force.py','PROOF.md','packet.py','results.json')]
    return {'schema':'oph.independent_native_force_review.v1','all_pass':True,
            'reviewed_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
            'source_energy_and_response_rebuilt_from_original_matrices':True,
            'exact_27_coordinate_gaussian_sign_and_normalization_verified':True,
            'full_nonlinear_remainder_positive_on_entire_declared_carrier':True,
            'independent_sparse_exact_witness_force_interval':[str(x) for x in attraction],
            'independent_force_decimal_interval':[float(x) for x in attraction],
            'actual_ground_counterexamples_independently_certified':checked,
            'mutation_controls':['wrong derivative off-diagonal sign rejected',
                                 'half interaction normalization rejected',
                                 'opposite charge yields repulsion',
                                 'stiffness rescaling changes force',
                                 'onsite mass screens response',
                                 'universal positive-energy geometry charge refuted'],
            'scope':'Native matter derivative plus declared finite spatial elastic extension; no native minimizer controller, physical G, or universal mass coupling.'}


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parents[1]/'native_force_review.json'
    target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: independent source derivative, finite nonlinear force and ground-relative obstruction')
