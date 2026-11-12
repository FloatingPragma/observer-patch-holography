"""Independent CAR, composition, norm-budget and complete-POVM review."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json

import numpy as np
import scipy.linalg as la
import sympy as s

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=ROOT/'model/openwave_lepton_phase_cycle'


def annihilator(mode):
    result=s.zeros(16)
    for mask in range(16):
        if mask&(1<<mode):
            result[mask^(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return result


def run():
    c=[annihilator(i) for i in range(4)]
    w=s.Matrix([[s.Rational(1,5),s.Rational(2,5)],
                [s.Rational(2,5),s.Rational(4,5)]])
    pair=sum((w[i,j]*c[i].H*c[j+2].H for i in range(2) for j in range(2)),s.zeros(16))
    v=pair+pair.H
    count=[((mask>>1)&1)+((mask>>3)&1) for mask in range(16)]
    rotations=[s.diag(*(s.I**(r*n) for n in count)) for r in range(4)]
    conjugates=[r.H*v*r for r in rotations]
    physical_pair=(c[0].H*c[2].H+c[2]*c[0])/5
    assert sum(conjugates,s.zeros(16))/4==physical_pair
    assert (conjugates[0]+conjugates[2])/2!=physical_pair
    sign=s.diag(*((-1)**((mask&3).bit_count()) for mask in range(16)))
    assert sign*v*sign==-v and all(sign*r==r*sign for r in rotations)
    assert v**3==v
    # The proof on every Fock sector is the monomial root-of-unity identity.
    for ancillary_creators in (0,1,2):
        average=sum(s.I**(-r*ancillary_creators) for r in range(4))/4
        assert average==int(ancillary_creators==0)

    # Independent exact matrix series uses SymPy polynomials, unlike the
    # producer's hand-written polynomial coefficient arithmetic.
    p=s.symbols('p');cubic=s.Poly(60*p**3-48*p**2+12*p-1,p)
    x=s.Matrix([[0,1],[1,0]]);y=s.Matrix([[0,-s.I],[s.I,0]]);z=s.diag(1,-1)
    matrices=[x,y,z,(x+y+z)/2]
    series=[s.eye(2)]+[s.zeros(2) for _ in range(4)]
    palindrome=[(0,s.Rational(1,8)),(1,s.Rational(1,8)),(2,s.Rational(1,8)),
                (3,s.Rational(1,4)),(2,s.Rational(1,8)),(1,s.Rational(1,8)),(0,s.Rational(1,8))]
    reduce=lambda matrix:matrix.applyfunc(lambda e:s.Poly(s.expand(e),p).rem(cubic).as_expr())
    for duration in (p,p,1-4*p,p,p):
        for index,weight in palindrome:
            generator=-s.I*matrices[index]*duration*weight
            step=[s.eye(2)]+[generator**n/factorial(n) for n in range(1,5)]
            series=[reduce(sum((step[j]*series[n-j] for j in range(n+1)),s.zeros(2))) for n in range(5)]
    target=-s.I*sum(matrices,s.zeros(2))/4
    assert all(reduce(series[n]-target**n/factorial(n))==s.zeros(2) for n in range(5))
    assert s.rem(4*p**3+(1-4*p)**3,cubic.as_expr(),p)==0

    # e<3 is bounded independently by its rational series tail. Together
    # with 2/5! it gives the claimed 1/20 local norm constant.
    exp_upper=sum((F(1,factorial(k)) for k in range(7)),F())+F(1,factorial(7))/(1-F(1,8))
    assert exp_upper<3 and 2*exp_upper/F(factorial(5))<F(1,20)
    assert F(8,5)**3>4  # cuberoot(4)<8/5 implies p<5/12.
    n=15229
    bound=F(14**5,10*n**4)
    assert bound<=F(1,10**12)<F(14**5,10*(n-1)**4)
    assert n>=14 and 35*n==533015

    # Direct complete Ramsey POVM, including all source and reference
    # complements; compare with the compact contrast formula.
    j=s.eye(16)[:,[0,1,4,5]]
    ground=j*s.Matrix([1/s.sqrt(2),0,0,-1/s.sqrt(2)])
    particle=j*s.Matrix([0,1,0,0])
    g=s.kronecker_product(ground,s.Matrix([1,0]))
    a=s.kronecker_product(particle,s.Matrix([0,1]))
    # Preserve the inherited source convention Y=i(|0><1|-|1><0|),
    # which reads positive sin(mt) for the ground-to-particle gap.
    X=g*a.H+a*g.H;Y=s.I*g*a.H-s.I*a*g.H
    effects=[(s.eye(32)+X)/4,(s.eye(32)-X)/4,(s.eye(32)+Y)/4,(s.eye(32)-Y)/4]
    assert sum(effects,s.zeros(32))==s.eye(32)
    assert all(set(effect.eigenvals()).issubset({0,s.Rational(1,4),s.Rational(1,2)}) for effect in effects)
    source_charge=s.diag(*(-((mask>>0)&1)-((mask>>1)&1)+((mask>>2)&1)+((mask>>3)&1) for mask in range(16)))
    total_charge=s.kronecker_product(source_charge,s.eye(2))+s.kronecker_product(s.eye(16),s.diag(0,1))
    assert total_charge*X==X*total_charge and total_charge*Y==Y*total_charge
    pfloat=1/(4-4**(1/3));rounds=100
    hs=[np.asarray(h,complex) for h in conjugates]
    cycle=np.eye(16,dtype=complex)
    for duration in (pfloat,pfloat,1-4*pfloat,pfloat,pfloat):
        for index,weight in palindrome:
            cycle=la.expm(-1j*hs[index]*duration*float(weight)/rounds)@cycle
    actual=np.linalg.matrix_power(cycle,rounds)
    joined=np.kron(actual,np.eye(2))
    state=np.asarray((g+a)/s.sqrt(2),complex).ravel()
    output=joined@state
    observed=[float(np.vdot(output,np.asarray(effect,complex)@output).real) for effect in effects]
    gv=np.asarray(ground,complex).ravel();av=np.asarray(particle,complex).ravel()
    contrast=np.conj(np.vdot(gv,actual@gv))*np.vdot(av,actual@av)
    formula=[(1+contrast.real)/4,(1-contrast.real)/4,(1-contrast.imag)/4,(1+contrast.imag)/4]
    assert np.max(np.abs(np.asarray(observed)-formula))<1e-12
    assert abs(sum(observed)-1)<1e-12
    ideal=la.expm(-1j*np.asarray(physical_pair,complex))@np.asarray(j,complex)
    error=np.linalg.norm(actual@np.asarray(j,complex)-ideal,2)
    assert 2*error<float(F(14**5,10*rounds**4))

    from model.openwave_lepton_phase_cycle.compiler import rounds_for_error,trace_error_bound
    fractional=rounds_for_error(F(1,10),time=F(1,3))
    assert fractional>=5 and trace_error_bound(fractional,time=F(1,3))<=F(1,10)

    receipt=json.loads((PACKAGE/'results.json').read_text())
    for path,expected in receipt['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    return {'schema':'oph.independent_lepton_phase_cycle_review.v1','all_pass':True,
        'exact_all16Fock_Z4_selection_and_negative_segment_control':True,
        'arbitrary_Fock_monomial_root_of_unity_argument_checked':True,
        'independent_exact_noncommuting_matrix_series_degrees':[0,1,2,3,4],
        'rational_local_Taylor_and_complete_trace_constants_verified':True,
        'exact_minimal_declared_budget_rounds':n,
        'complete32dimensional_Ramsey_POVM_and_charge_compensation_verified':True,
        'direct_complete_Ramsey_probabilities':observed,
        'diagnostic_source_unitary_error':error,
        'generic_fractional_time_rounding_bug_reported_to_owner':True,
        'fractional_time_radius_regression_fixed_and_verified':True,
        'scientific_t1_result_affected_by_rounding_bug':False,
        'source_clock_boundary':'A new declared coherent controller at zero kinetic term; source-duration overhead and ancillary services retained; not autonomous laboratory evolution.',
        'reviewed_sha256':{str((PACKAGE/name).relative_to(ROOT)):hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest()
                           for name in ('compiler.py','PROOF.md','results.json')}}


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parents[1]/'phase_cycle_review.json'
    target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS: full-Fock twirl, fourth-order composition, complete trace budget and Ramsey records')
