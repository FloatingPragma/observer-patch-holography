"""Independent metric reconstruction, clocks and complete field-record audit."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib
import json
import mpmath as mp
import numpy as np
import sympy as s
from scipy.linalg import eigh,expm
from scipy.sparse import diags,eye,kron
from scipy.sparse.linalg import cg

ROOT=Path(__file__).resolve().parents[3]
PACKAGE=ROOT/'model/openwave_metric_observables'


def run():
    saved=json.loads((PACKAGE/'results.json').read_text())
    for path,digest in saved['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    mp.mp.dps=85
    number=lambda x:mp.mpf(F(x).numerator)/F(x).denominator
    enclosed=lambda value,bounds:number(bounds[0])<=value<=number(bounds[1])
    radii=list(map(F,saved['metric']['transverse_slice_r']))
    metric=saved['metric']
    for row in metric['local_metric_rows']:
        r=s.Rational(row['r']);C=s.diag(*map(s.Rational,row['C']))
        g=s.diag(*map(s.Rational,row['g']))
        assert C==s.diag(-1/r,1/r,r*r,r*r)
        assert -g.det()==r*r and C==s.sqrt(-g.det())*g.inv()
        assert -g[0,0]/g[1,1]==1 and -g[0,0]/g[2,2]==r**3
    travel=mp.fsum(number(r)**(-mp.mpf('1.5')) for r in radii)
    assert enclosed(travel,metric['piecewise_metric_transverse_ray_time_interval'])
    assert enclosed(travel-31,metric['piecewise_metric_transverse_delay_interval'])

    # A separately assembled sparse positive Laplacian reproduces the dyadic
    # geometry witness used by the spectral producer, without using its FFT.
    n=31;one=eye(n,format='csr');line=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1],format='csr')
    L=kron(kron(line,one),one)+kron(kron(one,line),one)+kron(kron(one,one),line)
    charge=np.zeros(n**3);charge[np.ravel_multi_index((15,15,15),(n,n,n))]=-float(F(9,140000))
    q,info=cg(L,charge,rtol=1e-14,atol=0);assert info==0
    dyadic=np.rint(q.reshape((n,n,n))*2**48).astype(np.int64)
    assert hashlib.sha256(dyadic.astype('<i8').tobytes()).hexdigest()==metric['integer_geometry_sha256']

    # Rebuild M directly from the oscillator terms, not the metric-clock
    # producer or its source_matrices helper. All eight eigenlevels remain.
    from model.native_equilibration import vacuum_repair as native
    symbol=s.symbols('r',positive=True)
    terms=[s.Matrix(term.astype(int)) for term in native._terms(2)]
    M=sum((a*b for a,b in zip(native._coefficient_expressions(2,symbol),terms)),s.zeros(8))/s.sqrt(14)
    def spectrum(r):
        exact=M.subs(symbol,s.Rational(r)).applyfunc(s.simplify)
        matrix=mp.matrix([[mp.mpf(str(s.N(exact[i,j],90))) for j in range(8)] for i in range(8)])
        return list(mp.eigsy(matrix,eigvals_only=True))
    eigs=[spectrum(radii[j]) for j in (15,23)]
    lapse=number(radii[15]/radii[23]);ratios=[]
    for row in saved['native_clocks']['rows']:
        level=row['excited_level'];gaps=[ev[level]-ev[0] for ev in eigs]
        for gap,bounds,prob in zip(gaps,row['gap_over_sqrt14_intervals'],row['Ramsey_plus_Y_probability_intervals']):
            assert enclosed(gap,bounds)
            assert enclosed((1-mp.sin(gap))/2,prob)
        ratio=gaps[0]/gaps[1];ratios.append(ratio)
        assert enclosed(ratio,row['near_to_far_frequency_ratio_interval'])
        assert enclosed(ratio-lapse,row['ratio_minus_metric_lapse_interval'])
        assert abs(ratio-lapse)>mp.mpf('1e-8')
    assert all(abs(a-b)>mp.mpf('1e-8') for j,a in enumerate(ratios) for b in ratios[j+1:])

    # Canonical masses give a direct tridiagonal K. Independently exponentiate
    # its positive square root on all sites; retain vacuum and all31 records.
    wave=saved['spatial_scalar_wave'];position=np.arange(31)
    initial=np.exp(-(position-7)**2/8+0.7j*position);initial/=np.linalg.norm(initial)
    for run in wave['runs']:
        records=[]
        for name,rs in (('curved',radii),('flat_control',[F(1)]*31)):
            r=np.asarray(rs,float);extended=np.r_[1.,r,1.]
            a=(extended[:-1]**2+extended[1:]**2)/2
            K=np.diag(r*(a[:-1]+a[1:]))+np.diag(-a[1:-1]*np.sqrt(r[:-1]*r[1:]),1)+np.diag(-a[1:-1]*np.sqrt(r[:-1]*r[1:]),-1)
            values,vectors=eigh(K);assert values[0]>0 and values[-1]<4
            H=(vectors*np.sqrt(values))@vectors.T
            row=run[name];time=row['source_time'];psi=expm(-1j*time*H)@initial
            actual=np.r_[.5,.5*np.abs(psi)**2];records.append(actual)
            assert len(actual)==32 and abs(sum(actual)-1)<2e-15
            assert np.max(np.abs(actual-row['complete_vacuum_and_position_probabilities']))<2e-13
            assert np.max(np.abs(np.sqrt(values)-row['normal_frequencies']))<2e-14
            x=F(3*time,16);eps=x**37/F(factorial(37))/(1-x/38)
            assert F(row['analytic_evolution_operator_error_bound'])==16*eps*(1+eps)**15
        tv=np.abs(records[0]-records[1]).sum()/2
        assert abs(tv-run['complete_output_total_variation'])<2e-13 and tv>1e-7

    # Compute linearized curvature via Christoffel symbols, independently of
    # the producer's direct linearized-Ricci formula.
    t,x,y,z=s.symbols('t x y z');coords=(t,x,y,z);qf=s.Function('q')(x,y,z)
    eta=s.diag(-1,1,1,1);h=s.diag(-2*qf,2*qf,-qf,-qf)
    Gamma=[[[sum(eta[a,d]*(s.diff(h[d,c],coords[b])+s.diff(h[d,b],coords[c])-s.diff(h[b,c],coords[d])) for d in range(4))/2 for c in range(4)]for b in range(4)]for a in range(4)]
    Ric=s.Matrix(4,4,lambda a,b:sum(s.diff(Gamma[c][a][b],coords[c])-s.diff(Gamma[c][a][c],coords[b])for c in range(4)))
    G=s.simplify(Ric-eta*s.trace(eta*Ric)/2)
    expected=s.diff(qf,x,2)-(s.diff(qf,y,2)+s.diff(qf,z,2))/2
    assert s.simplify(G[0,0]-expected)==0
    assert expected.subs(qf,-1/s.sqrt(x*x+y*y+z*z)).doit().subs({x:3,y:0,z:0})==-s.Rational(1,9)
    point=(21,15,15);second=[]
    for axis in range(3):
        left=list(point);right=list(point);left[axis]-=1;right[axis]+=1
        second.append(F(int(dyadic[tuple(left)])-2*int(dyadic[point])+int(dyadic[tuple(right)]),2**48))
    einstein=saved['Einstein_premise_test']
    assert sum(second)==F(einstein['discrete_laplacian_q'])
    assert second[0]-(second[1]+second[2])/2==F(einstein['discrete_linearized_G00'])
    assert not einstein['same_source_vacuum_Einstein_premise_passed']
    assert not wave['analytic_Taylor_bound_includes_floating_roundoff']
    assert not wave['finite_group_delay_identified_with_metric_ray_time']
    assert not saved['native_clocks']['clock_comparison_refutes_equivalence_principle']
    assert not metric['full_nonlinear_geometry_minimizer_asserted']
    return {'schema':'oph.independent_metric_observables_review.v1','all_pass':True,
            'independent_sparse_geometry_reproduces_all_dyadic_integers':True,
            'high_precision_native_clock_eigenvalues_and_Ramsey_probabilities_inside_certified_intervals':14,
            'independent_complete_field_evolution_cases':6,
            'independent_Christoffel_curvature_and_nonzero_vacuum_residual':True,
            'scope':'Declared static geometry and separate finite scalar field. Fixed-coordinate native clocks fail a common lapse ratio; no same-proper-device or physical Einstein/light transfer is asserted.',
            'reviewed_sha256':{str((PACKAGE/name).relative_to(ROOT)):hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest()for name in ('geometry.py','clocks.py','waves.py','PROOF.md','results.json')}}


if __name__=='__main__':
    value=run();target=Path(__file__).resolve().parents[1]/'metric_review.json'
    target.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    print('PASS: independent metric, native clocks, complete scalar field and Einstein control')
