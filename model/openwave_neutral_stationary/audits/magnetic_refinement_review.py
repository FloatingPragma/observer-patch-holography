"""Independent spectral, analytic-tail and finite-time magnetic controls."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import itertools
import math
import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp
import mpmath as mp
from model.openwave_magnetic_refinement import certified_probe as source
from model.openwave_intrinsic_magnetic.solver import dirac,GAMMA,BETA,SIGMA
from model.openwave_intrinsic_magnetic.solver import geometry
from model.openwave_larmor_dynamics.dynamics import field_minimizer
from model.openwave_field_observables.lattice import curl,curl_adjoint

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_magnetic_refinement'


def exact_number(x):
    q=F(float(x));return mp.mpf(q.numerator)/q.denominator


def analytic_residual(vector,energy,n,a,B):
    """Direct analytic Peierls matrix, without the producer Taylor matvec."""
    with mp.workdps(95):
        a=mp.mpf(a.numerator)/a.denominator;B=mp.mpf(B.numerator)/B.denominator
        points=list(itertools.product(range(n),repeat=3));index={x:i for i,x in enumerate(points)}
        nn=len(points);shift=(n-1)//2
        values=[[mp.mpc(exact_number(vector[j*nn+i].real),exact_number(vector[j*nn+i].imag))
                 for i in range(nn)] for j in range(4)]
        e=exact_number(energy);res=[[mp.mpc(0) for _ in points] for _ in range(4)]
        for i,p in enumerate(points):
            xyz=[x-shift for x in p];mass=8+sum(x*x for x in xyz)*a*a/100+mp.mpf(3)/(2*a)
            for spin in range(4):res[spin][i]=((1 if spin<2 else -1)*mass-e)*values[spin][i]
            for axis in range(3):
                for step in (-1,1):
                    q=list(p);q[axis]+=step
                    if not 0<=q[axis]<n:continue
                    j=index[tuple(q)]
                    angle=B*a*a*(-xyz[1] if axis==0 else xyz[0] if axis==1 else 0)/2
                    phase=mp.exp(step*mp.j*angle)
                    for left in range(4):
                        for right in range(4):
                            gamma=complex(GAMMA[axis][left,right]);beta=complex(BETA[left,right])
                            if gamma==beta==0:continue
                            coefficient=(step*mp.j*mp.mpc(gamma)-mp.mpc(beta)/2)*phase/(2*a)
                            res[left][i]+=coefficient*values[right][j]
        norm2=sum(abs(z)**2 for row in values for z in row)
        raw=mp.sqrt(sum(abs(z)**2 for row in res for z in row))
        return {'raw':mp.nstr(raw,45),'normalized':mp.nstr(raw/mp.sqrt(norm2),45),
                'norm_squared':mp.nstr(norm2,45)}


def sturm_count(diagonal,off,x):
    # Sign variation of leading characteristic minors, independent of LDL pivots.
    minors=[F(1),diagonal[0]-x]
    for d in diagonal[1:]:minors.append((d-x)*minors[-1]-off*off*minors[-2])
    assert all(minors)
    return sum(a*b<0 for a,b in zip(minors,minors[1:]))


def run():
    receipt=json.loads((TARGET/'results.json').read_text())['same_refined_source_precession']
    n=17;a=F(1,3);B=F(1,100)
    diagonal=[F(1,2)/a+F(1,100)*(a*j)**2 for j in range(-8,9)];off=-F(1,2)/(2*a)
    for name in ('zero_probe_all_mode_gap','probe_all_mode_gap'):
        for j,interval in enumerate(receipt[name]['one_dimensional_lowest_two_exact_Sturm_intervals']):
            lo,hi=map(F,interval)
            assert sturm_count(diagonal,off,lo)==j and sturm_count(diagonal,off,hi)==j+1
    analytic=[]
    for filename,field in (('probe_zero.npz',F(0)),('probe_B001.npz',B)):
        with np.load(TARGET/'witnesses'/filename,allow_pickle=False) as z:
            for j in range(1 if not field else 2):
                independent=analytic_residual(z['vectors'][:,j],z['energies'][j],17,a,field)
                cert=receipt['zero_probe_exact_matrix_residual'] if not field else receipt['probe_exact_matrix_residuals'][j]
                r=F(cert['analytic_matrix_residual_upper'])
                with mp.workdps(95):
                    assert mp.mpf(independent['normalized'])<mp.mpf(r.numerator)/r.denominator
                analytic.append({'file':filename,'line':j,**independent,'certified_upper':str(r)})
    tails=[]
    with mp.workdps(100):
        for q in (F(1,225),F(-1,225),F(1,20),F(2,5)):
            z=mp.mpf(q.numerator)/q.denominator
            polynomial=sum((mp.j*z)**k/mp.factorial(k) for k in range(19))
            error=abs(mp.exp(mp.j*z)-polynomial);upper=2*abs(z)**19/mp.factorial(19)
            assert error<=upper
            tails.append({'angle':str(q),'actual_tail':mp.nstr(error,30),'analytic_upper':mp.nstr(upper,30)})
    dense=[]
    chirality=GAMMA[0]@GAMMA[1]@GAMMA[2]@BETA
    theta=la.block_diag(1j*SIGMA[1],1j*SIGMA[1]);I=np.eye(27)
    for field in (0.,.01):
        h,_=dirac(3,1.,9,np.zeros((9,9,9,3)),external_B=field);h=h.toarray()
        ch=np.kron(chirality,I);anti=np.kron(theta,I)
        assert la.norm(ch@h+h@ch)<1e-12
        if field==0:assert la.norm(anti@h.conj()@anti.conj().T-h)<1e-12
        eig=la.eigvalsh(h);absolute=np.sort(abs(eig));bound=source.scalar_mass_gap(3,F(1),F(str(field)))
        assert float(F(bound['fifth_absolute_D_eigenvalue_lower']))<=absolute[4]+1e-12
        if field==0:assert np.ptp(absolute[:4])<1e-12
        dense.append({'external_B':field,'first_five_absolute_levels':absolute[:5].tolist(),
                      'fifth_level_lower':bound['fifth_absolute_D_eigenvalue_lower']})
    # Independent nonlinear unitary example verifies the dimension-free
    # Duhamel inequality without replacing it by an exponential estimate.
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.]);h0=np.diag([0.,.7])
    initial=np.array([math.sqrt(.8),1j*math.sqrt(.2)]);delta=.03
    def rhs(t,psi):return -1j*(h0+delta*np.vdot(psi,sz@psi).real*sx)@psi
    times=np.linspace(0,8,33);solution=solve_ivp(rhs,(0,8),initial,t_eval=times,rtol=1e-11,atol=1e-12)
    maximum=0.
    for t,psi in zip(times,solution.y.T):
        exact=la.expm(-1j*h0*t)@initial;error=la.norm(psi-exact)
        assert abs(la.norm(psi)-1)<1e-9 and error<=delta*t+1e-10
        maximum=max(maximum,float(error))
    # Uniform rational bounds imply the phase normalization and positive enclosure.
    prop=F(receipt['complete_state_error_upper_at_duration']);spin=F(receipt['uniform_spin_phase_certificate']['uniform_frozen_transverse_spin_radius_lower'])
    assert prop<F(1,1000) and spin>F(49,100)
    # For the possibly nonunit normalized two-frequency vector, S+ has norm1.
    assert (2+prop)*prop/spin<F(1,2)
    assert F(4,3)*(2+prop)*prop/spin<4*prop/spin
    lo,hi=map(F,receipt['finite_g_interval_exact']);assert 0<lo<hi<2
    # Same weak-coupling member: reconstruct the stationary spin line and
    # current bond-by-bond, then reject zero field after scaling by kappa.
    with np.load(TARGET/'witnesses/probe_zero.npz',allow_pickle=False) as z:
        q,_=la.qr(z['vectors'][:,:2],mode='economic')
    nn=17**3;sz=np.repeat((.5,-.5,.5,-.5),nn)
    _,spinors=la.eigh(q.conj().T@(sz[:,None]*q));psi=q@spinors[:,-1];psi/=la.norm(psi)
    _,positions,edges=geometry(17,1/3,63);parts=psi.reshape(4,-1)
    current=np.zeros((63,63,63,3));dipole=np.zeros(3)
    for left,right,axis,address in edges:
        derivative=-(GAMMA[axis]+.5j*BETA)/2
        value=-2*np.vdot(parts[:,left],derivative@parts[:,right]).real
        current[address+(axis,)]=value
        dipole+=np.cross((positions[left]+positions[right])/2,np.eye(3)[axis])*value/2
    kappa=1e-14;field=field_minimizer(current,1/3,kappa)
    scaled=la.norm((curl_adjoint(curl(field,1/3),1/3)-kappa*current/(1/3)**3).reshape(-1))/kappa
    mutation=la.norm(current.reshape(-1)/(1/3)**3)
    stationary=receipt['same_member_intrinsic_stationary_witness']
    assert scaled<1e-8 and mutation>.03 and la.norm(field.reshape(-1))>0
    assert np.max(abs(dipole-stationary['actual_intrinsic_current_moment']))<1e-12
    paths=set(TARGET.glob('*.py'))|set(TARGET.glob('*.md'))|{TARGET/'results.json'}|set((TARGET/'witnesses').glob('*.npz'))
    paths.update(ROOT/name for name in json.loads((TARGET/'results.json').read_text())['input_sha256'])
    return {'schema':'oph.magnetic_refinement.independent_review.v1',
        'full_19652_row_analytic_residuals':analytic,'analytic_precision_digits':95,
        'independent_exact_characteristic_minor_Sturm_counts':8,'phase_Taylor_controls':tails,
        'full_small_signed_kernel_controls':dense,
        'nonlinear_Duhamel_control_times':33,'nonlinear_Duhamel_maximum_norm_error':maximum,
        'same_member_independent_static_current_moment':dipole.tolist(),
        'same_member_scaled_Maxwell_residual':float(scaled),
        'same_member_zero_field_mutation_scaled_residual':float(mutation),
        'finite_g_interval':[str(lo),str(hi)],
        'source_scope':'Declared finite positive Dirac kernel, confinement and unique quasistatic mean field. Static residuals do not identify a global matter minimizer; the dynamic g interval is finite-duration, with supplied mass/coupling and no continuum or arbitrary-reference nonlinear channel claim.',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        'verdict':'PASS'}


if __name__=='__main__':
    result=run();(HERE/'magnetic_refinement_review.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['verdict'])
