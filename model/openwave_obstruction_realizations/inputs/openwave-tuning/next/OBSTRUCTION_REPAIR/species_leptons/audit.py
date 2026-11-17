"""Independent bounded row-10/11 audit and a source-bound 3D beta repair.

Run from any directory. Writes only adjacent RESULT.json. The Wilson member
is explicitly new; the source's empirical hadron/neutral inputs are retained.
"""
from itertools import permutations, product
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
import mpmath as mp
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RER = ROOT / 'reverse-engineering-reality'
PUSH = ROOT / 'openwave-push'
PINS = {}


def consume(path):
    data = path.read_bytes()
    PINS[str(path.relative_to(ROOT))] = hashlib.sha256(data).hexdigest()
    return data


def neutral():
    old = json.loads(consume(RER / 'model/openwave_neutrino_oscillation/results.json'))
    words = json.loads(consume(RER / 'model/native_flavour_transport/results.json'))
    Y = s.Matrix(old['source_family_Y'])
    V = s.Matrix(words['primary']['transport'])
    clean = lambda M: M.applyfunc(s.simplify)
    assert clean(V.H * V - s.eye(3)) == s.zeros(3)
    M = clean(V.conjugate() * Y * V.H)
    assert clean(M - M.T) == s.zeros(3)
    y = np.array(Y.evalf(30), complex)
    v = np.array(V.evalf(30), complex)
    masses, F = np.linalg.eigh(y)
    C = F.T @ v @ F
    B = abs(C)**2
    takagi_error = np.max(abs((v@F).T @ np.array(M.evalf(30), complex) @ (v@F) - np.diag(masses)))
    assert takagi_error < 2e-13
    oldC = np.array(old['actual_charged_current_mixing_real']) + 1j*np.array(old['actual_charged_current_mixing_imag'])
    assert np.max(abs(oldC-oldC.T)) < 2e-13
    # Values are repository comparison inputs, not newly verified measurements.
    x, y23, z = .308, .470, .02215
    k = 2*np.sqrt(x*(1-x)*y23*(1-y23)*z)
    base = np.array([[(1-x)*(1-z), x*(1-z), z],
        [x*(1-y23)+(1-x)*y23*z, (1-x)*(1-y23)+x*y23*z, y23*(1-z)],
        [x*y23+(1-x)*(1-y23)*z, (1-x)*y23+x*(1-y23)*z, (1-y23)*(1-z)]])
    slope = np.array([[0,0,0],[k,-k,0],[-k,k,0]])
    permutation_rows = []
    for order in permutations(range(3)):
        d0 = base[:,order]-base[:,order].T
        d1 = slope[:,order]-slope[:,order].T
        norm = float(np.sum(d1*d1))
        c = -float(np.sum(d0*d1))/norm if norm else None
        best = min(1.,max(-1.,c)) if c is not None else 0.
        error = float(np.max(abs(d0+best*d1)))
        assert error > .01
        permutation_rows.append(dict(neutral_columns=order,required_cos_delta=c,
            minimum_max_antisymmetry_over_physical_delta=error))
    angles = dict(sin2_theta12=float(B[0,1]/(1-B[0,2])),
                  sin2_theta13=float(B[0,2]),sin2_theta23=float(B[1,2]/(1-B[0,2])))
    # Complete five-outcome instrument, reconstructed without importing producer.
    records=[]
    for momentum,time,initial in product((1.,3.,10.),(0.,1.,10.,100.),range(3)):
        E=np.hypot(momentum,masses); w=(E+momentum)/(2*E)
        K=np.diag(np.sqrt(w))@C.conj().T
        U=np.diag(np.exp(-1j*(E-E[0])*time)); out=K.conj().T@U@K[:,initial]
        p=list(abs(out)**2); success=float(np.vdot(K[:,initial],K[:,initial]).real)
        p += [1-success,success-sum(p)]
        assert min(p)>-1e-13 and abs(sum(p)-1)<2e-13
        records.append(dict(momentum=momentum,time=time,initial=initial,probabilities=p))
    return dict(single_phase_obstruction='C is complex symmetric, so |C|^2 is symmetric; excluded by the cited angle triple under all six neutral assignments.',
        comparisons=dict(sin2_theta12=x,sin2_theta13=z,sin2_theta23=y23),permutation_checks=permutation_rows,
        repair=dict(word=words['primary']['primary_word'],selection=words['predeclared_scan']['selection'],
            phase='-i',mass_coefficient='same native Y',neutral_coefficient='V* Y V-dagger',
            antisymmetry=float(np.max(abs(B-B.T))),takagi_error=float(takagi_error),angles=angles,
            masses_in_source_units=masses.tolist(),complete_records=records,
            independent_of_neutrino_comparison_inputs=True,observed_angles_matched=False,
            full_common_source_transfer_proved=False))


def lepton():
    consume(PUSH/'rowfull/ROW_11/r11lib.py')
    consume(PUSH/'rowfull/ROW_11/STATUS.md')
    mp.mp.dps=60
    b=mp.sqrt(mp.mpf(1)/2);theta=mp.mpf(2)/9
    def ratios(t,mod=1):
        masses=sorted((1+2*b*mod*mp.cos(t+2*mp.pi*k/3))**2 for k in range(3))
        return masses[1]/masses[0],masses[2]/masses[0]
    target=mp.mpf('206.7682827');sigma=mp.mpf('0.0000046')
    limit=ratios(theta)
    rows=[]
    for N in (2**16,2**20,2**24,2**28):
        zz=mp.exp(N*mp.log(1-theta+theta*mp.exp(1j/mp.mpf(N))))
        vals=ratios(mp.arg(zz),abs(zz))
        rows.append(dict(rounds=N,muon_ratio=mp.nstr(vals[0],35),tau_ratio=mp.nstr(vals[1],35),
            muon_residual_sigma=mp.nstr((vals[0]-target)/sigma,20)))
    derivative=mp.diff(lambda t:mp.log(ratios(t)[0]),theta)
    corrected=mp.findroot(lambda t:ratios(t)[0]-target,theta)
    return dict(limit_muon_ratio=mp.nstr(limit[0],40),limit_tau_ratio=mp.nstr(limit[1],40),
        muon_residual_sigma=mp.nstr((limit[0]-target)/sigma,30),rounds=rows,
        logarithmic_phase_sensitivity=mp.nstr(derivative,30),
        diagnostic_target_fitted_phase_shift=mp.nstr(corrected-theta,30),
        fitted_phase_is_not_a_solution=True,
        next_physical_obligation='Compute a mass-shell correction from the fixed interacting source; do not select phase or round count using the muon target.',
        solved=False)


def beta():
    path=PUSH/'rowfull/ROW_25/consumed/source_arrays.npz'
    consume(path); A=np.load(path)
    # All matrices below are the consumed 120-state source, not a fitted angular polynomial.
    W=A['W'];local=A['local'];P=A['pauli'];g=float(A['ga'])
    projected=np.array([[W[:,2:].conj().T@local[x,mu,0]@W[:,0] for mu in range(4)] for x in range(2)])
    current_error=max(np.max(abs(projected[:,0]-np.array([.5,0]))),
        max(np.max(abs(projected[:,mu]-.5*g*P[mu][:,0])) for mu in range(1,4)))
    assert current_error<2e-10
    Q=float(A['En']-A['Ep']);G=1/64;t=1/100000
    assert 1<Q<np.sqrt(8)
    assert 1.2752<g<1.2754 and np.max(A['neutral'])<1e-6
    rows=[]
    for L,w,chirality in product((4,8),(0.,1.),(1.,-1.,0.)):
        th=np.array(list(product(2*np.pi*np.arange(L)/L,repeat=3)))
        p=2*np.sin(th);p[abs(p)<1e-14]=0
        N=len(th);mass_e=1+w*np.sum(1-np.cos(th),axis=1)
        Ee=np.sqrt((p*p).sum(axis=1)+mass_e**2); ae=p/Ee[:,None]
        ve=(p*2*np.cos(th)+mass_e[:,None]*w*np.sin(th))/Ee[:,None]
        ve[abs(ve)<1e-13]=0
        form=np.cos((th[:,0,None]-th[None,:,0])/2)**2
        probability=np.zeros((N,N));vh_last=None
        for j,mu in enumerate(A['neutral']):
            mass_n=mu+w*np.sum(1-np.cos(th),axis=1)
            En=np.sqrt((p*p).sum(axis=1)+mass_n**2); bn=-p/En[:,None]
            vh=(p*2*np.cos(th)+mass_n[:,None]*w*np.sin(th))/En[:,None]
            vh[abs(vh)<1e-13]=0;vh_last=vh
            spin=(1+3*g*g+(1-g*g)*(ae@bn.T)+chirality*(2*g*(1-g)*ae[:,2,None]+2*g*(1+g)*bn[None,:,2]))/2
            assert spin.min()>-1e-12
            mismatch=Ee[:,None]+En[None,:]-Q
            probability+=(G*t/N)**2*abs(A['R'][0,j])**2*form*spin*np.sinc(mismatch*t/(2*np.pi))**2
        pe=probability.sum(axis=1);pn=probability.sum(axis=0);total=float(pe.sum())
        odd_e=float(pe@np.sign(ve[:,2]));odd_n=float(pn@np.sign(-vh_last[:,2]))
        delta=(54*G*t)**2/2; remainder=(8*G*t+delta)*delta
        # Complete original-input instrument: 1/120 preparation success, selected
        # event and its full Fock complement; this is not conditioning away failures.
        gt=F(1,64)*F(1,100000); tr=F(1,100000); glo=F('1.2752')
        exact_delta=(54*gt)**2/2
        exact_remainder=(8*gt+exact_delta)*exact_delta
        weight_error=gt**2*16*(14*tr)**2/12
        # A pair of staggered z momenta has nonnegative directional sum.
        # The z=+/-pi/2 planes alone give S/N >= 4/(7L).
        electron_lower=gt**2*glo*(glo-1)*F(2,7*L)
        neutral_lower=gt**2*glo*(glo+1)*F(2,7*L)
        e_margin=electron_lower-weight_error-exact_remainder
        n_margin=neutral_lower-weight_error-exact_remainder
        assert e_margin>0 and n_margin>0
        rows.append(dict(L=L,w=w,chirality=chirality,source_time=str(tr),selected_probability=total,
            electron_velocity_odd=odd_e,antineutrino_velocity_odd=odd_n,
            electron_conditional_velocity_asymmetry=odd_e/total,
            antineutrino_conditional_velocity_asymmetry=odd_n/total,
            full_Fock_Dyson_probability_error_upper=remainder,
            exact_Wilson_left_electron_negative_margin=str(e_margin/120) if w==1 and chirality==1 else None,
            exact_Wilson_left_antineutrino_positive_margin=str(n_margin/120) if w==1 and chirality==1 else None,
            electron_margin_over_omitted_sectors=abs(odd_e)-remainder,
            original_failed_preparation=119/120,
            original_selected_interval=[max(0,total-remainder)/120,min(1,total+remainder)/120],
            original_prepared_complement_interval=[max(0,1-total-remainder)/120,min(1,1-total+remainder)/120]))
        if w==0 or chirality==0:
            assert abs(odd_e)<1e-25 and abs(odd_n)<1e-25
        else:
            assert chirality*odd_e<0 and chirality*odd_n>0
            assert abs(odd_e)>remainder
    return dict(declared_Wilson_coefficient=1,source_axial=g,source_Q=Q,
        current_projection_residual=float(current_error),records=rows,
        three_dimensional_onshell_doubler_bounds=dict(electron_wrong_halfzone_min=np.sqrt(8),
            neutral_wrong_halfzone_plus_electron_min=3.,physical_Q=Q),
        meaning='3D group-velocity parity signal repaired on the Wilson-amended declared source at a controlled short time; strict signs have a rational analytic bound independent of the floating sums. Numerical probability intervals still carry unbounded floating roundoff and are diagnostics.',
        inherited_inputs=['proton/neutron/pion masses','axial coupling','neutral masses/mixing','electron unit'],
        missing=['source-derived neutron-proton gap','source-derived chirality and axial value','physical lifetime','transfer to shared common native source'],
        full_row_solved=False)


if __name__=='__main__':
    result=dict(schema='oph.obstruction_repair.species_leptons.v1',row10=neutral(),row11=lepton(),row25=beta())
    for p,h in PINS.items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h, 'source changed during run: '+p
    result['input_sha256']=PINS
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'row10_angles':result['row10']['repair']['angles'],
        'row11_muon_sigma':result['row11']['muon_residual_sigma'],
        'row25_records':len(result['row25']['records']),
        'row25_wilson':[r for r in result['row25']['records'] if r['w']==1 and r['chirality']==1]},indent=2))
