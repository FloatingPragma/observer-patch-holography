"""Independent source-clock signs, herald coherence and rigorous trig checks."""
from fractions import Fraction
from pathlib import Path
import hashlib,json
import numpy as np
import scipy.linalg as la
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_native_beat_clock'


def finite_controls():
    rng=np.random.default_rng(54972);cases=0;max_error=0.;min_contrast=1.
    size=11;p=np.zeros((size,size));p[:3,:3]=np.eye(3)
    q=np.zeros_like(p);q[5:8,5:8]=np.eye(3)
    a=np.zeros_like(p);a[5:8,:3]=np.eye(3)
    x=a+a.T;y=1j*(a-a.T);identity=np.eye(size)
    unitary=identity-p-q+(p+q+a-a.T)/np.sqrt(2)
    assert la.norm(unitary.conj().T@unitary-identity)<1e-14
    effects=[(identity+x)/4,(identity-x)/4,(identity+y)/4,(identity-y)/4]
    assert la.norm(sum(effects)-identity)<1e-14 and min(la.eigvalsh(e) .min() for e in effects)>-1e-15
    omega=4.69;h0=np.diag([0.,0.,0.,2.,3.,omega,omega,omega,omega+1,omega+2,11.])
    radius=.2
    for coupling in (.001,.01,.03):
        for repeat in range(6):
            t=np.zeros((size,size),complex)
            for starts,targets in ((range(3),range(3,5)),(range(5,8),range(8,10))):
                for i in starts:
                    for j in targets:
                        t[i,j]=rng.normal()+1j*rng.normal();t[j,i]=t[i,j].conjugate()
            t/=la.norm(t,2);h=h0+coupling*t
            ev,v=la.eigh(h);neutral=abs(ev)<radius;charged=abs(ev-omega)<radius
            en=v[:,neutral]@v[:,neutral].conj().T;ec=v[:,charged]@v[:,charged].conj().T
            herald=en+ec;eta=coupling/radius
            # A purification over an untouched two-level reference.
            raw=rng.normal(size=(3,2))+1j*rng.normal(size=(3,2));raw/=la.norm(raw)
            state=np.zeros((size,2),complex);state[:3]=raw;state=unitary@state
            accepted=herald@state;success=la.norm(accepted)**2;accepted/=np.sqrt(success)
            assert success>=1-eta**2-1e-13
            overlap=np.vdot(state.ravel(),accepted.ravel())
            trace=2*np.sqrt(max(0.,1-abs(overlap)**2));assert trace<=2*eta+1e-12
            initial_work=np.trace(state.conj().T@h@state).real
            assert abs(initial_work-omega/2)<1e-13
            # A common flag retains cross-reference coherence; arm-resolved
            # nonselective measurement instead deletes both quadratures.
            density=accepted@accepted.conj().T
            dephased=en@density@en+ec@density@ec
            assert abs(np.trace(dephased@x))+abs(np.trace(dephased@y))<1e-13
            min_contrast=min(min_contrast,float(abs(np.trace(density@x))))
            for time in (0.,.5,2.,20.):
                evolution=(v*np.exp(-1j*time*ev))@v.conj().T
                evolved=evolution@accepted
                observed=np.array([np.trace(evolved.conj().T@o@evolved).real for o in (x,y)])
                ideal=np.array([np.cos(omega*time),-np.sin(omega*time)])
                error=float(np.max(abs(observed-ideal)));max_error=max(max_error,error)
                assert error<=min(2.,2*eta+2*coupling*time)+1e-12
                probabilities=[np.trace(evolved.conj().T@e@evolved).real for e in effects]
                assert min(probabilities)>-1e-14 and abs(sum(probabilities)-1)<1e-13
            indices=np.where(abs(v.conj().T@x@v)>1e-8)
            frequencies=[]
            for i,j in zip(*indices):
                if charged[i] and neutral[j]:frequencies.append(ev[i]-ev[j])
            assert frequencies and min(frequencies)>=omega-2*radius and max(frequencies)<=omega+2*radius
            cases+=1
    # An exact partial-isometry between occupied families commutes with N_ref.
    d=s.zeros(7);d[4,2]=d[5,3]=1
    n=s.diag(0,0,1,1,1,1,2);mu,de,dl=s.symbols('mu de dl')
    h=s.diag(0,0,mu+dl,mu+dl,mu+de,mu+de,17)
    assert n*d==d*n and h*d-d*h==(de-dl)*d
    assert s.simplify(s.trace((d+d.T)*s.diag(1,0,0,0,0,0,0)))==0
    return {'purified_finite_source_examples':cases,'waveform_comparisons':4*cases,
        'largest_observed_waveform_error':max_error,'smallest_accepted_X_contrast':min_contrast,
        'common_flag_and_arm_dephasing_distinguished':True,
        'initial_energy_work_and_reference_cost_cancellation_checked':True}


def trig_controls(receipt):
    from model.openwave_native_beat_clock.clock import trig_box
    mp.mp.dps=110;checks=0;widest=mp.mpf(0)
    def m(value):
        f=Fraction(str(value));return mp.mpf(f.numerator)/f.denominator
    values=[Fraction(0),Fraction(1,3),Fraction(-10**8,3),Fraction(10**12,7)]
    for angle in values:
        width=Fraction(1,10**45);sine,cosine=trig_box((angle-width,angle+width))
        for point in (angle-width,angle,angle+width):
            assert m(sine[0])<=mp.sin(m(point))<=m(sine[1])
            assert m(cosine[0])<=mp.cos(m(point))<=m(cosine[1]);checks+=1
    # Recover the algebraic frequencies from the native integer polynomials.
    # The rounded receipt interval midpoint is not an exact frequency and
    # cannot substitute for it in a million-time, very narrow trig test.
    local=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json').read_text())['exact_local_channels']
    x=s.Symbol('x');energies={}
    for family in ('ec','L'):
        row=local[family];polynomial=s.Poly(s.sympify(row['complete_target_annihilator']),x)
        coefficients=[m(c) for c in polynomial.all_coeffs()]
        interval=row['complete_local_energy_records'][0]['energy_interval']
        guess=sum(map(m,interval))/2
        energies[family]=mp.findroot(lambda z:mp.polyval(coefficients,z),(guess-mp.mpf('1e-6'),guess+mp.mpf('1e-6')))
    delta=sum(map(m,receipt['frequency']['source_gap_after_reference_subtraction_interval']))/2
    guess=energies['ec']-delta
    ground=mp.findroot(lambda z:(z*(z-1))**3-313*(z*(z-1))**2+31708*z*(z-1)-1041468,
                       (guess-mp.mpf('1e-6'),guess+mp.mpf('1e-6')))
    frequencies=(4+energies['ec']-ground,energies['ec']-energies['L'])
    for section,frequency,omega in ((receipt['time_rows'],receipt['frequency']['total_beat_interval'],frequencies[0]),
        (receipt['native_interchannel_beat']['time_rows'],receipt['native_interchannel_beat']['frequency_interval'],frequencies[1])):
        assert m(frequency[0])<=omega<=m(frequency[1])
        for row in section:
            angle=omega*m(row['source_time'])
            for key,value in (('ideal_sine_interval',mp.sin(angle)),('ideal_cosine_interval',mp.cos(angle))):
                lo,hi=map(m,row[key]);assert lo<=value<=hi;widest=max(widest,hi-lo);checks+=1
            error=m(row['expectation_error_upper'])
            ideals=((1+mp.cos(angle))/4,(1-mp.cos(angle))/4,(1-mp.sin(angle))/4,(1+mp.sin(angle))/4)
            for interval,ideal in zip(row['complete_four_conditional_record_intervals'],ideals):
                lo,hi=map(m,interval);assert lo<=ideal<=hi
                assert lo>=0 and hi<=mp.mpf('.5')+mp.mpf('1e-29')
    return {'high_precision_interval_checks':checks,'widest_saved_trig_interval':str(widest),
        'negative_and_trillion_scale_phase_reduction_checked':True,
        'ideal_Y_sign_checked':True,'floating_reference_is_diagnostic_only':True}


def report():
    paths=[TARGET/name for name in ('clock.py','PROOF.md','results.json')]+[Path(__file__)]
    paths.append(ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json')
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt=json.loads((TARGET/'results.json').read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in receipt['source_sha256'].items())
    result={'schema':'oph.independent_native_beat_review.v1','input_sha256':pins,
        'finite_controls':finite_controls(),'trigonometric_intervals':trig_controls(receipt),
        'scope':'A declared one-time native-pole coherent service and one common full-H herald support the proved free beat. Native microscopic pulse compilation, passive self-starting, unique exact full-H frequency and a physical Compton dictionary are not inferred.',
        'verdict':'PASS_INDEPENDENT_NATIVE_BEAT_REVIEW'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in pins.items())
    return result


if __name__=='__main__':
    result=report();(HERE/'beat_clock_review.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(result['verdict'])
