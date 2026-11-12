"""Current normalization, full three-body phase space and finite-time records."""
from fractions import Fraction as F
from functools import lru_cache
import math
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from model.common_source import calibrated_matter as native
from model.openwave_weak_decay.decay import mass_data
from model.lepton_mass_compiler.compiler import electron_calibration

G=F(1,64)
B0=F(1,2000)
B1=F(1,250)
HUGE_BIN_WIDTH=F(1,2**60)
SIGMA=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
BAR=(np.eye(2),*(-z for z in SIGMA))
ETA=(1,-1,-1,-1)

@lru_cache(maxsize=1)
def source_data():
    mass,heavy,boxes=mass_data()
    lo=F(boxes[1][0])/F(boxes[2][1]);hi=F(boxes[1][1])/F(boxes[2][0])
    assert F(3,625)<lo<hi<F(49,10000)
    J=native.weak_raising_current(families=False)
    ell=native._LH_tensor()*s.Matrix([0,1])
    neutral=next(i for i in range(15)if ell[i]!=0)
    D=native._oriented_components('charged_lepton')[1]
    charged=next(iter(D.todok()))[0]
    generators=((J+J.H)/2,-s.I*(J-J.H)/2,(J*J.H-J.H*J)/2)
    # The two ordered family-current products supply this factor of two.
    coefficient=s.simplify(2*sum(T[charged,neutral]*T[neutral,charged]for T in generators))
    assert coefficient==1
    return dict(mass=mass,heavy=heavy,interval=[str(lo),str(hi)],coefficient=str(coefficient),
                neutral=neutral,charged=charged,source_effect_intervals=boxes)

def left_density(p):
    return p[0]*np.eye(2)-sum(p[j+1]*SIGMA[j]for j in range(3))

def minkowski(a,b):return float(a[0]*b[0]-np.dot(a[1:],b[1:]))

def kinematics(m,z,costheta,phi):
    """Rest parent M=1; z is the electron/parent-neutrino invariant mass squared."""
    q=(1-z)/2;pair_energy=(1+z)/2;beta=q/pair_energy;gamma=pair_energy/np.sqrt(z)
    k=(z-m*m)/(2*np.sqrt(z));estar=(z+m*m)/(2*np.sqrt(z));st=np.sqrt(1-costheta**2)
    pe=np.array([gamma*(estar-beta*k*costheta),k*st*np.cos(phi),k*st*np.sin(phi),gamma*(k*costheta-beta*estar)])
    pn=np.array([gamma*(k+beta*k*costheta),-k*st*np.cos(phi),-k*st*np.sin(phi),gamma*(-k*costheta-beta*k)])
    pa=np.array([q,0.,0.,q]);P=np.array([1.,0.,0.,0.])
    assert np.max(abs(pe+pn+pa-P))<1e-13
    return P,pe,pn,pa

def spin_trace(P,pe,pn,pa,parent_average=True):
    de,dn,da=map(left_density,(pe,pn,pa));dp=left_density(P)/(2 if parent_average else 1)
    value=sum(ETA[mu]*ETA[nu]*np.trace(de@BAR[mu]@da@BAR[nu])*
              np.trace(dn@BAR[mu]@dp@BAR[nu])for mu in range(4)for nu in range(4))
    assert abs(value.imag)<1e-12
    return float(value.real)

def spin_certificate():
    m=source_data()['mass'];rows=[]
    for z in(.07,.4,.91):
        for ct,phi in((-.73,.1),(.3,1.4),(.86,2.7)):
            P,pe,pn,pa=kinematics(m,z,ct,phi)
            direct=spin_trace(P,pe,pn,pa);invariant=8*minkowski(P,pa)*minkowski(pe,pn)
            assert abs(direct-invariant)<2e-14
            assert abs(direct-2*(1-z)*(z-m*m))<2e-14
            doubled=spin_trace(P,pe,pn,pa,False)
            assert abs(doubled-2*direct)<2e-14
            rows.append(dict(pair_s=z,costheta=ct,phi=phi,trace=direct,invariant=invariant,
                             omitted_parent_spin_average=doubled))
    return rows

@lru_cache(maxsize=1)
def exact_phase_space():
    z,a=s.symbols('z a',positive=True)
    integrand=s.expand((1-z)**2*(z-a)**2/z)
    antiderivative=s.integrate(integrand,z)
    integral=s.simplify(antiderivative.subs(z,1)-antiderivative.subs(z,a))
    f=1-8*a+8*a**3-a**4-12*a*a*s.log(a)
    assert s.simplify(12*integral-f)==0
    assert s.integrate(z*(1-z)**2,(z,0,1))==s.Rational(1,12)
    return {'phase_space_dPhi3_ds':'(M^2-s)(s-m^2)/(128*pi^3*M^2*s)',
            'unpolarized_native_amplitude_squared':'2*G^2*(M^2-s)*(s-m^2)',
            'integral_antiderivative':str(antiderivative),'integral':str(integral),
            'rate':'G^2*M^5*f(m/M)/(1536*pi^3)',
            'f':'1-8*r^2+8*r^6-r^8-12*r^4*log(r^2)',
            'parent_spin_average':'1/2','massless_neutral_helicities_each':1,
            'identical_particle_factor':'1: distinct neutral and antineutral families',
            'current_cross_family_coefficient':source_data()['coefficient']}

_X,_W=leggauss(24)
def phase_factor(r):
    if r<0 or r>1:raise ValueError('mass ratio must be in [0,1]')
    if r==0:return 1.
    if r==1:return 0.
    a=r*r
    if a<.5:return 1-8*a+8*a**3-a**4-12*a*a*np.log(a)
    w=(_X+1)/2
    return float(12*(1-a)**5*np.sum((_W/2)*w*w*(1-w)**2/(a+(1-a)*w)))

def density(E,m=None):
    if m is None:m=source_data()['mass']
    if E<=m or E>=2:return 0.
    return E**5*phase_factor(m/E)*(2-E)**4/(3072*np.pi**4)

def kernel_bounds():
    """Exact rational derivative and Fourier-tail bounds, uniform in native mass."""
    m=F(49,10000);lo=F(3,625);D=F(3072*3**4)
    # e^7 exceeds 2/m_lo using a positive rational partial sum.
    assert sum((F(7)**k/F(math.factorial(k))for k in range(16)),F(0))>2/lo
    S0=F(32,3)
    S1=32+64*m*m+9*m**5+384*m**4
    S2=80+96*m*m+179*m**4
    S3=160+96*m*m+52*m**3
    A0=F(256,315)/D
    A2=(16*S2+64*S1+48*S0)/D
    A3=(16*S3+96*S2+144*S1+48*S0)/D
    assert 4*A0*A2<B0**2
    assert F(27,8)*A0*A3*A3<B1**3
    # Positive integral gives f <= 1. Negative terms alone give the lower bound.
    f_lower=1-8*m*m-m**8
    rate_lower=f_lower/(1536*F(22,7)**3)
    assert rate_lower>F(1,50000)
    return {'mass_interval':[str(lo),str(m)],'log_2_over_mass_upper':'7',
            'S_derivative_L1_bounds':[str(v)for v in(S0,S1,S2,S3)],
            'rho_L1_upper':str(A0),'rho_second_derivative_L1_upper':str(A2),
            'rho_third_derivative_L1_upper':str(A3),
            'correlation_L1_upper':str(B0),'time_weighted_correlation_L1_upper':str(B1),
            'Gamma_over_G_squared_lower':str(rate_lower),
            'B0_squared_margin':str(B0**2-4*A0*A2),
            'B1_cubed_margin':str(B1**3-F(27,8)*A0*A3*A3)}

def survival_bounds():
    bounds=kernel_bounds();m=source_data()['mass'];rate=float(G*G)*phase_factor(m)/(1536*np.pi**3)
    rows=[]
    for k in(1,2,3,4):
        delta=G*G*B1*(1+25*k)
        # t < k*50000/G^2. A vector error t*h/2 gives probability error <= t*h.
        finite_probability=F(k*50000)*HUGE_BIN_WIDTH/(G*G)
        error=2*delta+finite_probability;center=float(np.exp(-k))
        rows.append({'nominal_lifetimes':k,'source_time_diagnostic':k/rate,
                     'continuum_survival_center':center,'continuum_amplitude_error_upper':str(delta),
                     'finite_bin_survival_error_upper':str(finite_probability),
                     'total_survival_probability_error_upper':str(error),
                     'survival_interval_diagnostic':[center-float(error),center+float(error)],
                     'inclusive_emission_interval_diagnostic':[1-center-float(error),1-center+float(error)]})
    bins=math.ceil(float(2-m)/float(HUGE_BIN_WIDTH))
    return {'native_G':'1/64 in middle-mass units; dimension energy^(-2)',
            'Gamma_in_middle_mass_units':rate,'f_native':phase_factor(m),
            'reservoir_uv_amplitude':'u(E)=(2-E)^2 for m<=E<=2, zero outside; M=1',
            'reservoir_density':'rho(E)=E^5*f(m/E)*(2-E)^4/(3072*pi^4)',
            'comparison_amplitude_bound':'G^2*B1*(1+G^2*B0*t)',
            'finite_bin_width_upper':str(HUGE_BIN_WIDTH),
            'finite_bin_count_diagnostic':bins,'large_bin_model_executed':False,
            'rows':rows,'rational_certificate':bounds}

def bin_data(n=64,m=None):
    if n<2:raise ValueError('at least two energy bins are required')
    if m is None:m=source_data()['mass']
    edges=np.linspace(m,2,n+1);centers=(edges[:-1]+edges[1:])/2
    weights=np.array([quad(density,a,b,args=(m,),epsabs=1e-20,epsrel=1e-11)[0]
                      for a,b in zip(edges[:-1],edges[1:])])
    assert np.all(weights>0)
    return centers,weights

def finite_records(n=64,time=1000.,coupling=float(G)):
    energies,weights=bin_data(n);H=np.diag(np.r_[1.,energies]);H[0,1:]=H[1:,0]=coupling*np.sqrt(weights)
    eigenvalues,U=np.linalg.eigh(H);state=U@(np.exp(-1j*time*eigenvalues)*U[0].conj())
    probabilities=abs(state)**2;complete=float(probabilities.sum());assert abs(complete-1)<2e-13
    full_energy=float(np.vdot(state,H@state).real)
    dephased_energy=float(np.dot(probabilities,np.r_[1.,energies]));work=dephased_energy-1
    interference=2*coupling*np.real(np.conj(state[0])*np.dot(np.sqrt(weights),state[1:]))
    assert abs(full_energy-1)<2e-13 and abs(work+interference)<2e-13
    assert abs(work)<=abs(coupling)*np.sqrt(weights.sum())+2e-13
    return {'energy_bins':n,'time':time,'coupling':coupling,'parent_survival':float(probabilities[0]),
            'all_emission_bin_probabilities':probabilities[1:].tolist(),'inclusive_emission':float(probabilities[1:].sum()),
            'bin_energies':energies.tolist(),'integrated_bin_spectral_weights':weights.tolist(),
            'complete_probability_sum':complete,'full_energy_before_readout':full_energy,
            'mean_energy_after_complete_record_readout':dephased_energy,'readout_energy_work':work,
            'energy_work_upper':float(abs(coupling)*np.sqrt(weights.sum())),
            'numerical_integration_scope':'double-precision small-carrier diagnostic; rational large-carrier theorem is separate'}

def certificate():
    source=source_data();spin=spin_certificate();phase=exact_phase_space();survival=survival_bounds()
    rows=[finite_records(n,t)for n in(32,64,128)for t in(1.,1000.,10000.)]
    zero=finite_records(32,1000.,0.);assert zero['inclusive_emission']==0
    calibration=electron_calibration();M=float(calibration['rows'][1]['electron_anchored_GeV'])
    rate=survival['Gamma_in_middle_mass_units']
    return {'schema':'oph.inclusive_positive_frequency_native_weak_decay.v1',
            'architecture':'full positive-frequency native charged-current middle-to-light plus two distinct massless neutrals, its Hermitian reverse, smooth UV amplitude, all spin/angular/energy channels; free spectators retained',
            'full_fock_invariant_sector':'N_parent+N_light=N_parent+N_neutral=N_parent+N_antineutral=1; one parent or three daughters, no selected two-state replacement',
            'reference':'daughter vacuum and all complementary native free quasiparticles vacant relative to the filled free sea; declared reference, not an interacting-ground claim',
            'source':source,'native_current_spin_traces':spin,'exact_phase_space':phase,
            'controlled_finite_time_survival':survival,'executed_finite_complete_records':rows,'G_zero_control':zero,
            'conditional_units':{'electron_anchor_is_empirical':True,'native_phase_2_over_9_is_architecture':True,
                'middle_mass_GeV':M,'width_GeV':rate*M,'G_GeV_inverse_squared':float(G)/(M*M),
                'time_unit':'hbar / middle_mass; no observed weak clock or muon lifetime supplied',
                'calibration_source':calibration},
            'readout_work_uniform_upper':'G*sqrt(256/(315*3072*3^4)) in middle-mass units',
            'arbitrary_reference_domain':'unitary and complete record channel on the exact admitted invariant sector, with spectator identity; the norm finite-bin bound is stable under tensoring any external reference',
            'scope':'positive-frequency contact member, not equivalence to the unfiltered270-CAR contact; finite-bin theorem and continuum comparison do not assert an irreversible finite unitary, interacting vacuum, observed lifetime, or electroweak parameter fit',
            'verdict':'PASS_INCLUSIVE_NATIVE_CURRENT_DECAY_WITH_FINITE_SURVIVAL_BOUNDS'}
