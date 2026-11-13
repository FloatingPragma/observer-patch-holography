"""Native forward-current contraction and a finite complete medium instrument."""
from functools import lru_cache
from itertools import product
from fractions import Fraction as F
import numpy as np
import sympy as sy
from scipy.linalg import expm
from model.common_source import calibrated_matter as internal
from model.openwave_neutrino_calibration.calibration import calibration
from model.openwave_baryon_weak_formfactors.native import PAULI
from model.openwave_lepton_widths.widths import units

@lru_cache(None)
def native_forward():
    # Every coefficient is a contraction of the literal fifteen-state current.
    raising=np.array(internal.weak_raising_current(families=False)).astype(complex)
    ts=((raising+raising.conj().T)/2,(raising-raising.conj().T)/(2j),(raising@raising.conj().T-raising.conj().T@raising)/2)
    small=[t[np.ix_([10,11],[10,11])] for t in ts]
    J=[(mu,np.kron(PAULI[mu],np.kron(t,np.eye(3)))) for mu in range(4) for t in small]
    def vec(spin,weak,family):
        v=np.zeros(12,complex);v[(spin*2+weak)*3+family]=1;return v
    nu=[vec(s,0,f) for f,s in product(range(3),range(2))]
    electron=[vec(s,1,0) for s in range(2)]
    M=np.zeros((12,12),complex)
    # Two ordered current products, fermionic exchange, and the squared
    # 1/sqrt(2) active weight of a rest electron give total prefactor one.
    for i,s,j,r in product(range(6),range(2),range(6),range(2)):
        M[2*i+s,2*j+r]=sum((1 if mu==0 else -1)*((nu[i].conj()@a@nu[j])*(electron[s].conj()@a@electron[r])-(nu[i].conj()@a@electron[r])*(electron[s].conj()@a@nu[j])) for mu,a in J)
    exact=sy.Matrix(M).applyfunc(sy.nsimplify);x=sy.Symbol('x')
    assert exact==exact.H and sy.factor(exact.charpoly(x).as_expr())==x**9*(x-1)*(x+1)**2
    mean=np.einsum('isjs->ij',M.reshape(6,2,6,2))/2
    expected=np.diag([F(1,4)]*2+[F(-1,4)]*4).astype(float)
    assert np.max(abs(mean-expected))<1e-14
    return {'matrix':M,'mean':mean,'exact_polynomial':str(x**9*(x-1)*(x+1)**2),'exact_matrix':[[str(v) for v in row] for row in exact.tolist()]}


def source(momentum_eV,antineutrino=False):
    c=calibration();m=c['physical_masses_eV'];p=float(momentum_eV);E=np.hypot(p,m)
    # Stable chirality weights, including the suppressed helicity.
    small=m/np.sqrt(2*E*(E+p));large=np.sqrt((E+p)/(2*E));u=np.column_stack((small,large)).reshape(6)
    v=np.column_stack((large,small)).reshape(6);R=c['R']
    Q=np.kron(R,np.eye(2))@np.diag(u);D=np.kron(Q,np.eye(2));M=D.conj().T@native_forward()['matrix']@D
    mean=np.einsum('isjs->ij',M.reshape(6,2,6,2))/2
    if antineutrino:M=-M.conj();mean=-mean.conj();R=R.conj()
    gaps=(m*m-m[0]**2)/(E+E[0]);free=np.repeat(gaps,2)
    assert np.linalg.norm(M-M.conj().T)<1e-14 and np.linalg.norm(mean-mean.conj().T)<1e-14
    assert np.linalg.norm(M,2)<=1+1e-14 and np.linalg.norm(mean,2)<=.25+1e-14
    return {'M':M,'mean':mean,'free':free,'u':u,'v':v,'R':R,'masses':m,'p':p,'Q':Q}


def unit_dictionary(rho_g_cm3=2.8,electron_fraction=.5,length_km=1300.):
    # The mass-density scenario uses one gram per mole of nucleon number as
    # a declared material convention. Density/composition are apparatus inputs.
    NA=F('6.02214076e23');hc=F('1.973269804e-5')
    ne=float(F(str(rho_g_cm3))*F(str(electron_fraction))*NA)
    density=ne*float(hc)**3;time=length_km*1e5/float(hc)
    G=units()['native_Weyl_G_GeV_minus2']*1e-18
    return {'electron_density_cm_minus3':ne,'electron_density_eV3':density,'length_km':length_km,'time_eV_inverse':time,'Weyl_G_eV_minus2':G,'charged_current_potential_eV':G*density/2,'native_common_neutral_current_shift_eV':-G*density/4,'integrated_native_contact':G*density*time,'hbar_c_eV_cm':float(hc),'density_g_cm3':rho_g_cm3,'electron_fraction':electron_fraction,'density_role':'declared homogeneous electron-medium scenario; not an Earth density fit'}


def complete_record(rho,d,initial_flavour=1):
    active=[]
    for f in range(3):
        K=np.zeros((2,6),complex)
        for i,s in product(range(3),range(2)):K[s,2*i+s]=d['R'][f,i]*d['u'][2*i+s]
        active.append(float(np.trace(K@rho@K.conj().T).real))
    sterile=float(np.sum(np.diag(rho).real*d['v']**2))
    failure=float(sum(abs(d['R'][initial_flavour,i])**2*d['v'][2*i+1]**2 for i in range(3)))
    assert min(active+[sterile,failure])>-1e-12
    assert abs(sum(active)+sterile+failure-1)<5e-10
    return {'e':active[0],'mu':active[1],'tau':active[2],'sterile_quantum_output':sterile,'failed_preparation':failure}


def input_state(d,flavour=1):
    psi=np.zeros(6,complex)
    for i in range(3):psi[2*i+1]=d['R'][flavour,i].conjugate()*d['u'][2*i+1]
    return np.outer(psi,psi.conj())


def collision_channel(d,free_phase,contact_phase,N):
    H=np.kron(np.diag(free_phase),np.eye(2))+contact_phase*d['M']
    W=expm(-1j*H/N).reshape(6,2,6,2)
    Ks=[W[:,s,:,r]/np.sqrt(2) for s,r in product(range(2),repeat=2)]
    assert np.linalg.norm(sum(K.conj().T@K for K in Ks)-np.eye(6))<2e-14
    S=sum(np.kron(K.conj(),K) for K in Ks)
    return np.linalg.matrix_power(S,N),Ks


def run(momentum_GeV,length_km,rho_g_cm3=2.8,N=2**14,antineutrino=False):
    d=source(momentum_GeV*1e9,antineutrino);u=unit_dictionary(rho_g_cm3,.5,length_km)
    phase=d['free']*u['time_eV_inverse'];alpha=u['integrated_native_contact'];initial=input_state(d)
    C,Ks=collision_channel(d,phase,alpha,N);rho=(C@initial.reshape(36,order='F')).reshape(6,6,order='F')
    Hbar=np.diag(phase)+alpha*d['mean'];U=expm(-1j*Hbar);mean=U@initial@U.conj().T
    uv=expm(-1j*np.diag(phase));vacuum=uv@initial@uv.conj().T
    records=complete_record(rho,d);average=complete_record(mean,d);vac=complete_record(vacuum,d)
    B=1.25*abs(alpha);bound=2*B*B/N*np.exp(2*B/N)
    observed=.5*sum(abs(records[k]-average[k]) for k in records)
    # matrix_power accumulation is numerical evidence; the analytic channel
    # inequality below is independent of these floating diagnostics.
    assert observed<bound+1e-9
    # Independent standard active three-flavour matter law, derived after
    # taking the complete relativistic/chirality comparison specified in proof.
    R=d['R'];m=d['masses'];sign=-1 if antineutrino else 1
    Hstd=np.diag((m*m-m[0]**2)/(2*d['p'])*u['time_eV_inverse'])+sign*alpha/2*np.outer(R[0].conj(),R[0])
    std=R@expm(-1j*Hstd)@R.conj().T;stdprob=abs(std[:,1])**2
    active_diff=max(abs(average[k]-stdprob[j]) for j,k in enumerate(('e','mu','tau')))
    chirality=(max(m)/(2*d['p']))**2;comparison=20*chirality*(1+max(abs(phase))+abs(alpha))
    assert active_diff<comparison+2e-14
    dephased=sum(K@mean@K.conj().T for K in Ks)
    return {'momentum_GeV':momentum_GeV,'antineutrino':antineutrino,'units':u,'collisions':N,'complete_finite_medium_records':records,'complete_coherent_mean_records':average,'complete_vacuum_records':vac,'finite_medium_vs_mean_TV_diagnostic':observed,'uniform_complete_channel_TV_bound':bound,'standard_active_matter_probabilities_diagnostic':stdprob.tolist(),'mean_vs_active_matter_max_error_diagnostic':active_diff,'relativistic_chirality_record_comparison_bound':comparison,'actual_native_forward_norm':float(np.linalg.norm(d['M'],2)),'full_hilbert_dimension_active_plus_sterile_helicities':6,'electron_spin_record_paths':'4^N Kraus words retained as a complete composed instrument; evaluation uses its36-dimensional superoperator','source_work':{'per_collision_external_switching_absolute_upper_eV':2*units()['native_Weyl_G_GeV_minus2']*1e-18*u['electron_density_eV3'],'mean_neutral_free_energy_change_eV':float(np.trace(np.diag(d['free'])@(rho-initial)).real),'mean_source_energy_is_conserved_between_switches':True}}


def certificate():
    native=native_forward();rows=[run(E,L,rho,N,anti) for E,L,rho in ((.6,295.,2.6),(2.,810.,2.8),(2.5,1300.,2.8)) for N in (2**10,2**14) for anti in (False,True)]
    u=unit_dictionary();c=calibration();dm=float(c['mass_squared_eV2_exact'][2]-c['mass_squared_eV2_exact'][0]);cos2=1-2*float(c['oscillation_inputs']['sin2_theta13']);resonance=dm*cos2/(2*u['charged_current_potential_eV'])/1e9
    return {'schema':'oph.native_neutrino_electron_medium.v1','verdict':'PASS_NATIVE_MATTER_POTENTIAL_AND_COMPLETE_COLLISION_RECORDS','native_forward_current_exact_matrix':native['exact_matrix'],'native_forward_characteristic_polynomial':native['exact_polynomial'],'exact_forward_current_norm':'1','exact_unpolarized_mean':'(P_e/2-I/4) tensor I_spin, before massive active-band compression','physical_charged_current_potential':'G_Weyl*n_e/2=sqrt(2)*GF*n_e','same_mass_tensor':'all three calibrated Dirac masses, full family frame, both helicities and their sterile readout complements','finite_medium_executions':rows,'two_flavour_atmospheric_resonance_scale_GeV':resonance,'resonance_scope':'diagnostic two-flavour scale of the derived potential; full records use all three flavours and the calibrated CP phase','medium_source_contract':'Native positive-frequency forward-current port architecture. Exact number-conserving quartic forward terms retain the complete6x2 admitted sector and its Fock complement. This is not an unproved projection equivalence for the unfiltered momentum-changing contact source.','apparatus':'fresh unpolarized rest-electron spin ports; every spin/environment word retained. Density, composition, flight length and collision count are declared apparatus inputs.','empirical_ancestry':'six neutral oscillation inputs and chosen lightest mass; GF weak normalization with muon ancestry; no matter-oscillation probability or resonance is fitted','physical_beam_event_likelihood_or_autonomous_released_neutrino_claimed':False,'scope':'A derived same-current matter potential and controlled finite complete observer instruments, with physical unit dictionary; no measured matter-event fit or full unfiltered scattering theorem.'}
